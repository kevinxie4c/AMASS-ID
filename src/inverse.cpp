#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <chrono>
#include <Eigen/Core>
#include <getopt.h>
#include <dart/dart.hpp>
#include <Python.h>
#include "SimCharacter.h"
#include "PyUtil.h"
#include "IOUtil.h"
#include "c_butterworth/c_butterworth_types.h"
#include "c_butterworth/c_butterworth.h"
#include "c_butterworth/rt_nonfinite.h"

#define ALGLIB	0
#define MOSEK	1
#define OPTIMIZOR   MOSEK

#if OPTIMIZOR == ALGLIB
#include "ext/optimization.h"
#else
#include <mosek.h>
#endif

using namespace std;
using namespace Eigen;
using namespace dart::dynamics;
using namespace dart::simulation;
#if OPTIMIZOR == ALGLIB
using namespace alglib;
#endif

void printUsage(char * prgname)
{
    cout << "usage: " << prgname << " [options] pose_file contact_file" << endl;
    cout << endl;
    cout << "options:" << endl;
    cout << "-j --char_file=string" << endl;
    cout << "-f --frame_time=double" << endl;
    cout << "-c --cutoff_freq=double" << endl;
    cout << "-g --ground_offset=double" << endl;
    cout << "-o --outdir=string" << endl;
}

void mosekOK(MSKrescodee r)
{
    if (r != MSK_RES_OK)
    {
	char symname[MSK_MAX_STR_LEN];
	char desc[MSK_MAX_STR_LEN];

	MSK_getcodedesc(r, symname, desc);
	cout << "Error " << symname << " - '" << desc << "'\n" << endl;
	exit(0);
    }
}

static void MSKAPI printstr(void *handle, const char str[])
{
    cout << str << endl;
}


int main(int argc, char* argv[])
{
    string jsonFilename = "data/character.json";
    string poseFilename;
    string contactNodesFilename;
    string outdir = "output";
    double frameTime = 1.0 / 120.0, cutoffFreq = 2, groundOffset = 0;

    while (1)
    {
	int c;
	static struct option long_options[] =
	{
	    { "char_file", required_argument, NULL, 'j' },
	    { "frame_time", required_argument, NULL, 'f' },
	    { "cutoff_freq", required_argument, NULL, 'c' },
	    { "ground_offset", required_argument, NULL, 'g' },
	    { "outdir", required_argument, NULL, 'o' },
	    { 0, 0, 0, 0 }
	};
	int option_index = 0;

	c = getopt_long(argc, argv, "j:f:c:g:o:", long_options, &option_index);
	if (c == -1)
        break;

	switch (c)
	{
	    case 'j':
		jsonFilename = optarg;
		break;
	    case 'f':
		frameTime = stod(optarg);
		break;
	    case 'c':
		cutoffFreq = stod(optarg);
		break;
	    case 'g':
		groundOffset = stod(optarg);
		break;
	    case 'o':
		outdir = optarg;
		break;
	    default:
		printUsage(argv[0]);
		exit(0);
	}
    }

    if (optind + 2 == argc)
    {
	poseFilename = argv[optind];
	contactNodesFilename = argv[optind + 1];
    }
    else
    {
	printUsage(argv[0]);
	exit(0);
    }

    SimCharacter character(jsonFilename);
    WorldPtr world = World::create();
    SkeletonPtr &skeleton = character.skeleton;

    std::vector<dart::dynamics::BodyNode*> bns = skeleton->getBodyNodes();

    Py_Initialize();
    //PySys_SetArgv(argc, argv);
    PyRun_SimpleString("import sys\n"
		       "sys.path.append(\".\")\n"
                       "print('python version:', sys.version)\n");

    PyObject *pDict = load_npz(poseFilename);
    //print_py_obj(pDict);
    PyObject *npa;
    npa = PyMapping_GetItemString(pDict, "trans");
    //print_py_obj(npa);
    Eigen::MatrixXd trans = npa2mat(npa);
    //Py_DECREF(npa);
    //std::cout << trans.row(0) << std::endl;
    npa = PyMapping_GetItemString(pDict, "poses");
    Eigen::MatrixXd poses_orig = npa2mat(npa);
    //Py_DECREF(npa);
    //std::cout << poses_orig.row(0) << std::endl;
    Py_DECREF(pDict);

    std::vector<size_t> indices{0, 3, 12, 21, 30, 6, 15, 24, 33, 9, 18, 27, 36, 45, 39, 48, 54, 60, 42, 51, 57, 63};
    Eigen::MatrixXd poses(poses_orig.rows(), indices.size() * 3);
    for (size_t i = 0; i < indices.size(); ++i)
	poses.middleCols(i * 3, 3) = poses_orig.middleCols(indices[i], 3);

    Eigen::MatrixXd position_mat = Eigen::MatrixXd(poses.rows(), poses.cols() + 3);
    position_mat.leftCols(3) = poses.leftCols(3);
    position_mat.middleCols(3, 3) = trans;
    position_mat.rightCols(poses.cols() - 3) = poses.rightCols(poses.cols() - 3);

    world->addSkeleton(skeleton);
    vector<VectorXd> positions;
    vector<VectorXd> velocities;
    vector<VectorXd> accelerations;
    //skeleton->setGravity(Vector3d(0, -9.8, 0));
    skeleton->setGravity(Vector3d(0, 0, -9.8));
    for (size_t i = 0; i < position_mat.rows(); ++i)
	positions.push_back(position_mat.row(i));
    for (size_t i = 0; i < positions.size() - 1; ++i)
	velocities.push_back(skeleton->getPositionDifferences(positions[i + 1], positions[i]) / frameTime);

    coder::array<double, 1> x, y;
    vector<VectorXd> tmp(velocities);
    x.set_size(velocities.size());
    size_t ndof = skeleton->getDofs().size();
    for (size_t i = 0; i < ndof; ++i)
    { for (int j = 0; j < x.size(0); ++j)
	    x[j] = velocities[j][i];
	c_butterworth(x, cutoffFreq, 100, y);
	for (int j = 0; j < x.size(0); ++j)
	    velocities[j][i] = y[j];
    }
    c_butterworth_terminate();

    for (size_t i = 0; i < velocities.size() - 1; ++i)
	accelerations.push_back(skeleton->getVelocityDifferences(velocities[i + 1], velocities[i]) / frameTime);

    vector<vector<string>> contactNodes;
    vector<vector<Vector3d>> contactPoints;
    ifstream input(contactNodesFilename);
    string line;
    while (getline(input, line))
    {
	stringstream strStream(line);
	vector<string> nodeList;
	vector<Vector3d> pointList;
	string s;
	size_t n;
	while (strStream >> s >> n)
	{
	    double x, y, z;
	    for (size_t i = 0; i < n; ++i)
	    {
		strStream >> x >> y >> z;
		nodeList.push_back(s);
		pointList.push_back(Vector3d(x, y, z));
	    }
	}
	contactNodes.push_back(nodeList);
	contactPoints.push_back(pointList);
    }
    input.close();

    ofstream pout(outdir + "/positions.txt");
    ofstream fout(outdir + "/forces.txt");
    ofstream vout(outdir + "/velocities.txt");
    ofstream aout(outdir + "/accelerations.txt");
    ofstream cfout(outdir + "/contact_forces.txt");
    ofstream eout(outdir + "/errors.txt");

    double mu = 1;
    double reg = 10;
    MSKenv_t env = NULL;
    mosekOK(MSK_makeenv(&env, NULL));
    
    //for (size_t i = 0; i < accelerations.size(); ++i)
    for (size_t i = 0; i < 500; ++i)
    {
	cout << "frame " << i << endl;
	skeleton->setPositions(positions[i]);
	skeleton->setVelocities(velocities[i]);
	skeleton->setAccelerations(accelerations[i]);
	//VectorXd v = VectorXd::Zero(ndof);
	//skeleton->setPositions(v);
	//skeleton->setVelocities(v);
	//skeleton->setAccelerations(v);

	MatrixXd M = skeleton->getMassMatrix();
	MatrixXd C = skeleton->getCoriolisAndGravityForces();
	//std::cout << C << std::endl;
	//MatrixXd G = skeleton->getGravityForces();
	//std::cout << G << std::endl;

	// TODO: confirm the following computation of dL is correct
	//Eigen::Vector3d dL = Eigen::Vector3d::Zero(); // rate of angular momentum at center of mass
	//Eigen::Vector3d com = skeleton->getCOM();
	//for (const BodyNode *bn: skeleton->getBodyNodes())
	//{
	//    dL += (bn->getCOM() - com).cross(bn->getMass() * bn->getCOMLinearAcceleration());

	//    Eigen::Matrix3d R = bn->getTransform().rotation();
	//    Eigen::Matrix3d I = bn->getInertia().getMoment(); // moment of inertia in local frame
	//    Eigen::Vector3d omiga = bn->getCOMSpatialVelocity().head(3); // angular velocity in local frame
	//    Eigen::Vector3d dOmiga = bn->getCOMSpatialAcceleration().head(3); // angular acceleration in local frame

	//    dL += R * (I * dOmiga + omiga.cross(I * omiga));
	//}

	VectorXd qddot = accelerations[i];

	
	if (contactNodes[i].size() > 0)
	{
	    std::chrono::steady_clock::time_point time_p = std::chrono::steady_clock::now();
	    size_t m = ndof;
	    size_t n = contactNodes[i].size() * 5;
	    /*
	     * Optimize contact force lambda by minimizing
	     *	    ||(M * qddot + C - Jt * B1 * B2 * lambda)_1:6||
	     * By computing the L2 norm, this becomes a QP problem:
	     *	    min 1/2 x^T * A * x - b^T * x
	     * where
	     *	    A = A1^T * A1
	     *	    b = a^T * A1
	     *	    A1 = Jt * B1 * B2
	     *	    a = M * qddot + C
	     */
	    MatrixXd emA(m, n);	// A1
	    vector<MatrixXd> B2list;
	    for (size_t j = 0; j < contactNodes[i].size(); ++j)
	    {
		const BodyNode *bn = skeleton->getBodyNode(contactNodes[i][j]);
		const Vector3d &point = contactPoints[i][j];
		Isometry3d transform = bn->getWorldTransform();
		MatrixXd Jt = skeleton->getWorldJacobian(bn, transform.inverse() * point).transpose();
		MatrixXd B1(6, 3);
		//B1.topRows(3) = dart::math::makeSkewSymmetric(point);
		//B1.topRows(3) = dart::math::makeSkewSymmetric(transform.inverse() * point);
		//B1.topRows(3) = dart::math::makeSkewSymmetric(point - transform.translation());
		B1.topRows(3) = Matrix3d::Zero();
		B1.bottomRows(3) = Matrix3d::Identity();
		Vector3d normal = Vector3d::UnitZ();
		Vector3d tangent1 = Vector3d::UnitX();
		Vector3d tangent2 = Vector3d::UnitY();
		MatrixXd B2(3, 5);
		B2.col(0) = normal;
		B2.col(1) = tangent1;
		B2.col(2) = tangent2;
		B2.col(3) = -tangent1;
		B2.col(4) = -tangent2;
		B2list.push_back(B2);
		MatrixXd mat = Jt * B1 * B2;
		emA.middleCols(j * 5, 5) = mat;
	    }
	    //qddot = VectorXd::Zero(ndof);
	    //qddot[5] = -9.8;
	    VectorXd a = (M * qddot + C);
	    //std::cout << "M" << std::endl;
	    //std::cout << M << std::endl;
	    //std::cout << "qddot" << std::endl;
	    //std::cout << qddot << std::endl;
	    //std::cout << "M * qddot" << std::endl;
	    //std::cout << M * qddot << std::endl;
	    //std::cout << "C" << std::endl;
	    //std::cout << C << std::endl;
	    //std::cout << "a" << std::endl;
	    //std::cout << a << std::endl;
	    VectorXd a6 = a.topRows(6);
	    MatrixXd emA6 = emA.topRows(6);
	    MatrixXd AtA = emA6.transpose() * emA6;
	    MatrixXd atA = a6.transpose() * emA6;
	    AtA += MatrixXd::Identity(n, n) * reg / n;
	    cout << "matrix construction (Eigen): " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - time_p).count() << " ms" << endl;
	    time_p = std::chrono::steady_clock::now();

#if OPTIMIZOR == ALGLIB
	    real_2d_array A;
	    //sparsematrix A;
	    real_1d_array b;
	    real_1d_array s;
	    A.setlength(n, n);
	    //sparsecreate(n, n, 0, A);
	    b.setlength(n);
	    s.setlength(n);
	    for (size_t i = 0; i < n; ++i)
	    {
		for (size_t j = 0; j < n; ++j)
		    A[i][j] = AtA(i, j);
		//for (size_t j = i; j < n; ++j)
		//    if (abs(AtA(i, j)) > 1e-8)
		//	sparseset(A, i, j, AtA(i, j));
		b[i] = -atA(0, i);
		s[i] = 1;
	    }
	    //real_2d_array c;
	    sparsematrix c;
	    integer_1d_array ct;
	    //c.setlength(4 * contactNodes[i].size(), n + 1);
	    sparsecreate(4 * contactNodes[i].size(), n + 1, 0, c);
	    ct.setlength(4 * contactNodes[i].size());
	    //for (size_t j = 0; j < 4 * contactNodes[i].size(); ++j)
	    //    for (size_t k = 0; k < n + 1; ++k)
	    //    c[j][k] = 0;
	    for (size_t j = 0; j < contactNodes[i].size(); ++j)
	    {
		for (size_t k = 0; k < 4; ++k)
		{
		    //c[4 * j + k][5 * j + 0] = mu;
		    //c[4 * j + k][5 * j + k + 1] = -1;
		    sparseset(c, 4 * j + k, 5 * j + 0, mu);
		    sparseset(c, 4 * j + k, 5 * j + k + 1, -1);

		    ct[4 * j + k] = 1;
		}
	    }
	    real_1d_array bndl;
	    real_1d_array bndu;
	    bndl.setlength(n);
	    bndu.setlength(n);
	    for (size_t i = 0; i < n; ++i)
	    {
		bndl[i] = 0;
		bndu[i] = 1e5;
	    }
	    real_1d_array x;
	    minqpstate state;
	    minqpreport rep;
	    minqpcreate(n, state);
	    minqpsetquadraticterm(state, A);
	    //minqpsetquadratictermsparse(state, A, true);
	    minqpsetlinearterm(state, b);
	    //cout << "c" << endl;
	    //cout << c.tostring(2) << endl;
	    //cout << "ct" << endl;
	    //cout << ct.tostring() << endl;
	    //minqpsetlc(state, c, ct);
	    minqpsetlcsparse(state, c, ct, 4 * contactNodes[i].size());
	    minqpsetbc(state, bndl, bndu);
	    minqpsetscale(state, s);
	    minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
	    cout << "matrix construction (Alglib): " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - time_p).count() << " ms" << endl;
	    minqpoptimize(state);
	    minqpresults(state, x, rep);

#else
	    MSKtask_t task = NULL;
	    size_t num_var = n;
	    size_t num_con = 4 * contactNodes[i].size();
	    mosekOK(MSK_maketask(env, num_con, num_var, &task));
	    mosekOK(MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr));
	    mosekOK(MSK_appendcons(task, num_con));
	    mosekOK(MSK_appendvars(task, num_var));

	    // set up lower triangular part of "A" ("Q" in Mosek)
	    vector<MSKint32t> qsubi, qsubj;
	    vector<double> qval;
	    for (size_t i = 0; i < n; ++i)
	    {
	        for (size_t j = 0; j <= i; ++j)
	        {
	            if (abs(AtA(i, j)) > 1e-8)
	            {
	        	qsubi.push_back(i);
	        	qsubj.push_back(j);
	        	qval.push_back(AtA(i, j));
	        	//mosekOK(MSK_putqobjij(task, i, j, AtA(i, j)));
	            }
	        }
	    }
	    mosekOK(MSK_putqobj(task, qval.size(), qsubi.data(), qsubj.data(), qval.data()));
	    for (size_t j = 0; j < num_var; ++j)
	    {
		mosekOK(MSK_putcj(task, j, -atA(0, j))); // linear term "b" ("c" in Mosek)
		mosekOK(MSK_putvarbound(task, j, MSK_BK_LO, 0.0, +MSK_INFINITY));
	    }
	    // set up constraints
	    for (size_t j = 0; j < contactNodes[i].size(); ++j)
	    {
		for (size_t k = 0; k < 4; ++k)
		{
		    mosekOK(MSK_putaij(task, 4 * j + k, 5 * j + 0, mu));
		    mosekOK(MSK_putaij(task, 4 * j + k, 5 * j + k + 1, -1));
		    mosekOK(MSK_putconbound(task, 4 * j + k, MSK_BK_LO, 0.0, +MSK_INFINITY));
		}
	    }
	    cout << "matrix construction (Mosek): " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - time_p).count() << " ms" << endl;
	    MSKrescodee trmcode;
	    mosekOK(MSK_optimizetrm(task, &trmcode));
	    MSK_solutionsummary(task, MSK_STREAM_MSG);
	    MSKsolstae solsta;
	    MSK_getsolsta(task, MSK_SOL_ITR, &solsta);
	    vector<double> x(num_var);
	    switch (solsta)
	    {
		case MSK_SOL_STA_OPTIMAL:
		    MSK_getxx(task, MSK_SOL_ITR, x.data());
		    break;

		case MSK_SOL_STA_DUAL_INFEAS_CER:
		case MSK_SOL_STA_PRIM_INFEAS_CER:
		    cout << "Primal or dual infeasibility certificate found." << endl;
		    break;

		case MSK_SOL_STA_UNKNOWN:
		    cout << "The status of the solution could not be determined. Termination code: " << trmcode << endl;
		    break;

		default:
		    cout << "Other solution status." << endl;
		    break;
	    }
	    MSK_deletetask(&task);

#endif
	    VectorXd lambda(n);
	    for (size_t i = 0; i < n; ++i)
		lambda[i] = x[i];

	    //cout << "lambda" << endl;
	    //cout << lambda << endl;
	    //cout << "error1 = " << (lambda.transpose() * AtA * lambda - atA * lambda * 2 + a6.transpose() * a6).cwiseSqrt() << endl; 
	    VectorXd Q = a - emA * lambda;
	    //cout << "error2 = " << Q.head(6).norm() << endl; 
	    eout << Q.head(6).norm() << endl; 
	    fout << Q.transpose() << endl;
	    for (size_t j = 0; j < contactNodes[i].size(); ++j)
	    {
		VectorXd lamb = lambda.segment(j * 5, 5);
		Vector3d f = B2list[j] * lamb;
		cfout << contactPoints[i][j].transpose() << " " << f.transpose() << " ";
	    }
	    cfout << endl;
	}
	else
	{
	    VectorXd Q = M * qddot + C;
	    //Q.head(6) = VectorXd::Zero(6);
	    fout << Q.transpose() << endl;
	    cfout << endl;
	}

	pout << positions[i].transpose() << endl;
	vout << velocities[i].transpose() << endl;
	aout << accelerations[i].transpose() << endl;
    }
    MSK_deleteenv(&env);
    pout.close();
    fout.close();
    vout.close();
    aout.close();
    cfout.close();
    eout.close();

    return 0;
}
