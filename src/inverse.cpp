#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <chrono>
#include <cstring>
#include <cassert>
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
#define FULL_SOLVE
#define CONIC_OPTIMIZATION

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
    cout << "\t-j, --char_file=string" << endl;
    cout << "\t-f, --frame_time=double" << endl;
    cout << "\t-c, --cutoff_freq=double" << endl;
    cout << "\t-u, --mu=double" << endl;
    cout << "\t-s, --step_length=int" << endl;
    cout << "\t-r, --regularization=double" << endl;
    cout << "\t-n, --num_frame=int" << endl;
    cout << "\t-o, --outdir=string" << endl;
    cout << "\t-A, --start_frame=int" << endl;
    cout << "\t-E, --end_frame=int" << endl;
    cout << "\t-F, --filter_type=none|position|velocity" << endl;
    cout << "\t-R, --root_rot_weight" << endl;
    cout << "\t-T, --root_trans_weight" << endl;
    cout << "\t-S, --use_sim_state" << endl;
}

#if OPTIMIZOR == MOSEK
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
#endif


int main(int argc, char* argv[])
{
    string jsonFilename = "data/character.json";
    string poseFilename;
    string contactNodesFilename;
    string outdir = "output";
    double frameTime = 1.0 / 120.0, cutoffFreq = 2, mu = 1, reg = 0.0001;
    size_t numFrame = 1;
    size_t stepLength = 1;
    size_t startFrame = 0, endFrame = 0;
    bool endFrameSet = false;
    int filterType = 2;	// 0: none; 1: position; 2: velocity.
    bool useSimState = false;
    double rootRotWeight = 0;
    double rootTransWeight = 0;

    while (1)
    {
	int c;
	static struct option long_options[] =
	{
	    { "char_file", required_argument, NULL, 'j' },
	    { "frame_time", required_argument, NULL, 'f' },
	    { "cutoff_freq", required_argument, NULL, 'c' },
	    { "mu", required_argument, NULL, 'u' },
	    { "step_length", required_argument, NULL, 's' },
	    { "regularization", required_argument, NULL, 'r' },
	    { "num_frame", required_argument, NULL, 'n' },
	    { "outdir", required_argument, NULL, 'o' },
	    { "start_frame", required_argument, NULL, 'A' },
	    { "end_frame", required_argument, NULL, 'E' },
	    { "filter_type", required_argument, NULL, 'F' },
	    { "root_rot_weight", required_argument, NULL, 'R' },
	    { "root_trans_weight", required_argument, NULL, 'T' },
	    { "use_sim_state", no_argument, NULL, 0 },
	    { 0, 0, 0, 0 }
	};
	int option_index = 0;

	c = getopt_long(argc, argv, "j:f:c:u:s:r:n:o:A:E:F:R:T:S", long_options, &option_index);
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
	    case 'u':
		mu = stod(optarg);
		break;
	    case 's':
		stepLength = stoi(optarg);
		break;
	    case 'r':
		reg = stod(optarg);
		break;
	    case 'n':
		numFrame = stoi(optarg);
		break;
	    case 'o':
		outdir = optarg;
		break;
	    case 'A':
		startFrame = stoi(optarg);
		break;
	    case 'E':
		endFrame = stoi(optarg);
		endFrameSet = true;
		break;
	    case 'F':
		if (strcmp(optarg, "none") == 0)
		    filterType = 0;
		else if (strcmp(optarg, "position") == 0)
		    filterType = 1;
		else if (strcmp(optarg, "velocity") == 0)
		    filterType = 2;
		else
		{
		    cout << "unknown filter_type: " << optarg << endl;
		    exit(0);
		}
		break;
	    case 'R':
		rootRotWeight = stod(optarg);
		break;
	    case 'T':
		rootTransWeight = stod(optarg);
		break;
	    case 'S':
		useSimState = true;
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
    //skeleton->setGravity(Vector3d(0, -9.8, 0));   // y-axis points up
    skeleton->setGravity(Vector3d(0, 0, -9.8));	    // z-axis points up

    for (size_t i = 0; i < position_mat.rows(); ++i)
	positions.push_back(position_mat.row(i));

    size_t ndof = skeleton->getDofs().size();

    if (filterType == 1)
    {
	coder::array<double, 1> x, y;
	vector<VectorXd> tmp(positions);
	x.set_size(positions.size());
	for (size_t i = 0; i < ndof; ++i)
	{ for (int j = 0; j < x.size(0); ++j)
		x[j] = positions[j][i];
	    c_butterworth(x, cutoffFreq, 100, y);
	    for (int j = 0; j < x.size(0); ++j)
		positions[j][i] = y[j];
	}
	c_butterworth_terminate();
    }

    for (size_t i = 0; i < positions.size() - stepLength; ++i)
	velocities.push_back(skeleton->getPositionDifferences(positions[i + stepLength], positions[i]) / (stepLength * frameTime));

    if (filterType == 2)
    {
	coder::array<double, 1> x, y;
	vector<VectorXd> tmp(velocities);
	x.set_size(velocities.size());
	for (size_t i = 0; i < ndof; ++i)
	{ for (int j = 0; j < x.size(0); ++j)
		x[j] = velocities[j][i];
	    c_butterworth(x, cutoffFreq, 100, y);
	    for (int j = 0; j < x.size(0); ++j)
		velocities[j][i] = y[j];
	}
	c_butterworth_terminate();
    }

    for (size_t i = 0; i < velocities.size() - stepLength; ++i)
	accelerations.push_back(skeleton->getVelocityDifferences(velocities[i + stepLength], velocities[i]) / (stepLength * frameTime));

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

    /* Jacobian test
    SimCharacter acrobot("Acrobot.json");
    acrobot.skeleton->setPositions(Vector2d(1.57079632679, 0));
    Vector2d dq = Vector2d(1, 2);
    acrobot.skeleton->setVelocities(dq);
    BodyNodePtr bn = acrobot.skeleton->getBodyNode(1);
    MatrixXd J = acrobot.skeleton->getWorldJacobian(bn, Vector3d(0, -2, 0));
    cout << "J * dq" << endl;
    cout << J * dq << endl;
    cout << "sv" << endl;
    cout << bn->getSpatialVelocity() << endl;
    cout << "T.t" << endl;
    cout << bn->getWorldTransform().translation() << endl;
    cout << "T.R" << endl;
    cout << bn->getWorldTransform().linear() << endl;
    return 0;
    */

    SkeletonPtr kin_skeleton = skeleton->cloneSkeleton();
    size_t f_start = startFrame, f_end = endFrameSet ? endFrame : accelerations.size();
    vector<vector<Vector3d>> transformedContactPoints(contactPoints.size());
#ifdef FULL_SOLVE
    VectorXd pos = positions[f_start];
    VectorXd vel = velocities[f_start];
    VectorXd vel_n;
    VectorXd vel_hat;
    MatrixXd W = MatrixXd::Identity(ndof, ndof);
    if (useSimState)
    {
	for (int i = 0; i < 3; ++i)
	{
	    W(i, i) = rootRotWeight;
	    W(3 + i, 3 + i) = rootTransWeight;
	}
    }
#endif

#if OPTIMIZOR == MOSEK
    MSKenv_t env = NULL;
    mosekOK(MSK_makeenv(&env, NULL));
#endif
    
    //for (size_t i = 0; i < accelerations.size(); ++i)
    for (size_t i = f_start; i < f_end; i += stepLength)
    {
	cout << "frame " << i << endl;
	kin_skeleton->setPositions(positions[i]);
#ifdef FULL_SOLVE
	if (!useSimState)
	{
	    pos = positions[i];
	    vel = velocities[i];
	    vel_hat = velocities[i + 1];
	}
	else
	{
	    vel_hat = skeleton->getPositionDifferences(positions[i + stepLength], pos) / (stepLength * frameTime);
	    //skeleton->setVelocities(vel);
	    //skeleton->integratePositions(frameTime);
	    //VectorXd pos_n = skeleton->getPositions();
	    //vel_hat = skeleton->getPositionDifferences(positions[i + 2], pos_n) / frameTime;
	}
	skeleton->setPositions(pos);
	skeleton->setVelocities(vel);
#else
	skeleton->setPositions(positions[i]);
	skeleton->setVelocities(velocities[i]);
	//skeleton->setAccelerations(accelerations[i]);
	//VectorXd v = VectorXd::Zero(ndof);
	//skeleton->setPositions(v);
	//skeleton->setVelocities(v);
	//skeleton->setAccelerations(v);
#endif

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

#ifdef FULL_SOLVE
	std::chrono::steady_clock::time_point time_p = std::chrono::steady_clock::now();
	if (numFrame == 1)
	//if (false)
	{
	    size_t m = ndof;
	    //size_t n = contactNodes[i].size() * 5 + ndof - 6;
	    size_t n_contacts = contactNodes[i].size() * 5;
	    size_t n = contactNodes[i].size() * 5 + ndof;
	    double h = stepLength * frameTime;
	    //double h = 1;
	    MatrixXd H(m, n);
	    MatrixXd hMinv = h * M.inverse();
	    vector<MatrixXd> B2list;
	    MatrixXd JtB(m, contactNodes[i].size() * 5);
	    if (contactNodes[i].size() > 0)
	    {
		//MatrixXd JtB(m, contactNodes[i].size() * 5);
		for (size_t j = 0; j < contactNodes[i].size(); ++j)
		{
		    const BodyNode *bn = skeleton->getBodyNode(contactNodes[i][j]);
		    const Vector3d &point = contactPoints[i][j];
		    //Isometry3d transform = bn->getWorldTransform();
		    Isometry3d transform_kin = kin_skeleton->getBodyNode(contactNodes[i][j])->getWorldTransform();
		    Isometry3d transform_sim = skeleton->getBodyNode(contactNodes[i][j])->getWorldTransform();
		    transformedContactPoints[i].push_back(transform_sim * transform_kin.inverse() * point);
		    //MatrixXd Jt = skeleton->getWorldJacobian(bn, transform.inverse() * point).transpose();
		    MatrixXd Jt = skeleton->getWorldJacobian(bn).transpose();
		    MatrixXd B1(6, 3);
		    //B1.topRows(3) = dart::math::makeSkewSymmetric(point);
		    //B1.topRows(3) = dart::math::makeSkewSymmetric(transform.inverse() * point);
		    B1.topRows(3) = dart::math::makeSkewSymmetric(point - transform_kin.translation());
		    //B1.topRows(3) = Matrix3d::Zero();
		    B1.bottomRows(3) = Matrix3d::Identity();
		    //Vector3d normal = Vector3d::UnitY();	// y-axis points up
		    //Vector3d tangent1 = Vector3d::UnitX();
		    //Vector3d tangent2 = Vector3d::UnitZ();
		    Vector3d normal = Vector3d::UnitZ();	// z-axis points up
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
		    JtB.middleCols(j * 5, 5) = mat;
		}
		H.leftCols(contactNodes[i].size() * 5) = hMinv * JtB;
	    }
	    //H.rightCols(ndof - 6) = hMinv.rightCols(ndof - 6);
	    H.rightCols(ndof) = hMinv.rightCols(ndof);
	    VectorXd hMinvC = hMinv * C;
	    VectorXd d = vel - vel_hat - hMinvC;

	    MSKtask_t task = NULL;
#ifdef CONIC_OPTIMIZATION
	    size_t num_var = n + m + (n_contacts + 6) + 1; // x = (lambda, xi_1 xi_2, t)
	    size_t num_con = 4 * contactNodes[i].size() + m + (n_contacts + 6);
	    mosekOK(MSK_maketask(env, num_con, num_var, &task));
	    mosekOK(MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr));
	    mosekOK(MSK_appendcons(task, num_con));
	    mosekOK(MSK_appendvars(task, num_var));

	    // set up linear term in the objective
	    for (size_t j = 0; j < num_var - 1; ++j)
		mosekOK(MSK_putcj(task, j, 0.0));
	    mosekOK(MSK_putcj(task, num_var - 1, 1.0));

	    // box constraints
	    for (size_t j = 0; j < contactNodes[i].size() * 5; ++j)
		mosekOK(MSK_putvarbound(task, j, MSK_BK_LO, 0.0, +MSK_INFINITY));
	    for (size_t j = contactNodes[i].size() * 5; j < num_var; ++j)
		mosekOK(MSK_putvarbound(task, j, MSK_BK_FR, -MSK_INFINITY, +MSK_INFINITY));

	    // linear constraints
	    vector<MSKint32t> asubi, asubj;
	    vector<double> aval;
	    for (size_t j = 0; j < contactNodes[i].size(); ++j)
	    {
		for (size_t k = 0; k < 4; ++k)
		{
		    //mosekOK(MSK_putaij(task, 4 * j + k, 5 * j + 0, mu));
		    asubi.push_back(4 * j + k);
		    asubj.push_back(5 * j + 0);
		    aval.push_back(mu);
		    //mosekOK(MSK_putaij(task, 4 * j + k, 5 * j + k + 1, -1));
		    asubi.push_back(4 * j + k);
		    asubj.push_back(5 * j + k + 1);
		    aval.push_back(-1.0);
		    mosekOK(MSK_putconbound(task, 4 * j + k, MSK_BK_LO, 0.0, +MSK_INFINITY));
		}
	    }
	    for (size_t j = 0; j < m; ++j)
	    {
		for (size_t k = 0; k < n; ++k)
		    if (abs(H(j, k)) > 1e-8)
		    {
			//mosekOK(MSK_putaij(task, 4 * contactNodes[i].size() + j, k, -H(j, k)));
			asubi.push_back(4 * contactNodes[i].size() + j);
			asubj.push_back(k);
			aval.push_back(-H(j, k));
		    }
		//mosekOK(MSK_putaij(task, 4 * contactNodes[i].size() + j, n + j, 1.0));
		asubi.push_back(4 * contactNodes[i].size() + j);
		asubj.push_back(n + j);
		aval.push_back(1.0);
		mosekOK(MSK_putconbound(task, 4 * contactNodes[i].size() + j, MSK_BK_FX, d[j], d[j]));
	    }
	    for (size_t j = 0; j < n_contacts + 6; ++j)
	    {
		asubi.push_back(4 * contactNodes[i].size() + m + j);
		asubj.push_back(j);
		if (j < n_contacts)
		    aval.push_back(sqrt(reg / n_contacts));
		else if (j < n_contacts + 3)
		    aval.push_back(sqrt(rootRotWeight));
		else
		    aval.push_back(sqrt(rootTransWeight));
		asubi.push_back(4 * contactNodes[i].size() + m + j);
		asubj.push_back(n + m + j);
		aval.push_back(-1.0);
		mosekOK(MSK_putconbound(task, 4 * contactNodes[i].size() + m + j, MSK_BK_FX, 0.0, 0.0));
	    }
	    mosekOK(MSK_putaijlist(task, aval.size(), asubi.data(), asubj.data(), aval.data()));

	    // conic constraints
	    vector<MSKint32t> csub(m + (n_contacts + 6) + 1);
	    csub[0] = n + m + (n_contacts + 6); // t
	    for (size_t i = 0; i < m + (n_contacts + 6); ++i)
		csub[i + 1] = n + i;
	    mosekOK(MSK_appendcone(task, MSK_CT_QUAD, 0.0, csub.size(), csub.data()));
#else	// CONIC_OPTIMIZATION
	    MatrixXd A = H.transpose() * W * H;
	    //A += MatrixXd::Identity(n, n) * reg / n;    // make it positive definite and evenly "spread" the forces
	    for (size_t j = 0; j < contactNodes[i].size() * 5; ++j)
		A(j, j) += reg / (contactNodes[i].size() * 5);
	    VectorXd b = H.transpose() * W * d;

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
		    //if (abs(A(i, j)) > 1e-8)
		    {
			qsubi.push_back(i);
			qsubj.push_back(j);
			qval.push_back(A(i, j));
			//mosekOK(MSK_putqobjij(task, i, j, A(i, j)));
		    }
		}
	    }
	    mosekOK(MSK_putqobj(task, qval.size(), qsubi.data(), qsubj.data(), qval.data()));
	    for (size_t j = 0; j < num_var; ++j)
		mosekOK(MSK_putcj(task, j, b[j])); // linear term "b" ("c" in Mosek)
	    for (size_t j = 0; j < contactNodes[i].size() * 5; ++j)
		mosekOK(MSK_putvarbound(task, j, MSK_BK_LO, 0.0, +MSK_INFINITY));
	    for (size_t j = contactNodes[i].size() * 5; j < num_var; ++j)
		mosekOK(MSK_putvarbound(task, j, MSK_BK_FR, -MSK_INFINITY, +MSK_INFINITY));
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
#endif	// CONIC_OPTIMIZATION
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

	    VectorXd lambda(n);
	    for (size_t i = 0; i < n; ++i)
		lambda[i] = x[i];
	    //VectorXd xi(m);
	    //for (size_t i = 0; i < m; ++i)
	    //    xi[i] = x[n + i];
	    //cout << "xi " << xi.transpose() << endl;
	    //cout << "t " << x[n + m] << endl;;
	    VectorXd Q = VectorXd::Zero(ndof);
	    //Q.tail(ndof - 6) = lambda.tail(ndof - 6);
	    Q = lambda.tail(ndof);
	    fout << Q.transpose() << endl;
	    for (size_t j = 0; j < contactNodes[i].size(); ++j)
	    {
		VectorXd lamb = lambda.segment(j * 5, 5);
		Vector3d f = B2list[j] * lamb;
		cfout << transformedContactPoints[i][j].transpose() << " " << f.transpose() << " ";
	    }
	    cfout << endl;
	    vel_n = vel - hMinvC + H * lambda;
	    VectorXd err = vel_n - vel_hat;
	    eout << err.norm() << endl;
	    //eout << err.head(6).norm() << endl;
	    cout << "err1 " << err.norm() << endl;
	    //cout << "err2 " << (lambda.transpose() * A * lambda + 2 * b.transpose() * lambda + d.transpose() * d).cwiseSqrt() << endl;
	    cout << "err3 " << (H * lambda + vel - hMinvC - vel_hat).norm() << endl;
	    VectorXd acc = M.colPivHouseholderQr().solve(JtB * lambda.head(contactNodes[i].size() * 5) + Q - C);
	    cout << "err4 " << (velocities[i] + acc * h - velocities[i + 1]).norm() << endl;
	    cout << "err5 " << (velocities[i] + accelerations[i] * h - velocities[i + 1]).norm() << endl;
	    //cout << "acc:\n" << acc << endl;
	    //cout << "acceleration:\n" << accelerations[i] << endl;
	    //cout << "acc_diff " << (acc - accelerations[i]).norm() << endl;
	    //VectorXd _Q = M * (accelerations[i]) + C;
	    //VectorXd _lambda = VectorXd::Zero(n);
	    //_lambda.tail(ndof) = _Q;
	    //cout << "err6 " << (_lambda.transpose() * A * _lambda + 2 * b.transpose() * _lambda + d.transpose() * d).cwiseSqrt() << endl;
	    //cout << "lambda\n" << lambda << endl;
	    //cout << "_lambda\n" << _lambda << endl;
	    //cout << "A\n" << A << endl;
	    //cout << "b\n" << b << endl;
	    //cout << "d\n" << d << endl;
	}
	else
	{
	    vector<MatrixXd> aH;
	    vector<VectorXd> ahMinvC;
	    vector<MatrixXd> B2list;
	    MatrixXd JtB(ndof, contactNodes[i].size() * 5);
	    size_t num_lambda = 0, num_contacts = 0;
	    double h = stepLength;
	    for (size_t j = 0; j < numFrame; ++j)
	    {
		size_t m = ndof;
		size_t n = contactNodes[i + j].size() * 5 + ndof - 6;
		num_lambda += n;
		MatrixXd H(m, n);
		MatrixXd hMinv = h * M.inverse();
		MatrixXd _JtB(m, contactNodes[i + j].size() * 5);
		if (j == 0)
		{
		    skeleton->setPositions(pos);
		    skeleton->setVelocities(vel);
		}
		else
		{
		    skeleton->setPositions(positions[i + j]);
		    skeleton->setVelocities(velocities[i + j]);
		}
		kin_skeleton->setPositions(positions[i + j]);
		kin_skeleton->setVelocities(velocities[i + j]);
		M = skeleton->getMassMatrix();
		C = skeleton->getCoriolisAndGravityForces();
		if (contactNodes[i + j].size() > 0)
		{
		    for (size_t k = 0; k < contactNodes[i + j].size(); ++k)
		    {
			const BodyNode *bn = skeleton->getBodyNode(contactNodes[i + j][k]);
			const Vector3d &point = contactPoints[i + j][k];
			Isometry3d transform_kin = kin_skeleton->getBodyNode(contactNodes[i + j][k])->getWorldTransform();
			Isometry3d transform_sim = skeleton->getBodyNode(contactNodes[i + j][j])->getWorldTransform();
			transformedContactPoints[i].push_back(transform_sim * transform_kin.inverse() * point);
			MatrixXd Jt = skeleton->getWorldJacobian(bn).transpose();
			MatrixXd B1(6, 3);
			B1.topRows(3) = dart::math::makeSkewSymmetric(point - transform_kin.translation());
			B1.bottomRows(3) = Matrix3d::Identity();
			//Vector3d normal = Vector3d::UnitY();	// y-axis points up
			//Vector3d tangent1 = Vector3d::UnitX();
			//Vector3d tangent2 = Vector3d::UnitZ();
			Vector3d normal = Vector3d::UnitZ();	// z-axis points up
			Vector3d tangent1 = Vector3d::UnitX();
			Vector3d tangent2 = Vector3d::UnitY();
			MatrixXd B2(3, 5);
			B2.col(0) = normal;
			B2.col(1) = tangent1;
			B2.col(2) = tangent2;
			B2.col(3) = -tangent1;
			B2.col(4) = -tangent2;
			MatrixXd mat = Jt * B1 * B2;
			_JtB.middleCols(k * 5, 5) = mat;
			if (j == 0)
			{
			    B2list.push_back(B2);
			    JtB = _JtB;
			}
			
		    }
		    H.leftCols(contactNodes[i + j].size() * 5) = hMinv * _JtB;
		    num_contacts += contactNodes[i + j].size();
		}
		H.rightCols(ndof - 6) = hMinv.rightCols(ndof - 6);
		aH.push_back(H);
		VectorXd hMinvC = hMinv * C;
		assert(hMinvC.rows() == ndof);
		ahMinvC.push_back(hMinvC);
	    }
	    // (lambda, xi_1, xi_2, q, t)
	    size_t num_var = num_lambda + ndof * numFrame + num_lambda + ndof * numFrame + 1;
	    size_t num_con = 4 * num_contacts + ndof * numFrame + num_lambda + ndof * numFrame;
	    MSKtask_t task = NULL;
	    mosekOK(MSK_maketask(env, num_con, num_var, &task));
	    mosekOK(MSK_linkfunctotaskstream(task, MSK_STREAM_LOG, NULL, printstr));
	    mosekOK(MSK_appendcons(task, num_con));
	    mosekOK(MSK_appendvars(task, num_var));

	    // set up linear term in the objective
	    for (size_t j = 0; j < num_var - 1; ++j)
		mosekOK(MSK_putcj(task, j, 0.0));
	    mosekOK(MSK_putcj(task, num_var - 1, 1.0));

	    // box constraints
	    size_t ii = 0;
	    for (size_t j = 0; j < numFrame; ++j)
	    {
		for (size_t k = 0; k < contactNodes[i + j].size() * 5; ++k)
		    mosekOK(MSK_putvarbound(task, ii++, MSK_BK_LO, 0.0, +MSK_INFINITY));
		for (size_t k = 0; k < ndof - 6; ++k)
		    mosekOK(MSK_putvarbound(task, ii++, MSK_BK_FR, -MSK_INFINITY, +MSK_INFINITY));
	    }
	    for (size_t j = ii; j < num_var; ++j)
		mosekOK(MSK_putvarbound(task, j, MSK_BK_FR, -MSK_INFINITY, +MSK_INFINITY));

	    // linear constraints
	    vector<MSKint32t> asubi, asubj;
	    vector<double> aval;
	    ii = 0;
	    size_t kk = 0;
	    for (size_t j = 0; j < numFrame; ++j)
	    {
		for (size_t jj = 0; jj < contactNodes[i + j].size(); ++jj)
		{
		    for (size_t k = 0; k < 4; ++k)
		    {
			asubi.push_back(kk);
			asubj.push_back(ii + 5 * jj + 0);
			aval.push_back(mu);
			asubi.push_back(kk);
			asubj.push_back(ii + 5 * jj + k + 1);
			aval.push_back(-1.0);
			mosekOK(MSK_putconbound(task, kk, MSK_BK_LO, 0.0, +MSK_INFINITY));
			++kk;
		    }
		}
		ii += 5 * contactNodes[i + j].size() + ndof - 6;
	    }
	    assert(kk == 4 * num_contacts);
	    size_t base = kk;
	    ii = 0;
	    for (size_t j = 0; j < numFrame; ++j)
	    {
		for (size_t k = 0; k < ndof; ++k)
		{
		    asubi.push_back(base + ii);
		    asubj.push_back(num_lambda + ii);
		    aval.push_back(-1.0);
		    asubi.push_back(base + ii);
		    asubj.push_back(num_lambda + ndof * numFrame + num_lambda + ii);
		    aval.push_back(1.0);
		    mosekOK(MSK_putconbound(task, base + ii, MSK_BK_FX, velocities[i + j][k], velocities[i + j][k]));
		    ++ii;
		}
	    }
	    assert(ii == ndof * numFrame);
	    base += ii;
	    for (size_t j = 0; j < num_lambda; ++j)
	    {
		asubi.push_back(base + j);
		asubj.push_back(j);
		aval.push_back(sqrt(reg / num_lambda));
		asubi.push_back(base + j);
		asubj.push_back(num_lambda + ndof * numFrame + j);
		aval.push_back(-1.0);
		mosekOK(MSK_putconbound(task, base + j, MSK_BK_FX, 0.0, 0.0));
	    }
	    base += num_lambda;
	    size_t base_k = 0;
	    size_t q_base = num_lambda + ndof * numFrame + num_lambda;
	    for (size_t j = 0; j < numFrame; ++j)
	    {
		MatrixXd &H = aH[j];
		VectorXd &hMinvC = ahMinvC[j];
		for (size_t jj = 0; jj < H.rows(); ++jj)
		    for (size_t kk = 0; kk < H.cols(); ++kk)
			if (abs(H(jj, kk)) > 1e-8)
			{
			    asubi.push_back(base + jj);
			    asubj.push_back(base_k + kk);
			    aval.push_back(-H(jj, kk));
			}
		for (size_t k = 0; k < ndof; ++k)
		{
		    asubi.push_back(base + k);
		    asubj.push_back(q_base + j * ndof + k);
		    aval.push_back(1.0);
		}
		if (j == 0)
		{
		    for (size_t k = 0; k < ndof; ++k)
			mosekOK(MSK_putconbound(task, base + k, MSK_BK_FX, vel[k] - hMinvC[k], vel[k] - hMinvC[k]));
		}
		else
		{
		    for (size_t k = 0; k < ndof; ++k)
		    {
			asubi.push_back(base + k);
			asubj.push_back(q_base + (j - 1) * ndof + k);
			aval.push_back(-1.0);
		    }
		    for (size_t k = 0; k < ndof; ++k)
			mosekOK(MSK_putconbound(task, base + k, MSK_BK_FX, -hMinvC[k], - hMinvC[k]));
		}
		base += ndof;
		base_k += H.cols();
	    }
	    mosekOK(MSK_putaijlist(task, aval.size(), asubi.data(), asubj.data(), aval.data()));

	    // conic constraints
	    vector<MSKint32t> csub(ndof * numFrame + num_lambda + 1);
	    csub[0] = num_lambda + ndof * numFrame + num_lambda + ndof * numFrame; // t
	    for (size_t i = 0; i < ndof * numFrame + num_lambda; ++i)
		csub[i + 1] = num_lambda + i;
	    mosekOK(MSK_appendcone(task, MSK_CT_QUAD, 0.0, csub.size(), csub.data()));

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

	    VectorXd lambda(contactNodes[i].size() * 5 + ndof - 6);
	    for (size_t i = 0; i < contactNodes[i].size() * 5 + ndof - 6; ++i)
		lambda[i] = x[i];
	    VectorXd Q = VectorXd::Zero(ndof);
	    Q.tail(ndof - 6) = lambda.tail(ndof - 6);
	    fout << Q.transpose() << endl;
	    for (size_t j = 0; j < contactNodes[i].size(); ++j)
	    {
		VectorXd lamb = lambda.segment(j * 5, 5);
		Vector3d f = B2list[j] * lamb;
		cfout << transformedContactPoints[i][j].transpose() << " " << f.transpose() << " ";
	    }
	    cfout << endl;
	    MatrixXd H = aH[0];
	    VectorXd hMinvC = ahMinvC[0];
	    vel_n = vel - hMinvC + H * lambda;
	    VectorXd err = vel_n - vel_hat;
	    eout << err.norm() << endl;
	    cout << "err1 " << err.norm() << endl;
	    cout << "err3 " << (H * lambda + vel - hMinvC - vel_hat).norm() << endl;
	    VectorXd acc = M.colPivHouseholderQr().solve(JtB * lambda.head(contactNodes[i].size() * 5) + Q - C);
	    cout << "err4 " << (velocities[i] + acc * h - velocities[i + 1]).norm() << endl;
	    cout << "err5 " << (velocities[i] + accelerations[i] * h - velocities[i + 1]).norm() << endl;
	}

#else

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
		//Isometry3d transform = bn->getWorldTransform();
		Isometry3d transform = kin_skeleton->getBodyNode(contactNodes[i][j])->getWorldTransform();
		//MatrixXd Jt = skeleton->getWorldJacobian(bn, transform.inverse() * point).transpose();
		MatrixXd Jt = skeleton->getWorldJacobian(bn).transpose();
		MatrixXd B1(6, 3);
		//B1.topRows(3) = dart::math::makeSkewSymmetric(point);
		//B1.topRows(3) = dart::math::makeSkewSymmetric(transform.inverse() * point);
		B1.topRows(3) = dart::math::makeSkewSymmetric(point - transform.translation());
		//B1.topRows(3) = Matrix3d::Zero();
		B1.bottomRows(3) = Matrix3d::Identity();
		//Vector3d normal = Vector3d::UnitY();	// y-axis points up
		//Vector3d tangent1 = Vector3d::UnitX();
		//Vector3d tangent2 = Vector3d::UnitZ();
		Vector3d normal = Vector3d::UnitZ();	// z-axis points up
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
		//for (size_t j = 0; j < n; ++j)
		//    A[i][j] = AtA(i, j);
		//for (size_t j = i; j < n; ++j)
		//    if (abs(AtA(i, j)) > 1e-8)
		//	sparseset(A, i, j, AtA(i, j));
		b[i] = -atA(0, i);
		s[i] = 1;
	    }
	    A.setcontent(n, n, AtA.data());
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

	    //Q.head(6) = VectorXd::Zero(6);
	    //acc = M.colPivHouseholderQr().solve(Q + emA * lambda - C);
	}
	else
	{
	    VectorXd Q = M * qddot + C;
	    //Q.head(6) = VectorXd::Zero(6);
	    eout << Q.head(6).norm() << endl; 
	    fout << Q.transpose() << endl;
	    cfout << endl;

	    //Q.head(6) = VectorXd::Zero(6);
	    //acc = M.colPivHouseholderQr().solve(Q - C);
	}
#endif

	/*
	pout << positions[i].transpose() << endl;
	vout << velocities[i].transpose() << endl;
	*/
#ifdef FULL_SOLVE
	pout << pos.transpose() << endl;
	vout << vel.transpose() << endl;
	aout << accelerations[i].transpose() << endl;

	vel = vel_n;
	skeleton->setVelocities(vel);
	skeleton->integratePositions(stepLength * frameTime);
	pos = skeleton->getPositions();
	//vel = vel_n;
#else
	pout << positions[i].transpose() << endl;
	vout << velocities[i].transpose() << endl;
	aout << accelerations[i].transpose() << endl;
#endif
    }
#if OPTIMIZOR == MOSEK
    MSK_deleteenv(&env);
#endif
    pout.close();
    fout.close();
    vout.close();
    aout.close();
    cfout.close();
    eout.close();

    return 0;
}
