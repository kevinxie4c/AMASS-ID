#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <Eigen/Core>
#include <dart/dart.hpp>
#include <Python.h>
#include "SimCharacter.h"
#include "PyUtil.h"
#include "IOUtil.h"
#include "ext/optimization.h"
#include "c_butterworth/c_butterworth_types.h"
#include "c_butterworth/c_butterworth.h"
#include "c_butterworth/rt_nonfinite.h"

using namespace std;
using namespace Eigen;
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace alglib;

int main(int argc, char* argv[])
{
    string jsonFilename;
    string posFilename;
    string contactNodesFilename;
    double frameTime, cutoffFreq, groundOffset;
    if (argc == 7)
    {
	jsonFilename = argv[1];
	posFilename = argv[2];
	contactNodesFilename = argv[3];
	frameTime = stod(argv[4]);
	cutoffFreq = stod(argv[5]);
	groundOffset = stod(argv[6]);
    }
    else
    {
	cout << "usage: " << argv[0] << " json_file pos_file contact_file frame_time cutoff_freq ground_offset" << endl;
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

    PyObject *pDict = load_npz(posFilename);
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
	if (strStream >> s >> n)
	{
	    double x, y, z;
	    for (size_t i = 0; i < n; ++i)
	    {
		strStream >> x >> y >> z;
		nodeList.push_back(s);
		pointList.push_back(Vector3d(x, y, z));
	    }
	    contactNodes.push_back(nodeList);
	    contactPoints.push_back(pointList);
	}
    }
    input.close();

    ofstream pout("positions.txt");
    ofstream fout("forces.txt");
    ofstream vout("velocities.txt");
    ofstream aout("accelerations.txt");
    ofstream cfout("contact_forces.txt");
    ofstream eout("errors.txt");

    double mu = 1;
    double reg = 10;
    
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
	Eigen::Vector3d dL = Eigen::Vector3d::Zero(); // rate of angular momentum at center of mass
	Eigen::Vector3d com = skeleton->getCOM();
	for (const BodyNode *bn: skeleton->getBodyNodes())
	{
	    dL += (bn->getCOM() - com).cross(bn->getMass() * bn->getCOMLinearAcceleration());

	    Eigen::Matrix3d R = bn->getTransform().rotation();
	    Eigen::Matrix3d I = bn->getInertia().getMoment(); // moment of inertia in local frame
	    Eigen::Vector3d omiga = bn->getCOMSpatialVelocity().head(3); // angular velocity in local frame
	    Eigen::Vector3d dOmiga = bn->getCOMSpatialAcceleration().head(3); // angular acceleration in local frame

	    dL += R * (I * dOmiga + omiga.cross(I * omiga));
	}

	VectorXd qddot = accelerations[i];

	
	if (contactNodes[i].size() > 0)
	{
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
		B1.topRows(3) = dart::math::makeSkewSymmetric(point - transform.translation());
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
	    real_2d_array A;
	    //sparsematrix A;
	    real_1d_array b;
	    real_1d_array s;
	    A.setlength(n, n);
	    //sparsecreate(n, n, 0, A);
	    b.setlength(n);
	    s.setlength(n);
	    AtA += MatrixXd::Identity(n, n) * reg / n;
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
	    minqpoptimize(state);
	    minqpresults(state, x, rep);

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
    pout.close();
    fout.close();
    vout.close();
    aout.close();
    cfout.close();
    eout.close();

    return 0;
}
