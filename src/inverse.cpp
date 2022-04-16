#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <Eigen/Core>
#include "dart/dart.hpp"
#include "BVHData.h"
#include "ext/optimization.h"
#include "c_butterworth/c_butterworth_types.h"
#include "c_butterworth/c_butterworth.h"
#include "c_butterworth/rt_nonfinite.h"

using namespace std;
using namespace Eigen;
using namespace dart::dynamics;
using namespace dart::simulation;
using namespace alglib;

std::vector<double> readListFrom(const std::string &filename)
{
    std::ifstream input(filename);
    std::vector<double> list;
    double d;
    while (input >> d)
	list.push_back(d);
    input.close();
    return list;
}

std::vector<Eigen::VectorXd> readVectorXdListFrom(const std::string &filename)
{
    std::ifstream input(filename);
    std::vector<Eigen::VectorXd> result;
    std::string line;
    while (std::getline(input, line))
    {
	std::stringstream strStream(line);
	std::vector<double> list;
	double d;
	while (strStream >> d)
	    list.push_back(d);
	Eigen::VectorXd v(list.size());
	for (size_t i = 0; i < list.size(); ++i)
	    v[i] = list[i];
	result.push_back(v);
    }
    input.close();
    return result;
}

int main(int argc, char* argv[])
{
    string bvhFilename;
    string posFilename;
    string massFilename;
    string comFilename;
    string contactNodesFilename;
    double frameTime, cutoffFreq;
    if (argc == 8)
    {
	bvhFilename = argv[1];
	posFilename = argv[2];
	massFilename = argv[3];
	comFilename = argv[4];
	contactNodesFilename = argv[5];
	frameTime = stod(argv[6]);
	cutoffFreq = stod(argv[7]);
    }
    else
    {
	cout << "usage: " << argv[0] << " bvh_file pos_file mass_file com_file contact_file frame_time cutoff_freq" << endl;
	exit(0);
    }
    BVHData bvh;
    bvh.loadBVH(bvhFilename, "", "", 1);
    WorldPtr world = World::create();
    SkeletonPtr &skeleton = bvh.skeleton;

    vector<double> massList = readListFrom(massFilename);
    vector<VectorXd> comList = readVectorXdListFrom(comFilename);
    std::vector<dart::dynamics::BodyNode*> bns = bvh.skeleton->getBodyNodes();
    for (size_t i = 0; i < bns.size(); ++i)
    {
	bns[i]->setMass(massList[i]);
	bns[i]->setLocalCOM(comList[i]);
	cout << bns[i]->getName() << endl;
	cout << bns[i]->getMass() << endl;
	cout << bns[i]->getLocalCOM() << endl;
    }

    world->addSkeleton(skeleton);
    vector<VectorXd> positions;
    vector<VectorXd> velocities;
    vector<VectorXd> accelerations;
    //skeleton->setGravity(Vector3d(0, -9.8, 0));
    skeleton->setGravity(Vector3d(0, 0, -9.8));
    positions = readVectorXdListFrom(posFilename);
    for (size_t i = 0; i < positions.size() - 1; ++i)
	velocities.push_back(bvh.skeleton->getPositionDifferences(positions[i + 1], positions[i]) / frameTime);

    ofstream frameOut("bvh-frame.txt");
    for (VectorXd &v: bvh.frameToEulerAngle(positions))
	frameOut << v.transpose() << endl;
    frameOut.close();

    coder::array<double, 1> x, y;
    vector<VectorXd> tmp(velocities);
    x.set_size(velocities.size());
    size_t ndof = skeleton->getDofs().size();
    for (size_t i = 0; i < ndof; ++i)
    {
	for (int j = 0; j < x.size(0); ++j)
	    x[j] = velocities[j][i];
	c_butterworth(x, cutoffFreq, 100, y);
	for (int j = 0; j < x.size(0); ++j)
	    velocities[j][i] = y[j];
    }
    c_butterworth_terminate();

    for (size_t i = 0; i < velocities.size() - 1; ++i)
	accelerations.push_back(bvh.skeleton->getVelocityDifferences(velocities[i + 1], velocities[i]) / frameTime);

    vector<vector<string>> contactNodes;
    ifstream input(contactNodesFilename);
    string line;
    while (getline(input, line))
    {
	stringstream strStream(line);
	vector<string> list;
	string s;
	while (strStream >> s)
	    list.push_back(s);
	contactNodes.push_back(list);
    }
    input.close();

    ofstream pout("positions.txt");
    ofstream fout("forces.txt");
    ofstream vout("velocities.txt");
    ofstream aout("accelerations.txt");

    for (size_t i = 0; i < accelerations.size(); ++i)
    {
	cout << "frame " << i << endl;
	skeleton->setPositions(positions[i]);
	skeleton->setVelocities(velocities[i]);
	skeleton->setAccelerations(accelerations[i]);
	MatrixXd M = skeleton->getMassMatrix();
	MatrixXd C = skeleton->getCoriolisAndGravityForces();

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

	if (contactNodes[i].size() == 1)
	{
	    MatrixXd J = skeleton->getWorldJacobian(skeleton->getBodyNode(contactNodes[i][0]), skeleton->getBodyNode(contactNodes[i][0])->getLocalCOM());
	    MatrixXd Jt = J.transpose();
	    VectorXd w = Jt.topRows(6).inverse() *  (M * qddot + C).topRows(6);
	    VectorXd Q = M * qddot + C - Jt * w;
	    fout << Q.transpose() << endl;
	    cout << w.transpose() << endl;
	}
	else if (contactNodes[i].size() > 1)
	{
	    vector<MatrixXd> Jt;
	    vector<VectorXd> w;
	    for (string &name: contactNodes[i])
	    {
		Jt.push_back(skeleton->getWorldJacobian(skeleton->getBodyNode(name), skeleton->getBodyNode(name)->getLocalCOM()).transpose());
		w.push_back(VectorXd(6));
	    }
	    size_t n = contactNodes[i].size();
	    real_2d_array A;
	    real_1d_array b;
	    real_1d_array s;
	    A.setlength(6 * n, 6 * n);
	    b.setlength(6 * n);
	    s.setlength(6 * n);
	    for (size_t i = 0; i < 6 * n; ++i)
		for (size_t j = 0; j < 6 * n; ++j)
		    A[i][j] = 0;
	    for (size_t i = 0; i < 6 * n; ++i)
	    {
		A[i][i] = 1;
		b[i] = 0;
		s[i] = 1;
	    }
	    VectorXd d = M * qddot + C;
	    //cout << "d:" << endl << d << endl;
	    real_2d_array c;
	    integer_1d_array ct;
	    c.setlength(6, 6 * n + 1);
	    ct.setlength(6);
	    for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < 6; ++j)
		    for (size_t k = 0; k < 6; ++k)
			c[j][6 * i + k] = Jt[i](j, k);
	    for (size_t i = 0; i < 6; ++i)
	    {
		c[i][6 * n] = d[i];
		ct[i] = 0;
	    }
	    real_1d_array x;
	    minqpstate state;
	    minqpreport rep;
	    minqpcreate(6 * n, state);
	    minqpsetquadraticterm(state, A);
	    minqpsetlinearterm(state, b);
	    minqpsetlc(state, c, ct);
	    minqpsetscale(state, s);
	    minqpsetalgobleic(state, 0.0, 0.0, 0.0, 0);
	    minqpoptimize(state);
	    minqpresults(state, x, rep);
	    //cout << "A:" << endl << A.tostring(1).c_str() << endl;
	    //cout << "b:" << endl << b.tostring(1).c_str() << endl;
	    //cout << "s:" << endl << s.tostring(1).c_str() << endl;
	    //cout << "c:" << endl << c.tostring(1).c_str() << endl;
	    //cout << "ct:" << endl << ct.tostring().c_str() << endl;
	    //cout << "x:" << endl << x.tostring(1).c_str() << endl;
	    for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < 6; ++j)
		    w[i][j] = x[6 * i + j];
	    VectorXd Q = M * qddot + C;
	    for (size_t i = 0; i < n; ++i)
		Q -= Jt[i] * w[i];
	    fout << Q.transpose() << endl;
	    for (size_t i = 0; i < n; ++i)
		cout << w[i].transpose() << endl;
	}
	else
	{
	    VectorXd Q = M * qddot + C;
	    Q.head(6) = VectorXd::Zero(6);
	    fout << Q.transpose() << endl;
	}

	pout << positions[i].transpose() << endl;
	vout << velocities[i].transpose() << endl;
	aout << accelerations[i].transpose() << endl;
    }
    pout.close();
    fout.close();
    vout.close();
    aout.close();

    return 0;
}
