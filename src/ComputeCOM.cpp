#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <Eigen/Core>
#include <getopt.h>
#include <dart/dart.hpp>
#include <Python.h>
#include "SimCharacter.h"
#include "PyUtil.h"
#include "IOUtil.h"

using namespace std;
using namespace Eigen;
using namespace dart::dynamics;
using namespace dart::simulation;

void printUsage(char * prgname)
{
    cout << "usage: " << prgname << " [options] pose_file output_file" << endl;
    cout << endl;
    cout << "options:" << endl;
    cout << "\t-j, --char_file=string" << endl;
    cout << "\t-z, --npz" << endl;
}

int main(int argc, char* argv[])
{
    string jsonFilename = "data/character.json";
    bool is_npz = false;
    string poseFilename;
    string outName;

    while (1)
    {
	int c;
	static struct option long_options[] =
	{
	    { "char_file", required_argument, NULL, 'j' },
	    { "npz", no_argument, NULL, 0 },
	    { 0, 0, 0, 0 }
	};
	int option_index = 0;
	c = getopt_long(argc, argv, "j:o:z", long_options, &option_index);
	if (c == -1)
	    break;

	switch (c)
	{
	    case 'j':
		jsonFilename = optarg;
		break;
	    case 'z':
		is_npz = true;
		break;
	    default:
		printUsage(argv[0]);
		exit(0);
	}
    }

    if (optind + 2 == argc)
    {
	poseFilename = argv[optind];
	outName = argv[optind + 1];
    }
    else
    {
	printUsage(argv[0]);
	exit(0);
    }

    SimCharacter character(jsonFilename);
    WorldPtr world = World::create();
    SkeletonPtr &skeleton = character.skeleton;
    world->addSkeleton(skeleton);
    //skeleton->setGravity(Vector3d(0, -9.8, 0));   // y-axis points up
    skeleton->setGravity(Vector3d(0, 0, -9.8));	    // z-axis points up

    vector<VectorXd> positions;
    vector<VectorXd> velocities;
    vector<VectorXd> accelerations;
    Eigen::MatrixXd position_mat;

    if (is_npz)
    {
	Py_Initialize();
	PyRun_SimpleString("import sys\n"
		"sys.path.append(\".\")\n"
		"print('python version:', sys.version)\n");

	PyObject *pDict = load_npz(poseFilename);
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

	 position_mat = Eigen::MatrixXd(poses.rows(), poses.cols() + 3);
	position_mat.leftCols(3) = poses.leftCols(3);
	position_mat.middleCols(3, 3) = trans;
	position_mat.rightCols(poses.cols() - 3) = poses.rightCols(poses.cols() - 3);

    }
    else
    {
	position_mat = readMatrixXFrom(poseFilename);
    }

    for (size_t i = 0; i < position_mat.rows(); ++i)
	positions.push_back(position_mat.row(i));
    size_t f_start = 0, f_end = positions.size() - 2;
    ofstream com_out(outName);
    
    for (size_t i = f_start; i < f_end; ++i)
    {
	skeleton->setPositions(positions[i]);
	com_out << skeleton->getCOM().transpose() << endl;
    }

    return 0;
}
