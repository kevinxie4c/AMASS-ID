#include <iostream>
#include <fstream>
#include <vector>
#include <nlohmann/json.hpp>
#include <string>
#include <dart/dart.hpp>
#include <dart/gui/osg/osg.hpp>
#include <dart/external/imgui/imgui.h>
#include <osgUtil/SmoothingVisitor>
#include <Python.h>
#include "SimCharacter.h"
#include "PyUtil.h"
#include "IOUtil.h"

int frame = 0;
Eigen::MatrixXd positions;
SimCharacter character("data/character.json");
osg::ref_ptr<osg::Geometry> mesh;
std::vector<Eigen::MatrixXd> vertices;
std::vector<Eigen::VectorXd> pointIndices;
std::vector<Eigen::VectorXd> contacts;
osg::ref_ptr<osg::Group> group;

class TestWidget : public dart::gui::osg::ImGuiWidget
{
public:

    TestWidget(dart::gui::osg::ImGuiViewer* viewer, dart::simulation::WorldPtr world) : mViewer(viewer), mWorld(world)
    {
    }

    void render() override
    {
	ImGui::SetNextWindowPos(ImVec2(10, 20));
	ImGui::SetNextWindowSize(ImVec2(2000, 200));
	ImGui::Begin("Control");
	//ImGui::SliderInt("frame", &frame, 0, positions.rows() - 1);
	ImGui::SliderInt("frame", &frame, 0, 499);
	character.skeleton->setPositions(positions.row(frame));
	if (frame < vertices.size())
	{
	    osg::Vec3Array* vertexArray = static_cast<osg::Vec3Array*>(mesh->getVertexArray());
	    osg::Vec4Array* colorArray = static_cast<osg::Vec4Array*>(mesh->getColorArray());
	    for (int i = 0; i < vertices[0].rows(); ++i)
	    {
		Eigen::VectorXd v = vertices[frame].row(i);
		(*vertexArray)[i] = osg::Vec3(v[0], v[1], v[2]);
		(*colorArray)[i] = osg::Vec4(1.0, 1.0, 1.0, 1.0);
	    }
	    for (int i = 0; i < pointIndices[frame].rows(); ++i)
	    {
		(*colorArray)[pointIndices[frame][i]] = osg::Vec4(1.0, 0.0, 0.0, 1.0);
	    }
	    group->removeChildren(0, group->getNumChildren());
	    for (size_t i = 0; i < contacts[frame].rows(); i += 6)
	    {
		Eigen::VectorXd v = contacts[frame];
		osg::ref_ptr<osg::ShapeDrawable> shape = new osg::ShapeDrawable;
		shape->setShape(new osg::Sphere(osg::Vec3(v[i + 0], v[i + 1], v[i + 2]), 0.02));
		shape->setColor(osg::Vec4(0.0, 1.0, 0.0, 1.0));
		group->addChild(shape);
	    }
	    mesh->dirtyDisplayList();
	    mesh->dirtyBound();
	}
	ImGui::End();
    }

private:

    osg::ref_ptr<dart::gui::osg::ImGuiViewer> mViewer;
    dart::simulation::WorldPtr mWorld;

};


int main(int argc, char *argv[])
{
    char *fname = argv[1];
    //Py_SetProgramName(argv[0]);
    Py_Initialize();
    //PySys_SetArgv(argc, argv);
    PyRun_SimpleString("import sys\n"
		       "sys.path.append(\".\")\n"
                       "print('python version:', sys.version)\n");

    PyObject *pDict = load_npz(fname);
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

    std::cout << "loading mesh..." << std::endl;
    const char *module = "amass_mesh";
    const char *func = "AMASSmesh";
    PyObject *pName, *pModule, *pFunc, *pValue, *pArgs;
    pName = PyUnicode_FromString(module);
    pModule = PyImport_Import(pName);
    Py_DECREF(pName);
    if (pModule != NULL)
    {
	pFunc = PyObject_GetAttrString(pModule, func);
	if (pFunc && PyCallable_Check(pFunc))
	{
	    pValue = PyUnicode_FromString(fname);
	    pArgs = PyTuple_New(1);
	    PyTuple_SetItem(pArgs, 0, pValue);
	    pValue = PyObject_CallObject(pFunc, pArgs);
	    Py_DECREF(pArgs);
	    Py_DECREF(pFunc);
	    Py_DECREF(pModule);
	}
	else
	{
	    if (PyErr_Occurred())
                PyErr_Print();
            fprintf(stderr, "failed to call function \"%s\"\n", func);
	    return 0;
	}
	Py_XDECREF(pFunc);
        Py_DECREF(pModule);
    }
    else
    {
	PyErr_Print();
        fprintf(stderr, "failed to load \"%s\"\n", module);
	return 0;
    }
    PyObject *pFaces, *pVertices;
    pFaces = PyObject_GetAttrString(pValue, "faces");
    pVertices = PyObject_GetAttrString(pValue, "vertices");
    Eigen::MatrixXd faces = npa2mat(pFaces);
    vertices = npa2mat_vector(pVertices);
    //std::cout << faces.row(0) << std::endl;
    //std::cout << vertices[0].row(0) << std::endl;
    Py_DECREF(pFaces);
    Py_DECREF(pVertices);
    Py_DECREF(pValue);
    std::cout << "mesh loaded" << std::endl;
    /*
    pDict = load_npz("data/0005_Walking001_mesh.npz");
    npa = PyMapping_GetItemString(pDict, "faces");
    Eigen::MatrixXd faces = npa2mat(npa);
    npa = PyMapping_GetItemString(pDict, "vertices");
    std::vector<Eigen::MatrixXd> vertices = npa2mat_vector(npa);
    Py_DECREF(pDict);
    std::cout << faces.row(0) << std::endl;
    std::cout << vertices[0].row(0) << std::endl;
    */

    std::vector<size_t> indices{0, 3, 12, 21, 30, 6, 15, 24, 33, 9, 18, 27, 36, 45, 39, 48, 54, 60, 42, 51, 57, 63};
    Eigen::MatrixXd poses(poses_orig.rows(), indices.size() * 3);
    for (size_t i = 0; i < indices.size(); ++i)
	poses.middleCols(i * 3, 3) = poses_orig.middleCols(indices[i], 3);

    positions = Eigen::MatrixXd(poses.rows(), poses.cols() + 3);
    positions.leftCols(3) = poses.leftCols(3);
    positions.middleCols(3, 3) = trans;
    positions.rightCols(poses.cols() - 3) = poses.rightCols(poses.cols() - 3);

    pointIndices = readVectorXdListFrom("contact-point-indices.txt");

    for (size_t i = 0; i < character.skeleton->getNumDofs(); ++i)
	std::cout << i << "\t" << character.skeleton->getDof(i)->getName() << std::endl;
    dart::simulation::WorldPtr world = dart::simulation::World::create();
    world->addSkeleton(character.skeleton);
    character.skeleton->setPositions(positions.row(frame));
    osg::ref_ptr<dart::gui::osg::ImGuiViewer> viewer = new dart::gui::osg::ImGuiViewer(osg::Vec4(0.1, 0.1, 0.1, 1.0));
    osg::ref_ptr<dart::gui::osg::WorldNode> worldNode = new dart::gui::osg::RealTimeWorldNode(world);

    osg::ref_ptr<osg::Vec3Array> vertexArray = new osg::Vec3Array();
    osg::ref_ptr<osg::Vec4Array> colorArray = new osg::Vec4Array();
    osg::ref_ptr<osg::DrawElementsUInt> indexArray = new osg::DrawElementsUInt(GL_TRIANGLES);
    for (int i = 0; i < vertices[0].rows(); ++i)
    {
	Eigen::VectorXd v = vertices[0].row(i);
	vertexArray->push_back(osg::Vec3(v[0], v[1], v[2]));
	colorArray->push_back(osg::Vec4(1.0, 1.0, 1.0, 1.0));
    }
    for (int i = 0; i < pointIndices[0].rows(); ++i)
    {
	(*colorArray)[pointIndices[0][i]] = osg::Vec4(1.0, 0.0, 0.0, 1.0);
    }
    for (int i = 0; i < faces.rows(); ++i)
    {
	Eigen::VectorXd v = faces.row(i);
	indexArray->push_back(v[0]);
	indexArray->push_back(v[1]);
	indexArray->push_back(v[2]);
    }
    mesh = new osg::Geometry();
    mesh->setVertexArray(vertexArray);
    mesh->setColorArray(colorArray);
    mesh->setColorBinding(osg::Geometry::BIND_PER_VERTEX);
    mesh->addPrimitiveSet(indexArray);
    osgUtil::SmoothingVisitor::smooth(*mesh);
    mesh->setDataVariance(osg::Object::DYNAMIC);
    osg::ref_ptr<osg::Geode> geoNode = new osg::Geode();
    geoNode->addDrawable(mesh);
    worldNode->addChild(geoNode);

    contacts = readVectorXdListFrom("contact_forces.txt");
    group = new osg::Group();
    for (size_t i = 0; i < contacts[0].rows(); i += 6)
    {
	Eigen::VectorXd v = contacts[0];
	osg::ref_ptr<osg::ShapeDrawable> shape = new osg::ShapeDrawable;
	shape->setShape(new osg::Sphere(osg::Vec3(v[i + 0], v[i + 1], v[i + 2]), 0.02));
	shape->setColor(osg::Vec4(0.0, 1.0, 0.0, 1.0));
	group->addChild(shape);
    }
    worldNode->addChild(group);

    viewer->addWorldNode(worldNode);
    viewer->getImGuiHandler()->addWidget(std::make_shared<TestWidget>(viewer, world));

    osg::Vec3d eye(0, 0, 5);
    osg::Vec3d center(0, 0, 0);
    osg::Vec3d up(0, 1, 0);
    viewer->getCameraManipulator()->setHomePosition(eye, center, up);
    viewer->run();

    if (Py_FinalizeEx() < 0) {
        return 120;
    }

    return 0;
}
