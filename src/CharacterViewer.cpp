#include <iostream>
#include <fstream>
#include <vector>
#include <nlohmann/json.hpp>
#include <string>
#include <getopt.h>
#include <dart/dart.hpp>
#include <dart/gui/osg/osg.hpp>
#include <dart/external/imgui/imgui.h>
#include <osgUtil/SmoothingVisitor>
#include <Python.h>
#include <osg/BlendFunc>
#include <osg/CullFace>
#include <osg/FrontFace>
#include "SimCharacter.h"
#include "PyUtil.h"
#include "IOUtil.h"

double forceScale = 0.003;
double epsilon = 0.01;
float meshAlpha = 0.3;
bool inverted = false;
bool transparent = true;
bool playing = false;
bool showWorldFrame = false;

Eigen::MatrixXd positions;
dart::dynamics::SkeletonPtr skeleton = nullptr;
osg::ref_ptr<osg::Geometry> mesh;
std::vector<Eigen::MatrixXd> vertices;
std::vector<Eigen::VectorXd> pointIndices;
std::vector<Eigen::VectorXd> contacts;
osg::ref_ptr<osg::Group> group;
osg::StateSet *stateSet;
size_t startFrame = 0, endFrame = 999;
int frame = startFrame;

osg::ref_ptr<osg::Node> makeArrow(const osg::Vec3 &start, const osg::Vec3 &end, const osg::Vec4 &color = osg::Vec4(0.0, 1.0, 0.0, 1.0))
{
    osg::Vec3 dir = end - start;
    osg::ref_ptr<osg::MatrixTransform> transform = new osg::MatrixTransform();
    osg::ref_ptr<osg::Geode> node = new osg::Geode();
    osg::ref_ptr<osg::ShapeDrawable> cylinder = new osg::ShapeDrawable();
    osg::ref_ptr<osg::Shape> cylinderShape = new osg::Cylinder(osg::Vec3(0, 0, dir.length() / 2), 0.002, dir.length());
    cylinder->setShape(cylinderShape.get());
    cylinder->setColor(color);
    osg::ref_ptr<osg::ShapeDrawable> cone = new osg::ShapeDrawable();
    osg::ref_ptr<osg::Shape> coneShape = new osg::Cone(osg::Vec3(0, 0, dir.length()), 0.004, 0.008);
    cone->setShape(coneShape.get());
    cone->setColor(color);
    node->addDrawable(cylinder.get());
    node->addDrawable(cone.get());
    dir.normalize();
    //node->setMatrix(osg::Matrix::translate(start) * osg::Matrix::rotate(osg::Vec3(0, 0, 1), dir));
    //std::cout << "dir " << dir[0] << " " << dir[1] << " " << dir[2] << std::endl;
    //std::cout << "start " << start[0] << " " << start[1] << " " << start[2] << std::endl;
    transform->addChild(node.get());
    transform->setMatrix(osg::Matrix::rotate(osg::Vec3(0, 0, 1), dir) * osg::Matrix::translate(start));
    return transform;
}

class CustomCallBack: public osg::NodeCallback
{
public:
    virtual void operator()(osg::Node* node, osg::NodeVisitor* nv )
    {
	osg::Group *group = (osg::Group*)node;
	group->removeChildren(0, group->getNumChildren());
	if (showWorldFrame)
	{
	    group->addChild(makeArrow(osg::Vec3(0, 0, 0), osg::Vec3(1, 0, 0), osg::Vec4(1, 0, 0, 1)));
	    group->addChild(makeArrow(osg::Vec3(0, 0, 0), osg::Vec3(0, 1, 0), osg::Vec4(0, 1, 0, 1)));
	    group->addChild(makeArrow(osg::Vec3(0, 0, 0), osg::Vec3(0, 0, 1), osg::Vec4(0, 0, 1, 1)));
	}
	const Eigen::VectorXd &v = contacts[frame];
	for (size_t i = 0; i < v.rows(); i += 6)
	{
	    //osg::ref_ptr<osg::ShapeDrawable> shape = new osg::ShapeDrawable();
	    //shape->setShape(new osg::Sphere(osg::Vec3(v[i + 0], v[i + 1], v[i + 2]), 0.02));
	    //shape->setColor(osg::Vec4(0.0, 1.0, 0.0, 1.0));
	    //group->addChild(shape);
	    //shape = new osg::ShapeDrawable();
	    //shape->setShape(new osg::Sphere(osg::Vec3(v[i + 0] + forceScale * v[i + 3], v[i + 1] + forceScale * v[i + 4], v[i + 2] + forceScale * v[i + 5]), 0.02));
	    //shape->setColor(osg::Vec4(0.0, 1.0, 0.0, 1.0));
	    //group->addChild(shape);
	    osg::Vec3 point = osg::Vec3(v[i + 0], v[i + 1], v[i + 2]);
	    osg::Vec3 force = osg::Vec3(v[i + 3], v[i + 4], v[i + 5]);
	    //std::cout << "point " << point[0] << " " << point[1] << " " << point[2] << std::endl;
	    double sign = inverted ? -1 : 1;
	    if (force.length() > epsilon)
		group->addChild(makeArrow(point, point + force * sign * forceScale).get());
	}
    }
};

class TestWidget: public dart::gui::osg::ImGuiWidget
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
	ImGui::SliderInt("frame", &frame, startFrame, endFrame - 1);
	ImGui::Checkbox("Inverted contact force", &inverted);
	ImGui::Checkbox("Transparent mesh", &transparent);
	ImGui::Checkbox("Show world frame", &showWorldFrame);
	ImGui::SliderFloat("Alpha", &meshAlpha, 0.0f, 1.0f);
	ImGui::End();
	if (showWorldFrame)
	{
	}
	if (transparent)
	{
	    stateSet->setMode(GL_DEPTH_TEST, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	    stateSet->setMode(GL_BLEND, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
	}
	else
	{
	    stateSet->setMode(GL_DEPTH_TEST, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
	    stateSet->setMode(GL_BLEND, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	}
	if (skeleton != nullptr)
	    skeleton->setPositions(positions.row(frame));
	if (frame < vertices.size())
	{
	    osg::Vec3Array* vertexArray = static_cast<osg::Vec3Array*>(mesh->getVertexArray());
	    osg::Vec4Array* colorArray = static_cast<osg::Vec4Array*>(mesh->getColorArray());
	    for (int i = 0; i < vertices[0].rows(); ++i)
	    {
		const Eigen::VectorXd &v = vertices[frame].row(i);
		(*vertexArray)[i] = osg::Vec3(v[0], v[1], v[2]);
		(*colorArray)[i] = osg::Vec4(1.0, 1.0, 1.0, meshAlpha);
	    }
	    for (int i = 0; i < pointIndices[frame].rows(); ++i)
	    {
		(*colorArray)[pointIndices[frame][i]] = osg::Vec4(1.0, 0.0, 0.0, meshAlpha);
	    }
	    mesh->dirtyDisplayList();
	    mesh->dirtyBound();
	    /*
	    group->removeChildren(0, group->getNumChildren());
	    for (size_t i = 0; i < contacts[frame].rows(); i += 6)
	    {
		const Eigen::VectorXd &v = contacts[frame];
		//osg::ref_ptr<osg::ShapeDrawable> shape = new osg::ShapeDrawable();
		//shape->setShape(new osg::Sphere(osg::Vec3(v[i + 0], v[i + 1], v[i + 2]), 0.02));
		//shape->setColor(osg::Vec4(0.0, 1.0, 0.0, 1.0));
		//group->addChild(shape);
		//shape = new osg::ShapeDrawable();
		//shape->setShape(new osg::Sphere(osg::Vec3(v[i + 0] + forceScale * v[i + 3], v[i + 1] + forceScale * v[i + 4], v[i + 2] + forceScale * v[i + 5]), 0.02));
		//shape->setColor(osg::Vec4(0.0, 1.0, 0.0, 1.0));
		//group->addChild(shape);
		osg::Vec3 point = osg::Vec3(v[i + 0], v[i + 1], v[i + 2]);
		osg::Vec3 force = osg::Vec3(v[i + 3], v[i + 4], v[i + 5]);
		//std::cout << "point " << point[0] << " " << point[1] << " " << point[2] << std::endl;
		if (force.length() > epsilon)
		    group->addChild(makeArrow(point, point + force * forceScale).get());
	    }
	    */
	}
    }

private:

    osg::ref_ptr<dart::gui::osg::ImGuiViewer> mViewer;
    dart::simulation::WorldPtr mWorld;

};

class CustomEventHandler: public osgGA::GUIEventHandler
{
public:

    CustomEventHandler(osg::ref_ptr<dart::gui::osg::ImGuiViewer> viewer) : mViewer(viewer) {}

    virtual bool handle(const osgGA::GUIEventAdapter& ea, osgGA::GUIActionAdapter&) override
    {
	switch (ea.getEventType())
	{
	    case osgGA::GUIEventAdapter::KEYDOWN:
		switch (ea.getKey())
		{
		    case 'r':
			if (mViewer->isRecording())
			    mViewer->pauseRecording();
			else
			    mViewer->record("output");
			return true;
		    case ' ':
			playing = !playing;
			return true;
		}
		break;
	    case osgGA::GUIEventAdapter::FRAME:
		if (playing && frame < endFrame)
		    ++frame;
		return true;
	}
	return false;
    }

    osg::ref_ptr<dart::gui::osg::ImGuiViewer> mViewer;

};

void printUsage(char * prgname)
{
    std::cout << "usage: " << prgname << " AMASS_pose_file contact_point_indices_file contact_forces_file" << std::endl;
    std::cout << std::endl;
    std::cout << "options:" << std::endl;
    std::cout << "\t-j, --char_file=string" << std::endl;
    std::cout << "\t-p, --pose_file=string" << std::endl;
    std::cout << "\t-A, --start_frame=int" << std::endl;
    std::cout << "\t-E, --end_frame=int" << std::endl;
}

int main(int argc, char *argv[])
{
    std::string jsonFilename = "data/character.json";
    std::string amassFilename;
    std::string poseFilename;
    std::string cpiFilename;
    std::string cfFilename;

    while (1)
    {
	int c;
	static struct option long_options[] =
	{
	    { "char_file", required_argument, NULL, 'j' },
	    { "pose_file", required_argument, NULL, 'p' },
	    { "start_frame", required_argument, NULL, 'A' },
	    { "end_frame", required_argument, NULL, 'E' },
	    { 0, 0, 0, 0 }
	};
	int option_index = 0;

	c = getopt_long(argc, argv, "j:p:A:E:", long_options, &option_index);
	if (c == -1)
        break;

	switch (c)
	{
	    case 'j':
		jsonFilename = optarg;
		break;
	    case 'p':
		poseFilename = optarg;
		break;
	    case 'A':
		startFrame = std::stoi(optarg);
		break;
	    case 'E':
		endFrame = std::stoi(optarg);
		break;
	    default:
		printUsage(argv[0]);
		exit(0);
	}
    }

    if (optind + 3 == argc)
    {
	amassFilename = argv[optind];
	cpiFilename = argv[optind + 1];
	cfFilename = argv[optind + 2];
    }
    else
    {
	printUsage(argv[0]);
	exit(0);
    }

    frame = startFrame;

    dart::simulation::WorldPtr world = dart::simulation::World::create();
    world->setGravity(Eigen::Vector3d::Zero());
    osg::ref_ptr<dart::gui::osg::ImGuiViewer> viewer = new dart::gui::osg::ImGuiViewer(osg::Vec4(0.1, 0.1, 0.1, 1.0));
    osg::ref_ptr<dart::gui::osg::WorldNode> worldNode = new dart::gui::osg::RealTimeWorldNode(world);

    pointIndices = readVectorXdListFrom(cpiFilename);

    //Py_SetProgramName(argv[0]);
    Py_Initialize();
    //PySys_SetArgv(argc, argv);
    PyRun_SimpleString("import sys\n"
		       "sys.path.append(\".\")\n"
		       "print('python version:', sys.version)\n");

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
	    pValue = PyUnicode_FromString(amassFilename.c_str());
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

    /*
    PyObject *pDict = load_npz(amassFilename);
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

    positions = Eigen::MatrixXd(poses.rows(), poses.cols() + 3);
    positions.leftCols(3) = poses.leftCols(3);
    positions.middleCols(3, 3) = trans;
    positions.rightCols(poses.cols() - 3) = poses.rightCols(poses.cols() - 3);
    */

    osg::ref_ptr<osg::Vec3Array> vertexArray = new osg::Vec3Array();
    osg::ref_ptr<osg::Vec4Array> colorArray = new osg::Vec4Array();
    osg::ref_ptr<osg::DrawElementsUInt> indexArray = new osg::DrawElementsUInt(GL_TRIANGLES);
    for (int i = 0; i < vertices[0].rows(); ++i)
    {
	const Eigen::VectorXd &v = vertices[0].row(i);
	vertexArray->push_back(osg::Vec3(v[0], v[1], v[2]));
	colorArray->push_back(osg::Vec4(1.0, 1.0, 1.0, meshAlpha));
    }
    for (int i = 0; i < pointIndices[0].rows(); ++i)
    {
	(*colorArray)[pointIndices[0][i]] = osg::Vec4(1.0, 0.0, 0.0, meshAlpha);
    }
    for (int i = 0; i < faces.rows(); ++i)
    {
	const Eigen::VectorXd &v = faces.row(i);
	indexArray->push_back(v[0]);
	indexArray->push_back(v[1]);
	indexArray->push_back(v[2]);
    }
    mesh = new osg::Geometry();
    mesh->setVertexArray(vertexArray.get());
    mesh->setColorArray(colorArray.get());
    mesh->setColorBinding(osg::Geometry::BIND_PER_VERTEX);
    mesh->addPrimitiveSet(indexArray.get());
    osgUtil::SmoothingVisitor::smooth(*mesh);
    mesh->setDataVariance(osg::Object::DYNAMIC);
    osg::ref_ptr<osg::Geode> geoNode = new osg::Geode();
    geoNode->addDrawable(mesh.get());

    //osg::ref_ptr<osg::StateSet> stateSet = geoNode->getOrCreateStateSet();
    stateSet = geoNode->getOrCreateStateSet();
    if (transparent)
    {
	stateSet->setMode(GL_DEPTH_TEST, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
	stateSet->setMode(GL_BLEND, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
    }
    else
    {
	stateSet->setMode(GL_DEPTH_TEST, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
	stateSet->setMode(GL_BLEND, osg::StateAttribute::OFF | osg::StateAttribute::OVERRIDE);
    }
    stateSet->setMode(GL_CULL_FACE, osg::StateAttribute::ON | osg::StateAttribute::OVERRIDE);
    osg::ref_ptr<osg::BlendFunc> blendFunc = new osg::BlendFunc(osg::BlendFunc::SRC_ALPHA, osg::BlendFunc::ONE_MINUS_SRC_ALPHA);
    stateSet->setAttributeAndModes(blendFunc.get());
    osg::ref_ptr<osg::CullFace> cullFace = new osg::CullFace(osg::CullFace::BACK);
    stateSet->setAttributeAndModes(cullFace.get());
    osg::ref_ptr<osg::FrontFace> frontFace = new osg::FrontFace(osg::FrontFace::COUNTER_CLOCKWISE);
    stateSet->setAttributeAndModes(frontFace.get());
    stateSet->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);
    worldNode->addChild(geoNode.get());

    if (!poseFilename.empty())
    {
	positions = readMatrixXFrom(poseFilename);
	SimCharacter character(jsonFilename);
	skeleton = character.skeleton;
	world->addSkeleton(character.skeleton);
	character.skeleton->setPositions(positions.row(frame));
    }

    /*
    for (size_t i = 0; i < character.skeleton->getNumDofs(); ++i)
	std::cout << i << "\t" << character.skeleton->getDof(i)->getName() << std::endl;
    */

    contacts = readVectorXdListFrom(cfFilename);
    group = new osg::Group();
    for (size_t i = 0; i < contacts[0].rows(); i += 6)
    {
	const Eigen::VectorXd &v = contacts[0];
	//osg::ref_ptr<osg::ShapeDrawable> shape = new osg::ShapeDrawable();
	//shape->setShape(new osg::Sphere(osg::Vec3(v[i + 0], v[i + 1], v[i + 2]), 0.02));
	//shape->setColor(osg::Vec4(0.0, 1.0, 0.0, 1.0));
	//group->addChild(shape);
	//shape = new osg::ShapeDrawable();
	//shape->setShape(new osg::Sphere(osg::Vec3(v[i + 0] + forceScale * v[i + 3], v[i + 1] + forceScale * v[i + 4], v[i + 2] + forceScale * v[i + 5]), 0.02));
	//shape->setColor(osg::Vec4(0.0, 1.0, 0.0, 1.0));
	//group->addChild(shape);
	osg::Vec3 point = osg::Vec3(v[i + 0], v[i + 1], v[i + 2]);
	osg::Vec3 force = osg::Vec3(v[i + 3], v[i + 4], v[i + 5]);
	//force = osg::Vec3(0, 0, 50);
	//std::cout << "point " << point[0] << " " << point[1] << " " << point[2] << std::endl;
	double sign = inverted ? -1 : 1;
	if (force.length() > epsilon)
	    group->addChild(makeArrow(point, point + force * sign * forceScale).get());
    }
    group->setUpdateCallback(new CustomCallBack());
    worldNode->addChild(group.get());

    viewer->addWorldNode(worldNode.get());
    viewer->getImGuiHandler()->addWidget(std::make_shared<TestWidget>(viewer, world));

    osg::Vec3d eye(0, 0, 5);
    osg::Vec3d center(0, 0, 0);
    osg::Vec3d up(0, 1, 0);
    viewer->getCameraManipulator()->setHomePosition(eye, center, up);
    viewer->addEventHandler(new CustomEventHandler(viewer));
    //viewer->record("output");
    viewer->run();

    if (Py_FinalizeEx() < 0) {
        return 120;
    }

    return 0;
}
