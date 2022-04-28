#include <iostream>
#include <fstream>
#include <vector>
#include <nlohmann/json.hpp>
#include <string>
#include <dart/dart.hpp>
#include <dart/gui/osg/osg.hpp>
#include <dart/external/imgui/imgui.h>
#include <Python.h>
#include "SimCharacter.h"
#include "PyUtil.h"

int frame = 0;
Eigen::MatrixXd positions;
SimCharacter character("data/character.json");

class TestWidget : public dart::gui::osg::ImGuiWidget
{
public:

    TestWidget(dart::gui::osg::ImGuiViewer* viewer, dart::simulation::WorldPtr world) : mViewer(viewer), mWorld(world)
    {
    }

    void render() override
    {
	ImGui::SetNextWindowPos(ImVec2(10, 20));
	ImGui::SetNextWindowSize(ImVec2(1000, 200));
	ImGui::Begin("Control");
	ImGui::SliderInt("frame", &frame, 0, 500);
	character.skeleton->setPositions(positions.row(frame));
	ImGui::End();
    }

private:

    osg::ref_ptr<dart::gui::osg::ImGuiViewer> mViewer;
    dart::simulation::WorldPtr mWorld;

};


int main()
{
    Py_Initialize();
    PyRun_SimpleString("import sys\n"
                       "print('python version:', sys.version)\n");

    PyObject* pDict = load_npz("data/0005_Walking001_poses.npz");
    //print_py_obj(pDict);
    PyObject* npa;
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

    dart::simulation::WorldPtr world = dart::simulation::World::create();
    world->addSkeleton(character.skeleton);
    character.skeleton->setPositions(positions.row(frame));
    osg::ref_ptr<dart::gui::osg::ImGuiViewer> viewer = new dart::gui::osg::ImGuiViewer();
    osg::ref_ptr<dart::gui::osg::WorldNode> worldNode = new dart::gui::osg::RealTimeWorldNode(world);
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
