//
// Created by Lind Xiao on 6/6/19.
//
#include "CutMesh.h"
#include "cxxopts.hpp"
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <math.h>
#include <memory>
#include <utility>
// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;
OTMapping::CutMesh CM;
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
    using namespace Eigen;
    using namespace std;
    if(key != 'i') {
        CM.plot_CutMesh(viewer, key);
    }
    return false;
}

using namespace std;
using namespace OTMapping;
int main(int argc, char *argv[]){
    using namespace std;
    using namespace Eigen;
    if (argc < 2) {
        cout << "Usage cutmeshtest_bin --file mesh.obj --n sample_num --cut" << endl;
        exit(0);
    }
    cxxopts::Options options("CutMesh", "One line description of MyProgram");
    options.add_options()
            ("f, file", "File name", cxxopts::value<std::string>())
            ("n, num","Sample number", cxxopts::value<int>())
            ("c, cut","cut the mesh", cxxopts::value<bool>())
            ;
    auto args = options.parse(argc, argv);
    // Load a mesh in OBJ format
    igl::readOBJ(args["file"].as<std::string>(), V, F);
    igl::opengl::glfw::Viewer viewer;
    int sample_num=args["n"].as<int>();
    auto func = [](Eigen::Vector3d x)->double{return std::sin(x[0]+x[1]+x[2]);};
    CM.set_initial(V, F, sample_num , func);
    std::cout << "holy shit1" << std::endl;
    if(args["cut"].as<bool>()) {
        Eigen::Vector3d p0(0, 0, 0);
        Eigen::Vector3d p1(0, 0, 0);
        Eigen::Vector3d n1(0, 0, 1);
        Eigen::Vector3d n0(0, 1, 0.5);
        CM.cut_with(p0, n0, p1, n1);
    }
    else{
        CM.separate_cube_faces();
    }
    CM.perturb(3, 0.1, 0.02);
    CM.plot_CutMesh(viewer,'i');
    viewer.callback_key_down = &key_down;
    std::cout << "holy shit4" << std::endl;
    viewer.launch();
}
