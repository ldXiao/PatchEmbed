//
// Created by Lind Xiao on 6/6/19.
//
#include "CutMesh.h"
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <math.h>
// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;


using namespace std;
int main(int argc, char *argv[]){
    using namespace std;
    using namespace Eigen;
    if (argc != 2) {
        cout << "Usage cutmeshtest_bin mesh.obj" << endl;
        exit(0);
    }
    std::cout << "holy shit" << std::endl;
    // Load a mesh in OBJ format
    igl::readOBJ(argv[1], V, F);
    igl::opengl::glfw::Viewer viewer;
    CutMesh CM;
    auto func = [](Eigen::Vector3d x)->double{return std::sin(x[0]+x[1]+x[2]);};
    CM.set_initial(V, F, 1000, func);
    CM.plot_CutMesh(viewer,"aha");
    viewer.launch();
}
