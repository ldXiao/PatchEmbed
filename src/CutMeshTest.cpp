//
// Created by Lind Xiao on 6/6/19.
//
#include "CutMesh.h"
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
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
    std::cout << "hello shit" << std::endl;
    // Load a mesh in OBJ format
    igl::readOBJ(argv[1], V, F);
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.launch();
}
