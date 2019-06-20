//
// Created by Lind Xiao on 6/6/19.
//
#include "CutMesh.h"
#include <cxxopts.hpp>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <math.h>
#include <memory>
#include <fstream>
#include <utility>
#include <nlohmann/json.hpp>
// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;
OTMapping::CutMesh CM;
using json = nlohmann::json;
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
            ("j, json","json storing parameters", cxxopts::value<std::string>())
            ;
    auto args = options.parse(argc, argv);
    // Load a mesh in OBJ format
    igl::opengl::glfw::Viewer viewer;
    json param_json;
    {
        std::ifstream temp(args["json"].as<std::string>());
        param_json= json::parse(temp);
        std::cout<< param_json<<std::endl;
    }
    CM.set_initial_from_json(param_json);
    CM.plot_CutMesh(viewer,'i');
    viewer.callback_key_down = &key_down;
    std::cout << "holy shit4" << std::endl;
    viewer.launch();
}
