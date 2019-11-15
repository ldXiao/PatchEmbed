//
// Created by Lind Xiao on 6/11/19.
//

#ifndef OTMAPPING_QUADRATIC_OT_H
#define OTMAPPING_QUADRATIC_OT_H
#include <cxxopts.hpp>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include <vector>
#include <math.h>
#include <memory>
#include <fstream>
#include <utility>
#include <nlohmann/json.hpp>
#include <igl/copyleft/cgal/mesh_to_polyhedron.h>
#include <igl/copyleft/cgal/polyhedron_to_mesh.h>
#include "graphcut_cgal.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd C;
Eigen::MatrixXd L;
using json = nlohmann::json;
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
    igl::read_triangle_mesh(param_json["mesh_file"], V, F);
    int cluster_num = param_json["cluster_num"];
    double smoothing_lambda = param_json["smoothing_lambda"];
    sdf_segmentation(V, F, L, cluster_num, smoothing_lambda);
    std::vector<double> v= sdf_computation(V, F);
//    Eigen::VectorXd C1;
    Eigen::VectorXd C1 = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(v.data(), v.size());
    std::cout << C1 << std::endl;
    std::cout << "------------------------" << std::endl;
    for(auto it = v.begin(); it!= v.end(); ++it){
        std::cout << *it << std::endl;
    }
    viewer.data().set_mesh(V,F);
    viewer.data().set_colors(L);
    viewer.launch();
}

#endif //OTMAPPING_QUADRATIC_OT_H
