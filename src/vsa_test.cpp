//
// Created by Lind Xiao on 7/23/19.
//

#include "VSA_cgal.hpp"
#include <igl/read_triangle_mesh.h>
#include <iostream>
#include <Eigen/Core>
#include <map>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/jet.h>
#include <cxxopts.hpp>
#include <nlohmann/json.hpp>
#include "CutGraph.h"
igl::opengl::glfw::Viewer viewer;
using json = nlohmann::json;
int main(int argc, char *argv[]){
    using namespace std;
    if (argc < 2) {
        cout << "Usage cutmeshtest_bin --file mesh.obj --n sample_num --cut" << endl;
        exit(0);
    }
    cxxopts::Options options("CutMesh", "One line description of MyProgram");
    options.add_options()
            ("j, json","json storing parameters", cxxopts::value<std::string>())
            ;
    auto args = options.parse(argc, argv);
    json param_json;
    {
        std::ifstream temp(args["json"].as<std::string>());
        param_json= json::parse(temp);
        std::cout<< param_json<<std::endl;
    }
    Eigen::MatrixXd V;
    Eigen::MatrixXi F, fp;
    OTMapping::face_tuple M;
    std::map<OTMapping::face_tuple, int> map;
    int count;
    igl::read_triangle_mesh(param_json["mesh_file"], V, F);
    std::cout << F.rows()<< std::endl;
    OTMapping::vsa_compute(V, F, param_json["proxy_num"], fp, count);
    Eigen::MatrixXd C;
    CutGraph cg;
    cg._set_Vertices(V);
    cg._set_Edges_from_KNN(10);
    igl::jet(fp, 1,param_json["proxy_num"], C);
    viewer.data().set_mesh(V,F);
    viewer.data().set_colors(C);
    viewer.launch();
    return 0;
}