//
// Created by Lind Xiao on 6/9/19.
//
#include <iostream>
#include "CutMesh.h"
#include "CutGraph.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/random_points_on_mesh.h>
#include <CGAL/internal/Surface_mesh_segmentation/Alpha_expansion_graph_cut.h>
#include <igl/read_triangle_mesh.h>
#include <igl/jet.h>
#include <cxxopts.hpp>
#include <Eigen/Core>
#include "VSA_cgal.hpp"
#include "graphcut_cgal.h"
#include <nlohmann/json.hpp>
using json = nlohmann::json;
igl::opengl::glfw::Viewer viewer;
Eigen::MatrixXd V0, V1; // mesh vertices
Eigen::MatrixXi F0, F1; // mesh faces
Eigen::MatrixXd S0, S1; // sample on two meshes
Eigen::MatrixXi I0, I1; // Sample source index to faces
Eigen::MatrixXd C0, C1; // Sample color for mesh faces
Eigen::MatrixXd FC0, FC1;
Eigen::MatrixXi SL0, SL1, SLf; // Sample Label
Eigen::MatrixXd heads;
Eigen::MatrixXd tails;
//                    igl::slice(this->SamplePerturb,this->SkeletonIndices1,1)
int proxy_num;
int label_num;
int point_size=7;
CutGraph CtGrph;
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
    using namespace Eigen;
    using namespace std;
    switch(key){
        case '0':
            viewer.data().clear();
            viewer.data().point_size = point_size;
            viewer.data().set_mesh(V0, F0);
            viewer.data().add_points(S0, C0);
            viewer.data().set_colors(FC0);
            break;
        case '1':
            viewer.data().clear();
            viewer.data().point_size = point_size;
            viewer.data().set_mesh(V1,F1);
            viewer.data().set_colors(Eigen::RowVector3d(1,1,0));
//            std::cout << S1 << std::endl;
            igl::jet(SL1, 1, proxy_num, C1);
            viewer.data().add_points(S1, C1);
            heads = igl::slice(S1, CtGrph.Edges.col(0),1);
            tails = igl::slice(S1, CtGrph.Edges.col(1),1);
//                    igl::slice(this->SamplePerturb,this->SkeletonIndices1,1)
            viewer.data().add_edges(heads, tails, Eigen::RowVector3d(1,1,1));
            break;
        case '2':
            viewer.data().clear();
            viewer.data().point_size = point_size;
            viewer.data().set_mesh(V1,F1);
            viewer.data().set_colors(Eigen::RowVector3d(1,1,1));
            CtGrph._set_ProbabilityMatrix(label_num, SL1);
            SLf = SL1;
            OTMapping::Alpha_expansion_graph_cut_boykov_kolmogorov_Eigen(
                    CtGrph.Edges,
                    CtGrph.EdgeWeights,
                    CtGrph.ProbabilityMatrix,
                    SLf
                    );
            igl::jet(SLf,1, label_num, C1);
            std::cout << SLf << std::endl;
            viewer.data().add_points(S1,C1);
            heads = igl::slice(S1, CtGrph.Edges.col(0),1);
            tails = igl::slice(S1, CtGrph.Edges.col(1),1);
//                    igl::slice(this->SamplePerturb,this->SkeletonIndices1,1)
//            viewer.data().add_edges(heads, tails, Eigen::RowVector3d(0,0,0));
            break;
    }

    return false;
}

int main(int argc, char *argv[]) {
    using namespace std;
    unsigned int sample_num;
    if (argc < 2) {
        cout << "Usage cutmeshtest_bin --file mesh.obj --n sample_num --cut" << endl;
        exit(0);
    }
    cxxopts::Options options("CutGraph", "One line description of MyProgram");
    options.add_options()
            ("j, json", "json storing parameters", cxxopts::value<std::string>());
    auto args = options.parse(argc, argv);
    // Load a mesh in OBJ format
    json param_json;
    {
        std::ifstream temp(args["json"].as<std::string>());
        param_json = json::parse(temp);
        std::cout << param_json << std::endl;
    }
    igl::read_triangle_mesh(param_json["mesh_file0"], V0, F0);
    igl::read_triangle_mesh(param_json["mesh_file1"], V1, F1);
    sample_num = param_json["sample_num"];
    proxy_num = param_json["proxy_num"];
    {
        Eigen::SparseMatrix<double> B0, B1;
        igl::random_points_on_mesh(sample_num, V0, F0, B0, I0);
        igl::random_points_on_mesh(sample_num, V1, F1, B1, I1);
        S0 = B0 * V0;
        S1 = B1 * V1;
    }

    {
        // use vsa to generate face color and sample color for V0, C0, S0
        int count;
        Eigen::MatrixXi FL; // face_label for mesh 0
        OTMapping::vsa_compute(V0, F0, param_json["proxy_num"], FL, count);
        label_num = FL.maxCoeff()+1;
        igl::jet(FL, 1, param_json["proxy_num"], FC0);
//        generate_sample_color(FC0, I0,param_json["sample_num"], C0);
        generate_sample_label(FL, I0, param_json["sample_num"], SL0);
        igl::jet(SL0, 1, param_json["proxy_num"], C0);
        NN_sample_label_transport(S0, S1, SL0, SL1);
        igl::jet(SL1, 1, param_json["proxy_num"], C1);
    }
    std::cout << "ha"<< std::endl;
//    CtGrph._set_Vertices(S0);

    CtGrph.initialize(
            S1,
            param_json["KNN_valance"],
            label_num,
            param_json["lambda"]
            );

    viewer.callback_key_down = &key_down;
    viewer.launch();
}
