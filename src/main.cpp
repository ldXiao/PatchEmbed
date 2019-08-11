//
// Created by Lind Xiao on 6/9/19.
//
#include <iostream>
#include "CutMesh.h"
#include "CutGraph.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/random_points_on_mesh.h>
#include <igl/read_triangle_mesh.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/jet.h>
#include <igl/edges.h>
#include <igl/loop.h>
#include <igl/upsample.h>
#include <cxxopts.hpp>
#include <Eigen/Core>
#include <map>
#include <vector>
#include "VSA_cgal.hpp"
#include "graphcut_cgal.h"
#include <nlohmann/json.hpp>
using json = nlohmann::json;
using namespace OTMapping;
igl::opengl::glfw::Viewer viewer;
Eigen::MatrixXd V0, V1, V1_1; // mesh vertices
Eigen::MatrixXi F0, F1, F1_1; // mesh faces
Eigen::MatrixXd S0, S1, S0_1, S1_1; // sample on two meshes
Eigen::MatrixXi I0, I1, I0_1, I1_1; // Sample source index to faces
Eigen::MatrixXd C0, C1, C0_1, C1_1; // Sample color for mesh faces
Eigen::MatrixXd FC0, FC1, FC1_1;
Eigen::MatrixXi SL0, SL1, SLf, SL0_1 ,SL1_1; // Sample Label
Eigen::MatrixXd heads;
Eigen::MatrixXd tails;
//                    igl::slice(this->SamplePerturb,this->SkeletonIndices1,1)
int proxy_num;
int label_num;
int point_size=7;
double lambda_refine;
CutGraph CtGrph, CtGrph1;
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
        case '1': {
            // show  mesh knn graph edgweight and nn transported labels
            viewer.data().clear();
            viewer.data().point_size = point_size;
            viewer.data().set_mesh(V1, F1);
            viewer.data().set_colors(Eigen::RowVector3d(1, 1, 0));
            igl::jet(SL1, 0, proxy_num-1, C1);
            Eigen::MatrixXd C3;
            igl::jet(CtGrph.EdgeWeights, 0, 2, C3);
            viewer.data().add_points(S1, C1);
            heads = igl::slice(S1, CtGrph.Edges.col(0), 1);
            tails = igl::slice(S1, CtGrph.Edges.col(1), 1);
            viewer.data().add_edges(heads, tails, C3);
            break;
        }
        case '2': {
            // show graph cut on knn graph
            viewer.data().clear();
            viewer.data().point_size = point_size;
            viewer.data().set_mesh(V1, F1);
            viewer.data().set_colors(Eigen::RowVector3d(1, 1, 1));
            CtGrph._set_ProbabilityMatrix(label_num, SL1);
            SLf = SL1;
            std::cout << "initial label" << SLf << std::endl;
            std::cout << CtGrph.EdgeWeights << std::endl;
            OTMapping::Alpha_expansion_graph_cut_boykov_kolmogorov_Eigen(
                    CtGrph.Edges,
                    CtGrph.EdgeWeights,
                    CtGrph.ProbabilityMatrix,
                    SLf
            );
            igl::jet(SLf, 0, SLf.maxCoeff(), C1);
            std::cout << SLf << std::endl;
            viewer.data().add_points(S1, C1);
            heads = igl::slice(S1, CtGrph.Edges.col(0), 1);
            tails = igl::slice(S1, CtGrph.Edges.col(1), 1);
            viewer.data().add_edges(heads, tails, Eigen::RowVector3d(0,0,0));
            break;
        }

        case '3': {
            viewer.data().clear();
            viewer.data().point_size = point_size;
            viewer.data().set_mesh(V1, F1);
            Eigen::MatrixXi FL1_mod2;
            Eigen::MatrixXd probmat2;
            NN_sample_label_vote_face_label(label_num, I1, SL1, F1, FL1_mod2, probmat2);
            igl::jet(FL1_mod2, 0, label_num-1, FC1);
            viewer.data().set_colors(FC1);
            break;
        }
        case '4': {
            viewer.data().clear();
            viewer.data().point_size = point_size;
            viewer.data().set_mesh(V1, F1);
            Eigen::MatrixXi FL1_mod;
            Eigen::MatrixXd probmat;
            NN_sample_label_vote_face_label(label_num, I1, SL1, F1, FL1_mod, probmat);
            OTMapping::refine_labels_graph_cut(V1, F1, probmat.transpose(), FL1_mod, lambda_refine);
            igl::jet(FL1_mod, 0, label_num-1, FC1);
            viewer.data().set_colors(FC1);
            break;
        }
        case '5': {
            // show vertices based initilization
            viewer.data().clear();
            viewer.data().point_size = point_size;
            viewer.data().set_mesh(V1_1, F1_1);
            viewer.data().set_colors(Eigen::RowVector3d(1, 1, 1));
            igl::jet(SL1_1, 0, SL1_1.maxCoeff(), C1_1);
            viewer.data().add_points(S1_1, C1_1);
            break;
        }
        case '6': {
            // vertices-based graph cut
            viewer.data().clear();
            viewer.data().point_size = point_size;
            viewer.data().set_mesh(V1_1, F1_1);
            viewer.data().set_colors(Eigen::RowVector3d(1, 1, 1));
            CtGrph1._set_ProbabilityMatrix(label_num, SL1_1);
            SLf = SL1_1;
            OTMapping::Alpha_expansion_graph_cut_boykov_kolmogorov_Eigen(
                    CtGrph1.Edges,
                    CtGrph1.EdgeWeights,
                    CtGrph1.ProbabilityMatrix,
                    SLf
            );
            Eigen::MatrixXd C5;
            igl::jet(SLf, 0, SLf.maxCoeff(), C5);
            viewer.data().add_points(S1_1, C5);
            break;
        }
        case '7':{
            // show
            viewer.data().clear();
            viewer.data().point_size = point_size;
            Eigen::MatrixXi FL1_1_mod;
            Eigen::MatrixXd probmat;
            viewer.data().set_mesh(V1_1, F1_1);
            vertex_label_vote_face_label(label_num, SL1_1, F1_1, FL1_1_mod, probmat);
            igl::jet(FL1_1_mod, 0, label_num-1, FC1_1);
            viewer.data().set_colors(FC1_1);
            break;
        }

        case '8':{
            // show
            viewer.data().clear();
            viewer.data().point_size = point_size;
            viewer.data().set_mesh(V1_1, F1_1);
            Eigen::MatrixXi FL1_1_mod;
            Eigen::MatrixXd probmat;
            vertex_label_vote_face_label(label_num, SL1_1, F1_1, FL1_1_mod, probmat);
            OTMapping::refine_labels_graph_cut(V1_1, F1_1, probmat.transpose(), FL1_1_mod, lambda_refine);
            igl::jet(FL1_1_mod, 0, label_num-1, FC1_1);
            viewer.data().set_colors(FC1_1);
            break;
        }
        case '9':{
            viewer.data().clear();
            viewer.data().point_size = point_size;
            viewer.data().set_mesh(V1_1, F1_1);
            Eigen::MatrixXi FL1_1_mod;
            Eigen::MatrixXd probmat;
            vertex_label_vote_face_label(label_num, SL1_1, F1_1, FL1_1_mod, probmat);
            Eigen::MatrixXd FW;
            Eigen::MatrixXi Edg;
            set_Face_Edges(F1_1, Edg);
            set_EdgeWeight(lambda_refine, FL1_1_mod, Edg, FW);
            std::vector<std::pair<std::size_t, std::size_t> > a_edges;
            std::vector<double> edge_weights;
//            OTMapping::calculate_and_log_normalize_dihedral_angles(
//                    V1_1,
//                    F1_1,
//                    lambda_refine,
//                    a_edges,
//                    edge_weights);
            OTMapping::Alpha_expansion_graph_cut_boykov_kolmogorov_Eigen(
                    Edg,
                    FW,
                    probmat,
                    FL1_1_mod
            );
            igl::jet(FL1_1_mod, 0, label_num-1, FC1_1);
            viewer.data().set_colors(FC1_1);
            break;
        }
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
    std::cout << "h1"<< std::endl;
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
    std::cout << "h2"<< std::endl;
    Eigen::MatrixXi FL0; // face_label for mesh 0
    {
        // use vsa to generate face color and sample color for V0, C0, S0
        int count;

        OTMapping::vsa_compute(V0, F0, param_json["proxy_num"], FL0, count);
        label_num = FL0.maxCoeff()+1;
        igl::jet(FL0, 0, label_num-1, FC0);
//        generate_sample_color(FC0, I0,param_json["sample_num"], C0);
        generate_sample_label(FL0, I0, param_json["sample_num"], SL0);
        igl::jet(SL0, 0, label_num-1, C0);

        NN_sample_label_transport(S0, S1, SL0, SL1);
        igl::jet(SL1, 0, label_num-1, C1);
    }
    std::cout << "h3"<< std::endl;
    {
        std::string upsample_method = param_json["upsample"];
        int num_subdiv = param_json["num_subdiv"];
        if(upsample_method != "loop"){
            igl::upsample(V1,F1, V1_1, F1_1, num_subdiv);
            std::cout << "ha1"<< std::endl;
            int nsample_num = V1_1.rows();
            S1_1 = V1_1;
            Eigen::SparseMatrix<double> B0_1;
            igl::random_points_on_mesh(nsample_num, V0, F0, B0_1, I0_1);
            std::cout << "ha2"<< std::endl;
            S0_1 = B0_1 * V0;
            generate_sample_label(FL0, I0_1, S0_1.rows(), SL0_1);
            std::cout << "ha3"<< std::endl;
            igl::jet(SL0_1, 0, label_num-1, C0_1);
            NN_sample_label_transport(S0_1, S1_1, SL0_1, SL1_1);
            std::cout << "ha4"<< std::endl;
            igl::jet(SL1_1, 0, label_num-1, C1_1);
            std::cout << "ha5"<< std::endl;
        }
    }
    std::cout << "ha"<< std::endl;
//    CtGrph._set_Vertices(S0);

    CtGrph.initialize(
            S1,
            param_json["KNN_valance"],
            label_num,
            param_json["lambda"],
            SL1
            );
    CtGrph1.initialize2(
            S1_1,
            F1_1,
            label_num,
            param_json["lambda"],
            SL1_1
            );
    lambda_refine = param_json["lambda"];
    viewer.callback_key_down = &key_down;
    viewer.launch();
}
