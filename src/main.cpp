//
// Created by Lind Xiao on 7/30/19.
//
#include <map>
#include <vector>
#include <iostream>
#include <pybind11/embed.h> // everything needed for embedding
#include <pybind11/pybind11.h>
#include <igl/read_triangle_mesh.h>
#include <igl/readDMAT.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_face_normals.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/embree/line_mesh_intersection.h>
#include <imgui/imgui.h>
#include <Eigen/Core>
#include "bcclean.h"
#include "graphcut_cgal.h"
#include <igl/upsample.h>
#include <igl/random_points_on_mesh.h>
#include <igl/jet.h>
#include <igl/boundary_facets.h>
#include <igl/readMSH.h>
#include <igl/remove_unreferenced.h>
#include <nanoflann.hpp>
#include <cxxopts.hpp>
#include <nlohmann/json.hpp>
namespace py = pybind11;
using json = nlohmann::json;
igl::opengl::glfw::Viewer viewer;

Eigen::MatrixXd V_bad, V_tet,V_good, V_good_refine; // mesh vertices
Eigen::MatrixXi F_bad, T_tet,F_good, F_good_refine; // mesh faces
Eigen::MatrixXd S_bad; // sample on two meshes
Eigen::MatrixXi I0, I1, I0_1, I1_1; // Sample source index to faces
Eigen::MatrixXd C_bad, C_good_cut, C_good_refine; // Sample color for mesh faces
Eigen::MatrixXi SL_bad, VL_good_refine;
Eigen::MatrixXi FL_bad, FL_good_refine, FL_good_cut;

Eigen::MatrixXd prob_mat;
Eigen::VectorXi b;
// constrained face id
int label_num;
double lambda_refine=1;
unsigned int file_idx=2;
int num_subdiv = 1;
double stop_energy=10;

std::map<int, std::vector<int> > patch_dict_bad;
// Currently selected face
int selected;



void update_states(unsigned int file_index){
    using namespace std;

    py::module sys = py::module::import("sys");
//    py::module torch = py::module::import("torch");
    py::module os = py::module::import("os");
    sys.attr("path").attr("append")("../src");
    py::module utlis = py::module::import("utlis");
    std::cout << "2"<< std::endl;
    // handle and generate all the necessary files
    std::string bad_mesh_file, good_mesh_file, tet_file, face_label_dmat, face_label_yml;
    py::object py_bad_mesh_file =
            utlis.attr("kth_existing_file")("../data/ABC/obj/abc_0000_obj_v00",file_index);
    bad_mesh_file = py_bad_mesh_file.cast<string>();
    py::object py_tet_file = utlis.attr("tetrahedralize_bad_mesh")(bad_mesh_file,"", stop_energy);
    std::cout << "3"<< std::endl;
    tet_file = py_tet_file.cast<string>();
    py::object py_face_label_yml = utlis.attr("kth_existing_file")("../data/ABC/feat/abc_0000_feat_v00",file_index);
    face_label_yml = py_face_label_yml.cast<string>();
    py::object py_face_label_dmat = utlis.attr("parse_feat")(face_label_yml);
    face_label_dmat = py_face_label_dmat.cast<string>();
    std::cout << "4"<< std::endl;
    // read the files into eigen matrices
    igl::read_triangle_mesh(bad_mesh_file, V_bad, F_bad);
    std::cout << "5"<< std::endl;
    igl::readMSH(tet_file, V_tet, T_tet);
    {
        Eigen::MatrixXi F_good_temp;
        igl::boundary_facets(T_tet, F_good_temp);
        Eigen::VectorXi J, K;
        igl::remove_unreferenced(V_tet, F_good_temp, V_good, F_good, J, K);
    }
    std::cout << "6"<< std::endl;
    try {
        igl::upsample(V_good, F_good, V_good_refine, F_good_refine, num_subdiv);
    } catch(...){
        std::cout << "failed to upsample by" << num_subdiv << "subdivisions" <<std::endl;
        V_good_refine = V_good;
        F_good_refine = F_good;
    }

    std::cout << "7"<< std::endl;
    igl::readDMAT(face_label_dmat, FL_bad);
    bcclean::build_patch_dict(FL_bad, patch_dict_bad);
    std::cout << "8"<< std::endl;
    label_num = FL_bad.maxCoeff()+1;
    int sample_num = V_good_refine.rows();
    std::cout << sample_num<< std::endl;
    {
        // handle all the labels
        Eigen::MatrixXi I_bad;
        Eigen::SparseMatrix<double> B_bad;
        igl::random_points_on_mesh(sample_num, V_bad, F_bad, B_bad, I_bad);
        S_bad = B_bad * V_bad;
        // generate random sample on bad mesh;
        std::cout << "9"<< std::endl;
        bcclean::generate_sample_label(FL_bad, I_bad, sample_num, SL_bad);
    }

    igl::jet(FL_bad, 0, label_num-1, C_bad);

}

void update_NN(){
    // generate the sample label on bad mesh;
    std::cout << "10"<< std::endl;
    bcclean::NN_sample_label_transport(S_bad, V_good_refine, SL_bad, VL_good_refine);
    // nearest neighbor transport sample labels to vertices;
    std::cout << "11"<< std::endl;
    bcclean::vertex_label_vote_face_label(label_num, VL_good_refine, F_good_refine, FL_good_refine, prob_mat);
    std::cout << "12"<< std::endl;
    igl::jet(FL_good_refine, 0, label_num-1, C_good_refine);
}

void update_intersection(){
    bcclean::LM_intersection_label_transport(
    V_bad,
    F_bad,
    FL_bad,
    V_good_refine,
    F_good_refine,
    VL_good_refine);
    bcclean::vertex_label_vote_face_label(label_num, VL_good_refine, F_good_refine, FL_good_refine, prob_mat);
    std::cout << "12"<< std::endl;
    igl::jet(FL_good_refine, 0, label_num-1, C_good_refine);
}

void refine_cuts(){
    FL_good_cut = FL_good_refine;
    bcclean::refine_labels_graph_cut(V_good_refine,F_good_refine, prob_mat.transpose(), FL_good_cut, lambda_refine);
    std::cout << "13"<< std::endl;
    igl::jet(FL_good_cut, 0, label_num-1, C_good_cut);
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
    using namespace Eigen;
    using namespace std;
    switch(key) {
        case '1':{
            viewer.data().clear();
            viewer.data().set_mesh(V_bad, F_bad);
            viewer.data().set_colors(C_bad);
            break;
        }
        case '2': {
            viewer.data().clear();
            update_NN();
            viewer.data().set_mesh(V_good_refine, F_good_refine);
            viewer.data().set_colors(C_good_refine);
            break;

        }
        case '3':{
            viewer.data().clear();
            update_intersection();
            viewer.data().set_mesh(V_good_refine, F_good_refine);
            viewer.data().set_colors(C_good_refine);
            break;
        }
        case '4':{
            viewer.data().clear();
            refine_cuts();
            viewer.data().set_mesh(V_good_refine, F_good_refine);
            viewer.data().set_colors(C_good_cut);
            break;
        }
        case '[': {
            if (file_idx == 0){
                file_idx = 99;
            }
            else (file_idx = (file_idx-1)% 100);
            num_subdiv = 1;
            update_states(file_idx);
            std::cout <<"new idx"<< file_idx << std::endl;
            break;
        }
        case ']': {
            file_idx = (file_idx+1)% 100;
            num_subdiv = 1;
            update_states(file_idx);
            std::cout <<"new idx"<< file_idx << std::endl;
            break;
        }
    }

    return false;
}

bool mouse_down(igl::opengl::glfw::Viewer &viewer, int, int) {
    int fid_ray;
    Eigen::Vector3f bary;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
                                 viewer.core().proj, viewer.core().viewport, V_bad, F_bad, fid_ray, bary)) {
        bool found = false;
        for (int i = 0; i < b.size(); ++i) {
            if (b(i) == fid_ray) {
                found = true;
                selected = i;
            }
        }

        if (!found) {
            b.conservativeResize(b.size() + 1);
            b(b.size() - 1) = fid_ray;
        }


        return true;
    }
    return false;
};


int main(int argc, char *argv[]) {
    using namespace std;
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
    lambda_refine=param_json["lambda_refine"];
    file_idx=param_json["file_idx"];
    num_subdiv = param_json["num_subdiv"];
    stop_energy= param_json["stop_energy"];
    py::scoped_interpreter guard{};
    std::cout << "1"<< std::endl;
    update_states(file_idx);
    viewer.callback_key_down = &key_down;
    viewer.launch();
}