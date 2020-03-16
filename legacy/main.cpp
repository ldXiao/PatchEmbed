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
#include <igl/writeDMAT.h>
#include <igl/writeOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/per_face_normals.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/embree/line_mesh_intersection.h>
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>
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
Eigen::MatrixXi I_bad;
Eigen::MatrixXd C_bad, C_good_cut, C_good_refine; // Sample color for mesh faces
Eigen::MatrixXi SL_bad, VL_good_refine;
Eigen::MatrixXi FT_bad, FL_bad, FL_good_refine, FL_good_cut;
// FT_bad store the total feature tags, 
//while FL_bad store the patch label after manually merging
Eigen::MatrixXd prob_mat;
bool clickable= true;
bool invert_normal= false;
unsigned int current_show_mode ='1';
// only set true when it is bad mesh



int label_num, patch_num, sample_num;
double lambda_refine=1;
unsigned int file_idx=2;
int num_subdiv = 1;
double stop_energy=10;

int selected_label = 4;
int selected_face_index = -1;
int selected_patch_index = -1; 
bool selected; // activated when mouse click a patch
std::map<int, std::vector<int> > patch_dict_bad; //  to re-tag patches
std::vector<std::string> patch_dict_bad_labels; // for use in menu combo


void modular_patch_label(
    const Eigen::MatrixXi & FT_bad, 
    const int &label_number, 
    Eigen::MatrixXi & FL_bad){
        FL_bad = FT_bad;
    for(int i=0; i< FL_bad.rows(); ++i){
        FL_bad(i,0)= FT_bad(i,0) % label_number;
    }
}

void update_states(std::string data_root){
    using namespace std;

    py::module sys = py::module::import("sys");
    py::module os = py::module::import("os");
    sys.attr("path").attr("append")("../src");
    py::module utlis = py::module::import("utlis");
    // handle and generate all the necessary files
    std::string bad_mesh_file, good_mesh_file, tet_file, face_label_dmat, face_label_yml;
    py::object py_bad_mesh_file =
            utlis.attr("get_bad_mesh")(data_root);
    bad_mesh_file = py_bad_mesh_file.cast<string>();
    std::cout << "tetrahedralizing..."<< std::endl;
    py::object py_tet_file = utlis.attr("tetrahedralize_bad_mesh")(bad_mesh_file,"", stop_energy);
    tet_file = py_tet_file.cast<string>();
    py::object py_face_label_yml = utlis.attr("get_feat_file")(data_root);
    face_label_yml = py_face_label_yml.cast<string>();
    py::object py_face_label_dmat = utlis.attr("parse_feat")(face_label_yml);
    face_label_dmat = py_face_label_dmat.cast<string>();
    // read the files into eigen matrices
    std::cout << "reading bad mesh..."<< std::endl;
    igl::read_triangle_mesh(bad_mesh_file, V_bad, F_bad);
    std::cout << "reading tet mesh..."<< std::endl;
    igl::readMSH(tet_file, V_tet, T_tet);
    {
        std::cout << "generating surface good mesh..."<< std::endl;
        Eigen::MatrixXi F_good_temp;
        igl::boundary_facets(T_tet, F_good_temp);
        Eigen::VectorXi J, K;
        igl::remove_unreferenced(V_tet, F_good_temp, V_good, F_good, J, K);
    }
    try {
        igl::upsample(V_good, F_good, V_good_refine, F_good_refine, num_subdiv);
    } catch(...){
        std::cout << "failed to upsample by" << num_subdiv << "subdivisions" <<std::endl;
        std::cout << "use original instead..."<< std::endl;
        V_good_refine = V_good;
        F_good_refine = F_good;
    }

    igl::readDMAT(face_label_dmat, FT_bad);
    bcclean::build_patch_dict(FT_bad, patch_dict_bad);
    patch_dict_bad_labels.clear();
    for(int i=0; i < label_num; ++i) {
        patch_dict_bad_labels.push_back(std::to_string(i));
    }

    FL_bad = FT_bad;
    patch_num = FT_bad.maxCoeff()+1;
    for(int i=0; i< FL_bad.rows(); ++i){
        FL_bad(i,0)= FT_bad(i,0) % label_num;
    }
    label_num = FL_bad.maxCoeff()+1;
    sample_num = V_good_refine.rows();
    std::cout << sample_num<< std::endl;
    {
        // handle all the labels
        Eigen::SparseMatrix<double> B_bad;
        igl::random_points_on_mesh(sample_num, V_bad, F_bad, B_bad, I_bad);
        S_bad = B_bad * V_bad;
        // generate random sample on bad mesh;
        bcclean::generate_sample_label(FL_bad, I_bad, sample_num, SL_bad);
    }

    igl::jet(FL_bad, 0, label_num-1, C_bad);
}



void update_NN(){
    // generate the sample label on bad mesh;
    bcclean::generate_sample_label(FL_bad, I_bad, sample_num, SL_bad);
    bcclean::NN_sample_label_transport(S_bad, V_good_refine, SL_bad, VL_good_refine);
    // nearest neighbor transport sample labels to vertices;
    bcclean::vertex_label_vote_face_label(label_num, VL_good_refine, F_good_refine, FL_good_refine, prob_mat);
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
    igl::jet(FL_good_refine, 0, label_num-1, C_good_refine);
}

void refine_cuts(){
    FL_good_cut = FL_good_refine;
    bcclean::refine_labels_graph_cut(V_good_refine,F_good_refine, prob_mat.transpose(), FL_good_cut, lambda_refine);
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
    }
    viewer.data().invert_normals = invert_normal;

    return false;
}

void plot_selected_patch_on_bad_mesh(
    igl::opengl::glfw::Viewer &viewer,
    const Eigen::MatrixXd & V_bad, 
    const Eigen::MatrixXi & F_bad, 
    const Eigen::MatrixXi & FL_bad, 
    std::map<int, std::vector<int> >  & patch_dict_bad, 
    int selected_patch_idx){
        Eigen::MatrixXd temp_C_bad; // only work when patch selected
        igl::jet(FL_bad,0, label_num-1, temp_C_bad);
        auto v = patch_dict_bad[selected_patch_idx];
        for(int f: v){
            temp_C_bad.row(f) = Eigen::RowVector3d(1,1,1);
        }
        viewer.data().clear();
        viewer.data().set_mesh(V_bad, F_bad);
        viewer.data().set_colors(temp_C_bad);
}

bool mouse_down(igl::opengl::glfw::Viewer &viewer, int, int) {
    int fid_ray;
    Eigen::Vector3f bary;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    if(clickable){
        // mouse_down only work on bad_mesh
        if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
                                 viewer.core().proj, viewer.core().viewport, V_bad, F_bad, fid_ray, bary)) 
        {
            selected_face_index = fid_ray;
            selected_patch_index = FT_bad(selected_face_index,0);
            selected_label = FL_bad(selected_face_index,0);
            plot_selected_patch_on_bad_mesh(viewer, V_bad, F_bad, FL_bad, patch_dict_bad, selected_patch_index);
            selected = true;   // activate selected state
            return true;
        }
    }
    return false;
};

void reset_label_num(int new_label_num){
    label_num = new_label_num;
    modular_patch_label(FT_bad, label_num, FL_bad);
}


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
        std::cout <<"the json parameters are" << param_json << std::endl;
    }
    lambda_refine=param_json["lambda_refine"];
    num_subdiv = param_json["num_subdiv"];
    stop_energy= param_json["stop_energy"];
    label_num = param_json["label_num"];
    std::string data_root = param_json["data_root"];
    py::scoped_interpreter guard{};
    update_states(data_root);
    viewer.callback_key_down = &key_down;
    viewer.callback_mouse_down = &mouse_down;
    // link key_down function
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);
//
    menu.callback_draw_viewer_menu = [&]() {
        if (ImGui::CollapsingHeader("States parameters", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::InputDouble("lambda_refine", &lambda_refine);
            ImGui::InputInt("num_subdiv", &num_subdiv);
            ImGui::InputDouble("stop_energy", &stop_energy);
            ImGui::InputInt("label_num", &label_num);
            if (ImGui::Button("update states", ImVec2(0,0))){
                update_states(data_root);
            }
        }
        if(ImGui::CollapsingHeader("Label manipulation", ImGuiTreeNodeFlags_DefaultOpen)){
            int selected_label_num = std::distance(
                patch_dict_bad_labels.begin(), 
                std::find(patch_dict_bad_labels.begin(), patch_dict_bad_labels.end(), std::to_string(selected_label))
                );
            if (ImGui::Combo("Label patch as", &selected_label_num, patch_dict_bad_labels)){
                if(selected){
                    selected_label = std::stoi(patch_dict_bad_labels[selected_label_num]);
                    for (auto f:patch_dict_bad[selected_patch_index]){
                        FL_bad(f,0) = selected_label; 
                    }
                    igl::jet(FL_bad, 0, label_num-1, C_bad);
                    key_down(viewer, '1',0);
                    selected = false;
                }
            }
        }
        if(ImGui::CollapsingHeader("Visualization choices", ImGuiTreeNodeFlags_DefaultOpen)){
            if(ImGui::Button("flip normal", ImVec2(0,0))){
                invert_normal = ! invert_normal;
                key_down(viewer, current_show_mode, 0);
            }
            if (ImGui::Button("show bad mesh", ImVec2(-1,0))) {
                key_down(viewer, '1', 0);
                current_show_mode = '1';
                clickable = true;
            }
            if (ImGui::Button("show NN on good mesh", ImVec2(-2,0))){
                key_down(viewer, '2', 0);
                current_show_mode = '2';
                clickable = false;
            }
            if (ImGui::Button("show projection on good mesh", ImVec2(-3,0))){
                key_down(viewer, '3', 0);
                current_show_mode = '3';
                clickable = false;
            }
            if(ImGui::Button("graph cut", ImVec2(-4,0))){
                key_down(viewer, '4', 0);
                current_show_mode = '4';
                clickable = false;
            }
        }
        if(ImGui::CollapsingHeader("Outputs", ImGuiTreeNodeFlags_DefaultOpen)){
            if (ImGui::Button("Export good mesh", ImVec2(-1,0))) {
                igl::writeOBJ(data_root+"/good_mesh.obj", V_good_refine, F_good_refine);
            }
            if (ImGui::Button("Export feat file(dmat)", ImVec2(-2,0))){
                igl::writeDMAT(data_root+"/good_feat.dmat", FL_good_refine);
            }
        }
    };
    viewer.launch();
}