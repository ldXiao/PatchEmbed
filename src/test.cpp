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
#include <igl/opengl/ViewerData.h>
#include <igl/per_face_normals.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/embree/line_mesh_intersection.h>
#include <imgui/imgui.h>
#include <igl/jet.h>
#include <igl/writeOBJ.h>
#include <igl/upsample.h>
#include <imgui/imgui_internal.h>
#include <Eigen/Core>
#include "bcclean.h"
#include "edge.h"
#include "graphcut_cgal.h"
#include "patch.h"
#include "Match_Maker.h"
#include <igl/upsample.h>
#include <igl/random_points_on_mesh.h>
#include <igl/jet.h>
#include <igl/boundary_facets.h>
#include <igl/readMSH.h>
#include <igl/remove_unreferenced.h>
#include <igl/harmonic.h>
#include <nanoflann.hpp>
#include <cxxopts.hpp>
#include <nlohmann/json.hpp>
#include <unordered_map>

Eigen::MatrixXi F_bad, F_good;
Eigen::MatrixXd V_bad, V_good;
Eigen::MatrixXi FL_bad, VL_good;
Eigen::VectorXi II, JJ;
Eigen::VectorXi FL_good;
Eigen::MatrixXd prob_mat;
int label_num;
igl::opengl::glfw::Viewer viewer;
Eigen::MatrixXd V_uv;
Eigen::MatrixXd V_i, bnd_uv;
Eigen::MatrixXi F_i;
Eigen::VectorXi bnd;
int lb_in=3;
bcclean::pair_map<std::pair<int,int>,std::vector<int>>  pair_edge_list_dict_bad, pair_edge_list_dict_good;
std::vector<bcclean::edge> edge_list_bad, edge_list_good;
std::unordered_map<int, std::vector<int> > patch_edge_bad, patch_edge_good;
std::unordered_map<int, std::vector<int> > vertices_label_list_bad, vertices_label_list_good;
Eigen::MatrixXd C_bad, C_good;
bool visual_bad=false;
bool iden = false;
int upsp = 0;
std::vector<int> node_list_extern;
void project_face_labels(
    const Eigen::MatrixXd &V_bad, 
    const Eigen::MatrixXi &F_bad, 
    const Eigen::MatrixXi &FL_bad,
    const Eigen::MatrixXd &V_good,
    const Eigen::MatrixXi &F_good,
    Eigen::MatrixXi & FL_good,
    Eigen::MatrixXd & prob_mat){
        int label_num = 0;
        std::map<int, int> count_dict;
        for(int i =0; i< FL_bad.rows(); ++i){
            auto it = count_dict.find(FL_bad(i,0));
            if(it == count_dict.end()){
                label_num+=1;
                count_dict[FL_bad(i,0)]=1;
            }
        }
        Eigen::MatrixXi VL_good, FL_good_temp;
        std::cout << "r1"<< std::endl;
        // bcclean::LM_intersection_label_transport(V_bad,F_bad,FL_bad,V_good,F_good,VL_good);
        std::cout << "r2"<< std::endl;
        // bcclean::Barycenter_intersection_label_transport(V_bad,F_bad,FL_bad,V_good,F_good,FL_good);
        // bcclean::vertex_label_vote_face_label(label_num, VL_good, F_good, FL_good, prob_mat);
        bcclean::refine_proj_vote(V_bad, F_bad, FL_bad, V_good, F_good, label_num, 2, FL_good, prob_mat);
        std::cout << "r3"<< std::endl;
}
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
    using namespace Eigen;
    using namespace std;
    switch(key) {
        case ']':{
            if(visual_bad){
                viewer.data().clear();
                // igl::jet(FL_bad.unaryExpr([](const int x) { return x%8; }), 0, 8, C_bad);
                // igl::jet(FL_good.unaryExpr([](const int x) { return x%8; }), 0, 8, C_good);
                viewer.data().show_overlay_depth = false;
                viewer.data().set_mesh(V_good, F_good);
                // viewer.data().set_colors(C_bad);
                
                for(auto p: edge_list_bad){
                bcclean::plot_edge(viewer, V_bad, FL_bad, p);
                }
                visual_bad = !visual_bad;
            } else {
                viewer.data().clear();
                viewer.data().show_overlay_depth = false;
                viewer.data().set_mesh(V_good, F_good);
                viewer.data().set_colors(C_good);
                
                for(auto p: edge_list_good){
                bcclean::plot_edge(viewer, V_good, FL_good, p);
                }
                visual_bad = !visual_bad;
            }
            break;
        }
    case '=':{
        Eigen::MatrixXd V_good_copy = V_good;
        Eigen::MatrixXi F_good_copy = F_good;
        Eigen::VectorXi FL_good_copy = FL_good;
        bcclean::MatchMaker::trace_and_label(bcclean::patch::Vbase, bcclean::patch::Fbase, bcclean::patch::FL_mod, V_good_copy, F_good_copy, FL_good_copy);
        viewer.data().clear();
        viewer.data().set_mesh(V_good_copy, F_good_copy);
        Eigen::MatrixXd head, tail;
        igl::slice(V_good_copy, II, 1, head);
        igl::slice(V_good_copy, JJ, 1,tail);
        viewer.data().add_edges(head, tail , Eigen::RowVector3d(1,0,0));
        break;
    }
            
    }
    return false;
}
std::string data_root;

using json = nlohmann::json;
int main(int argc, char *argv[]){
    using namespace std;
    namespace py = pybind11;
    if (argc < 2) {
        std::cout << "Usage cutmeshtest_bin --file mesh.obj --n sample_num --cut" << std::endl;
        exit(0);
    }
    
    /*-----------------------------------------
    for json
    */
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
        data_root = param_json["data_root"];
        iden = param_json["iden"];
        upsp = param_json["upsp"];
    }
    /*-----------------------------------------
    for pybind
    */

    py::scoped_interpreter guard{};
    py::module sys = py::module::import("sys");
    py::module os = py::module::import("os");
    sys.attr("path").attr("append")("../src");
    py::module utlis = py::module::import("utlis");
    // handle and generate all the necessary files
    std::string bad_mesh_file, good_mesh_file, face_label_dmat, face_label_yml;
    py::object py_bad_mesh_file =
            utlis.attr("get_bad_mesh")(data_root);
    bad_mesh_file = py_bad_mesh_file.cast<string>();
    // py::object py_face_label_yml = utlis.attr("get_feat_file")(data_root);
    // face_label_yml = py_face_label_yml.cast<string>();
    // py::object py_face_label_dmat = utlis.attr("parse_feat")(face_label_yml);
    // face_label_dmat = py_face_label_dmat.cast<string>();
    good_mesh_file = data_root + "/"+"good.mesh__sf.obj";
    face_label_dmat = data_root + "/"+ "feat.dmat";
    igl::read_triangle_mesh(bad_mesh_file, V_bad, F_bad);
    if(iden)
    {
    igl::read_triangle_mesh(bad_mesh_file, V_good, F_good);
    }
    else
    {
        igl::read_triangle_mesh(good_mesh_file, V_good, F_good);
    }
    if(upsp> 0)
    {
        igl::upsample(V_good, F_good, upsp);
    }
    igl::readDMAT(face_label_dmat, FL_bad);
    label_num = FL_bad.maxCoeff()+1;
    bcclean::patch pat;
    bcclean::patch::SetStatics(V_bad, F_bad, FL_bad, label_num);
    std::map<int, bcclean::patch> patch_dict;
    bcclean::CollectPatches();
    viewer.data().set_mesh(V_good, F_good);
    viewer.callback_key_pressed = key_down;
    viewer.launch();
    // bcclean::trace_and_label(bcclean::patch::Vbase, bcclean::patch::Fbase, bcclean::patch::FL_mod, V_good, F_good, FL_good, dots, viewer);
    
    igl::jet(bcclean::patch::FL_mod,0, bcclean::patch::total_label_num-1, C_bad);
    igl::writeDMAT(data_root +"/test.dmat", bcclean::patch::FL_mod );
    igl::writeOBJ(data_root + "/test.obj", V_good, F_good);
    // viewer.data().set_colors(C_bad);
    
    return 0;
}