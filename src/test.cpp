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

Eigen::MatrixXi F_bad, F_good;
Eigen::MatrixXd V_bad, V_good;
Eigen::MatrixXi FL_bad, VL_good, FL_good;
Eigen::MatrixXd prob_mat;
int label_num;

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
        bcclean::LM_intersection_label_transport(V_bad,F_bad,FL_bad,V_good,F_good,VL_good);
        std::cout << "r2"<< std::endl;
        // bcclean::Barycenter_intersection_label_transport(V_bad,F_bad,FL_bad,V_good,F_good,FL_good);
        bcclean::vertex_label_vote_face_label(label_num, VL_good, F_good, FL_good, prob_mat);
        std::cout << "r3"<< std::endl;
}

int main(){
    
    igl::read_triangle_mesh("../data/test.obj", V_bad, F_bad);
    igl::read_triangle_mesh("../data/test_good.obj", V_good, F_good);
    igl::readDMAT("../data/test.dmat", FL_bad);
    label_num = FL_bad.maxCoeff()+1;
    bcclean::LM_intersection_label_transport(
    V_bad,
    F_bad,
    FL_bad,
    V_good,
    F_good,
    VL_good);
    bcclean::vertex_label_vote_face_label(label_num, VL_good, F_good, FL_good, prob_mat);
    bcclean::refine_labels_graph_cut(V_good, F_good, prob_mat.transpose(), FL_good, 1);
    std::cout <<"finishe"<<std::endl;
    std::vector<std::vector<bcclean::node>> vec = bcclean::build_label_nodes_list(V_good, F_good, FL_good);
    for(auto p: vec){
        std::cout << p.size() <<std::endl;
    }
}