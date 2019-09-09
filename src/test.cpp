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
#include "mapping_patch.h"
#include "graphcut_cgal.h"
#include "mapping_patch.h"
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


Eigen::MatrixXi F_bad, F_good;
Eigen::MatrixXd V_bad, V_good;
Eigen::MatrixXi FL_bad, VL_good, FL_good;
Eigen::MatrixXd prob_mat;
int label_num;
igl::opengl::glfw::Viewer viewer;
Eigen::MatrixXd V_uv;
Eigen::MatrixXd V_i, bnd_uv;
Eigen::MatrixXi F_i;
Eigen::VectorXi bnd;
int lb_in=3;
std::vector<std::vector<bcclean::node>> vec;
std::vector<bcclean::node> ordered_nodes;
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
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier) {
    using namespace Eigen;
    using namespace std;
    switch(key) {
        case ']':{
            viewer.data().clear();
            lb_in+=1;
            lb_in %= 14;
            bcclean::mapping_patch mp;
            
            if(vec[lb_in].size()>=3){
                bcclean::extract_label_patch_mesh(V_good, F_good, FL_good, lb_in,V_i, F_i);
                
                
                if(!mp.build_patch(V_i, F_i,vec[lb_in],lb_in)){return false;}
                viewer.data().set_mesh(mp._V_uv, mp._F_uv);
                for(int nail:mp._nails){
                    bcclean::node nd = mp._nails_nodes_dict[nail];
                    viewer.data().add_points(nd._position, Eigen::RowVector3d(1,1,1));
                }
            }
            break;
        }
        
    }
    return false;
}

int main(){
    
    igl::read_triangle_mesh("../data/2/00000006_d4fe04f0f5f84b52bd4f10e4_trimesh_001.obj", V_bad, F_bad);
    igl::read_triangle_mesh("../data/2/bench.mesh__sf.obj", V_good, F_good);
    igl::readDMAT("../data/2/feat.dmat", FL_bad);
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
    vec = bcclean::build_label_nodes_list(V_good, F_good, FL_good);
    
    
    for(auto p: vec){
        std::cout << p.size() <<std::endl;
    }
    bcclean::extract_label_patch_mesh(V_good, F_good, FL_good, lb_in,V_i, F_i);
    std::cout<<"extracted"<<std::endl;
    std::vector<std::vector<bcclean::node>> vec1 = vec;
    std::vector<bcclean::node> list = vec[lb_in];
    std::cout<< list.size()<< std::endl;
    // bcclean::map_vertices_to_regular_polygon(V_i, F_i, list, bnd, bnd_uv, ordered_nodes);
    bcclean::mapping_patch mp;
    mp.build_patch(V_i, F_i,list,lb_in);
    // igl::harmonic(V_i,F_i,bnd, bnd_uv, 1, V_uv);
    viewer.callback_key_down = &key_down;
    
    viewer.data().set_mesh(mp._V_uv, mp._F_uv);
    for(bcclean::node nd:list){
        viewer.data().add_points(nd._position, Eigen::RowVector3d(1,1,1));
    }
    viewer.launch();
}