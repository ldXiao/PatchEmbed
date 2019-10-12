#include "patch.h"
#include "planar_cut_simply_connect.h"
#include "patch_cut_relabel.h"
#include <unordered_map>
#include <igl/boundary_loop.h>
#include <igl/remove_unreferenced.h>
namespace bcclean{
    // all the static member must be defined out side the class scope
    // because they are intrinsically not part of thec class
    Eigen::MatrixXd patch::Vbase;
    Eigen::MatrixXi patch::Fbase;
    std::vector<node> patch::node_list; // list of all nodes
    std::vector<edge> patch::edge_list; // list of all edges;
    int patch::total_label_num;
    Eigen::VectorXi patch::FL; // a copy of initial face labels
    Eigen::VectorXi patch::FL_mod; // Face labels after modification

    void patch::SetStatics(const Eigen::MatrixXd & Vbase_in, const Eigen::MatrixXi & Fbase_in, const Eigen::VectorXi & FL_in, size_t total_label_num_in){
        Vbase = Vbase_in;
        Fbase = Fbase_in;
        FL = FL_in;
        FL_mod = FL;
        total_label_num =total_label_num_in;
        std::unordered_map<int, std::vector<int> > patch_edge_dict;
        bcclean::build_edge_list(Vbase, Fbase, FL, total_label_num, edge_list, patch_edge_dict);
        // initialize edge_list
        std::unordered_map<int, std::vector<int>> _vertex_label_list_dict;
        bcclean::build_vertex_label_list_dict(Fbase, FL, total_label_num, _vertex_label_list_dict);
        for(auto item: _vertex_label_list_dict){
            if(item.second.size()>2){
                // it is a node
                node nd;
                nd.initialize(total_label_num, item.first, Vbase.row(item.first), item.second);
                node_list.push_back(nd);
            }
        }
    }

    void CollectPatches(){
        // this function must be called after SetSattics
        // check all set otherwise exit with failure
        if(patch::total_label_num<=0){
            std::cout  << "Static members not initialized abort" <<std::endl;
            exit(EXIT_FAILURE);
        }
        std::cout << "1" << std::endl;
        // else do 
        std::map<int, int> label_count_dict;
        std::map<int, Eigen::MatrixXi> label_faces_dict;
        std::map<int, Eigen::VectorXi> label_FI_dict;
        for(int lb = 0; lb < patch::total_label_num; ++lb){
            label_count_dict[lb] =0 ;
            label_faces_dict[lb] = Eigen::MatrixXi::Zero(patch::Fbase.rows(),3);
            label_FI_dict[lb] = Eigen::VectorXi::Zero(patch::Fbase.rows());
        }
        std::cout << "2" << std::endl;
        for(int fidx =0 ; fidx < patch::FL.rows(); ++fidx){
            int lb = patch::FL(fidx);
            label_faces_dict[lb].row(label_count_dict[lb])=patch::Fbase.row(fidx);
            label_FI_dict[lb](label_count_dict[lb])= fidx;
            label_count_dict[lb]+=1;
        }
        std::cout << "3" << std::endl;
        for(int lb = 0; lb < patch::total_label_num; ++lb){
            label_faces_dict[lb].conservativeResize(label_count_dict[lb],3);
            label_FI_dict[lb].conservativeResize(label_count_dict[lb]);
        }
        std::cout << "4" << std::endl;
        int total_label_num_dummy =patch::total_label_num;
        for(int lb =0; lb < patch::total_label_num; ++ lb){
            patch pat;
            std::cout  <<"lb" << lb << std::endl;
            Eigen::MatrixXi Fraw;
            Eigen::MatrixXd Vraw;
            Eigen::VectorXi I;
            igl::remove_unreferenced(patch::Vbase, label_faces_dict[lb], Vraw, Fraw, I);

            // check Vraw, Fraw is a manifold
            // check the number of boundary loops of the manifold
            // if more than one it is not simply connected.
            pat.Vraw = Vraw;
            pat.Fraw = Fraw;
            pat.FI = label_FI_dict[lb];
            pat.label = lb;
            pat.VI = I; // VI()
            std::vector<std::vector<int > > boundary_loops;
            igl::boundary_loop(Fraw, boundary_loops);

            
            if(boundary_loops.size()>1){
                std::vector<bool> VCuts;
                std::cout << "ob" << lb << std::endl;
                planar_cut_simply_connect(Vraw, Fraw, boundary_loops, VCuts);
                std::cout << "ob1" << lb << std::endl;
                Eigen::VectorXi FL_mod_copy =patch::FL_mod;
                patch_cut_relabel(Fraw, label_FI_dict[lb], VCuts, FL_mod_copy, patch::FL_mod, total_label_num_dummy);
            }
            //
        }
        patch::total_label_num = total_label_num_dummy;
    }
}