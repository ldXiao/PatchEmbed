#include "patch.h"
#include "planar_cut_simply_connect.h"
#include "patch_cut_relabel.h"
#include "non_vertex_manifold_relabel.h"
#include <unordered_map>
#include <igl/boundary_loop.h>
#include <igl/remove_unreferenced.h>
#include <igl/is_vertex_manifold.h>
namespace bcclean{
    // all the static member must be defined out side the class scope
    // because they are intrinsically not part of thec class
    Eigen::MatrixXd patch::Vbase;
    Eigen::MatrixXi patch::Fbase;
    Eigen::MatrixXd patch::V_mod;
    Eigen::MatrixXi patch::F_mod;
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
        std::unordered_map<int, std::vector<bool> > patch_edge_direc_dict;
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
        // when collecting the patches , it could split the edges whose two nodes are on boundaries while itself is not part of boundary
        if(patch::total_label_num<=0){
            std::cout << "Static members not initialized, abort" <<std::endl;
            exit(EXIT_FAILURE);
        }
        // else do 
        std::map<int, int> label_count_dict; // count how many faces are there in a label patch recorded for resizing only
        std::map<int, Eigen::MatrixXi> label_faces_dict;
        std::map<int, Eigen::VectorXi> label_FI_dict;
        for(int lb = 0; lb < patch::total_label_num; ++lb){
            label_count_dict[lb] =0;
            label_faces_dict[lb] = Eigen::MatrixXi::Zero(patch::Fbase.rows(),3);
            label_FI_dict[lb] = Eigen::VectorXi::Zero(patch::Fbase.rows());
        }
        for(int fidx =0 ; fidx < patch::FL.rows(); ++fidx){
            int lb = patch::FL(fidx);
            label_faces_dict[lb].row(label_count_dict[lb])=patch::Fbase.row(fidx);
            label_FI_dict[lb](label_count_dict[lb])= fidx;
            label_count_dict[lb]+=1;
        }
        for(int lb = 0; lb < patch::total_label_num; ++lb){
            label_faces_dict[lb].conservativeResize(label_count_dict[lb],3);
            label_FI_dict[lb].conservativeResize(label_count_dict[lb]);
        }
        int total_label_num_dummy = patch::total_label_num;
        for(int lb =0; lb < patch::total_label_num; ++ lb){
            
            std::cout  <<"lb" << lb << std::endl;
            Eigen::MatrixXi Fraw;
            Eigen::MatrixXd Vraw;
            Eigen::VectorXi I, J;
            igl::remove_unreferenced(patch::Vbase, label_faces_dict[lb], Vraw, Fraw, I, J);
            // check is is non-manifold


            // if  Vraw Fraw is a manifold, cut open the loops
            Eigen::VectorXi NonManifoldVertices;
            if(igl::is_vertex_manifold(Fraw, NonManifoldVertices))
            {
                patch pat;
                // check Vraw, Fraw is a manifold
                // check the number of boundary loops of the manifold
                // if more than one it is not simply connected.
                pat.Vraw = Vraw;
                pat.Fraw = Fraw;
                pat.FI = label_FI_dict[lb];
                pat.label = lb;
                pat.VI = J; // VI()
                std::vector<std::vector<int > > boundary_loops;
                igl::boundary_loop(Fraw, boundary_loops);

                
                if(boundary_loops.size()>1){
                    std::vector<bool> VCuts;
                    std::vector<std::vector<bool> > TCuts;
                    std::cout << "label" << lb << "has " << boundary_loops.size() <<"loops"<< std::endl;
                    for(auto p : boundary_loops){
                        for(auto vi: p){
                            std::cout << vi << " ";
                        }
                        // std::cout << "loop" << std::endl;
                    }
                    planar_cut_simply_connect(Vraw, Fraw, patch::Vbase, patch::Fbase, pat.VI, pat.FI, patch::FL_mod, boundary_loops, VCuts, TCuts);
                    // std::cout << "ob1" << lb << std::endl;
                    Eigen::VectorXi FL_mod_copy =patch::FL_mod;
                    patch_cut_relabel(Fraw, pat.FI, VCuts, TCuts, FL_mod_copy, patch::FL_mod, total_label_num_dummy);
                }
            //
            }
            else{
                std::cout << "patch "<< lb << " is non-manfold, first level relabelling" << std::endl;
                std::vector<int> NMV;
                for(int j =0 ; j < NonManifoldVertices.rows(); j++)
                {
                    if(NonManifoldVertices(j)==0)
                    {
                        NMV.push_back(j);
                    }
                }
                Eigen::VectorXi FL_mod_copy =patch::FL_mod;
                std::map<int, Eigen::VectorXi> subpatchFI_dict;
                Prepocess::non_vertex_manifold_relabel(Vraw, Fraw, label_FI_dict[lb], NMV, FL_mod_copy, patch::FL_mod, total_label_num_dummy, subpatchFI_dict);
                for(auto item: subpatchFI_dict)
                {
                    int lb = item.first;
                    Eigen::VectorXi FI_lb = item.second;
                    Eigen::MatrixXi F_lb = Eigen::MatrixXi::Constant(Fraw.rows(),3, -1);
                    int count = 0;
                    for(int f =0 ; f < FI_lb.rows(); ++f)
                    {
                        F_lb.row(f) = patch::Fbase.row(FI_lb(f));
                        count +=1;
                    }
                    patch pat;
                    // check Vraw, Fraw is a manifold
                    // check the number of boundary loops of the manifold
                    // if more than one it is not simply connected
                    pat.FI = FI_lb;
                    pat.label = lb;
                    {
                        Eigen::VectorXi B;
                        igl::remove_unreferenced(patch::Vbase, F_lb, pat.Vraw, pat.Fraw, pat.FI, pat.VI);
                        assert(igl::is_vertex_manifold(pat.Fraw,B));
                    }
                    pat.VI = J; 
                    std::vector<std::vector<int> > boundary_loops;
                    igl::boundary_loop(pat.Fraw, boundary_loops);
                    if(boundary_loops.size()>1)
                    {
                        std::vector<bool> VCuts;
                        std::vector<std::vector<bool> > TCuts;
                        std::cout << "subpatch" << lb << "has " << boundary_loops.size() <<"loops"<< std::endl;
                        for(auto p : boundary_loops){
                            for(auto vi: p){
                                std::cout << vi << " ";
                            }
                            // std::cout << "loop" << std::endl;
                        }
                        planar_cut_simply_connect(pat.Vraw, pat.Fraw, patch::Vbase, patch::Fbase, pat.VI, pat.FI, patch::FL_mod, boundary_loops, VCuts, TCuts);
                        // std::cout << "ob1" << lb << std::endl;
                        Eigen::VectorXi FL_mod_copy =patch::FL_mod;
                        patch_cut_relabel(pat.Fraw, pat.FI, VCuts, TCuts, FL_mod_copy, patch::FL_mod, total_label_num_dummy);
                    }

                }

            }
        }
        std::cout << "patch relabel finished" << std::endl;
        patch::total_label_num = total_label_num_dummy;
    }
}