//
// Created by Lind Xiao on 7/25/19.
//

#include "bcclean.h"
#include "kdtree_NN_Eigen.hpp"
#include <Eigen/Core>
#include <iostream>
#include <map>
#include <tuple>
#include <igl/edges.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/sparse.h>
#include <igl/embree/line_mesh_intersection.h>
#include <igl/per_vertex_normals.h>
#include <igl/harmonic.h>
#include <igl/boundary_loop.h>
#include <igl/remove_unreferenced.h>
#include <random>
#include <vector>
#include <map>
#include <cmath>

#define PI 3.141592653

namespace bcclean {

    void generate_sample_color(
            const Eigen::MatrixXd &FC,
            const Eigen::MatrixXi &source_index,
            const int &sample_num,
            Eigen::MatrixXd &SC) {
        SC.resize(sample_num, 3);
        for (int i = 0; i < sample_num; ++i) {
            SC.row(i) = FC.row(source_index(i, 0));
        }

    }

    void generate_sample_label(
            const Eigen::MatrixXi &FL,
            const Eigen::MatrixXi &source_index,
            const int &sample_num,
            Eigen::MatrixXi &SL
    ) {
        SL.resize(sample_num, 3);
        for (int i = 0; i < sample_num; ++i) {
            SL(i, 0) = FL(source_index(i, 0), 0);
        }
    }

    void NN_sample_label_transport(
            const Eigen::MatrixXd &S0,
            const Eigen::MatrixXd &S1,
            const Eigen::MatrixXi &SL0,
            Eigen::MatrixXi &SL1) {
        int sample_num = S1.rows();
        SL1.resize(SL0.rows(), 1);
        bcclean::kd_tree_Eigen<double> kdt(S0.cols(),std::cref(S0),10);
        kdt.index->buildIndex();
        for (unsigned int i = 0; i < sample_num; ++i) {
            int min_idx = bcclean::kd_tree_NN_Eigen<double>(kdt, S1.row(i));
            SL1(i, 0) = SL0(min_idx, 0);
        }
    }


    void Barycenter_intersection_label_transport(
                const Eigen::MatrixXd &V0,
                const Eigen::MatrixXi & F0,
                const Eigen::MatrixXi & FL0,
                const Eigen::MatrixXd & V1,
                const Eigen::MatrixXi & F1,
                Eigen::MatrixXi & FL1){
            Eigen::MatrixXd N1;
            igl::per_face_normals(V1, F1, N1);
            Eigen::MatrixXd BaryCenters = Eigen::MatrixXd::Constant(F1.rows(), 3, 0);
            for(int j =0; j< F1.rows();++j){
                int ii, jj, kk;
                ii = F1(j,0);
                jj = F1(j,1);
                kk = F1(j,2);
                BaryCenters.row(j)= (V1.row(ii)+ V1.row(jj)+ V1.row(kk))/3;
            }
            Eigen::MatrixXd R1=igl::embree::line_mesh_intersection(BaryCenters, N1, V0, F0);
            FL1= Eigen::MatrixXi::Constant(F1.rows(), 1, -1);
            for(int i =0; i < F1.rows(); ++i){
                FL1(i,0)= FL0(std::round(R1(i,0)),0);
            }
        }

    void LM_intersection_label_transport(
            const Eigen::MatrixXd &V0,
            const Eigen::MatrixXi & F0,
            const Eigen::MatrixXi & FL0,
            const Eigen::MatrixXd & V1,
            const Eigen::MatrixXi & F1,
            Eigen::MatrixXi & VL1){
        Eigen::MatrixXd N1;
        std::cout << "r3.1"<< std::endl;
        igl::per_vertex_normals(V1, F1, N1);
        std::cout << "r4"<< std::endl;
        Eigen::MatrixXd R1=igl::embree::line_mesh_intersection(V1, N1, V0, F0);
        std::cout << "r5"<< std::endl;
        VL1= Eigen::MatrixXi::Constant(V1.rows(), 1, -1);
        int short_mem = 0;
        for(int i =0; i < VL1.rows(); ++i){
            int j = std::round(R1(i,0));
            if(j==-1){VL1(i,0) = short_mem;}
            else{VL1(i,0)= FL0(j,0); short_mem = FL0(j,0);}
        }
    }

    void construct_face_sample_dictionary(
            const Eigen::MatrixXi &I1,
            std::map<int, std::vector<int> > &dict) {
        for (int i = 0; i < I1.rows(); ++i) {
            auto find = dict.find(I1(i, 0));
            if (find != dict.end()) {
                dict[I1(i, 0)].push_back(i);
            } else {
                dict[I1(i, 0)] = std::vector<int>();
            }
        }
    }

    void NN_sample_label_vote_face_label(
            const int label_num,
            const Eigen::MatrixXi &I1,
            const Eigen::MatrixXi &SL1,
            const Eigen::MatrixXi &F1,
            Eigen::MatrixXi &FL1,
            Eigen::MatrixXd &probability_matrix
    ) {
        std::map<int, std::vector<int> > dict;
        construct_face_sample_dictionary(I1, dict);
        FL1 = Eigen::MatrixXi::Zero(F1.rows(), 1);
        std::cout << dict.size() << "size" << std::endl;
        probability_matrix = Eigen::MatrixXd::Constant(F1.rows(), label_num, 0.5);
        for (int fidx = 0; fidx < FL1.rows(); ++fidx) {
            auto find = dict.find(fidx);
            if (find != dict.end()) {
                auto vect = dict[fidx];
                int size = dict[fidx].size();
                int max_num_label = 0;
                int max_label_num = 0;
                for (auto &sp: vect) {
                    int lbl = SL1(sp, 0);
                    probability_matrix(fidx, lbl) = probability_matrix(fidx, lbl) + 1;
                    if (probability_matrix(fidx, lbl) > max_label_num) {
                        max_label_num = probability_matrix(fidx, lbl);
                        max_num_label = lbl;
                    }
                }
                double sum = probability_matrix.row(fidx).sum();
                probability_matrix.row(fidx) /= sum;
                FL1(fidx, 0) = max_num_label;
            } else {
                probability_matrix.row(fidx) = Eigen::MatrixXd::Constant(1, label_num, 1 / label_num);
                FL1(fidx, 0) = fidx % label_num;
            }
        }
    }

    void vertex_label_vote_face_label(
            const int label_num,
            const Eigen::MatrixXi &VL,
            const Eigen::MatrixXi &F,
            Eigen::MatrixXi &FL,
            Eigen::MatrixXd &prob_mat
    ) {
        FL = Eigen::MatrixXi::Zero(F.rows(), 1);
        prob_mat = Eigen::MatrixXd::Constant(FL.rows(), label_num, 0.1 / label_num);
        std::map<int, int> dict;
        for (int fidx = 0; fidx < F.rows(); ++fidx) {
            // int v0, v1, v2, l0, l1, l2;
            dict.clear();
            for (int i = 0; i < 3; ++i) {
                int vi = F(fidx, i);
                int li = VL(vi, 0);
                auto find = dict.find(li);
                if (find != dict.end()) {
                    dict[li] += 1;
                } else {
                    dict[li] = 1;
                }
            }
            std::cout << "built dict" <<std::endl;
            if (dict.size() == 3) {
                FL(fidx, 0) = dict.begin()->first;
                auto it = dict.begin();
                int l0 = it->first;
                it++;
                int l1 = it->first;
                it++;
                int l2 = it->first;

                prob_mat(fidx, l0) = 0.3;
                prob_mat(fidx, l1) = 0.3;
                prob_mat(fidx, l2) = 0.3;

            }
            std::cout << "b1" <<std::endl;
            if (dict.size() == 2) {
                auto it = dict.begin();
                int l0 = it->first;
                it++;
                int l1 = it->first;
                if (dict[l0] > dict[l1]) {
                    FL(fidx, 0) = l0;
                    prob_mat(fidx, l0) = 0.6;
                    prob_mat(fidx, l1) = 0.3;
                } else {
                    FL(fidx, 0) = l1;
                    prob_mat(fidx, l1) = 0.6;
                    prob_mat(fidx, l0) = 0.3;
                }
            }
            std::cout << "b2" <<std::endl;
            if (dict.size() == 1) {
                FL(fidx, 0) = dict.begin()->first;
                int l0 = dict.begin()->first;
                prob_mat(fidx, l0) = 0.9;
            }
        }


    }

    void
    set_EdgeWeight(const double &lambda, const Eigen::MatrixXi &SL, const Eigen::MatrixXi &E, Eigen::MatrixXd &EW) {
        EW = Eigen::MatrixXd::Constant(E.rows(), 1, lambda);
        std::map<int, int> boundary_vertices;
        for (int i = 0; i < E.rows(); ++i) {
            int j = E(i, 0);
            int k = E(i, 1);
            if (SL(j, 0) == SL(k, 0)) {
                continue;
            } else {
                EW(i, 0) = 0.01 * lambda;
            }
        }
    }

    void set_Face_Edges(const Eigen::MatrixXi &F, Eigen::MatrixXi &E) {
        Eigen::MatrixXi TT;
        igl::triangle_triangle_adjacency(F, TT);
        E = Eigen::MatrixXi::Constant(3 * F.rows(), 2, -1);
        int count = 0;
        for (int i = 0; i < TT.rows(); i++) {
            for (int j = 0; j < 3; j++) {
                if (TT(i, j) == -1) continue;
                int of = TT(i, j);
                E(count, 0) = i;
                E(count, 1) = of;
                count += 1;
            }
        }
        E.conservativeResize(count, 2);
    }

    void normalize_mesh(Eigen::MatrixXd & V_bad, Eigen::MatrixXd & V_good){
        Eigen::RowVector3d bb_min, bb_max, bb_center;
        bb_min = V_bad.colwise().minCoeff();
        bb_max = V_bad.colwise().maxCoeff();
        bb_center = (bb_min + bb_max)/2;
        double diag_lenth = (bb_max-bb_min).norm();
        V_bad.rowwise() -= bb_center;
        V_bad * 5/ diag_lenth;
        V_good.rowwise() -= bb_center;
        V_good * 5/ diag_lenth;
    }

    void build_patch_dict(const Eigen::MatrixXi &FL, std::map<int, std::vector<int> > & patch_dict){
        patch_dict.clear();
        for(int fidx =0 ; fidx < FL.rows(); fidx++){
            int patch_idx = FL(fidx,0);
//            std::vector<int> chunk;
            auto it = patch_dict.find(patch_idx);
            if(it == patch_dict.end()){
                std::vector<int> chunk;
                chunk.push_back(fidx);
                patch_dict[patch_idx]= chunk;
            }
            else{
                patch_dict[patch_idx].push_back(fidx);
            }
        }
    }

    bool node::initialize(
        const int total_label_num, 
        const Eigen::MatrixXd & position, 
        const std::vector<int> labels){
            _total_label_num = total_label_num;
            _occupied_label_num = 0;
            for(int i =0; i< _total_label_num; ++i){
                _label_occupy_dict[i]=0;
            }
            for(auto lb: labels){
                if(lb<_total_label_num && lb>-1){
                    _label_occupy_dict[lb]=1;
                }
                else{
                    std::cout<< "label out of range"<<std::endl;
                    return false;
                }
            }
            if(position.rows() == 1 and position.cols()==3){
                    this->_position = position;
                    return true;
            }
            else{
                    std::cout<< "position should be initialized as rowvector3d" <<std::endl;
                    return false;
            }
    };
    bool node::of_same_type(const node & b){
            if(_total_label_num != b._total_label_num){
                return false;
            }
            if(_occupied_label_num != b._occupied_label_num){
                return false;
            }
            for(int i =0; i< _total_label_num; ++i){
               if( _label_occupy_dict[i]!=b._label_occupy_dict.at(i)){
                   return false;
               }
            }
            return true;
    };
    bool node::at_same_position(const Eigen::MatrixXd& position){
        double tol = 10e-6;
        if ((position.row(0)-_position.row(0)).norm()< tol){
            return true;
        }
        else{
            return false;
        }
    }

    std::vector<std::vector<node>> build_label_nodes_list(
        const Eigen::MatrixXd &V, 
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXi &FL){
            // return lb->vector<nodes> where the nodes are unodered for each patch
            std::map<int, std::vector<int>> count_dict;
            // vertx_index -> label_list
            for(int fidx=0; fidx< F.rows(); ++fidx){
                int lb = FL(fidx,0);
                for(int i=0; i <3; ++i){
                    int vi = F(fidx, i);
                    auto & vect = count_dict[vi];
                    auto it = std::find (vect.begin(), vect.end(), lb);
                    if(it == vect.end()){
                        std::cout<<"push0"<<vi<<","<<lb<<std::endl;
                        vect.push_back(lb);
                        std::cout<<vect.size()<<std::endl;
                    }
                }
            }
            int total_label_num = FL.maxCoeff()+1;
            std::vector<std::vector<node>> result(total_label_num);
            for(auto & [vi, vect]: count_dict){
                std::sort(vect.begin(), vect.end());
                if(vect.size()>2){
                    Eigen::RowVector3d position = V.row(vi);
                    node nd;
                    nd.initialize(total_label_num, position, vect);
                    for(int lb: vect){
                        result[lb].push_back(nd);
                        std::cout<< "push"<<std::endl;
                    }
                }
            }
            return result;
        };

        void extract_label_patch_mesh(
            const Eigen::MatrixXd& V, 
            const Eigen::MatrixXi& F, 
            const Eigen::MatrixXi&FL, 
            const int lb_in, 
            Eigen::MatrixXd& V_i, 
            Eigen::MatrixXi& F_i){
                Eigen::MatrixXi F_l = Eigen::MatrixXi::Constant(F.rows(),3, 0);
                int count = 0;
                for(int fidx=0; fidx< FL.rows(); ++fidx){
                    int lb = FL(fidx,0);
                    if(lb==lb_in){
                        F_l.row(count)=F.row(fidx);
                        count += 1;
                    }
                }
                F_l.conservativeResize(count,3);
                {
                    Eigen::MatrixXi I;
                    igl::remove_unreferenced(V, F_l, V_i, F_i, I);
                }
        }

        void map_vertices_to_regular_polygon(
            const Eigen::MatrixXd &V, 
            const Eigen::MatrixXi & F, 
            std::vector<node> & nodes, 
            Eigen::VectorXi & bnd,
            Eigen::MatrixXd & bnd_uv,
            std::vector<node>& ordered_nodes){
                ordered_nodes.clear();
                igl::boundary_loop(F, bnd);
                bnd_uv = Eigen::MatrixXd::Constant(bnd.rows(),3,0);
                std::vector<int> nails;
                std::map<int, int> loop_patch_dict;
                int count = -1;
                for(int bnd_idx=0; bnd_idx<bnd.rows();++bnd_idx){
                    for(node & nd:nodes){
                        if(nd.at_same_position(V.row(bnd(bnd_idx)))){
                            ordered_nodes.push_back(nd);
                            nails.push_back(bnd_idx);
                            count +=1;
                        }
                    }
                    loop_patch_dict[bnd_idx]=count;
                }

                int node_num = ordered_nodes.size();
                for(auto item: loop_patch_dict){
                    if(item.second<0) item.second=node_num-1;
                    // fix the patch label
                }

                // arc length parametrization for each patch
                std::map<int, double> patch_arc_dict;
                double sum = 0;
                int nail_count = 0;
                int start = nails[0];
                for(int shift =0 ; shift< bnd.rows();++shift){
                    int curr_idx_raw =(start+shift);
                    int curr_idx = (start+shift) % bnd.rows();
                    int next_idx = (start+shift +1) % bnd.rows();
                    double len = (V.row(bnd(next_idx,0))-V.row(bnd(curr_idx,0))).norm();
                    patch_arc_dict[curr_idx]= sum;
                    if(loop_patch_dict[curr_idx]==loop_patch_dict[next_idx]){
                        sum += len;
                    }
                    if(loop_patch_dict[curr_idx]!=loop_patch_dict[next_idx]){
                        sum += len;
                        // get the full arc length of a patch
                        int curr_nail = nails[nail_count];
                        std::cout<< curr_nail <<"in"<<nails.size()<< std::endl;
                        std::cout<< bnd.rows() <<"boundary loops"<< std::endl;
                        for(int shuttle= curr_nail; shuttle< curr_idx_raw+1; ++shuttle){
                            double ratio = patch_arc_dict[shuttle%bnd.rows()]/sum;
                        
                            std::cout<< shuttle%bnd.rows()<<"  for"<<patch_arc_dict[shuttle%bnd.rows()] <<"/"<< sum <<"="<<ratio<<std::endl;
                            patch_arc_dict[shuttle%bnd.rows()]=ratio;
                        }
                        // quotient by the full arc length to get ratio
                        nail_count+=1;
                        sum = 0;
                    }
                }
                std::cout<< patch_arc_dict.size() <<"total"<< std::endl;
                for(auto item: patch_arc_dict){
                    std::cout << item.first <<".."<<item.second <<std::endl;
                }
                start = nails[0];
                nail_count = 0;
                for(int shift =0 ; shift< bnd.rows();++shift){
                    int curr_idx = (start+shift) % bnd.rows();
                    int next_idx = (start+shift+1) % bnd.rows();
                    double x, y;
                    double ratio = patch_arc_dict[curr_idx];
                    double curr_theta = nail_count *2 * PI/ node_num;
                    double next_theta = (nail_count+1)*2 * PI/node_num;
                    x = (1-ratio)* std::cos(curr_theta) + ratio * std::cos(next_theta);
                    y = (1-ratio) * std::sin(curr_theta) + ratio * std::sin(next_theta);
                    bnd_uv.row(curr_idx) = Eigen::RowVector3d(x,y,0);
                    if(loop_patch_dict[curr_idx]!=loop_patch_dict[next_idx]){
                        nail_count+=1;
                    }
                }
        }
}