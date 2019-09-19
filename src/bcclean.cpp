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
#include <igl/upsample.h>
#include <random>
#include <vector>
#include <map>
#include <cmath>



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

    void refine_proj_vote(
        const Eigen::MatrixXd &V0,
            const Eigen::MatrixXi & F0,
            const Eigen::MatrixXi & FL0,
            const Eigen::MatrixXd & V1,
            const Eigen::MatrixXi & F1,
            const int label_num,
            const int subdiv,
            Eigen::MatrixXi & FL1,
            Eigen::MatrixXd & prob_mat
    ){
        try{
            Eigen::MatrixXd NV;
            Eigen::MatrixXi NF, NFL, NVL;
            // FI is used to trace
            igl::upsample(V1, F1, NV, NF, subdiv);
            //if subdiv==1 fi face in F1 will be corresponds to 4 * fi +0, 4 * fi +1,4 * fi +2, 4 * fi +3 faces in NF
            // reversely for general nf in NF, it corresponds to face nf % pow(4, subdiv) in F1
            size_t modu = std::pow(4, subdiv);
            LM_intersection_label_transport(V0, F0, FL0, NV, NF, NVL);
            {
                Eigen::MatrixXd prob_mat_temp;
                vertex_label_vote_face_label(label_num, NVL, NF, NFL, prob_mat_temp);
            }
            std::vector< std::map<int, int> > F_label_dict;
            for(int nfidx=0; nfidx < NF.rows(); ++ nfidx){
                int fidx = (int) nfidx / modu;
                if(F_label_dict.size()< fidx+1){
                    F_label_dict.push_back(std::map<int, int>());
                } else{
                    int nfl = NFL(nfidx,0);
                    auto it =F_label_dict[fidx].find(nfl);
                    if(it == F_label_dict[fidx].end()){
                        F_label_dict[fidx][nfl]=1;
                    } else{
                        F_label_dict[fidx][nfl]+=1;
                    }
                }
            }
            int fid =0;
            prob_mat = Eigen::MatrixXd::Constant(F1.rows(), label_num, 0.05 / label_num);
            for(auto mp: F_label_dict){
                int target_lb = -1;
                int target_lb_count = 0;
                for(auto item : mp){
                    if(item.second >target_lb_count){
                        target_lb = item.first;
                        target_lb_count= item.second;
                    }
                    prob_mat(fid, item.first) = 0.95 * (double(item.second)/double(modu));
                    // std::cout << "prob_mat"<<0.95 * (item.second / modu)<< std::endl;
                }
                fid+=1;
            }
        } catch(...){
            std::cout << "upsample fails" <<std::endl;
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
}