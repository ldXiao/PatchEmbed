//
// Created by Lind Xiao on 7/25/19.
//

#include "CutGraph.h"
#include <Eigen/Core>
#include <iostream>
#include <map>
#include <tuple>
#include <igl/opengl/glfw/Viewer.h>
#include <random>
#include <vector>
void CutGraph::_set_Vertices(const Eigen::MatrixXd & S) {
    int num = S.rows();
    this->_vertices.resize(num);
    this->Vertices = S;
    {
        int count = 0;
        for (auto it = this->_vertices.begin(); it != this->_vertices.end(); ++it) {
            *it = std::make_tuple(S(count,0),S(count,1), S(count,2));
            count += 1;
        }
    }
}

void CutGraph::_set_Edges_from_KNN(int m) {
    int count = 0;
    int esize = 2*m * this->_vertices.size();
    this->_edges.resize(esize);
    this->Edges= Eigen::MatrixXi::Zero(esize,2);
    std::map<std::pair<int, int>, int > map;
    for(int j=0; j < this->Vertices.rows(); ++j){
        std::vector<std::pair<int, double> > j_stack;
        for (int k = 0; k < this->Vertices.rows(); ++k) {
            if(k==j) continue;
            double norm_jk = (this->Vertices.row(k) - this->Vertices.row(j)).norm();
            if(j_stack.size()==0 and k!=j){
                j_stack.push_back(std::make_pair(k, norm_jk));
            } else {
                auto it = j_stack.begin();
                while (it->second > norm_jk and it!= j_stack.end()) {
                    it++;
                    // loop to find a position to insert
                }
                if (it == j_stack.begin()) {
                    // if does not find within
                    if (j_stack.size() < m) {
                        // if the stack is not full
                        j_stack.insert(it, std::make_pair(k, norm_jk));
                    }
                } else {
                    // find a position to insert including j_stack.end()
                    j_stack.insert(it, std::make_pair(k, norm_jk));
                    if (j_stack.size() > m) {
                        j_stack.erase(j_stack.begin(), j_stack.begin() + 1);
                    }
                }
            }
        }
        // inserted the detected m-nearest-pair in to the this->Edges

        for(auto jt = j_stack.begin(); jt!= j_stack.end(); jt++){
            int k = jt->first;
            auto find0 = map.find(std::make_pair(j,k));
            auto find1 = map.find(std::make_pair(k,j));
            if( find0 == map.end() and find1 == map.end()) {
                // insert each edge in both directions
                this->_edges[count] = std::make_pair(j, k);
                this->_edges[count+1] = std::make_pair(k, j);
                Eigen::RowVector2i temp;
                temp << j,k;
                Eigen::RowVector2i temp1;
                temp1 << k,j;
                this->Edges.row(count) = temp;
                this->Edges.row(count+1)=temp;
                count += 2;
                map[std::make_pair(j,k)]=1;
                map[std::make_pair(k,j)]=1;
            }
        }
    }
    this->_edges.resize(count);
    this->Edges.conservativeResize(count,2);
}

void CutGraph::_set_EdgeWeights(double lambda){
    this->EdgeWeights = Eigen::MatrixXd::Constant(this->Edges.rows(),1, 0);
    for(int i =0; i< this->Edges.rows();++i){
        int j = this->Edges(i,0);
        int k = this->Edges(i,1);
        if(this->Labels(j,0)==this->Labels(k,0)){
            this->EdgeWeights(i,0)= lambda;
        }
    }
}

void CutGraph::_set_ProbabilityMatrix(int label_num, Eigen::MatrixXi observed_labels) {
    this->ProbabilityMatrix = Eigen::MatrixXd::Zero(this->sample_num, this->label_num);
    std::mt19937 gen;
    gen.seed(1002);
    std::uniform_real_distribution<> dis(0, 1);
    for(int i = 0; i< this->sample_num; ++i){
        int lbl = observed_labels(i,0);
        double p = dis(gen);
        this->ProbabilityMatrix.row(i) = Eigen::MatrixXd::Constant(1,this->label_num, (0)/ this->label_num);
        this->ProbabilityMatrix(i,lbl) = 1.0;
    }

    std::cout<< "adaesfhod9aewhfcla" << this->ProbabilityMatrix;
}


void generate_sample_color(
        const Eigen::MatrixXd &FC,
        const Eigen::MatrixXi & source_index,
        const int & sample_num,
        Eigen::MatrixXd & SC) {
    SC.resize(sample_num, 3);
    for(int i =0; i < sample_num;++i){
        SC.row(i)=FC.row(source_index(i,0));
    }

}

void generate_sample_label(
        const Eigen::MatrixXi & FL,
        const Eigen::MatrixXi & source_index,
        const int & sample_num,
        Eigen::MatrixXi & SL
        ){
    SL.resize(sample_num, 3);
    for(int i =0; i < sample_num;++i){
        SL(i,0)=FL(source_index(i,0),0);
    }
}

void NN_sample_label_transport(
        const Eigen::MatrixXd &S0,
        const Eigen::MatrixXd &S1,
        const Eigen::MatrixXi &SL0,
        Eigen::MatrixXi &SL1){
    int sample_num = S1.rows();
    Eigen::VectorXi Nearest(S1.rows());
    SL1.resize(SL0.rows(),1);
    for(unsigned int i=0; i< sample_num; ++i){
        int min_idx = 0;
        double min_dist = (S0.row(min_idx)-S1.row(i)).norm();
        for(unsigned int j=0; j< sample_num; ++j){
            double cur_dist=(S0.row(j)-S1.row(i)).norm();
            if(cur_dist < min_dist){
                min_dist = cur_dist;
                min_idx = j;
            }
        }
        SL1(i,0) = SL0(min_idx,0);
    }
}

void construct_face_sample_dictionary(
        const  Eigen::MatrixXi & I1,
        std::map<int, std::vector<int> > & dict){
    for(int i =0; i< I1.rows(); ++i){
        auto find = dict.find(I1(i,0));
        if(find != dict.end()){
            dict[I1(i,0)].push_back(i);
        } else{
            dict[I1(i,0)] = std::vector<int>();
        }
    }
}
void NN_sample_label_vote_face_label(
        const int label_num,
        const Eigen::MatrixXi & I1,
        const Eigen::MatrixXi & SL1,
        const Eigen::MatrixXi & F1,
        Eigen::MatrixXi & FL1,
        Eigen::MatrixXd & probability_matrix
        ){
    std::map<int, std::vector<int> > dict;
    construct_face_sample_dictionary(I1, dict);
    FL1 = Eigen::MatrixXi::Zero(F1.rows(), 1);
    std::cout  << dict.size() << "size" << std::endl;
    probability_matrix = Eigen::MatrixXd::Constant(F1.rows(),label_num, 1e-5);
    for(int fidx=0; fidx< FL1.rows(); ++fidx){
        auto find = dict.find(fidx);
        if(find != dict.end()) {
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

            probability_matrix.row(fidx) /= size;
            FL1(fidx, 0) = max_num_label;
        } else{
            probability_matrix.row(fidx) = Eigen::MatrixXd::Constant(1,label_num, 1/label_num);
            FL1(fidx, 0) = fidx % label_num;
        }
    }
}


