//
// Created by Lind Xiao on 7/25/19.
//

#ifndef OTMAPPING_CUTGRAPH_H
#define OTMAPPING_CUTGRAPH_H

#include <vector>
#include <map>
#include <Eigen/Core>
#include <tuple>
class CutGraph {
public:
    std::vector<std::tuple<double, double, double > > _vertices;
    std::vector<std::pair<std::size_t, std::size_t> > _edges;
    Eigen::MatrixXd EdgeWeights;
    Eigen::MatrixXd ProbabilityMatrix;
    Eigen::MatrixXd Vertices;
    Eigen::MatrixXi Edges;
    Eigen::MatrixXi Labels;
    int sample_num;
    int label_num;
    double Lambda;
    void _set_Vertices(const Eigen::MatrixXd &);
    void _set_Edges_from_KNN(int k);
    void _set_EdgeWeights(double lambda);
    void _set_ProbabilityMatrix(int label_num,Eigen::MatrixXi observed_labels);
    void initialize(const Eigen::MatrixXd&S, int knn_valance,int LabelNum, double lambda, const Eigen::MatrixXi & labels){
        this->_set_Vertices(S);
        this->_set_Edges_from_KNN(knn_valance);
        this->Lambda = lambda;
        this->Labels = labels;
        this->_set_EdgeWeights(Lambda);
        this->sample_num = S.rows();
        this->label_num = LabelNum;
    }
};

// helper function
void generate_sample_color(
        const Eigen::MatrixXd &FC,
        const Eigen::MatrixXi &source_index,
        const int & sample_num,
        Eigen::MatrixXd & SC);

void generate_sample_label(
        const Eigen::MatrixXi & FL,
        const Eigen::MatrixXi & source_index,
        const int & sample_num,
        Eigen::MatrixXi & SL
);

void NN_sample_label_transport(
        const Eigen::MatrixXd &S0,
        const Eigen::MatrixXd &S1,
        const Eigen::MatrixXi &SL0,
        Eigen::MatrixXi &SL1);

void construct_face_sample_dictionary(
        const  Eigen::MatrixXi & I1,
        std::map<int, std::vector<int> > & dict);


void NN_sample_label_vote_face_label(
        const int label_num,
        const Eigen::MatrixXi & I1,
        const Eigen::MatrixXi & SL1,
        const Eigen::MatrixXi & F1,
        Eigen::MatrixXi & FL1,
        Eigen::MatrixXd & probability_matrix
);

#endif //OTMAPPING_CUTGRAPH_H
