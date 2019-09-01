//
// Created by Lind Xiao on 7/25/19.
//

#ifndef OTMAPPING_CUTGRAPH_H
#define OTMAPPING_CUTGRAPH_H

#include <vector>
#include <map>
#include <Eigen/Core>
#include <tuple>
namespace  bcclean {
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

    void LM_intersection_label_transport(
            const Eigen::MatrixXd &V0,
            const Eigen::MatrixXi & F0,
            const Eigen::MatrixXi & FL0,
            const Eigen::MatrixXd & V1,
            const Eigen::MatrixXi & F1,
            Eigen::MatrixXi & VL1);

    void construct_face_sample_dictionary(
            const Eigen::MatrixXi &I1,
            std::map<int, std::vector<int> > &dict);

    void construct_face_nearby_sample_dictionary(
            const Eigen::MatrixXi &F,
            const std::map<int, std::vector<int> > &ETT,
            std::map<int, std::vector<int> > &dict);

    void Barycenter_intersection_label_transport(
        const Eigen::MatrixXd &V0,
        const Eigen::MatrixXi & F0,
        const Eigen::MatrixXi & FL0,
        const Eigen::MatrixXd & V1,
        const Eigen::MatrixXi & F1,
        Eigen::MatrixXi & FL1);

    void NN_sample_label_vote_face_label(
            const int label_num,
            const Eigen::MatrixXi &I1,
            const Eigen::MatrixXi &SL1,
            const Eigen::MatrixXi &F1,
            Eigen::MatrixXi &FL1,
            Eigen::MatrixXd &probability_matrix
    );

    void vertex_label_vote_face_label(
            const int label_num,
            const Eigen::MatrixXi &VL,
            const Eigen::MatrixXi &F,
            Eigen::MatrixXi &FL,
            Eigen::MatrixXd &prob_mat
    );

    void set_EdgeWeight(const double &lambda, const Eigen::MatrixXi &SL, const Eigen::MatrixXi &E, Eigen::MatrixXd &EW);

    void set_Face_Edges(const Eigen::MatrixXi &F, Eigen::MatrixXi &E);

    void normalize_mesh(Eigen::MatrixXd & V_bad, Eigen::MatrixXd & V_good);

    void build_patch_dict(const Eigen::MatrixXi &FL, std::map<int, std::vector<int> > & patch_dict);

}

#endif //OTMAPPING_CUTGRAPH_H
