//
// Created by lind on 6/12/19.
//

#ifndef OTMAPPING_SKELETON_H
#define OTMAPPING_SKELETON_H
#include <Eigen/Core>
#include <Eigen/Sparse>

class Skeleton {
public:
    Eigen::MatrixXd Nodes;
    Eigen::MatrixXi Edges;
    Eigen::SparseMatrix<double> ElasticityTensor;
    void set_initial(
            const Eigen::MatrixXd & Sample,
            const Eigen::MatrixXi & MeshFaces,
            const Eigen::VectorXi & SampleSourceIndices);


};


#endif //OTMAPPING_SKELETON_H
