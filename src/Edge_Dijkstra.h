#ifndef BCCLEAN_EDGE_DIJKSTRA_H
#define BCCLEAN_EDGE_DIJKSTRA_H
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "edge.h"
namespace bcclean {
namespace Trace{
    struct comparator{
        bool operator()(const std::pair<double, int>  &lhs, const std::pair<double, int>& rhs ){
            return lhs.first > rhs.first;
        }
    };
    void setWeight(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        const Eigen::MatrixXd & V_bad,
        const Eigen::MatrixXi & F_bad,
        const edge & edg,
        Eigen::SparseMatrix<double> & Weights
    );

    void Edge_Dijkstra(
        const std::vector<std::vector<int> > & VV,
        const int source,
        const int target,
        const Eigen::SparseMatrix<double> & Weights,
        std::vector<int> & path
    );

    void setWeight1(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        const std::vector<Eigen::RowVector3d> & V_bad,
        const edge & edg,
        Eigen::SparseMatrix<double> & Weights
    );
}
}
#endif
