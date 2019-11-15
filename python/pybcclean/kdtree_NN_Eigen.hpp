//
// Created by Lind Xiao on 8/13/19.
//

#ifndef OTMAPPING_KDTREE_NN_EIGEN_H
#define OTMAPPING_KDTREE_NN_EIGEN_H

#include <nanoflann.hpp>
#include <cstdlib>
#include <ctime>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Core>
using namespace Eigen;
using namespace std;
using namespace nanoflann;

namespace bcclean {
    template<typename num_t>
    using kd_tree_Eigen=KDTreeEigenMatrixAdaptor <Eigen::Matrix<num_t, Dynamic, Dynamic>>;


    template<typename num_t>
    int kd_tree_NN_Eigen(const kd_tree_Eigen<num_t> &mat_index, const Eigen::Matrix<num_t, 1, 3> & query_pt) {
        // do a knn search
        const size_t num_results = 1;
        vector <size_t> ret_indexes(num_results);
        vector <num_t> out_dists_sqr(num_results);

        nanoflann::KNNResultSet <num_t> resultSet(num_results);

        resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
        mat_index.index->findNeighbors(resultSet, &query_pt[0],
                                       nanoflann::SearchParams(10));
        return ret_indexes[0];
    }
}
#endif //OTMAPPING_KDTREE_NN_EIGEN_H
