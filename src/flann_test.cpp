//
// Created by Lind Xiao on 8/13/19.
//
/*************************************************************************/

#include <nanoflann.hpp>

#include <cstdlib>
#include <ctime>
#include <iostream>

#include <Eigen/Dense>
#include "kdtree_NN_Eigen.hpp"


int main(int argc, char **argv) {
    // Randomize Seed
    srand(static_cast<unsigned int>(time(nullptr)));

    Eigen::MatrixXd mat(5,3);
    mat << 1,0,0,
            2,0,0,
            3,0,0,
            4,0,0,
            1,1,1;
    Eigen::RowVector3d v(1,1,2);
    OTMapping::kd_tree_Eigen<double> kdt(mat.cols(),std::cref(mat),10);
    kdt.index->buildIndex();
    std::cout << OTMapping::kd_tree_NN_Eigen(kdt, v)<<std::endl;

    return 0;
}