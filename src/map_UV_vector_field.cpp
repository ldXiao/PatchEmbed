//
// Copyright (C) 2019 Lind Xiao on 4/4/19 <lx471@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "map_UV_vector_field.h"
#include <Eigen/Core>
#include <igl/local_basis.h>
#include <iostream>

bool map_vectors(
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXd& V_uv,       // Vertices of uv mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXd& R_uv,       // vector field on uv mesh
        Eigen::MatrixXd& R                 // return result vector field on 3d mesh
){
    try {
        // resize R to make enough space
        R.resize(R_uv.rows(), 3);
        // construct the local basis of UV mesh
        Eigen::MatrixXd B1_uv, B2_uv, B3_uv;
        igl::local_basis(V_uv, F, B1_uv, B2_uv, B3_uv);
        // construct the local basis of 3d mesh
        Eigen::MatrixXd B1, B2, B3;
        igl::local_basis(V, F, B1, B2, B3);

        // get the entries in different directions
        for (int i = 0; i < R_uv.rows(); ++i) {
            double x = R_uv.row(i) * B1_uv.row(i).transpose();
            double y = R_uv.row(i) * B2_uv.row(i).transpose();
            double z = R_uv.row(i) * B3_uv.row(i).transpose();

            R.row(i) = B1.row(i) * x + B2.row(i) * y + B3.row(i) * z;
        }
    }
    catch(...) {
        std::cerr << "Please check that the model imported indeed has a boundary" << std::endl;
        return false;
    }

    return true;
}

//bool map_3D_vec_UV(
//        const Eigen::MatrixXd& V,          // Vertices of the mesh
//        const Eigen::MatrixXd& V_uv,       // Vertices of uv mesh
//        const Eigen::MatrixXi& F,          // Faces
//        const Eigen::MatrixXd& R,       // vector field on 3d mesh
//        Eigen::MatrixXd& R_uv               // vector field on uv mesh
//){
//    try {
//        // resize R_uv to make enough space
//        R_uv.resize(R.rows(), 3);
//        // construct the local basis of UV mesh
//        Eigen::MatrixXd B1_uv, B2_uv, B3_uv;
//        igl::local_basis(V_uv, F, B1_uv, B2_uv, B3_uv);
//        // construct the local basis of 3d mesh
//        Eigen::MatrixXd B1, B2, B3;
//        igl::local_basis(V, F, B1, B2, B3);
//
//        // get the entries in different directions
//        for (int i = 0; i < R_uv.rows(); ++i) {
//            double x = R_uv.row(i) * B1_uv.row(i).transpose();
//            double y = R_uv.row(i) * B2_uv.row(i).transpose();
//            double z = R_uv.row(i) * B3_uv.row(i).transpose();
//
//            R_uv.row(i) = B1.row(i) * x + B2.row(i) * y + B3.row(i) * z;
//        }
//    }
//    catch(...) {
//        std::cerr << "Please check that the model imported indeed has a boundary" << std::endl;
//        return false;
//    }
//
//    return true;
//}
//




