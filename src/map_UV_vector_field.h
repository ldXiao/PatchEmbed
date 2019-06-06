//
// Copyright (C) 2019 Lind Xiao on 4/4/19 <lx471@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef NROSY_FIELD_MAP_UV_VECTOR_FIELD_H
#define NROSY_FIELD_MAP_UV_VECTOR_FIELD_H
#include <Eigen/Core>


bool map_vectors(
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXd& V_uv,       // Vertices of uv mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXd& R_uv,       // vector field on uv mesh
        Eigen::MatrixXd& R                 // return result vector field on 3d mesh
);


bool map_3D_vec_UV(
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXd& V_uv,       // Vertices of uv mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXd& R,       // vector field on 3d mesh
        Eigen::MatrixXd& R_uv               // vector field on uv mesh
);

#endif //NROSY_FIELD_MAP_UV_VECTOR_FIELD_H
