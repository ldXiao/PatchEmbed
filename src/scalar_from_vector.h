//
// Copyright (C) 2019 Lind Xiao on 4/4/19 <lx471@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef NROSY_FIELD_SCALAR_FROM_VECTOR_H
#define NROSY_FIELD_SCALAR_FROM_VECTOR_H

#endif //NROSY_FIELD_SCALAR_FROM_VECTOR_H

#include <Eigen/Core>

Eigen::MatrixXd scalar_from_vector(
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXi& TT,         // Adjacency triangle-triangle
        const Eigen::MatrixXd& R           // Constructed vector field
        );

Eigen::MatrixXd gradient_scalar(
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXd& s
        );