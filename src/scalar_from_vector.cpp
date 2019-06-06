//
// Copyright (C) 2019 Lind Xiao on 4/4/19 <lx471@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include "scalar_from_vector.h"
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <iostream>
using namespace std;
typedef Eigen::SparseMatrix<double> SparseMatrixXd;

std::vector<Eigen::Triplet<double>> to_triplets(Eigen::SparseMatrix<double> & M){
    std::vector<Eigen::Triplet<double>> v;
    for(int i = 0; i < M.outerSize(); i++)
        for(typename Eigen::SparseMatrix<double>::InnerIterator it(M,i); it; ++it)
            v.emplace_back(it.row(),it.col(),it.value());
    return v;
}


Eigen::MatrixXd scalar_from_vector(
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXi& TT,         // Adjacency triangle-triangle
        const Eigen::MatrixXd& R           // Constructed vector field
){
    // minimize A_f ((G s)_f - u_f)^2 = (s^t * G^t - u^t) diag(A) (G s -u)
    // = s^t G^t A G s -u^t AG s- s^t G^t A u + u^t A u
    // = s^t K s + s^t b  +c
    // K = G^t A G, b= -2 G^t A u , c =u^t A u
    // solve 2 K s + b =0 , G^t A G s = G^t A u
    SparseMatrixXd G(F.rows() * 3, V.rows());
    igl::grad(V, F, G);
    // After looking around in the documents, we know G is if formated as (f_x, f_y, f_z.....)
    SparseMatrixXd A(F.rows() * 3, F.rows() * 3);
    std::vector<Eigen::Triplet<double> > t_A;
    Eigen::MatrixXd areas(F.rows(),1);
    igl::doublearea(V, F, areas);
    for(int i = 0; i < F.rows(); ++i){
        double area = areas(i,0);
        t_A.push_back(Eigen::Triplet<double>(i , i, area));
        t_A.push_back(Eigen::Triplet<double>(i + F.rows(), i + F.rows(), area));
        t_A.push_back(Eigen::Triplet<double>(i + 2 * F.rows(), i + F.rows() * 2, area));
    }
    A.setFromTriplets(t_A.begin(), t_A.end());
    // calculate the extra constraints
    SparseMatrixXd GAG = G.adjoint() * A * G;
    SparseMatrixXd K_extra(1, GAG.rows() *3);
    K_extra.coeffRef(0,1) = 1;

    // convert the GAG back into triplet
    std::vector<Eigen::Triplet<double> > t_K = to_triplets(GAG);
    SparseMatrixXd K(GAG.rows() + 1, GAG.rows());
    t_K.push_back(Eigen::Triplet<double>(K.rows()-1, 0, 1.0 ));
    K.setFromTriplets(t_K.begin(), t_K.end());

    // import R with new constraint
    Eigen::MatrixXd u(3 * F.rows(), 1);
    for(int i =0; i < F.rows(); ++i){
        u(i,0)=R(i,0);
        u(i + F.rows() , 0) = R(i,1);
        u(i + F.rows()*2, 0) = R(i, 2);
    }

    Eigen::MatrixXd GAu = G.adjoint() * A * u;
    Eigen::MatrixXd b(K.rows(), 1);
    b << GAu, 1;
    Eigen::SimplicialLDLT<SparseMatrixXd> solver;
    solver.compute(K.adjoint() * K);
    assert(solver.info()==Eigen::Success);
    Eigen::MatrixXd s = solver.solve(K.adjoint()*b);
    assert(solver.info()==Eigen::Success);
    s(0,0)=1;
//    Eigen::SparseLU<SparseMatrixXd> solver2;
//    solver2.compute(GAG);
//    assert(solver2.info()==Eigen::Success);
//    Eigen::MatrixXd ss = solver2.solve(GAu);
//    assert(solver2.info()==Eigen::Success);
//    cout << ss - s << endl;
    return s;
}


Eigen::MatrixXd gradient_scalar(
        const Eigen::MatrixXd& V,          // Vertices of the mesh
        const Eigen::MatrixXi& F,          // Faces
        const Eigen::MatrixXd& s
){
    SparseMatrixXd G(F.rows() * 3, V.rows());
    igl::grad(V, F, G);
//    cout <<(G*s)<<endl;
    Eigen::MatrixXd Gs = Eigen::Map<const Eigen::MatrixXd>((G*s).eval().data(),F.rows(),3);
//    cout << (Gs)<<endl;
    return Gs;
}
