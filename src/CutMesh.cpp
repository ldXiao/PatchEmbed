//
// Created by Lind Xiao on 6/6/19.
//

#include "CutMesh.h"
#include <igl/random_points_on_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/random_points_on_mesh.h>
#include <igl/uniformly_sample_two_manifold.h>
#include <Eigen/Core>
#include <functional>

void CutMesh::plot_CutMesh(igl::opengl::glfw::Viewer &viewer, std::string options){
    viewer.data().clear();
    viewer.data().set_mesh(this->Vertices, this->Faces);

    viewer.data().add_points(this->SampleInitial, Eigen::RowVector3d(0, 0, 1));
}

void CutMesh::set_initial(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, const int n, std::function<double(Eigen::Vector3d)> func) {
    this->Vertices = V;
    std::cout << this->SampleInitial<<std::endl;
    this->Faces = F;
    {
        Eigen::MatrixXi LocalI;
//        igl::random_points_on_mesh(n, V, F, this->SampleInitial, LocalI);
        igl::uniformly_sample_two_manifold(V, F, n,0.1,this->SampleInitial);
    }
    this->SampleNum = SampleInitial.rows();
    std::cout << this->SampleInitial<<std::endl;
    this->SampleVals.resize(this->SampleNum);
    for(int i =0 ; i < this->SampleNum; ++i){
        SampleVals[i]=func(SampleInitial.row(i));
    }
}