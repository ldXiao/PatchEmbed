//
// Created by Lind Xiao on 6/6/19.
//

#include "CutMesh.h"
#include <igl/random_points_on_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/cgal/intersect_with_half_space.h>
#include <igl/random_points_on_mesh.h>
#include <igl/uniformly_sample_two_manifold.h>
#include <Eigen/Core>
#include <memory.h>
#include <utility>
#include <functional>
#include <tuple>
#include <random>
#include <map>
using namespace OTMapping;
using MeshPair = std::tuple<Eigen::MatrixXd, Eigen::MatrixXi>;
using SampleTriplet = std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>;
MeshPair two_half_spaces_intersect(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::Vector3d &p0,
        const Eigen::Vector3d &n0,
        const Eigen::Vector3d &p1,
        const Eigen::Vector3d &n1){
    Eigen::MatrixXd V0, V1;
    Eigen::MatrixXi F0, F1;
    Eigen::VectorXi J0, J1;
    /* intersect twice and try to retrace where each fase comes from*/
    igl::copyleft::cgal::intersect_with_half_space(V, F, p0, n0, V0, F0,J0);
    igl::copyleft::cgal::intersect_with_half_space(V0, F0, p1, n1, V1, F1,J1);
    Eigen::MatrixXd VD;
    Eigen::MatrixXi FD(F1.rows(),3);
    int count = 0;
    for(int f =0; f < F1.rows(); ++f){
        if(J1(f) < F0.rows() and J0(J1(f))<F.rows()){
            FD.row(count)=F1.row(f);
            count+=1;
        }
    }
    FD.conservativeResize(count, 3);
    {
        Eigen::VectorXi Local;
        igl::remove_unreferenced(V1, Eigen::MatrixXi(FD), VD, FD, Local);
    }
    return std::make_tuple(VD,FD);
}

void CutMesh::plot_CutMesh(igl::opengl::glfw::Viewer &viewer, unsigned char options) {
    std::map<int, Eigen::RowVector3d> colors;
    viewer.data().point_size = 5;
    switch(options) {
        case '1':
            viewer.data().clear();
            viewer.data().add_points(this->SampleInitial, Eigen::RowVector3d(0, 0, 1));
            viewer.data().set_mesh(this->Vertices, this->Faces);
        case '2':
            viewer.data().clear();
            viewer.data().add_points(this->SampleInitial, Eigen::RowVector3d(0, 0, 1));
            for(unsigned int i =0; i < this->ComponentsVertices.size(); ++i) {
                viewer.append_mesh();
                viewer.data().set_mesh(*(this->ComponentsVertices[i]), *(this->ComponentsFaces[i]));
                colors.emplace(viewer.data().id, 0.5*Eigen::RowVector3d::Random().array() + 0.5);
            }
            for(auto &data : viewer.data_list){
                data.set_colors(colors[data.id]);
            }
    }
}

void CutMesh::set_initial(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const int n,
                          std::function<double(Eigen::Vector3d)> func) {
    this->Vertices = V;

    this->Faces = F;
    {
        Eigen::MatrixXi LocalI;
        Eigen::SparseMatrix<double> B;
        igl::random_points_on_mesh(n, V, F, B, LocalI);
        this->SampleInitial = B * V;
    }
    this->SampleNum = SampleInitial.rows();
    this->SampleVals.resize(this->SampleNum);
    for (int i = 0; i < this->SampleNum; ++i) {
        SampleVals[i] = func(SampleInitial.row(i));
    }
}

void CutMesh::cut_with(const Eigen::Vector3d &p0,
        const Eigen::Vector3d &n0,
        const Eigen::Vector3d &p1,
        const Eigen::Vector3d &n1) {
    Eigen::MatrixXd VA, VB, VC, VD;
    Eigen::MatrixXi FA, FB, FC, FD;
    MeshPair PA = two_half_spaces_intersect(this->Vertices, this->Faces, p0, n0, p1, n1);
    MeshPair PB = two_half_spaces_intersect(this->Vertices, this->Faces, p0,-n0, p1, n1);
    MeshPair PC = two_half_spaces_intersect(this->Vertices, this->Faces, p0, n0, p1, -n1);
    MeshPair PD = two_half_spaces_intersect(this->Vertices, this->Faces, p0, -n0, p1, -n1);
    std::cout << std::get<0>(PA) <<std::endl;
    std::unique_ptr<Eigen::MatrixXd >
            ptrVA(new Eigen::MatrixXd()),
            ptrVB(new Eigen::MatrixXd()),
            ptrVC(new Eigen::MatrixXd()),
            ptrVD(new Eigen::MatrixXd());
    *ptrVA = std::get<0>(PA);
    *ptrVB = std::get<0>(PB);
    *ptrVC = std::get<0>(PC);
    *ptrVD = std::get<0>(PD);

    std::unique_ptr<Eigen::MatrixXi >
            ptrFA(new Eigen::MatrixXi()), ptrFB(new Eigen::MatrixXi()),
            ptrFC(new Eigen::MatrixXi()), ptrFD(new Eigen::MatrixXi());
    *ptrFA = std::get<1>(PA);
    *ptrFB = std::get<1>(PB);
    *ptrFC = std::get<1>(PC);
    *ptrFD = std::get<1>(PD);
    this->ComponentsVertices.push_back(std::move(ptrVA));
    this->ComponentsVertices.push_back(std::move(ptrVB));
    this->ComponentsVertices.push_back(std::move(ptrVC));
    this->ComponentsVertices.push_back(std::move(ptrVD));

    this->ComponentsFaces.push_back(std::move(ptrFA));
    this->ComponentsFaces.push_back(std::move(ptrFB));
    this->ComponentsFaces.push_back(std::move(ptrFC));
    this->ComponentsFaces.push_back(std::move(ptrFD));
}

void CutMesh::perturb(const int seed, double max_shift, double min_shift) {
    unsigned int num_components = this->ComponentsVertices.size();
    std::mt19937 gen;
    gen.seed(seed);
    std::uniform_real_distribution<> dis(max_shift, min_shift);
    std::map<int, Eigen::Vector3i> signs;
    signs.emplace(0, Eigen::Vector3i(-1,-1,1));
    signs.emplace(1, Eigen::Vector3i(-1,1,1));
    signs.emplace(2, Eigen::Vector3i(1,1,-1));
    signs.emplace(3, Eigen::Vector3i(1,-1,1));
    for(unsigned int i = 0; i < num_components; ++i){
        Eigen::RowVector3d shift;
        for(unsigned int j=0; j < 3; ++j){
            shift(j) = signs[i](j)*dis(gen);
        }
        for(unsigned int k =0; k < this->ComponentsVertices[i]->rows();k++){
            this->ComponentsVertices[i]->row(k)+= shift;
        }
    }
}