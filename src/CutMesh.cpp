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
#include "Sinkhorn.hpp"
using namespace OTMapping;
using MeshPair = std::tuple<Eigen::MatrixXd, Eigen::MatrixXi>;
using SamplePair = std::tuple<Eigen::MatrixXd, Eigen::VectorXi>;
using SampleTriplet = std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>;
SamplePair two_half_spaces_sample_select(
        const Eigen::MatrixXd & S,
        const Eigen::Vector3d &p0,
        const Eigen::Vector3d &n0,
        const Eigen::Vector3d &p1,
        const Eigen::Vector3d &n1
        ){
    Eigen::VectorXi Indices(S.rows());
    Eigen::RowVector3d P0= p0.transpose();
    Eigen::RowVector3d P1= p1.transpose();
    Eigen::MatrixXd S_sub(S.rows(),3);
    int count = 0;
    for(unsigned int i=0; i < S.rows(); ++i){
        if((S.row(i)-P0).dot(n0) < 0 and (S.row(i)-P1).dot(n1) < 0){
            Indices(count) = i;
            S_sub.row(count) = S.row(i);
            count+=1;
        }
    }
    Indices.conservativeResize(count);
    S_sub.conservativeResize(count,3);
    return std::make_tuple(S_sub, Indices);
}

MeshPair two_half_spaces_intersect(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::Vector3d &p0,
        const Eigen::Vector3d &n0,
        const Eigen::Vector3d &p1,
        const Eigen::Vector3d &n1){
    // Take a mesh and take boolean operation with two halfspaces return the resulting mesh Only the surface;
    //input: V Vertices
    //input: F Faces
    //input: p0,p1 points
    //input: n0, n1 normal vectors to describe the halfspace
    // retun std::tuple<Eigen::MatrixXd VD, Eigen::MatrixXi FD>;
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
    //update the viewer data with CutMesh data and options
    std::map<int, Eigen::RowVector3d> component_colors;
    component_colors.emplace(1, Eigen::RowVector3d(1,1,1));
    component_colors.emplace(2, Eigen::RowVector3d(0,1,0));
    component_colors.emplace(3, Eigen::RowVector3d(1,1,0));
    component_colors.emplace(4, Eigen::RowVector3d(1,0,1));
    std::map<int, Eigen::RowVector3d> component_sample_colors;
    component_sample_colors.emplace(1, Eigen::RowVector3d(1,0,0));
    component_sample_colors.emplace(2, Eigen::RowVector3d(0,1,0));
    component_sample_colors.emplace(3, Eigen::RowVector3d(1,1,0));
    component_sample_colors.emplace(4, Eigen::RowVector3d(1,0,1));

    Eigen::MatrixXd SamplePerturbColor(this->SampleNum, 3);
    for(unsigned int i =0; i < this->ComponentsSample.size(); ++i) {
        for(unsigned int j=0; j < this->ComponentsSample[i]->rows();++j){
            SamplePerturbColor.row(this->ComponentsSampleIndices[i]->coeff(j,0))
            = component_sample_colors[i+1];
        }
    }
    switch(options) {
        case 'i':
            viewer.init();
            viewer.data().point_size = 11;
            for(unsigned int i =0; i < this->ComponentsVertices.size(); ++i) {
                viewer.data().point_size = 11;
                viewer.append_mesh();
            }
            break;
        case '1':
//            viewer.data().clear();
            for(auto &data : viewer.data_list){
                data.clear();
            }
            std::cout << "called 1" <<std::endl;
            viewer.selected_data_index = 0;
            viewer.data().add_points(this->SampleInitial, Eigen::RowVector3d(0, 0, 1));
            viewer.data().set_mesh(this->Vertices, this->Faces);
            break;
        case '2':
            for(auto &data : viewer.data_list){
                data.clear();
            }
            viewer.selected_data_index = 0;
            viewer.data().add_points(this->SampleInitial, Eigen::RowVector3d(0, 0, 1));
            for(unsigned int i =0; i < this->ComponentsVertices.size(); ++i) {
                viewer.data_list[i+1].set_mesh(*(this->ComponentsVertices[i]), *(this->ComponentsFaces[i]));
                viewer.data_list[i+1].add_points(*(this->ComponentsSample[i]), Eigen::RowVector3d(1, 0, 0));
            }
            for(auto &data : viewer.data_list){
                data.set_colors(component_colors[data.id]);
                data.point_size = 11;
            }
            break;
        case '3':
            for(auto &data : viewer.data_list){
                data.clear();
            }
            viewer.selected_data_index = 0;
            viewer.data().add_points(this->SampleInitial, Eigen::RowVector3d(0, 0, 1));
            viewer.data().add_points(this->SamplePerturb, SamplePerturbColor);
            for(auto &data : viewer.data_list){
//                data.set_colors(component_colors[data.id]);
                data.point_size = 11;
                data.set_colors(component_colors[data.id]);
            }
            break;
        case '4':
            for(auto &data : viewer.data_list){
                data.clear();
                data.set_colors(Eigen::RowVector3d(1,1,1));
            }
            for(unsigned int i=0; i < this->ComponentsVertices.size() ;++i){
                Eigen::MatrixXd shift_back_vertices = *(this->ComponentsVertices[i]);
                for(unsigned int j=0; j < shift_back_vertices.rows();++j){
                    shift_back_vertices.row(j) -= *(this->Shifts[i]);
                }
                viewer.data_list[i+1].set_mesh(shift_back_vertices,*(this->ComponentsFaces[i]));
                viewer.data_list[i+1].set_colors(Eigen::RowVector3d(1,1,1));
            }
            viewer.selected_data_index = 0;
            viewer.data().point_size= 11;
            this->to_nearest();
            viewer.data().add_points(this->SampleInitial, this->TransportPlan*SamplePerturbColor);
            std::cout << "the color error for nearst neighbor ="
            << color_error(this->TransportPlan*SamplePerturbColor,SamplePerturbColor)<< std::endl;
            break;
        case '5':
            for(auto &data : viewer.data_list){
                data.clear();
                data.set_colors(Eigen::RowVector3d(1,1,1));
            }
            for(unsigned int i=0; i < this->ComponentsVertices.size() ;++i){
                Eigen::MatrixXd shift_back_vertices = *(this->ComponentsVertices[i]);
                for(unsigned int j=0; j < shift_back_vertices.rows();++j){
                    shift_back_vertices.row(j) -= *(this->Shifts[i]);
                }
                viewer.data_list[i+1].set_mesh(shift_back_vertices,*(this->ComponentsFaces[i]));
                viewer.data_list[i+1].set_colors(Eigen::RowVector3d(1,1,1));
            }
            this->compute_CostMatrix(this->SamplePerturb,this->SampleInitial,'a');
            this->Sinkhorn(1e-15,1e-15, 10);
            viewer.selected_data_index = 0;
            viewer.data().point_size= 11;
            viewer.data().add_points(this->SampleInitial, this->TransportPlan *SamplePerturbColor);
            std::cout << "the color error for L2-norm cost ="
                      << color_error(this->TransportPlan*SamplePerturbColor,SamplePerturbColor)<< std::endl;
    }
}

void CutMesh::set_initial(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const int n,
                          std::function<double(Eigen::Vector3d)> func) {
    //Set the initial mesh and samples
    this->Vertices = V;

    this->Faces = F;
    {
        Eigen::MatrixXi LocalI;
        Eigen::SparseMatrix<double> B;
        igl::random_points_on_mesh(n, V, F, B, LocalI);
        this->SampleInitial = B * V;
    }
    this->SampleNum = SampleInitial.rows();
    this->SamplePerturb = Eigen::MatrixXd::Zero(this->SampleNum, 3);
    this->SampleVals.resize(this->SampleNum);
    for (int i = 0; i < this->SampleNum; ++i) {
        SampleVals[i] = func(SampleInitial.row(i));
    }
}

void CutMesh::cut_with(const Eigen::Vector3d &p0,
        const Eigen::Vector3d &n0,
        const Eigen::Vector3d &p1,
        const Eigen::Vector3d &n1) {
    // slice with two planes to get 4 pieces of mesh
    int signs[4][2] = {{1,1},{-1,1},{1,-1},{-1,-1}};
    int sum = 0;
    for(unsigned int i=0; i < 4; ++i) {
        MeshPair PA = two_half_spaces_intersect(this->Vertices, this->Faces, p0, signs[i][0]*n0, p1, signs[i][1]*n1);
        SamplePair SA = two_half_spaces_sample_select(this->SampleInitial, p0, signs[i][0]*n0, p1, signs[i][1]*n1);
//        std::cout << std::get<0>(PA) << std::endl;
        std::unique_ptr<Eigen::MatrixXd> ptrVA(new Eigen::MatrixXd()), ptrSA(new Eigen::MatrixXd());
        *ptrVA = std::get<0>(PA);
        *ptrSA = std::get<0>(SA);
        std::unique_ptr<Eigen::MatrixXi> ptrFA(new Eigen::MatrixXi());
        std::unique_ptr<Eigen::VectorXi> ptrSIA(new Eigen::VectorXi());
        *ptrFA = std::get<1>(PA);
        *ptrSIA = std::get<1>(SA);
        sum += ptrSA->rows();
        this->ComponentsVertices.push_back(std::move(ptrVA));
        this->ComponentsFaces.push_back(std::move(ptrFA));
        this->ComponentsSample.push_back(std::move(ptrSA));
        this->ComponentsSampleIndices.push_back(std::move(ptrSIA));
    }
    try{
        if(sum!= this->SampleNum){
            throw sum;
        }
    }
    catch (int e){
        std::cout << "Exception accured, divided smaple num " <<
        sum  << "!=" << "totol sample num" << this->SampleNum << std::endl;
    }
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
        std::unique_ptr<Eigen::RowVector3d > shift_ptr(new Eigen::RowVector3d());
        *shift_ptr = shift;
        this->Shifts.push_back(std::move(shift_ptr));
        for(unsigned int k =0; k < this->ComponentsVertices[i]->rows();k++){
            this->ComponentsVertices[i]->row(k)+= shift;
        }
        for(unsigned int k =0; k < this->ComponentsSample[i]->rows();k++){
            this->ComponentsSample[i]->row(k)+= shift;
            this->SamplePerturb.row(this->ComponentsSampleIndices[i]->coeff(k,0))=this->ComponentsSample[i]->row(k);
        }
    }
//    std::cout <<this->SamplePerturb-this->SampleInitial << std::endl;
}

double color_error(Eigen::MatrixXd C0, Eigen::MatrixXd C1){
    return ((C0-C1).rowwise().norm()).sum();
}

Eigen::MatrixXd CutMesh::to_nearest(){
    this->TransportPlan = Eigen::MatrixXd::Zero(this->SampleNum,this->SampleNum);
    Eigen::MatrixXd SampleNearest(this->SampleNum, 3);
    for(unsigned int i=0; i< this->SampleNum; ++i){
        int min_idx = 0;
        double min_dist = (this->SamplePerturb.row(min_idx)-this->SampleInitial.row(i)).norm();
        for(unsigned int j=1; j< this->SampleNum; ++j){
            double cur_dist=(this->SamplePerturb.row(j)-this->SampleInitial.row(i)).norm();
            if(cur_dist < min_dist){
                min_dist = cur_dist;
                min_idx = j;

            }
        }
        this->TransportPlan(i, min_idx) = 1;
        SampleNearest.row(i) = this->SamplePerturb.row(min_idx);
    }
    return SampleNearest;
}

void CutMesh::compute_CostMatrix(Eigen::MatrixXd source, Eigen::MatrixXd target, char option) {
    try {
        if (source.cols() != 3 or target.cols() != 3) {
            throw "Unimplemented Error";
        }
        if (source.rows() != target.rows()){
            throw "Unequal Error";
        }
        int n = this->SampleNum;
        this->CostMatrix.resize(n, n);
        switch(option) {
            case 'a':
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        this->CostMatrix(i, j) = (this->SampleInitial.row(i) - this->SamplePerturb.row(j)).lpNorm<Eigen::Infinity>();
                    }
                }
            case 'b':
                for (int i = 0; i < target.rows(); ++i) {
                    for (int j = 0; j < source.rows(); ++j) {
                        this->CostMatrix(i, j) = (this->SampleInitial.row(i) - this->SamplePerturb.row(j)).lpNorm<Eigen::Infinity>();
                    }
                }

        }
    }
    catch(std::string a){
        std::cout << a << std::endl;
    }
}

void CutMesh::Sinkhorn(double eps, double stop_thresh, int max_iters){
    /*
    Compute the Sinkhorn divergence between two sum of dirac delta distributions, U, and V.
            This implementation is numerically stable with float32.
    :param eps: The reciprocal of the sinkhorn regularization parameter
    :param max_iters: The maximum number of Sinkhorn iterations
    :param stop_thresh: Stop if the change in iterates is below this value
    :return:
    */
    try {
        if (this->CostMatrix.rows()< this->SampleInitial.rows()) {
            throw "CostMatrix not initialized";
        }
//        if (this->TransportPlan.rows() < this->SampleInitial.rows()){
//            throw "TransportPlan not initialized";
//        }
    }
    catch(std::string a){
    }
    int n = this->SampleNum;
    this->TransportPlan = sinkhorn(
            Eigen::VectorXd::Constant(n,1),
            Eigen::VectorXd::Constant(n,1),
            this->CostMatrix, eps, max_iters, stop_thresh);
    std::cout<<" Totoal Cost for this TransportPlan=" << (this->TransportPlan.array() * this->CostMatrix.array()).sum()
    << '\n'<<" Totoal Cost for identity="
    << (Eigen::MatrixXd::Identity(n,n).array() * this->CostMatrix.array()).sum() <<std::endl;
}
