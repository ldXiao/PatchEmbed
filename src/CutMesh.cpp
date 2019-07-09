//
// Created by Lind Xiao on 6/6/19.
//

#include "CutMesh.h"
#include "Sinkhorn.hpp"
#include "Kruskal.hpp"
#include <igl/random_points_on_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/copyleft/cgal/intersect_with_half_space.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/readOBJ.h>
#include <Eigen/Core>
#include <memory.h>
#include <utility>
#include <functional>
#include <tuple>
#include <random>
#include <map>
#include <nlohmann/json.hpp>
#include "ExpressionValue.hpp"
#include "tinyexpr.h"
using namespace OTMapping;
using MeshPair = std::tuple<Eigen::MatrixXd, Eigen::MatrixXi>;
using MeshTriple = std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>;
using SamplePair = std::tuple<Eigen::MatrixXd, Eigen::VectorXi>;
using SampleTriplet = std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>;
SamplePair two_half_spaces_sample_select(
        const Eigen::MatrixXd & S,
        const Eigen::Vector3d &p0,
        const Eigen::Vector3d &n0,
        const Eigen::Vector3d &p1,
        const Eigen::Vector3d &n1){
    /* As the mesh is cutted by two plane,
     * this function would divide the generated
     * sample on the mesh and give back the subset
     * of the sample that are above both of the two planes
     * return std::pair(Sample_subset, subset_indices)
     * */
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

SamplePair components_sample_select(
        const Eigen::MatrixXd & S,
        const Eigen::VectorXi & SamplesSourceFaces,
        const Eigen::VectorXi & FacesSourceComponents,
        const int & component_idx
        ){
    Eigen::VectorXi Indices(S.rows());
    Eigen::MatrixXd S_sub(S.rows(),3);
    int count = 0;
    for(unsigned int i=0; i < S.rows(); ++i){
        if(FacesSourceComponents(SamplesSourceFaces(i))==component_idx){
            Indices(count) = i;
            S_sub.row(count) = S.row(i);
            count += 1;
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



MeshTriple cube_surface_select(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        const int & idx,
        const int & entry,
        const double &val){
    /*work for the test case of this cube.obj only
     * separate the six if the triangle's $(entry)
     * projection is contain is close to val within
     * some tolerance, select the face
     * return this subset of mesh together with
     * the face indices in the original mesh*/
    Eigen::MatrixXi FD(F.rows(),3);
    Eigen::MatrixXd VD;
    Eigen::VectorXi Indices(F.rows());
    int count = 0;
    for(int f =0; f < F.rows(); ++f){
        int a,b,c;
        a = F(f,0);
        b = F(f,1);
        c = F(f,2);
        if(std::abs(V(a,entry)-val)<0.05 and std::abs(V(b,entry)-val)<0.05 and std::abs(V(c,entry)-val)<0.05){
            FD.row(count)=F.row(f);
            Indices(count)=f;
            count+=1;
        }
    }
    FD.conservativeResize(count, 3);
    Indices.conservativeResize(count);
    {
        Eigen::VectorXi Local;
        igl::remove_unreferenced(V, Eigen::MatrixXi(FD), VD, FD, Local);
    }

    return std::make_tuple(VD, FD, Indices);
}

void CutMesh::plot_CutMesh(igl::opengl::glfw::Viewer &viewer, unsigned char options) {
    //update the viewer data with CutMesh data and options
    std::map<int, Eigen::RowVector3d> component_colors;
    component_colors.emplace(1, Eigen::RowVector3d(1,1,1));
    component_colors.emplace(2, Eigen::RowVector3d(0,1,0));
    component_colors.emplace(3, Eigen::RowVector3d(1,1,0));
    component_colors.emplace(4, Eigen::RowVector3d(1,0,1));
    component_colors.emplace(5, Eigen::RowVector3d(0,0,1));
    component_colors.emplace(6, Eigen::RowVector3d(0,1,1));
    std::map<int, Eigen::RowVector3d> component_sample_colors;
    component_sample_colors.emplace(1, Eigen::RowVector3d(1,0,0));
    component_sample_colors.emplace(2, Eigen::RowVector3d(0,1,0));
    component_sample_colors.emplace(3, Eigen::RowVector3d(1,1,1));
    component_sample_colors.emplace(4, Eigen::RowVector3d(1,0,1));
    component_sample_colors.emplace(5, Eigen::RowVector3d(0,0,1));
    component_sample_colors.emplace(6, Eigen::RowVector3d(0,1,1));

    Eigen::MatrixXd SamplePerturbColor(this->SampleNum, 3);
    for(unsigned int i =0; i < this->ComponentsSample.size(); ++i) {
        for(unsigned int j=0; j < this->ComponentsSample[i]->rows();++j){
            SamplePerturbColor.row(this->ComponentsSampleIndices[i]->coeff(j,0))
            = component_sample_colors[i+1];
        }
    }
    int point_size = this->point_size;
    switch(options) {
        case 'i':
            viewer.init();
            viewer.data().point_size = point_size;
            for(unsigned int i =0; i < this->ComponentsVertices.size(); ++i) {
                viewer.data().point_size = point_size;
                viewer.append_mesh();
            }
            break;
        case '1':
//            viewer.data().clear();
            for(auto &data : viewer.data_list){
                data.clear();
            }
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
                data.point_size = point_size;
            }
            break;
        case '3':
            for(auto &data : viewer.data_list){
                data.clear();
            }
            viewer.selected_data_index = 0;
            viewer.data().add_points(this->SampleInitial, Eigen::RowVector3d(0, 0, 1));
            viewer.data().add_points(this->SamplePerturb, SamplePerturbColor);
            viewer.data().add_edges(
                    igl::slice(this->SamplePerturb,this->SkeletonIndices0,1),
                    igl::slice(this->SamplePerturb,this->SkeletonIndices1,1),
                    Eigen::RowVector3d(0,0,0));
            for(auto &data : viewer.data_list){
//                data.set_colors(component_colors[data.id]);
                data.point_size = point_size;
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
            viewer.data().point_size= point_size;
            this->to_nearest();
            viewer.data().add_points(this->SampleInitial, this->TransportPlan*SamplePerturbColor);
            std::cout << "the color error for nearst neighbor ="
            << color_error(this->TransportPlan*SamplePerturbColor,SamplePerturbColor)<< std::endl;

            std::cout << "the function error in L2-norm cost ="
                      << color_error(this->TransportPlan* this->SampleVals, this->SampleVals)<<std::endl;
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
            this->Sinkhorn();
            viewer.selected_data_index = 0;
            viewer.data().point_size= point_size;
            viewer.data().add_points(this->SampleInitial, this->TransportPlan*SamplePerturbColor);
            std::cout << "the color error for L2-norm cost ="
                      << color_error(this->TransportPlan*SamplePerturbColor,SamplePerturbColor)<< std::endl;
//            std::cout << this->TransportPlan*SamplePerturbColor << std::endl;
            std::cout << "the function error in L2-norm cost ="
                      << color_error(this->TransportPlan* this->SampleVals, this->SampleVals)<<std::endl;
            break;
        case '6':
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
            this->VarSinkhorn();
            viewer.selected_data_index = 0;
            viewer.data().point_size= point_size;
            viewer.data().add_points(this->SampleInitial, this->TransportPlan*SamplePerturbColor);
            std::cout << "the color error for L2-norm cost ="
                      << color_error(this->TransportPlan*SamplePerturbColor,SamplePerturbColor)<< std::endl;
//            std::cout << this->TransportPlan*SamplePerturbColor << std::endl;

            std::cout << "the function error in L2-norm cost ="
                      << color_error(this->TransportPlan* this->SampleVals, this->SampleVals)<<std::endl;
            break;
        case '7':
            for(auto &data : viewer.data_list){
                data.clear();
            }
            this->_components_union();
            viewer.selected_data_index = 0;
            viewer.data().add_points(this->SamplePerturb, Eigen::RowVector3d(0, 0, 1));
//            viewer.data().set_mesh(this->TotalVerticesPerturb, this->TotalFacesPerturb);
            for(unsigned int i =0; i < this->ComponentsVertices.size(); ++i) {
                viewer.data_list[i+1].set_mesh(*(this->ComponentsVertices[i]), *(this->ComponentsFaces[i]));
                viewer.data_list[i+1].add_points(*(this->ComponentsSample[i]), Eigen::RowVector3d(1, 0, 0));
            }
            for(auto &data : viewer.data_list){
                data.set_colors(component_colors[data.id]);
                data.point_size = point_size;
            }

    }
}

void CutMesh::set_initial_from_json(const nlohmann::json & params){
    /*set the initial parameters from a json*/
    if(params.find("CutMesh_params")!= params.end()){
        auto CutMesh_params = params["CutMesh_params"];
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        igl::readOBJ(CutMesh_params["mesh_file"],V, F);
        ExpressionValue funcexpr;
        funcexpr.init(CutMesh_params["funcexpr"]);
        auto func = [&funcexpr](Eigen::Vector3d x)->double{return funcexpr(x[0],x[1],x[2]);};
        this->set_initial(V,F,CutMesh_params["sample_num"],func);
        if(CutMesh_params["cut"]){
            Eigen::Vector3d p0(0, 0, 0);
            Eigen::Vector3d p1(0, 0, 0);
            Eigen::Vector3d n1(0, 0, 1);
            Eigen::Vector3d n0(0, 1, 0.5);
            this->cut_with(p0, n0, p1, n1);
        }
        else{
            this->separate_cube_faces();
        }
        this->point_size = CutMesh_params["point_size"];
        this->perturb(CutMesh_params["random_seed"],CutMesh_params["perturb_range"],0.02, CutMesh_params["two_samples"]);
        if(params.find("Sinkhorn_params")!= params.end()) {
            auto Sinkhorn_params = params["Sinkhorn_params"];
            this->set_sinkhorn_const(
                    Sinkhorn_params["eps"],
                    Sinkhorn_params["threshold"],
                    Sinkhorn_params["max_iter"]);
            this->round = Sinkhorn_params["round"];
            this->compute_WeightMatrix(Sinkhorn_params["sigma"]);
        }
        this->build_graph(CutMesh_params["skeleton_range"],CutMesh_params["graph_valence"]);
        if(CutMesh_params["build_tree"]) this->build_tree();
//        this->build_Quad();

    }
}

void CutMesh::set_initial(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const int n,
                          std::function<double(Eigen::Vector3d)> func) {
    /*Set the initial mesh and samples*/
    this->Vertices = V;

    this->Faces = F;
    {
        Eigen::MatrixXi LocalI;
        Eigen::SparseMatrix<double> B;
        igl::random_points_on_mesh(n, V, F, B, LocalI);
        this->SampleInitial = B * V;
        this->SampleSourceFace = LocalI;
    }
    this->SampleNum = SampleInitial.rows();
    this->SamplePerturb = Eigen::MatrixXd::Zero(this->SampleNum, 3);
    this->SampleVals.resize(this->SampleNum);
    for (int i = 0; i < this->SampleNum; ++i) {
        SampleVals[i] = func(SampleInitial.row(i));
    }
}

void CutMesh::set_sinkhorn_const(const double eps, const double threshold, const int maxiter){
    this->SinkhornEps = eps;
    this->SinkhornThreshold = threshold;
    this->SinkhornMaxIter= maxiter;
}

void CutMesh::build_graph(double range, int m) {
    /*build  a m-nearst-enighbor graph */
    unsigned int num_components = this->ComponentsSample.size();
    std::vector<Eigen::Triplet<double> >Elastic_JK;
    Eigen::MatrixXd pair_Dense=Eigen::MatrixXd::Zero(this->SampleNum,this->SampleNum);
    this->SkeletonIndices0.resize(std::pow(this->SampleNum,2)/2,1);
    this->SkeletonIndices1.resize(std::pow(this->SampleNum,2)/2,1);
    int count = 0;
    for(int i =0; i< num_components; ++i){
        for(int j=0; j < this->ComponentsSampleIndices[i]->rows(); ++j){
            std::vector<std::pair<int, double> > j_stack;
            int idj = this->ComponentsSampleIndices[i]->coeff(j, 0);
            for (int k = 0; k < this->ComponentsSampleIndices[i]->rows(); ++k) {
                if(k==j) continue;
                int idk = this->ComponentsSampleIndices[i]->coeff(k, 0);
                double norm_jk = (this->SamplePerturb.row(idk) - this->SamplePerturb.row(idj)).norm();
                if(j_stack.size()==0 and k!=j){
                    j_stack.push_back(std::make_pair(idk, norm_jk));
                } else {
                    auto it = j_stack.begin();
                    while (it->second > norm_jk and it!= j_stack.end()) {
                        it++;
                        // loop to find a position to insert
                    }
                    if (it == j_stack.begin()) {
                        // if does not find within
                        if (j_stack.size() < m) {
                            // if the stack is not full
                            j_stack.insert(it, std::make_pair(idk, norm_jk));
                        }
                    } else {
                        // find a position to insert including j_stack.end()
                        j_stack.insert(it, std::make_pair(idk, norm_jk));
                        if (j_stack.size() > m) {
//                            j_stack.insert(it, std::make_pair(idk, norm_jk));
                            j_stack.erase(j_stack.begin(), j_stack.begin() + 1);
                        }
                    }
                }
            }
            // inserted the detected m-nearest-pair in to the Dense_matrix

            for(auto jt = j_stack.begin(); jt!= j_stack.end(); jt++){
                int idk = jt->first;
                std::cout << idj <<"," <<idk <<", " <<j_stack.size() << std::endl;
                pair_Dense(std::min(idj,idk), std::max(idj,idk))=jt->second;
            }
        }

        for(int j=0; j< pair_Dense.rows(); ++j){
            for(int k=0;k< pair_Dense.cols(); ++k){
                if(pair_Dense(j,k)!=0){
                    double norm_jk = pair_Dense(j,k);
                    Elastic_JK.push_back(Eigen::Triplet<double>(j, k, norm_jk));
                    (this->SkeletonIndices0)(count) = j;
                    (this->SkeletonIndices1)(count) = k;
                    count += 1;
                }
            }
        }
    }
    this->SkeletonIndices0.conservativeResize(count,1);
    this->SkeletonIndices1.conservativeResize(count,1);
    Eigen::SparseMatrix<double> temp(this->SampleNum,this->SampleNum);
    temp.setFromTriplets(Elastic_JK.begin(), Elastic_JK.end());
    this->ElasticTensor = temp;
}

void CutMesh::build_tree() {
    /* use kruskal minspaning tree algorithm to simplify
     * the knn grpah if necessary */
    this->ElasticTensor = Kruskal_MST(this->ElasticTensor);
    Eigen::VectorXi temp0(this->ElasticTensor.nonZeros());
    Eigen::VectorXi temp1(this->ElasticTensor.nonZeros());
    int count = 0;
    for (int k = 0; k <this->ElasticTensor.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(this->ElasticTensor, k); it; ++it) {
            temp0(count)= it.row();
            temp1(count)=it.col();
            count+=1;
        }
    }
    this->SkeletonIndices0 = temp0;
    this->SkeletonIndices1 = temp1;
    std::cout << "after min_spaning tree, the number of skeleton is"<< count << std::endl;
}

void CutMesh::build_Quad(){
    try {
        if (this->SkeletonIndices0.rows() < 1 or this->SkeletonIndices1.rows() < 1) {
            throw(this->SkeletonIndices0.rows());
        }
        std::vector<Eigen::Triplet<double> > Quad_JK;
        const int num  = this->SampleNum;
        const int sknum = this->SkeletonIndices0.rows();
        for(int i=0; i< sknum; ++i){
            int j = this->SkeletonIndices0.coeff(i,0);
            int k = this->SkeletonIndices1.coeff(i,0);
            std::cout << i << std::endl;
            for(int a=0; a < num; ++a) {
                for (int b = 0; b < num; ++b) {
                    double M_jakb = (
                            this->SamplePerturb.row(j)
                            -this->SamplePerturb.row(k)
                            -this->SampleInitial.row(a)
                            +this->SampleInitial.row(b)).squaredNorm();
//                    Quad_JK.push_back(Eigen::Triplet<double>(j* num+a, k* num+b, M_jakb));
//                    Quad_JK.push_back(Eigen::Triplet<double>(k* num+b, j* num+a, M_jakb));
                }
            }
            Eigen::SparseMatrix<double> temp(num * num, num * num);
            temp.setFromTriplets(Quad_JK.begin(), Quad_JK.end());
            this->QuadraticCostMatrix = temp;
        }

    } catch(int a){
        std::cout << "skeleton uninitialized" << a << " is too small" << std::endl;
    }
}

void CutMesh::separate_cube_faces(){
    std::map<int, std::pair<int, double> > level;
    level.emplace(0,std::make_pair(0, 0.5));
    level.emplace(1,std::make_pair(0, -0.5));
    level.emplace(2,std::make_pair(1, 0.5));
    level.emplace(3,std::make_pair(1, -0.5));
    level.emplace(4,std::make_pair(2, 0.5));
    level.emplace(5,std::make_pair(2, -0.5));
    int sum = 0;
    try {
        if (this->ComponentsVertices.size() != 0) {
            throw this->ComponentsVertices.size() != 0;
        }
        for(unsigned int idx=0; idx < 6; ++idx) {
            int entry = std::get<0>(level[idx]);
            double val = std::get<1>(level[idx]);

            MeshTriple TA = cube_surface_select(this->Vertices, this->Faces, idx, entry, val);
            Eigen::MatrixXd S_sub(this->SampleNum,3);
            Eigen::VectorXi S_idx(this->SampleNum);
            Eigen::VectorXi IdA = std::get<2>(TA);
            Eigen::MatrixXi FA=std::get<1>(TA);
            int count = 0;
            for(unsigned int i=0; i< this->SampleNum; ++i){
                bool found_in_face= false;

                for(unsigned int j= 0; j < IdA.rows();++j){
                    if(this->SampleSourceFace(i)== IdA(j)){
                        found_in_face = true;
                    }
                }
                if(found_in_face){
                    S_sub.row(count)=this->SampleInitial.row(i);
                    S_idx(count) = i;
                    count+=1;
                }
            }
            S_sub.conservativeResize(count,3);
            S_idx.conservativeResize(count);

            std::unique_ptr<Eigen::MatrixXd> ptrVA(new Eigen::MatrixXd()), ptrSA(new Eigen::MatrixXd());
            *ptrVA = std::get<0>(TA);
            *ptrSA = S_sub;
            std::unique_ptr<Eigen::MatrixXi> ptrFA(new Eigen::MatrixXi());
            std::unique_ptr<Eigen::VectorXi> ptrSIA(new Eigen::VectorXi());
            *ptrFA = std::get<1>(TA);
            *ptrSIA = S_idx;
            sum += ptrSA->rows();
            this->ComponentsVertices.push_back(std::move(ptrVA));
            this->ComponentsFaces.push_back(std::move(ptrFA));
            this->ComponentsSample.push_back(std::move(ptrSA));
            this->ComponentsSampleIndices.push_back(std::move(ptrSIA));
        }
        if(sum!= this->SampleNum){
            throw sum;
        }
    }
    catch (bool b){
        std::cout << "Exception accurred, do not separate twice"<<std::endl;
    }
    catch (int sum){
        std::cout << "Exception accurred, divided Sample num " <<
                  sum  << "!=" << "total sample num" << this->SampleNum << std::endl;
    }
}
void CutMesh::cut_with(const Eigen::Vector3d &p0,
        const Eigen::Vector3d &n0,
        const Eigen::Vector3d &p1,
        const Eigen::Vector3d &n1) {
    // slice with two planes to get 4 pieces of mesh
    int signs[4][2] = {{1,1},{-1,1},{1,-1},{-1,-1}};
    int sum = 0;
    try{
        if(this->ComponentsVertices.size()!=0){
            throw this->ComponentsVertices.size()!=0;
        }
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

        if(sum!= this->SampleNum){
            throw sum;
        }
    }
    catch (int e){
        std::cout << "Exception accurred, divided Sample num " <<
        sum  << "!=" << "totol sample num" << this->SampleNum << std::endl;
    }
    catch (bool b){
        std::cout << "Exception accurred, do not cut twice"<<std::endl;
    }
}

void CutMesh::perturb(const int seed, double max_shift, double min_shift, bool two_samples) {
    /*perturb the mesh compoenents with random shift vectors
     * and generate the perturbed samples*/
    unsigned int num_components = this->ComponentsVertices.size();
    std::mt19937 gen;
    gen.seed(seed);
    std::uniform_real_distribution<> dis(max_shift, min_shift);
    std::map<int, Eigen::Vector3i> signs;
    signs.emplace(0, Eigen::Vector3i(-1,-1,1));
    signs.emplace(1, Eigen::Vector3i(-1,1,1));
    signs.emplace(2, Eigen::Vector3i(1,1,-1));
    signs.emplace(3, Eigen::Vector3i(1,-1,1));
    signs.emplace(4, Eigen::Vector3i(1,-1,-1));
    signs.emplace(5, Eigen::Vector3i(-1,-1,-1));
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
            if(not two_samples) {
                this->SamplePerturb.row(this->ComponentsSampleIndices[i]->coeff(k, 0)) = this->ComponentsSample[i]->row(
                        k);
            } else{
//                this->SamplePerturb.row(this->ComponentsSampleIndices[i]->coeff(k, 0)) =
            }
        }
    }
}

void CutMesh::_components_union(){
    int NF=0;
    int NV =0;
    for(int i = 0; i< this->ComponentsFaces.size(); ++i){
        NF += this->ComponentsFaces[i]->rows();
        NV += this->ComponentsVertices[i]->rows();
    }
    this->TotalFacesPerturb= Eigen::MatrixXi::Zero(NF,3);
    this->TotalVerticesPerturb=Eigen::MatrixXd::Zero(NV,3);
    Eigen::VectorXi FacesSourceComponents = Eigen::VectorXi::Zero(NF);
    int countF = 0;
    int shiftF = 0;
    int countV = 0;
    for(int i =0 ; i < this->ComponentsFaces.size(); ++i){
        for(int j = 0; j < this->ComponentsFaces[i]->rows();++j){
            this->TotalFacesPerturb.row(countF) = this->ComponentsFaces[i]->row(j)
                    + Eigen::RowVectorXi::Constant(3, shiftF);
            FacesSourceComponents(countF)= i;
            countF += 1;
        }
        for(int k =0; k < this->ComponentsVertices[i]->rows();++k){
            this->TotalVerticesPerturb.row(countV) = this->ComponentsVertices[i]->row(k);
            countV += 1;
        }
        shiftF += this->ComponentsVertices[i]->rows();
    }
    for(int i =0 ; i< 4;++i){
        std::cout << this->ComponentsSample[i]->rows() << std::endl;
    }
    this->ComponentsSample.resize(0);
    {
        Eigen::MatrixXi LocalI;
        Eigen::SparseMatrix<double> B;
        igl::random_points_on_mesh(this->SampleNum, this->TotalVerticesPerturb, this->TotalFacesPerturb, B, LocalI);
        this->SamplePerturb = B * this->TotalVerticesPerturb;
        this->SampleSourceFace = LocalI;

        for (int i = 0; i < this->ComponentsFaces.size(); ++i) {
            std::unique_ptr<Eigen::MatrixXd> ptr_SA(new Eigen::MatrixXd());
            std::unique_ptr<Eigen::VectorXi > ptr_SIA(new Eigen::VectorXi());
            SamplePair SA = components_sample_select(
                    this->SamplePerturb,
                    this->SampleSourceFace,
                    FacesSourceComponents, i);
            *(ptr_SA) = std::get<0>(SA);
            *(ptr_SIA) = std::get<1>(SA);
            this->ComponentsSample[i]= std::move(ptr_SA);
            std::cout << this->ComponentsSample[i]->rows() << std::endl;
            this->ComponentsSampleIndices[i]= std::move(ptr_SIA);
        }
    }
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

void CutMesh::compute_CostMatrix(const Eigen::MatrixXd & source, const Eigen::MatrixXd & target, char option) {
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
                std::cout << this->SampleInitial << std::endl;
                std::cout << "-------------------" << std::endl;
                std::cout << this->SamplePerturb <<std::endl;
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        this->CostMatrix(i, j) = (this->SampleInitial.row(i) - this->SamplePerturb.row(j)).norm();
                    }
                }
                break;
            case 'b':
                for (int i = 0; i < target.rows(); ++i) {
                    for (int j = 0; j < source.rows(); ++j) {
                        this->CostMatrix(i, j) = (this->SampleInitial.row(i) - this->SamplePerturb.row(j)).lpNorm<Eigen::Infinity>();
                    }
                }
                break;
            case 'c':
                this->locate_Centers(this->WeightMatrixPerturb, this->TransportPlan, this->SamplePerturb, this->CentersPerturb, 'o');
                this->locate_Centers(this->WeightMatrixInitial, this->TransportPlan.transpose(), this->SampleInitial, this->CentersInitial, 'o');
                int num = this->SampleNum;
                Eigen::MatrixXd D_centersInitial(num, num);
                Eigen::MatrixXd D_centersPerturb(num, num);
                // store the squared distance of the samples from the centers
                for(unsigned int i =0; i < num; ++i){
                    for(unsigned int j =0 ; j < num; ++j){
                        D_centersInitial(i,j) = (this->CentersInitial.row(j)-this->SamplePerturb.row(i)).squaredNorm();
                        D_centersPerturb(i,j) = (this->CentersPerturb.row(j)-this->SampleInitial.row(i)).squaredNorm();
                    }
                }
//                std::cout << D_centerPerturb << std::endl;
//                std::cout << this->WeightMatrixPerturb;
                this->CostMatrix =  D_centersPerturb * this->WeightMatrixPerturb
                        + D_centersInitial * this->WeightMatrixInitial;
                break;
        }
    }
    catch(std::string a){
        std::cout << a << std::endl;
    }
}

void CutMesh::Sinkhorn(){
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
        if (this->SinkhornEps==0 or this->SinkhornMaxIter==0){
            throw "Sinkhorn constant not initialized";
        }
    }
    catch(std::string a){
        std::cout << "Exception accurred"<<
        " "<< a<<std::endl;
    }
    int n = this->SampleNum;
    this->TransportPlan = sinkhorn(
            Eigen::VectorXd::Constant(n,1),
            Eigen::VectorXd::Constant(n,1),
            this->CostMatrix, this->SinkhornEps, this->SinkhornMaxIter, this->SinkhornThreshold);
    if(this->round){
        for(unsigned int i=0; i< this->SampleNum; ++i){
            int max_idx = 0;
            double max_val = this->TransportPlan(i,0);
            this->TransportPlan(i,0)=1;
            for(unsigned int j=1; j< this->SampleNum; ++j){
                double cur_val = this->TransportPlan(i,j);
                if(cur_val > max_val){
                    max_val = cur_val;
                    this->TransportPlan(i, max_idx) = 0;
                    this->TransportPlan(i,j) = 1;
                    max_idx = j;
                }
                else{
                    this->TransportPlan(i,j) = 0;
                }
            }
        }
    }
    std::cout<<" Total Cost for this TransportPlan=" << (this->TransportPlan.array() * this->CostMatrix.array()).sum()
    << '\n'<<" Totoal Cost for identity="
    << (Eigen::MatrixXd::Identity(n,n).array() * this->CostMatrix.array()).sum() <<std::endl;
}

void CutMesh::VarSinkhorn(){
    // this is the two loop Sinkhorn described in the paper;
    // for the first step initialize the centers by the wasserstein metric
    int n = this->SampleNum;
//    std::cout << this->CostMatrix<<std::endl;
    this->TransportPlan = sinkhorn(
            Eigen::VectorXd::Constant(n,1),
            Eigen::VectorXd::Constant(n,1),
            this->CostMatrix, this->SinkhornEps, this->SinkhornMaxIter, this->SinkhornThreshold);
//    std::cout << -1 << "------" << this->CostMatrix.block(0,0,10,10);
    for(int i =0 ; i < 6 ; ++i){
        this->compute_CostMatrix(this->SamplePerturb, this->SampleInitial,'c');
//        std::cout << i << "------" << this->CostMatrix.block(0,0,10,10);
        this->TransportPlan = sinkhorn(
            Eigen::VectorXd::Constant(n,1),
            Eigen::VectorXd::Constant(n,1),
            this->CostMatrix, this->SinkhornEps, (i+1) * 20, this->SinkhornThreshold);

    }

    if(this->round){
        for(unsigned int i=0; i< this->SampleNum; ++i){
            int max_idx = 0;
            double max_val = this->TransportPlan(i,0);
            this->TransportPlan(i,0)=1;
            for(unsigned int j=1; j< this->SampleNum; ++j){
                double cur_val = this->TransportPlan(i,j);
                if(cur_val > max_val){
                    max_val = cur_val;
                    this->TransportPlan(i, max_idx) = 0;
                    this->TransportPlan(i,j) = 1;
                    max_idx = j;
                }
                else{
                    this->TransportPlan(i,j) = 0;
                }
            }
        }
    }
    std::cout<<" Totoal Cost for this TransportPlan=" << (this->TransportPlan.array() * this->CostMatrix.array()).sum()
    << '\n'<<" Totoal Cost for identity="
    << (Eigen::MatrixXd::Identity(n,n).array() * this->CostMatrix.array()).sum() <<std::endl;

}

void CutMesh::NewtonSinkhorn(){}

void CutMesh::locate_Centers(
    const Eigen::SparseMatrix<double> & weightmatrix,
    const Eigen::MatrixXd & transportplan,
    const Eigen::MatrixXd & sample,
    Eigen::MatrixXd & centers,
    char options){
    switch(options) {
        case 'o': // stand for original formula in the paper
            centers = weightmatrix * transportplan * sample;
            break;
    }

}




void CutMesh::compute_WeightMatrix(double sigma) {
    // TODO add face adjacency support.
    // this is a temporary solution to compute the wieght matrix;
    // because the sample on the bad mesh does not know the face indices of faces
    // we temporarily build it on each compoenent, later we will build the weight matrix
    // based on the adjacency of faces on bad mesh.
    // In practical usage, the the sample on bad mesh is generated randomly
    // rather then (currently for evaluation reason) generated on good mesh and then cut
    // Compute the weightmatrix W_ij=W_{x_i}(x_j)
    // it start with a dense matrix and then transfer into a sparsematrix
    // does not start with triplets because we need to normalize each row;
    std::vector<Eigen::Triplet<double> > Wxx_triplets;
    Eigen::MatrixXd Wxx_Dense=Eigen::MatrixXd::Zero(this->SampleNum,this->SampleNum);
    unsigned int num_components = this->ComponentsSample.size();
    for(int i =0; i< num_components; ++i) {
        for (int j = 0; j < this->ComponentsSampleIndices[i]->rows() - 1; ++j) {
            std::cout << j << std::endl;
            for (int k = j + 1; k < this->ComponentsSampleIndices[i]->rows(); ++k) {
                int idj = this->ComponentsSampleIndices[i]->coeff(j,0);
                int idk = this->ComponentsSampleIndices[i]->coeff(k,0);
                double norm_jk = (this->SamplePerturb.row(idk)-this->SamplePerturb.row(idj)).norm();
                double exp_jk = std::exp(-std::pow(norm_jk/sigma,2));
                if(exp_jk>0.5){
                    Wxx_Dense(idj,idk)=exp_jk;
                    Wxx_Dense(idk,idj)=exp_jk;
                }
            }
        }
    }
    for(int i =0; i< Wxx_Dense.rows();++i){
        double su = Wxx_Dense.row(i).sum();
        if(su > 0) {
            Wxx_Dense.row(i) /= Wxx_Dense.row(i).sum();
        }
    }
    for(int i =0; i< Wxx_Dense.rows();++i){
        for(int j =0 ; j< Wxx_Dense.cols();++j){
            if(Wxx_Dense.coeff(i,j)!=0){
                Wxx_triplets.push_back(Eigen::Triplet<double>(i,j,Wxx_Dense.coeff(i,j)));
            }
        }
    }

    Eigen::SparseMatrix<double> temp(this->SampleNum,this->SampleNum);
    temp.setFromTriplets(Wxx_triplets.begin(),Wxx_triplets.end());
    this->WeightMatrixPerturb = temp;
    weight_matrix(sigma,
            this->SampleInitial,
            this->Faces,
            this->SampleSourceFace,
            this->WeightMatrixInitial
            );
}

// helper function to build weightmatrix
void weight_matrix(
        double sigma,
        const Eigen::MatrixXd  & sample,
        const Eigen::MatrixXi & faces,
        const Eigen::MatrixXi & samplesource_indice,
        Eigen::SparseMatrix<double> & wmtx){
    //TODO add face adjacency support.
    int num = sample.rows();
    std::vector<Eigen::Triplet<double> > Wyy_triplets;
    Eigen::MatrixXd Wyy_Dense=Eigen::MatrixXd::Zero(num,num);
    for (int j = 0; j < num-1; ++j) {
        for (int k = j + 1; k < num; ++k) {
            double norm_jk = (sample.row(k)-sample.row(j)).norm();
            double exp_jk = std::exp(-std::pow(norm_jk/sigma,2));
            if(exp_jk>0.5){
                Wyy_Dense(j,k)=exp_jk;
                Wyy_Dense(k,j)=exp_jk;
            }
        }
    }
    for(int i =0; i< Wyy_Dense.rows();++i){
        Wyy_Dense.row(i) /= Wyy_Dense.row(i).sum();
    }
    for(int i =0; i< Wyy_Dense.rows();++i){
        for(int j =0 ; j< Wyy_Dense.cols();++j){
            if(Wyy_Dense.coeff(i,j)!=0){
                Wyy_triplets.push_back(Eigen::Triplet<double>(i,j,Wyy_Dense.coeff(i,j)));
            }
        }
    }
    Eigen::SparseMatrix<double> temp(num,num);
    temp.setFromTriplets(Wyy_triplets.begin(),Wyy_triplets.end());
    wmtx= temp;

}
