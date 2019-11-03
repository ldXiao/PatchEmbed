#include "mapping_patch.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <iostream>
#include <map>
#include <vector>
#include <iostream>
#include <igl/read_triangle_mesh.h>
#include <igl/readDMAT.h>
#include <igl/opengl/glfw/Viewer.h>
#include "bcclean.h"
#include "graphcut_cgal.h"
#include <igl/readMSH.h>
#include <igl/remove_unreferenced.h>
#include <igl/predicates/predicates.h>
#include <igl/vertex_triangle_adjacency.h>
#include "edge.h"
int main(){
    bcclean::node a1,a2,a3, b1, b2, b3, b4;
    Eigen::RowVector3d x(1,0,0);
    Eigen::RowVector3d y(0,1,0);
    Eigen::RowVector3d z(-1,0,0);
    Eigen::RowVector3d w(0,-1,0);
    Eigen::RowVector3d c(0,0,0);
    Eigen::MatrixXd V(5,3);
    V << 1, 0, 0,
         0, 1, 0,
         -1, 0, 0,
         0, -1, 0,
         0, 0, 0;
    Eigen::MatrixXi F(4,3);
    F << 4, 0, 1,
         4, 1, 2,
         4, 2, 3,
         4, 3, 0;
    std::vector<int> l1, l2, l3, l4;
    l1.push_back(0);
    l1.push_back(1);
    l1.push_back(2);
    l2.push_back(0);
    l2.push_back(2);
    l2.push_back(3);
    l3.push_back(0);
    l3.push_back(3);
    l3.push_back(4);
    l4.push_back(0);
    l4.push_back(4);
    l4.push_back(1);
    a1.initialize(5, x, l1);
    a2.initialize(5, y, l2);
    a3.initialize(5, z, l3);
    b1.initialize(5, y, l1);
    b2.initialize(5, z, l2);
    b3.initialize(5, x, l3);
    b4.initialize(5, w, l4);


    std::vector<bcclean::node> n1 {a1,a2,a3};
    std::vector<bcclean::node> n2 {b4, b3, b2, b1};
    std::vector<bcclean::node> n3 {b3, b2, b1, b4};
    bcclean::mapping_patch mp;
    bcclean::mapping_patch mp1;
    mp.build_patch(V, F, n2, 0);
    mp1.build_patch(V, F, n3, 0);
    std::map<int, int> mapping;
    bcclean::cyc_flip_mapping(n2, n3, mapping);
    mp1.adjust_fit_target(n2);
    Eigen::SparseMatrix<double> Bmp(V.rows(), V.rows());
    Eigen::VectorXi FImp;
    

    for(auto item : mp._edge_arc_ratio_list){
        std::cout << item.first<<", "<<item.second<<std::endl;
    }
    std::cout << mp._bnd_uv <<std::endl;
     Eigen::MatrixXd V1(5,2);
    V1 << 0.25, 0.25,
         0.25, -0.25,
         -0.25, 0.25,
         -0.25, -0.25,
         0, 0.1;
    Eigen::Vector2d x1;
    x1 << 0, 0.1;
    Eigen::Vector2d x2;
    x2 << -1, 0;
    Eigen::Vector2d x3;
    x3 << 0, 1;
    Eigen::Vector2d x4;
    x4 << 1, 0;
    switch(igl::predicates::incircle(x1, x2, x3, x4)){
        case igl::predicates::Orientation::INSIDE: std::cout << "inside" << std::endl; break;
        case igl::predicates::Orientation::OUTSIDE: std::cout << "outside" << std::endl; break;
    }
    std::vector<Eigen::Triplet<double> > tp;
    Eigen::SparseMatrix<double> B(5,5);
    
    // Eigen::VectorXi FI;
    // bcclean::project_check(mp, V1, tp, FI);
    // B.setFromTriplets(tp.begin(), tp.end());
    // std::cout<< FI << std::endl;
    // std::cout<< "---------------" << std::endl;
    // std::cout<< mp._V_uv << std::endl;
    // std::cout<< mp._F_uv << std::endl;
    // std::cout<< "---------------" << std::endl;
    // std::cout<< B << std::endl;
    
    // std::vector<std::vector<int> > VF;
    // std::vector<std::vector<int> > VFi;
    // igl::vertex_triangle_adjacency(mp._V_raw, mp._F_raw, VF, VFi);
    // int ida =0;
    // std::cout << ida <<std::endl;
    // if(bcclean::inject_identity_uv(mp, mp1, Bmp, FImp)){
    //     std::cout<< "---------------" << std::endl;
    //     std::cout << Bmp <<std::endl;
    //     std::cout<< "---------------" << std::endl;
    //     std::cout << FImp << std::endl;
    // }
    // std::cout<< "---------------" << std::endl;
    bcclean::pair_map<std::pair<int,int>, bcclean::node> aada; 
    aada[std::pair<int,int>(1,2)]=a1;
    std::cout <<aada[std::pair<int,int>(1,2)]._position <<std::endl;
}