#include "mapping_patch.h"
#include <Eigen/Core>
#include <iostream>
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
    bcclean::mapping_patch mp;
    mp.build_patch(V, F, n2, 0);
    for(auto item : mp._edge_arc_ratio_list){
        std::cout << item.first<<", "<<item.second<<std::endl;
    }
    std::cout << mp._bnd_uv <<std::endl;
}