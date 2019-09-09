#include "mapping_patch.h"
#include <Eigen/Core>
#include <iostream>
int main(){
    bcclean::node a1,a2,a3, b1, b2, b3;
    Eigen::RowVector3d x(1,0,0);
    Eigen::RowVector3d y(0,1,0);
    Eigen::RowVector3d z(0,0,1);
    std::vector<int> l1, l2, l3;
    l1.push_back(0);
    l1.push_back(1);
    l1.push_back(2);
    l2.push_back(0);
    l2.push_back(2);
    l2.push_back(3);
    l3.push_back(0);
    l3.push_back(1);
    l3.push_back(3);
    a1.initialize(4, x, l1);
    a2.initialize(4, y, l2);
    a3.initialize(4, z, l3);
    b1.initialize(4, y, l1);
    b2.initialize(4, z, l2);
    b3.initialize(4, x, l3);

    std::vector<bcclean::node> n1 {a1,a2,a3};
    std::vector<bcclean::node> n2 {b3,b2, b1};
    std::map<int, int> mapping;
    std::cout<<bcclean::cyc_flip_mapping(n1, n2, mapping)<<std::endl;
    for(auto item:mapping){
        std::cout<< item.first<<"->"<<item.second<<std::endl;
    }

}