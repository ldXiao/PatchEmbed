//
// Created by Lind Xiao on 7/25/19.
//

#include "CutGraph.h"
#include <Eigen/Core>
#include <iostream>
#include <map>
#include <tuple>
#include <igl/opengl/glfw/Viewer.h>
void CutGraph::_set_Vertices(Eigen::MatrixXd V) {
    int num = V.rows();
    this->_vertices.resize(num);
    this->Vertices = V;
    {
        int count = 0;
        for (auto it = this->_vertices.begin(); it != this->_vertices.end(); ++it) {
            *it = std::make_tuple(V(count,0),V(count,1), V(count,2));
            count += 1;
        }
    }
}

void CutGraph::_set_Edges_from_KNN(int m) {
    int count = 0;
    int esize = m * this->_vertices.size();
    this->_edges.resize(esize);
    this->Edges= Eigen::MatrixXi::Zero(esize,2);
    std::map<std::pair<int, int>, int > map;
    for(int j=0; j < this->Vertices.rows(); ++j){
        std::vector<std::pair<int, double> > j_stack;
        for (int k = 0; k < this->Vertices.rows(); ++k) {
            if(k==j) continue;
            double norm_jk = (this->Vertices.row(k) - this->Vertices.row(j)).norm();
            if(j_stack.size()==0 and k!=j){
                j_stack.push_back(std::make_pair(k, norm_jk));
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
                        j_stack.insert(it, std::make_pair(k, norm_jk));
                    }
                } else {
                    // find a position to insert including j_stack.end()
                    j_stack.insert(it, std::make_pair(k, norm_jk));
                    if (j_stack.size() > m) {
                        j_stack.erase(j_stack.begin(), j_stack.begin() + 1);
                    }
                }
            }
        }
        // inserted the detected m-nearest-pair in to the this->Edges

        for(auto jt = j_stack.begin(); jt!= j_stack.end(); jt++){
            int k = jt->first;
            auto find0 = map.find(std::make_pair(j,k));
            auto find1 = map.find(std::make_pair(k,j));
            if( find0 == map.end() and find1 == map.end()) {
                // insert each edge only once
                this->_edges[count] = std::make_pair(j, k);
                Eigen::RowVector2i temp;
                temp << j,k;
                this->Edges.row(count) = temp;
                count += 1;
                map[std::make_pair(j,k)]=1;
            }
        }
    }
    this->_edges.resize(count);
    this->Edges.conservativeResize(count,2);
}


