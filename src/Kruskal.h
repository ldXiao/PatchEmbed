//
// Created by Lind Xiao on 6/20/19.
//
#ifndef BCCLEAN_KRUSKAL_H
#define BCCLEAN_KRUSKAL_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <algorithm>
#include <iostream>
#include <vector>
#include <unordered_map>
#include "Kruskal.h"
namespace bcclean{

namespace Algo{
class DisjointSet{
    std::unordered_map<int, int> root;
    //constructor
public:
    void initialize(const std::vector<int> & universe){
        for(int i: universe){
            root[i]=i;
        }
    }

    int find_root(int k){
        while(root[k]!=k){
            k = root[k];
        }
        return k;
    }

    void make_union(int i, int j){
        int x = find_root(i);
        int y = find_root(j);

        root[x]=y;
    }
};

std::vector<Eigen::Triplet<double>> to_triplets(Eigen::SparseMatrix<double> & M);

std::vector<int> get_vertices(const std::vector<std::pair<int,std::pair<int, int> > > frame_graph);

std::vector<int> get_vertices(const Eigen::SparseMatrix<double> & Graph);

Eigen::SparseMatrix<double> Kruskal_MST(Eigen::SparseMatrix<double> Graph);

std::vector<std::pair<int, std::pair<int, int> > > Kruskal_MST(const std::vector<std::pair<int, std::pair<int, int> > > & frame_graph);

std::vector<int> MST_BFS(const std::vector<std::pair<int, std::pair<int, int> > > & frame_MST);
std::vector<int> Graph_BFS(const std::vector<std::pair<int, std::pair<int, int> > > & frame_graph, const int root);
std::vector<int> Constriant_Graph_BFS(
    const std::vector<std::pair<int, std::pair<int, int> > > & frame_graph,
    const std::map<int, std::vector<int> > & graph_node_dict,
    const std::vector<int> & constriant_list,
    const int root
);
}
}




#endif