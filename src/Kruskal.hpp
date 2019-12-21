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

std::vector<Eigen::Triplet<double>> to_triplets(Eigen::SparseMatrix<double> & M){
    std::vector<Eigen::Triplet<double>> v;
    for(int i = 0; i < M.outerSize(); i++)
        for(typename Eigen::SparseMatrix<double>::InnerIterator it(M,i); it; ++it)
            v.emplace_back(it.row(),it.col(),it.value());
    return v;
}

std::vector<int> get_vertices(const std::vector<std::pair<int,std::pair<int, int> > > frame_graph)

std::vector<int> get_vertices(const Eigen::SparseMatrix<double> & Graph){
    std::unordered_map<int,int> set;
    std::vector<int> universe;
    for(int i = 0; i < Graph.outerSize(); i++) {
        for (typename Eigen::SparseMatrix<double>::InnerIterator it(Graph, i); it; ++it) {
            if(set.find(it.row())==set.end()){
                set[it.row()]=1;
                universe.push_back(it.row());
            }
            if(set.find(it.col())==set.end()){
                set[it.col()]=1;
                universe.push_back(it.col());
            }
        }
    }
    return universe;
}

Eigen::SparseMatrix<double> Kruskal_MST(Eigen::SparseMatrix<double> Graph){
    //this graph input is sparse matrix, it is by default free from parallel edges
    //we can further make sure the the diagonal terms are all zero in practical use
    auto Graph_triplets = to_triplets(Graph);
    std::vector<Eigen::Triplet<double> > Tree_triplets;
    std::sort(Graph_triplets.begin(), Graph_triplets.end(),
              [](const Eigen::Triplet<double> & a, const Eigen::Triplet<double> & b) -> bool
              {
                  return a.value() < b.value();
              });
    std::vector<int> universe = get_vertices(Graph);
    DisjointSet DS;
    DS.initialize(universe);
    for(auto it=Graph_triplets.begin(); it!= Graph_triplets.end();++it){
        int x = it->row();
        int y = it->col();
        if(DS.find_root(x)!= DS.find_root(y)){
            DS.make_union(x,y);
            Tree_triplets.push_back(*it);
        }
    }
    Eigen::SparseMatrix<double> Tree(Graph.rows(),Graph.cols());
    Tree.setFromTriplets(Tree_triplets.begin(), Tree_triplets.end());
    return Tree;
}

#endif

