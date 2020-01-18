
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <algorithm>
#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include "Kruskal.h"
namespace bcclean{

namespace Algo{

std::vector<Eigen::Triplet<double>> to_triplets(Eigen::SparseMatrix<double> & M){
    std::vector<Eigen::Triplet<double>> v;
    for(int i = 0; i < M.outerSize(); i++)
        for(typename Eigen::SparseMatrix<double>::InnerIterator it(M,i); it; ++it)
            v.emplace_back(it.row(),it.col(),it.value());
    return v;
}

std::vector<int> get_vertices(const std::vector<std::pair<int,std::pair<int, int> > > frame_graph)
{
    std::vector<int> universe;
    std::unordered_map<int, int> set;
    for(auto edg: frame_graph)
    {
        if(set.find(edg.second.first)== set.end())
        {
            universe.push_back(edg.second.first);
            set[edg.second.first] = edg.first;
        }
        if(set.find(edg.second.second)== set.end())
        {
            universe.push_back(edg.second.second);
            set[edg.second.second] = edg.first;
        }
    }    
    return universe;
}

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

std::vector<std::pair<int, std::pair<int, int> > > Kruskal_MST(const std::vector<std::pair<int, std::pair<int, int> > > & frame_graph)
{
    std::vector<std::pair<int, std::pair<int, int> > > frame_tree;
    std::vector<std::pair<int, std::pair<int, int> > > frame_graph_copy = frame_graph;
    std::sort(frame_graph_copy.begin(), frame_graph_copy.end(),
        [](const std::pair<int, std::pair<int, int> > & a, const std::pair<int, std::pair<int, int> > & b)
        {
            return a.first < b.first;
        });
    std::vector<int> universe = get_vertices(frame_graph_copy);
    DisjointSet DS;
    DS.initialize(universe);
    for(auto edg: frame_graph_copy)
    {
        int x = edg.second.first;
        int y = edg.second.second;
        if(DS.find_root(x)!= DS.find_root(y))
        {
            DS.make_union(x,y);
            frame_tree.push_back(edg);
        }
    }
    return frame_tree;
}

std::vector<int> Graph_BFS(const std::vector<std::pair<int, std::pair<int, int> > > & frame_graph, const int root)
{
    std::vector<int> res;
    bool find_root = false;
    std::map<int, bool> vertices_visited;
    std::vector<int> vertices = get_vertices(frame_graph);
    for(auto  vv: vertices)
    {
        vertices_visited[vv]= false; // initialize vertices_visited
        if(vv == root) find_root = true;
    }

    if(! find_root)
    {
        exit(EXIT_FAILURE);
    }
    std::queue<int> visit_queue;
    visit_queue.push(root);
    vertices_visited[root]= true;
    while(visit_queue.size()!=0)
    {
        int cur_v = visit_queue.front();
        bool neighbor_visited = true;
        for(auto item: frame_graph)
        {
            int other= -1;
            if(item.second.first==cur_v || item.second.second == cur_v)
            {
                if(item.second.first== cur_v)
                {
                    other = item.second.second;
                    if(vertices_visited[other])
                    {
                        continue;
                    }
                    else
                    {
                        vertices_visited[other] = true;
                        visit_queue.push(other);
                    }
                    
                }
                else 
                {
                    other = item.second.first;
                    if(vertices_visited[other])
                    {
                        continue;
                    }
                    else
                    {
                        vertices_visited[other] = true;
                        visit_queue.push(other);
                    }
                }
            }
        }
        res.push_back(cur_v);
        visit_queue.pop();
    }
    return res;
}

void remove_vertices(
    std::vector<std::pair<int, std::pair<int, int> > > & frame_graph,
    const std::vector<int> & remove_list)
{
    auto remove_cond=[remove_list](const std::pair<int, std::pair<int, int> > edg)-> bool
    {
        return (std::find(remove_list.begin(), remove_list.end(), edg.second.first) != remove_list.end()) || (std::find(remove_list.begin(), remove_list.end(), edg.second.second) != remove_list.end());
    };
    frame_graph.erase(std::remove_if(frame_graph.begin(), frame_graph.end(), remove_cond), frame_graph.end());
    return;
}

void addInOrder(std::vector<int> & res, std::vector<int> & constriant_list, std::map<int, std::vector<int> > & graph_node_dict,  const std::vector<int> & forest_list)
{
    
    for(auto p: forest_list)
    {
        res.push_back(p);
        std::vector<int> node_p = graph_node_dict.at(p);
        if(node_p.size()>0)
        {
            std::vector<int> remove_constr;
            for(auto cc : constriant_list)
            {
                std::vector<int> & node_cc = graph_node_dict[cc];
                for(auto nd: node_p)
                {
                    node_cc.erase(std::remove(node_cc.begin(), node_cc.end(), nd), node_cc.end());
                }
                if(node_cc.size()==0)
                {
                    remove_constr.push_back(cc);
                    res.push_back(cc);
                }
            }
            if(remove_constr.size()>0)
            {
                for(auto ccr: remove_constr)
                {
                    constriant_list.erase(std::remove(constriant_list.begin(), constriant_list.end(), ccr), constriant_list.end());
                }
            }
        }
    }
}

std::vector<int> Constriant_Graph_BFS(
    const std::vector<std::pair<int, std::pair<int, int> > > & frame_graph,
    const std::map<int, std::vector<int> > & graph_node_dict,
    const std::vector<int> & constriant_list,
    const int root
)
{
    // inputs:
    // frame_graph (edge_idx ,(patch_idx, patch_idx)) of patch connectivity
    // graph_node_dict {patch_idx: [node_idx,....]} of graph_node_dict
    // constriant_list [patch_idx]
    // int root stating patch
    int forest_root = root;
    std::vector<int> res;
    std::vector<int> forest_list;
    std::vector<std::pair<int, std::pair<int, int> > > frame_graph_copy = frame_graph;
    std::vector<int> constriant_list_copy = constriant_list;
    std::map<int, std::vector<int> > graph_node_dict_copy = graph_node_dict;
    remove_vertices(frame_graph_copy, constriant_list);
    std::vector<int> reference = Graph_BFS(frame_graph,root);
    for(auto cc: constriant_list)
    {
        reference.erase(std::remove(reference.begin(),reference.end(), cc), reference.end());
    }
    while(reference.size()>0)
    {
        forest_root = reference[0];
        forest_list = Graph_BFS(frame_graph_copy, forest_root);
        for(auto ff: forest_list)
        {
            reference.erase(std::remove(reference.begin(),reference.end(), ff), reference.end());
        }
        remove_vertices(frame_graph_copy, forest_list);
        addInOrder(res, constriant_list_copy, graph_node_dict_copy, forest_list);
    }
    if(constriant_list_copy.size()>0)
    {
        for(auto cc: constriant_list_copy)
        {
            res.push_back(cc);
        }
    }
    return res;

}

std::vector<int> MST_BFS(const std::vector<std::pair<int, std::pair<int, int> > > & frame_MST)
{
    std::vector<int> res;
    for(auto item: frame_MST){
        if(std::find(res.begin(), res.end(),item.second.first)== res.end())
        {
            res.push_back(item.second.first);
        }
        if(std::find(res.begin(), res.end(), item.second.second)== res.end())
        {
            res.push_back(item.second.second);
        }
    }
    return res;
}
}
}


