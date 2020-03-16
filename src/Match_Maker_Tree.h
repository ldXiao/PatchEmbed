#ifndef BCCLEAN_MATCH_MAKER_TREE_H
#define BCCLEAN_MATCH_MAKER_TREE_H 
#include <Eigen/Core>
#include "patch.h"
#include "node.h"
#include "edge.h"
#include "kdtree_NN_Eigen.hpp"
#include <unordered_map>
#include <map>
#include <cstdlib>
#include <ctime>
#include <igl/dijkstra.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include <algorithm>
#include <functional>
#include <cmath>
#include "params.h"
#include "CellularGraph.h"
#include "TraceComplex.h"
namespace bcclean{
namespace MatchMaker{

    int _locate_seed_face(
        const CellularGraph & cg, 
        const TraceComplex & tc, 
        const int patch_idx);

    void silence_vertices(std::vector<std::vector<int > > & VV, const std::vector<int> & silent_indices);

    void update_local_sector(
        const std::vector<std::vector<int> > & VV, 
        const Eigen::MatrixXi & F,
        const std::map<int , std::map<int, bool> > & node_edge_visit_dict,
        const std::map<int, std::vector<int> > & node_edge_dict,
        const std::vector<std::vector<int> > & TEdges,
        const std::map<int, std::vector<int> > & CC_node_face_dict, 
        const int& source,
        const int& target,
        const int& cur_edge,
        std::vector<std::vector<int> > & VV_temp
    );


    bool MatchMakerTree(
        const CellularGraph & cg,
        Eigen::MatrixXd & V_good,
        Eigen::MatrixXi & F_good,
        Eigen::VectorXi & FL_good,
        const params param
    ); 
}
}

#endif // BCCLEAN_MATCH_MAKER_TREE_H
