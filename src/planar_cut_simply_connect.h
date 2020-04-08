#ifndef BCCLEAN_PLANARCUT_H
#define BCCLEAN_PLANARCUT_H
#include <Eigen/Core>
#include <vector>
#include <igl/bfs.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Sparse>
#include <algorithm>
#include "MatchMakerTree.h"
#include <nlohmann/json.hpp>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#include <igl/matrix_to_list.h>
#include <igl/extract_manifold_patches.h>
namespace bcclean{
    void planar_cut_simply_connect(
        Eigen::MatrixXd & V, // Vraw
        Eigen::MatrixXi & F,  // Fraw
        Eigen::MatrixXd & Vbase,
        Eigen::MatrixXi & Fbase, 
        Eigen::VectorXi & VI, // vertex mapping Vraw -> Vbase
        Eigen::VectorXi & FI, // face mapping Fraw -> Fbase
        Eigen::VectorXi & FL, // face label on Fbase
        const std::vector<std::vector<int> > & boundary_loops, 
        std::vector<bool> & VCuts, 
        std::vector<std::vector<bool> > & TCuts);

   

}
#endif