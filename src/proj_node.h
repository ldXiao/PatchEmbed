
#ifndef BCCLEAN_PROJ_NODE_H
#define BCCLEAN_PROJ_NODE_H
#include <vector>
#include <Eigen/Core>
#include <map>
#include "CellularGraph.h"
#include "TraceComplex.h"
namespace bcclean {
    void proj_node(
        const Eigen::MatrixXd & Vbad,
        const Eigen::MatrixXi & Fbad,
        const std::vector<int> & node_list_bad,
        Eigen::MatrixXd & Vgood,
        Eigen::MatrixXi & Fgood,
        std::map<int, int> & node_map
    );
    // try to project the nodes into good mesh and modify the good mesh correspondingly
    // if it can not find intersection use nearest neighbor instead 
    // always gaurantee that node_map is injective
    void proj_node_loop(
        const Eigen::MatrixXd & Vbad,
        const Eigen::MatrixXi & Fbad,
        const int & node_bad, // indices into Vbad
        const std::vector<int> & node_list_good, // nodes to exclude
        Eigen::MatrixXi & TT_good, // connectivity info
        std::vector<std::vector<int> > & VV_good,// connectivity info
        std::vector<std::vector<int> > & TEdges_good,
        std::vector<std::vector<int> > & VEdges_good, //edge vertices to exclude
        Eigen::MatrixXd & Vgood,
        Eigen::MatrixXi & Fgood,
        Eigen::VectorXi & FL_good,
        int & node_image
    );


    void proj_node(
        const CellularGraph & cg,
        const std::vector<int> & node_list_bad,
        Eigen::MatrixXd & Vgood,
        Eigen::MatrixXi & Fgood,
        std::map<int, int> & node_map
    );
    // try to project the nodes into good mesh and modify the good mesh correspondingly
    // if it can not find intersection use nearest neighbor instead 
    // always gaurantee that node_map is injective
    void proj_node_loop(
        const CellularGraph & cg,
        const int & node_bad, // indices into Vbad
        const std::vector<int> & node_list_good, // nodes to exclude
        Eigen::MatrixXi & TT_good, // connectivity info
        std::vector<std::vector<int> > & VV_good,// connectivity info
        std::vector<std::vector<int> > & TEdges_good,
        std::vector<std::vector<int> > & VEdges_good, //edge vertices to exclude
        Eigen::MatrixXd & Vgood,
        Eigen::MatrixXi & Fgood,
        Eigen::VectorXi & FL_good,
        int & node_image
    );

    void proj_node_loop(
        const CellularGraph & cg,
        const int & node_bad,
        MatchMaker::TraceComplex & tc,
        int & node_image
    );

    
}
#endif //BCCLEAN_PROJ_NODE_H