#ifndef BCCLEAN_BTDIFF_H
#define BCCLEAN_BTDIFF_H
#include <Eigen/Dense>
#include <map>
#include <unordered_map>
#include <vector>
#include "edge.h"
namespace bcclean{
namespace MatchMaker{
    bool backtrack_diff(
        const Eigen::MatrixXd & V_good,
        const Eigen::MatrixXd & V_bad,
        const int pidx,
        const std::unordered_map<int, std::vector<int> > & patch_edge_dict,
        const std::vector<edge> & edge_list,
        const std::unordered_map<int, std::vector<bool> > & patch_edge_direction_dict,
        const std::map<int, std::vector<int> > & edge_path_map,
        const double threshold
    );
}
}
#endif 