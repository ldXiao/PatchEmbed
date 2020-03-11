#include "backtrack_diff.h"
namespace bcclean{
namespace MatchMaker{
    double path_len(const Eigen::MatrixXd & V, 
                    const std::vector<int> & path
    )
    {
        double res=0;
        for(int p=0; p < path.size()-1; p++)
        {
            double edgeln = (V.row(path[p])- V.row(path[p+1])).norm();
            res += edgeln;
        }
        return res;
    }
    bool backtrack_diff(
        const Eigen::MatrixXd & V_good,
        const Eigen::MatrixXd & V_bad,
        const int pidx,
        const std::unordered_map<int, std::vector<int> > & patch_edge_dict,
        const std::vector<edge> & edge_list,
        const std::unordered_map<int, std::vector<bool> > & patch_edge_direction_dict,
        const std::map<int, std::vector<int> > & edge_path_map,
        const double threshold
    )
    {
        // the first step is to check wether the  orientation has been inverted.
        // postone 
        // the second step is to measure the max length difference
        std::vector<int> edge_indices = patch_edge_dict.at(pidx);
        double max_rel_dif= -1;
        double total_len_good =0;
        double total_len_bad = 0;
        for(auto edg_idx: edge_indices)
        {
            std::vector<int> path_good = edge_path_map.at(edg_idx);
            std::vector<int> path_bad = edge_list.at(edg_idx)._edge_vertices;
            total_len_good += path_len(V_good, path_good);
            total_len_bad += path_len(V_bad, path_bad);
        }
         max_rel_dif = std::max(max_rel_dif, std::abs(total_len_good/total_len_bad -1));
        if(max_rel_dif > threshold){
            return false;
        }
        return true;
    }

    bool backtrack_diff(
        const Eigen::MatrixXd & V_good,
        const CellularGraph & cg,
        const int pidx,
        const std::map<int, std::vector<int> > & edge_path_map,
        const double threshold
    ){
        // the first step is to check wether the  orientation has been inverted.
        // postone 
        // the second step is to measure the max length difference
        std::vector<int> edge_indices = cg._patch_edge_dict.at(pidx);
        double max_rel_dif= -1;
        double total_len_good =0;
        double total_len_bad = 0;
        Eigen::MatrixXd V_bad = Eigen::MatrixXd::Constant(cg._vertices.size(),3, -1);
        for(int vidx= 0; vidx< V_bad.rows(); ++vidx)
        {
            V_bad.row(vidx) = cg._vertices.at(vidx);
        }
        for(auto edg_idx: edge_indices)
        {
            std::vector<int> path_good = edge_path_map.at(edg_idx);
            std::vector<int> path_bad = cg._edge_list.at(edg_idx)._edge_vertices;
            total_len_good += path_len(V_good, path_good);
            total_len_bad += path_len(V_bad, path_bad);
        }
         max_rel_dif = std::max(max_rel_dif, std::abs(total_len_good/total_len_bad -1));
        if(max_rel_dif > threshold){
            return false;
        }
        return true;
    }
}
}