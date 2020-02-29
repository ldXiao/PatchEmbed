#ifndef BCCLEAN_NV_MANIFOLD_RELABEL_H
#define BCCLEAN_NV_MANIFOLD_RELABEL_H
#include <Eigen/Dense>
#include <vector>
#include <map>
namespace bcclean{
namespace Prepocess{
    void non_vertex_manifold_relabel(
        const Eigen::MatrixXd & Vraw,
        const Eigen::MatrixXi & Fraw, 
        const Eigen::VectorXi & FI,
        const std::vector<int> & NMV,
        const Eigen::VectorXi & FL, 
        const int cur_lb,
        Eigen::VectorXi & FL_mod, 
        int & total_label_num,
        std::map<int, Eigen::VectorXi> & subpatches
    );
}
}
#endif 