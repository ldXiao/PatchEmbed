#ifndef BCCLEAN_NV_MANIFOLD_RELABEL_H
#define BCCLEAN_NV_MANIFOLF_RELABEL_H
#include <Eigen/Dense>
#include <vector>
namespace bcclean{
namespace Prepocess{
    void non_vertex_manifold_relabel(
        const Eigen::MatrixXd & Vraw,
        const Eigen::MatrixXi & Fraw, 
        const Eigen::VectorXi & FI,
        const std::vector<int> & NMV,
        const Eigen::VectorXi & FL, 
        Eigen::VectorXi & FL_mod, 
        int & total_label_num
    );
}
}
#endif 