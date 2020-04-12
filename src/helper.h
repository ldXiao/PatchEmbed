#pragma once
#include <Eigen/Dense>
#include <vector>
namespace bcclean{
namespace Helper{
    void to_list(const Eigen::MatrixX3i & M, std::vector<Eigen::RowVector3i> & L);    
    void to_list(const Eigen::MatrixX3d & M, std::vector<Eigen::RowVector3d> & L);
    void to_matrix(const std::vector<Eigen::RowVector3i> & L, Eigen::MatrixXi & M);
    void to_matrix(const std::vector<Eigen::RowVector3d> & L, Eigen::MatrixXd & M);

}
}