#include "helper.h"
namespace bcclean{
namespace Helper{
    void to_list(const Eigen::MatrixX3i & M, std::vector<Eigen::RowVector3i> & L)
    {
        L.resize(M.rows());
        for(int fidx =0 ; fidx < M.rows(); ++fidx)
        {
            L[fidx] = M.row(fidx);
        }
        return;
    }
    void to_list(const Eigen::MatrixX3d & M, std::vector<Eigen::RowVector3d> & L)
    {
        L.resize(M.rows());
        for(int vidx =0; vidx< M.rows(); ++vidx)
        {
            L[vidx]= M.row(vidx);
        }
        return;
    }

    void to_matrix(const std::vector<Eigen::RowVector3i> & L, Eigen::MatrixXi & M)
    {
        M = Eigen::MatrixXi::Constant(L.size(), 3, 0);
        for(int fidx =0; fidx < L.size(); ++fidx)
        {
            M.row(fidx) = L[fidx];
        }
        return;
    }
    void to_matrix(const std::vector<Eigen::RowVector3d> & L, Eigen::MatrixXd & M)
    {
        M = Eigen::MatrixXd::Constant(L.size(), 3, 0);
        for(int vidx=0; vidx < L.size();++vidx)
        {
            M.row(vidx) = L[vidx];
        }
        return;
    }
}
}