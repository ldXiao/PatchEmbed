#include <Eigen/Dense>
#include <igl/barycenter.h>
#include <igl/writeDMAT.h>
#include "../kdtree_NN_Eigen.hpp"

void nn_transfer(const Eigen::MatrixXd & Vs,
                const Eigen::MatrixXi & Fs,
                const Eigen::VectorXi & FLs,
                const Eigen::MatrixXd & Vt,
                const Eigen::MatrixXi & Ft,
                Eigen::VectorXi& FLt,
                Eigen::VectorXi &VL)
{
    Eigen::MatrixXd BCs, BCt;
    igl::barycenter(Vs, Fs, BCs);
    igl::barycenter(Vt, Ft, BCt);
    bcclean::kd_tree_Eigen<double> KDTs(BCs.cols(),std::cref(BCs),10);
    FLt.resize(Ft.rows());
    KDTs.index->buildIndex();
    for(int fidx = 0 ; fidx < Ft.rows(); fidx++)
    {
        Eigen::RowVector3d query = BCt.row(fidx) ;
        int sfidx= bcclean::kd_tree_NN_Eigen(KDTs, query);
        FLt(fidx) = FLs(sfidx);
    }
    VL.resize(Vt.rows());
    for(int vidx=0; vidx < Vt.rows(); ++vidx)
    {
        Eigen::RowVector3d query = Vt.row(vidx);
        int sfidx = bcclean::kd_tree_NN_Eigen(KDTs,query);
        VL(vidx)= FLs(sfidx);
    }
    
    // igl::writeOBJ("../../blenders/nnfl.dmat",FLt);
    return;
}