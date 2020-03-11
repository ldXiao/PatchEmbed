#include "orientation_check.h"
#include "kdtree_NN_Eigen.hpp"
#include <igl/per_face_normals.h>
#include <igl/embree/line_mesh_intersection.h>
#include <igl/bfs_orient.h>
namespace bcclean {
    bool orientation_check(
        const Eigen::MatrixXd & VA,
        const Eigen::MatrixXi & FA,
        const Eigen::MatrixXd & VB,
        const Eigen::MatrixXi & FB
    )
    {
        Eigen::MatrixXd BaryCenterA = Eigen::MatrixXd::Constant(FA.rows(), 3, 0);
        Eigen::MatrixXd BaryCenterB = Eigen::MatrixXd::Constant(FB.rows(), 3, 0);
        for(int j =0; j< FA.rows();++j){
            int ii, jj, kk;
            ii = FA(j,0);
            jj = FA(j,1);
            kk = FA(j,2);
            BaryCenterA.row(j)= (VA.row(ii)+ VA.row(jj)+ VA.row(kk))/3;
        }
        for(int j =0; j< FB.rows();++j){
            int ii, jj, kk;
            ii = FB(j,0);
            jj = FB(j,1);
            kk = FB(j,2);
            BaryCenterB.row(j)= (VB.row(ii)+ VB.row(jj)+ VB.row(kk))/3;
        }

        Eigen::MatrixXd NA,NB, NAB;
        igl::per_face_normals(VA, FA, NA);
        igl::per_face_normals(VB,FB, NB);
        kd_tree_Eigen<double> kdt(BaryCenterA.cols(),std::cref(BaryCenterA),10);
        kdt.index->buildIndex();
        NAB = Eigen::MatrixXd::Constant(FB.rows(),3, 0);
        for(int fidx = 0; fidx< FB.rows(); fidx++)
        {
            Eigen::RowVector3d query = BaryCenterB.row(fidx);
            int fjdxA =kd_tree_NN_Eigen(kdt, query);
            if(fjdxA != -1)
            {
                NAB.row(fidx)= NA.row(fjdxA);
            }
        }
        double vote = 0;
        for(int fidx =0 ; fidx< FB.rows(); fidx++)
        {
            double contrib = NAB.row(fidx).dot(NB.row(fidx));
            vote += contrib;
        }
        return (vote > 0);
    }
    void flip_orientation_ifnecessary
    (
        const Eigen::MatrixXd & VA,
        const Eigen::MatrixXi & FA,
        const Eigen::MatrixXd & VB,
        Eigen::MatrixXi & FB
    )
    {
        if(! orientation_check(VA, FA, VB, FB))
        {
            for(int fidx =0 ; fidx < FB.rows(); fidx++)
            {
                int v1,v2;
                v1 = FB(fidx,1);
                v2 = FB(fidx, 2);
                FB(fidx,1)=v2;
                FB(fidx,2)=v1;
            }
        }
        return;
    }
}