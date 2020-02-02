#include "orientation_check.h"
#include <igl/per_face_normals.h>
#include <igl/embree/line_mesh_intersection.h>
namespace bcclean {
    bool orientation_check(
        const Eigen::MatrixXd & VA,
        const Eigen::MatrixXi & FA,
        const Eigen::MatrixXd & VB,
        const Eigen::MatrixXi & FB
    )
    {
        Eigen::MatrixXd BaryCenterA = Eigen::MatrixXd::Constant(FA.rows(), 3, 0);
        for(int j =0; j< FA.rows();++j){
            int ii, jj, kk;
            ii = FA(j,0);
            jj = FA(j,1);
            kk = FA(j,2);
            BaryCenterA.row(j)= (VA.row(ii)+ VA.row(jj)+ VA.row(kk))/3;
        }
        Eigen::MatrixXd NA,NB, NAB;
        igl::per_face_normals(VA, FA, NA);
        igl::per_face_normals(VB,FB, NB);
        NAB = Eigen::MatrixXd::Constant(FA.rows(),3, 0);
        Eigen::MatrixXd RR = igl::embree::line_mesh_intersection(BaryCenterA, NA, VB, FB);
        for(int fidx = 0; fidx< FA.rows(); fidx++)
        {
            int fjdxB = std::round(RR(fidx, 0));
            if(fjdxB != -1)
            {
                NAB.row(fidx)= NB.row(fjdxB);
            }
        }
        double vote = 0;
        for(int fidx =0 ; fidx< FA.rows(); fidx++)
        {
            double contrib = NAB.row(fidx).dot(NA.row(fidx));
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