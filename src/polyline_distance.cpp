#include "polyline_distance.h"
#include "kdtree_NN_Eigen.hpp"
namespace bcclean{
namespace Eval{
    double hausdorff1d(
        const Eigen::MatrixXd & VA,
        const std::vector<int> & pathA,
        const Eigen::MatrixXd & VB,
        const std::vector<int> & pathB
    ){
        // increase the density by 30
        int upratio = 30;
        int target_numA = upratio * (pathA.size()-1) + 1;
        int target_numB = upratio * (pathB.size()-1) +  1;
        Eigen::MatrixXd sampleA(target_numA, 3);
        Eigen::MatrixXd sampleB(target_numB, 3);
        for(int jj =0 ; jj < target_numA ; ++jj)
        {
            int batch = (int)(jj / upratio);
            int remain = jj % upratio;
            if(remain == 0) 
            {
                sampleA.row(jj) = VA.row(pathA[batch]);
            }
            else
            {
                double lambda = (double)(remain) / (double)(upratio);
                sampleA.row(jj) = (1-lambda)* VA.row(pathA[batch]) + lambda* VA.row(pathA[batch+1]);
            }
        }
        kd_tree_Eigen<double> sampleA_kdt(sampleA.cols(),std::cref(sampleA),10);
        sampleA_kdt.index->buildIndex();
        for(int jj =0 ; jj < target_numB ; ++jj)
        {
            int batch = (int)(jj / upratio);
            int remain = jj % upratio;
            if(remain == 0) 
            {
                sampleB.row(jj) = VB.row(pathB[batch]);
            }
            else
            {
                double lambda = (double)(remain) / (double)(upratio);
                sampleB.row(jj) = (1-lambda)* VB.row(pathB[batch]) + lambda* VB.row(pathB[batch+1]);
            }
        }
        kd_tree_Eigen<double> sampleB_kdt(sampleB.cols(),std::cref(sampleB),10);
        sampleB_kdt.index->buildIndex();
        double maxminAB= -1;
        double maxminBA= -1;
        for(int aj =0 ; aj < sampleA.rows(); aj++)
        {
            Eigen::RowVector3d jqueryA= sampleA.row(aj);
            int sampleB_idx= kd_tree_NN_Eigen(sampleB_kdt, jqueryA); 
            
            maxminAB = std::max((jqueryA - sampleB.row(sampleB_idx)).norm(), maxminAB);
        }
        for(int bj =0 ; bj < sampleB.rows(); bj++)
        {
            Eigen::RowVector3d jqueryB= sampleB.row(bj);
            int sampleA_idx= kd_tree_NN_Eigen(sampleA_kdt, jqueryB); 
            maxminBA = std::max((jqueryB - sampleA.row(sampleA_idx)).norm(), maxminBA);
        }
        return std::max(maxminAB, maxminBA);
        
    }

    double haursdorff1d(
        const std::vecotr<Eigen::RowVector3d> & VA,
        const std::vector<int> & pathA,
        const Eigen::MatrixXd & VB,
        const std::vector<int> & pathB
    )
    {
        // increase the density by 30
        int upratio = 30;
        int target_numA = upratio * (pathA.size()-1) + 1;
        int target_numB = upratio * (pathB.size()-1) +  1;
        Eigen::MatrixXd sampleA(target_numA, 3);
        Eigen::MatrixXd sampleB(target_numB, 3);
        for(int jj =0 ; jj < target_numA ; ++jj)
        {
            int batch = (int)(jj / upratio);
            int remain = jj % upratio;
            if(remain == 0) 
            {
                sampleA.row(jj) = VA[pathA[batch]];
            }
            else
            {
                double lambda = (double)(remain) / (double)(upratio);
                sampleA.row(jj) = (1-lambda)* VA[pathA[batch]] + lambda* VA[pathA[batch+1]];
            }
        }
        kd_tree_Eigen<double> sampleA_kdt(sampleA.cols(),std::cref(sampleA),10);
        sampleA_kdt.index->buildIndex();
        for(int jj =0 ; jj < target_numB ; ++jj)
        {
            int batch = (int)(jj / upratio);
            int remain = jj % upratio;
            if(remain == 0) 
            {
                sampleB.row(jj) = VB.row(pathB[batch]);
            }
            else
            {
                double lambda = (double)(remain) / (double)(upratio);
                sampleB.row(jj) = (1-lambda)* VB.row(pathB[batch]) + lambda* VB.row(pathB[batch+1]);
            }
        }
        kd_tree_Eigen<double> sampleB_kdt(sampleB.cols(),std::cref(sampleB),10);
        sampleB_kdt.index->buildIndex();
        double maxminAB= -1;
        double maxminBA= -1;
        for(int aj =0 ; aj < sampleA.rows(); aj++)
        {
            Eigen::RowVector3d jqueryA= sampleA.row(aj);
            int sampleB_idx= kd_tree_NN_Eigen(sampleB_kdt, jqueryA); 
            
            maxminAB = std::max((jqueryA - sampleB.row(sampleB_idx)).norm(), maxminAB);
        }
        for(int bj =0 ; bj < sampleB.rows(); bj++)
        {
            Eigen::RowVector3d jqueryB= sampleB.row(bj);
            int sampleA_idx= kd_tree_NN_Eigen(sampleA_kdt, jqueryB); 
            maxminBA = std::max((jqueryB - sampleA.row(sampleA_idx)).norm(), maxminBA);
        }
        return std::max(maxminAB, maxminBA);
    }

    double single_sample_trial(
        const Eigen::MatrixXd & VA,
        const std::vector<int> & pathA,
        const Eigen::MatrixXd & VB,
        const int & idxB
    )
    {
        // increase the density by 30
        int upratio = 30;
        int target_numA = upratio * (pathA.size()-1) + 1;
        Eigen::MatrixXd sampleA(target_numA, 3);
        for(int jj =0 ; jj < target_numA ; ++jj)
        {
            int batch = (int)(jj / upratio);
            int remain = jj % upratio;
            if(remain == 0) 
            {
                sampleA.row(jj) = VA.row(pathA[batch]);
            }
            else
            {
                double lambda = (double)(remain) / (double)(upratio);
                sampleA.row(jj) = (1-lambda)* VA.row(pathA[batch]) + lambda* VA.row(pathA[batch+1]);
            }
        }
        kd_tree_Eigen<double> sampleA_kdt(sampleA.cols(),std::cref(sampleA),10);
        Eigen::RowVector3d jqueryB= VB.row(idxB);
        int sampleA_idx= kd_tree_NN_Eigen(sampleA_kdt, jqueryB); 
        return (jqueryB- sampleA.row(sampleA_idx)).norm();
        
    }
}
}