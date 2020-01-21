#ifndef BCCLEAN_HAUSDORFF_H
#define BCCLEAN_HAUSDORFF_H
#include <Eigen/Dense>
namespace bcclean{
namespace Eval{
        double hausdorff1d(
            const Eigen::MatrixXd & VA,
            const std::vector<int> & pathA,
            const Eigen::MatrixXd & VB,
            const std::vector<int> & pathB
        );  

        
        double single_sample_trial(
            const Eigen::MatrixXd & VA,
            const std::vector<int> & pathA,
            const Eigen::MatrixXd & VB,
            const int & idxB
        );   
}
}
#endif