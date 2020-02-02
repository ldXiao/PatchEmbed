#include <Eigen/Dense>

namespace bcclean{
    bool orientation_check(
        const Eigen::MatrixXd & VA,
        const Eigen::MatrixXi & FA,
        const Eigen::MatrixXd & VB,
        const Eigen::MatrixXi & FB
        );
    
    void flip_orientation_ifnecessary
    (
        const Eigen::MatrixXd & VA,
        const Eigen::MatrixXi & FA,
        const Eigen::MatrixXd & VB,
        Eigen::MatrixXi & FB
    );
     
}