#ifndef BCCLEAN_DEGENERATE_H
#define BCCLEAN_DEGENERATE_H
#include <Eigen/Dense>
namespace bcclean {
    void degenerate_clean(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        Eigen::VectorXi & FL,
        const double threshold);
}
#endif