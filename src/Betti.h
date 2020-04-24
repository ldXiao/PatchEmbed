#ifndef BCCLEAN_BETTI_H
#define BCCLEAN_BETTI_H
#include <igl/edges.h>
namespace bcclean{
int Betti(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F);
}
#endif