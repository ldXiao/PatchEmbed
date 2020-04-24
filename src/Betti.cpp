#include "Betti.h"
namespace bcclean {
int Betti(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F)
    {
        Eigen::MatrixXi E;
        igl::edges(F,E);
        return V.rows() - E.rows() + F.rows();
    }
}