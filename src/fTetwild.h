#ifndef CUSTOM_FTETWILD_H
#define CUSTOM_FTETWILD_H


#include <Eigen/Dense>
namespace bcclean{
    namespace Tet{
        int fTetwild(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F,const double edge_len_r,const int stop_eng,Eigen::MatrixXd & VS, Eigen::MatrixXi & FS);
    }
}
#endif