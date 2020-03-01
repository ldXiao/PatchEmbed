#ifndef BCCLEAN_RECURCUT_H
#define BCCLEAN_RECURCUT_H
#include <Eigen/Dense>
#include <map>
namespace bcclean{
namespace Preprocess{
    void recursive_cut(
        Eigen::MatrixXd & Vbase,
        Eigen::MatrixXi & Fbase, 
        Eigen::VectorXi & FLbase, // face label on Fbas
        std::map<int, Eigen::MatrixXi> & label_faces_dict,
        std::map<int, Eigen::VectorXi> & label_FI_dict,
        const int lb
    );
}
}
#endif //