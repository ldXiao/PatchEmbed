#include "patch_cut_relabel.h"
namespace bcclean{
    void patch_cut_relabel(const Eigen::MatrixXi & Fraw, const Eigen::VectorXi FI, const std::vector<bool> & VCuts, const Eigen::VectorXi & FL, Eigen::VectorXi & FL_mod, int & total_label_num){
        assert(total_label_num == FL.maxCoeff()+1);
        FL_mod = FL;
        Eigen::VectorXi PatchL;

    }
}