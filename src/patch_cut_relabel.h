#include <Eigen/Core>
#include <vector>
#include <map>
namespace bcclean{
    void patch_cut_relabel(const Eigen::MatrixXi & Fraw, const Eigen::VectorXi & FI, const std::vector<bool> & VCuts, const std::vector<std::vector<bool> > & TCuts, const Eigen::VectorXi & FL, Eigen::VectorXi & FL_mod, int & total_label_num);
    void patch_cut_relabel(
        const Eigen::MatrixXi & Fraw, 
        const Eigen::VectorXi & FI, 
        const std::vector<bool> & VCuts, 
        const std::vector<std::vector<bool> > & TCuts, 
        const Eigen::VectorXi & FL, 
        Eigen::VectorXi & FL_mod, 
        std::map<int, Eigen::VectorXi> & subpatchFI_dict);
}