#ifndef BCCLEAN_NODE_H
#define BCCLEAN_NODE_H
#include <Eigen/Core>
#include <vector>
#include <map>
namespace bcclean{
    class node {
        public:
        Eigen::MatrixXd _position;
        int _total_label_num;
        int _occupied_label_num;
        std::map<int, int> _label_occupy_dict;
        bool of_same_type(const node & b);
        bool at_same_position(const Eigen::MatrixXd & position);
        bool initialize(const int total_label_num, const Eigen::MatrixXd & position, const std::vector<int> labels);
    };
}
#endif //BCCLEAN_NODE_H