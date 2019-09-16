#ifndef BCCLEAN_EDGE_H
#define BCCLEAN_EDGE_H
#include <Eigen/Core>
#include <vector>
#include <map>
#include <utility>
#include <unordered_map>
#include "node.h"

namespace bcclean{
    struct pairhash {
        public:
        template <typename T, typename U>
        std::size_t operator()(const std::pair<T, U> &x) const
        {
            return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
        }
    };
    template<class U, class T>
    using pair_map = std::unordered_map<U, T, pairhash>;

    enum class Edge_Compare_Result {IDENTICAL, SAME_TYPE, NON_MATCH};
    class edge {
        public:
        int head;
        int tail;
        const int total_label_num;
        std::pair<int, int> _label_pair; // lab1, lab2, wheere lab1 < lab2;
        std::vector<int> _edge_vertices;
        bool matched;
        Edge_Compare_Result compare_edge(edge b);
    };
    bool build_edge_dict(const Eigen::MatrixXd& V, Eigen::MatrixXi&F, const  Eigen::VectorXi & FL, pair_map<std::pair<int,int>, std::vector<edge> > & edge_dict);
}

#endif // BCCLEAN_EDGE_H