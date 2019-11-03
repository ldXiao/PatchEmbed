#ifndef BCCLEAN_EDGE_H
#define BCCLEAN_EDGE_H
#include <Eigen/Core>
#include <vector>
#include <map>
#include <utility>
#include <unordered_map>
#include "node.h"
#include <igl/opengl/glfw/Viewer.h>

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
        int total_label_num;
        std::pair<int, int> _label_pair; // lab1, lab2, where lab1 < lab2;
        std::vector<int> _edge_vertices;
        bool matched = true;
        Edge_Compare_Result _compare_edge(edge b);
    };
    bool build_vertex_label_list_dict(const Eigen::MatrixXi&F, const Eigen::MatrixXi & FL, const size_t total_label_num, std::unordered_map<int, std::vector<int>> & vertex_label_list_dict);
    bool build_edge_list(const Eigen::MatrixXd& V, const Eigen::MatrixXi&F, const Eigen::MatrixXi & FL, const size_t total_label_num, std::vector<edge> & edge_list, std::unordered_map<int, std::vector<int> > & patch_edge_dict);
    void plot_edge(igl::opengl::glfw::Viewer & viewer, const Eigen::MatrixXd & V, const Eigen::VectorXi &FL, const edge & edg);

}

#endif // BCCLEAN_EDGE_H