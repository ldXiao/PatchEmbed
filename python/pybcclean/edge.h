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
        int total_label_num;
        std::pair<int, int> _label_pair; // lab1, lab2, where lab1 < lab2;
        std::vector<int> _edge_vertices;
        bool matched = true;
        Edge_Compare_Result _compare_edge(edge b);
    };
    bool build_edge_dict(const Eigen::MatrixXd& V, Eigen::MatrixXi&F, const  Eigen::VectorXi & FL, const  size_t total_label_num, pair_map<std::pair<int,int>, std::vector<edge> > & edge_dict);
    bool build_vertex_label_list_dict(const Eigen::MatrixXi&F, const Eigen::MatrixXi & FL, const size_t total_label_num, std::unordered_map<int, std::vector<int>> & vertex_label_list_dict);
    bool build_edge_list(const Eigen::MatrixXd& V, const Eigen::MatrixXi&F, const Eigen::MatrixXi & FL, const size_t total_label_num, std::vector<edge> & edge_list, std::unordered_map<int, std::vector<int> > & patch_edge_dict);
    bool simple_match_two_edges(const edge & edg0, const edge & edg1);
    void match_patch_edges(const int & lb0, const int & lb1, const std::unordered_map<int, std::vector<int> > & patch_edge,  std::vector<edge> & edge_list);
    void build_pair_edge_list(const std::vector<edge> & edge_list, pair_map<std::pair<int,int>,std::vector<int>> & pair_edge_list_dict);
    void match_pair_edge_list_dicts(
        const pair_map<std::pair<int,int>, std::vector<int> > & pair_edge_list_dict0, 
        const pair_map<std::pair<int,int>, std::vector<int> > & pair_edge_list_dict1 ,
        std::vector<edge> & edge_list0,
        std::vector<edge> & edge_list1
    );
}

#endif // BCCLEAN_EDGE_H