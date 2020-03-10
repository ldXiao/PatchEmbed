#ifndef BCCLEAN_TRACECOMPLEX_H
#define BCCLEAN_TRACECOMPLEX_H
#include <Eigen/Dense>
#include <vector>
#include <map>
namespace bcclean{
namespace MatchMaker{
    class TraceComplex{
        public:
        std::vector<std::vector<int> > VV;
        std::vector<std::vector<int> > VF;
        std::vector<std::vector<int> > VEdges;
        std::vector<std::vector<int> > TEdges;
        Eigen::MatrixXi TT;
        Eigen::MatrixXi F;
        Eigen::MatrixXd V;
        Eigen::VectorXi FL;
        std::map<int, std::vector<int> > edge_path_map;
        std::map<int, int> node_image_map;
    };
}
}
#endif