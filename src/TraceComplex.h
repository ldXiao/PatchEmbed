#ifndef BCCLEAN_TRACECOMPLEX_H
#define BCCLEAN_TRACECOMPLEX_H
#include <Eigen/Dense>
#include <vector>
#include <map>
#include <utility>
namespace bcclean{
namespace MatchMaker{
    class TraceComplex{
        public:
        std::vector<std::vector<int> > _VV;
        std::vector<std::vector<int> > _VF;
        std::vector<std::vector<int> > _VEdges;
        std::vector<std::vector<int> > _TEdges;
        Eigen::MatrixXi _TT;
        Eigen::MatrixXi _F;
        Eigen::MatrixXd _V;
        Eigen::VectorXi _FL;
        std::map<int, std::vector<int> > _edge_path_map;
        std::map<int, int> _node_image_map;
        std::vector<int> _total_silence_list;
        std::vector<int> _node_list;
        void initialize(const Eigen::MatrixXd & Vin, const Eigen::MatrixXi & Fin);
        void insert_update(const Eigen::MatrixXd & baryentry);
        bool split_detect(std::pair<int,int> & split);
        void split_update(const std::pair<int, int> & split);
    
    };
}
}
#endif