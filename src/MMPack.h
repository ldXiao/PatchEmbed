#ifndef BCCLEAN_MMPACK_H
#define BCCLEAN_MMPACK_H
#include <vector>
#include <Eigen/Dense>
namespace bcclean {
    class MMPack{
        public:
            Eigen::MatrixXd V;
            Eigen::MatrixXi F;
            Eigen::VectorXi FL;
            Eigen::MatrixXi TT;
            std::vector<std::vector<int> > VV;
            std::vector<std::vector<int> > VF;
            std::vector<std::vector<int> > VEdges;
            std::vector<std::vector<int> > TEdges;
    };
}
#endif // BCCLEAN_MMPACK_H