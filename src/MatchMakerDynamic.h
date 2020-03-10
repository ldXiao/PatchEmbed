#ifndef BCCLEAN_MATCHMAKERDYNAMIC_H
#define BCCLEAN_MATCHMAKERDYNAMIC_H
#include "CellularGraph.h"
#include "params.h"
#include <Eigen/Dense>
namespace bcclean{
namespace MatchMaker{
    bool MatchMakerDynamic(
        const CellularGraph & cg,
        const params & param,
        Eigen::MatrixXd & V_good,
        Eigen::MatrixXi & F_good,
        Eigen::VectorXi & FL_good
    );
}
}


#endif// 