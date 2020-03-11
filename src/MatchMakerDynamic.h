#ifndef BCCLEAN_MATCHMAKERDYNAMIC_H
#define BCCLEAN_MATCHMAKERDYNAMIC_H
/* this BTCMM is the abbreviation of  backtracking cellular match maker */
#include "Match_Maker_Tree.h"
#include "proj_node.h"
#include "CellularGraph.h"
#include "params.h"
#include <igl/slice.h>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#include <utility>
#include <nlohmann/json.hpp>
#include <Eigen/Core>
namespace bcclean{
namespace MatchMaker{
    bool BTCMM1(
        const CellularGraph & cg,
        Eigen::MatrixXd & V_good,
        Eigen::MatrixXi & F_good,
        Eigen::VectorXi & FL_good,
        const params params
    );
}
}


#endif// 