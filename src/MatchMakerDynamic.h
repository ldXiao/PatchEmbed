#ifndef BCCLEAN_MATCHMAKERDYNAMIC_H
#define BCCLEAN_MATCHMAKERDYNAMIC_H
/* this BTCMM is the abbreviation of  backtracking cellular match maker */
#include "MatchMakerTree.h"
#include "proj_node.h"
#include "CellularGraph.h"
#include "params.h"
#include <igl/slice.h>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#include <utility>
#include <nlohmann/json.hpp>
#include <Eigen/Core>
#include <spdlog/spdlog.h>
#include <memory>
namespace bcclean{
namespace MatchMaker{
    bool MatchMakerPatch(
        const CellularGraph & cg,
        Eigen::MatrixXd & V_good,
        Eigen::MatrixXi & F_good,
        Eigen::VectorXi & FL_good,
        const params params,
        std::shared_ptr<spdlog::logger> logger
    );
}
}


#endif// 