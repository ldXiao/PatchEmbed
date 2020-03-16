#include "Match_Maker_Tree.h"
#include "proj_node.h"
#include "params.h"
#include <igl/slice.h>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#include <utility>
#include <nlohmann/json.hpp>
#include <Eigen/Core>
namespace bcclean{
namespace MatchMaker{
    bool trace_and_label_loop(
        const Eigen::MatrixXd & V_bad,
        const Eigen::MatrixXi & F_bad,
        const Eigen::VectorXi & FL_bad,
        Eigen::MatrixXd & V_good,
        Eigen::MatrixXi & F_good,
        Eigen::VectorXi & FL_good,
        const params params
    );
}
}