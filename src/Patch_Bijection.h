#include <Eigen/Dense>
#include <igl/harmonic.h>
#include "patch.h"
#include "CellularGraph.h"
#include <vector>
namespace bcclean {
namespace Bijection{
   void mapping2polygon(
        const CellularGraph & cg,
        const int pidx,
        Eigen::MatrixXd & V_uv,
        Eigen::MatrixXi & F_uv,
        std::vector<int> nodesp_uv,
        Eigen::VectorXi & VI, 
        Eigen::VectorXi & FI
    );

    void BijGlobal(
        const CellularGraph & cga,
        const CellularGraph & cgb,
        Eigen::MatrixXd & M_a2b // #Va * 4 matrix
    );
}
}