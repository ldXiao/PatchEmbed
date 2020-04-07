//
// Created by Zhongshi Jiang on 5/2/17.
//

#ifndef SCAFFOLD_TEST_TETRAHEDRAL_REFINEMENT_H
#define SCAFFOLD_TEST_TETRAHEDRAL_REFINEMENT_H

#include <vector>
#include <igl/igl_inline.h>
#include <functional>
#include <set>
#include <Eigen/Core>

namespace igl{ namespace dev {
bool tet_tuple_edge_split(int ti, int fi, int ei, bool ai,
                          std::vector<Eigen::RowVector3d> &V,
                          std::vector<Eigen::RowVector4i> &T,
                          std::vector<Eigen::RowVector4i> &TT,
                          std::vector<Eigen::RowVector4i> &TTif,
                          std::vector<Eigen::Matrix<int, 4, 3>> &TTie,
                          std::vector<int> &new_tets_id);
}
}
#ifndef IGL_STATIC_LIBRARY
#include "tet_split.cpp"
#endif
#endif //SCAFFOLD_TEST_TETRAHEDRAL_REFINEMENT_H
