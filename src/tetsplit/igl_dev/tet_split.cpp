//
// Created by Zhongshi Jiang on 5/2/17.
//

#include "tet_split.h"

#include <spdlog/spdlog.h>

#include "retain_tetrahedral_adjacency.h"
#include "tetrahedron_tetrahedron_adjacency.h"
#include "tetrahedron_tuple.h"

namespace igl::dev {

bool tet_tuple_edge_split(int ti, int fi, int ei, bool ai,
                          std::vector<Eigen::RowVector3d> &V,
                          std::vector<Eigen::RowVector4i> &T,
                          std::vector<Eigen::RowVector4i> &TT,
                          std::vector<Eigen::RowVector4i> &TTif,
                          std::vector<Eigen::Matrix<int, 4, 3>> &TTie,
                          std::vector<int> &new_tets_id) {
  ai = true;
  int vert_a = igl::dev::tet_tuple_get_vert(ti, fi, ei, ai, T, TT, TTif, TTie);
  int vert_b = igl::dev::tet_tuple_get_vert(ti, fi, ei, !ai, T, TT, TTif, TTie);

  // get the edge pointing to another vertex.
  // swtich e/f
  // igl::dev::tet_tuple_switch_edge(ti, fi, ei, ai, T, TT, TTif, TTie);
  // igl::dev::tet_tuple_switch_face(ti, fi, ei, ai, T, TT, TTif, TTie);

  // loop through the one ring of tets around an edge.
  int f = fi, t = ti, e = ei;
  spdlog::trace("t {} f {} e {}", t, f, e);
  spdlog::trace("va {} vb {}", vert_a, vert_b);
  std::set<int> one_ring_tets;

  do {
    one_ring_tets.insert(t);
    // swtch f/t/f/e
    igl::dev::tet_tuple_switch_face(t, f, e, ai, T, TT, TTif, TTie);
    if (igl::dev::tet_tuple_is_on_boundary(t, f, e, ai, T, TT, TTif, TTie))break;
    igl::dev::tet_tuple_switch_tet(t, f, e, ai, T, TT, TTif, TTie);
    spdlog::trace("t {} f {} e {}", t, f, e);
  } while (t != ti);

  // attempt addition
  V.push_back((V[vert_a] + V[vert_b]) / 2);
  int vert_m = V.size() - 1;

  std::vector<Eigen::RowVector4i> new_tets;
  new_tets.reserve(2 * one_ring_tets.size());
  for (auto old_v : {vert_a, vert_b})
    for (auto s : one_ring_tets) {
      auto ts = T[s];
      for (auto j : {0, 1, 2, 3})
        if (ts[j] == old_v) {
          ts[j] = vert_m;
          break;
        }
      new_tets.push_back(ts);
    }

  std::set<int> influence_id;
  for (auto t : one_ring_tets)
    for (auto f : {0, 1, 2, 3}) influence_id.insert(TT[t][f]);
  influence_id.erase(-1);

  auto &delete_tets = one_ring_tets;
  new_tets_id.clear();
  new_tets_id.insert(new_tets_id.end(), delete_tets.begin(), delete_tets.end());
  for (int i = 0, loop_num = new_tets.size() - delete_tets.size(); i < loop_num;
       i++)
    new_tets_id.push_back(i + T.size());
  if (delete_tets.size() > new_tets.size()) new_tets_id.resize(new_tets.size());

  // retain_connectivity.
  std::set<int> surround_id;
  std::set_difference(influence_id.begin(), influence_id.end(),
                      one_ring_tets.begin(), one_ring_tets.end(),
                      std::inserter(surround_id, surround_id.end()));
  retain_tetrahedral_adjacency(one_ring_tets, surround_id, new_tets, T, TT,
                               TTif, TTie);

  return true;
}
}  // namespace igl::dev