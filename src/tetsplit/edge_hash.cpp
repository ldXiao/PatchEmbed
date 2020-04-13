#include "edge_hash.hpp"

#include "igl_dev/tetrahedron_tetrahedron_adjacency.h"
#include <spdlog/spdlog.h>
void make_hash(const std::vector<Vec3d>& V, const std::vector<Vec4i>& T,
               const std::vector<Vec4i>& TT, const std::vector<int>& newtets, double eps,
               edge2tuple& hashtable) {
  auto& FF = igl::dev::tetrahedron_local_FF;
  for (auto t:newtets) {
    for (int f = 0; f < 4; f++) {
      if (TT[t][f] != -1) continue;  // only hash boundary edges.
      for (int e = 0; e < 3; e++) {
        auto v0 = T[t][FF(f, e)], v1 = T[t][FF(f, (e + 1) % 3)];
        std::array<double, 6> key{V[v0][0], V[v0][1], V[v0][2],
                                  V[v1][0], V[v1][1], V[v1][2]};
        for (int k=0;k<6;k++) {
          key[k] = std::round(key[k] / eps);
        };
        // spdlog::trace("key in {} {} {} {} {} {}", key[0],key[1],key[2],key[3],key[4], key[5]);
        hashtable[key] = {t, f, e, 1};
      }
    }
  }
}