#ifndef BCCLEAN_TETSPLIT_EDGE_HASH_HPP
#define BCCLEAN_TETSPLIT_EDGE_HASH_HPP
#include <Eigen/Core>
#include <unordered_map>

#include "common.hpp"

// https://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine
template <typename T>
struct arr_hash : std::unary_function<T, size_t> {
  std::size_t operator()(std::array<T, 6> const& arr) const {
    size_t seed = 0;
    for (auto a : arr)
      seed ^= std::hash<T>()(a) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  }
};

using edge2tuple = std::unordered_map<std::array<double, 6>, std::array<int, 4>,
                                      arr_hash<double>>;

void make_hash(const std::vector<Vec3d>& V, const std::vector<Vec4i>& T,
               const std::vector<Vec4i>& TT, const std::vector<int>& newtets, double eps, edge2tuple& hashtable);

#endif