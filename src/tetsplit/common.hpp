#ifndef BCCLEAN_TETSPLIT_COMMON_HPP
#define BCCLEAN_TETSPLIT_COMMON_HPP

#include <Eigen/Core>
using RowMatd = Eigen::Matrix<double, -1, -1, Eigen::RowMajor>;
using RowMati = Eigen::Matrix<int, -1, -1, Eigen::RowMajor>;
using Vec3d = Eigen::RowVector3d;
using Vec3i = Eigen::RowVector3i;
using Vec4i = Eigen::RowVector4i;

inline constexpr auto vec2eigen = [](const auto& vec, auto& mat) {
  mat.resize(vec.size(), vec[0].size());
  for (int i = 0; i < mat.rows(); i++) {
    for (int j = 0; j < mat.cols(); j++) mat(i, j) = vec[i][j];
  }
};

constexpr auto eigen2vec = [](const auto& mat, auto& vec) {
  vec.resize(mat.rows());
  for (int i = 0; i < vec.size(); i++)
    for (int j = 0; j < mat.cols(); j++) vec[i][j] = mat(i, j);
};

#endif