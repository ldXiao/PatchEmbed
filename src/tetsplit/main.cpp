#include <igl/boundary_facets.h>
#include <igl/readMESH.h>
#include <igl/writeMESH.h>
#include "igl_dev/tet_split.h"
#include "igl_dev/tetrahedron_tetrahedron_adjacency.h"
#include <spdlog/spdlog.h>
#define TEST_FILE "../data/0/good.mesh"

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

void reader(std::string filename, std::vector<Vec3d>& V, std::vector<Vec3i>& F,
            std::vector<Vec4i>& T, std::vector<int>& FT,
            std::vector<int>& FTi) {
  // we also need the map from F to T, and the corner info.
  RowMatd mV;
  RowMati mT;
  RowMati mF;
  igl::readMESH(filename, mV, mT, mF);
  Eigen::VectorXi mFT, mFTi;
  if (mF.rows() == 0) igl::boundary_facets(mT, mF, mFT, mFTi);
  spdlog::info("Vert{} Face{} Tets{}", mV.rows(), mF.rows(), mT.rows());

  eigen2vec(mV, V);
  eigen2vec(mF, F);
  eigen2vec(mT, T);
  FT.resize(mFT.size()); for (int i=0; i<mFT.size(); i++) FT[i] = mFT[i];
  FTi.resize(mFTi.size()); for (int i=0; i<mFTi.size(); i++) FTi[i] = mFTi[i];
}

void writer(std::string filename, const std::vector<Vec3d>& V,
            const std::vector<Vec4i>& T) {
  RowMatd mV;
  RowMati mT, mF;
  vec2eigen(V, mV);
  vec2eigen(T, mT);
  igl::writeMESH(filename, mV, mT, mF);
}

int main() {
  spdlog::set_level(spdlog::level::trace);

  std::vector<Vec3d> V;
  std::vector<Vec3i> F;
  std::vector<Vec4i> T, TT, TTif;
  std::vector<Eigen::Matrix<int, 4, 3>> TTie;
  std::vector<int> FT, FTi;
  reader(TEST_FILE, V, F, T, FT, FTi);
  igl::dev::tetrahedron_tetrahedron_adjacency(T, TT, TTif, TTie);

  std::vector<int> newtets;
  int fid = 10, eid = 0;

  int ti = FT[fid], fi = FTi[fid];
  igl::dev::tet_tuple_edge_split(ti, fi, eid, true, V, T, TT, TTif, TTie,
                                 newtets);
  writer("temp.mesh", V, T);
  spdlog::info("Vert{}  Tets{}", V.size(), T.size());
  return 0;
}