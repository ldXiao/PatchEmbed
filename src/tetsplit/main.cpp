#include <igl/boundary_facets.h>
#include <igl/readDMAT.h>
#include <igl/readMESH.h>
#include <igl/writeMESH.h>
#include <spdlog/spdlog.h>
#include <numeric>
#include "common.hpp"
#include "edge_hash.hpp"
#include "igl_dev/tet_split.h"
#include "igl_dev/tetrahedron_tetrahedron_adjacency.h"

#define TEST_FILE "../data/2/CC0/tet.mesh"
#define SPLIT_FILE "../data/2/CC0/splits_record.dmat"

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
  FT.resize(mFT.size());
  for (int i = 0; i < mFT.size(); i++) FT[i] = mFT[i];
  FTi.resize(mFTi.size());
  for (int i = 0; i < mFTi.size(); i++) FTi[i] = mFTi[i];
}

void writer(std::string filename, const std::vector<Vec3d>& V,
            const std::vector<Vec4i>& T) {
  RowMatd mV;
  RowMati mT, mF;
  vec2eigen(V, mV);
  vec2eigen(T, mT);
  igl::writeMESH(filename, mV, mT, mF);
}

void record_reader(std::string filename, Eigen::MatrixXd& record_mat) {
  igl::readDMAT(SPLIT_FILE, record_mat);
  spdlog::trace("Read DMAT, {} {}", record_mat.rows(), record_mat.cols());
  double sum = record_mat.col(0).sum();
  spdlog::trace("sum first row = {}", sum);
}

int main() {
  spdlog::set_level(spdlog::level::trace);

  std::vector<Vec3d> V;
  std::vector<Vec3i> F;
  std::vector<Vec4i> T, TT, TTif;
  std::vector<Eigen::Matrix<int, 4, 3>> TTie;
  std::vector<int> FT, FTi;

  Eigen::MatrixXd record_mat;
  record_reader(SPLIT_FILE, record_mat);
  reader(TEST_FILE, V, F, T, FT, FTi);
  igl::dev::tetrahedron_tetrahedron_adjacency(T, TT, TTif, TTie);

  // create hashing structure.
  edge2tuple hashtable;
  double eps = 1e-10;
  std::vector<int> alltets(T.size());
  std::iota(alltets.begin(), alltets.end(), 0);
  make_hash(V, T, TT, alltets, eps, hashtable);

  for (int i = 0; i < record_mat.rows(); i++) {
    if (record_mat(i, 0) == 0) {  // edge split
      std::array<double, 6> edge_pair;
      for (int k = 0; k < 6; k++)
        edge_pair[k] = std::round(record_mat(i, k + 1) / eps);
      Vec3d newvert(record_mat(i, 7), record_mat(i, 8), record_mat(i, 9));
      // auto key = edge_pair;
      // spdlog::info("query {} {} {} {} {} {}",
      // key[0],key[1],key[2],key[3],key[4], key[5]);
      auto query = hashtable.find(edge_pair);
      if (query == hashtable.end()) {
        spdlog::error("Not found in table: i {}", i);
        exit(1);
      }
      auto [ti, fi, ei, _] = query->second;
      std::vector<int> newtets;
      igl::dev::tet_tuple_edge_split(ti, fi, ei, true, V, T, TT, TTif, TTie,
                                     newtets);
      V.back() = newvert;
      make_hash(V,T, TT, newtets, eps, hashtable);
    }
    else if(record_mat(i,0)==1){// vertices insert
      spdlog::error("unsupported yet");
      return -1;
    }
  }
  writer("temp.mesh", V, T);
  spdlog::info("Vert{}  Tets{}", V.size(), T.size());
  return 0;
}