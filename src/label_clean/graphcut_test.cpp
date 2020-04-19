//
// Created by Lind Xiao on 6/11/19.
//

#include <cxxopts.hpp>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/writeOBJ.h>
#include <igl/read_triangle_mesh.h>
#include <igl/jet.h>
#include <vector>

#include <memory>

#include <utility>
#include "../graphcut_cgal.h"
#include "nn.hpp"

Eigen::MatrixXd Vs, Vt, prob_mat;
Eigen::MatrixXi Ft, Fs;
Eigen::VectorXi FLs, FLt, FLv, VL;

void vertex_label_vote_face_label(
            const int label_num,
            const Eigen::VectorXi &VL,
            const Eigen::MatrixXi &F,
            Eigen::VectorXi &FL,
            Eigen::MatrixXd &prob_mat
    ) {
        FL = Eigen::VectorXi::Constant(F.rows(),0);
        prob_mat = Eigen::MatrixXd::Constant(FL.rows(), label_num, 0.1 / label_num);
        std::map<int, int> dict;
        for (int fidx = 0; fidx < F.rows(); ++fidx) {
            // int v0, v1, v2, l0, l1, l2;
            dict.clear();
            for (int i = 0; i < 3; ++i) {
                int vi = F(fidx, i);
                int li = VL(vi);
                auto find = dict.find(li);
                if (find != dict.end()) {
                    dict[li] += 1;
                } else {
                    dict[li] = 1;
                }
            }
            std::cout << "built dict" <<std::endl;
            if (dict.size() == 3) {
                FL(fidx, 0) = dict.begin()->first;
                auto it = dict.begin();
                int l0 = it->first;
                it++;
                int l1 = it->first;
                it++;
                int l2 = it->first;

                prob_mat(fidx, l0) = 0.3;
                prob_mat(fidx, l1) = 0.3;
                prob_mat(fidx, l2) = 0.3;

            }
            std::cout << "b1" <<std::endl;
            if (dict.size() == 2) {
                auto it = dict.begin();
                int l0 = it->first;
                it++;
                int l1 = it->first;
                if (dict[l0] > dict[l1]) {
                    FL(fidx, 0) = l0;
                    prob_mat(fidx, l0) = 0.6;
                    prob_mat(fidx, l1) = 0.3;
                } else {
                    FL(fidx, 0) = l1;
                    prob_mat(fidx, l1) = 0.6;
                    prob_mat(fidx, l0) = 0.3;
                }
            }
            std::cout << "b2" <<std::endl;
            if (dict.size() == 1) {
                FL(fidx, 0) = dict.begin()->first;
                int l0 = dict.begin()->first;
                prob_mat(fidx, l0) = 0.9;
            }
        }


    }


int main(int argc, char *argv[]){
using namespace std;
    using namespace Eigen;
    if (argc < 2) {
        cout << "Usage cutmeshtest_bin --file mesh.obj --n sample_num --cut" << endl;
        exit(0);
    }
    cxxopts::Options options("CutMesh", "One line description of MyProgram");
    options.add_options()
            ("m, mesh", "",cxxopts::value<std::string>())
            ("l, label", "",cxxopts::value<std::string>())
            ("t, target","", cxxopts::value<std::string>())
            ;
    auto args = options.parse(argc, argv);
    // Load a mesh in OBJ format
    

    
    igl::read_triangle_mesh(args["mesh"].as<std::string>(), Vs, Fs);
    igl::readDMAT(args["label"].as<std::string>(), FLs);
    igl::read_triangle_mesh(args["target"].as<std::string>(), Vt, Ft);
    nn_transfer(Vs, Fs, FLs, Vt,Ft, FLt, VL);
    igl::writeOBJ("../../../blenders/nn.obj", Vt, Ft);
    igl::writeDMAT("../../../blenders/nnfl.dmat", FLt);
    int label_num = FLs.maxCoeff()+1;
    vertex_label_vote_face_label(label_num,VL,Ft,FLt,prob_mat);
    igl::writeDMAT("../../../blenders/vote.dmat", FLt);
    Eigen::MatrixXi FLt_copy = FLt;
    bcclean::refine_labels_graph_cut(Vt,Ft, prob_mat.transpose(), FLt_copy,1);
    igl::writeDMAT("../../../blenders/gcgl.dmat", FLt_copy);
    return 0;

}

