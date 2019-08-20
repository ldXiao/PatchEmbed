#include <pybind11/pybind11.h>
#include <map>
#include "bcclean.h"
#include "graphcut_cgal.h"
#include <pybind11/eigen.h>

int add(int i, int j) {
    return i + j;
}

void project_face_labels(
    const Eigen::MatrixXd &V_bad, 
    const Eigen::MatrixXi &F_bad, 
    const Eigen::MatrixXi &FL_bad,
    const Eigen::MatrixXd &V_good,
    const Eigen::MatrixXi &F_good,
    Eigen::MatrixXi & FL_good,
    Eigen::MatrixXd & prob_mat){
        int label_num = 0;
        std::map<int, int> count_dict;
        for(int i =0; i< FL_bad.rows(); ++i){
            auto it = count_dict.find(FL_bad(i,0));
            if(it == count_dict.end()){
                label_num+=1;
                count_dict[FL_bad(i,0)]=1;
            }
        }
        Eigen::MatrixXi VL_good;
        bcclean::LM_intersection_label_transport(V_bad,F_bad,FL_bad,V_good,F_good,VL_good);
        bcclean::vertex_label_vote_face_label(label_num, VL_good, F_good, FL_good, prob_mat);
}

void refine_labels(
    const Eigen::MatrixXd &V_good,
    const Eigen::MatrixXi &F_good,
    const Eigen::MatrixXd &prob_mat,
    Eigen::MatrixXi & FL_good,
    double lambda_refine){
    bcclean::refine_labels_graph_cut(V_good,F_good, prob_mat.transpose(), FL_good, lambda_refine);
}

namespace py = pybind11;

PYBIND11_MODULE(pybcclean, m) {
    m.doc() = R"pbdoc(
        Pybind11 gen pybcclean 
        -----------------------

        .. currentmodule:: pybcclean

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";


    m.def("add_any", [](py::EigenDRef<Eigen::MatrixXd> x, int r, int c, double v) { x(r,c) += v; });
    m.def(
        "project_face_labels",
        [](
        py::EigenDRef<Eigen::MatrixXd> V_bad, 
        py::EigenDRef<Eigen::MatrixXi> F_bad,
        py::EigenDRef<Eigen::MatrixXi> FL_bad,
        py::EigenDRef<Eigen::MatrixXd> V_good, 
        py::EigenDRef<Eigen::MatrixXi> F_good
        ) {
            Eigen::MatrixXd prob_mat;
            Eigen::MatrixXi FL_good;
            project_face_labels(V_bad, F_bad, FL_bad, V_good, F_good, FL_good, prob_mat);
            return py::make_tuple(prob_mat, FL_good);
        },
        R"pbdoc(
        project face labels onto good_mesh

        V_bad: numpy array of #V_bad x 3
        F_bad: numpy array of #F_bad x 3,
        FL_bad: numpy array #F x 1, contain #label 
        V_bad: numpy array of #V_good x 3,
        F_bad: numpy array of #F_good x 3,
        return FL_good, # F_good x 1, contain #label
                prob_mat # F_good x #labels
        )pbdoc"
    );
    m.def("refine_labels",
    [](
        py::EigenDRef<Eigen::MatrixXd> V_good, 
        py::EigenDRef<Eigen::MatrixXi> F_good,
        py::EigenDRef<Eigen::MatrixXd> prob_mat,
        Eigen::MatrixXi& FL_good,
        double lambda_refine
    ){
        Eigen::MatrixXi FL_good_cut= FL_good;
        bcclean::refine_labels_graph_cut(V_good,F_good, prob_mat.transpose(), FL_good_cut,lambda_refine);
        return FL_good_cut;
    },
    R"pbdoc(
        project face labels onto good_mesh

        V_good: numpy array of #V_good x 3,
        F_good: numpy array of #F_good x 3,
        FL_good: numpy array of #F_good x1,
        prob_mat: numpy array of #F_good x #label,
        lambda_refine: a positive constant to control graph_cut
        return FL_good_cut, # F_good x 1, contain #label
        )pbdoc"
    );
    

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
