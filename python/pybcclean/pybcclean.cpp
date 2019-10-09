#include <pybind11/pybind11.h>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <vector>
#include <utility>
#include "bcclean.h"
#include "edge.h"
#include "graphcut_cgal.h"
#include <igl/boundary_loop.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>

void project_face_labels(
    const Eigen::MatrixXd &V_bad, 
    const Eigen::MatrixXi &F_bad, 
    const Eigen::MatrixXi &FL_bad,
    const Eigen::MatrixXd &V_good,
    const Eigen::MatrixXi &F_good,
    const int label_num,
    const int subdiv,
    Eigen::MatrixXi & FL_good,
    Eigen::MatrixXd & prob_mat){
        // int label_num = 0;
        // std::map<int, int> count_dict;
        // for(int i =0; i< FL_bad.rows(); ++i){
        //     auto it = count_dict.find(FL_bad(i,0));
        //     if(it == count_dict.end()){
        //         label_num+=1;
        //         count_dict[FL_bad(i,0)]=1;
        //     }
        // }
        Eigen::MatrixXi VL_good, FL_good_temp;
        std::cout << "r1"<< std::endl;
        // bcclean::LM_intersection_label_transport(V_bad,F_bad,FL_bad,V_good,F_good,VL_good);
        // std::cout << "r2"<< std::endl;
        // // bcclean::Barycenter_intersection_label_transport(V_bad,F_bad,FL_bad,V_good,F_good,FL_good);
        // bcclean::vertex_label_vote_face_label(label_num, VL_good, F_good, FL_good, prob_mat);
        bcclean::refine_proj_vote(V_bad,F_bad,FL_bad, V_good, F_good, label_num, subdiv,FL_good, prob_mat);
        std::cout << "r2"<< std::endl;
}

void refine_labels(
    const Eigen::MatrixXd &V_good,
    const Eigen::MatrixXi &F_good,
    const Eigen::MatrixXd &prob_mat,
    Eigen::MatrixXi & FL_good,
    double lambda_refine){
    bcclean::refine_labels_graph_cut(V_good,F_good, prob_mat.transpose(), FL_good, lambda_refine);
}

class Node{
    public:
        int vidx;
        std::vector<int> labels;
        void initialize(const int vidx_, const std::vector<int> & labels_){
            vidx = vidx_;
            labels = labels_;
            std::sort(labels.begin(), labels.end());
        }
};

class Edge{
    public:
        int head;
        int tail;
        std::vector<int> vertices_list;
        std::pair<int, int> label_pair;
        void initialize(const int head_, const int tail_, const std::pair<int,int> & label_pair_, std::vector<int> vertices_list_){
            head = head_;
            tail = tail_;
            label_pair = label_pair_;
            vertices_list = vertices_list_;
        }
};

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
    m.def(
        "project_face_labels",
        [](
        py::EigenDRef<Eigen::MatrixXd> V_bad, 
        py::EigenDRef<Eigen::MatrixXi> F_bad,
        const Eigen::MatrixXi FL_bad,
        int label_num,
        int subdiv,
        py::EigenDRef<Eigen::MatrixXd> V_good, 
        py::EigenDRef<Eigen::MatrixXi> F_good
        ) {
            Eigen::MatrixXd prob_mat;
            Eigen::MatrixXi FL_good;
            project_face_labels(V_bad, F_bad, FL_bad, V_good, F_good, label_num, subdiv, FL_good, prob_mat);
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
        bcclean::refine_labels_graph_cut(V_good, F_good, prob_mat.transpose(), FL_good_cut,lambda_refine);
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

    m.def("reorient",
        [](
           py::EigenDRef<Eigen::MatrixXd> V, 
           py::EigenDRef<Eigen::MatrixXi> F 
        ){
            Eigen::MatrixXi FF;
            Eigen::VectorXi I;
            bcclean::reorient(V, F, FF, I);
            return py::make_tuple(FF, I);
        }
    );

    py::class_<Node>(m, "Node")
        .def(py::init<>())
        .def_readwrite("vidx", &Node::vidx)
        .def_readwrite("labels", &Node::labels);

    py::class_<Edge>(m, "Edge")
        .def(py::init<>())
        .def_readwrite("head", &Edge::head)
        .def_readwrite("tail", &Edge::tail)
        .def_readwrite("label_pair", &Edge::label_pair)
        .def_readwrite("vertices_list", &Edge::vertices_list);
    
    m.def("build_edge_list",
    [](const Eigen::MatrixXd V, const Eigen::MatrixXi F, const Eigen::MatrixXi & FL, const int total_label_num){
        std::vector<bcclean::edge> edge_list;
        std::vector<Edge> Edge_list;
        std::unordered_map<int, std::vector<int> > patch_edge_dict;
        bcclean::build_edge_list(V, F, FL, total_label_num, edge_list, patch_edge_dict);
        
        for(auto edg: edge_list){
            Edge EDG;
            EDG.initialize(edg.head, edg.tail, edg._label_pair, edg._edge_vertices);
            Edge_list.push_back(EDG);
        }
        return py::make_tuple(Edge_list, patch_edge_dict)
        ;},
        R"pbdoc(
        detect return a list of nodes
        V: numpy array of #V x 3,
        F: numpy array of #F x 3,
        FL: numpy array of #F x1,
        label_num: total number of labels
        return  tuple(edge_list: list(Edge), patch_edge_dict: dict{label:list(edge_indx)})
        )pbdoc"
    );

    m.def("build_node_list",
    [](const Eigen::MatrixXi F, const Eigen::MatrixXi & FL, const int total_label_num){
        std::unordered_map<int, std::vector<int> > vertex_label_list_dict;
        std::vector<Node> node_list;
        bcclean::build_vertex_label_list_dict(F, FL, total_label_num, vertex_label_list_dict);
        for(auto item: vertex_label_list_dict){
            if(item.second.size()>2){
                Node nd;
                nd.initialize(item.first,item.second);
                node_list.push_back(nd);
            }
        }
        return node_list;
    },
    R"pbdoc(
        detect return a list of nodes
        F: numpy array of #F x 3,
        FL: numpy array of #F x1,
        return  list(Node), list of Node
        )pbdoc"
    );

    m.def("boundary_loops",
    [](const Eigen::MatrixXi & F){
        std::vector<std::vector<int> > L;
        igl::boundary_loop(F, L);
        return L;
    });
    

    

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
