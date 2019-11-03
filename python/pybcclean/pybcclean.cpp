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
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>

 void splits_detect(
        const Eigen::MatrixXi & F,
        const Eigen::MatrixXi & TT, // triangle-triangle adjacency
        const std::vector<bool> & VCuts, // indicate whether a vertex is on the boundary
        const std::vector<std::vector<bool> > & TCuts, // indicate whether wich edge of a face is on the boundary
        std::map<std::pair<int, int>, int > & splits
    )
    {
        // splits stores pairs of adjacent vertices indices that need to be splited
        // indiced into F
        splits.clear();
        // loop over all faces
        for(int fidx=0; fidx < F.rows(); ++fidx){
            for(int edge_idx=0; edge_idx < 3; ++ edge_idx){
                int v0 = F(fidx, edge_idx);
                int v1 = F(fidx, (edge_idx+1) % 3);
                if(VCuts[v0] && VCuts[v1] && !TCuts[fidx][edge_idx]){
                    int cofidx = TT(fidx, edge_idx);
                    assert(cofidx != -1);
                    if(cofidx != -1){
                        int vl = std::min(v0, v1);
                        int vg = std::max(v0, v1);
                        splits[std::make_pair(vl, vg)] = 1;
                    }
                    // make sure each pair of indices are only pushed onece
                }
            }
        }
        return;
        
    }

bool determin_adj_configure1(
        const Eigen::MatrixXi & F, 
        const std::vector<std::vector<int > > & VF, 
        const int uidx, 
        const int vidx,
        int & fupidx,
        int & fdownidx,
        int & vupidx,
        int & vdownidx)
    {
        std::vector<int> inter(VF[uidx].size()+ VF[vidx].size());
        auto it = std::set_intersection(VF[uidx].begin(), VF[uidx].end(), VF[vidx].begin(), VF[vidx].end(), inter.begin());
        inter.resize(it-inter.begin());
        if(inter.size()!=2) return false; // there should be always only two triangles
        
        // choose  the triangle that has positive orientation on u -> v
        // should update VF in the end of loop


        /* decide up and down faces vertices
            vup
            /   \
           / fup \
          u ----- v
           \fdown/
            \   /
            vdown

            both triangle oriented outward the screen
            */
        

        for(auto trig: inter)
        {
            for(int edgepos=0; edgepos<3; ++ edgepos)
            {
                if(F(trig, edgepos) == uidx && F(trig,(edgepos+1)%3) == vidx){
                    fupidx = trig;
                    for(auto other:inter){
                        if(other!= fupidx) fdownidx = other;
                    }
                    // decide the up and down triangles


                    vupidx = F(trig, (edgepos+2)%3);
                    for(int downedge=0; downedge < 3;++downedge)
                    {
                        if(F(fdownidx, downedge)!= uidx && F(fdownidx, downedge)!= vidx) vdownidx = F(fdownidx, downedge);
                    }
                }
            }
        }
        return true;
    }

 void splits_update(
        const std::map<std::pair<int, int>, int> & splits,
        Eigen::MatrixXd & Vraw, // raw mesh
        Eigen::MatrixXi & Fraw,
        std::vector<bool> & VCuts, // indicate whether a vertex is on the boundary
        std::vector<std::vector<bool> > & TCuts,
        std::vector<std::vector<int> > & VV // adjacency list on Fraw, posibly some ofthem has been silenced
    )
    {
        /*
        splits stores pairs of adjacent triangle indices and corresponding vertives that need to be splited

        this funciton will correspondingly  raw mesh also the face mapping
        
        everytime split two faces and create four smaller triangles add two of them to the end of Fraw Fbase

        also add the face mapping and cut info into the end of FI, VCuts, TCuts
        */

        const int task_num = splits.size();
        // each split task will increace face num by 2 and vertices num by 1
        const int nF_num_raw = task_num * 2 + Fraw.rows();
        const int nV_num_raw = task_num + Vraw.rows();

        // preallocate all the memory
        Eigen::MatrixXd nVraw = Eigen::MatrixXd::Zero(nV_num_raw , 3); // raw mesh
        nVraw.block(0, 0, Vraw.rows(), 3) = Vraw;
        Eigen::MatrixXi nFraw = Eigen::MatrixXi::Zero(nF_num_raw, 3);
        nFraw.block(0, 0, Fraw.rows(), 3) = Fraw;


        std::vector<bool> nVCuts(nV_num_raw); // indicate whether a vertex is on the boundary
        std::vector<std::vector<bool> > nTCuts(nF_num_raw);
        VV.resize(nV_num_raw);
        for(int count =0; count <nV_num_raw; ++ count)
        {
            if(count < VCuts.size()) nVCuts[count] = VCuts[count];
            else nVCuts[count] = false;
        }
        for(int fcount = 0; fcount  < nF_num_raw; ++fcount)
        {   
            if(fcount < TCuts.size()) nTCuts[fcount] = TCuts[fcount];
            else nTCuts[fcount] = {false, false, false};
        }
        std::vector<std::vector<int > > nVF_raw;// #V list of lists of incident faces (adjacency list)
        // the naming indicates that they are to be updated by hand in each loop
        {
            std::vector<std::vector<int > > II;  //  #V list of lists of index of incidence within incident faces listed in VF, local variable do not care about them
            igl::vertex_triangle_adjacency(Vraw, Fraw, nVF_raw, II);
        }
        // initialize VFs

        // because splits might be connected directly, this process has to be done one by one
        // after  spliting each edge, the info in splits needs not to be updated because
        // splitting does not change initial indices of existing vertices
        int task_count = 0;
        for(auto item: splits){
            int uidx_raw = item.first.first;
            int vidx_raw = item.first.second;
            int widx_raw = Vraw.rows() + task_count;
            Eigen::RowVector3d  wpos = (Vraw.row(uidx_raw)+Vraw.row(vidx_raw))/2;
            nVraw.row(widx_raw) = wpos;
            // they are two representations of an identical vertex
            







            // decide the adj configurations
            /* decide up and down faces vertices
                vup
                /   \
                / fup \
            u ----- v
                \fdown/
                \   /
                vdown

                both triangle oriented outward the screen
                */

            int fupidx_raw, fdownidx_raw;
            int vupidx_raw, vdownidx_raw;
            determin_adj_configure1(
                nFraw, nVF_raw, uidx_raw, vidx_raw,
                fupidx_raw, fdownidx_raw, vupidx_raw, vdownidx_raw
            );
            // decide the indices of existing triangles and vertces

            


            // create new triangles

            /* 
                    vup
                    / | \
                    /  |  \
                    / f1|f0 \
                /    |    \
                u-----w ----v
                \    |    /
                    \ f2| f3/
                    \  |  /
                    \ | /
                    vdown
            */
            int f0idx_raw, f1idx_raw, f2idx_raw, f3idx_raw;
            f0idx_raw = fupidx_raw;
            f2idx_raw = fdownidx_raw;
            f1idx_raw = Fraw.rows() + 2 * task_count;
            f3idx_raw = f1idx_raw + 1;
            // update TCuts before updating nFraw
            //  poseone updates of nFraw after TCuts because Tcuts relies on Fraw
            // update nVCuts
            // only one vertices is added
            nVCuts[widx_raw] = false;



            // update nTCuts
            std::map<std::pair<int,int> , bool> nTCuts_record;
            // store cut info of all 5 edges of fupidx_raw, and fdownidx_raw in nTCuts resulted from last loop
            for(int edgepos =0 ; edgepos < 3 ; ++edgepos){
                for(auto updown: {fupidx_raw, fdownidx_raw}){
                    int uu = nFraw(updown, edgepos);
                    int vv = nFraw(updown, (edgepos+1)%3);
                    int vl = std::min(uu,vv);
                    int vg = std::max(uu, vv);
                    nTCuts_record[std::make_pair(vl, vg)] = nTCuts[updown][edgepos];
                }
            }
            std::pair<int, int> key; 
            key = std::make_pair(std::min(vidx_raw, vupidx_raw), std::max(vidx_raw, vupidx_raw));
            nTCuts[f0idx_raw] = {false, nTCuts_record[key], false};

            key = std::make_pair(std::min(uidx_raw, vupidx_raw), std::max(uidx_raw, vupidx_raw));
            nTCuts[f1idx_raw] = {false, nTCuts_record[key], false};

            key = std::make_pair(std::min(uidx_raw, vdownidx_raw), std::max(uidx_raw, vdownidx_raw));
            nTCuts[f2idx_raw] = {false, nTCuts_record[key], false};

            key = std::make_pair(std::min(vidx_raw, vdownidx_raw), std::max(vidx_raw, vdownidx_raw));
            nTCuts[f3idx_raw] = {false, nTCuts_record[key], false};
            


            // update nFraw the faces with correct orientation
            nFraw.row(f0idx_raw) = Eigen::RowVector3i(widx_raw, vidx_raw, vupidx_raw);
            nFraw.row(f1idx_raw) = Eigen::RowVector3i(widx_raw, vupidx_raw, uidx_raw);
            nFraw.row(f2idx_raw) = Eigen::RowVector3i(widx_raw, uidx_raw, vdownidx_raw);
            nFraw.row(f3idx_raw) = Eigen::RowVector3i(widx_raw, vdownidx_raw, vidx_raw);
            
        


            
            // update nVFs VV the connectivitity info
            // update nVF_raw;
            /*  only have to deal with f1 f3 becase they are new faces
                    vup
                    / | \
                    /  |  \
                    / f1|f0 \
                /    |    \
                u-----w ----v
                \    |    /
                    \ f2| f3/
                    \  |  /
                    \ | /
                    vdown
            */
            nVF_raw.resize(nV_num_raw);
            
            nVF_raw[vidx_raw].push_back(f3idx_raw);
            nVF_raw[vidx_raw].erase(std::remove(nVF_raw[vidx_raw].begin(), nVF_raw[vidx_raw].end(), f2idx_raw), nVF_raw[vidx_raw].end()); 

            nVF_raw[vupidx_raw].push_back(f1idx_raw);

            nVF_raw[uidx_raw].push_back(f1idx_raw);
            nVF_raw[uidx_raw].erase(std::remove(nVF_raw[uidx_raw].begin(), nVF_raw[uidx_raw].end(), f0idx_raw), nVF_raw[uidx_raw].end());

            nVF_raw[vdownidx_raw].push_back(f3idx_raw);
            nVF_raw[widx_raw] = {vidx_raw, vupidx_raw, uidx_raw, vdownidx_raw};



            for(auto idx_raw: {vidx_raw, vupidx_raw, uidx_raw, vdownidx_raw}){
                if(VV[idx_raw].size()!=0){
                    // idx_raw has not been silenced add connection
                    VV[idx_raw].push_back(widx_raw);
                    VV[widx_raw].push_back(idx_raw);
                }
            }

            // also silence the connection between vidx_raw and uidx_raw
            // if one of them are silenced no need to further remove the other
            VV[vidx_raw].erase(std::remove(VV[vidx_raw].begin(), VV[vidx_raw].end(), uidx_raw), VV[vidx_raw].end()); 
            VV[uidx_raw].erase(std::remove(VV[uidx_raw].begin(), VV[uidx_raw].end(), vidx_raw), VV[uidx_raw].end());

            // the above two lines will do nothing if they find nothing

            task_count += 1; // finishes one task
        }


        Vraw = nVraw;
        Fraw = nFraw;
        VCuts = nVCuts;
        TCuts = nTCuts;
    }

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
        detect return a list of edges 
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

    m.def("splits_detect",
    []
        (
            const Eigen::MatrixXi & F,
            const Eigen::MatrixXi & TT,
            const std::vector<bool> & VCuts, // indicate whether a vertex is on the boundary
            const std::vector<std::vector<bool> > & TCuts){
                std::map<std::pair<int, int>, int > splits;
                splits_detect(F, TT, VCuts, TCuts, splits);
                std::vector<std::pair<int, int> > vsplits;
                for(auto item: splits){
                    vsplits.push_back(item.first);
                }
                return vsplits;
            }
    );

    m.def("splits_update",
    [](const std::vector<std::pair<int, int> > & vsplits,
        Eigen::MatrixXd  Vraw, // raw mesh
        Eigen::MatrixXi  Fraw,
        std::vector<bool>  VCuts, // indicate whether a vertex is on the boundary
        std::vector<std::vector<bool> >  TCuts,
        std::vector<std::vector<int> >  VV){
            std::map<std::pair<int, int>, int > splits;
                for(auto item: vsplits){
                    splits[item]=1;
                } 
            splits_update(splits, Vraw, Fraw, VCuts, TCuts, VV);
            return py::make_tuple(Vraw, Fraw, VCuts, TCuts, VV);

        }
    );
    

    

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
