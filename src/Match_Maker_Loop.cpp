#include "Match_Maker_Loop.h"
#include "Match_Maker.h"
namespace bcclean{
namespace MatchMaker{
    using json = nlohmann::json;
    void trace_and_label_loop(
        const Eigen::MatrixXd & V_bad,
        const Eigen::MatrixXi & F_bad,
        const Eigen::VectorXi & FL_bad,
        Eigen::MatrixXd & V_good,
        Eigen::MatrixXi & F_good,
        Eigen::VectorXi & FL_good
    )
    {
        // datas dump to file for debug
        std::map<int, int> edge_order_map; // store and maintain the order of added edges {order: edge_dx}
        std::map<int, std::vector<int> > edge_path_map; // {edge_idx, path}
        int run_count = 0;
        // Randomize Seed
        srand(static_cast<unsigned int>(time(nullptr)));
        int total_label_num = FL_bad.maxCoeff()+1;


        // PART 0 GET THE FRAME ON BAD MESH
        std::vector<bcclean::edge> edge_list;
        std::unordered_map<int, std::vector<int> > patch_edge_dict;
        std::unordered_map<int, std::vector<bool> > patch_edge_direction_dict;
        // label -> list(edge_idx) map indices into edge_list
        build_edge_list_loop(V_bad, F_bad, FL_bad, total_label_num, edge_list, patch_edge_dict, patch_edge_direction_dict);

        // we can assume that  all nodes have valance more than 3
        std::vector<int> node_list_bad;
        {
            std::unordered_map<int, std::vector<int> > vertex_label_list_dict;
            build_vertex_label_list_dict(F_bad, FL_bad, total_label_num, vertex_label_list_dict);
            for(auto item: vertex_label_list_dict)
            {
                if(item.second.size()>2)
                {
                    node_list_bad.push_back(item.first);
                }
            }
        }

        
        std::map<int, std::vector<int> > node_edgepool_dict;
        {
            int edg_idx =0;
            for(auto edg:edge_list)
            {
                int head= edg.head;
                int tail=edg.tail;
                assert(head != tail); // both are nodes
                assert(edg._edge_vertices.size()>1);
                for(auto ht:{head, tail})
                {
                    auto it = std::find(node_edgepool_dict[ht].begin(), node_edgepool_dict[ht].end(), edg_idx);
                    if(it == node_edgepool_dict[ht].end())
                    {
                        // push it into list
                        node_edgepool_dict[ht].push_back(edg_idx);
                    }
                }
                edg_idx +=1;
            }
        }

        // the resulting node_edgepool_dict will be non-orderd but will be used later for counter clockwise ordered dictionary
        

        // PART 0.1 get a counter clock orderd map for node_idx -> edge_idx going out of node_idx
        std::map<int, std::vector<int> > node_faces_dict_bad;
        CC_faces_per_node(V_bad, F_bad, node_list_bad, node_faces_dict_bad);
        // the  faces in node_faces_dict_bad[nd] is ordered Counter clockwise;
        for(auto item: node_faces_dict_bad)
        {
            printf("\n");
            printf("node  %d\n", item.first);
            printf("face CC: ");
            for(auto x: item.second){
                printf("%d, ", x);
            }
        }
        

        // also order the node_edge_dict;
        // firsly find th CC ordered out-Vertex of each nd
        std::map<int, std::vector<int> > node_outv_dict_bad;
        for(auto item: node_faces_dict_bad)
        {
            int nd = item.first;
            node_outv_dict_bad[item.first] = {};
            for(auto fidx : item.second)
            {
                std::vector<int> pool = {F_bad(fidx,0), F_bad(fidx,1), F_bad(fidx,2)};
                int bench = -1;
                int ii = 0;
                for(auto v_f: pool)
                {
                    if(v_f == nd)
                    {
                        bench = ii;
                        break;
                    }
                    ii += 1;
                }
                int outv = pool[(bench+1)%3];
                node_outv_dict_bad[nd].push_back(outv);
            }
        }

        std::map<int, std::vector<int> > node_edge_dict;
        for(auto item: node_outv_dict_bad)
        {
            int nd = item.first;
            node_edge_dict[nd] = std::vector<int>();
            std::vector<int> edge_pool = node_edgepool_dict[nd];
            for(auto outv : node_outv_dict_bad[nd])
            {
                for(auto edg_idx: edge_pool)
                {
                    edge  edg = edge_list[edg_idx];
                    auto it = std::find(edg._edge_vertices.begin(), edg._edge_vertices.end(), outv);
                    auto jt = std::find(node_edge_dict[nd].begin(), node_edge_dict[nd].end(), edg_idx);
                    if(it != edg._edge_vertices.end() and jt == node_edge_dict[nd].end())
                    {
                        node_edge_dict[nd].push_back(edg_idx);
                    }

                }
            }
        }

        std::vector<int> sorted_node_list_bad = node_list_bad;
        std::map<int , std::vector<bool> > node_edge_visit_dict;
        std::map<int, bool> node_visit_dict;
        std::sort(sorted_node_list_bad.begin(), sorted_node_list_bad.end(), [node_edge_dict](int nda, int ndb){return (node_edge_dict.at(nda).size() > node_edge_dict.at(ndb).size());});
        for(auto nd: sorted_node_list_bad)
        {
            for(auto q: node_edge_dict[nd])
            {
                node_edge_visit_dict[nd].push_back(false);
            }
            node_visit_dict[nd] = false;
        }


        // we will use only adjacency list VV, bool lists VEdges, TEdges
        // initializations
        std::vector<std::vector<int> > VV_good, VF_good;
        igl::adjacency_list(F_good, VV_good);
        {
            std::vector<std::vector<int> >  VFi_good;
            igl::vertex_triangle_adjacency(V_good, F_good, VF_good, VFi_good);

        }
        std::vector<std::vector<int> > VEdges_good(V_good.rows());
        std::vector<std::vector<int> > TEdges_good(F_good.rows());
        for(int count =0; count <V_good.rows(); ++ count)
        {
            VEdges_good[count] = std::vector<int>();
        }
        for(int fcount = 0; fcount  < F_good.rows(); ++fcount)
        {   
            TEdges_good[fcount] = {-1, -1, -1};
        }
        // start with any of the patches and trace the loop of the disk
        for(auto item: patch_edge_dict)
        {
            int lb = item.first;
            std::vector<bool> direction_list = patch_edge_direction_dict[lb];
            std::vector<int> loop_good; // store the vertex indices of the loop in good_mesh
            int local_edg_count =0;
            for(auto edg_idx:item.second)
            {
                int ndhead_bad, ndtail_bad;
                if(patch_edge_direction_dict[lb][local_edg_count])
                {
                    ndhead_bad = patch_edge_dict[lb][0];
                    ndtail_bad = patch_edge_dict[lb].back();
                    // when the direction is correct
                } else{
                    ndtail_bad = patch_edge_dict[lb][0];
                    ndhead_bad = patch_edge_dict[lb].back();
                }
                local_edg_count += 1;
            }
            // project the nodehead_bad onto the good mesh outside current FL_good region and FL_good will be dynamically
        }
        
    }
}
}