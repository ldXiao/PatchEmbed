#include "Match_Maker_Loop.h"
#include "Match_Maker_Tree.h"
#include "Kruskal.h"
namespace bcclean{
namespace MatchMaker{
    using json = nlohmann::json;

    void _build_dual_frame_graph(
        const std::vector<bcclean::edge> edge_list,
        std::vector<std::pair<int, std::pair<int, int> > > & dual_frame_graph
    )
    {
        dual_frame_graph.clear();
        int edge_idx = 0;
        for(auto edge: edge_list)
        {
            dual_frame_graph.push_back(std::make_pair(edge_idx, edge._label_pair));
            edge_idx++;
        }
        return;
    }

    void trace_for_edge_loop(
        const Eigen::MatrixXd & V_bad,
        const Eigen::MatrixXi & F_bad,
        const std::vector<edge> & edge_list,
        const std::vector<int> & node_list_good,
        const std::map<int, int> & node_image_dict,
        const std::map<int, std::vector<int> > & node_edge_dict,
        const int edge_idx,
        Eigen::MatrixXd & V_good,
        Eigen::MatrixXi & F_good,
        std::vector<std::vector<int> > & VV_good,
        std::vector<std::vector<int> > & VEdges_good,
        std::vector<std::vector<int> > & TEdges_good,
        std::vector<int> & total_silence_list,
        std::map<int, std::map<int, bool> > & node_edge_visit_dict,
        std::map<int, std::vector<int> > & edge_path_map,
        json & path_json
    )
    {
        
        Eigen::MatrixXi TT_good;
        std::vector<std::vector<int> > VF_good;
        igl::triangle_triangle_adjacency(F_good, TT_good);
        {
            std::vector<std::vector<int> > VFi_good;
            igl::vertex_triangle_adjacency(V_good, F_good, VF_good, VFi_good);
        }
        edge edg = edge_list[edge_idx];
        int target =-1;
        int target_bad = -1;
        int source_bad = -1;
        int source = node_image_dict.at(edg.head);     
        source_bad = edg.head;
        target = node_image_dict.at(edg.tail);
        target_bad = edg.tail;
        assert(target != -1 && source != -1); 
        std::vector<double> Weights;
        // update local sector
        // (a) create CC_node_face_list
        std::map<int, std::vector<int> > CC_node_face_dict;
        std::vector<std::vector<int > > VV_temp;
        CC_faces_per_node(V_good, F_good, {source, target}, CC_node_face_dict);
        {
            std::map<int, std::map<int, bool> > node_edge_visit_dict_temp;
            std::map<int, std::vector<int> > node_edge_dict_temp;
            for(auto item: node_edge_dict)
            {
                node_edge_dict_temp[node_image_dict.at(item.first)] = node_edge_dict.at(item.first);
                node_edge_visit_dict_temp[node_image_dict.at(item.first)] = node_edge_visit_dict.at(item.first);
            }
            update_local_sector(
                VV_good, 
                F_good, 
                node_edge_visit_dict_temp, 
                node_edge_dict_temp,
                TEdges_good,
                CC_node_face_dict,
                source,
                target,
                edge_idx,
                VV_temp);
        }
        setWeights(V_good, V_bad, edg, 10, 1,  Weights);
        // the Weights is vertex based

        // dijkstra_trace(....,VEdges, TEdges);
        std::vector<int> path;
        dijkstra_trace(VV_temp, source, target, Weights, path);
        std::vector<int> path_records(path.size()-2);
        std::printf("for mst edge %d, find a path:\n",edge_idx);
        for(auto rec : path)
        {
            std::cout << rec<<", ";
        }
        std::cout << "\n";
        for(int p =0 ; p < path.size()-2; ++p)
        {
            path_records[p] = path[p+1];
        }
        // path update VV
        silence_vertices(VV_good,path_records);
        for(auto rec: path_records)
        {
            total_silence_list.push_back(rec);
        }
        //path update VEdges
        for(auto vidx:path_records){
            VEdges_good[vidx].push_back(edge_idx);
        }

        // path updates TEdges
        // setf the triangle edges in cuts to be true
        for(int rc_idx=0; rc_idx < path.size()-1; ++rc_idx){
            int uidx = path[rc_idx];
            int vidx = path[(rc_idx+1)%path.size()];
            std::vector<int> inter(VF_good[uidx].size()+ VF_good[vidx].size());
            auto it = std::set_intersection(VF_good[uidx].begin(), VF_good[uidx].end(), VF_good[vidx].begin(), VF_good[vidx].end(), inter.begin());
            inter.resize(it-inter.begin());
            // there should be only one comman adjacent triangle for boundary vertices
            for(auto trg: inter){
                for(int edgpos =0; edgpos < 3 ; ++edgpos){
                    int uuidx = F_good(trg, edgpos);
                    int vvidx = F_good(trg, (edgpos+1)% 3);
                    if(uuidx == uidx && vvidx == vidx){
                        TEdges_good[trg][edgpos] = edge_idx;//1457-1459 1587 1588 2189 2371 2373 2374 2716 2742 1794 2712 3046 3047
                        // break;
                    }
                    if(uuidx == vidx && vvidx == uidx){
                        TEdges_good[trg][edgpos] = edge_idx;
                        // break;
                    }
                }
            }
        }

        // split_detect
        std::pair<int, int> split;
        while(split_detect(F_good, TT_good, node_list_good,VEdges_good, TEdges_good, split))
        {
            

            // splits_update
            splits_update(split, V_good, F_good, VEdges_good, TEdges_good, VV_good);
            igl::triangle_triangle_adjacency(F_good, TT_good);
        }
        silence_vertices(VV_good, total_silence_list);

        edge_path_map[edge_idx] = path;

        path_json[std::to_string(edge_idx)] = path;   
        
        std::ofstream file;
        file.open("../debug_paths.json");
        file << path_json;

        igl::writeOBJ("../debug_mesh.obj", V_good, F_good);

        // update visit_dict or loop condition update
        node_edge_visit_dict[target_bad][edge_idx]=true;
        node_edge_visit_dict[source_bad][edge_idx]=true;
    
    }

    void trace_and_label_loop(
        const Eigen::MatrixXd & V_bad,
        const Eigen::MatrixXi & F_bad,
        const Eigen::VectorXi & FL_bad,
        Eigen::MatrixXd & V_good,
        Eigen::MatrixXi & F_good,
        Eigen::VectorXi & FL_good,
        int kk
    )
    {
        // datas dump to file for debug
        std::map<int, int> edge_order_map; // store and maintain the order of added edges {order: edge_dx}
        std::map<int, std::vector<int> > edge_path_map; // {edge_idx, path}
        // Randomize Seed
        srand(static_cast<unsigned int>(time(nullptr)));
        int total_label_num = FL_bad.maxCoeff()+1;
        // PART 0 GET THE FRAME ON BAD MESH
        std::vector<bcclean::edge> edge_list;
        std::unordered_map<int, std::vector<int> > patch_edge_dict;
        std::unordered_map<int, std::vector<bool> > patch_edge_direction_dict;
        build_edge_list_loop(V_bad, F_bad, FL_bad, total_label_num, edge_list, patch_edge_dict, patch_edge_direction_dict);
        FL_good = Eigen::VectorXi::Constant(F_good.rows(), -1);
        std::vector<std::pair<int, std::pair<int, int> > > dual_frame_graph;
        _build_dual_frame_graph(edge_list, dual_frame_graph);
        std::vector<std::pair<int, std::pair<int, int> > > dual_frame_MST
         = Algo::Kruskal_MST(dual_frame_graph);
        std::vector<int> patch_order = Algo::MST_BFS(dual_frame_MST);
        


        for(int patch_idx: patch_order)
        {
            
        }
    }
}
}