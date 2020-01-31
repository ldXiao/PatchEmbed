#include "Match_Maker_Loop.h"
#include "Match_Maker_Tree.h"
#include "Kruskal.h"
#include "polyline_distance.h"
#include "loop_colorize.h"
#include "params.h"
#include <igl/Timer.h>
#include <igl/facet_components.h>
#include <climits>
#include <igl/bounding_box_diagonal.h>
namespace bcclean{
namespace MatchMaker{
    using json = nlohmann::json;

    bool check_problematic_patch(
        const Eigen::MatrixXd & V,
        const std::unordered_map<int, std::vector<int> > & patch_edge_dict,
        const std::map<int, std::vector<int> > & patch_node_dict,
        const std::vector<edge> & edge_list,
        const int patch_idx,
        const double & bound
    )
    {
        double mindis = std::numeric_limits<double>::max();
        std::vector<int> node_list_p = patch_node_dict.at(patch_idx);
        std::vector<int> edge_list_p = patch_edge_dict.at(patch_idx);
        // case 0 if there is only 2 nodes in this patch(loop)
        // we directly compare the hausdorff distance between the two edges
        assert(edge_list_p.size()>1);
        if(edge_list_p.size()==2)
        {
            int edgidxA = edge_list_p.at(0);
            int edgidxB = edge_list_p.at(1);
            std::vector<int> pathA = edge_list.at(edgidxA)._edge_vertices;
            std::vector<int> pathB = edge_list.at(edgidxB)._edge_vertices;
            mindis = Eval::hausdorff1d(V, pathA, V, pathB);
            if(mindis < bound)
            {
                return true;
            }
        }
        else{
            
            for(auto ndidx: node_list_p)
            {
                std::vector<int> nadj_edges; // nonadjacent_edges of the node
                for(auto edgidx:edge_list_p)
                {
                    const edge &  edg = edge_list.at(edgidx);
                    if(edg.head != ndidx && edg.tail  != ndidx)
                    {
                        mindis = std::min(mindis, Eval::single_sample_trial(V, edg._edge_vertices, V, ndidx));
                    }
                }
            }
            if(mindis < bound)
            {
                return true;
            }
            else return false;

        }
    }


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
        Eigen::VectorXi & FL_good,
        std::vector<std::vector<int> > & VV_good,
        std::vector<std::vector<int> > & VEdges_good,
        std::vector<std::vector<int> > & TEdges_good,
        std::vector<int> & total_silence_list,
        std::map<int, std::map<int, bool> > & node_edge_visit_dict,
        std::map<int, std::vector<int> > & edge_path_map,
        json & path_json,
        bool debug
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
                if(node_image_dict.find(item.first)!=node_image_dict.end()){
                    node_edge_dict_temp[node_image_dict.at(item.first)] = node_edge_dict.at(item.first);
                    node_edge_visit_dict_temp[node_image_dict.at(item.first)] = node_edge_visit_dict.at(item.first);
                }
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
        if(debug)
        {
            Eigen::VectorXd source_target=Eigen::VectorXd::Constant(6,0);
            for(int xx: {0,1,2})
            {
                source_target(xx)=V_good(source,xx);
                source_target(xx+3)=V_good(target,xx);
            }
            
            igl::writeDMAT("../dbginfo/source_target.dmat",source_target);
        }
        assert(path.size()>=2);
        if(path.size()<2)
        {
            std::cout << "for path" << edge_idx << std::endl;
            std::cout << "start" << source  << "target" << target << std::endl;
            std::cout << "path of size only "<< path.size() << std::endl;
            exit(EXIT_FAILURE);
        }
        std::vector<int> path_records(path.size()-2);
        std::printf("for edge %d, find a path:\n",edge_idx);
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
        igl::triangle_triangle_adjacency(F_good, TT_good);
        while(split_detect(F_good, TT_good, node_list_good,VEdges_good, TEdges_good, split))
        {
            // splits_update
            splits_update(split, V_good, F_good, FL_good, VEdges_good, TEdges_good, VV_good);
            igl::triangle_triangle_adjacency(F_good, TT_good);
        }
        silence_vertices(VV_good, total_silence_list);

        edge_path_map[edge_idx] = path;
        if(debug)
        {
            path_json[std::to_string(edge_idx)] = path;   
            std::ofstream file;
            file.open("../dbginfo/debug_paths.json");
            file << path_json;
            igl::writeOBJ("../dbginfo/debug_mesh.obj", V_good, F_good);
        }
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
        const params params
    )
    {
        // datas dump to file for debug
        std::map<int, int> edge_order_map; // store and maintain the order of added edges {order: edge_dx}
        std::map<int, std::vector<int> > edge_path_map; // {edge_idx, path}
        int total_label_num = FL_bad.maxCoeff()+1;
        // PART 0 GET THE FRAME ON BAD MESH
        std::vector<bcclean::edge> edge_list;
        std::unordered_map<int, std::vector<int> > patch_edge_dict;
        std::unordered_map<int, std::vector<bool> > patch_edge_direction_dict;
        if(params.debug)
        {
            igl::writeOBJ("../dbginfo/debug_mesh_bad.obj", V_bad, F_bad);
            igl::writeDMAT("../dbginfo/FL_bad.dmat", FL_bad);
        }
        int largest_patch;
        build_edge_list_loop(V_bad, F_bad, FL_bad, total_label_num, edge_list, patch_edge_dict, patch_edge_direction_dict,largest_patch);
        {
            json path_json_bad;
            int edg_idx =0;
            for(auto edg: edge_list)
            {
                path_json_bad[std::to_string(edg_idx)]= edg._edge_vertices;
                edg_idx += 1;
            }
            std::ofstream file;
            file.open("../dbginfo/debug_path_bad.json");
            file << path_json_bad;
        }
        FL_good = Eigen::VectorXi::Constant(F_good.rows(), -1);
        std::vector<std::pair<int, std::pair<int, int> > > dual_frame_graph;
        _build_dual_frame_graph(edge_list, dual_frame_graph);
        // std::vector<std::pair<int, std::pair<int, int> > > dual_frame_MST
        //  = Algo::Kruskal_MST(dual_frame_graph);
        // std::vector<int> patch_order = Algo::MST_BFS(dual_frame_MST);
        std::vector<int> patch_order = Algo::Graph_BFS(dual_frame_graph, largest_patch);
        std::map<int, std::vector<int> > patch_node_dict;
        std::vector<int> problematic_patches;
        {
        for(auto item: patch_edge_dict)
        {
            int ptidx = item.first;
            patch_node_dict[ptidx] = std::vector<int>();
            for(auto edgidx: patch_edge_dict[ptidx])
            {
                int head = edge_list[edgidx].head;
                int tail = edge_list[edgidx].tail;
                if(std::find(patch_node_dict[ptidx].begin(), patch_node_dict[ptidx].end(), head)== patch_node_dict[ptidx].end())
                {
                    patch_node_dict[ptidx].push_back(head);
                }
                if(std::find(patch_node_dict[ptidx].begin(), patch_node_dict[ptidx].end(), tail)== patch_node_dict[ptidx].end())
                {
                    patch_node_dict[ptidx].push_back(tail);
                }
            }
        }
        
        }
        
        double bound = params.guard_len_r * igl::bounding_box_diagonal(V_bad);
        for(auto item: patch_edge_dict)
        {
            int ptidx = item.first;
            if(check_problematic_patch(V_bad, patch_edge_dict, patch_node_dict, edge_list, ptidx, bound)){
                problematic_patches.push_back(ptidx);
            }
        }

        std::vector<int> patch_order_adv = Algo::Constriant_Graph_BFS(dual_frame_graph, patch_node_dict, problematic_patches, largest_patch);
        // we can assume that  all nodes have valance more than or equal to 3


        std::vector<int> node_list_bad;
        _gen_node_list(F_bad, FL_bad, total_label_num, node_list_bad);

        std::map<int, std::vector<int> > node_edge_dict;
        _gen_node_CCedges_dict(V_bad, F_bad, edge_list, node_list_bad, node_edge_dict);


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
        // start with the nodes with largest valance and deal with the edge starting with this node in counter clock order
        std::vector<int> total_silence_list; // store only the vertices on the path interior (head tail excluded)


        json path_json;

        std::map<int, int> node_image_dict;
        std::vector<int> node_list_good; // potnetially replace all use with node_image_dict
        std::map<int , std::map<int, bool> > node_edge_visit_dict;
        for(auto nd: node_list_bad)
        {
            for(auto q: node_edge_dict[nd])
            {
                node_edge_visit_dict[nd][q]=false;
            }
        }
        std::map<int, bool> edge_visit_dict;
        int kkk = 0;
        for(auto edg: edge_list )
        {
            edge_visit_dict[kkk]= false;
            kkk++;
        }

        for(auto patch_idx : patch_order_adv)
        {
            if(params.debug)
            {
                igl::writeDMAT("../dbginfo/cur_patch.dmat", Eigen::VectorXi::Constant(1,patch_idx));
            }
            for(auto edge_idx: patch_edge_dict[patch_idx])
            {
                if(edge_visit_dict[edge_idx])
                {
                    continue;
                }
                edge_visit_dict[edge_idx]=true;
                int source_bad = edge_list.at(edge_idx).head;
                int target_bad = edge_list.at(edge_idx).tail;
                if(params.debug)
                {
                    Eigen::VectorXd source_target_bad=Eigen::VectorXd::Constant(6,0);
                    for(int xx: {0,1,2})
                    {
                        source_target_bad(xx)=V_bad(source_bad,xx);
                        source_target_bad(xx+3)=V_bad(target_bad,xx);
                    }
                    
                    igl::writeDMAT("../dbginfo/source_target_bad.dmat",source_target_bad);
                }
                // find the corresponding source and target on goodmesh
                // project the source_bad and target_bad onto the good_mesh outside any colored patch
                int source = -1;
                int target = -1;

                Eigen::MatrixXi TT_good;
                igl::triangle_triangle_adjacency(F_good, TT_good);
                //update the node_list_good and node_image_dict;

                if(node_image_dict.find(source_bad)==node_image_dict.end())
                {
                    proj_node_loop(
                        V_bad,
                        F_bad,
                        source_bad,
                        node_list_good,
                        TT_good,
                        VV_good,
                        TEdges_good,
                        VEdges_good,
                        V_good,
                        F_good,
                        FL_good,
                        source
                    );
                    assert(source!= -1);
                    node_image_dict[source_bad]=source;
                    node_list_good.push_back(source);
                    igl::triangle_triangle_adjacency(F_good, TT_good);
                    std::pair<int, int> split;
                    while(split_detect(F_good, TT_good, node_list_good,VEdges_good, TEdges_good, split))
                    {
                        // splits_update
                        splits_update(split, V_good, F_good, FL_good, VEdges_good, TEdges_good, VV_good);
                        igl::triangle_triangle_adjacency(F_good, TT_good);
                    }
                }
                if(node_image_dict.find(target_bad)==node_image_dict.end())
                {
                    proj_node_loop(
                        V_bad,
                        F_bad,
                        target_bad,
                        node_list_good,
                        TT_good,
                        VV_good,
                        TEdges_good,
                        VEdges_good,
                        V_good,
                        F_good,
                        FL_good,
                        target
                    );
                    assert(target!= -1);
                    node_image_dict[target_bad]=target;
                    node_list_good.push_back(target);
                    igl::triangle_triangle_adjacency(F_good, TT_good);
                    std::pair<int, int> split;
                    while(split_detect(F_good, TT_good, node_list_good,VEdges_good, TEdges_good, split))
                    {
                        // splits_update
                        splits_update(split, V_good, F_good, FL_good, VEdges_good, TEdges_good, VV_good);
                        igl::triangle_triangle_adjacency(F_good, TT_good);
                    }
                }
                
                
                //
                trace_for_edge_loop(
                    V_bad,
                    F_bad, 
                    edge_list, 
                    node_list_good, 
                    node_image_dict, 
                    node_edge_dict, 
                    edge_idx, 
                    V_good, 
                    F_good, 
                    FL_good,
                    VV_good, 
                    VEdges_good, 
                    TEdges_good, 
                    total_silence_list,
                    node_edge_visit_dict, 
                    edge_path_map, 
                    path_json,
                    params.debug);
                
            }
            // one loop patch finished
            // colorize FL_good with label patch_idx
            int stem_edge = patch_edge_dict[patch_idx][0];
            int v0, v1;
            v0 = edge_path_map[stem_edge][0];
            v1 = edge_path_map[stem_edge][1];
            std::vector<std::vector<int> > VF_good;
            {
                std::vector<std::vector<int> > VFi_good;
                igl::vertex_triangle_adjacency(V_good, F_good, VF_good, VFi_good);
            }
            bool directionCC = patch_edge_direction_dict[patch_idx][0];
            std::vector<int> inter(VF_good[v0].size()+ VF_good[v1].size());
            auto it = std::set_intersection(VF_good[v0].begin(), VF_good[v0].end(), VF_good[v1].begin(), VF_good[v1].end(), inter.begin());
            inter.resize(it-inter.begin());
            assert(inter.size()==2);
            int ffa = inter[0];
            int ffb = inter[1];
            bool ffa_coline = false; 
            for(int j: {0,1,2})
            {
                if(F_good(ffa,j)==v0 && F_good(ffa,(j+1)% 3)==v1)
                {
                    ffa_coline=true;
                    break;
                }
            }
            int ff_in;
            if(ffa_coline != directionCC)
            {
                ff_in = ffa;
            }
            else
            {
                ff_in = ffb;
            }
            loop_colorize(V_good, F_good, TEdges_good, ff_in, patch_idx, FL_good);
            if(params.debug)
            {
                igl::writeDMAT("../dbginfo/FL_loop.dmat", FL_good);
            }
        }
    }

}
}