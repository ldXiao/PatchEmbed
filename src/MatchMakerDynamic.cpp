#include "Match_Maker_Tree.h"
#include "Kruskal.h"
#include "backtrack_diff.h"
#include "polyline_distance.h"
#include "loop_colorize.h"
#include "params.h"
#include <igl/Timer.h>
#include <climits>
#include <list>
#include <deque>
#include <igl/bounding_box_diagonal.h>
#include <igl/barycentric_to_global.h>
#include <spdlog/common.h>
#include <algorithm>
#include "Edge_Dijkstra.h"
#include "CellularGraph.h"
#include "MatchMakerDynamic.h"
#include "TraceComplex.h"
#include "TransferCellGraph.h"
#include "Patch_Bijection.h"
/* BTCMM means backtracking cellular matchmaker */
namespace bcclean{
namespace MatchMaker{
    using json = nlohmann::json;

 

    void build_dual_frame_graph(
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


    bool BTCMM1_for_edge(
        const CellularGraph & cg,
        TraceComplex & tc,
        const int edge_idx,
        std::map<int, std::map<int, bool> > & node_edge_visit_dict,
        json & path_json,
        const params param
    )
    {
        
        edge edg = cg._edge_list[edge_idx];
        int target =-1;
        int target_bad = -1;
        int source_bad = -1;
        int source = tc._node_image_map.at(edg.head);     
        source_bad = edg.head;
        target = tc._node_image_map.at(edg.tail);
        target_bad = edg.tail;
        assert(target != -1 && source != -1); 
        // update local sector
        // (a) create CC_node_face_list
        std::map<int, std::vector<int> > CC_node_face_dict;
        std::vector<std::vector<int > > VV_temp,VV_temp1;
        CCfaces_per_node(tc._V,tc._F, tc._VF,{source,target},CC_node_face_dict);
        {
            std::map<int, std::map<int, bool> > node_edge_visit_dict_temp;
            std::map<int, std::vector<int> > node_edge_dict_temp;
            for(auto item: cg._node_edge_dict)
            {
                if(tc._node_image_map.find(item.first)!=tc._node_image_map.end()){
                    node_edge_dict_temp[tc._node_image_map.at(item.first)] = cg._node_edge_dict.at(item.first);
                    node_edge_visit_dict_temp[tc._node_image_map.at(item.first)] = node_edge_visit_dict.at(item.first);
                }
            }
            update_local_sector(
                tc._VV,
                tc._F,
                node_edge_visit_dict_temp, 
                node_edge_dict_temp,
                tc._TEdges,
                CC_node_face_dict,
                source,
                target,
                edge_idx,
                VV_temp1
            );
        }
        // the Weights is edge based
        Eigen::SparseMatrix<double> SpWeight;
        Trace::setWeight1(tc._V, tc._F, cg._vertices,edg,SpWeight);
        std::vector<int> path;
        Trace::Edge_Dijkstra(VV_temp1, target, source, SpWeight,path);// because the path will be inverted, we switch the position of source and target
        // dijkstra_trace(VV_temp, source, target, Weights, path);
        if(param.debug)
        {
            Eigen::VectorXd source_target=Eigen::VectorXd::Constant(6,0);
            for(int xx: {0,1,2})
            {
                source_target(xx)=tc._V(source,xx);
                source_target(xx+3)=tc._V(target,xx);
            }
            
            igl::writeDMAT(param.data_root+"/source_target.dmat",source_target);
        }
        // assert(path.size()>=2);
        if(path.size()<2)
        {
            std::cout << "for path" << edge_idx << std::endl;
            std::cout << "start" << source  << "target" << target << std::endl;
            std::cout << "path of size only "<< path.size() << std::endl;
            return false;
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
        silence_vertices(tc._VV, path_records);
        for(auto rec: path_records)
        {
            tc._total_silence_list.push_back(rec);
        }
        //path update VEdges
        for(auto vidx:path_records){
            tc._VEdges[vidx].push_back(edge_idx);
        }

        // path updates TEdges
        // set the triangle edges in cuts to be true
        for(int rc_idx=0; rc_idx < path.size()-1; ++rc_idx){
            int uidx = path[rc_idx];
            int vidx = path[(rc_idx+1)%path.size()];
            std::vector<int> vtr = tc._VF[vidx];
            std::vector<int> utr = tc._VF[uidx];
            std::sort(vtr.begin(), vtr.end());
            std::sort(utr.begin(), utr.end());
            std::vector<int> inter(vtr.size()+utr.size());
            auto it = std::set_intersection(vtr.begin(), vtr.end(), utr.begin(), utr.end(), inter.begin());
            inter.resize(it-inter.begin());
            // there should be only one comman adjacent triangle for boundary vertices
            for(auto trg: inter){
                for(int edgpos =0; edgpos < 3 ; ++edgpos){
                    int uuidx = tc._F(trg, edgpos);
                    int vvidx = tc._F(trg, (edgpos+1)% 3);
                    if(uuidx == uidx && vvidx == vidx){
                        tc._TEdges[trg][edgpos] = edge_idx;
                    }
                    if(uuidx == vidx && vvidx == uidx){
                        tc._TEdges[trg][edgpos] = edge_idx;
                    }
                }
            }
        }

        // splits_detect

        std::vector<std::pair<int, int> > splits;
        tc.splits_detect(splits);
        for(auto split : splits)
        {
            // splits_update
            tc.split_update(split);
        }
        silence_vertices(tc._VV, tc._total_silence_list);
        tc._edge_path_map[edge_idx] = path;
        if(param.debug)
        {
            path_json[std::to_string(edge_idx)] = path;   
            std::ofstream file;
            file.open(param.data_root+"/debug_paths.json");
            file << path_json;
            igl::writeOBJ(param.data_root+"/debug_mesh.obj", tc._V, tc._F);
            igl::writeDMAT(param.data_root+"/FL_loop.dmat", tc._FL);
        }
        // update visit_dict or loop condition update
        node_edge_visit_dict[target_bad][edge_idx]=true;
        node_edge_visit_dict[source_bad][edge_idx]=true;
        return true;
    }






    bool BTCMM1(
        const CellularGraph & cg,
        Eigen::MatrixXd & V_good,
        Eigen::MatrixXi & F_good,
        Eigen::VectorXi & FL_good,
        const params param
    )
    {


        // datas dump to file for debug
        std::map<int, int> edge_order_map; // store and maintain the order of added edges {order: edge_dx}
        int total_label_num = cg.label_num;
        // PART 0 GET THE FRAME ON BAD MESH
        std::vector<bcclean::edge> edge_list = cg._edge_list;
        // std::unordered_map<int, std::vector<int> > patch_edge_dict;
        // std::unordered_map<int, std::vector<bool> > patch_edge_direction_dict;
        if(param.debug)
        {
            igl::writeOBJ(param.data_root+"/debug_mesh_bad.obj", cg.V, cg.F);
            igl::writeDMAT(param.data_root+"/FL_bad.dmat",cg.FL);
        }
        int largest_patch = cg.root_cell;
        // build_edge_list_loop(V_bad, F_bad, FL_bad, total_label_num, edge_list, patch_edge_dict, patch_edge_direction_dict,largest_patch);
        {
            json path_json_bad;
            int edg_idx =0;
            for(auto edg: cg._edge_list)
            {   


                std::vector<int> path_raw;
                for(auto vidx: edg._edge_vertices)
                {
                    path_raw.push_back(cg._ivmap.at(vidx));
                }
                path_json_bad[std::to_string(edg_idx)]= path_raw;
                edg_idx += 1;
            }
            std::ofstream file;
            file.open(param.data_root+"/debug_path_bad.json");
            file << path_json_bad;
        }
        std::vector<std::pair<int, std::pair<int, int> > > dual_frame_graph;

        build_dual_frame_graph(cg._edge_list, dual_frame_graph);

        std::vector<int> patch_order = Algo::Graph_BFS(dual_frame_graph, largest_patch);


        TraceComplex tc;
        tc.initialize(V_good,F_good);

        json path_json;


        std::map<int, std::map<int, bool> >  node_edge_visit_dict;
        for(auto nd: cg._nodes)
        {
            for(auto q: cg._node_edge_dict.at(nd))
            {
                node_edge_visit_dict[nd][q]=false;
            }
        }
        std::map<int, bool> edge_visit_dict;
        int kkk = 0;
        for(auto edg: cg._edge_list)
        {
            edge_visit_dict[kkk]= false;
            kkk++;
        }
        std::list<int> patch_queue(patch_order.begin(), patch_order.end());
        std::list<int> recycle;
        int switch_count  = 0;
        double bcthreshold = param.backtrack_threshold;
        double arthreshold = param.area_threshold;
        while(! patch_queue.empty())
        {
            bool curpatch_succ = true;  
            /* copy part */
            TraceComplex tc_copy = tc;
            
            


            std::map<int, bool> edge_visit_dict_copy = edge_visit_dict;
            std::map<int , std::map<int, bool> > node_edge_visit_dict_copy = node_edge_visit_dict; 
            
             
            json path_json_copy = path_json;

            // copy everything in advance , if backtrack happens 
            // replace the origin object with the copies



            int patch_idx = patch_queue.front();
            patch_queue.pop_front();

            if(param.debug)
            {
                igl::writeDMAT(param.data_root+"/cur_patch.dmat", Eigen::VectorXi::Constant(1,patch_idx));

            }
            bool all_edge_traced = true;
            for(auto edge_idx: cg._patch_edge_dict.at(patch_idx))
            {
                if(edge_visit_dict[edge_idx])
                {
                    continue;
                }
                all_edge_traced = false;

                edge_visit_dict[edge_idx] = true;
                int source_bad = cg._edge_list.at(edge_idx).head;
                int target_bad = cg._edge_list.at(edge_idx).tail;
                if(param.debug)
                {
                    Eigen::VectorXd source_target_bad = Eigen::VectorXd::Constant(6,0);
                    for(int xx: {0,1,2})
                    {
                        source_target_bad(xx)=cg._vertices[source_bad](xx);
                        source_target_bad(xx+3)=cg._vertices[target_bad](xx);
                    }
                    
                    igl::writeDMAT(param.data_root+"/source_target_bad.dmat",source_target_bad);
                }
                //
                //
                int source = -1;
                int target  = -1;

                
                
                if(tc._node_image_map.find(source_bad)==tc._node_image_map.end())
                {
                    proj_node_loop(cg, source_bad, tc, source); 
                    // updated the tc._node_listd and tc._node_image_map;
                    
                    assert(source!= -1);
                    
                    

                    std::vector<std::pair<int, int> > splits;
                    tc.splits_detect(splits);
                    for(auto split : splits)
                    {
                        // splits_update
                        tc.split_update(split);
                    }
                    
                   
                }
                if(tc._node_image_map.find(target_bad)==tc._node_image_map.end())
                {
        
                    proj_node_loop(cg, target_bad, tc, target);
                    assert(target!= -1);
                    std::vector<std::pair<int, int> > splits;
                    tc.splits_detect(splits);
                    for(auto split : splits)
                    {
                        // splits_update
                        tc.split_update(split);
                    }

                }




                //
                if(!BTCMM1_for_edge(
                    cg,
                    tc,
                    edge_idx, 
                    node_edge_visit_dict, 
                    path_json,
                    param)) 
                    {
                        // roll back current patch
                        curpatch_succ = false;
                        break;
                    }
            
            }
            
            bool withinthreshold= false; 
            if(curpatch_succ){
                // if the current is traced successfully, we go on to compare the difference
                int ff_in = _locate_seed_face(cg, tc, patch_idx);
                double newarea = loop_colorize(tc._V, tc._F, tc._TEdges, ff_in, patch_idx, tc._FL);
                double target_area = cg._patch_area_dict.at(patch_idx);
                bool area_withinthreshold = (std::abs(newarea/target_area) < arthreshold)|| (switch_count > 2);
                withinthreshold=backtrack_diff(
                    tc._V,
                    cg,
                    patch_idx,
                    tc._edge_path_map,
                    bcthreshold
                ) && area_withinthreshold;
            }
            if(!withinthreshold && !(all_edge_traced))
            {

                // if all_edge_traced, we don't care about the thresholds
                // the traced patch error is larget than the backtrack_threshold
                // abort the result in this loop
                recycle.push_back(patch_idx);
                // reverse copy everything
                std::cout << "patch" << patch_idx << "postponed" << std::endl;

                /* copy part*/
                /* copy part */
                tc = tc_copy;
                edge_visit_dict = edge_visit_dict_copy;
                node_edge_visit_dict = node_edge_visit_dict_copy;
                path_json = path_json_copy;
            }
            else
            {
                if(param.debug)
                {
                    igl::writeOBJ(param.data_root+"/debug_mesh.obj", tc._V, tc._F);
                    igl::writeDMAT(param.data_root+"/FL_loop.dmat", tc._FL);
                }

                // loop over patch_queue and  recycle to find patches that where all edges has been traced
                // and place them at the front of patch_queue
                std::vector<int> q_relocate_list;
                std::vector<int> r_relocate_list;
                for(auto pidx: patch_queue)
                {
                    bool p_finished = true;
                    for(auto eidx: cg._patch_edge_dict.at(pidx))
                    {
                        if(!edge_visit_dict.at(eidx))
                        {
                            p_finished = false;
                            break;
                        }
                    }
                    if(p_finished)
                    {
                        q_relocate_list.push_back(pidx);
                    }
                }
               for(auto pidx: recycle)
                {
                    bool p_finished = true;
                    for(auto eidx: cg._patch_edge_dict.at(pidx))
                    {
                        if(!edge_visit_dict.at(eidx))
                        {
                            p_finished = false;
                            break;
                        }
                    }
                    if(p_finished)
                    {
                        r_relocate_list.push_back(pidx);
                    }
                } 
                for(auto pidx: q_relocate_list)
                {
                    patch_queue.remove(pidx);
                }
                for(auto pidx: r_relocate_list)
                {
                    recycle.remove(pidx);
                }
                for(auto pidx: q_relocate_list)
                {
                    patch_queue.push_front(pidx);
                }
                for(auto pidx: r_relocate_list)
                {
                    patch_queue.push_front(pidx);
                }

            }
            if(patch_queue.empty() && !recycle.empty())
            {
                // if it is still empty after relocation 
                // switch recycle and patch_queue
                if(switch_count < 8)
                {
                    std::list<int> temp = recycle;
                    recycle = patch_queue;
                    patch_queue = temp;
                    bcthreshold = 1.5 * bcthreshold;
                    arthreshold = 1.5 * arthreshold;
                    switch_count +=1;
                }
            }

        }
        F_good = tc._F;
        V_good = tc._V;
        FL_good = tc._FL;
        if(param.debug)
        { 
            std::ofstream file;
            file.open(param.data_root+"/debug_paths.json");
            file << path_json;
            igl::writeOBJ(param.data_root+"/debug_mesh.obj", tc._V, tc._F);
            igl::writeDMAT(param.data_root+"/FL_loop.dmat", tc._FL);
            std::ofstream split_file;
            split_file.open(param.data_root+"/splits_record.txt");
            for(auto & vec: tc._splits_record)
            {
                if(vec.size()==9)
                {
                    // it is splits
                    split_file << 0;
                    split_file << " ";
                }
                else
                {
                    // it is an insert
                    split_file << 1;
                    split_file << " ";
                }
                
                for(auto & val: vec)
                {
                    split_file << val;
                    split_file << " ";
                }
                split_file << "\n";
            }
        }
        if(recycle.empty()){
            CellularGraph cgt;
            Bijection::TransferCellGraph(cg, tc, cgt);
            Eigen::MatrixXd M_s2t;
            Bijection::BijGlobal(cg,cgt, M_s2t);
            Eigen::MatrixXd Vmap=igl::barycentric_to_global(cgt.V, cgt.F, M_s2t);
            igl::writeOBJ(param.data_root+"/map.obj", Vmap, cg.F);
            return true;
        }
        else return false;

    }

}
}