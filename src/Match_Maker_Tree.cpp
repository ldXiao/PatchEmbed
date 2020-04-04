#include "proj_node.h"
#include "loop_colorize.h"
#include "Kruskal.h"
#include <igl/slice.h>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#include <utility>
#include <nlohmann/json.hpp>
#include "Match_Maker_Tree.h"
#include "params.h"
#include "CellularGraph.h"
#include "TraceComplex.h"
#include "Edge_Dijkstra.h"
#include "Patch_Bijection.h"
#include "TransferCellGraph.h"
#include <igl/barycentric_to_global.h>
namespace bcclean {
namespace MatchMaker{
    using json = nlohmann::json;
    int _locate_seed_face(
        const CellularGraph & cg, 
        const TraceComplex & tc, 
        const int patch_idx)
    {
        // one loop patch finished
        // colorize tc._FL with label patch_idx
        int stem_edge = cg._patch_edge_dict.at(patch_idx)[0];
        int v0, v1;
        v0 = tc._edge_path_map.at(stem_edge)[0];
        v1 = tc._edge_path_map.at(stem_edge)[1];

        bool directionCC = cg._patch_edge_direction_dict.at(patch_idx)[0];
        
        std::vector<int> vtr0 = tc._VF[v0];
        std::vector<int> vtr1 = tc._VF[v1];
        std::sort(vtr0.begin(), vtr0.end());
        std::sort(vtr1.begin(), vtr1.end());
        std::vector<int> inter(vtr0.size()+vtr1.size());
        auto it = std::set_intersection(vtr0.begin(), vtr0.end(), vtr1.begin(), vtr1.end(), inter.begin());
        inter.resize(it-inter.begin());
        assert(inter.size()==2);
        int ffa = inter[0];
        int ffb = inter[1];
        bool ffa_coline = false; 
        for(int j: {0,1,2})
        {
            if(tc._F(ffa,j)==v0 && tc._F(ffa,(j+1)% 3)==v1)
            {
                ffa_coline=true;
                break;
            }
        }
        int ff_in;
        if(ffa_coline == directionCC)
        {
            ff_in = ffa;
        }
        else
        {
            ff_in = ffb;
        }
        return ff_in;
    }

    void _build_frame_graph(
        const std::vector<bcclean::edge> & edge_list,
        std::vector<std::pair<int,std::pair<int,int> > > & frame_graph
    )
    {
        /*
            modify in place a frame_graph of edges, with identical index with edge_list and head and tail
            will be used for generating a min spanning tree
         */

        frame_graph.clear();
        int count = 0;
        for(auto edg: edge_list)
        {
            frame_graph.push_back(std::make_pair(count, std::make_pair(edg.head, edg.tail)));
            count++;
        }
        return;
    }

    void sector_silence_vertices(std::vector<std::vector<int> > & VV, const int source, const std::vector<int> & silent_indices)
    {
        // only remove the connection between source and silent_indices
        for(auto index : silent_indices)
        {
            
            std::vector<int> & adjs = VV[index];
            adjs.erase(std::remove(adjs.begin(), adjs.end(), source), adjs.end()); // major solwing part
            
            std::vector<int> & scadjs = VV[source];
            scadjs.erase(std::remove(scadjs.begin(), scadjs.end(), index), scadjs.end());
            // for(auto & adjs: VV)
            // {
            //     adjs.erase(std::remove(adjs.begin(), adjs.end(), index), adjs.end()); // major solwing part
            // }
        }

    }

    void silence_vertices(std::vector<std::vector<int > > & VV, const std::vector<int> & silent_indices)
    {
        for(auto index : silent_indices)
        {
            for(auto vidx: VV[index])
            {
                std::vector<int> & adjs = VV[vidx];
                adjs.erase(std::remove(adjs.begin(), adjs.end(), index), adjs.end()); // major solwing part
            }
            VV[index].clear();
            // for(auto & adjs: VV)
            // {
            //     adjs.erase(std::remove(adjs.begin(), adjs.end(), index), adjs.end()); // major solwing part
            // }
        }
    }





    void setWeights(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXd & V_bad,
        const edge & edg,
        const int & upratio,
        const int & power,
        std::vector<double> & Weights
    )
    {
        Weights.clear();
        Weights.resize(V.rows());
        std::vector<int> Elist = edg._edge_vertices;
        // arc length upsample of
        int target_num = upratio * (Elist.size()-1) + 1;
        Eigen::MatrixXd sample(target_num, 3);
        for(int jj =0 ; jj < target_num ; ++jj)
        {
            int batch = (int)(jj / upratio);
            int remain = jj % upratio;
            if(remain == 0) 
            {
                sample.row(jj) = V_bad.row(Elist[batch]);
            }
            else
            {
                double lambda = (double)(remain) / (double)(upratio);
                sample.row(jj) = (1-lambda)* V_bad.row(Elist[batch]) + lambda* V_bad.row(Elist[batch+1]);
            }
        }
        kd_tree_Eigen<double> sample_kdt(sample.cols(),std::cref(sample),10);
        sample_kdt.index->buildIndex();
        for(int j =0; j < V.rows(); ++j){
            Eigen::RowVector3d jquery= V.row(j);
            int sample_idx= kd_tree_NN_Eigen(sample_kdt, jquery);
            Weights[j] = pow((jquery - sample.row(sample_idx)).norm(), power);
        }
        
        
    }

    void sector_silence_list(
        const Eigen::MatrixXi & F,
        const std::map<int, std::vector<int> > & node_edge_dict,
        const std::map<int, std::map<int, bool> > & node_edge_visit_dict,
        const std::vector<std::vector<int> > & TEdges,
        const std::map<int, std::vector<int> > & CC_node_face_dict, 
        const int& source,
        const int& cur_edge,
        std::vector<int> & temp_silence_list
    )
    {
        temp_silence_list.clear();
        /*
                assume we are dealling with edge e0
                for a source assume  Counter clockwisely there are edges (ei, ebf, ek, e0, ebh) and ebf, and ebh are already set. We need to silence the faces in the sector spanned by ebh - ebf
            */

            // step 1 determin ebf and ebh
            int ebf, ebh; 
            // ebf is the right most edge in CC order to the left of e0
            // ebh is the left most edge in CC order to the right of e0
            // step 1.1 determine the position of e0
            int e0pos =0;
            for(auto ex : node_edge_dict.at(source))
            {
                if(ex != cur_edge)
                {
                    e0pos += 1;
                } 
                else
                {
                    break;
                }   
            }

            // step 1.2 determine the position of ebf
            
            int ebfpos = e0pos;
            bool find_ebf = false;
            int CCsize  = node_edge_dict.at(source).size();
            for(int ii =0; ii < CCsize; ++ii)
            {
                ebfpos =(e0pos - ii + CCsize) % CCsize;
                int ebftemp = node_edge_dict.at(source)[ebfpos];
                if(node_edge_visit_dict.at(source).at(ebftemp))
                {
                    find_ebf = true;
                    break;
                }
            }
            // similarly determine the position of ebh
            int ebhpos = e0pos;
            bool find_ebh = false;
            for(int ii =0; ii < CCsize; ++ii)
            {
                ebhpos =(e0pos + ii) % CCsize;
                int ebhtemp = node_edge_dict.at(source).at(ebhpos);
                if(node_edge_visit_dict.at(source).at(ebhtemp))
                {
                    find_ebh = true;
                    break;
                }
            }
            if( !find_ebf || !find_ebh)
            {
                return; // no other occupied edges
            }
            if(ebfpos == ebhpos) {
                return; // only one occupied edges
            }
            ebf = node_edge_dict.at(source)[ebfpos];
            ebh = node_edge_dict.at(source)[ebhpos];
            
             // currently occupied edge less than one nothing to go further

            // the sector is spanned CC wisely from ebh to ebf
            int sector_start = -1;
            for(auto fidx: CC_node_face_dict.at(source))
            {
                sector_start += 1;
                int e_pos = -1;
                for(auto v_pos: {0, 1, 2})
                {
                    if(F(fidx, v_pos) == source)
                    {
                        e_pos = v_pos;
                        break;
                    }
                }
                // check the out edge of the triangle fidx in TEdges
                // e_pos = v_pos
                if(TEdges[fidx][e_pos]== ebh) // sector start with ebh
                {
                    break;
                }
            }

            // current setor_start store the starting index of faces in CC_node_face_dict[source]
            for(int ff =0 ; ff < CC_node_face_dict.at(source).size() ; ++ff)
            {
                int secface = (sector_start + ff) %CC_node_face_dict.at(source).size() ;
                int e_pos = -1;
                for(auto v_pos: {0, 1, 2})
                {
                    if(F(CC_node_face_dict.at(source)[secface], v_pos) == source)
                    {
                        e_pos = v_pos;
                        break;
                    }
                }
                if(TEdges[CC_node_face_dict.at(source)[secface]][(e_pos+2)%3]== ebf) // sector end with ebf
                {
                    break;
                }
           
                temp_silence_list.push_back(F(CC_node_face_dict.at(source)[secface], (e_pos+2)%3));
                
                
            }
    }

    void update_local_sector(
        const std::vector<std::vector<int> > & VV, 
        const Eigen::MatrixXi & F,
        const std::map<int , std::map<int, bool> > & node_edge_visit_dict,
        const std::map<int, std::vector<int> > & node_edge_dict,
        const std::vector<std::vector<int> > & TEdges,
        const std::map<int, std::vector<int> > & CC_node_face_dict, 
        const int& source,
        const int& target,
        const int& cur_edge,
        std::vector<std::vector<int> > & VV_temp
    )
    {
        VV_temp = VV; // copy VV into VV_temp

        {
            // set silence all other nodes in VV_temp
            std::vector<int> other_node_list;
            for(auto item: node_edge_dict)
            {
                if(item.first != source && item.first != target)
                {
                    other_node_list.push_back(item.first);
                }
            }

            silence_vertices(VV_temp, other_node_list);
        }
        
        // deal with the sector of source first
        std::vector<int> temp_silence_list;
        sector_silence_list(
            F,
            node_edge_dict, 
            node_edge_visit_dict, 
            TEdges, 
            CC_node_face_dict, 
            source, 
            cur_edge, 
            temp_silence_list);
        // if target is in temp_silence_list do not silence it
        temp_silence_list.erase(std::remove(temp_silence_list.begin(), temp_silence_list.end(), target), temp_silence_list.end());
        sector_silence_vertices(VV_temp,source, temp_silence_list);
        // deal with the sctor of target
        sector_silence_list(
            F,
            node_edge_dict, 
            node_edge_visit_dict, 
            TEdges, 
            CC_node_face_dict, 
            target,
            cur_edge, 
            temp_silence_list);
        temp_silence_list.erase(std::remove(temp_silence_list.begin(), temp_silence_list.end(), source), temp_silence_list.end());
        sector_silence_vertices(VV_temp, target, temp_silence_list);
    }


   

    bool trace_for_edge(
        const CellularGraph & cg,
        TraceComplex & tc,
        const int edge_idx,
        std::map<int, std::map<int, bool> > & node_edge_visit_dict,
        json  & path_json,
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
        Trace::Edge_Dijkstra(VV_temp1, target, source, SpWeight,path); // switch the position of source and target to make sure path[0]= source and path[-1] = target
        // dijkstra_trace(VV_temp, source, target, Weights, path);
        if(param.debug)
        {
            Eigen::VectorXd source_target=Eigen::VectorXd::Constant(6,0);
            for(int xx: {0,1,2})
            {
                source_target(xx)=tc._V(source,xx);
                source_target(xx+3)=tc._V(target,xx);
            }
            
            igl::writeDMAT(param.data_root+"/source_target_"+param.tracing+".dmat",source_target);
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
            file.open(param.data_root+"/debug_paths_"+param.tracing+".json");
            file << path_json;
            igl::writeOBJ(param.data_root+"/debug_mesh_"+param.tracing+".obj", tc._V, tc._F);
            igl::writeDMAT(param.data_root+"/FL_"+param.tracing+".dmat", tc._FL);
        }
        // update visit_dict or loop condition update
        node_edge_visit_dict[target_bad][edge_idx]=true;
        node_edge_visit_dict[source_bad][edge_idx]=true;
        return true;
    }

    bool MatchMakerTree(
        const CellularGraph & cg,
        Eigen::MatrixXd & V_good,
        Eigen::MatrixXi & F_good,
        Eigen::VectorXi & FL_good,
        const params param
    )
    {
        int total_label_num = cg.label_num;
        
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
        std::vector<std::pair<int, std::pair<int, int> > > frame_graph, frame_MST;
        _build_frame_graph(cg._edge_list,frame_graph);
        frame_MST = Algo::Kruskal_MST(frame_graph);
        json mst_json; // store the min spanning tree
        for(auto item: frame_MST)
        {
            mst_json[std::to_string(item.first)] = cg._edge_list[item.first]._edge_vertices;   
        }

        TraceComplex tc;
        tc.initialize(V_good, F_good);
        for(auto edg: cg._edge_list)
        {
            int source_bad = edg.head;
            int target_bad = edg.tail;
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

            int source = -1;
                int target  = -1;

                
                
                if(tc._node_image_map.find(source_bad)==tc._node_image_map.end())
                {
                    proj_node_loop(cg, source_bad, tc, source); 
                    // updated the tc._node_listd and tc._node_image_map;
                    assert(source!= -1);               
                }
                if(tc._node_image_map.find(target_bad)==tc._node_image_map.end())
                {
        
                    proj_node_loop(cg, target_bad, tc, target);
                    assert(target!= -1);
                }

        }

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

        // trace for mst
        for(auto frame_edge: frame_MST)
        {
            int edge_idx = frame_edge.first;

            if(!trace_for_edge(
                cg,
                tc,
                edge_idx,
                node_edge_visit_dict,
                path_json,
                param 
            )) return false;
        }
        std::cout << "tracing mst done" << std::endl;
        // trace of remaining edges

        for(auto item: cg._node_edge_dict)
        {
            // Main loop for tracing
            int nd = item.first;
            for (auto edge_idx: cg._node_edge_dict.at(nd)){
                if(node_edge_visit_dict[nd][edge_idx])
                {
                    continue;
                }
                if(!trace_for_edge(
                    cg,
                    tc,
                    edge_idx,
                    node_edge_visit_dict,
                    path_json,
                    param 
                )) return false;
            }
        }

        for(auto item: cg._patch_edge_dict)
        {
            int patch_idx = item.first;

            int ff_in = _locate_seed_face(cg, tc, patch_idx);
            loop_colorize(tc._V, tc._F, tc._TEdges,ff_in, patch_idx, tc._FL);
        }
        if(param.debug)
        {
            igl::writeOBJ(param.data_root+"/debug_mesh_tree.obj", tc._V, tc._F);
            igl::writeDMAT(param.data_root+"/FL_tree.dmat", tc._FL);
            CellularGraph cgt;
            Bijection::TransferCellGraph(cg, tc, cgt);
            Eigen::MatrixXd M_s2t;
            Bijection::BijGlobal(cg,cgt, M_s2t);
            Eigen::MatrixXd Vmap=igl::barycentric_to_global(cgt.V, cgt.F, M_s2t);
            igl::writeOBJ(param.data_root+"/map.obj", Vmap, cg.F);
        }
        F_good = tc._F;
        V_good = tc._V;
        FL_good = tc._FL;
        return true;

    }





}
}
