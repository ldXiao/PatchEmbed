#include "edge.h"
#include <Eigen/Core>
#include <unordered_map>
#include <iostream>
#include <igl/remove_unreferenced.h>
#include <igl/boundary_loop.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/list_to_matrix.h>
#include <algorithm>
#include <utility>
namespace bcclean{
    bool _match_two_loops(std::vector<int> loop_a, std::vector<int> loop_b){
        if(loop_a.size()!= loop_b.size()){return false;}
        int n = loop_a.size();
        int head_a = loop_a[0];
        int start_b = 0;
        bool found = false;
        bool reversed = true;
        bool match_vertices = true;
        for(auto b_num: loop_b){
            if(b_num == head_a){
                found = true;
                break;
            }
            start_b+=1;
        }
        if(!found){return false;}
        if(n==1){return true;}
        {
            int neck_a = loop_a[1];
            int neck_b1 = loop_b[(start_b+1)% n];
            int neck_b2 = loop_b[(start_b-1+n)%n];
            if(neck_a == neck_b1){
                reversed = false;
            } else if(neck_a == neck_b2){
                reversed = true;
            } else {
                return false;
            }
        }
        
        for(int offset =0; offset <n ; ++offset){
            int na = loop_a[offset];
            int nb;
            {
            if(!reversed){nb = loop_b[(start_b+offset)%n];}
            else{nb = loop_b[(start_b-offset+n)%n];}
            }
            if(na !=nb){
                match_vertices= false;
                break;
            }
        }
        return match_vertices;

    }

    Edge_Compare_Result edge::_compare_edge(edge b){
        if(b._label_pair == _label_pair){
            if(b.head == -1 && b.tail == -1){
                // b is a loop
                if(head == -1 && head == -1){
                    // this edge is also a loop
                    // continue to compare
                    if(b._edge_vertices.size()!= _edge_vertices.size()){
                        return Edge_Compare_Result::SAME_TYPE;
                    } else{
                        if(_match_two_loops(_edge_vertices, b._edge_vertices)){return Edge_Compare_Result::IDENTICAL;}
                        else return Edge_Compare_Result::SAME_TYPE;
                    }
                } else return Edge_Compare_Result::NON_MATCH;
            }
            bool match_node = (b.head == head) && (b.tail == tail);
            bool anti_match_node = (b.head == tail) && (b.tail == head);
            if(match_node || anti_match_node){
                bool match_vertices = true;
                int start = (match_node) ? 0 : _edge_vertices.size()-1;
                int increment = (match_node) ? 1 : -1;
                int idx = start;
                if(b._edge_vertices.size() != _edge_vertices.size()){
                    return Edge_Compare_Result::SAME_TYPE;
                }
                for(auto p: b._edge_vertices){
                    if(p == _edge_vertices[idx]){
                        idx += increment;
                    } else{
                        match_vertices = false;
                        break;
                    }
                }
                if(match_vertices){
                    return Edge_Compare_Result::IDENTICAL;
                } else {
                    return Edge_Compare_Result::SAME_TYPE;
                }
                
            } else {
                return Edge_Compare_Result::SAME_TYPE;
            }
        } else{
            return Edge_Compare_Result::NON_MATCH;
        }
    }

    bool extract_patch_mesh(
            const Eigen::MatrixXd& V, 
            const Eigen::MatrixXi& F, 
            const Eigen::MatrixXi&FL, 
            const int lb_i, 
            Eigen::MatrixXd& V_i, 
            Eigen::MatrixXi& F_i,
            Eigen::VectorXi& I_i){
                Eigen::MatrixXi F_l = Eigen::MatrixXi::Constant(F.rows(),3, 0);
                int count = 0;
                for(int fidx=0; fidx< FL.rows(); ++fidx){
                    int lb = FL(fidx,0);
                    if(lb==lb_i){
                        F_l.row(count)=F.row(fidx);
                        count += 1;
                    }
                }
                if(count == 0){return false;}
                F_l.conservativeResize(count,3);
                {   
                    Eigen::VectorXi J_i;
                    igl::remove_unreferenced(V, F_l, V_i, F_i,J_i, I_i);
                }
                return true;
        }
    
    bool build_vertex_label_list_dict(const Eigen::MatrixXi&F, const Eigen::MatrixXi & FL, const size_t total_label_num, std::unordered_map<int, std::vector<int>> & vertex_label_list_dict){
        // vertex_label_list_dict store all the vertex_indx- label-list(ordered)
        vertex_label_list_dict.clear();
        
        for(int fidx=0; fidx< F.rows(); ++fidx){
            int lb = FL(fidx,0);
            for(int i=0; i <3; ++i){
                int vi = F(fidx, i);
                auto jt = vertex_label_list_dict.find(vi);
                if(jt == vertex_label_list_dict.end()){
                    // std::cout << "pushed"<<std::endl;
                    vertex_label_list_dict[vi]=std::vector<int>();
                    // vertex_label_list_dict[vi].clear();
                }
                auto & vect = vertex_label_list_dict[vi];
                auto it = std::find (vect.begin(), vect.end(), lb);
                if(it == vect.end()){
                    vect.push_back(lb);
                }
            }
        }
        // for(auto & vect: vertex_label_list_dict){
        //     std::sort(vect.second.begin(),vect.second.end());
        //     for(auto q: vect.second){
        //         std::cout << q <<" ";
        //     }
        //     std::cout  << "---------"<<std::endl;
        // }
        return true;
    }

    std::pair<int, int> _adjacent2vertex2labels(const Eigen::MatrixXi & FL, const std::vector<std::vector<int>> & VF, int va, int vb){
        // VF should have been initliazed by igl::vertex_triangle_adjacency(V, F, VF,VFi)
        // va, vb are two vertices
        // find the two triangle that share the two vertices
        std::vector<int> fas = VF.at(va);
        std::vector<int> fbs = VF.at(vb);
        std::vector<int> inter(10);
        std::vector<int>::iterator it;
        it=std::set_intersection (fas.begin(),fas.end(), fbs.begin(), fbs.end(), inter.begin());
        inter.resize(it-inter.begin());
        try{
            if(inter.size()==2){
                inter[0] = FL(inter[0],0);
                inter[1] = FL(inter[1],0);
                std::sort(inter.begin(), inter.end());
                return std::pair<int, int>(inter[0], inter[1]);
            } else {
                throw;
            }
        }
        catch(...){std::cout << "Exception occurred"; exit(EXIT_FAILURE);}
    }

    bool build_edge_list(const Eigen::MatrixXd& V, const Eigen::MatrixXi&F, const Eigen::MatrixXi & FL, const size_t total_label_num, std::vector<edge> & edge_list, std::unordered_map<int, std::vector<int> > & patch_edge_dict){
        // edge list store all the edge uniquely in a vector
        // patch_edge_dict store the patch label -> list of edges in loop order
        // return lb->vector<nodes> where the nodes are unodered for each patch
        std::unordered_map<int, std::vector<int>> vertex_label_list_dict;
        edge_list.clear();
        patch_edge_dict.clear();
        std::vector<std::vector<int>> VF; // vertex - face index list
        {
            std::vector<std::vector<int>> VFi;
            igl::vertex_triangle_adjacency(V, F, VF,VFi);
        }
        // reinitialize
        // vertx_index -> label_list
        if(! build_vertex_label_list_dict(F, FL, total_label_num, vertex_label_list_dict)){
            return false;
        }
        for(int lb_i=0; lb_i < total_label_num; ++lb_i){

            patch_edge_dict[lb_i] = std::vector<int>();
            // loop for each patch
            Eigen::MatrixXd V_i;
            Eigen::MatrixXi F_i;
            Eigen::VectorXi I_i;
            extract_patch_mesh(V, F, FL, lb_i, V_i, F_i, I_i);
            std::vector<std::vector<int> > L;
            igl::boundary_loop(F_i, L);
            if(L.size()!=1){
                // input mesh is not isomorphic to disk
                std::cout << "not simply connected" << std::endl;
                std::cout << L.size()<< " loops detected, treated separately" << std::endl;
            }
            for(auto loop : L){
                // for each loop detect the first node to set as a starting head of first edge
                int start = -1;
                int count = 0;
                for(auto sss: loop){
                                    std::cout << I_i(sss) << ' ';
                                }
                
                std::cout << "-----loop"<<std::endl;
                if(loop.size()>2){
                std::vector<edge> patch_local_edges;
                for(auto v_idx:loop){
                    int v_idx_raw = I_i(v_idx);
                    if(vertex_label_list_dict[v_idx_raw].size() > 2){
                        //detected a node and set it to be start;
                        start = count;
                        break;
                    }
                    count +=1;
                }
                // after the above loop start is either -1 or the index on loop
                if(start == -1){
                    // this loop has no node
                    edge edg;
                    edg.head = -1;
                    edg.tail = -1;
                    edg.total_label_num = total_label_num;
                    std::vector<int> edge_vertices;
                    for(auto p: loop){edge_vertices.push_back(I_i(p));}
                    // p is an index on V_i use I_i to source back to V
                    edg._edge_vertices = edge_vertices;
                    int va  = edge_vertices[0];
                    int vb  = edge_vertices[1];
                    edg._label_pair = _adjacent2vertex2labels(FL, VF, va, vb);
                    patch_local_edges.push_back(edg);
                } else {
                    // start node detect, loop over all edges to push into pathc_local_edges
                    std::vector<int> curr_path = {I_i(loop[start])};
                    for(int offset =1; offset < loop.size()+1;++offset){
                        int loop_idx  = (start + offset) % loop.size();
                        int v_idx = loop[loop_idx];
                        int v_idx_raw = I_i(v_idx);
                        if(vertex_label_list_dict[v_idx_raw].size() > 2){
                            //detected a node and set it to be tail curr_path;
                            curr_path.push_back(v_idx_raw);
                            // a full path is detected , constuct the corresponding edge;
                            {
                                edge edg;
                                edg.head =curr_path[0];  
                                edg.tail = v_idx_raw;
                                edg.total_label_num = total_label_num;
                                edg._edge_vertices = curr_path;
                                int va  = edg._edge_vertices[0];
                                int vb  = edg._edge_vertices[1];
                                for(auto sss: curr_path){
                                    std::cout << sss << ' ';
                                }
                                std::cout << "-----"<<std::endl;
                                edg._label_pair = _adjacent2vertex2labels(FL, VF, va, vb);
                                patch_local_edges.push_back(edg);
                            }
                            curr_path.clear();
                        }
                        curr_path.push_back(v_idx_raw);
                    }
                }
                
                std::vector<bool> inserts;
                std::vector<int> local_indices_list; // store the final edge indices in this patch
                int count_local = 0;
                for(edge & edg: patch_local_edges){
                    bool insert = true;
                    int count_glob =0;
                    int final_idx = -1;
                    for(edge & edg1: edge_list){
                        if(edg._compare_edge(edg1)==Edge_Compare_Result::IDENTICAL){
                            insert = false;
                            final_idx = count_glob;
                            break;
                        }
                        count_glob +=1;
                    }
                    if(final_idx==-1){
                        final_idx = count_local+ edge_list.size();
                    }
                    inserts.push_back(insert);
                    local_indices_list.push_back(final_idx);
                    if(insert){count_local+=1;}
                }
                for(int i=0; i < inserts.size();++i){
                    if(inserts[i]){
                        edge_list.push_back(patch_local_edges[i]);
                    }
                    patch_edge_dict[lb_i].push_back(local_indices_list[i]);
                }
            }
            }
        }
        return true;
    }

    bool build_edge_list_loop(const Eigen::MatrixXd& V, const Eigen::MatrixXi&F, const Eigen::MatrixXi & FL, const size_t total_label_num, std::vector<edge> & edge_list, std::unordered_map<int, std::vector<int> > & patch_edge_dict,
    std::unordered_map<int, std::vector<bool> > & patch_edgedirec_dict, int & largest_patch){
        // edge list store all the edge uniquely in a vector
        // patch_edge_dict store the patch label -> list of edges in loop order
        // return lb->vector<nodes> where the nodes are unodered for each patch
        std::unordered_map<int, std::vector<int>> vertex_label_list_dict;
        edge_list.clear();
        patch_edge_dict.clear();
        patch_edgedirec_dict.clear();
        std::vector<std::vector<int>> VF; // vertex - face index list
        {
            std::vector<std::vector<int>> VFi;
            igl::vertex_triangle_adjacency(V, F, VF,VFi);
        }
        // reinitialize
        // vertx_index -> label_list
        if(! build_vertex_label_list_dict(F, FL, total_label_num, vertex_label_list_dict)){
            return false;
        }
        int largest_patch_idx = 0;
        int largest_patch_traingle_count=-1;
        for(int lb_i=0; lb_i < total_label_num; ++lb_i){
            
            patch_edge_dict[lb_i] = std::vector<int>();
            // loop for each patch
            Eigen::MatrixXd V_i;
            Eigen::MatrixXi F_i;
            Eigen::VectorXi I_i;
            extract_patch_mesh(V, F, FL, lb_i, V_i, F_i, I_i);
            if(F_i.rows()> largest_patch_traingle_count)
            {
                largest_patch_idx = lb_i;
                largest_patch_traingle_count = F_i.rows();
            }
            std::vector<std::vector<int> > L;
            igl::boundary_loop(F_i, L);
            if(L.size()!=1){
                // input mesh is not isomorphic to disk
                std::cout << "not simply connected" << std::endl;
                std::cout << L.size()<< " loops detected, treated separately" << std::endl;
            }
            for(auto loop : L){
                // for each loop detect the first node to set as a starting head of first edge
                int start = -1;
                int count = 0;
                for(auto sss: loop){
                                    std::cout << I_i(sss) << ' ';
                                }
                
                std::cout << "-----loop"<<std::endl;
                if(loop.size()>2){
                std::vector<edge> patch_local_edges;
                for(auto v_idx:loop){
                    int v_idx_raw = I_i(v_idx);
                    if(vertex_label_list_dict[v_idx_raw].size() > 2){
                        //detected a node and set it to be start;
                        start = count;
                        break;
                    }
                    count +=1;
                }
                // after the above loop start is either -1 or the index on loop
                if(start == -1){
                    // this loop has no node
                    edge edg;
                    edg.head = -1;
                    edg.tail = -1;
                    edg.total_label_num = total_label_num;
                    std::vector<int> edge_vertices;
                    for(auto p: loop){edge_vertices.push_back(I_i(p));}
                    // p is an index on V_i use I_i to source back to V
                    edg._edge_vertices = edge_vertices;
                    int va  = edge_vertices[0];
                    int vb  = edge_vertices[1];
                    edg._label_pair = _adjacent2vertex2labels(FL, VF, va, vb);
                    patch_local_edges.push_back(edg);
                } else {
                    // start node detect, loop over all edges to push into pathc_local_edges
                    std::vector<int> curr_path = {I_i(loop[start])};
                    for(int offset =1; offset < loop.size()+1;++offset){
                        int loop_idx  = (start + offset) % loop.size();
                        int v_idx = loop[loop_idx];
                        int v_idx_raw = I_i(v_idx);
                        if(vertex_label_list_dict[v_idx_raw].size() > 2){
                            //detected a node and set it to be tail curr_path;
                            curr_path.push_back(v_idx_raw);
                            // a full path is detected , constuct the corresponding edge;
                            {
                                edge edg;
                                edg.head =curr_path[0];  
                                edg.tail = v_idx_raw;
                                edg.total_label_num = total_label_num;
                                edg._edge_vertices = curr_path;
                                int va  = edg._edge_vertices[0];
                                int vb  = edg._edge_vertices[1];
                                for(auto sss: curr_path){
                                    std::cout << sss << ' ';
                                }
                                std::cout << "-----"<<std::endl;
                                edg._label_pair = _adjacent2vertex2labels(FL, VF, va, vb);
                                patch_local_edges.push_back(edg);
                            }
                            curr_path.clear();
                        }
                        curr_path.push_back(v_idx_raw);
                    }
                }
                
                std::vector<bool> inserts;
                std::vector<bool> local_direc_list; // store whether the direction of final edge is correct in this path
                std::vector<int> local_indices_list; // store the final edge indices in this patch
                int count_local = 0;
                for(edge & edg: patch_local_edges){
                    bool insert = true;
                    bool correc_direc = true;
                    int count_glob =0;
                    int final_idx = -1;
                    for(edge & edg1: edge_list){
                        if(edg._compare_edge(edg1)==Edge_Compare_Result::IDENTICAL){
                            insert = false;
                            if(edg.head != edg1.head)
                            {
                                correc_direc = false;
                            }
                            final_idx = count_glob;
                            break;
                        }
                        count_glob +=1;
                    }
                    if(final_idx==-1){
                        final_idx = count_local+ edge_list.size();
                    }
                    inserts.push_back(insert);
                    local_direc_list.push_back(correc_direc);
                    local_indices_list.push_back(final_idx);
                    if(insert){count_local+=1;}
                }
                for(int i=0; i < inserts.size();++i){
                    if(inserts[i]){
                        edge_list.push_back(patch_local_edges[i]);
                    }
                    patch_edge_dict[lb_i].push_back(local_indices_list[i]);
                    patch_edgedirec_dict[lb_i].push_back(local_direc_list[i]);
                }
            }
            }
        }
        largest_patch = largest_patch_idx;
        return true;
    }



    
}