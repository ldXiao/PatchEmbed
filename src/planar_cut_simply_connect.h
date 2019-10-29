#include <Eigen/Core>
#include <vector>
#include <igl/bfs.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <Eigen/Sparse>
#include <algorithm>
namespace bcclean{
    void silence_vertices(std::vector<std::vector<int > > & VV, const std::vector<int> & silent_indices){
        for(auto index : silent_indices){
            VV[index].clear();
            for(auto & adjs: VV){
                adjs.erase(std::remove(adjs.begin(), adjs.end(), index), adjs.end());
            }
        }
    }

    void planar_cut_simply_connect(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, const std::vector<std::vector<int> > & boundary_loops, std::vector<bool> & VCuts, std::vector<std::vector<bool> > & TCuts){
        // V #V x 3 vertices of a nonsimply connected mesh
        // F #F x 3 faces of the mesh
        //  std::vector<int> & VLoopLabels, #V x 1 indicate whether a vertice is on a loop if not label = -1 if on label = loop index ( loop index starting from 0)
        // return 
        //Vcuts # V x 1 indicate a vertice is on the cut or boundary.
        // TCuts # #F * 3 indicate an edge in triangle is in cut or boundary


        // resize VCuts
        std::vector<int> VLoops(V.rows());
        assert(boundary_loops.size()>1);
        VCuts.resize(V.rows());
        TCuts.resize(F.rows());
        for(int count =0; count <V.rows(); ++ count){
            VCuts[count]= false;
        }
        for(int fcount = 0; fcount  < F.rows(); ++fcount){
            TCuts[fcount] = {false, false, false};
        }

        std::vector<std::vector<int > > VF;// #V list of lists of incident faces (adjacency list)
        std::vector<std::vector<int > > VFi;  //  #V list of lists of index of incidence within incident faces listed in VF
        igl::vertex_triangle_adjacency(V, F, VF, VFi);
        for(auto loop: boundary_loops){
            for(int loop_idx=0; loop_idx < loop.size(); ++loop_idx){
                int uidx = loop[loop_idx];
                int vidx = loop[(loop_idx+1)%loop.size()];
                std::vector<int> inter;
                std::set_intersection(VF[uidx].begin(), VF[uidx].end(), VF[vidx].begin(), VF[vidx].end(), std::back_inserter(inter));
                // there should be only one comman adjacent triangle for boundary vertices
                for(auto trg: inter){
                    for(int edgpos =0; edgpos < 3 ; ++edgpos){
                        int uuidx = F(trg, edgpos);
                        int vvidx = F(trg, (edgpos+1)%3);
                        if(uuidx == uidx && vvidx == vidx){
                            TCuts[trg][edgpos] = true;
                            break;
                        }
                        if(uuidx == vidx && vvidx == uidx){
                            TCuts[trg][edgpos] = true;
                            break;
                        }
                    }
                }
            }
        }
        // initialize VLoops and VCuts;
        int loop_num = 0;
        for(auto loop: boundary_loops){
            for(auto vidx: loop){
                VCuts[vidx]=true;
                VLoops[vidx]=loop_num;
            }
            loop_num += 1;
        }

        std::vector<std::vector<int> > VV;
        igl::adjacency_list(F, VV); // initialize the adjacency list;

        // start from loop 0
        // start with the first by default
        std::map<int, int> loop_visit_count;
        std::vector<int> remaining_loops;
        for(int i =0; i < boundary_loops.size();++i){
            remaining_loops.push_back(i);
            loop_visit_count[i]=0;
        }
        int curloop = 0;
        int root=boundary_loops[curloop][0];
        loop_visit_count[curloop]=1;
        while(remaining_loops.size() != 0){
            std::vector<int > D, P;
            igl::bfs(VV, root, D, P); // D is the spanning tree in the discover order and P store the predessors, all contain indices into VV
            int tail;
            for(auto vidx: D){
                if(VCuts[vidx]){
                    if(VLoops[vidx] ==0 && remaining_loops.size()>1){
                        continue; 
                    }
                    if(VLoops[vidx]== -1){continue;}
                    if(loop_visit_count[VLoops[vidx]] == 2){
                        continue;
                    }
                    tail = vidx;
                    break;
                } else{
                    continue;
                }
            }
            curloop = VLoops[tail];
            loop_visit_count[curloop] += 1;
            int cur  = tail;
            std::vector<int> records;
            while(cur != -1){
                // root not reach, go on back search
                records.push_back(cur);
                VCuts[cur] = true;
                cur = P[cur];
            }
            // remove the indices in records in the adjacency info
            silence_vertices(VV, records);

            // set the triangle edges in cut to be true
            for(int rc_idx=0; rc_idx < records.size()-1; ++rc_idx){
                int uidx = records[rc_idx];
                int vidx = records[(rc_idx+1)%records.size()];
                std::vector<int> inter;
                std::set_intersection(VF[uidx].begin(), VF[uidx].end(), VF[vidx].begin(), VF[vidx].end(), std::back_inserter(inter));
                // there should be only one comman adjacent triangle for boundary vertices
                for(auto trg: inter){
                    for(int edgpos =0; edgpos < 3 ; ++edgpos){
                        int uuidx = F(trg, edgpos);
                        int vvidx = F(trg, (edgpos+1)% 3);
                        if(uuidx == uidx && vvidx == vidx){
                            TCuts[trg][edgpos] = true;
                            break;
                        }
                        if(uuidx == vidx && vvidx == uidx){
                            TCuts[trg][edgpos] = true;
                            break;
                        }
                    }
                }
            }
            // find next root;
            int mid = int(boundary_loops[curloop].size() / 2);
            int tail_idx=-1;
            for(auto vidx :boundary_loops[curloop]){
                tail_idx+=1;
                if(vidx == tail){
                    break;
                }
            }
            root = boundary_loops[curloop][(tail_idx+mid) %  boundary_loops[curloop].size()];
            // visited the loop twice, pop it out from remaining_loops
            if(curloop != 0){loop_visit_count[curloop] +=1;}
            if(loop_visit_count[curloop]>=2){
                remaining_loops.erase(std::remove(remaining_loops.begin(), remaining_loops.end(),curloop), remaining_loops.end());
            }
        }
    }
}