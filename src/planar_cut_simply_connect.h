#include <Eigen/Core>
#include <vector>
#include <igl/bfs.h>
#include <igl/adjacency_list.h>
#include <Eigen/Sparse>
#include <algorithm>
namespace bcclean{
    void silence_vertices(std::vector<std::vector<int > > & VV, const std::vector<int> & silent_indices){
        for(auto index : silent_indices){
            VV[index].clear();
            for(auto & adjs: VV){
                std::remove(adjs.begin(), adjs.end(), index);
            }
        }
    }

    void planar_cut_simply_connect(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, const std::vector<std::vector<int> > & boundary_loops, std::vector<bool> & VCuts){
        // V #V x 3 vertices of a nonsimply connected mesh
        // F #F x 3 faces of the mesh
        //  std::vector<int> & VLoopLabels, #V x 1 indicate whether a vertice is on a loop if not label = -1 if on label = loop index ( loop index starting from 0)
        // return Vcuts # V x 1 indicate a vertice is on the cut or boundary.


        // resize VCuts
        std::vector<int> VLoops(V.rows());
        assert(boundary_loops.size()>1);
        VCuts.resize(V.rows());
        for(int count =0; count <V.rows(); ++ count){
            VCuts[count]= false;
        }



        // initialize VLoops and VCuts;
        int loop_idx = 0;
        for(auto loop: boundary_loops){
            for(auto vidx: loop){
                VCuts[vidx]=true;
                VLoops[vidx]=loop_idx;
            }
            loop_idx += 1;
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

            // find next root;
            int mid = boundary_loops[curloop].size() % 2;
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
                std::remove(remaining_loops.begin(), remaining_loops.end(),curloop);
            }
        }
    }
}