#include "non_vertex_manifold_relabel.h"
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <queue>
#include <list>
#include <map>
namespace bcclean {
namespace Prepocess{
    void non_vertex_manifold_relabel(
        const Eigen::MatrixXd & Vraw,
        const Eigen::MatrixXi & Fraw, 
        const Eigen::VectorXi & FI,
        const std::vector<int> & NMV,
        const Eigen::VectorXi & FL, 
        const int cur_lb,
        Eigen::VectorXi & FL_mod, 
        int & total_label_num,
        std::map<int, Eigen::VectorXi> & subpatches
    )
    {
        subpatches.clear();
        std::map<int, int> subpatch_face_count;
        Eigen::MatrixXi TT;
        Eigen::VectorXi FL_patch= Eigen::VectorXi::Constant(Fraw.rows(), -1);
        igl::triangle_triangle_adjacency(Fraw, TT);
        std::vector<std::vector<int> > VF, VFi;
        igl::vertex_triangle_adjacency(Vraw, Fraw, VF, VFi);
        std::vector<int> branch_starts;
        for( auto probv : NMV)
        {
            std::list<int> Tpool; 
            for(auto elem : VF.at(probv))
            {
                Tpool.push_back(elem);
            }
            while(Tpool.size()>0)
            {
                int root  = Tpool.front();
                Tpool.remove(root);
                branch_starts.push_back(root);
                std::queue<int> search_queue;
                search_queue.push(root);
                while(!search_queue.empty())
                {
                    int cur = search_queue.front();
                    search_queue.pop();
                    std::list<int> remove_list;
                    for(auto cand: Tpool)
                    {
                        if((TT.row(cand).array()==cur).any())
                        {
                            remove_list.push_back(cand);
                        }
                    }
                    for(auto elem: remove_list)
                    {
                        Tpool.remove(elem);
                        search_queue.push(elem);
                    }
                }
            }
        }
        std::vector<std::queue<int> > SQs;
        std::vector<int> nLB;
        int count = 0;
        for(auto ss: branch_starts)
        {
            std::queue<int> sq;
            sq.push(ss);
            SQs.push_back(sq);
            if(! nLB.empty()){
                nLB.push_back(count + total_label_num - 1);
                FL_patch(ss) =count + total_label_num - 1;
            }
            else{
                nLB.push_back(cur_lb);
                FL_patch(ss) = cur_lb;
            }
            count +=1;
        }
        total_label_num = total_label_num + count-1;

        bool finished= false;
        while( ! finished)
        {
            finished = true;
            for(int j = 0; j < SQs.size(); ++j)
            {
                std::queue<int> & sq = SQs[j];
                int lb = nLB[j];
                if(! sq.empty())
                {
                    finished = false;
                    int cur = sq.front();
                    sq.pop();
                    for(auto j : {0,1,2})
                    {
                        int adj = TT(cur, j);
                        if(adj != -1 && FL_patch(adj)== -1)
                        {
                            std::vector<int> neck;
                            for(auto k : {0,1,2})
                            {
                                int test = TT(adj,k);
                                if(test != -1 && test != cur)
                                {
                                    neck.push_back(test);
                                }
                            }
                            if(neck.size()==2)
                            {
                                int lb0 = FL_patch(neck[0]);
                                int lb1 = FL_patch(neck[1]);
                                if(lb1 == lb0)
                                {
                                    if(lb1 != lb && lb1 != -1)
                                    {
                                        continue;
                                    }
                                }
                            }
                            sq.push(adj);
                            FL_patch(adj) = lb;
                        }
                    }
                }
            }
            // check finish;
        }
        for(auto lb: nLB)
        {
            subpatches[lb] = Eigen::VectorXi::Constant(Fraw.rows(), 0);
            subpatch_face_count[lb]= 0;
        }
        for(int fidx = 0; fidx < Fraw.rows(); ++ fidx)
        {
            FL_mod(FI(fidx)) = FL_patch(fidx);
            int lb = FL_patch(fidx);
            subpatches[lb](subpatch_face_count[lb]) = FI(fidx);
            subpatch_face_count[lb]+=1;
        }
        for(auto lb: nLB)

        {
            subpatches[lb].conservativeResize(subpatch_face_count[lb]);
        }

        return;
    }
}
}