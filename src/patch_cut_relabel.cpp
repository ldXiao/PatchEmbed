#include "patch_cut_relabel.h"
#include <igl/triangle_triangle_adjacency.h>
#include <queue>
namespace bcclean{
    void patch_cut_relabel(const Eigen::MatrixXi & Fraw, const Eigen::VectorXi FI, const std::vector<bool> & VCuts, const Eigen::VectorXi & FL, Eigen::VectorXi & FL_mod, int & total_label_num){
        assert(total_label_num == FL.maxCoeff()+1);
        FL_mod = FL;
        int ini_label = FL(FI(0));
        Eigen::VectorXi PatchL = Eigen::VectorXi::Zero(Fraw.rows());
        for(int frawidx =0 ; frawidx < Fraw.rows(); ++frawidx){
            PatchL(frawidx) = -1;
        }
        Eigen::MatrixXi TT, TTi;
        igl::triangle_triangle_adjacency(Fraw, TT, TTi);
        // start with face 0 use TT and TTi info to get connected components
        // traverse all faces connected to face 0 and do not cross cuts label them to be 0;
        int CC = 0; // number of connected components
        int root_face = 0;
        bool find_empty = true;
        std::queue<int> search_queue;
        search_queue.push(root_face);
        while(find_empty){
            while(search_queue.size()!=0){
                int cur_face = search_queue.front();
                search_queue.pop(); // remove head
                
                PatchL(cur_face) = CC;
                Eigen::RowVector3i adjs = TT.row(cur_face);
                for(int j =0 ; j <3 ; ++j){
                    int face_j = adjs(j);
                    int vidx, vidy;
                    vidx = Fraw(cur_face, j);
                    vidy  = Fraw(cur_face, (j+1)%3);
                    if(!(VCuts[vidx] && VCuts[vidy])){
                        // the edge is not in VCuts
                        if(PatchL(face_j)== -1){
                            // not visited before;
                            search_queue.push(face_j);
                        }
                    }
                }
            }
            find_empty = false;
            for(int i =0; i < PatchL.rows();++i){
                if(PatchL(i)== -1){
                    root_face = i;
                    find_empty = true;
                    search_queue.push(root_face);
                    CC += 1;
                    break; // found an empty face to label
                }
            }
        }
        for(int frawidx =0 ; frawidx < Fraw.rows(); ++frawidx){
            if(PatchL(frawidx)==0){
                FL_mod(FI(frawidx))= ini_label;
            } else {
                FL_mod(FI(frawidx)) = total_label_num + PatchL(frawidx); 
            }
        }
        total_label_num+= CC;
    }
}