#include "patch_cut_relabel.h"
#include <igl/triangle_triangle_adjacency.h>
#include <queue>
namespace bcclean{
    void patch_cut_relabel(const Eigen::MatrixXi & Fraw, const Eigen::VectorXi & FI, const std::vector<bool> & VCuts, const std::vector<std::vector<bool> > & TCuts, const Eigen::VectorXi & FL, Eigen::VectorXi & FL_mod, int & total_label_num)
    {
        assert(total_label_num == FL.maxCoeff()+1);
        FL_mod = FL;
        int ini_label = FL(FI(0));
        Eigen::VectorXi PatchL = Eigen::VectorXi::Constant(Fraw.rows(), -1);
        // for(int frawidx =0 ; frawidx < Fraw.rows(); ++frawidx){
        //     PatchL(frawidx) = -1;
        // }
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
                    if(face_j == -1) {
                        continue;
                    }
                    // int vidx, vidy;
                    // vidx = Fraw(cur_face, j);
                    // vidy  = Fraw(cur_face, (j+1)%3);
                    if(!TCuts[cur_face][j]){
                        // the edge is not in VCuts
                        if(PatchL(face_j)== -1){
                            // not visited before;
                            PatchL(face_j) = CC;
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
                int newlb = total_label_num + PatchL(frawidx) -1;
                FL_mod(FI(frawidx)) = newlb; 
            }
        }
        
        total_label_num+= CC;
        std::cout << "patch " << ini_label << "cut into " << CC+1 << "pieces" << std::endl;
        // std::cout << PatchL << std::endl;
    }

    void patch_cut_relabel(
        const Eigen::MatrixXi & Fraw, 
        const Eigen::VectorXi & FI, 
        const std::vector<bool> & VCuts, 
        const std::vector<std::vector<bool> > & TCuts, 
        const Eigen::VectorXi & FL, 
        Eigen::VectorXi & FL_mod, 
        std::map<int, Eigen::VectorXi> & subpatchFI_dict)
    {
        int total_label_num = FL.maxCoeff()+1;
        FL_mod = FL;
        int ini_label = FL(FI(0));
        Eigen::VectorXi PatchL = Eigen::VectorXi::Constant(Fraw.rows(), -1);
        // for(int frawidx =0 ; frawidx < Fraw.rows(); ++frawidx){
        //     PatchL(frawidx) = -1;
        // }
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
                    if(face_j == -1) {
                        continue;
                    }
                    // int vidx, vidy;
                    // vidx = Fraw(cur_face, j);
                    // vidy  = Fraw(cur_face, (j+1)%3);
                    if(!TCuts[cur_face][j]){
                        // the edge is not in VCuts
                        if(PatchL(face_j)== -1){
                            // not visited before;
                            PatchL(face_j) = CC;
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
        subpatchFI_dict.clear();
        std::map<int,int> subpatchFcount_dict;
        for(int i = 0; i < CC+1; ++i)
        {
            if(i == 0)
            {
                subpatchFI_dict[ini_label] = Eigen::VectorXi::Constant(Fraw.rows(),-1);
                subpatchFcount_dict[ini_label]=0;
            }
            else
            {
                subpatchFI_dict[total_label_num+i-1] = Eigen::VectorXi::Constant(Fraw.rows(),-1);
                subpatchFcount_dict[total_label_num+i-1] = 0;
            }
        }
        for(int frawidx =0 ; frawidx < Fraw.rows(); ++frawidx){
            if(PatchL(frawidx)==0){
                FL_mod(FI(frawidx))= ini_label;
                subpatchFI_dict[ini_label](subpatchFcount_dict[ini_label]) = FI(frawidx);
                subpatchFcount_dict[ini_label] +=1;
            } else {
                int newlb = total_label_num + PatchL(frawidx) -1;
                FL_mod(FI(frawidx)) = newlb;
                subpatchFI_dict[newlb](subpatchFcount_dict[newlb]) = FI(frawidx);
                subpatchFcount_dict[newlb] += 1; 
            }
        }
        for(auto item : subpatchFcount_dict)
        {
            int lb = item.first;
            int size = item.second;
            subpatchFI_dict[lb].conservativeResize(size);
        }
        
        total_label_num+= CC;
        std::cout << "patch " << ini_label << "cut into " << CC+1 << "pieces" << std::endl;
        // std::cout << PatchL << std::endl;
    }
}