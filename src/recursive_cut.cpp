#include "recursive_cut.h"
#include "non_vertex_manifold_relabel.h"
#include "planar_cut_simply_connect.h"
#include "patch_cut_relabel.h"
#include <igl/matrix_to_list.h>
#include <igl/remove_unreferenced.h>
#include <igl/is_vertex_manifold.h>
#include <igl/boundary_loop.h>
#include <igl/writeOBJ.h>
namespace bcclean{
namespace Preprocess{
    void recursive_cut(
        Eigen::MatrixXd & Vbase,
        Eigen::MatrixXi & Fbase, 
        Eigen::VectorXi & FLbase, // face label on Fbas
        std::map<int, Eigen::MatrixXi> & label_faces_dict,
        std::map<int, Eigen::VectorXi> & label_FI_dict,
        const int lb
    )
    {
        std::cout  <<"analyzing patch:" << lb << std::endl;
        Eigen::MatrixXi Fraw;
        Eigen::MatrixXd Vraw;
        Eigen::VectorXi I, J;
        igl::remove_unreferenced(Vbase, label_faces_dict[lb], Vraw, Fraw, I, J);
        // first step remove nonmanifold vertices
        Eigen::VectorXi NonManifoldVertices;
        bool is_manifold = igl::is_vertex_manifold(Fraw, NonManifoldVertices);
        std::vector<std::vector<int > > boundary_loops;
        igl::boundary_loop(Fraw, boundary_loops);
        int total_label_num = label_FI_dict.size();
        if(!is_manifold)
        {
            // relabel around on each nonmaniofld branch
            std::cout << "patch "<< lb << " is non-manfold, first level relabelling" << std::endl;
            std::vector<int> NMV;
            // igl::writeOBJ("nm.obj", Vraw, Fraw);
            for(int j =0 ; j < NonManifoldVertices.rows(); j++)
            {
                if(NonManifoldVertices(j)==0)
                {
                    NMV.push_back(j);
                }
            }
            std::map<int, Eigen::VectorXi> subpatchFI_dict;
            Eigen::VectorXi FLbase_copy;
            Prepocess::non_vertex_manifold_relabel(Vraw, Fraw, label_FI_dict[lb], NMV, FLbase_copy,lb, FLbase, total_label_num, subpatchFI_dict);
            for(auto item: subpatchFI_dict)
            {
                int lb_sub = item.first;
                Eigen::VectorXi FI_sub = item.second;
                Eigen::MatrixXi F_sub = Eigen::MatrixXi::Constant(Fraw.rows(),3, -1);
                Eigen::MatrixXd V_sub;
                int count = 0;
                for(int f =0 ; f < FI_sub.rows(); ++f)
                {
                    F_sub.row(f) = Fbase.row(FI_sub(f));
                    count +=1;
                }
                F_sub.conservativeResize(count,3);
                // check Vraw, Fraw is a manifold
                // check the number of boundary loops of the manifold
                // if more than one it is not simply connected

                // {
                //     Eigen::MatrixXi F_sub_copy = F_sub;
                //     Eigen::VectorXi I,J;
                //     igl::remove_unreferenced(Vbase, F_sub_copy, V_sub, F_sub, I,J);
                // }
                label_faces_dict[lb_sub]= F_sub;
                label_FI_dict[lb_sub] = FI_sub;
                recursive_cut(Vbase, Fbase, FLbase, label_faces_dict, label_FI_dict, lb_sub);
            }

        }
        else if(boundary_loops.size()>1)
        {
            std::vector<bool> VCuts;
            std::vector<std::vector<bool> > TCuts;
            std::cout << "patch" << lb << "has " << boundary_loops.size() <<"loops"<< std::endl;
            Eigen::VectorXi I, J;
            Eigen::VectorXi VI, FI;
            FI = label_FI_dict[lb];
            igl::remove_unreferenced(Vbase, label_faces_dict[lb], Vraw, Fraw, I, J);
            igl::writeOBJ("sub.obj", Vraw, Fraw);
            VI = J;
            std::vector<int> VIII;
            VIII = igl::matrix_to_list(VI);
            
            planar_cut_simply_connect(Vraw, Fraw, Vbase, Fbase, VI, FI, FLbase, boundary_loops, VCuts, TCuts);
            // std::cout << "ob1" << lb << std::endl;
            Eigen::VectorXi FLbase_copy = FLbase;
            std::map<int, Eigen::VectorXi > subpatchFI_dict;
            patch_cut_relabel(Fraw, FI, VCuts, TCuts, FLbase_copy, FLbase, subpatchFI_dict);
            for(auto item: subpatchFI_dict)
            {
                int lb_sub = item.first;
                label_FI_dict[lb_sub] = item.second;
            }
            
            for(auto item: subpatchFI_dict)
            {
                int lb_sub = item.first;
                Eigen::VectorXi FI_sub = item.second;
                Eigen::MatrixXi F_sub = Eigen::MatrixXi::Constant(Fraw.rows(),3, -1);
                Eigen::MatrixXd V_sub;
                std::vector<int> FIII = igl::matrix_to_list(FI_sub);
                int count = 0;
                for(int f =0 ; f < FI_sub.rows(); ++f)
                {
                    F_sub.row(f) = Fbase.row(FI_sub(f));
                    count +=1;
                }
                F_sub.conservativeResize(count,3);
                // check Vraw, Fraw is a manifold
                // check the number of boundary loops of the manifold
                // if more than one it is not simply connected

                // {
                //     Eigen::MatrixXi F_sub_copy = F_sub;
                //     Eigen::VectorXi I,J;
                //     igl::remove_unreferenced(Vbase, F_sub_copy, V_sub, F_sub, I,J);
                //     igl::writeOBJ("sub.obj", V_sub, F_sub);
                // }
                label_faces_dict[lb_sub]= F_sub;
                recursive_cut(Vbase, Fbase, FLbase, label_faces_dict, label_FI_dict, lb_sub);
            }
            
        }
        else
        {
            std::cout << "patch" << lb << "is a disk now" << std::endl;  
            return;  
        }
    }
}
}