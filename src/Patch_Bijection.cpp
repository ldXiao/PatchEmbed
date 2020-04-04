#include <Eigen/Dense>
#include <igl/harmonic.h>
#include <igl/remove_unreferenced.h>
#include <igl/embree/line_mesh_intersection.h>
#include "patch.h"
#include "CellularGraph.h"
#include <igl/writeOBJ.h>
#include <igl/boundary_loop.h>
#include <igl/barycentric_to_global.h>
#include <igl/list_to_matrix.h>
#include <list>
#include <vector>
namespace bcclean {
namespace Bijection{
    const double PI = 3.14159265358979323846;
    bool extract_patch_mesh(
        const Eigen::MatrixXd& V, 
        const Eigen::MatrixXi& F, 
        const Eigen::MatrixXi&FL, 
        const int lb_i, 
        Eigen::MatrixXd& V_i, 
        Eigen::MatrixXi& F_i,
        Eigen::VectorXi& I_i,
        Eigen::VectorXi& J_i)
    {
        Eigen::MatrixXi F_l = Eigen::MatrixXi::Constant(F.rows(),3, 0);
        J_i = Eigen::VectorXi::Constant(F.rows(),0);
        int count = 0;
        for(int fidx=0; fidx< FL.rows(); ++fidx){
            int lb = FL(fidx,0);
            if(lb==lb_i){
                F_l.row(count)=F.row(fidx);
                J_i(count)= fidx;
                count += 1;
            }
        }
        if(count == 0){return false;}
        F_l.conservativeResize(count,3);
        J_i.conservativeResize(count);
        {   
            Eigen::VectorXi J_j;
            igl::remove_unreferenced(V, F_l, V_i, F_i,J_j, I_i);
        }
        return true;
    }

    double path_len(const std::vector<Eigen::RowVector3d> & V, 
                    const std::vector<int> & path
    )
    {
        double res=0;
        for(int p=0; p < path.size()-1; p++)
        {
            double edgeln = (V.at(path[p])- V.at(path[p+1])).norm();
            res += edgeln;
        }
        return res;
    }

    void set_bnd_uv(
        const Eigen::MatrixXd & V,
        const std::vector<int> & bnd,
        const std::map<int, bool> & bnd_nodes_dict,
        const std::vector<double> & edge_len_dict,
        Eigen::MatrixXd & bnd_uv
    )
    {
        bnd_uv = Eigen::MatrixXd::Constant(bnd.size(),3,-1);
        int cur_arc = 0;
        double shift = 0;
        int node_num = edge_len_dict.size();
        for(int idx =0 ; idx< bnd.size();++idx){
            int curr_idx = bnd[(idx) % bnd.size()];// indicis into V
            int next_idx = bnd[(idx+1) % bnd.size()];
            double x, y;
            double e_len = (V.row(curr_idx)-V.row(next_idx)).norm();
            double ratio = shift/ edge_len_dict.at(cur_arc);
            shift+= e_len;
            double curr_theta = cur_arc *2 * PI/ node_num;
            double next_theta = (cur_arc+1)*2 * PI/node_num;
            x = (1-ratio)* std::cos(curr_theta) + ratio * std::cos(next_theta);
            y = (1-ratio) * std::sin(curr_theta) + ratio * std::sin(next_theta);
            bnd_uv.row(idx) = Eigen::RowVector3d(x,y,0);
            if(bnd_nodes_dict.at(next_idx)){
                cur_arc+=1;
                shift = 0;
            }
        }
        return;  
    }

    void mapping2polygon(
        const CellularGraph & cg,
        const int pidx,
        Eigen::MatrixXd & V_uv,
        Eigen::MatrixXi & F_uv,
        std::vector<int> nodesp_uv,
        Eigen::VectorXi & VI, 
        Eigen::VectorXi & FI
    ){
        Eigen::MatrixXd Vp;
        Eigen::MatrixXi Fp;
        std::map<int, int> invIp;  
        extract_patch_mesh(cg.V, cg.F, cg.FL, pidx, Vp, Fp, VI, FI);
        for(int vidx =0 ; vidx < VI.rows(); vidx++)
        {
            int vidx_raw = VI(vidx);
            invIp[vidx_raw]= vidx;
        }
        std::vector<int> bnd, bnd_cg;
        std::map<int, bool> bnd_nodes_dict;
        std::vector<double> edge_len_list;
       
        
        std::list<int> edge_queue;
        for(int edgidx : cg._patch_edge_dict.at(pidx))
        {
            edge_queue.push_back(edgidx);
        }
        int start_edg = edge_queue.front();
        edge_queue.pop_front();
        bnd_cg = cg._edge_list.at(start_edg)._edge_vertices;
        edge_len_list.push_back(path_len(cg._vertices,cg._edge_list.at(start_edg)._edge_vertices));
        std::vector< std::vector<int>> loops;
        igl::boundary_loop(Fp, loops);
        for(int i = 0; i < bnd_cg.size(); i ++)
        {
            int vidx_cg = bnd_cg.at(i);
            bnd_nodes_dict[invIp.at(cg._ivmap.at(vidx_cg))] = (i ==0 || i == bnd_cg.size()-1);
        }
        
        while(! edge_queue.empty())
        {
            int end_raw = bnd_cg[bnd_cg.size()-1];
            int find_idx = -1; 
            for(int edgidx: edge_queue)
            {
                edge edg = cg._edge_list[edgidx];
                if(edg.head == bnd_cg.at(bnd_cg.size()-1))
                {
                    find_idx = edgidx;
                    bnd_nodes_dict[invIp.at(cg._ivmap.at(edg.head))] = true;
                    for(int i = 1; i< edg._edge_vertices.size(); ++i)
                    {
                        int vidx_cg = edg._edge_vertices.at(i);
                        bnd_cg.push_back(vidx_cg);
                        bnd_nodes_dict[invIp.at(cg._ivmap.at(vidx_cg))] = (i == edg._edge_vertices.size()-1);
                    }
                    edge_len_list.push_back(path_len(cg._vertices, edg._edge_vertices));
                    break;
                }
                else if(edg.tail == bnd_cg.at(bnd_cg.size()-1))
                {
                    find_idx = edgidx;
                    bnd_nodes_dict[invIp.at(cg._ivmap.at(edg.tail))] = true;
                    for(int i = edg._edge_vertices.size()-2; i > -1 ; i--)
                    {
                         int vidx_cg = edg._edge_vertices.at(i);
                        bnd_cg.push_back(vidx_cg);
                        bnd_nodes_dict[invIp.at(cg._ivmap.at(vidx_cg))] = (i == 0);
                    }
                    edge_len_list.push_back(path_len(cg._vertices, edg._edge_vertices));
                    break;
                }
            }
            assert(find_idx != -1);
            edge_queue.remove(find_idx);
        }
        bnd_cg.resize(bnd_cg.size()-1);
        for(auto vidx_cg: bnd_cg)
        {
            bnd.push_back(invIp.at(cg._ivmap.at(vidx_cg)));
        }

        
        
        Eigen::MatrixXd bnd_uv;
        Eigen::VectorXi bndp;
        igl::list_to_matrix(bnd, bndp);
        set_bnd_uv(Vp,bnd, bnd_nodes_dict, edge_len_list, bnd_uv);

        igl::harmonic(Vp, Fp, bndp, bnd_uv, 1, V_uv);
        F_uv = Fp;
    }
    void BijLocal(
       const Eigen::MatrixXd & Va_uv,
       const Eigen::MatrixXi & Fa,
       const Eigen::VectorXi & FIa,
       const Eigen::MatrixXd & Vb_uv,
       const Eigen::MatrixXi & Fb,
       const Eigen::VectorXi & FIb,
       Eigen::MatrixXd & M_a2b // #Va_uv * 4 matrix
    )
    {
        Eigen::MatrixXd Na = Eigen::MatrixXd::Constant(Va_uv.rows(), 3, 0);
        Eigen::MatrixXd Va_shift = Va_uv;
        Eigen::MatrixXd Vb_large = Vb_uv * 1.00000000001;
        for(int i =0 ; i< Na.rows(); ++i)
        {
            Na(i,2) = -1;
            Va_shift(i,2) = 1;
        }

        M_a2b = igl::embree::line_mesh_intersection(Va_shift, Na, Vb_large, Fb);
        for(int idx = 0; idx < Va_uv.rows(); ++idx)
        {
            int fidx = std::round(M_a2b(idx,0));
            M_a2b(idx,0) = FIb(fidx);
        }
        return;
    }

    void BijGlobal(
        const CellularGraph & cga,
        const CellularGraph & cgb,
        Eigen::MatrixXd & M_a2b // #Va * 4 matrix
    ){
        M_a2b = Eigen::MatrixXd::Constant(cga.V.rows(), 3,0);
        for(auto item: cga._patch_edge_dict)
        {
            Eigen::MatrixXd Va_uv, Vb_uv,Phi_a2b;
            Eigen::MatrixXi Fa_uv, Fb_uv;
            Eigen::VectorXi FIa, FIb, VIa, VIb;
            std::map<int, int> invVIa, invVIb;
            std::vector<int> nodesa_uv, nodesb_uv;
            int pidx = item.first;
            
            mapping2polygon(cgb, pidx, Vb_uv, Fb_uv, nodesb_uv, VIb, FIb);
            igl::writeOBJ("../dbginfo/patchb.obj", Vb_uv, Fb_uv);
            mapping2polygon(cga, pidx, Va_uv, Fa_uv, nodesa_uv, VIa, FIa);
            igl::writeOBJ("../dbginfo/patcha.obj", Va_uv, Fa_uv);
            BijLocal(Va_uv, Fa_uv, FIa, Vb_uv, Fb_uv, FIb, Phi_a2b);
            // Eigen::MatrixXd Vmap = igl::barycentric_to_global(Vb_uv, Fb_uv, Phi_a2b);
            // igl::writeOBJ("../dbginfo/map1.obj", Vmap, Fa_uv);
            for(int idx =0 ; idx< Phi_a2b.rows(); ++idx)
            {
                int fidx_raw = std::round(Phi_a2b(idx,0));
                M_a2b.row(VIa(idx)) = Phi_a2b.row(idx);
            }
        }
        return;
    }
}
}