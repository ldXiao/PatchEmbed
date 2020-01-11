#include "proj_node.h"
#include <igl/embree/line_mesh_intersection.h>
#include <igl/per_vertex_normals.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "kdtree_NN_Eigen.hpp"
#include <igl/barycentric_to_global.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/barycenter.h>
#include <queue>
namespace bcclean {

    int insertV_baryCord(
        const std::vector<int> & node_list,
        Eigen::MatrixXd & baryentry,
        Eigen::MatrixXd & V,
        Eigen::MatrixXi & F,
        Eigen::VectorXi & FL,
        Eigen::MatrixXi & TT,
        std::vector<std::vector<int> > & VV,
        std::vector<std::vector<int> > & TEdges,
        std::vector<std::vector<int> > & VEdges
    )
    {
        int fidx = std::round(baryentry(0,0));
        assert(fidx != -1);
        assert(baryentry.rows()==1);
        int v0 = F(fidx,0);
        int v1 = F(fidx,1);
        int v2 = F(fidx,2);
        int nvidx = -1;
        double u,v,w, eps;
        u = baryentry(0,1);
        v = baryentry(0,2);
        w = 1 - u - v;
        eps = 10e-2;
        std::vector<double> bcs = {w, u, v};
        int vcount =0;
        int contribute_count = 0;
        std::map<int, bool> contribute_dict;
        for(auto bc: bcs)
        {
            if(bc > eps)
            {
                contribute_count +=1;
                contribute_dict[vcount] = true;
            }
            else
            {
                contribute_dict[vcount] = false;
            }
            vcount +=1;
        }
        
        // else if(contribute_count == 2)
        // {
        //     // new node should be places on edge
        //     int non_contrib_v = 0;
        //     for(auto item: contribute_dict)
        //     {
        //         if(!item.second)
        //         {
        //             break;
        //         }
        //         non_contrib_v +=1;
        //     }

        // }
        bool addnew=true;
        if(contribute_count <3)
        {
            // only one vertices contribute to the node
            int cc  =0;
            for(auto item: contribute_dict)
            {
                if(item.second)
                {
                    break;
                }
                cc +=1;
            }
            nvidx = F(fidx, cc);
            if(std::find(node_list.begin(), node_list.end(),nvidx)!=node_list.end())
            {
                addnew = true;
            }
            else if(VEdges.empty())
            {
                addnew = false;
            }
            else if(VEdges[nvidx].empty())
            {
                addnew = false;
            }
            else
            {
                addnew = true;
            }
        }
        if(addnew)
        {
            Eigen::MatrixXd nVentry = igl::barycentric_to_global(V, F, baryentry);
            Eigen::MatrixXd nV = Eigen::MatrixXd::Zero(V.rows()+1, 3);
            nV.block(0,0, V.rows(), 3) = V;
            nV.row(V.rows()) = nVentry;
            nvidx = V.rows();
            Eigen::MatrixXi nF = Eigen::MatrixXi::Zero(F.rows()+2,3);
            FL.conservativeResize(F.rows()+2);
            
            
            nF.block(0,0,F.rows(),3) = F;
            /*
                v0

                nv
            v1      v2
            */
            // choose v0 v1 nv to inherit initial face
            nF.row(fidx) = Eigen::RowVector3i(v0, v1, nvidx);
            nF.row(F.rows()) = Eigen::RowVector3i(v1, v2, nvidx);
            nF.row(F.rows()+1) = Eigen::RowVector3i(v2, v0, nvidx);
            FL(F.rows()) = FL(fidx);
            FL(F.rows()+1) = FL(fidx);
            if(VEdges.size()!=0)
            {
                TT.conservativeResize(F.rows()+2,3);
                int nb0 = TT(fidx,0);
                int nb1 = TT(fidx,1);
                int nb2 = TT(fidx,2);
                int nf0 = fidx;
                int nf1 = F.rows();
                int nf2 = F.rows()+1;
                TT.row(nf0) = Eigen::RowVector3i(nb0, nf1, nf2);
                TT.row(nf1) = Eigen::RowVector3i(nb1, nf2, nf0);
                TT.row(nf2) = Eigen::RowVector3i(nb2, nf0, nf1);
                for(auto j: {0,1,2})
                {
                    if(TT(nb0,j)==fidx){TT(nb0,j) = nf0;}
                    if(TT(nb1,j)==fidx){TT(nb1,j) = nf1;}
                    if(TT(nb2,j)==fidx){TT(nb2,j) = nf2;}
                } 
                VV.push_back({v0,v1,v2});
                VV[v0].push_back(nvidx);
                VV[v1].push_back(nvidx);
                VV[v2].push_back(nvidx);
                VEdges.push_back(std::vector<int>());
                TEdges.push_back({TEdges[fidx][1],-1,-1});
                TEdges.push_back({TEdges[fidx][2],-1,-1});
                TEdges[fidx]={TEdges[fidx][0],-1,-1};
                
            }

            F = nF;
            
            V = nV;
        
        }

       return nvidx;
    }

    void proj_node(
        const Eigen::MatrixXd & Vbad,
        const Eigen::MatrixXi & Fbad,
        const std::vector<int> & node_list_bad, // indices into Vbad
        Eigen::MatrixXd & Vgood,
        Eigen::MatrixXi & Fgood,
        std::map<int, int> & node_map
    )
    // try to project the nodes onto the good mesh faces and use the nearest vertex of that triangle.
    // if the vertex has been taken, use the second closest vertex
    {
        node_map.clear();
        std::vector<int> node_list_good;
        Eigen::MatrixXd node_normal;
        Eigen::MatrixXd node_v;
        {
            Eigen::MatrixXd N;
            igl::per_vertex_normals(Vbad, Fbad, N);
            //
            // slice the N into node_normal;
            node_normal=Eigen::MatrixXd::Zero(N.rows(), 3);
            node_v=Eigen::MatrixXd::Zero(N.rows(), 3);
            int nd_count = 0;
            for(auto nd: node_list_bad)
            {
                node_normal.row(nd_count) = N.row(nd);
                node_v.row(nd_count) = Vbad.row(nd);
                nd_count += 1;
            }
            node_normal.conservativeResize(nd_count,3);
            node_v.conservativeResize(nd_count, 3);
        }
        Eigen::MatrixXd Vgood_copy = Vgood;

        Eigen::MatrixXd RR = igl::embree::line_mesh_intersection(node_v, node_normal, Vgood, Fgood);
        std::vector<int> node_losers;
        // if R(j,0)==-1, it means that normal does not intersect with good mesh
        for(int jj = 0; jj < RR.rows(); ++ jj)
        {
            if(std::round(RR(jj,0))!= -1)
            {
                Eigen::VectorXi FL_temp; // dummy parameter
                Eigen::MatrixXi TT_temp; // dummy parameter 
                std::vector<std::vector<int> > VV_temp;
                std::vector<std::vector<int> > VEdges_temp;
                std::vector<std::vector<int> > TEdges_temp;
                Eigen::MatrixXd bc = RR.row(jj);
                int nvidx = insertV_baryCord(node_list_good, bc, Vgood, Fgood, FL_temp, TT_temp, VV_temp, TEdges_temp, VEdges_temp);
                node_map[node_list_bad[jj]] = nvidx;
                node_list_good.push_back(nvidx);
            } 
            else
            {
                node_losers.push_back(node_list_bad[jj]);
            }
        }

        if(node_losers.size()>0)
        {
            std::map<int, int> img_records;
            kd_tree_Eigen<double> kdt(Vgood_copy.cols(),std::cref(Vgood_copy),10);
            kdt.index->buildIndex();
            for(auto nd: node_losers)
            {
                Eigen::RowVector3d querypt = Vbad.row(nd);
                int img = kd_tree_NN_Eigen( kdt, querypt);
                assert(img_records.find(img)== img_records.end());
                img_records[img] = 1;
                node_map[nd] = img;
                // TODO remove detect images from kdtree and makesure all the resulting nodes are different
            }
        }

    }

    void proj_node_loop(
        const Eigen::MatrixXd & Vbad,
        const Eigen::MatrixXi & Fbad,
        const int & node_bad, // input indices into Vbad
        const std::vector<int> & node_list_good, // nodes to exclude
        Eigen::MatrixXi & TT_good, // connectivity info
        std::vector<std::vector<int> > & VV_good,// connectivity info
        std::vector<std::vector<int> > & TEdges_good,
        std::vector<std::vector<int> > & VEdges_good, //edge vertices to exclude
        Eigen::MatrixXd & Vgood,
        Eigen::MatrixXi & Fgood,
        Eigen::VectorXi & FL_good,
        int & node_image
    )
    {
        node_image = -1;
        Eigen::MatrixXd node_normal;
        Eigen::MatrixXd node_v;
        {
            Eigen::MatrixXd N;
            igl::per_vertex_normals(Vbad, Fbad, N);
            //
            // slice the N into node_normal;
            node_normal=Eigen::MatrixXd::Zero(1, 3);
            node_v=Eigen::MatrixXd::Zero(1, 3);


            node_normal.row(0) = N.row(node_bad);
            node_v.row(0) = Vbad.row(node_bad);
        }
        Eigen::MatrixXd Vgood_copy = Vgood;
        Eigen::MatrixXd RR = igl::embree::line_mesh_intersection(node_v, node_normal, Vgood, Fgood);
        std::vector<int> node_losers;

        Eigen::MatrixXd bc = RR.row(0);
        int fidx = std::round(bc(0,0));
        std::vector<int> visit_list;
        // build a kdtree 
        Eigen::MatrixXd Centers;
        igl::barycenter(Vgood, Fgood, Centers);
        kd_tree_Eigen<double> kdt(Centers.cols(),std::cref(Centers),10);
        kdt.index->buildIndex();
        Eigen::RowVector3d query = Vbad.row(node_bad);
        int nnidx= kd_tree_NN_Eigen(kdt, query);
        double d0 = (Vbad.row(node_bad)-Vgood.row(Fgood(nnidx,0))).norm();
        double d1 = (Vbad.row(node_bad)-Vgood.row(Fgood(nnidx,1))).norm();
        double d2 = (Vbad.row(node_bad)-Vgood.row(Fgood(nnidx,2))).norm();
        double dnnbc = (Vbad.row(node_bad)-Centers.row(nnidx)).norm();
        Eigen::MatrixXd nnbc = bc;
        nnbc(0,0)= nnidx;
        double reverse_sum = (1.1/d0)+(1.1/d1)+(1.1/d2);
        nnbc(0,1)= ((1.1)/d1)/(reverse_sum);
        nnbc(0,2) =(1.1/d2)/reverse_sum;
        if(fidx == -1)
        {
            fidx = nnidx;
            bc = nnbc;
        }
        else 
        {
            double dproj = (Vbad.row(node_bad)-Centers.row(fidx)).norm();
            if(dproj > 3 * dnnbc)
            {
                // proj error is too large
                // choose nn 
                fidx = nnidx;
                bc = nnbc;
            }
        }
        if(FL_good(fidx)!= -1)
        {
            // triangle occupied
            // search for nonoccupied triangle and choose the nearset non-boundary vertice
            // we can assume for each triangle at least one vertice is not on boundary
            int target_face=-1;
            std::queue<int> search_queue;
            search_queue.push(fidx);
            while(search_queue.size()!=0){
                int cur_face = search_queue.front();
                search_queue.pop(); // remove head
                visit_list.push_back(cur_face);
                Eigen::RowVector3i adjs = TT_good.row(cur_face);
                for(int j =0 ; j <3 ; ++j){
                    int face_j = adjs(j);
                    if(face_j == -1) {
                        continue;
                    }
                    if(FL_good(face_j)!= -1)
                    {
                        if(std::find(visit_list.begin(), visit_list.end(), face_j)== visit_list.end())
                        {
                            // not visited before
                            search_queue.push(face_j);
                        }   
                    }
                    else
                    {
                        target_face = face_j;
                    }
                }
                if(target_face!= -1){break;}
            }

            // we can assume at least one of the vertices in target_face is not on Cut
            for(auto vv: {0,1,2})
            {
                if(
                    VEdges_good[Fgood(target_face,vv)].empty() 
                    && 
                    (std::find(node_list_good.begin(), node_list_good.end(),Fgood(target_face,vv))==node_list_good.end())
                )
                {
                    node_image = Fgood(target_face,vv);
                    break;
                }
            }
            assert(node_image!= -1);
            return;
        }   
        else
        {
            // triangle nonoccupied
            int nvidx = insertV_baryCord(node_list_good, bc, Vgood, Fgood, FL_good, TT_good, VV_good, TEdges_good, VEdges_good);
            // node_list_return.push_back(nvidx);
            node_image = nvidx;
        }
        
        
        
    }
}