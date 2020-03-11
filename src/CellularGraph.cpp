#include "CellularGraph.h"
#include <igl/per_vertex_normals.h>
#include <igl/vertex_triangle_adjacency.h>
namespace bcclean{

    void CCfaces_per_node(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        const std::vector<int> & node_list,
        std::map<int, std::vector<int> > & node_faces_dict
    )
    {
        // calculate node_faces_dict
        // node_list give a list of nodes indexed into Vertex
        // where key is node index and value is a std::vector<int> of face indices in Counterclock loop direction.
        node_faces_dict.clear();
        std::vector<std::vector<int> > VF, VFi;
        igl::vertex_triangle_adjacency(V, F, VF, VFi);
       
        for(auto nd: node_list)
        {
            std::vector<int> pool = VF[nd];
            int head = pool.back();
            pool.pop_back();
            node_faces_dict[nd].push_back(head);

            while(pool.size()>0)
            {
                int fidx = node_faces_dict[nd].back();
                int nd_idx = 0;

                for(int j =0; j < 3; ++j)
                {
                    // fi_list.push_back(F(fidx, j));
                    if(F(fidx, j) == nd){
                        nd_idx = j;
                        break;
                    }
                }
                // remove the next vertex in fi_list
                // u---> nd --> v  remove v here account for only the incomming edge
                std::vector<int> fi_list = {nd, F(fidx, (nd_idx + 2) %3)};
                std::sort(fi_list.begin(), fi_list.end());
                int k = 0;
                for(auto fkdx : pool)
                {
                    std::vector<int> fk_list = {F(fkdx, 0), F(fkdx, 1), F(fkdx, 2)};
                    std::sort(fk_list.begin(), fk_list.end());
                    std::vector<int> inter;
                    std::set_intersection(fk_list.begin(), fk_list.end(), fi_list.begin(), fi_list.end(), std::back_inserter(inter));
                    // inter.resize(it - inter.begin());
                    if(inter.size()==2){
                        // find the next triangle
                        // push it into node_faces_dict[nd]
                        // pop it from pool
                        node_faces_dict[nd].push_back(fkdx);
                        pool.erase(pool.begin()+k);
                    }
                    k+=1;
                }                
            }
        }
        return;

    }

    void _gen_node_CCedges_dict(
        const Eigen::MatrixXd & V_bad,
        const Eigen::MatrixXi & F_bad,
        const std::vector<edge> & edge_list,
        const std::vector<int> & node_list_bad,
        std::map<int, std::vector<int> > &  node_edge_dict   
    )
    {
        std::map<int, std::vector<int> > node_edgepool_dict;
        {
        int edg_idx =0;
        for(auto edg:edge_list)
        {
            int head= edg.head;
            int tail=edg.tail;
            assert(head != tail); // both are nodes
            assert(edg._edge_vertices.size()>1);
            for(auto ht:{head, tail})
            {
                auto it = std::find(node_edgepool_dict[ht].begin(), node_edgepool_dict[ht].end(), edg_idx);
                if(it == node_edgepool_dict[ht].end())
                {
                    // push it into list
                    node_edgepool_dict[ht].push_back(edg_idx);
                }
            }
            edg_idx +=1;
        }
        }

        // the resulting node_edgepool_dict will be non-orderd but will be used later for counter clockwise ordered dictionary
        

        // PART 0.1 get a counter clock orderd map for node_idx -> edge_idx going out of node_idx
        std::map<int, std::vector<int> > node_faces_dict_bad;
        CCfaces_per_node(V_bad, F_bad, node_list_bad, node_faces_dict_bad);
        // the  faces in node_faces_dict_bad[nd] is ordered Counter clockwise;
        for(auto item: node_faces_dict_bad)
        {
            printf("\n");
            printf("node  %d\n", item.first);
            printf("face CC: ");
            for(auto x: item.second){
                printf("%d, ", x);
            }
            printf("\n");
        }
        
        

        // also order the node_edge_dict;
        // firsly find th CC ordered out-Vertex of each nd
        std::map<int, std::vector<int> > node_outv_dict_bad;
        for(auto item: node_faces_dict_bad)
        {
            int nd = item.first;
            node_outv_dict_bad[item.first] = {};
            for(auto fidx : item.second)
            {
                std::vector<int> pool = {F_bad(fidx,0), F_bad(fidx,1), F_bad(fidx,2)};
                int bench = -1;
                int ii = 0;
                for(auto v_f: pool)
                {
                    if(v_f == nd)
                    {
                        bench = ii;
                        break;
                    }
                    ii += 1;
                }
                int outv = pool[(bench+1)%3];
                node_outv_dict_bad[nd].push_back(outv);
            }
        }

        
        for(auto item: node_outv_dict_bad)
        {
            int nd = item.first;
            node_edge_dict[nd] = std::vector<int>();
            std::vector<int> edge_pool = node_edgepool_dict[nd];
            for(auto outv : node_outv_dict_bad[nd])
            {
                for(auto edg_idx: edge_pool)
                {
                    edge  edg = edge_list[edg_idx];
                    assert(edg._edge_vertices.size()>1);
                    int second = edg._edge_vertices[1];
                    int secondlast = edg._edge_vertices[edg._edge_vertices.size()-2];
                    if(second == outv || secondlast == outv)
                    {
                        node_edge_dict[nd].push_back(edg_idx);
                    }

                }
            }
        }
    }

    CellularGraph CellularGraph::GenCellularGraph(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        const Eigen::VectorXi & FL
    )
    {
        CellularGraph cg;
        cg.V =  V;
        cg.F = F;
        cg.FL = FL; 
        cg.label_num = FL.maxCoeff()+1;
        std::vector<std::vector<int>> VF; // vertex - face index list
        {
            std::vector<std::vector<int>> VFi;
            igl::vertex_triangle_adjacency(V, F, VF,VFi);
        }
        Eigen::MatrixXd N;
        igl::per_vertex_normals(V, F, N);
        std::unordered_map<int, std::vector<int> > vertex_label_list_dict;
        build_vertex_label_list_dict(F, FL, cg.label_num, vertex_label_list_dict);
        cg._vertices.clear();
        cg._normals.clear();
        cg._nodes.clear();
        cg._edge_list.clear();
        int count = 0;
        std::map<int, int> vmap;
        std::vector<int> node_list_raw;
        for(auto item: vertex_label_list_dict)
        {
            int ss = item.second.size();
            int vidx = item.first;
            if(ss>1)
            {
                cg._vertices.push_back(V.row(vidx));
                cg._normals.push_back(N.row(vidx));
                if(ss >2)
                {
                    cg._nodes.push_back(count);
                    node_list_raw.push_back(vidx);
                }
                vmap[vidx]= count;
                count += 1;       
            }
        }
        std::vector<edge> edge_list_raw;
        build_edge_list_loop(V, F, FL, cg.label_num, edge_list_raw,cg._patch_edge_dict, cg._patch_edge_direction_dict, cg.root_cell);
        
        for(auto edg_raw: edge_list_raw)
        {
            edge edg = edg_raw;
            edg._edge_vertices.clear();
            for(auto vidx_raw: edg_raw._edge_vertices)
            {
                edg._edge_vertices.push_back(vmap.at(vidx_raw));
            }
            edg.head = vmap.at(edg_raw.head);
            edg.tail = vmap.at(edg_raw.tail);
            cg._edge_list.push_back(edg);
        }

        _gen_node_CCedges_dict(V, F, edge_list_raw, node_list_raw, cg._node_edge_dict);

        return cg;
    }

   


    
}