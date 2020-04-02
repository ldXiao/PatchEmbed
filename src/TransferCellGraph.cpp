#include "CellularGraph.h"
#include "TraceComplex.h"
#include "edge.h"
#include <igl/per_vertex_normals.h>
#include <igl/remove_unreferenced.h>
namespace bcclean{
namespace Bijection{
    
    
    void TransferCellGraph(
        const  CellularGraph & cg,
        const  MatchMaker::TraceComplex & tc,
        CellularGraph & cgt
    )
    {
        cgt.V = tc._V;
        cgt.F = tc._F;
        cgt.FL = tc._FL;
        cgt.root_cell = cg.root_cell;
        cgt.label_num = cg.label_num;
        cgt._vertices.clear();
        cgt._normals.clear();
        cgt._nodes.clear(); 
        cgt._edge_list.clear();
        cgt._patch_edge_dict = cg._patch_edge_dict; // copy
        cgt._patch_edge_direction_dict = cg._patch_edge_direction_dict; // copy
        cgt._node_edge_dict = cg._node_edge_dict; // copy
        cgt._vmap.clear();
        cgt._ivmap.clear();
        cgt._edge_list.resize(cg._edge_list.size());
        // cgt._patch_area_dict; preserve
        int count = 0; // indices into cgt._vertices
        std::map<int, int> node_record_raw;
        for(auto item: tc._edge_path_map)
        {
            edge nedg, edg;
            int edge_idx = item.first;
            std::vector<int> path_raw = item.second;
            edg  = cg._edge_list.at(edge_idx);
            nedg.head = tc._node_image_map.at(edg.head);
            nedg.tail = tc._node_image_map.at(edg.tail);
            nedg.matched = true;
            nedg.total_label_num = cgt.label_num;
            nedg._label_pair = edg._label_pair;
            nedg._edge_vertices.clear();
            nedg._edge_vertices = path_raw;
            int pos = 0;
            Eigen::MatrixXd N;
            igl::per_vertex_normals(tc._V, tc._F, N);
            for(int vidx_raw: path_raw)
            {
                if(pos == 0 || pos == path_raw.size()-1){
                    // verify if it is a node
                    if(node_record_raw.find(vidx_raw)== node_record_raw.end())
                    {
                        cgt._vertices.push_back(tc._V.row(vidx_raw));
                        cgt._nodes.push_back(count);
                        cgt._normals.push_back(N.row(vidx_raw));
                        cgt._vmap[vidx_raw] = count;
                        cgt._ivmap[count] = vidx_raw;
                        node_record_raw[vidx_raw] = count;
                        count += 1;
                    }
                }
                else{
                    cgt._vertices.push_back(tc._V.row(vidx_raw));
                    cgt._nodes.push_back(count);
                    cgt._normals.push_back(N.row(vidx_raw));
                    cgt._vmap[vidx_raw] = count;
                    cgt._ivmap[count] = vidx_raw;
                    count += 1; 
                }
                pos += 1;
            }
            cgt._edge_list[edge_idx]=nedg;
        }
    }
}
}