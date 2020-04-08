#include "CellularGraph.h"
#include "TraceComplex.h"
#include "edge.h"
#include <igl/per_vertex_normals.h>
#include <igl/remove_unreferenced.h>
#include <igl/matrix_to_list.h>
#include <igl/list_to_matrix.h>
namespace bcclean{
namespace Bijection{
    
    
    void TransferCellGraph(
        const  CellularGraph & cg,
        const  MatchMaker::TraceComplex & tc,
        CellularGraph & cgt
    )
    {
        cgt.V = tc._V;
        cgt.F = Eigen::MatrixXi::Constant(tc._F.size(),3,0);
        for(int fidx= 0; fidx< tc._F.size(); ++fidx)
        {
            cgt.F.row(fidx) = tc._F[fidx];
        }
        igl::list_to_matrix(tc._FL, cgt.FL);
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
        cgt._edge_list.clear();
        // cgt._patch_area_dict; preserve
        int edge_idx = 0; // indices into cgt._vertices
        std::map<int, int> node_record_raw;
        Eigen::MatrixXd N;
        igl::per_vertex_normals(cgt.V, cgt.F, N);
        int count = 0;
        for(auto edg: cg._edge_list)
        {
            edge nedg;
            nedg.matched = true;
            nedg.total_label_num = cgt.label_num;
            nedg._label_pair = edg._label_pair;
            nedg._edge_vertices.clear();
            std::vector<int> path_raw = tc._edge_path_map.at(edge_idx);
            assert(path_raw.size()>1);
            for(int vidx_raw:path_raw)
            {
                int vidx_cg = -1;
                if(node_record_raw.find(vidx_raw)== node_record_raw.end())
                {
                    node_record_raw[vidx_raw]= count;
                    cgt._vertices.push_back(tc._V.row(vidx_raw));
                    count += 1;
                }
                vidx_cg = node_record_raw[vidx_raw];
                
                cgt._nodes.push_back(vidx_cg);
                cgt._normals.push_back(N.row(vidx_raw));
                cgt._vmap[vidx_raw] = vidx_cg;
                cgt._ivmap[vidx_cg] = vidx_raw;
                nedg._edge_vertices.push_back(vidx_cg);
            }
            nedg.head = nedg._edge_vertices[0];
            nedg.tail = nedg._edge_vertices[nedg._edge_vertices.size()-1];
            cgt._edge_list.push_back(nedg);
            edge_idx += 1;
        }
    
    }
}
}