#include "Match_Maker.h"
#include "proj_node.h"
#include "loop_colorize.h"
#include <igl/slice.h>
#include <igl/writeOBJ.h>
#include <igl/writeDMAT.h>
#include <utility>
#include <nlohmann/json.hpp>

namespace bcclean {
namespace MatchMaker{
    using json = nlohmann::json;

    void silence_vertices1(std::vector<std::vector<int > > & VV, const std::vector<int> & silent_indices)
    {
        for(auto index : silent_indices)
        {
            VV[index].clear();
            for(auto & adjs: VV)
            {
                adjs.erase(std::remove(adjs.begin(), adjs.end(), index), adjs.end());
            }
        }
    }

    void CC_faces_per_node(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        const std::vector<int> & node_list,
        std::map<int, std::vector<int> > & node_faces_dict
    )
    {
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



    void setWeights(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXd & V_bad,
        const edge & edg,
        const int & upratio,
        const int & power,
        std::vector<double> & Weights
    )
    {
        Weights.clear();
        Weights.resize(V.rows());
        std::vector<int> Elist = edg._edge_vertices;
        // arc length upsample of
        int target_num = upratio * (Elist.size()-1) + 1;
        Eigen::MatrixXd sample(target_num, 3);
        for(int jj =0 ; jj < target_num ; ++jj)
        {
            int batch = (int)(jj / upratio);
            int remain = jj % upratio;
            if(remain == 0) 
            {
                sample.row(jj) = V_bad.row(Elist[batch]);
            }
            else
            {
                double lambda = (double)(remain) / (double)(upratio);
                sample.row(jj) = (1-lambda)* V_bad.row(Elist[batch]) + lambda* V_bad.row(Elist[batch+1]);
            }
        }
        kd_tree_Eigen<double> sample_kdt(sample.cols(),std::cref(sample),10);
        sample_kdt.index->buildIndex();
        for(int j =0; j < V.rows(); ++j){
            Eigen::RowVector3d jquery= V.row(j);
            int sample_idx= kd_tree_NN_Eigen(sample_kdt, jquery);
            Weights[j] = pow((jquery - sample.row(sample_idx)).norm(), power);
        }
        
        
    }

    void sector_silence_list(
        const Eigen::MatrixXi & F,
        const std::map<int, std::vector<int> > & node_edge_dict,
        const std::map<int, std::vector<bool> > & node_edge_visit_dict,
        const std::vector<std::vector<int> > & TEdges,
        const std::map<int, std::vector<int> > & CC_node_face_dict, 
        const int& source,
        const int& cur_edge,
        std::vector<int> & temp_silence_list
    )
    {
        /*
                assume we are dealling with edge e0
                for a source assume  Counter clockwisely there are edges (ei, ebf, ek, e0, ebh) and ebf, and ebh are already set. We need to silence the faces in the sector spanned by ebh - ebf
            */

            // step 1 determin ebf and ebh
            int ebf, ebh; 
            // ebf is the right most edge in CC order to the left of e0
            // ebh is the left most edge in CC order to the right of e0
            // step 1.1 determine the position of e0
            int e0pos =0;
            for(auto ex : node_edge_dict.at(source))
            {
                if(ex != cur_edge)
                {
                    e0pos += 1;
                } 
                else
                {
                    break;
                }   
            }

            // step 1.2 determine the position of ebf
            
            int ebfpos = e0pos;
            bool find_ebf = false;
            int CCsize  = node_edge_visit_dict.at(source).size();
            for(int ii =0; ii < CCsize; ++ii)
            {
                ebfpos =(e0pos - ii + CCsize) % CCsize;
                if(node_edge_visit_dict.at(source)[ebfpos])
                {
                    find_ebf = true;
                    break;
                }
            }
            // similarly determine the position of ebh
            int ebhpos = e0pos;
            bool find_ebh = false;
            for(int ii =0; ii < CCsize; ++ii)
            {
                ebhpos =(e0pos + ii) % CCsize;
                if(node_edge_visit_dict.at(source)[ebhpos])
                {
                    find_ebh = true;
                    break;
                }
            }
            if( !find_ebf || !find_ebh)
            {
                return; // no other occupied edges
            }
            if(ebfpos == ebhpos) {
                return; // only one occupied edges
            }
            ebf = node_edge_dict.at(source)[ebfpos];
            ebh = node_edge_dict.at(source)[ebhpos];
            
             // currently occupied edge less than one nothing to go further

            // the sector is spanned CC wisely from ebh to ebf
            int sector_start = -1;
            for(auto fidx: CC_node_face_dict.at(source))
            {
                sector_start += 1;
                int e_pos = -1;
                for(auto v_pos: {0, 1, 2})
                {
                    if(F(fidx, v_pos) == source)
                    {
                        e_pos = v_pos;
                        break;
                    }
                }
                // check the out edge of the triangle fidx in TEdges
                // e_pos = v_pos
                if(TEdges[fidx][e_pos]== ebh) // sector start with ebh
                {
                    break;
                }
            }

            // current setor_start store the starting index of faces in CC_node_face_dict[source]
            for(int ff =0 ; ff < CC_node_face_dict.at(source).size() ; ++ff)
            {
                int secface = (sector_start + ff) %CC_node_face_dict.at(source).size() ;
                int e_pos = -1;
                for(auto v_pos: {0, 1, 2})
                {
                    if(F(CC_node_face_dict.at(source)[secface], v_pos) == source)
                    {
                        e_pos = v_pos;
                        break;
                    }
                }
                if(TEdges[CC_node_face_dict.at(source)[secface]][(e_pos+2)%3]== ebf) // sector end with ebf
                {
                    break;
                }
                 
                    
                temp_silence_list.push_back(F(CC_node_face_dict.at(source)[secface], (e_pos+2)%3));
                
                
            }
    }

    void update_local_sector(
        const std::vector<std::vector<int> > & VV, 
        const Eigen::MatrixXi & F,
        const std::map<int , std::vector<bool> > & node_edge_visit_dict,
        const std::map<int, std::vector<int> > & node_edge_dict,
        const std::vector<std::vector<int> > & TEdges,
        const std::map<int, std::vector<int> > & CC_node_face_dict, 
        const int& source,
        const int& target,
        const int& cur_edge,
        std::vector<std::vector<int> > & VV_temp
    )
    {
        VV_temp = VV; // copy VV into VV_temp

        {
            // set silence all other nodes in VV_temp
            std::vector<int> node_list;
            for(auto item: node_edge_dict)
            {
                if(item.first != source && item.first != target)
                {
                    node_list.push_back(item.first);
                }
            }

            silence_vertices1(VV_temp, node_list);
        }
        
        // deal with the sector of source first
        std::vector<int> temp_silence_list;
        sector_silence_list(
            F,
            node_edge_dict, 
            node_edge_visit_dict, 
            TEdges, 
            CC_node_face_dict, 
            source, 
            cur_edge, 
            temp_silence_list);

        // deal with the sctor of target
        sector_silence_list(
            F,
            node_edge_dict, 
            node_edge_visit_dict, 
            TEdges, 
            CC_node_face_dict, 
            target,
            cur_edge, 
            temp_silence_list);
        silence_vertices1(VV_temp, temp_silence_list);
    }


    bool dijkstra_trace(
        const std::vector<std::vector<int> > & VV,
        const int & source,
        const int & target,
        const std::vector<double> & Weights,
        std::vector<int> & path
    )
    {
        Eigen::VectorXi previous;
        Eigen::VectorXd min_dist;
        std::set<int> target_set = {target};

        igl::dijkstra(source, target_set, VV, Weights, min_dist, previous);


        int cur = target;
        int count = 0;
        while(cur != -1)
        {
            count +=1 ;
            path.push_back(cur);
            cur = previous(cur);
        }
        if(path[count - 1]!= source){
            // path finding failed 
            path.clear();
        }
    }


    bool split_detect(
        const Eigen::MatrixXi & F,
        const Eigen::MatrixXi & TT, // triangle-triangle adjacency
        const std::vector<int> & node_list,
        const std::vector<std::vector<int> > & VEdges, // indicate whether a vertex is on the boundary
        const std::vector<std::vector<int> > & TEdges, // indicate whether wich edge of a face is on the boundary
        std::pair<int, int> & splits
    )
    {
        // splits stores pairs of adjacent vertices indices that need to be splited
        // indiced into F

        // loop over all faces
        bool add = false;
        for(int fidx=0; fidx < F.rows(); ++fidx){
            for(int edge_idx=0; edge_idx < 3; ++ edge_idx){
                int v0 = F(fidx, edge_idx);
                int v1 = F(fidx, (edge_idx+1) % 3);
                bool nodefind0 = (std::find(node_list.begin(), node_list.end(), v0)!= node_list.end());
                bool nodefind1 = (std::find(node_list.begin(), node_list.end(), v1)!= node_list.end());

                // if at both them are node and the connecting edge is not on the cut split them
                if(nodefind0 && nodefind1 && (TEdges[fidx][edge_idx]== -1))
                {
                    int cofidx = TT(fidx, edge_idx);
                    assert(cofidx != -1);
                    if(cofidx != -1){
                        add = true;
                    }
                }
                // if one of them is node and the other is on cut
                else if((nodefind0 && VEdges[v1].size() !=0) ||(nodefind1 && VEdges[v0].size()!=0))
                {
                    if(TEdges[fidx][edge_idx]== -1)
                    {
                        int cofidx = TT(fidx, edge_idx);
                        assert(cofidx != -1);
                        if(cofidx != -1){
                            add = true;
                        }
                    }
                }
                // ff both of them are on cut and the line connecting them are not on cut
                else if((VEdges[v0].size() !=0) && (VEdges[v1].size() !=0) && (TEdges[fidx][edge_idx]== -1)){
                    int cofidx = TT(fidx, edge_idx);
                    assert(cofidx != -1);
                    if(cofidx != -1){
                        add = true;
                    }
                    // make sure each pair of indices are only pushed once
                }
                if(add)
                {
                    int vl = std::min(v0, v1);
                    int vg = std::max(v0, v1);
                    splits=std::make_pair(vl, vg);
                    return true;
                }
            }
        }
        return false;
        
    }


    bool determin_adj_configure1(
        const Eigen::MatrixXi & F, 
        const std::vector<std::vector<int > > & VF, 
        const int uidx, 
        const int vidx,
        int & fupidx,
        int & fdownidx,
        int & vupidx,
        int & vdownidx)
    {
        std::vector<int> inter(VF[uidx].size()+ VF[vidx].size());
        auto it = std::set_intersection(VF[uidx].begin(), VF[uidx].end(), VF[vidx].begin(), VF[vidx].end(), inter.begin());
        inter.resize(it-inter.begin());
        if(inter.size()!=2) return false; // there should be always only two triangles
        
        // choose  the triangle that has positive orientation on u -> v
        // should update VF in the end of loop


        /* decide up and down faces vertices
            vup
            /   \
            / fup \
        u ----- v
            \fdown/
            \   /
            vdown

            both triangle oriented outward the screen
            */
        

        for(auto trig: inter)
        {
            for(int edgepos=0; edgepos<3; ++ edgepos)
            {
                if(F(trig, edgepos) == uidx && F(trig,(edgepos+1)%3) == vidx){
                    fupidx = trig;
                    for(auto other:inter){
                        if(other!= fupidx) fdownidx = other;
                    }
                    // decide the up and down triangles


                    vupidx = F(trig, (edgepos+2)%3);
                    for(int downedge=0; downedge < 3;++downedge)
                    {
                        if(F(fdownidx, downedge)!= uidx && F(fdownidx, downedge)!= vidx) vdownidx = F(fdownidx, downedge);
                    }
                }
            }
        }
        return true;
    }

    void splits_update(
        const std::pair<int, int> & splits,
        Eigen::MatrixXd & Vraw, // raw mesh
        Eigen::MatrixXi & Fraw,
        std::vector<std::vector<int> > & VEdges, // indicate whether a vertex is on the boundary
        std::vector<std::vector<int> > & TEdges,
        std::vector<std::vector<int> > & VV // adjacency list on Fraw, posibly some ofthem has been silenced
    )
    {
        /*
        splits stores pairs of adjacent triangle indices and corresponding vertives that need to be splited

        this funciton will correspondingly  raw mesh also the face mapping
        
        everytime split two faces and create four smaller triangles add two of them to the end of Fraw Fbase

        also add the face mapping and cut info into the end of FI, VEdges, TEdges
        */

        const int task_num = 1;
        // each split task will increace face num by 2 and vertices num by 1
        const int nF_num_raw = task_num * 2 + Fraw.rows();
        const int nV_num_raw = task_num + Vraw.rows();

        // preallocate all the memory
        Eigen::MatrixXd nVraw = Eigen::MatrixXd::Zero(nV_num_raw , 3); // raw mesh
        nVraw.block(0, 0, Vraw.rows(), 3) = Vraw;
        Eigen::MatrixXi nFraw = Eigen::MatrixXi::Zero(nF_num_raw, 3);
        nFraw.block(0, 0, Fraw.rows(), 3) = Fraw;


        std::vector<std::vector<int> > nVEdges(nV_num_raw); // indicate whether a vertex is on the boundary
        std::vector<std::vector<int> > nTEdges(nF_num_raw);
        VV.resize(nV_num_raw);
        for(int count =0; count <nV_num_raw; ++ count)
        {
            if(count < VEdges.size()) nVEdges[count] = VEdges[count];
            else nVEdges[count] = std::vector<int>();
        }
        for(int fcount = 0; fcount  < nF_num_raw; ++fcount)
        {   
            if(fcount < TEdges.size()) nTEdges[fcount] = TEdges[fcount];
            else nTEdges[fcount] = {-1, -1, -1};
        }
        std::vector<std::vector<int > > nVF_raw;// #V list of lists of incident faces (adjacency list)
        // the naming indicates that they are to be updated by hand in each loop
        {
            std::vector<std::vector<int > > II;  //  #V list of lists of index of incidence within incident faces listed in VF, local variable do not care about them
            igl::vertex_triangle_adjacency(Vraw, Fraw, nVF_raw, II);
        }
        // initialize VFs

        // because splits might be connected directly, this process has to be done one by one
        // after  spliting each edge, the info in splits needs not to be updated because
        // splitting does not change initial indices of existing vertices
        int task_count = 0;
        int uidx_raw = splits.first;
        int vidx_raw = splits.second;
        int widx_raw = Vraw.rows() + task_count;
        Eigen::RowVector3d  wpos = (Vraw.row(uidx_raw)+Vraw.row(vidx_raw))/2;
        nVraw.row(widx_raw) = wpos;
        // they are two representations of an identical vertex
        

        // decide the adj configurations
        /* decide up and down faces vertices
            vup
            /   \
          / fup \
        u ----- v
          \fdown/
            \   /
            vdown

            both triangle oriented outward the screen
            */

        int fupidx_raw, fdownidx_raw;
        int vupidx_raw, vdownidx_raw;
        determin_adj_configure1(
            nFraw, nVF_raw, uidx_raw, vidx_raw,
            fupidx_raw, fdownidx_raw, vupidx_raw, vdownidx_raw
        );
        // decide the indices of existing triangles and vertices

        


        // create new triangles

        /* 
                vup
                / | \
                /  |  \
              / f1|f0 \
            /    |    \
            u-----w ----v
            \    |    /
             \ f2| f3/
              \  |  /
                \ | /
                vdown
        */
        int f0idx_raw, f1idx_raw, f2idx_raw, f3idx_raw;
        f0idx_raw = fupidx_raw;
        f2idx_raw = fdownidx_raw;
        f1idx_raw = Fraw.rows();
        f3idx_raw = f1idx_raw + 1;
        
        nVEdges[widx_raw] = std::vector<int>();



        // update nTEdges
        std::map<std::pair<int,int> , int> nTEdges_record;
        // store cut info of all 5 edges of fupidx_raw, and fdownidx_raw in nTEdges resulted from last loop
        for(int edgepos =0 ; edgepos < 3 ; ++edgepos){
            for(auto updown: {fupidx_raw, fdownidx_raw}){
                int uu = nFraw(updown, edgepos);
                int vv = nFraw(updown, (edgepos+1)%3);
                int vl = std::min(uu,vv);
                int vg = std::max(uu, vv);
                nTEdges_record[std::make_pair(vl, vg)] = nTEdges[updown][edgepos];
            }
        }
        std::pair<int, int> key; 
        key = std::make_pair(std::min(vidx_raw, vupidx_raw), std::max(vidx_raw, vupidx_raw));
        nTEdges[f0idx_raw] = {-1, nTEdges_record[key], -1};

        key = std::make_pair(std::min(uidx_raw, vupidx_raw), std::max(uidx_raw, vupidx_raw));
        nTEdges[f1idx_raw] = {-1, nTEdges_record[key], -1};

        key = std::make_pair(std::min(uidx_raw, vdownidx_raw), std::max(uidx_raw, vdownidx_raw));
        nTEdges[f2idx_raw] = {-1, nTEdges_record[key], -1};

        key = std::make_pair(std::min(vidx_raw, vdownidx_raw), std::max(vidx_raw, vdownidx_raw));
        nTEdges[f3idx_raw] = {-1, nTEdges_record[key], -1};
        


        // update nFraw the faces with correct orientation
        nFraw.row(f0idx_raw) = Eigen::RowVector3i(widx_raw, vidx_raw, vupidx_raw);
        nFraw.row(f1idx_raw) = Eigen::RowVector3i(widx_raw, vupidx_raw, uidx_raw);
        nFraw.row(f2idx_raw) = Eigen::RowVector3i(widx_raw, uidx_raw, vdownidx_raw);
        nFraw.row(f3idx_raw) = Eigen::RowVector3i(widx_raw, vdownidx_raw, vidx_raw);
        
    


        
        // update nVFs VV the connectivitity info
        // update nVF_raw;
        /*  only have to deal with f1 f3 becase they are new faces
                vup
               / | \
              /  |  \
             / f1|f0 \
            /    |    \
            u-----w ----v
            \    |    /
             \ f2| f3/
              \  |  /
                \ | /
                vdown
        */
        nVF_raw.resize(nV_num_raw);
        
        nVF_raw[vidx_raw].push_back(f3idx_raw);
        nVF_raw[vidx_raw].erase(std::remove(nVF_raw[vidx_raw].begin(), nVF_raw[vidx_raw].end(), f2idx_raw), nVF_raw[vidx_raw].end()); 

        nVF_raw[vupidx_raw].push_back(f1idx_raw);

        nVF_raw[uidx_raw].push_back(f1idx_raw);
        nVF_raw[uidx_raw].erase(std::remove(nVF_raw[uidx_raw].begin(), nVF_raw[uidx_raw].end(), f0idx_raw), nVF_raw[uidx_raw].end());

        nVF_raw[vdownidx_raw].push_back(f3idx_raw);
        nVF_raw[widx_raw] = {vidx_raw, vupidx_raw, uidx_raw, vdownidx_raw};



        for(auto idx_raw: {vidx_raw, vupidx_raw, uidx_raw, vdownidx_raw}){

            // idx_raw has not been silenced add connection
            VV[idx_raw].push_back(widx_raw);
            VV[widx_raw].push_back(idx_raw);
            // add all coonection into VV and silence the cuts afterwards
        }

        // also silence the connection between vidx_raw and uidx_raw
        // if one of them are silenced no need to further remove the other
        VV[vidx_raw].erase(std::remove(VV[vidx_raw].begin(), VV[vidx_raw].end(), uidx_raw), VV[vidx_raw].end()); 
        VV[uidx_raw].erase(std::remove(VV[uidx_raw].begin(), VV[uidx_raw].end(), vidx_raw), VV[uidx_raw].end());

        // the above two lines will do nothing if they find nothing




        Vraw = nVraw;
        Fraw = nFraw;
        VEdges = nVEdges;
        TEdges = nTEdges;
    }

    void trace_and_label(
        const Eigen::MatrixXd & V_bad,
        const Eigen::MatrixXi & F_bad,
        const Eigen::VectorXi & FL_bad,
        Eigen::MatrixXd & V_good,
        Eigen::MatrixXi & F_good,
        Eigen::VectorXi & FL_good,
        int & kk
    )
    {
        std::map<int, int> edge_order_map; // store and maintain the order of added edges {order: edge_dx}
        std::map<int, std::vector<int> > edge_path_map; // {edge_idx, path}
        int run_count = 0;
        // Randomize Seed
        srand(static_cast<unsigned int>(time(nullptr)));
        int total_label_num = FL_bad.maxCoeff()+1;


        // PART 0 GET THE FRAME ON BAD MESH
        std::vector<bcclean::edge> edge_list;
        std::unordered_map<int, std::vector<int> > patch_edge_dict;
        std::unordered_map<int, std::vector<bool> > patch_edge_direction_dict;
        // label -> list(edge_idx) map indices into edge_list
        build_edge_list_loop(V_bad, F_bad, FL_bad, total_label_num, edge_list, patch_edge_dict, patch_edge_direction_dict);

        // we can assume that  all nodes have valance more than 3
        std::vector<int> node_list_bad;
        {
            std::unordered_map<int, std::vector<int> > vertex_label_list_dict;
            build_vertex_label_list_dict(F_bad, FL_bad, total_label_num, vertex_label_list_dict);
            for(auto item: vertex_label_list_dict)
            {
                if(item.second.size()>2)
                {
                    node_list_bad.push_back(item.first);
                }
            }
        }

        
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
        CC_faces_per_node(V_bad, F_bad, node_list_bad, node_faces_dict_bad);
        // the  faces in node_faces_dict_bad[nd] is ordered Counter clockwise;
        for(auto item: node_faces_dict_bad)
        {
            printf("\n");
            printf("node  %d\n", item.first);
            printf("face CC: ");
            for(auto x: item.second){
                printf("%d, ", x);
            }
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

        std::map<int, std::vector<int> > node_edge_dict;
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
                    auto it = std::find(edg._edge_vertices.begin(), edg._edge_vertices.end(), outv);
                    auto jt = std::find(node_edge_dict[nd].begin(), node_edge_dict[nd].end(), edg_idx);
                    if(it != edg._edge_vertices.end() and jt == node_edge_dict[nd].end())
                    {
                        node_edge_dict[nd].push_back(edg_idx);
                    }

                }
            }
        }

        // PART 1 finds the Nearest Neighbor of nodes on good mesh
        // build the kdtree;
        kd_tree_Eigen<double> kdt(V_good.cols(),std::cref(V_good),10);
        kdt.index->buildIndex();
        std::map<int, int> node_imge_dict;
        std::vector<int> node_list_good;
        bcclean::proj_node(V_bad, F_bad, node_list_bad, V_good, F_good, node_imge_dict);
        
        for(auto item:node_imge_dict)
        {
            node_list_good.push_back(item.second);
        }
        std::vector<int> sorted_node_list_bad = node_list_bad;
        std::map<int , std::vector<bool> > node_edge_visit_dict;
        std::sort(sorted_node_list_bad.begin(), sorted_node_list_bad.end(), [node_edge_dict](int nda, int ndb){return (node_edge_dict.at(nda).size() > node_edge_dict.at(ndb).size());});
        // shift sorted_node_list by kk
        std::vector<int> sorted_node_list_bad_copy = sorted_node_list_bad;
        int sss = sorted_node_list_bad.size();
        for(int jj =0; jj < sss; ++jj)
        {
            sorted_node_list_bad[jj] = sorted_node_list_bad_copy[(jj+kk)%sss];
        }
        // remove all the above one finish the true implementation
        for(auto nd: sorted_node_list_bad)
        {
            for(auto q: node_edge_dict[nd])
            {
                node_edge_visit_dict[nd].push_back(false);
            }
        }


        // we will use only adjacency list VV, bool lists VEdges, TEdges
        // initializations
        std::vector<std::vector<int> > VV_good, VF_good;
        igl::adjacency_list(F_good, VV_good);
        {
            std::vector<std::vector<int> >  VFi_good;
            igl::vertex_triangle_adjacency(V_good, F_good, VF_good, VFi_good);

        }
        std::vector<std::vector<int> > VEdges_good(V_good.rows());
        std::vector<std::vector<int> > TEdges_good(F_good.rows());
        for(int count =0; count <V_good.rows(); ++ count)
        {
            VEdges_good[count] = std::vector<int>();
        }
        for(int fcount = 0; fcount  < F_good.rows(); ++fcount)
        {   
            TEdges_good[fcount] = {-1, -1, -1};
        }
        // start with the nodes with largest valance and deal with the edge starting with this node in counter clock order
        std::vector<int> total_silence_list;
        for(auto nd: sorted_node_list_bad)
        {
            // Main loop for tracing
            int edg_nd = -1; // the indices of target edge in node_edge_dict[nd]

            for (auto edge_idx: node_edge_dict[nd]){
                edg_nd +=1;
                if(node_edge_visit_dict[nd][edg_nd])
                {
                    continue;
                }
                run_count+=1;

                Eigen::MatrixXi TT_good;
                igl::triangle_triangle_adjacency(F_good, TT_good);
                {
                    std::vector<std::vector<int> > VFi_good;
                    igl::vertex_triangle_adjacency(V_good, F_good, VF_good, VFi_good);
                }
                edge edg = edge_list[edge_idx];
                int target =-1;
                int target_bad = -1;
                int source_bad = -1;
                // if(edg.head == nd){target_bad = edg.tail; target = node_imge_dict[target_bad];}
                // if(edg.tail == nd){target_bad= edg.head; target= node_imge_dict[edg.head];}
                // int source = node_imge_dict[nd];
                int source = node_imge_dict[edg.head];
                source_bad = edg.head;
                target = node_imge_dict[edg.tail];
                target_bad = edg.tail;
                assert(target != -1 && source != -1); 

                std::vector<double> Weights;
                // update local sector
                // (a) create CC_node_face_list
                std::map<int, std::vector<int> > CC_node_face_dict;
                std::vector<std::vector<int > > VV_temp;
                CC_faces_per_node(V_good, F_good, {source, target}, CC_node_face_dict);
                {
                    std::map<int, std::vector<bool> > node_edge_visit_dict_temp;
                    std::map<int, std::vector<int> > node_edge_dict_temp;
                    for(auto item: node_edge_dict)
                    {
                        node_edge_dict_temp[node_imge_dict.at(item.first)] = node_edge_dict.at(item.first);
                        node_edge_visit_dict_temp[node_imge_dict.at(item.first)] = node_edge_visit_dict.at(item.first);
                    }
                    update_local_sector(
                        VV_good, 
                        F_good, 
                        node_edge_visit_dict_temp, 
                        node_edge_dict_temp,
                        TEdges_good,
                        CC_node_face_dict,
                        source,
                        target,
                        edge_idx,
                        VV_temp);
                }
                setWeights(V_good, V_bad, edg, 10, 1,  Weights);
                // the Weights is vertex based

                
                // dijkstra_trace(....,VEdges, TEdges);
                std::vector<int> path;
                dijkstra_trace(VV_temp, source, target, Weights, path);
                std::vector<int> path_records(path.size()-2);
                
                for(int p =0 ; p < path.size()-2; ++p)
                {
                    path_records[p] = path[p+1];
                }


                // path update VV
                silence_vertices1(VV_good,path_records);
                for(auto rec: path_records)
                {
                    total_silence_list.push_back(rec);
                }
                //path update VEdges
                for(auto vidx:path_records){
                    VEdges_good[vidx].push_back(edge_idx);
                }

                // path updates TEdges
                // setf the triangle edges in cuts to be true
                for(int rc_idx=0; rc_idx < path.size()-1; ++rc_idx){
                    int uidx = path[rc_idx];
                    int vidx = path[(rc_idx+1)%path.size()];
                    std::vector<int> inter(VF_good[uidx].size()+ VF_good[vidx].size());
                    auto it = std::set_intersection(VF_good[uidx].begin(), VF_good[uidx].end(), VF_good[vidx].begin(), VF_good[vidx].end(), inter.begin());
                    inter.resize(it-inter.begin());
                    // there should be only one comman adjacent triangle for boundary vertices
                    for(auto trg: inter){
                        for(int edgpos =0; edgpos < 3 ; ++edgpos){
                            int uuidx = F_good(trg, edgpos);
                            int vvidx = F_good(trg, (edgpos+1)% 3);
                            if(uuidx == uidx && vvidx == vidx){
                                TEdges_good[trg][edgpos] = edge_idx;//1457-1459 1587 1588 2189 2371 2373 2374 2716 2742 1794 2712 3046 3047
                                // break;
                            }
                            if(uuidx == vidx && vvidx == uidx){
                                TEdges_good[trg][edgpos] = edge_idx;
                                // break;
                            }
                        }
                    }
                }



                // split_detect
                std::pair<int, int> splits;
                while(split_detect(F_good, TT_good, node_list_good,VEdges_good, TEdges_good, splits))
                {
                    

                    // splits_update
                    splits_update(splits, V_good, F_good, VEdges_good, TEdges_good, VV_good);
                    igl::triangle_triangle_adjacency(F_good, TT_good);
                }
                silence_vertices1(VV_good, total_silence_list);

                json path_json;
                edge_order_map[run_count] = edge_idx;
                edge_path_map[edge_idx] = path;
                for(auto item: edge_order_map)
                {
                    path_json[std::to_string(item.second)] = edge_path_map[item.second];   
                }
                std::ofstream file;
                file.open("../debug_paths.json");
                file << path_json;

                igl::writeOBJ("../debug_mesh.obj", V_good, F_good);

                // update visit_dict or loop condition update
                node_edge_visit_dict[nd][edg_nd] = true;
                int pos =0;
                for(auto tedge_idx: node_edge_dict[target_bad]){
                    if(tedge_idx == edge_idx)
                    {
                        node_edge_visit_dict[target_bad][pos]=true;
                    }
                    pos +=1;
                }
                pos = 0;
                for(auto sedge_idx: node_edge_dict[source_bad]){
                    if(sedge_idx == edge_idx)
                    {
                        node_edge_visit_dict[source_bad][pos]=true;
                    }
                    pos +=1;
                } 

            }

            silence_vertices1(VV_good,{node_imge_dict[nd]});
            // total_silence_list.push_back(nd);
        }
        std::vector<std::vector<int> > finalVF, finalVFi;
        FL_good = Eigen::VectorXi::Constant(F_good.rows(),-1);
        igl::vertex_triangle_adjacency(V_good, F_good, finalVF, finalVFi);
        for(auto item: patch_edge_dict)
        {
            int lb = item.first;
            int edge_idx = item.second[0];
            bool correct_direction = patch_edge_direction_dict[lb][0];
            int va = edge_path_map[edge_idx][0];
            int vb = edge_path_map[edge_idx][1];
            std::vector<int> fas = finalVF.at(va);
            std::vector<int> fbs = finalVF.at(vb);
            std::vector<int> inter(10);
            std::vector<int>::iterator it;
            it=std::set_intersection (fas.begin(),fas.end(), fbs.begin(), fbs.end(), inter.begin());
            inter.resize(it-inter.begin());
            int fa = inter[0];
            int fb = inter[1];
            bool fa_on_left = false;
            int seed_face = -1;
            for(auto feidx:{0,1,2})
            {
                if(F_good(fa,feidx) == va && F_good(fa, (feidx +1)% 3) == vb)
                {
                    fa_on_left = true;
                    break;
                }
            }
            if(correct_direction != fa_on_left)
            {
                seed_face = fa;
            }
            else 
            {
                seed_face = fb;
            }
            loop_colorize(V_good, F_good, TEdges_good, seed_face, lb, FL_good);
        }
        igl::writeDMAT("../FL_final.dmat", FL_good);
        
        
    }
}
}
