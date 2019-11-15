#include <Eigen/Core>
#include <vector>
#include <igl/bfs.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include <Eigen/Sparse>
#include <algorithm>
namespace bcclean{
    void silence_vertices(std::vector<std::vector<int > > & VV, const std::vector<int> & silent_indices){
        for(auto index : silent_indices){
            VV[index].clear();
            for(auto & adjs: VV){
                adjs.erase(std::remove(adjs.begin(), adjs.end(), index), adjs.end());
            }
        }
    }

    bool determin_adj_configure(
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

     void restricted_splits_detect(
        const Eigen::MatrixXi & Fraw,
        const Eigen::MatrixXi & TT, // triangle-triangle adjacency
        const std::vector<bool> & VCuts, // indicate whether a vertex is on the boundary
        const std::vector<std::vector<bool> > & TCuts, // indicate whether wich edge of a face is on the boundary
        std::map<std::pair<int, int>, int > & splits
    )
    {
        // splits stores pairs of adjacent triangle indices that need to be splited
        // indiced into Fraw
        splits.clear();
        // loop over all faces
        for(int fidx=0; fidx < Fraw.rows(); ++fidx){
            for(int edge_idx=0; edge_idx < 3; ++ edge_idx){
                int v0 = Fraw(fidx, edge_idx);
                int v1 = Fraw(fidx, (edge_idx+1) % 3);
                if(VCuts[v0] && VCuts[v1] && !TCuts[fidx][edge_idx]){
                    int cofidx = TT(fidx, edge_idx);
                    assert(cofidx != -1);
                    if(cofidx != -1){
                        int vl = std::min(v0, v1);
                        int vg = std::max(v0, v1);
                        splits[std::make_pair(vl, vg)] = 1;
                    }
                    // make sure each pair of indices are only pushed onece
                }
            }
        }
        return;
        
    }

    void restricted_splits_update(
        const std::map<std::pair<int, int>, int> & splits,
        Eigen::MatrixXd & Vraw, // raw mesh
        Eigen::MatrixXi & Fraw,
        Eigen::MatrixXd & Vbase, // base mesh
        Eigen::MatrixXi & Fbase,
        Eigen::VectorXi & VI, // vertex mapping Vraw -> Vbase
        Eigen::VectorXi & FI, // face mapping Fraw -> Fbase
        Eigen::VectorXi & FL, // face label on Fbase
        std::vector<bool> & VCuts, // indicate whether a vertex is on the boundary
        std::vector<std::vector<bool> > & TCuts,
        std::vector<std::vector<int> > & VV // adjacency list on Fraw, posibly some ofthem has been silenced
    )
    {
        /*
        splits stores pairs of adjacent triangle indices and corresponding vertives that need to be splited

        this funciton will correspondingly modily base and raw mesh also the face mapping
        add new vertices to the tail of Vraw and Vbase
        
        everytime split two faces and create four smaller triangles add two of them to the end of Fraw Fbase

        also add the face mapping and cut info into the end of FI, VCuts, TCuts
        */

        const int task_num = splits.size();
        // each split task will increace face num by 2 and vertices num by 1
        const int nF_num_raw = task_num * 2 + Fraw.rows();
        const int nF_num_base = task_num * 2 + Fbase.rows();
        const int nV_num_raw = task_num + Vraw.rows();
        const int nV_num_base = task_num + Vbase.rows();

        // preallocate all the memory
        Eigen::MatrixXd nVraw = Eigen::MatrixXd::Zero(nV_num_raw , 3); // raw mesh
        nVraw.block(0, 0, Vraw.rows(), 3) = Vraw;
        Eigen::MatrixXi nFraw = Eigen::MatrixXi::Zero(nF_num_raw, 3);
        nFraw.block(0, 0, Fraw.rows(), 3) = Fraw;
        Eigen::MatrixXd nVbase = Eigen::MatrixXd::Zero(nV_num_base,3); // base mesh
        nVbase.block(0, 0, Vbase.rows(), 3) = Vbase;
        Eigen::MatrixXi nFbase= Eigen::MatrixXi::Zero(nF_num_base, 3);
        nFbase.block(0, 0, Fbase.rows(), 3) = Fbase;
        Eigen::VectorXi nFI = Eigen::VectorXi::Zero(nF_num_raw); // face mapping
        nFI.block(0,0, FI.rows(),1) = FI;
        Eigen::VectorXi nVI = Eigen::VectorXi::Zero(nV_num_raw);
        nVI.block(0,0, VI.rows(),1) = VI;
        Eigen::VectorXi nFL = Eigen::VectorXi::Zero(nF_num_base); // face mapping
        nFL.block(0,0, FL.rows(),1) = FL;
        std::vector<bool> nVCuts(nV_num_raw); // indicate whether a vertex is on the boundary
        std::vector<std::vector<bool> > nTCuts(nF_num_raw);
        VV.resize(nV_num_raw);
        for(int count =0; count <nV_num_raw; ++ count)
        {
            if(count < VCuts.size()) nVCuts[count] = VCuts[count];
            else nVCuts[count] = false;
        }
        for(int fcount = 0; fcount  < nF_num_raw; ++fcount)
        {   
            if(fcount < TCuts.size()) nTCuts[fcount] = TCuts[fcount];
            else nTCuts[fcount] = {false, false, false};
        }
        std::vector<std::vector<int > > nVF_raw, nVF_base;// #V list of lists of incident faces (adjacency list)
        // the naming indicates that they are to be updated by hand in each loop
        {
            std::vector<std::vector<int > > II, JJ;  //  #V list of lists of index of incidence within incident faces listed in VF, local variable do not care about them
            igl::vertex_triangle_adjacency(Vraw, Fraw, nVF_raw, II);
            igl::vertex_triangle_adjacency(Vbase, Fbase, nVF_base, JJ);
        }
        // initialize VFs

        // because splits might be connected directly, this process has to be done one by one
        // after  spliting each edge, the info in splits needs not to be updated because
        // splitting does not change initial indices of existing vertices
        int task_count = 0;
        for(auto item: splits){
            int uidx_raw = item.first.first;
            int vidx_raw = item.first.second;
            int uidx_base = VI(uidx_raw); // the indices in splits can only be existing vertices
            int vidx_base = VI(vidx_raw);
            int widx_raw = Vraw.rows() + task_count;
            int widx_base = Vbase.rows() + task_count;
            Eigen::RowVector3d  wpos = (Vraw.row(uidx_raw)+Vraw.row(vidx_raw))/2;
            nVraw.row(widx_raw) = wpos;
            nVbase.row(widx_base) = wpos;
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
            determin_adj_configure(
                nFraw, nVF_raw, uidx_raw, vidx_raw,
                fupidx_raw, fdownidx_raw, vupidx_raw, vdownidx_raw
            );

            int fupidx_base, fdownidx_base;
            int vupidx_base, vdownidx_base;
            determin_adj_configure(
                nFbase, nVF_base, uidx_base, vidx_base,
                fupidx_base, fdownidx_base, vupidx_base, vdownidx_base
            );
            // decide the indices of existing triangles and vertces

            


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
            f1idx_raw = Fraw.rows() + 2 * task_count;
            f3idx_raw = f1idx_raw + 1;
            // update TCuts before updating nFraw
            //  poseone updates of nFraw after TCuts because Tcuts relies on Fraw
            // update nVCuts
            // only one vertices is added
            nVCuts[widx_raw] = false;



            // update nTCuts
            std::map<std::pair<int,int> , bool> nTCuts_record;
            // store cut info of all 5 edges of fupidx_raw, and fdownidx_raw in nTCuts resulted from last loop
            for(int edgepos =0 ; edgepos < 3 ; ++edgepos){
                for(auto updown: {fupidx_raw, fdownidx_raw}){
                    int uu = nFraw(updown, edgepos);
                    int vv = nFraw(updown, (edgepos+1)%3);
                    int vl = std::min(uu,vv);
                    int vg = std::max(uu, vv);
                    nTCuts_record[std::make_pair(vl, vg)] = nTCuts[updown][edgepos];
                }
            }
            std::pair<int, int> key; 
            key = std::make_pair(std::min(vidx_raw, vupidx_raw), std::max(vidx_raw, vupidx_raw));
            nTCuts[f0idx_raw] = {false, nTCuts_record[key], false};

            key = std::make_pair(std::min(uidx_raw, vupidx_raw), std::max(uidx_raw, vupidx_raw));
            nTCuts[f1idx_raw] = {false, nTCuts_record[key], false};

            key = std::make_pair(std::min(uidx_raw, vdownidx_raw), std::max(uidx_raw, vdownidx_raw));
            nTCuts[f2idx_raw] = {false, nTCuts_record[key], false};

            key = std::make_pair(std::min(vidx_raw, vdownidx_raw), std::max(vidx_raw, vdownidx_raw));
            nTCuts[f3idx_raw] = {false, nTCuts_record[key], false};
            


            // update nFraw the faces with correct orientation
            nFraw.row(f0idx_raw) = Eigen::RowVector3i(widx_raw, vidx_raw, vupidx_raw);
            nFraw.row(f1idx_raw) = Eigen::RowVector3i(widx_raw, vupidx_raw, uidx_raw);
            nFraw.row(f2idx_raw) = Eigen::RowVector3i(widx_raw, uidx_raw, vdownidx_raw);
            nFraw.row(f3idx_raw) = Eigen::RowVector3i(widx_raw, vdownidx_raw, vidx_raw);


            int f0idx_base, f1idx_base, f2idx_base, f3idx_base;
            f0idx_base = fupidx_base;
            f2idx_base = fdownidx_base;
            f1idx_base = Fbase.rows() + 2 * task_count;
            f3idx_base = f1idx_base + 1;
            // initialize the faces with correct orientation
            nFbase.row(f0idx_base) = Eigen::RowVector3i(widx_base, vidx_base, vupidx_base);
            nFbase.row(f1idx_base) = Eigen::RowVector3i(widx_base, vupidx_base, uidx_base);
            nFbase.row(f2idx_base) = Eigen::RowVector3i(widx_base, uidx_base, vdownidx_base);
            nFbase.row(f3idx_base) = Eigen::RowVector3i(widx_base, vdownidx_base, vidx_base);

            
            


            // updates on nFs are done

            // update the nVI nFI  nFL
            nVI(widx_raw) = widx_base;
            
            nFI(f0idx_raw) = f0idx_base;
            nFI(f1idx_raw) = f1idx_base;
            nFI(f2idx_raw) = f2idx_base;
            nFI(f3idx_raw) = f3idx_base;

            nFL(f0idx_base) = nFL(fupidx_base);
            nFL(f1idx_base) = nFL(fupidx_base);
            nFL(f2idx_base) = nFL(fdownidx_base);
            nFL(f3idx_base) = nFL(fdownidx_base);


            
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


            nVF_base.resize(nV_num_base);

            nVF_base[vidx_base].push_back(f3idx_base);
            nVF_base[vidx_base].erase(std::remove(nVF_base[vidx_base].begin(), nVF_base[vidx_base].end(), f2idx_base), nVF_base[vidx_base].end());


            nVF_base[vupidx_base].push_back(f1idx_base);

            nVF_base[uidx_base].push_back(f1idx_base);
            nVF_base[uidx_base].erase(std::remove(nVF_base[uidx_base].begin(), nVF_base[uidx_base].end(), f0idx_base), nVF_base[uidx_base].end());

            nVF_base[vdownidx_base].push_back(f3idx_base);
            nVF_base[widx_base] = {vidx_base, vupidx_base, uidx_base, vdownidx_base};

            for(auto idx_raw: {vidx_raw, vupidx_raw, uidx_raw, vdownidx_raw}){
                if(VV[idx_raw].size()!=0){
                    // idx_raw has not been silenced add connection
                    VV[idx_raw].push_back(widx_raw);
                    VV[widx_raw].push_back(idx_raw);
                }
            }

            // also silence the connection between vidx_raw and uidx_raw
            // if one of them are silenced no need to further remove the other
            VV[vidx_raw].erase(std::remove(VV[vidx_raw].begin(), VV[vidx_raw].end(), uidx_raw), VV[vidx_raw].end()); 
            VV[uidx_raw].erase(std::remove(VV[uidx_raw].begin(), VV[uidx_raw].end(), vidx_raw), VV[uidx_raw].end());

            // the above two lines will do nothing if they find nothing

            task_count += 1; // finishes one task
        }


        Vraw = nVraw;
        Fraw = nFraw;
        Vbase = nVbase;
        Fbase = nFbase;
        VI = nVI;
        FI = nFI;
        FL = nFL;
        VCuts = nVCuts;
        TCuts = nTCuts;
    }

    void planar_cut_simply_connect(
        Eigen::MatrixXd & V, // Vraw
        Eigen::MatrixXi & F,  // Fraw
        Eigen::MatrixXd & Vbase,
        Eigen::MatrixXi & Fbase, 
        Eigen::VectorXi & VI, // vertex mapping Vraw -> Vbase
        Eigen::VectorXi & FI, // face mapping Fraw -> Fbase
        Eigen::VectorXi & FL, // face label on Fbase
        const std::vector<std::vector<int> > & boundary_loops, 
        std::vector<bool> & VCuts, 
        std::vector<std::vector<bool> > & TCuts)
    {
        // V #V x 3 vertices of a nonsimply connected mesh
        // F #F x 3 faces of the mesh
        //  std::vector<int> & VLoopLabels, #V x 1 indicate whether a vertice is on a loop if not label = -1 if on label = loop index ( loop index starting from 0)
        // return 
        //Vcuts # V x 1 indicate a vertice is on the cut or boundary.
        // TCuts # #F * 3 indicate an edge in triangle is in cut or boundary


        // resize VCuts
        std::vector<int> VLoops(V.rows());
        assert(boundary_loops.size()>1);
        VCuts.resize(V.rows());
        TCuts.resize(F.rows());
        for(int count =0; count <V.rows(); ++ count)
        {
            VCuts[count]= false;
        }
        for(int fcount = 0; fcount  < F.rows(); ++fcount)
        {
            TCuts[fcount] = {false, false, false};
        }

        std::vector<std::vector<int > > VF;// #V list of lists of incident faces (adjacency list)
        std::vector<std::vector<int > > VFi;  //  #V list of lists of index of incidence within incident faces listed in VF
        igl::vertex_triangle_adjacency(V, F, VF, VFi);
        for(auto loop: boundary_loops){
            for(int loop_idx=0; loop_idx < loop.size(); ++loop_idx){
                int uidx = loop[loop_idx];
                int vidx = loop[(loop_idx+1)%loop.size()];
                std::vector<int> inter(VF[uidx].size()+ VF[vidx].size());
                auto it = std::set_intersection(VF[uidx].begin(), VF[uidx].end(), VF[vidx].begin(), VF[vidx].end(), inter.begin());
                inter.resize(it-inter.begin());
                // there should be only one comman adjacent triangle for boundary vertices
                for(auto trg: inter){
                    for(int edgpos =0; edgpos < 3 ; ++edgpos){
                        int uuidx = F(trg, edgpos);
                        int vvidx = F(trg, (edgpos+1)%3);
                        if(uuidx == uidx && vvidx == vidx){
                            TCuts[trg][edgpos] = true;
                        }
                        if(uuidx == vidx && vvidx == uidx){
                            TCuts[trg][edgpos] = true;
                        }
                    }
                }
            }
        }
        Eigen::MatrixXi TT;
        igl::triangle_triangle_adjacency(F, TT);
        
        // initialize VLoops and VCuts;
        int loop_num = 0;
        for(auto loop: boundary_loops){
            for(auto vidx: loop){
                VCuts[vidx]=true;
                VLoops[vidx]=loop_num;
            }
            loop_num += 1;
        }

        std::vector<std::vector<int> > VV;
        igl::adjacency_list(F, VV); // initialize the adjacency list;

        // start from loop 0
        // start with the first by default
        std::map<int, int> loop_visit_count;
        std::vector<int> remaining_loops;
        for(int i =0; i < boundary_loops.size();++i){
            remaining_loops.push_back(i);
            loop_visit_count[i]=0;
        }
        int curloop = 0;
        int root=boundary_loops[curloop][0];
        loop_visit_count[curloop]=1;



        std::map<std::pair<int, int>, int > splits;
        restricted_splits_detect(F, TT, VCuts, TCuts, splits);
        if(splits.size()>0)
        {
            restricted_splits_update(splits, V, F, Vbase, Fbase, VI, FI, FL, VCuts, TCuts, VV);
        }
        while(remaining_loops.size() != 0)
        {
            std::vector<int > D, P;
            igl::bfs(VV, root, D, P); // D is the spanning tree in the discover order and P store the predessors, all contain indices into VV
            int tail;
            for(auto vidx: D){
                if(VCuts[vidx]){
                    if(VLoops[vidx] ==0 && remaining_loops.size()>1){
                        continue; 
                    }
                    if(VLoops[vidx]== -1){continue;}
                    if(loop_visit_count[VLoops[vidx]] == 2){
                        continue;
                    }
                    tail = vidx;
                    break;
                } else{
                    continue;
                }
            }
            curloop = VLoops[tail];
            loop_visit_count[curloop] += 1;
            if(loop_visit_count[curloop]==2){
                silence_vertices(VV, boundary_loops[curloop]);
            }
            int cur  = tail;
            std::vector<int> records;
            igl::vertex_triangle_adjacency(V, F, VF, VFi); // update VF
            while(cur != -1){
                // root not reach, go on back search
                records.push_back(cur);
                VCuts[cur] = true;
                cur = P[cur];

            }
            // remove the indices in records in the adjacency info
            silence_vertices(VV, records);

            // set the triangle edges in cut to be true
            for(int rc_idx=0; rc_idx < records.size()-1; ++rc_idx){
                int uidx = records[rc_idx];
                int vidx = records[(rc_idx+1)%records.size()];
                std::vector<int> inter(VF[uidx].size()+ VF[vidx].size());
                auto it = std::set_intersection(VF[uidx].begin(), VF[uidx].end(), VF[vidx].begin(), VF[vidx].end(), inter.begin());
                inter.resize(it-inter.begin());
                // there should be only one comman adjacent triangle for boundary vertices
                for(auto trg: inter){
                    for(int edgpos =0; edgpos < 3 ; ++edgpos){
                        int uuidx = F(trg, edgpos);
                        int vvidx = F(trg, (edgpos+1)% 3);
                        if(uuidx == uidx && vvidx == vidx){
                            TCuts[trg][edgpos] = true;
                            // break;
                        }
                        if(uuidx == vidx && vvidx == uidx){
                            TCuts[trg][edgpos] = true;
                            // break;
                        }
                    }
                }
            }
            // find next root;
            int mid = int(boundary_loops[curloop].size() / 2);
            int tail_idx=-1;
            for(auto vidx :boundary_loops[curloop]){
                tail_idx+=1;
                if(vidx == tail){
                    break;
                }
            }
            root = boundary_loops[curloop][(tail_idx+mid) %  boundary_loops[curloop].size()];
            // visited the loop twice, pop it out from remaining_loops
            if(curloop != 0){loop_visit_count[curloop] +=1;}
            if(loop_visit_count[curloop]>=2){
                remaining_loops.erase(std::remove(remaining_loops.begin(), remaining_loops.end(),curloop), remaining_loops.end());
            }

            igl::triangle_triangle_adjacency(F, TT); 
            restricted_splits_detect(F, TT, VCuts, TCuts, splits);
            if(splits.size()>0){
                restricted_splits_update(splits, V, F, Vbase, Fbase, VI, FI, FL, VCuts, TCuts, VV);
            }
        }
    }

   

}