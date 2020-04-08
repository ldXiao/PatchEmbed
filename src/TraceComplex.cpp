#include "TraceComplex.h"
#include "helper.h"
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/barycentric_to_global.h>
#include <algorithm>
#include <igl/doublearea.h>
#include <igl/writeOBJ.h>
namespace bcclean
{
namespace MatchMaker
{
       bool determine_adj_configure(
        const std::vector<Eigen::RowVector3i> & F, 
        const std::vector<std::vector<int > > & VF, 
        const int uidx, 
        const int vidx,
        int & fupidx,
        int & fdownidx,
        int & vupidx,
        int & vdownidx)
    {
        std::vector<int> inter(VF[uidx].size()+ VF[vidx].size());
        std::vector<int> vtr = VF[vidx];
        std::vector<int> utr = VF[uidx];
        std::sort(vtr.begin(), vtr.end());
        std::sort(utr.begin(), utr.end());
        auto it = std::set_intersection(vtr.begin(), vtr.end(), utr.begin(), utr.end(), inter.begin());
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
                Eigen::RowVector3i left;
                left = F[trig];
                if(F[trig](edgepos) == uidx && F[trig]((edgepos+1)%3) == vidx){
                    fupidx = trig;
                    for(auto other:inter){
                        if(other!= fupidx) fdownidx = other;
                    }
                    // decide the up and down triangles


                    vupidx = F[trig]( (edgepos+2)%3);
                    for(int downedge=0; downedge < 3;++downedge)
                    {
                        if(F[fdownidx](downedge)!= uidx && F[fdownidx](downedge)!= vidx) vdownidx = F[fdownidx](downedge);
                    }
                }
            }
        }
        return true;
    }

    void TraceComplex::initialize(const Eigen::MatrixXd & Vin, const Eigen::MatrixXi & Fin)
    {
        Helper::to_list(Vin, this->_V);
        this->_F.resize(Fin.rows());
        Eigen::VectorXi DblA;
        igl::doublearea(Vin, Fin, DblA);
        Eigen::MatrixXi TT;
        igl::triangle_triangle_adjacency(Fin, TT);
        _TT.clear();
        for(int fidx=0; fidx < Fin.rows(); ++ fidx)
        {
            _TT.push_back(TT.row(fidx));
            _DblA.push_back(DblA(fidx));
            _FL.push_back(-1);
            _F[fidx]= Fin.row(fidx);
        }
        igl::adjacency_list(Fin, _VV);
        {
            std::vector<std::vector<int> >  VFi;
            igl::vertex_triangle_adjacency(Vin, Fin, _VF, VFi);
        }
        _VEdges.resize(Vin.rows());
        for(int count =0; count <Vin.rows(); ++ count)
        {
            _VEdges[count] = std::vector<int>();
        }
        _TEdges.resize(Fin.rows());
        for(int fcount = 0; fcount  < Fin.rows(); ++fcount)
        {   
            _TEdges[fcount] = {-1, -1, -1};
        }
    }

    void TraceComplex::insert_update(
        const Eigen::MatrixXd & baryentry
    )
    {
        int fidx = std::round(baryentry(0,0));
        int v0 = _F[fidx](0);
        int v1 = _F[fidx](1);
        int v2 = _F[fidx](2);
        assert(fidx != -1);
        assert(baryentry.rows()==1);
        Eigen::RowVector3d nVentry;
        {
            Eigen::MatrixXi F;
            Eigen::MatrixXd V;
            Helper::to_matrix(_F, F);
            Helper::to_matrix(_V, V);
            nVentry = igl::barycentric_to_global(V, F, baryentry).row(0);
        }
        std::vector<Eigen::RowVector3d> nV;
        nV = _V;
        nV.resize(_V.size()+1);
        nV[_V.size()] = nVentry;
        int nvidx = _V.size();
        _FL.resize(_F.size()+2);
        _DblA.resize(_F.size()+2);
        // update splits_record
        std::vector<double> record = {_V[v0](0), _V[v0](1), _V[v0](2), _V[v1](0), _V[v1](1), _V[v1](2), _V[v2](0),_V[v2](1),_V[v2](2), nVentry(0), nVentry(1), nVentry(2)};
        
        _splits_record.push_back(record);
        std::vector<Eigen::RowVector3i> nF = _F;
        nF.resize(_F.size()+2);
        /*
            v0

            nv
        v1      v2
        */
        // choose v0 v1 nv to inherit initial face
        nF[fidx] = Eigen::RowVector3i(v0, v1, nvidx);
        nF[_F.size()] = Eigen::RowVector3i(v1, v2, nvidx);
        nF[_F.size()+1] = Eigen::RowVector3i(v2, v0, nvidx);
        _FL[_F.size()] = _FL[fidx];
        _FL[_F.size()+1] = _FL[fidx];
        

        // update _DblA

        double a0, a1,a2;
        a0 = ((nV[v0]-nVentry).cross(nV[v1]-nVentry)).norm();
        a1 = ((nV[v1]-nVentry).cross(nV[v2]-nVentry)).norm();
        a2 = ((nV[v2]-nVentry).cross(nV[v0]-nVentry)).norm();
        _DblA[fidx] = a0;
        _DblA[_F.size()] = a1;
        _DblA[_F.size()+1] = a2;
        _FL[_F.size()] =-1;
        _FL[_F.size()+1] = -1;

        // update TT
        _TT.resize(_F.size()+2);
        int nb0 = _TT[fidx](0);
        int nb1 = _TT[fidx](1);
        int nb2 = _TT[fidx](2);
        int nf0 = fidx;
        int nf1 = _F.size();
        int nf2 = _F.size()+1;
        _TT[nf0] = Eigen::RowVector3i(nb0, nf1, nf2);
        _TT[nf1] = Eigen::RowVector3i(nb1, nf2, nf0);
        _TT[nf2] = Eigen::RowVector3i(nb2, nf0, nf1);
        for(auto j: {0,1,2})
        {
            if(_TT[nb0](j)==fidx)
            {
                _TT[nb0](j) = nf0;
            }
            if(_TT[nb1](j)==fidx)
            {
                _TT[nb1](j) = nf1;
            }
            if(_TT[nb2](j)==fidx)
            {
                _TT[nb2](j) = nf2;
            }
        } 

        // update VV
        _VV.push_back(std::vector<int>());
        for(auto vi: {v0, v1, v2})
        {
            if(_VEdges[vi].empty())
            {
                // if vi is not on cut
                _VV[nvidx].push_back(vi);
                _VV[vi].push_back(nvidx);
            }
        }
        // update VF
        _VF.push_back({nf0, nf1, nf2});
        _VF[v0].push_back(nf2);
        _VF[v1].push_back(nf1);
        _VF[v2].push_back(nf1);
        _VF[v2].push_back(nf2);
        _VF[v2].erase(std::remove(_VF[v2].begin(), _VF[v2].end(), nf0), _VF[v2].end()); 

        // update _VEdges _TEdges
        _VEdges.push_back(std::vector<int>());
        _TEdges.push_back({_TEdges[fidx][1],-1,-1});
        _TEdges.push_back({_TEdges[fidx][2],-1,-1});
        _TEdges[fidx]={_TEdges[fidx][0],-1,-1}; 


        _F = nF;
        _V = nV;
        // igl::writeOBJ("../dbginfo/operation"+std::to_string(this->operation_count)+".obj", _V, _F);
        std::ofstream split_file;
        split_file.open("../dbginfo/splits_record.txt");
        for(auto & vec: this->_splits_record)
        {
            if(vec.size()==9)
            {
                // it is splits
                split_file << 0;
                split_file << " ";
            }
            else
            {
                // it is an insert
                split_file << 1;
                split_file << " ";
            }
            
            for(auto & val: vec)
            {
                split_file << val;
                split_file << " ";
            }
            split_file << "\n";
        }
        this->operation_count+=1;
    }


    void TraceComplex::splits_detect(std::vector<std::pair<int, int> >& splits)
    {
        splits.clear();
        for(int fidx=0; fidx < _F.size(); ++fidx){
            for(int edge_idx=0; edge_idx < 3; ++ edge_idx){
                bool add = false;
                int v0 = _F[fidx]( edge_idx);
                int v1 = _F[fidx]( (edge_idx+1) % 3);
                bool nodefind0 = (std::find(_node_list.begin(), _node_list.end(), v0)!= _node_list.end());
                bool nodefind1 = (std::find(_node_list.begin(), _node_list.end(), v1)!= _node_list.end());

                // if at both them are node and the connecting edge is not on the cut split them
                if(nodefind0 && nodefind1 && (_TEdges[fidx][edge_idx]== -1))
                {
                    int cofidx = _TT[fidx]( edge_idx);
                    assert(cofidx != -1);
                    if(cofidx != -1){
                        add = true;
                    }
                }
                // if one of them is node and the other is on cut
                else if((nodefind0 && _VEdges[v1].size() !=0) ||(nodefind1 && _VEdges[v0].size()!=0))
                {
                    if(_TEdges[fidx][edge_idx]== -1)
                    {
                        int cofidx = _TT[fidx]( edge_idx);
                        assert(cofidx != -1);
                        if(cofidx != -1){
                            add = true;
                        }
                    }
                }
                // ff both of them are on cut and the line connecting them are not on cut
                else if((_VEdges[v0].size() !=0) && (_VEdges[v1].size() !=0) && (_TEdges[fidx][edge_idx]== -1)){
                    int cofidx = _TT[fidx]( edge_idx);
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
                    std::pair<int, int> temp_pair =std::make_pair(vl, vg);
                    if(std::find(splits.begin(), splits.end(), temp_pair) == splits.end())
                    {
                        splits.push_back(temp_pair);
                    }

                }
            }
        }

    }

    void TraceComplex::split_update(const std::pair<int,int> & split){
       /*
        split stores pairs of adjacent triangle indices and corresponding vertives that need to be splited

        this funciton will correspondingly  raw mesh also the face mapping
        
        everytime split two faces and create four smaller triangles add two of them to the end of _F Fbase

        also add the face mapping and cut info into the end of FI, _VEdges, _TEdges
        */

        // each split task will increace face num by 2 and vertices num by 1
        const int nF_num = 2 + _F.size();
        const int nV_num = 1 + _V.size();
        const int oldF_num = _F.size();
        const int oldV_num = _V.size();
        // preallocate all the memory
        std::vector<Eigen::RowVector3d> nV = _V; // raw mesh
        nV.resize(nV_num);
        std::vector<Eigen::RowVector3i> nF = _F;
        nF.resize(_F.size()+2);
        


        // splitting does not change initial indices of existing vertices
        int uidx = split.first;
        int vidx = split.second;
        int widx = oldV_num;
        
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

        int fupidx, fdownidx;
        int vupidx, vdownidx;
        {
            bool nobug = determine_adj_configure(
                _F, _VF, uidx, vidx,
                fupidx, fdownidx, vupidx, vdownidx
            );
            assert(nobug);
        }
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
        int f0idx, f1idx, f2idx, f3idx;
        // choose the convention to preserve f0 and f2
        f0idx = fupidx;
        f2idx = fdownidx;
        f1idx = oldF_num;
        f3idx = oldF_num + 1;
        
        
        // update _VEdges
        _VEdges.resize(nV_num); // nV_num = oldV_num + 1
        _VEdges[widx] = std::vector<int>(); // widx == oldV_num



        // update TEdges and TT, these has to be updated before _F
        _TEdges.resize(nF_num);
        _TT.resize(nF_num);
        std::map<std::pair<int,int> , int> TEdges_record;
        std::map<std::pair<int, int>, int> TT_record;
        // store cut info of all 5 edges of fupidx, and fdownidx in _TEdges resulted from last loop before we update Tdges
        for(int edgepos =0 ; edgepos < 3 ; ++edgepos){
            for(auto updown: {fupidx, fdownidx}){
                int uu = _F[updown]( edgepos);
                int vv = _F[updown]( (edgepos+1)%3);
                int vl = std::min(uu,vv);
                int vg = std::max(uu, vv);
                TEdges_record[std::make_pair(vl, vg)] = _TEdges[updown][edgepos];
                TT_record[std::make_pair(vl, vg)] = _TT[updown]( edgepos);
            }
        }
        std::pair<int, int> key0, key1, key2, key3; 
        key0 = std::make_pair(std::min(vidx, vupidx), std::max(vidx, vupidx));
        _TEdges[f0idx] = {-1, TEdges_record[key0], -1};
        _TT[f0idx] = Eigen::RowVector3i(f3idx, TT_record[key0], f1idx);
        
        key1 = std::make_pair(std::min(uidx, vupidx), std::max(uidx, vupidx));
        _TEdges[f1idx] = {-1, TEdges_record[key1], -1};
        _TT[f1idx] = Eigen::RowVector3i(f0idx, TT_record[key0], f2idx);

        key2 = std::make_pair(std::min(uidx, vdownidx), std::max(uidx, vdownidx));
        _TEdges[f2idx] = {-1, TEdges_record[key2], -1};
        _TT[f2idx] = Eigen::RowVector3i(f1idx, TT_record[key2], f3idx);

        key3 = std::make_pair(std::min(vidx, vdownidx), std::max(vidx, vdownidx));
        _TEdges[f3idx] = {-1, TEdges_record[key3], -1};
        _TT[f3idx] = Eigen::RowVector3i(f2idx, TT_record[key3], f0idx);
        



    
        // update _FL;
        _FL.resize(nF_num);
        _FL[f1idx] =_FL[f0idx]; 
        _FL[f3idx] =_FL[f2idx];


         // update _F the faces with correct orientation
        nF[f0idx] = Eigen::RowVector3i(widx, vidx, vupidx);
        nF[f1idx] = Eigen::RowVector3i(widx, vupidx, uidx);
        nF[f2idx] = Eigen::RowVector3i(widx, uidx, vdownidx);
        nF[f3idx] = Eigen::RowVector3i(widx, vdownidx, vidx);
        _F = nF;
        // update _V
        Eigen::RowVector3d  wpos = (_V[uidx]+_V[vidx])/2;
        nV[widx] = wpos;
        _V = nV;
        std::vector<double> record = {_V[uidx](0), _V[uidx](1), _V[uidx](2), _V[vidx](0), _V[vidx](1), _V[vidx](2) ,_V[widx](0),_V[widx](1), _V[widx](2)};
        _splits_record.push_back(record);
        
        // iglw::writeOBJ("../dbginfo/operation"+std::to_string(this->operation_count)+".obj", _V, _F);
        std::ofstream split_file;
        split_file.open("../dbginfo/splits_record.txt");
        for(auto & vec: this->_splits_record)
        {
            if(vec.size()==9)
            {
                // it is splits
                split_file << 0;
                split_file << " ";
            }
            else
            {
                // it is an insert
                split_file << 1;
                split_file << " ";
            }
            
            for(auto & val: vec)
            {
                split_file << val;
                split_file << " ";
            }
            split_file << "\n";
        }
        this->operation_count++;
        
        // update _DblA;
        _DblA.resize(nF_num);
        _DblA[f0idx] = ((nV[vupidx]-wpos).cross(nV[vidx]-wpos)).norm();
        _DblA[f1idx] = ((nV[vupidx]-wpos).cross(nV[uidx]-wpos)).norm();
        _DblA[f2idx] = ((nV[vdownidx]-wpos).cross(nV[uidx]-wpos)).norm();
        _DblA[f3idx] = ((nV[vdownidx]-wpos).cross(nV[vidx]-wpos)).norm(); 
        
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
        // update VF
        _VF.resize(nV_num);
        _VF[vidx].push_back(f3idx);
        _VF[vidx].erase(std::remove(_VF[vidx].begin(), _VF[vidx].end(), f2idx), _VF[vidx].end()); 

        _VF[vupidx].push_back(f1idx);

        _VF[uidx].push_back(f1idx);
        _VF[uidx].erase(std::remove(_VF[uidx].begin(), _VF[uidx].end(), f0idx), _VF[uidx].end());

        _VF[vdownidx].push_back(f3idx);
        _VF[widx] = {f0idx, f1idx, f2idx,f3idx};


        // update VV
        _VV.resize(nV_num);
        for(auto idx: {vidx, vupidx, uidx, vdownidx}){

            // idx_raw has not been silenced add connection
            _VV[idx].push_back(widx);
            _VV[widx].push_back(idx);
            // add all connection into VV and silence the vertices in total_silence_list afterwards
        }

        // also silence the connection between vidx_raw and uidx_raw
        // if one of them is already silenced no need to further remove the other
        _VV[vidx].erase(std::remove(_VV[vidx].begin(), _VV[vidx].end(), uidx), _VV[vidx].end()); 
        _VV[uidx].erase(std::remove(_VV[uidx].begin(), _VV[uidx].end(), vidx), _VV[uidx].end());

        // the above two lines will do nothing if they find nothing
        
        
    }
}
} // namespace bcclean
