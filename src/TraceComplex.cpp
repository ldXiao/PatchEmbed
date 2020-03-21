#include "TraceComplex.h"
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/adjacency_list.h>
#include <igl/barycentric_to_global.h>
#include <algorithm>
#include <igl/doublearea.h>
namespace bcclean
{
namespace MatchMaker
{
       bool determine_adj_configure(
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
                left = F.row(trig);
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

    void TraceComplex::initialize(const Eigen::MatrixXd & Vin, const Eigen::MatrixXi & Fin)
    {
        this->_V = Vin;
        this->_F = Fin;
        this->_FL = Eigen::VectorXi::Constant(Fin.rows(),-1);
        igl::doublearea(Vin, Fin, _DblA);
        igl::triangle_triangle_adjacency(Fin, _TT);
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
        int v0 = _F(fidx,0);
        int v1 = _F(fidx,1);
        int v2 = _F(fidx,2);
        assert(fidx != -1);
        assert(baryentry.rows()==1);
        Eigen::RowVector3d nVentry = igl::barycentric_to_global(_V, _F, baryentry).row(0);
        Eigen::MatrixX3d nV = Eigen::MatrixXd::Zero(_V.rows()+1, 3);
        nV.block(0,0, _V.rows(), 3) = _V;
        nV.row(_V.rows()) = nVentry;
        int nvidx = _V.rows();
        Eigen::MatrixXi nF = Eigen::MatrixXi::Zero(_F.rows()+2,3);
        _FL.conservativeResize(_F.rows()+2);
        _DblA.conservativeResize(_F.rows()+2);
        // update splits_record
        std::vector<double> record = {_V(v0,0), _V(v0,1), _V(v0,2), _V(v1,0), _V(v1,1), _V(v1,2), _V(v2,0),_V(v2,1),_V(v2,2)};
        _splits_record.push_back(record);
        nF.block(0,0,_F.rows(),3) = _F;
        /*
            v0

            nv
        v1      v2
        */
        // choose v0 v1 nv to inherit initial face
        nF.row(fidx) = Eigen::RowVector3i(v0, v1, nvidx);
        nF.row(_F.rows()) = Eigen::RowVector3i(v1, v2, nvidx);
        nF.row(_F.rows()+1) = Eigen::RowVector3i(v2, v0, nvidx);
        _FL(_F.rows()) = _FL(fidx);
        _FL(_F.rows()+1) = _FL(fidx);
        

        // update _DblA

        double a0, a1,a2;
        a0 = ((nV.row(v0)-nVentry).cross(nV.row(v1)-nVentry)).norm();
        a1 = ((nV.row(v1)-nVentry).cross(nV.row(v2)-nVentry)).norm();
        a2 = ((nV.row(v2)-nVentry).cross(nV.row(v0)-nVentry)).norm();
        _DblA(fidx) = a0;
        _DblA(_F.rows()) = a1;
        _DblA(_F.rows()+1) = a2;

        // update TT
        _TT.conservativeResize(_F.rows()+2,3);
        int nb0 = _TT(fidx,0);
        int nb1 = _TT(fidx,1);
        int nb2 = _TT(fidx,2);
        int nf0 = fidx;
        int nf1 = _F.rows();
        int nf2 = _F.rows()+1;
        _TT.row(nf0) = Eigen::RowVector3i(nb0, nf1, nf2);
        _TT.row(nf1) = Eigen::RowVector3i(nb1, nf2, nf0);
        _TT.row(nf2) = Eigen::RowVector3i(nb2, nf0, nf1);
        for(auto j: {0,1,2})
        {
            if(_TT(nb0,j)==fidx){_TT(nb0,j) = nf0;}
            if(_TT(nb1,j)==fidx){_TT(nb1,j) = nf1;}
            if(_TT(nb2,j)==fidx){_TT(nb2,j) = nf2;}
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
        
    }


    void TraceComplex::splits_detect(std::vector<std::pair<int, int> >& splits)
    {
        splits.clear();
        for(int fidx=0; fidx < _F.rows(); ++fidx){
            for(int edge_idx=0; edge_idx < 3; ++ edge_idx){
                bool add = false;
                int v0 = _F(fidx, edge_idx);
                int v1 = _F(fidx, (edge_idx+1) % 3);
                bool nodefind0 = (std::find(_node_list.begin(), _node_list.end(), v0)!= _node_list.end());
                bool nodefind1 = (std::find(_node_list.begin(), _node_list.end(), v1)!= _node_list.end());

                // if at both them are node and the connecting edge is not on the cut split them
                if(nodefind0 && nodefind1 && (_TEdges[fidx][edge_idx]== -1))
                {
                    int cofidx = _TT(fidx, edge_idx);
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
                        int cofidx = _TT(fidx, edge_idx);
                        assert(cofidx != -1);
                        if(cofidx != -1){
                            add = true;
                        }
                    }
                }
                // ff both of them are on cut and the line connecting them are not on cut
                else if((_VEdges[v0].size() !=0) && (_VEdges[v1].size() !=0) && (_TEdges[fidx][edge_idx]== -1)){
                    int cofidx = _TT(fidx, edge_idx);
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
        const int nF_num = 2 + _F.rows();
        const int nV_num = 1 + _V.rows();
        const int oldF_num = _F.rows();
        const int oldV_num = _V.rows();
        // preallocate all the memory
        Eigen::MatrixX3d nV = Eigen::MatrixXd::Zero(nV_num , 3); // raw mesh
        nV.block(0, 0, oldV_num, 3) = _V;
        Eigen::MatrixXi nF = Eigen::MatrixXi::Zero(nF_num, 3);
        nF.block(0, 0, oldF_num, 3) = _F;
        


        // splitting does not change initial indices of existing vertices
        int uidx = split.first;
        int vidx = split.second;
        int widx = oldV_num;
        std::vector<double> record = {_V(uidx,0), _V(uidx,1), _V(uidx,2), _V(vidx,0), _V(vidx,1), _V(vidx,2) };
        _splits_record.push_back(record);
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
        _TT.conservativeResize(nF_num,3);
        std::map<std::pair<int,int> , int> TEdges_record;
        std::map<std::pair<int, int>, int> TT_record;
        // store cut info of all 5 edges of fupidx, and fdownidx in _TEdges resulted from last loop before we update Tdges
        for(int edgepos =0 ; edgepos < 3 ; ++edgepos){
            for(auto updown: {fupidx, fdownidx}){
                int uu = _F(updown, edgepos);
                int vv = _F(updown, (edgepos+1)%3);
                int vl = std::min(uu,vv);
                int vg = std::max(uu, vv);
                TEdges_record[std::make_pair(vl, vg)] = _TEdges[updown][edgepos];
                TT_record[std::make_pair(vl, vg)] = _TT(updown, edgepos);
            }
        }
        std::pair<int, int> key0, key1, key2, key3; 
        key0 = std::make_pair(std::min(vidx, vupidx), std::max(vidx, vupidx));
        _TEdges[f0idx] = {-1, TEdges_record[key0], -1};
        _TT.row(f0idx) = Eigen::RowVector3i(f3idx, TT_record[key0], f1idx);
        
        key1 = std::make_pair(std::min(uidx, vupidx), std::max(uidx, vupidx));
        _TEdges[f1idx] = {-1, TEdges_record[key1], -1};
        _TT.row(f1idx) = Eigen::RowVector3i(f0idx, TT_record[key0], f2idx);

        key2 = std::make_pair(std::min(uidx, vdownidx), std::max(uidx, vdownidx));
        _TEdges[f2idx] = {-1, TEdges_record[key2], -1};
        _TT.row(f2idx) = Eigen::RowVector3i(f1idx, TT_record[key2], f3idx);

        key3 = std::make_pair(std::min(vidx, vdownidx), std::max(vidx, vdownidx));
        _TEdges[f3idx] = {-1, TEdges_record[key3], -1};
        _TT.row(f3idx) = Eigen::RowVector3i(f2idx, TT_record[key3], f0idx);
        



    
        // update _FL;
        _FL.conservativeResize(nF_num);
        _FL(f1idx) =_FL(f0idx); 
        _FL(f3idx) =_FL(f2idx);


         // update _F the faces with correct orientation
        nF.row(f0idx) = Eigen::RowVector3i(widx, vidx, vupidx);
        nF.row(f1idx) = Eigen::RowVector3i(widx, vupidx, uidx);
        nF.row(f2idx) = Eigen::RowVector3i(widx, uidx, vdownidx);
        nF.row(f3idx) = Eigen::RowVector3i(widx, vdownidx, vidx);
        _F = nF;
        // update _V
        Eigen::RowVector3d  wpos = (_V.row(uidx)+_V.row(vidx))/2;
        nV.row(widx) = wpos;
        _V = nV;

        // update _DblA;
        _DblA.conservativeResize(nF_num);
        _DblA(f0idx) = ((nV.row(vupidx)-wpos).cross(nV.row(vidx)-wpos)).norm();
        _DblA(f1idx) = ((nV.row(vupidx)-wpos).cross(nV.row(uidx)-wpos)).norm();
        _DblA(f2idx) = ((nV.row(vdownidx)-wpos).cross(nV.row(uidx)-wpos)).norm();
        _DblA(f3idx) = ((nV.row(vdownidx)-wpos).cross(nV.row(vidx)-wpos)).norm(); 
        
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
