#include "mapping_patch.h"
#include <map>
#include <vector>
#include <algorithm>
#include <Eigen/Core>
#include <igl/remove_unreferenced.h>
#include <igl/boundary_loop.h>
#include <igl/boundary_facets.h>
#include <igl/list_to_matrix.h>
#include <igl/vertex_triangle_adjacency.h>
#include <iostream>

const double PI = 3.14159265358979323846;

namespace bcclean{
    void build_patch_dict(const Eigen::MatrixXi &FL, std::map<int, std::vector<int> > & patch_dict){
        patch_dict.clear();
        for(int fidx =0 ; fidx < FL.rows(); fidx++){
            int patch_idx = FL(fidx,0);
            auto it = patch_dict.find(patch_idx);
            if(it == patch_dict.end()){
                std::vector<int> chunk;
                chunk.push_back(fidx);
                patch_dict[patch_idx]= chunk;
            }
            else{
                patch_dict[patch_idx].push_back(fidx);
            }
        }
    }

    bool node::initialize(
        const int total_label_num, 
        const Eigen::MatrixXd & position, 
        const std::vector<int> labels){
            _total_label_num = total_label_num;
            _occupied_label_num = 0;
            for(int i =0; i< _total_label_num; ++i){
                _label_occupy_dict[i]=0;
            }
            for(auto lb: labels){
                if(lb<_total_label_num && lb>-1){
                    _label_occupy_dict[lb]=1;
                }
                else{
                    std::cout<< "label out of range"<<std::endl;
                    return false;
                }
            }
            if(position.rows() == 1 and position.cols()==3){
                    this->_position = position;
                    return true;
            }
            else{
                    std::cout<< "position should be initialized as rowvector3d" <<std::endl;
                    return false;
            }
    };
    
    bool node::of_same_type(const node & b){
            if(_total_label_num != b._total_label_num){
                return false;
            }
            if(_occupied_label_num != b._occupied_label_num){
                return false;
            }
            for(int i =0; i< _total_label_num; ++i){
               if( _label_occupy_dict[i]!=b._label_occupy_dict.at(i)){
                   return false;
               }
            }
            return true;
    };
    bool node::at_same_position(const Eigen::MatrixXd& position){
        double tol = 10e-6;
        if ((position.row(0)-_position.row(0)).norm()< tol){
            return true;
        }
        else{
            return false;
        }
    }

    std::vector<std::vector<node>> build_label_nodes_list(
        const Eigen::MatrixXd &V, 
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXi &FL){
            // return lb->vector<nodes> where the nodes are unodered for each patch
            std::map<int, std::vector<int>> count_dict;
            // vertx_index -> label_list
            for(int fidx=0; fidx< F.rows(); ++fidx){
                int lb = FL(fidx,0);
                for(int i=0; i <3; ++i){
                    int vi = F(fidx, i);
                    auto & vect = count_dict[vi];
                    auto it = std::find (vect.begin(), vect.end(), lb);
                    if(it == vect.end()){
                        std::cout<<"push0"<<vi<<","<<lb<<std::endl;
                        vect.push_back(lb);
                        std::cout<<vect.size()<<std::endl;
                    }
                }
            }
            int total_label_num = FL.maxCoeff()+1;
            std::vector<std::vector<node>> result(total_label_num);
            for(auto & [vi, vect]: count_dict){
                if(vect.size()>2){
                    Eigen::RowVector3d position = V.row(vi);
                    node nd;
                    nd.initialize(total_label_num, position, vect);
                    for(int lb: vect){
                        result[lb].push_back(nd);
                        std::cout<< "push"<<std::endl;
                    }
                }
            }
            return result;
        };

        void extract_label_patch_mesh(
            const Eigen::MatrixXd& V, 
            const Eigen::MatrixXi& F, 
            const Eigen::MatrixXi&FL, 
            const int lb_in, 
            Eigen::MatrixXd& V_i, 
            Eigen::MatrixXi& F_i){
                Eigen::MatrixXi F_l = Eigen::MatrixXi::Constant(F.rows(),3, 0);
                int count = 0;
                for(int fidx=0; fidx< FL.rows(); ++fidx){
                    int lb = FL(fidx,0);
                    if(lb==lb_in){
                        F_l.row(count)=F.row(fidx);
                        count += 1;
                    }
                }
                F_l.conservativeResize(count,3);
                {
                    Eigen::MatrixXi I;
                    igl::remove_unreferenced(V, F_l, V_i, F_i, I);
                }
        }

        void map_vertices_to_regular_polygon(
            const Eigen::MatrixXd &V, 
            const Eigen::MatrixXi & F, 
            std::vector<node> & nodes, 
            Eigen::VectorXi & bnd,
            Eigen::MatrixXd & bnd_uv,
            std::vector<node>& ordered_nodes){
                ordered_nodes.clear();
                igl::boundary_loop(F, bnd);
                bnd_uv = Eigen::MatrixXd::Constant(bnd.rows(),3,0);
                std::vector<int> nails;
                std::map<int, int> loop_patch_dict;
                int count = -1;
                for(int bnd_idx=0; bnd_idx<bnd.rows(); ++bnd_idx){
                    for(node & nd:nodes){
                        if(nd.at_same_position(V.row(bnd(bnd_idx)))){
                            ordered_nodes.push_back(nd);
                            nails.push_back(bnd_idx);
                            count +=1;
                        }
                    }
                    loop_patch_dict[bnd_idx]=count;
                }

                int node_num = ordered_nodes.size();
                for(auto item: loop_patch_dict){
                    if(item.second<0) item.second=node_num-1;
                    // fix the patch label
                }

                // arc length parametrization for each patch
                std::map<int, double> patch_arc_dict;
                double sum = 0;
                int nail_count = 0;
                int start = nails[0];
                for(int shift =0 ; shift< bnd.rows();++shift){
                    int curr_idx_raw =(start+shift);
                    int curr_idx = (start+shift) % bnd.rows();
                    int next_idx = (start+shift +1) % bnd.rows();
                    double len = (V.row(bnd(next_idx,0))-V.row(bnd(curr_idx,0))).norm();
                    patch_arc_dict[curr_idx]= sum;
                    if(loop_patch_dict[curr_idx]==loop_patch_dict[next_idx]){
                        sum += len;
                    }
                    if(loop_patch_dict[curr_idx]!=loop_patch_dict[next_idx]){
                        sum += len;
                        // get the full arc length of a patch
                        int curr_nail = nails[nail_count];
                        std::cout<< curr_nail <<"in"<<nails.size()<< std::endl;
                        std::cout<< bnd.rows() <<"boundary loops"<< std::endl;
                        for(int shuttle= curr_nail; shuttle< curr_idx_raw+1; ++shuttle){
                            double ratio = patch_arc_dict[shuttle%bnd.rows()]/sum;
                        
                            std::cout<< shuttle%bnd.rows()<<"  for"<<patch_arc_dict[shuttle%bnd.rows()] <<"/"<< sum <<"="<<ratio<<std::endl;
                            patch_arc_dict[shuttle%bnd.rows()]=ratio;
                        }
                        // quotient by the full arc length to get ratio
                        nail_count+=1;
                        sum = 0;
                    }
                }
                std::cout<< patch_arc_dict.size() <<"total"<< std::endl;
                for(auto item: patch_arc_dict){
                    std::cout << item.first <<".."<<item.second <<std::endl;
                }
                start = nails[0];
                nail_count = 0;
                for(int shift =0 ; shift< bnd.rows();++shift){
                    int curr_idx = (start+shift) % bnd.rows();
                    int next_idx = (start+shift+1) % bnd.rows();
                    double x, y;
                    double ratio = patch_arc_dict[curr_idx];
                    double curr_theta = nail_count *2 * PI/ node_num;
                    double next_theta = (nail_count+1)*2 * PI/node_num;
                    x = (1-ratio)* std::cos(curr_theta) + ratio * std::cos(next_theta);
                    y = (1-ratio) * std::sin(curr_theta) + ratio * std::sin(next_theta);
                    bnd_uv.row(curr_idx) = Eigen::RowVector3d(x,y,0);
                    if(loop_patch_dict[curr_idx]!=loop_patch_dict[next_idx]){
                        nail_count+=1;
                    }
                }
        }
    
    bool no_consecutive_on_loop(int loop_len, std::vector<int> peaks){
        if(peaks.size()<2){
            return true;
        }
        bool res = true;
        std::sort(peaks.begin(), peaks.end());
        
        for(int i =0; i<peaks.size();++i){
            int curr = peaks[i];
            int next = peaks[(i+1)%peaks.size()];
            if((next-curr)==1){
                res = false;
                break;
            }
            if((next+loop_len-curr)==1){
                res = false;
                break;
            }
        }
        return res;
    }

    bool split_ears(
        const Eigen::MatrixXd & V, 
        const Eigen::MatrixXi & F, 
        const Eigen::VectorXi & bnd,
        const std::vector<int> & nails,
        Eigen::MatrixXd &NV,
        Eigen::MatrixXi &NF,
        Eigen::VectorXi &I,
        Eigen::VectorXi &J){
            // TODO improve performace use 
            std::vector<int> degenerate_faces;
            std::vector<int> degenerate_boundary_peaks;
            for(int fidx=0; fidx< F.rows(); ++fidx){
                std::vector<int> count_list;
                for(int vinf; vinf <3; ++vinf){
                    for(int shuttle=0; shuttle < bnd.rows();++shuttle){
                        if(bnd[shuttle]==F(fidx,vinf)){
                            count_list.push_back(bnd[shuttle]);
                            break;
                        }
                    }
                }
                if(count_list.size()>2){
                    // degenerate face detected
                    int ear_peak = count_list[1];
                    auto it = std::find(nails.begin(), nails.end(), ear_peak);
                    if(it==nails.end()){
                        degenerate_faces.push_back(fidx);
                        // the count_list already store the order of triangle points in bnd push the middle one into the degenerate_boundary_peaks if it is not node
                        degenerate_boundary_peaks.push_back(ear_peak);
                    }
                }
            }
            if(degenerate_faces.size()==0){
                NV = V;
                NF = F;
                std::cout<< "no degenerate face detected"<<std::endl;
                return true;
            }
            else{
                // TODO add functions to deal with ear spliting
                return false;
            }
        }
    
    
    bool build_nails(const Eigen::MatrixXd & V, const Eigen::VectorXi & bnd, std::vector<node> & nodes, std::vector<int> & nails, std::map<int, node> & nails_node_dict){
        int count = 0;
        nails.clear();
        nails_node_dict.clear();
        for(int bnd_idx=0; bnd_idx<bnd.rows(); ++bnd_idx){
            for(node & nd:nodes){
                if(nd.at_same_position(V.row(bnd(bnd_idx)))){
                    nails.push_back(bnd_idx);       
                    nails_node_dict[bnd_idx]=nd;
                }
            }
        }
        if(nails.size()<3){
            std::cout << "too few nodes, unable to map"<<std::endl;
            return false;
        }
        return true;
    }
    
    bool cyc_flip_mapping(
        std::vector<node> & nodes, 
        std::vector<node> & target_nodes, 
        std::map<int,int> & mapping){
            /* tested */
            /* would find the permutation between nodes and target_nodes 
            which are product of cyclic permutation and reversing*/
            mapping.clear();
            if(nodes.size()!=target_nodes.size()){
                return false;
            }
            const int size = nodes.size();
            for(int shift=0; shift < size; ++shift){
                std::vector<int> match_list;
                for(int node_idx=0; node_idx<size; ++node_idx){
                    int target_node_idx = (node_idx+shift) % size;
                    if(nodes[node_idx].of_same_type(target_nodes[target_node_idx])){
                        match_list.push_back(target_node_idx);
                    } else{
                        break;
                    }
                }
                if(match_list.size()==size){
                    for(int i = 0; i< size; ++i){
                        // initialize positive cyc
                        mapping[i]=match_list[i];
                    }
                    return true;
                }
            }
            for(int shift=0; shift < size; ++shift){
                std::vector<int> match_list;
                for(int node_idx=0; node_idx<size; ++node_idx){
                    int target_node_idx = (size-node_idx+shift) % size;
                    // negative loop
                    if(nodes[node_idx].of_same_type(target_nodes[target_node_idx])){
                        match_list.push_back(target_node_idx);
                    } else{
                        break;
                    }
                }
                if(match_list.size()==size){
                    for(int i = 0; i< size; ++i){
                        mapping[i]=match_list[i];
                    }
                    return true;
                }
            }
            return false;
    }

    bool reordering_arcs(
        const Eigen::VectorXi & bnd, 
        const std::vector<int> & nails, 
        std::map<int, node> nails_node_dict, 
        std::vector<node> & target_nodes,
        std::vector<int> & ccw_ordered_nails){
            ccw_ordered_nails.clear();
            std::vector<node> nodes;
            for(auto nl:nails){
                nodes.push_back(nails_node_dict[nl]);
            }
            std::map<int,int> mapping;
            // cycl or flip mapping betweins nails if any
            if(cyc_flip_mapping(nodes, target_nodes, mapping)){
                std::cout << mapping.size()<<"mapping size"<<std::endl;
                ccw_ordered_nails.resize(mapping.size());
                for(int i = 0; i< mapping.size();++i){
                    ccw_ordered_nails[mapping[i]]=nails[i];
                }
                return true;
            }
            return false;
    }

    int _loop_next(int total, int curr, bool reverse){
        if(curr>=total or curr < 0 or total < 1){
            return -1;
        }
        if(reverse){
            if (curr==0){
                return total -1;
            }
            else{
                return curr - 1;
            }
        }
        else{
            return (curr +1) % total;
        }
    }
    void set_edge_arc_ratio_list(
        const Eigen::MatrixXd & V,
        const Eigen::VectorXi & bnd,
        const std::vector<int> & nails,
        std::vector<int> & ccw_ordered_nails,
        std::map<int, double> & edge_arc_ratio_list){
            edge_arc_ratio_list.clear();
            int first = ccw_ordered_nails[0];
            int second = ccw_ordered_nails[1];
            bool reverse= false;
            double sum = 0;
            if(second < first){reverse = true;}
            int curr_nail_idx = 0;
            int curr_idx = ccw_ordered_nails[curr_nail_idx];
            std::vector<int> arc_idices;
            for(int count =0 ; count < bnd.rows(); ++count){
                int next_nail_idx = _loop_next(ccw_ordered_nails.size(), curr_nail_idx, false);
                int next_nail = ccw_ordered_nails[next_nail_idx];
                int next_idx = _loop_next(bnd.rows(), curr_idx, reverse);
                double len = (V.row(bnd(next_idx,0))-V.row(bnd(curr_idx,0))).norm();
                edge_arc_ratio_list[curr_idx]=sum;
                arc_idices.push_back(curr_idx);
                curr_idx = next_idx;
                if(next_idx!= next_nail){
                    sum += len;
                }
                else{
                    sum += len;
                    curr_nail_idx = next_nail_idx;
                    for(auto it = arc_idices.begin(); it != arc_idices.end(); ++it){
                        int index = *it;
                        edge_arc_ratio_list[index]/= sum;
                        //quotient
                    }
                    arc_idices.clear();
                    sum = 0;
                }
            }
            
        }

    bool set_bnd_uv(
        const Eigen::VectorXi & bnd, 
        const std::vector<int> & ccw_ordered_nails,
        const std::map<int, double> & edge_arc_ratio_list,
        Eigen::MatrixXd & bnd_uv){
            int first = ccw_ordered_nails[0];
            int second = ccw_ordered_nails[1];
            bool reverse= false;
            if(second < first){reverse = true;}
            int start = ccw_ordered_nails[0];
            int curr_nail_idx = 0;
            int curr_idx = ccw_ordered_nails[curr_nail_idx];
            bnd_uv.resize(bnd.rows(), 3);
            for(int count =0 ; count < bnd.rows(); ++count){
                int next_nail_idx = _loop_next(ccw_ordered_nails.size(), curr_nail_idx, false);
                int next_nail = ccw_ordered_nails[next_nail_idx];
                int next_idx = _loop_next(bnd.rows(), curr_idx, reverse);
                double x, y;
                double ratio = edge_arc_ratio_list.at(curr_idx);
                double curr_theta = curr_nail_idx *2 * PI/ ccw_ordered_nails.size();
                double next_theta = next_nail_idx *2 * PI/ ccw_ordered_nails.size();
                x = (1-ratio)* std::cos(curr_theta) + ratio * std::cos(next_theta);
                y = (1-ratio) * std::sin(curr_theta) + ratio * std::sin(next_theta);
                bnd_uv.row(curr_idx) = Eigen::RowVector3d(x,y,0);
                curr_idx = next_idx;
                if(next_idx!= next_nail){
                    continue;
                }
                else{
                    curr_nail_idx = next_nail_idx;
                }
            }

        }

    bool mapping_patch::build_patch(
        const Eigen::MatrixXd &Vi, 
        const Eigen::MatrixXi & Fi, 
        std::vector<node> & nodes, int lb_in){
        std::cout << "1"<<std::endl;
        _V_raw = Vi;
        _F_raw = Fi;
        std::vector<std::vector<int> > L;
        igl::boundary_loop(_F_raw, L);
        std::cout << "2"<<std::endl;
        if(L.size()!=1){
            // input mesh is not isomorphic to disk
            return false;
        }
        igl::list_to_matrix(L[0], _bnd_raw);
        std::cout << "3"<<std::endl;
        if(! build_nails(_V_raw, _bnd_raw, nodes, _nails, _nails_nodes_dict)){
            std::cout << "not a canonical patch"<<std::endl;
            return false;
        }
        {
            Eigen::VectorXi I;
            Eigen::VectorXi J;
            if(!split_ears(_V_raw, _F_raw, _bnd_raw, _nails, _V_ndg, _F_ndg, I, J)){
                std::cout << "ears_detected"<<std::endl;
                return false;
            }
             _bnd_ndg = _bnd_raw;
             std::vector<node> target_nodes;
             for(auto nl:_nails){
                 target_nodes.push_back(_nails_nodes_dict[nl]);
             }
             reordering_arcs(_bnd_ndg, _nails, _nails_nodes_dict, target_nodes, _ccw_ordered_nails);
             // initilize the arc numbeirng with itself.
        }

        set_edge_arc_ratio_list(_V_ndg, _bnd_ndg, _nails, _ccw_ordered_nails, _edge_arc_ratio_list);
        set_bnd_uv(_bnd_ndg, _ccw_ordered_nails, _edge_arc_ratio_list, _bnd_uv);
        
    }
}