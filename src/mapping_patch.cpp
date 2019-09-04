#include "mapping_patch.h"
#include <map>
#include <vector>
#include <algorithm>
#include <Eigen/Core>
#include <igl/remove_unreferenced.h>
#include <igl/boundary_loop.h>
#include <igl/boundary_facets.h>
#define PI 3.141592653

namespace bcclean{
    void build_patch_dict(const Eigen::MatrixXi &FL, std::map<int, std::vector<int> > & patch_dict){
        patch_dict.clear();
        for(int fidx =0 ; fidx < FL.rows(); fidx++){
            int patch_idx = FL(fidx,0);
//            std::vector<int> chunk;
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
                std::sort(vect.begin(), vect.end());
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
}