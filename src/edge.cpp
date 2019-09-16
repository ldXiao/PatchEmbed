#include "edge.h"
#include <Eigen/Core>
#include <unordered_map>
#include <iostream>
namespace bcclean{
    Edge_Compare_Result edge::compare_edge(edge b){
        if(b._label_pair == _label_pair){
            bool match_node = (b.head == head) && (b.tail == tail);
            bool anti_match_node = (b.head == tail) && (b.tail == head);
            if(match_node || anti_match_node){
                bool match_vertices = true;
                int start = (match_node) ? 0 : _edge_vertices.size()-1;
                int increment = (match_node) ? 1 : -1;
                int idx = start;
                if(b._edge_vertices.size() != _edge_vertices.size()){
                    return Edge_Compare_Result::SAME_TYPE;
                }
                for(auto p: b._edge_vertices){
                    if(p == _edge_vertices[idx]){
                        idx += increment;
                    } else{
                        match_vertices = false;
                        break;
                    }
                }
                if(match_vertices){
                    return Edge_Compare_Result::IDENTICAL;
                } else {
                    return Edge_Compare_Result::SAME_TYPE;
                }
                
            } else {
                return Edge_Compare_Result::SAME_TYPE;
            }
        } else{
            return Edge_Compare_Result::NON_MATCH;
        }
    }

    bool build_edge_dict(const Eigen::MatrixXd& V, Eigen::MatrixXi&F, const  Eigen::MatrixXi & FL, pair_map<std::pair<int,int>, std::vector<edge> > & edge_dict){
        // return lb->vector<nodes> where the nodes are unodered for each patch
        std::unordered_map<int, std::vector<int>> count_dict;
        // vertx_index -> label_list
        for(int fidx=0; fidx< F.rows(); ++fidx){
            int lb = FL(fidx,0);
            for(int i=0; i <3; ++i){
                int vi = F(fidx, i);
                auto & vect = count_dict[vi];
                auto it = std::find (vect.begin(), vect.end(), lb);
                if(it == vect.end()){
                    vect.push_back(lb);
                }
            }
        }

        


    }
}