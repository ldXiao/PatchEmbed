#include "node.h"
#include <Eigen/Core>
#include <vector>
#include <map>
#include <iostream>
namespace bcclean{
        bool node::initialize(
        const int total_label_num,
        const int vidx, 
        const Eigen::MatrixXd & position, 
        const std::vector<int> labels){
            _total_label_num = total_label_num;
            _occupied_label_num = vidx;
            for(int i =0; i< _total_label_num; ++i){
                _label_occupy_dict[i]=0;
            }
            for(auto lb: labels){
                if(lb<_total_label_num && lb>-1){
                    _label_occupy_dict[lb]=1;
                }
                else{
                    std::cout << "label out of range" << std::endl;
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
    };
}