//
// Created by Lind Xiao on 6/8/19.
//
#include "Sinkhorn.hpp"
#include <Eigen/Core>
#include <iostream>
#include <nlohmann/json.hpp>
#include <vector>
#include "ExpressionValue.hpp"
int main(){
//    auto stable_log_sum = [](Eigen::ArrayXXd X){
//        double max_X = X.maxCoeff();
//        X = X - Eigen::ArrayXXd::Constant(X.rows(),X.cols(),max_X);
//
//        Eigen::ArrayXd ret = Eigen::log(Eigen::exp(X).rowwise().sum())+Eigen::ArrayXd::Constant(X.rows(), max_X);
//        return ret;
//    };
//    auto row_wise_shift = [](Eigen::ArrayXXd M, Eigen::ArrayXd x){
//        int m =M.cols();
//        for(int i = 0; i < M.rows(); ++i){
//            M.row(i)+= Eigen::ArrayXXd::Constant(1, m,x(i)) ;
//        }
//        return M;
//    };
    std::cout << "everything OK"<<std::endl;
    Eigen::MatrixXd K(3,3);
    K <<1,2,3,
        4,5,6,
        7,8,9;

    ExpressionValue jasd;
    std::string aasd = "sin(1)";
    jasd.init(aasd);
//    K << -1.09861, -694.246, -1099.71,
//            -1387.39, -1610.54, -1792.86,
//            -1947.01, -2080.54 -2198.32;
    Eigen::MatrixXd Y =sinkhorn(
            Eigen::VectorXd::Constant(3,1/double(3)),
            Eigen::VectorXd::Constant(3,1/double(3)),
            K, 1e-6, 10, 1e-5);
    std::cout << Y<<std::endl;
    std::vector<int> vec;
    for(int i=0; i< 10;++i){
        vec.push_back(i);
    }
    vec.insert(vec.end(), 123);
    for(auto it = vec.begin(); it!=vec.end();++it){
        std::cout << *it << std::endl;
    }
}

