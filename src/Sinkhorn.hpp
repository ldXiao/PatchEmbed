//
// Created by Lind Xiao on 6/8/19.
//
#include <Eigen/Core>
#include <iostream>
#ifndef OTMAPPING_SINKHORN_H
#define OTMAPPING_SINKHORN_H
auto stable_log_sum = [](Eigen::ArrayXXd X){
    Eigen::ArrayXd max_X = X.rowwise().maxCoeff();
    for(int i=0; i< X.cols();++i){
        X.col(i)-= max_X;
    }
    Eigen::ArrayXd ret = Eigen::log(Eigen::exp(X).rowwise().sum());
    for(int i=0; i< ret.cols();++i){
        ret.col(i)+=max_X;
    }
    return ret;
};
auto row_wise_shift = [](Eigen::ArrayXXd M, Eigen::ArrayXd x){
    int m =M.cols();
    for(int i = 0; i < M.rows(); ++i){
        M.col(i)+= Eigen::ArrayXXd::Constant(m, 1,x(i)) ;
    }
    return M;
};
auto col_wise_shift = [](Eigen::ArrayXXd M, Eigen::ArrayXd x){
    int n  = M.rows();
    for(int i = 0; i < M.cols(); ++i){
        M.row(i)+= Eigen::ArrayXXd::Constant(1, n,x(i));
    }
    return M;
};

Eigen::MatrixXd sinkhorn(
        const Eigen::VectorXd & source_distr,
        const Eigen::VectorXd & target_distr,
        const Eigen::MatrixXd & K,
        double eps,
        int max_iters,
        double stop_thresh){
    Eigen::MatrixXd log_K= Eigen::log(K.array()).matrix();
    Eigen::ArrayXd u = Eigen::ArrayXd::Zero(source_distr.rows());
    Eigen::ArrayXd v = eps * Eigen::log(target_distr.array());
    Eigen::ArrayXXd C = log_K.array();
//    std::cout<<"C" << C << std::endl;
//    std::cout<<"log(a) " << Eigen::log(target_distr.array()) << std::endl;
    Eigen::ArrayXXd C_t = log_K.transpose();
    Eigen::ArrayXXd summand_u;
    Eigen::ArrayXXd summand_v;
    for(int curr_iter = 0; curr_iter< max_iters; ++curr_iter){

        Eigen::ArrayXd u_prev=u;
        Eigen::ArrayXd v_prev= v;

        summand_u = row_wise_shift(-C, v)/eps;
        u= eps* (Eigen::log(target_distr.array())-stable_log_sum(summand_u));
//        std::cout<<"u " << u << std::endl;
        summand_v = row_wise_shift(-C_t, u)/eps;
        v= eps* (Eigen::log(source_distr.array())-stable_log_sum(summand_v));
        double err_u = Eigen::abs(u_prev - u).sum();
        double err_v = Eigen::abs(v_prev - v).sum();
        std::cout << curr_iter << ", "<< err_u << std::endl;
        if (err_u < stop_thresh and err_v < stop_thresh){
            break;
        }
    }
    Eigen::ArrayXXd log_P = row_wise_shift(col_wise_shift(-log_K.array(), u),v)/eps;
    Eigen::MatrixXd  P = Eigen::exp(log_P).matrix();
    return P;
}
#endif //OTMAPPING_SINKHORN_H
