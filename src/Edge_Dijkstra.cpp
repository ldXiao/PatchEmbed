#include "Edge_Dijkstra.h"
#include "polyline_distance.h"
#include "kdtree_NN_Eigen.hpp"
#include "helper.h"
#include <utility>
#include <vector>
#include <igl/edges.h>
#include <queue>
#include <limits>
namespace bcclean {
namespace Trace{
    void setWeight(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        const Eigen::MatrixXd & V_bad,
        const Eigen::MatrixXi & F_bad,
        const edge & edg,
        Eigen::SparseMatrix<double> & Weights
    )
    {
        std::vector<Eigen::Triplet<double>> ww;
        Eigen::MatrixXi Edges;
        igl::edges(F, Edges);
        Weights.resize(V.rows(), V.rows());
        int upratio = 10;
        int power=1;
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
        for(int eidx = 0 ; eidx < Edges.rows(); eidx++)
        {
            int x= Edges(eidx, 0);
            int y = Edges(eidx,1);
            Eigen::RowVector3d query = (V.row(x) + V.row(y)) * 0.5 ;
            int sample_idx= kd_tree_NN_Eigen(sample_kdt, query);
            double dis = pow((query - sample.row(sample_idx)).norm(), power);
            ww.emplace_back(std::min(x,y), std::max(x,y),dis);
        }
        Weights.setFromTriplets(ww.begin(), ww.end());
        return;
    }

    void Edge_Dijkstra(
        const std::vector<std::vector<int> > & VV,
        const int source,
        const int target,
        const Eigen::SparseMatrix<double> & Weights,
        std::vector<int> & path
    ){
        std::vector<double> dist(VV.size());
        std::vector<int> prev(VV.size());
        std::fill(dist.begin(), dist.end(), std::numeric_limits<double>::max());
        std::fill(prev.begin(), prev.end(), -1);
        dist[source] = 0;
        std::priority_queue<std::pair<double, int>, std::vector<std::pair<double,int> >, comparator > PQ;
        PQ.emplace(std::make_pair(0,source));
        while(!PQ.empty())
        {
            std::pair<double, int> item = PQ.top();
            PQ.pop();
            int u = item.second;
            for(auto v: VV.at(u)){
                double alt = dist[u] + Weights.coeff(std::min(u,v), std::max(u,v));
                if(alt < dist[v])
                {
                    dist[v] = alt;
                    prev[v] = u;
                    PQ.push(std::make_pair(alt,v));
                }
            }
        }
        int pp = target;
        // assert(pp!= -1);
        path.clear();
        while(pp!= -1)
        {
            path.push_back(pp);
            pp = prev[pp];
        }
        // assert(path.size()>=2);
        // assert(path[path.size()-1]== source);
        return;
    }

    void setWeight1(
        const std::vector<Eigen::RowVector3d> & V,
        const std::vector<Eigen::RowVector3i> & F,
        const std::vector<Eigen::RowVector3d> & V_bad,
        const edge & edg,
        Eigen::SparseMatrix<double> & Weights
    )
    {
        std::vector<Eigen::Triplet<double>> ww;
        Eigen::MatrixXi Edges;
        Eigen::MatrixXi F_matrix;
        Helper::to_matrix(F,F_matrix);
        igl::edges(F_matrix, Edges);
        Weights.resize(V.size(), V.size());
        int upratio = 10;
        int power=1;
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
                sample.row(jj) = V_bad[Elist[batch]];
            }
            else
            {
                double lambda = (double)(remain) / (double)(upratio);
                sample.row(jj) = (1-lambda)* V_bad[Elist[batch]] + lambda* V_bad[Elist[batch+1]];
            }
        }
        kd_tree_Eigen<double> sample_kdt(sample.cols(),std::cref(sample),10);
        sample_kdt.index->buildIndex();
        for(int eidx = 0 ; eidx < Edges.rows(); eidx++)
        {
            int x= Edges(eidx, 0);
            int y = Edges(eidx,1);
            Eigen::RowVector3d query = (V.at(x) + V.at(y)) * 0.5 ;
            int sample_idx= kd_tree_NN_Eigen(sample_kdt, query);
            double dis = pow((query - sample.row(sample_idx)).norm(), power);
            ww.emplace_back(std::min(x,y), std::max(x,y),dis);
        }
        Weights.setFromTriplets(ww.begin(), ww.end());
        return;
    }
}
}