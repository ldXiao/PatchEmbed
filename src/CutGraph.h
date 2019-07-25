//
// Created by Lind Xiao on 7/25/19.
//

#ifndef OTMAPPING_CUTGRAPH_H
#define OTMAPPING_CUTGRAPH_H

#include <vector>
#include <Eigen/Core>
#include <tuple>
class CutGraph {
public:
    std::vector<std::tuple<double, double, double > > _vertices;
    std::vector<std::pair<std::size_t, std::size_t> > _edges;
    std::vector<double> EdgeWeights;
    Eigen::MatrixXd Vertices;
    Eigen::MatrixXi Edges;
    void _set_Vertices(Eigen::MatrixXd);
    void _set_Edges_from_KNN(int k);
    void _compute_EdgeWeights();
};


#endif //OTMAPPING_CUTGRAPH_H
