//
// Created by Lind Xiao on 6/6/19.
//

#ifndef OTMAPPING_CUTMESH_H
#define OTMAPPING_CUTMESH_H


#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/StdVector>
#include <vector>
#include <functional>
#include <utility>
#include <tuple>
namespace OTMapping {
    using MeshPair = std::tuple<Eigen::MatrixXd, Eigen::MatrixXi>;
    using SampleTriplet = std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>;

    class CutMesh {
    public:
        Eigen::MatrixXd Vertices;
        Eigen::MatrixXi Faces;

        Eigen::MatrixXd SampleInitial;
        Eigen::MatrixXd SamplePost;
        Eigen::VectorXd SampleVals;
        Eigen::MatrixXd SampleNorms;
        int SampleNum;

        std::vector<std::unique_ptr<Eigen::MatrixXd> > ComponentsVertices;
        std::vector<std::unique_ptr<Eigen::MatrixXi> > ComponentsFaces;
        std::vector<std::unique_ptr<Eigen::MatrixXd> > ComponentSample;
        std::vector<std::unique_ptr<Eigen::VectorXd> > ComponentVals;
        std::vector<std::unique_ptr<Eigen::MatrixXd> > ComponentNorms;

        Eigen::MatrixXd TransportPlan;
        Eigen::MatrixXd CostMatrix;
        Eigen::VectorXd u;
        Eigen::VectorXd v;

        void plot_CutMesh(igl::opengl::glfw::Viewer &viewer, unsigned char options);

        void set_initial(const Eigen::MatrixXd &, const Eigen::MatrixXi &, const int,
                         std::function<double(Eigen::Vector3d)>);

        void cut_with(const Eigen::Vector3d & p0,
                const Eigen::Vector3d & n0,
                const Eigen::Vector3d & p1,
                const Eigen::Vector3d & n1);
        void perturb(const int seed, double max_shift,  double min_shift);

    private:
        Eigen::MatrixXd _sample0;
        Eigen::MatrixXd _sample_vals0;
        Eigen::MatrixXd _sample1;
        Eigen::MatrixXd _sample_vals1;
        Eigen::MatrixXd _sample2;
        Eigen::MatrixXd _sample_vals2;
        Eigen::MatrixXd _sample3;
        Eigen::MatrixXd _sample_vals3;

        Eigen::MatrixXd _vertices0;
        Eigen::MatrixXi _faces0;
        Eigen::MatrixXd _vertices1;
        Eigen::MatrixXi _faces1;
        Eigen::MatrixXd _vertices2;
        Eigen::MatrixXi _faces2;
        Eigen::MatrixXd _vertices3;
        Eigen::MatrixXi _faces3;

    };
}

#endif //OTMAPPING_CUTMESH_H
