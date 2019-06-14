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
#include <nlohmann/json.hpp>
namespace OTMapping {
    using MeshPair = std::tuple<Eigen::MatrixXd, Eigen::MatrixXi>;
    using SampleTriplet = std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>;

    class CutMesh {
    public:
        Eigen::MatrixXd Vertices;
        Eigen::MatrixXi Faces;

        Eigen::MatrixXd SampleInitial;
        Eigen::VectorXd SampleVals;
        Eigen::MatrixXd SamplePerturb;
        int SampleNum;
        int point_size;

        std::vector<std::unique_ptr<Eigen::MatrixXd> > ComponentsVertices;
        std::vector<std::unique_ptr<Eigen::RowVector3d> > Shifts;
        std::vector<std::unique_ptr<Eigen::MatrixXi> > ComponentsFaces;
        std::vector<std::unique_ptr<Eigen::MatrixXd> > ComponentsSample;
        std::vector<std::unique_ptr<Eigen::VectorXi> > ComponentsSampleIndices;
        // store the indices of samples in original SampleInitial


        

        Eigen::MatrixXd TransportPlan;
        Eigen::MatrixXd CostMatrix;
        Eigen::VectorXi SampleSourceFace;
        void plot_CutMesh(igl::opengl::glfw::Viewer &viewer, unsigned char options);

        void set_initial(const Eigen::MatrixXd &, const Eigen::MatrixXi &, const int,
                         std::function<double(Eigen::Vector3d)>);
        void set_sinkhorn_const(const double, const double, const int);

        void set_initial_from_json(const nlohmann::json &);
        void separate_cube_faces();

        void cut_with(const Eigen::Vector3d & p0,
                const Eigen::Vector3d & n0,
                const Eigen::Vector3d & p1,
                const Eigen::Vector3d & n1);

        void perturb(const int seed, double max_shift,  double min_shift);

        Eigen::MatrixXd to_nearest();

        void compute_CostMatrix(Eigen::MatrixXd, Eigen::MatrixXd,char options);

        void Sinkhorn();
        double SinkhornEps=0;
        double SinkhornThreshold=0;
        int SinkhornMaxIter=0;
        bool round = true;
    };
}
double color_error(Eigen::MatrixXd C0, Eigen::MatrixXd C1);

#endif //OTMAPPING_CUTMESH_H
