//
// Created by Lind Xiao on 6/6/19.
//

#ifndef OTMAPPING_CUTMESH_H
#define OTMAPPING_CUTMESH_H


#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/StdVector>
#include <Eigen/Sparse>
#include <vector>
#include <functional>
#include <utility>
#include <tuple>
#include <nlohmann/json.hpp>
// TODO add function to relocate the centers
// TODO add CostMatrix computation tool with centers
// TODO rename the perturb and initial blah..
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
        Eigen::MatrixXd Skeleton0;
        Eigen::MatrixXd Skeleton1;
        Eigen::VectorXi SkeletonIndices0;
        Eigen::VectorXi SkeletonIndices1;
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
        Eigen::SparseMatrix<double> WeightMatrixPerturb;
        Eigen::SparseMatrix<double> WeightMatrixInitial;
        Eigen::SparseMatrix<double> VarianceMatrixPerturb;
        Eigen::SparseMatrix<double> VarianceMatrixInitial;
        Eigen::MatrixXd CentersPerturb;
        Eigen::MatrixXd CentersInitial;
        Eigen::VectorXi SampleSourceFace;
        Eigen::SparseMatrix<double> ElasticTensor;
        Eigen::SparseMatrix<double> QuadraticCostMatrix;
        void plot_CutMesh(igl::opengl::glfw::Viewer &viewer, unsigned char options);
        void set_initial(const Eigen::MatrixXd &, const Eigen::MatrixXi &, const int,
                         std::function<double(Eigen::Vector3d)>);
        void set_sinkhorn_const(const double, const double, const int);
        void build_graph(double range, double l);
        void build_tree();
        void set_initial_from_json(const nlohmann::json &);
        void separate_cube_faces();

        void compute_WeightMatrix(double sigma);
        void locate_Centers(
            const Eigen::SparseMatrix<double> & weightmatrix,
            const Eigen::MatrixXd & transportplan, 
            const Eigen::MatrixXd & sample,
            Eigen::MatrixXd & centers, 
            char options);



        void cut_with(const Eigen::Vector3d & p0,
                const Eigen::Vector3d & n0,
                const Eigen::Vector3d & p1,
                const Eigen::Vector3d & n1);

        void perturb(const int seed, double max_shift,  double min_shift);

        Eigen::MatrixXd to_nearest();

        void compute_CostMatrix(const Eigen::MatrixXd &, const Eigen::MatrixXd &,char options);
        void Sinkhorn();
        void VarSinkhorn();
        double SinkhornEps=0;
        double SinkhornThreshold=0;
        int SinkhornMaxIter=0;
        bool round = true;
    };
}

// helper funciton declaritions
double color_error(Eigen::MatrixXd C0, Eigen::MatrixXd C1);
void weight_matrix(
        double sigma,
        const Eigen::MatrixXd  & sample,
        const Eigen::MatrixXi & faces,
        const Eigen::MatrixXi & samplesource_indice,
        Eigen::SparseMatrix<double> & wmtx);
#endif //OTMAPPING_CUTMESH_H
