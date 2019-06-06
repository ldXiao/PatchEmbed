//
// Created by Lind Xiao on 6/6/19.
//

#ifndef OTMAPPING_CUTMESH_H
#define OTMAPPING_CUTMESH_H


#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/StdVector>
#include <vector>
class CutMesh {
public:
    Eigen::MatrixXd Vertices;
    Eigen::MatrixXi Faces;

    Eigen::MatrixXd sampleInitial;
    Eigen::MatrixXd samplePost;
    Eigen::VectorXd sample_vals;
    int samples_nums;

    std::vector<std::unique_ptr<Eigen::MatrixXd> >  ComponentsVertices;
    std::vector<std::unique_ptr<Eigen::MatrixXi> > ComponentsFaces;
    std::vector<std::unique_ptr<Eigen::MatrixXd> > ComponentSample;
    std::vector<std::unique_ptr<Eigen::VectorXd> > ComponentVals;

    Eigen::MatrixXd TransportPlan;
    Eigen::VectorXd u;
    Eigen::VectorXd v;
    void SetInitial(const int, Eigen::MatrixXd, Eigen::MatrixXi);
    void CutWith(Eigen::Matrix3d p0,Eigen::Matrix3d n0, Eigen::Matrix3d p1, Eigen::Matrix3d n1);
    void PlotCutMesh(igl::opengl::glfw::Viewer &viewer, std::string options);

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


#endif //OTMAPPING_CUTMESH_H
