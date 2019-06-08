//
// Created by Lind Xiao on 6/6/19.
//
#include "CutMesh.h"
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <math.h>
#include <memory>
#include <utility>
// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;


using namespace std;
using namespace OTMapping;
class box{
public:
    double width;
    double height;
    double bridth;
    std::vector<int> chain;
    void setVal(double a, double b , double c){width=a; height=b; bridth =c;}
    double getVolum(){
        return width * bridth * bridth;
    }
};
int main(int argc, char *argv[]){
    using namespace std;
    using namespace Eigen;
    if (argc != 2) {
        cout << "Usage cutmeshtest_bin mesh.obj" << endl;
        exit(0);
    }
    std::cout << "holy shit0" << std::endl;
    // Load a mesh in OBJ format
    igl::readOBJ(argv[1], V, F);
    igl::opengl::glfw::Viewer viewer;
    CutMesh CM;
    auto func = [](Eigen::Vector3d x)->double{return std::sin(x[0]+x[1]+x[2]);};
    CM.set_initial(V, F, 2000, func);
//    CM.plot_CutMesh(viewer,"aha");
    box aboox;
    std::cout << "holy shit1" << std::endl;
    aboox.setVal(1,2,3);
//    std::unique_ptr<box> ab;
//    Eigen::MatrixXd* Mxptr;
//    *Mxptr=Eigen::MatrixXd::Zero(5,5);
//    Mxptr->resize(V.rows(),V.cols());
//    std::cout << Mxptr->row(0)<<std::endl;
    Eigen::MatrixXd Z =Eigen::MatrixXd::Zero(5,5);
    std::unique_ptr<Eigen::MatrixXd > Y(new Eigen::MatrixXd());
    *Y = V;
    Eigen::Vector3d p0(0,0,0);
    Eigen::Vector3d p1(0,0,0);
    Eigen::Vector3d n1(0,0,1);
    Eigen::Vector3d n0(0,1,0.5);
    std::cout << "holy shit3" << std::endl;
    CM.cut_with(p0,n0,p1,n1);
    CM.perturb(3, 0.1, 0.02);
    CM.plot_CutMesh(viewer,'2');
    std::cout << "holy shit4" << std::endl;
    viewer.launch();
}
