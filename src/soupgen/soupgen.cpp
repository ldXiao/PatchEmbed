#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <igl/writePLY.h>
#include <igl/jet.h>
#include <igl/per_vertex_normals.h>
void to_soup(const Eigen::MatrixXd & V, 
                const Eigen::MatrixXi &F, 
                const Eigen::VectorXi & FL,
                Eigen::MatrixXd & Vsoup,
                Eigen::MatrixXi & Fsoup,
                Eigen::MatrixXd & Nsoup,
                Eigen::MatrixXd & Csoup)
{
    Eigen::MatrixXd N, Nsoup,C;
    int label_num = FL.maxCoeff()+1;
    igl::jet(FL, 0, label_num-1, C);
    igl::per_vertex_normals(V,F,N);
    Vsoup.resize(F.rows() * 3, 3);
    Fsoup.resize(F.rows(),3);
    Nsoup.resize(Vsoup.rows());
    Csoup.resize(Vsoup.rows());
    for(int fidx = 0; fidx < Fsoup.rows() ;++ fidx)
    {
        for(int j : {0,1,2})
        {
            int vj = F(fidx,j);
            Vsoup.row(fidx *3 +j) = V.row(vj);
            Nsoup.row(fidx * 3 +j) = N.row(vj);
            Csoup.row(fidx*3+j) = C.row(fidx);
            Fsoup(fidx, j) = fidx *3 +j;
        }
    }
}
int main(int argc, char *argv[]){
    std::string objfile = argv[1];
    std::string dmatfile = argv[2];
    Eigen::MatrixXd V, Vsoup;
    Eigen::MatrixXi F, Fsoup, C;
    Eigen::VectorXi FL;
    igl::readDMAT(dmatfile, FL);
    igl::readOBJ(objfile, V, F);
    Eigen::MatrixXd Vsoup, Csoup, Nsoup;
    Eigen::MatrixXi Fsoup;
    to_soup(V, F, FL, Vsoup, Fsoup, Nsoup, Csoup);
    igl::writePLY("test.ply", Vsoup, Fsoup, Nsoup, Csoup,true);
}