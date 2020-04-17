#include <Eigen/Dense>
#include <iostream>
#include <string>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <igl/jet.h>
#include <igl/per_vertex_normals.h>
void writePlyColor(std::string name,const Eigen::MatrixXd & Vsoup, const Eigen::MatrixXi & Fsoup, const Eigen::MatrixXd & Nsoup, const Eigen::MatrixXi Csoup)
{
    std::ofstream file;
    file.open(name);
    file << "ply\n";
    file << "format ascii 1.0\n";
    file << "element vertex " <<Vsoup.rows() << "\n";
    file << "property double x\n";
    file <<"property double y\n";
    file << "property double z\n";
    file <<"property double nx\n";
    file <<"property double ny\n";
    file <<"property double nz\n";
    file <<"property uchar red\n";
    file <<"property uchar green\n";
    file <<"property uchar blue\n";
    file <<"element face " << Fsoup.rows() << "\n";
    file <<"property list uchar int vertex_indices\n";
    file <<"end_header\n";
    for(int vidx =0 ; vidx < Vsoup.rows(); ++vidx)
    {
        file << Vsoup(vidx,0) <<" ";
        file << Vsoup(vidx,1) <<" ";
        file << Vsoup(vidx,2) <<" ";
        file << Nsoup(vidx,0) <<" ";
        file << Nsoup(vidx,1) <<" ";
        file << Nsoup(vidx,2) <<" ";
        file << Csoup(vidx,0) <<" ";
        file << Csoup(vidx,1) <<" ";
        file << Csoup(vidx,2) <<"\n";
    }
    for(int fidx = 0; fidx < Fsoup.rows(); ++fidx)
    {
        file << 3 << " ";
        file << Fsoup(fidx,0) << " ";
        file << Fsoup(fidx,1) << " ";
        file << Fsoup(fidx,2) << "\n";
    }
}
void to_soup(const Eigen::MatrixXd & V, 
                const Eigen::MatrixXi &F, 
                const Eigen::VectorXi & FL,
                Eigen::MatrixXd & Vsoup,
                Eigen::MatrixXi & Fsoup,
                Eigen::MatrixXd & Nsoup,
                Eigen::MatrixXd & Csoup)
{
    Eigen::MatrixXd N,C;
    int label_num = FL.maxCoeff()+1;
    igl::jet(FL, 0, label_num-1, C);
    if(FL.minCoeff()==-1)
    {
        Eigen::VectorXi  NFL = FL.unaryExpr([](int vl){return (vl>-1)? 1: -1;});
        for(int fidx =0; fidx < FL.rows(); ++fidx){
            if(NFL(fidx)>0)
            {
                C.row(fidx)=Eigen::RowVector3d(0.6,0.6,0.95);
            }
            else{
                C.row(fidx)=Eigen::RowVector3d(0.95,0.3,0.23);
            }
        }

    }
    igl::per_vertex_normals(V,F,N);
    Vsoup.resize(F.rows() * 3, 3);
    Fsoup.resize(F.rows(),3);
    Nsoup.resize(Vsoup.rows(),3);
    Csoup.resize(Vsoup.rows(),3);
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
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXi FL;
    igl::readDMAT(dmatfile, FL);
    igl::readOBJ(objfile, V, F);
    Eigen::MatrixXd Vsoup, Csoup, Nsoup;
    Eigen::MatrixXi Fsoup;
    to_soup(V, F, FL, Vsoup, Fsoup, Nsoup, Csoup);
    Eigen::MatrixXi CsoupI;
    CsoupI.resizeLike(Csoup);
    for(int vidx = 0; vidx< Csoup.rows();++vidx)
    {
        for(int j : {0,1,2})
        {
            CsoupI(vidx,j) = std::round(Csoup(vidx,j)* 255);
        }
    }
    writePlyColor("../../blenders/test.ply", Vsoup, Fsoup, Nsoup, CsoupI);
}