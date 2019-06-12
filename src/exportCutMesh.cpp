//
// Created by lind on 6/12/19.
//
#include "CutMesh.h"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <cxxopts.hpp>
#include <Eigen/StdVector>
Eigen::MatrixXd V;
Eigen::MatrixXi F;
OTMapping::CutMesh CM;
using namespace std;
using namespace OTMapping;
int main(int argc, char *argv[]){
    using namespace std;
    using namespace Eigen;
    if (argc < 2) {
        cout << "Usage cutmeshtest_bin --infile mesh.obj --outfile cutmesh.obj" << endl;
        exit(0);
    }
    cxxopts::Options options("CutMesh", "One line description of MyProgram");
    options.add_options()
            ("i, infile", "File name", cxxopts::value<std::string>())
            ("o, outfile","Sample number", cxxopts::value<std::string>())
            ("c, cut", "cut or separate(for cube only)", cxxopts::value<bool>())
            ;
    auto args = options.parse(argc, argv);

    // Load a mesh in OBJ format
    igl::readOBJ(args["infile"].as<std::string>(), V, F);
    auto func = [](Eigen::Vector3d x)->double{return std::sin(x[0]+x[1]+x[2]);};
    CM.set_initial(V, F,10, func);
    std::cout << "holy shit1" << std::endl;
    if(args["cut"].as<bool>()) {
        Eigen::Vector3d p0(0, 0, 0);
        Eigen::Vector3d p1(0, 0, 0);
        Eigen::Vector3d n1(0, 0, 1);
        Eigen::Vector3d n0(0, 1, 0.5);
        CM.cut_with(p0, n0, p1, n1);
    }
    else{
        CM.separate_cube_faces();
        cout << "sepraate"<<endl;
    }
    std::cout << "holy shit2" << std::endl;
    CM.perturb(3, 0.1, 0.02);
    std::cout << CM.ComponentsVertices[0]->array() << std::endl;
    Eigen::MatrixXd Vout(CM.ComponentsVertices[0]->rows(),3);
    Vout<<CM.ComponentsVertices[0]->array();
    Eigen::MatrixXi Fout(CM.ComponentsFaces[0]->rows(),3);
    Fout<<(CM.ComponentsFaces[0]->array());
    std::cout << "holy shit3" << std::endl;
    int countv=Vout.rows();
    int countf = Fout.rows();
    for(int i =1; i< CM.ComponentsVertices.size(); ++i){
        // Concatenate (VA,FA) and (VB,FB) into (V,F)
        Eigen::MatrixXd Vtemp(countv+CM.ComponentsVertices[i]->rows(),CM.ComponentsVertices[i]->cols());
        Vtemp<<Vout,*(CM.ComponentsVertices[i]);
        Vout = Vtemp;
        Eigen::MatrixXi Ftemp(countf+CM.ComponentsFaces[i]->rows(),CM.ComponentsFaces[i]->cols());
        std::cout << Fout << std::endl;
        std::cout << "________________________________" << std::endl;
        Ftemp<<Fout,(CM.ComponentsFaces[i]->array()+countv);
        countf += CM.ComponentsFaces[i]->rows();
        countv += CM.ComponentsVertices[i]->rows();
        Fout = Ftemp;
    }
    std::cout << "holy shit4" << std::endl;
    igl::writeOBJ(args["outfile"].as<std::string>(),Vout,Fout);
    return 0;

}

