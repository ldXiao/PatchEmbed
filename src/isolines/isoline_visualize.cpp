#include <iostream>
#include <fstream>
#include <igl/copyleft/cgal/wire_mesh.h>
#include <igl/readOBJ.h>
#include <igl/list_to_matrix.h>
#include <nlohmann/json.hpp>
#include <igl/bounding_box_diagonal.h>
#include <igl/remove_unreferenced.h>
#include <string>
#include <igl/exact_geodesic.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
using json = nlohmann::json;
int main(int argc, char * argv[]){
    std::string obj_s = argv[1]; // obj_in
    double r = std::stod(argv[2]);
    
    Eigen::MatrixXd Vs, Vt;
    Eigen::MatrixXi Es, Ft;
    
    igl::readOBJ(obj_s, Vs, Es);
    double h = igl::bounding_box_diagonal(Vs) * r;
    Es.conservativeResize(Es.rows(),2);
    Eigen::VectorXi J;
    igl::copyleft::cgal::wire_mesh(Vs, Es, h,3,false,Vt,Ft,J);
    igl::writeOBJ("../../blenders/isolines/ioslines_c.obj",Vt,Ft);
    return 0;
}