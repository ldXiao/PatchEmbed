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
    std::string obj_t = argv[2];
    std::string mapping = argv[3];
    Eigen::MatrixXd Vs, Vt;
    Eigen::MatrixXi Fs, Ft;
    igl::readOBJ(obj_s, Vs, Fs);
    igl::readOBJ(obj_t, Vt, Ft);
    Eigen::MatrixXd M_t2s;
    
    // igl::exact_geodesic
    // igl::writeOBJ("../../blenders/wiremesh.obj",VW,FW);
    return 0;
}