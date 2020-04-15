#include <iostream>
#include <fstream>
#include <igl/copyleft/cgal/wire_mesh.h>
#include <igl/readOBJ.h>
#include <igl/list_to_matrix.h>
#include <nlohmann/json.hpp>
#include <igl/bounding_box_diagonal.h>
#include <igl/remove_unreferenced.h>
#include <string>
using json = nlohmann::json;
int main(int argc, char * argv[]){
    std::string obj_file = argv[1];
    std::ifstream path_stream(argv[2]);
    json path_json;
    double r = std::stod(argv[3]);
    path_json = json::parse(path_stream);
    Eigen::MatrixXd V, VV,VW;
    Eigen::MatrixXi WE, F,FW;
    igl::readOBJ(obj_file, V, F);
    std::vector<std::vector<int>> WE_list;
    for(auto item: path_json.items())
    {
        std::vector<int> path = path_json[item.key()].get<std::vector<int> >();
        for(int i = 0; i < path.size()-1; ++i){
            WE_list.push_back({path[i], path[i+1]});
        }
    }
    igl::list_to_matrix(WE_list, WE);
    Eigen::VectorXi I, J;
    igl::remove_unreferenced(V,WE,VV, WE, I, J);
    double h = igl::bounding_box_diagonal(V) * r;
    igl::copyleft::cgal::wire_mesh(VV, WE, h,3,false,VW,FW,J);
    igl::writeOBJ("../../blenders/wiremesh.obj",VW,FW);
    return 0;
}