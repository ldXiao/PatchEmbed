#include <yaml-cpp/yaml.h>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <string>
#include <unordered_map>
#include <igl/writeDMAT.h>
#include <Eigen/Dense>
namespace fs = std::filesystem;
int main(int argc, char *argv[])
{
    std::unordered_map<int, int> face_label_dict; 
    std::string yaml_file = argv[1];
    fs::path yaml_path = yaml_file;
    std::string par_dir = yaml_path.parent_path(); 
    YAML::Node conf = YAML::LoadFile(yaml_file);
    int lb =0;
    int count = 0;
    for( auto surf: conf["surfaces"])
    {
        for(auto fidx: surf["face_indices"])
        {
            face_label_dict[fidx.as<int>()] = lb;
            count +=1;
        }
        lb +=1;
    }
    Eigen::VectorXi fl = Eigen::VectorXi::Constant(count, 0);
    for(auto item: face_label_dict)
    {
        fl(item.first) = item.second;
    } 
    igl::writeDMAT(par_dir+"/feat.dmat", fl);
}