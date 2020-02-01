#include <regex>
#include <filesystem>
#include <map>
#include <string>
#include <igl/read_triangle_mesh.h>
#include <igl/readDMAT.h>
#include <igl/facet_components.h>
#include <igl/writeDMAT.h>
#include <igl/writeOBJ.h>
#include <cxxopts.hpp>
#include <Eigen/Dense>
#include <igl/remove_unreferenced.h>
#include <tuple>
#include "fTetwild.h"
using VFL = std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>;
namespace fs = std::filesystem;
void decomposeVFL(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F, const Eigen::VectorXi & FL, std::map<int, VFL> & vfls, std::map<int, std::map<int, int> > & cc_label_map)
{
    vfls.clear();
    cc_label_map.clear();
    Eigen::VectorXi FC;
    igl::facet_components(F, FC);
    int numCC = FC.maxCoeff()+1;
    std::map<int, int> CC_fnums;
    std::map<int, std::map<int, int> > CC_l_log;
    std::map<int, Eigen::MatrixXi> CC_faces;
    std::map<int, Eigen::VectorXi> CC_labels;
    for(int cc= 0; cc< numCC; ++cc)
    {
        CC_faces[cc] = Eigen::MatrixXi::Constant(F.rows(),3,0);
        CC_labels[cc] = Eigen::VectorXi::Constant(F.rows(),0);
        CC_fnums[cc]=0;
        CC_l_log[cc] = std::map<int,int>();
        cc_label_map[cc] = std::map<int, int>();
    }
    for(int fidx =0 ; fidx< F.rows(); ++fidx)
    {
        int cc = FC(fidx);
        CC_faces[cc].row(CC_fnums[cc]) = F.row(fidx);
        
        if(CC_l_log[cc].find(FL(fidx)) == CC_l_log[cc].end())
        {
            int cur_cc_lnum = CC_l_log[cc].size();
            CC_l_log[cc][FL(fidx)] = cur_cc_lnum;
            cc_label_map[cc][cur_cc_lnum] = FL(fidx);
        } 
        CC_labels[cc](CC_fnums[cc]) = CC_l_log[cc][FL(fidx)];
        CC_fnums[cc]+=1;
    }
    for(int cc= 0; cc< numCC; ++cc)
    {
        CC_faces[cc].conservativeResize(CC_fnums[cc],3);
        CC_labels[cc].conservativeResize(CC_fnums[cc]);
        Eigen::MatrixXi Fcc = CC_faces[cc];
        Eigen::MatrixXd Vcc;
        Eigen::VectorXi FLcc = CC_labels[cc];
        Eigen::VectorXi Icc;
        igl::remove_unreferenced(V,Fcc, Vcc, Fcc, Icc);
        vfls[cc] = std::make_tuple(Vcc, Fcc, FLcc);
    }
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cout << "Usage ftet_gen -w working_directory -e edge_len_r -s stop_energy" << std::endl;
        exit(0);
    }
    cxxopts::Options options("ftet_gen", "One line description of MyProgram");
    options.add_options()
            ("w, workdir", "working directory", cxxopts::value<std::string>())
            ("e, edgelen", "edge lenth relative", cxxopts::value<double>())
            ("s, stopeng", "stoping energy", cxxopts::value<double>())
            ;
    auto args = options.parse(argc, argv);
    double edge_len_r = args["edgelen"].as<double>();
    double stop_eng = args["stopeng"].as<double>();
    std::string data_root = args["workdir"].as<std::string>();
    std::string face_label_dmat = data_root + "/"+ "feat.dmat";
    std::string bad_mesh_file, good_mesh_file, face_label_yml;
    std::regex r(".*trimesh.*\\.obj");
    for (const auto & entry : std::filesystem::directory_iterator(data_root))
    {
        if(std::regex_match(entry.path().c_str(), r))
        {
            bad_mesh_file = entry.path().c_str();
            std::printf("got bad mesh: %s\n",bad_mesh_file.c_str());
        }
    }

    Eigen::MatrixXd V_bad, V_good;
    Eigen::MatrixXi F_bad, F_good;
    Eigen::VectorXi FL_bad, FL_good;
    igl::readDMAT(face_label_dmat, FL_bad);
    igl::read_triangle_mesh(bad_mesh_file, V_bad, F_bad);
    std::map<int, std::map<int, int> > ComponentsLabelMaps;
    std::map<int, VFL> vfls;
    decomposeVFL(V_bad, F_bad, FL_bad, vfls, ComponentsLabelMaps);
    for(int cc= 0; cc < vfls.size(); ++cc)
    {
        Eigen::MatrixXd CCV_bad = std::get<0>(vfls[cc]);
        Eigen::MatrixXi CCF_bad = std::get<1>(vfls[cc]);
        Eigen::VectorXi CCFL_bad = std::get<2>(vfls[cc]);
        int CClabel_num = CCFL_bad.maxCoeff()+1;
        Eigen::MatrixXd CCV_good;
        Eigen::MatrixXi CCF_good;
        Eigen::VectorXi CCFL_good;
        std::string output_file_bad, output_file_good, output_label_good;
        output_file_bad = data_root + "/"+ "CC"+std::to_string(cc)+"-bad.obj";
        output_file_good = data_root +"/"+"CC"+std::to_string(cc)+"edg_len_r"+std::to_string(edge_len_r)+"-good.obj";
        output_label_good = data_root +"/"+"CC"+std::to_string(cc)+"edg_len_r"+std::to_string(edge_len_r); +"-label.dmat";
        bool file_exists = false;
        for (const auto & entry : std::filesystem::directory_iterator(data_root))
        {
            if(entry.path().string()==output_file_good)
            {
                file_exists = true;
                std::printf("good mesh exists for components %d\n", cc);
            }
        }
        
        if(!file_exists){
            bcclean::Tet::fTetwild(CCV_bad, CCF_bad, edge_len_r, stop_eng, CCV_good, CCF_good);
            igl::writeOBJ(output_file_good, CCV_good, CCF_good);
        }
        else{
            continue;
        }
    }
}