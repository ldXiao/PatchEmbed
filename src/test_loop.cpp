#include <map>
#include <vector>
#include <iostream>


#include <regex>
#include <filesystem>
#include <string>
#include <igl/read_triangle_mesh.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/writeOBJ.h>
#include <igl/upsample.h>
#include <igl/facet_components.h>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>
#include <igl/all_edges.h>
#include <Eigen/Core>
#include "bcclean.h"
#include "edge.h"
#include "graphcut_cgal.h"
#include "patch.h"
#include "Match_Maker_Tree.h"
#include "Match_Maker_Loop.h"
#include "fTetwild.h"

#include <cxxopts.hpp>
#include <nlohmann/json.hpp>
#include <unordered_map>
#include <tuple>
using VFL = std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>;

int Betti(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F)
{
    Eigen::MatrixXi E;
    igl::all_edges(F,E);
    return V.rows() - E.rows() + F.rows();
}

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

int main(int argc, char *argv[]){
    using json = nlohmann::json;
    using namespace std;
    

    
    if (argc < 2) {
        std::cout << "Usage test_loop -j jsonfile" << std::endl;
        exit(0);
    }
    
    /*-----------------------------------------
    for json
    */
    cxxopts::Options options("Testloop", "One line description of MyProgram");
    options.add_options()
            ("j, json", "json storing parameters", cxxopts::value<std::string>());
    auto args = options.parse(argc, argv);
    // Load a mesh in OBJ format
    bool debug=true;
    bool iden = false;
    int upsp = 0;
    double edge_len_r=0.01;
    bool only_tet = false;
    std::string data_root, tracing;
    json param_json;
    int stop_eng = 10;
    {
        std::ifstream temp(args["json"].as<std::string>());
        param_json = json::parse(temp);
        std::cout <<"the json parameters are" << param_json << std::endl;
        data_root = param_json["data_root"];
        iden = param_json["iden"];
        upsp = param_json["upsp"];
        debug = param_json["debug"];
        only_tet = param_json["only_tet"];
        edge_len_r = param_json["edge_len_r"];
        tracing = param_json["tracing"];
        stop_eng = param_json["stop_eng"];
    }
    std::string bad_mesh_file, good_mesh_file, face_label_dmat, face_label_yml;
    std::regex r(".*trimesh.*\\.obj");
    for (const auto & entry : std::filesystem::directory_iterator(data_root))
    {
        if(std::regex_match(entry.path().c_str(), r))
        {
            bad_mesh_file = entry.path().c_str();
            std::printf("got bad mesh: %s\n",bad_mesh_file.c_str());
        }
    }
    
    good_mesh_file = data_root + "/"+"good.mesh__sf.obj";
    face_label_dmat = data_root + "/"+ "feat.dmat";
    Eigen::MatrixXd V_bad, V_good;
    Eigen::MatrixXi F_bad, F_good;
    Eigen::VectorXi FL_bad, FL_good;
    igl::read_triangle_mesh(bad_mesh_file, V_bad, F_bad);
    if(iden)
    {
        igl::read_triangle_mesh(bad_mesh_file, V_good, F_good);
    }
    else
    {
        igl::read_triangle_mesh(good_mesh_file, V_good, F_good);
    }
    if(upsp> 0)
    {
        igl::upsample(V_good, F_good, upsp);
    }
    igl::readDMAT(face_label_dmat, FL_bad);
    std::map<int, VFL> vfls;
    std::map<int, std::map<int, int> > ComponentsLabelMaps;
    decomposeVFL(V_bad, F_bad, FL_bad, vfls, ComponentsLabelMaps);
    for(int cc= 0; cc < vfls.size(); ++cc)
    {
        Eigen::MatrixXd CCV_bad = std::get<0>(vfls[cc]);
        Eigen::MatrixXi CCF_bad = std::get<1>(vfls[cc]);
        Eigen::VectorXi CCFL_bad = std::get<2>(vfls[cc]);
        int CClabel_num = CCFL_bad.maxCoeff()+1;
        bcclean::patch pat;
        bcclean::patch::SetStatics(CCV_bad, CCF_bad, CCFL_bad, CClabel_num);
        bcclean::CollectPatches();
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
        if(only_tet){
            bcclean::Tet::fTetwild(CCV_bad, CCF_bad, edge_len_r,stop_eng, CCV_good, CCF_good);
            igl::writeOBJ(output_file_good, CCV_good, CCF_good);
            break;
        }
        if(file_exists){
            igl::read_triangle_mesh(output_file_good, CCV_good, CCF_good);
        }
        else{
            bcclean::Tet::fTetwild(CCV_bad, CCF_bad, edge_len_r,stop_eng, CCV_good, CCF_good);
            igl::writeOBJ(output_file_good, CCV_good, CCF_good); 
        }
        if(tracing=="loop"){
            bcclean::MatchMaker::trace_and_label_loop(bcclean::patch::Vbase, bcclean::patch::Fbase, bcclean::patch::FL_mod, CCV_good, CCF_good, CCFL_good, debug);
        } else if (tracing == "tree")
        {
            bcclean::MatchMaker::trace_and_label(bcclean::patch::Vbase, bcclean::patch::Fbase, bcclean::patch::FL_mod, CCV_good, CCF_good, CCFL_good, debug); 
        }

    }
    // bcclean::Tet::fTetwild(V_good, F_good,0.01, V_good, F_good);
    // int label_num = FL_bad.maxCoeff()+1;
    // bcclean::patch pat;
    // bcclean::patch::SetStatics(V_bad, F_bad, FL_bad, label_num);
    // bcclean::CollectPatches();
    // Eigen::MatrixXd V_good_copy = V_good;
    // Eigen::MatrixXi F_good_copy = F_good;
    // Eigen::VectorXi FL_good_copy = FL_good;
    // bcclean::MatchMaker::trace_and_label_loop(bcclean::patch::Vbase, bcclean::patch::Fbase, bcclean::patch::FL_mod, V_good_copy, F_good_copy, FL_good_copy, debug); 
    // igl::writeDMAT(data_root +"/test.dmat", bcclean::patch::FL_mod );
    // igl::writeOBJ(data_root + "/test.obj", V_good, F_good);

    return 0;
}