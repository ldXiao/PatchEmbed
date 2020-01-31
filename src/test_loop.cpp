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
#include <igl/bounding_box_diagonal.h>
#include <igl/edges.h>
#include <Eigen/Core>
#include "bcclean.h"
#include "edge.h"
#include "graphcut_cgal.h"
#include "patch.h"
#include "Match_Maker_Tree.h"
#include "Match_Maker_Loop.h"
#include "fTetwild.h"
#include "degenerate_clean.h"
#include "polyline_distance.h"
#include "params.h"
#include <cxxopts.hpp>
#include <nlohmann/json.hpp>
#include <unordered_map>
#include <tuple>
using VFL = std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>;

int Betti(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F)
{
    Eigen::MatrixXi E;
    igl::edges(F,E);
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
    bool re_tet = false;
    std::string data_root, tracing;
    json param_json;
    bcclean::params param;
    int stop_eng = 10;
    {
        std::ifstream temp(args["json"].as<std::string>());
        param_json = json::parse(temp);
        std::cout <<"the json parameters are" << param_json << std::endl;
        param.data_root = param_json["data_root"];
        param.iden = param_json["iden"];
        param.upsp = param_json["upsp"];
        param.debug = param_json["debug"];
        param.guard_len_r = param_json["guard_len_r"];
        re_tet = param_json["re_tet"];
        param.edge_len_r = param_json["edge_len_r"];
        tracing = param_json["tracing"];
        param.stop_eng = param_json["stop_eng"];
        param.merge_threshold = param_json["merge_threshold"];
    }
    std::string bad_mesh_file, good_mesh_file, face_label_dmat, face_label_yml;
    std::regex r(".*trimesh.*\\.obj");
    for (const auto & entry : std::filesystem::directory_iterator(param.data_root))
    {
        if(std::regex_match(entry.path().c_str(), r))
        {
            bad_mesh_file = entry.path().c_str();
            std::printf("got bad mesh: %s\n",bad_mesh_file.c_str());
        }
    }
    
    good_mesh_file = param.data_root + "/"+"good.mesh__sf.obj";
    face_label_dmat = param.data_root + "/"+ "feat.dmat";
    Eigen::MatrixXd V_bad, V_good;
    Eigen::MatrixXi F_bad, F_good;
    Eigen::VectorXi FL_bad, FL_good;
    igl::read_triangle_mesh(bad_mesh_file, V_bad, F_bad);
    // if(iden)
    // {
    //     igl::read_triangle_mesh(bad_mesh_file, V_good, F_good);
    // }
    // else
    // {
    //     igl::read_triangle_mesh(good_mesh_file, V_good, F_good);
    // }

    igl::readDMAT(face_label_dmat, FL_bad);
    bcclean::degenerate_clean(V_bad, F_bad, FL_bad, param.merge_threshold);
    igl::writeDMAT("../dbginfo/dgfl.dmat", FL_bad);
    igl::writeOBJ("../dbginfo/dgfl.obj", V_bad, F_bad);
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
        output_file_bad = param.data_root + "/"+ "CC"+std::to_string(cc)+"-bad.obj";
        output_file_good = param.data_root +"/"+"CC"+std::to_string(cc)+"edg_len_r"+std::to_string(param.edge_len_r)+"-good.obj";
        output_label_good = param.data_root +"/"+"CC"+std::to_string(cc)+"edg_len_r"+std::to_string(param.edge_len_r); +"-label.dmat";
        bool file_exists = false;
        for (const auto & entry : std::filesystem::directory_iterator(param.data_root))
        {
            if(entry.path().string()==output_file_good)
            {
                file_exists = true;
                std::printf("good mesh exists for components %d\n", cc);
            }
        }
        
        if(!file_exists || re_tet){
            bcclean::Tet::fTetwild(CCV_bad, CCF_bad, param.edge_len_r,stop_eng, CCV_good, CCF_good);
            igl::writeOBJ(output_file_good, CCV_good, CCF_good);
            break;
        }
        else{
            igl::read_triangle_mesh(output_file_good, CCV_good, CCF_good); 
        }
        int betti_bad=Betti(CCV_bad, CCF_bad);
        int betti_good=Betti(CCV_good, CCF_good);
        if(betti_bad!= betti_good)
        {
            std::cout << "Inconsisitant topology, abort" << std::endl;
        }
        {
            Eigen::VectorXi CCFC;
            igl::facet_components(CCF_good,CCFC);
            if(CCFC.maxCoeff()>0)
            {
                std::cout << "Inconsisitant topology, abort" << std::endl;
                return EXIT_FAILURE; 
            }
        }
        if(param.upsp> 0)
        {
            igl::upsample(CCV_good, CCF_good, param.upsp);
        }
        if(tracing=="loop"){
            bcclean::MatchMaker::trace_and_label_loop(bcclean::patch::Vbase, bcclean::patch::Fbase, bcclean::patch::FL_mod, CCV_good, CCF_good, CCFL_good, param);
        } else if (tracing == "tree")
        {
            bcclean::MatchMaker::trace_and_label(bcclean::patch::Vbase, bcclean::patch::Fbase, bcclean::patch::FL_mod, CCV_good, CCF_good, CCFL_good, param.debug); 
        }
        if(param.debug)
        {
            json path_good, path_bad;
            path_good = json::parse("../dbginfo/debug_paths.json");
            path_bad  = json::parse("../dbginfo/debug_paths_bad.json");
            double max_error= -1;
            for(auto item: path_good.items())
            {
                double err;
                err = bcclean::Eval::hausdorff1d(CCV_bad, path_bad[item.key()], CCV_good, item.value());
                if(max_error< err)
                {
                    max_error = err;
                }
            }
            double dd = igl::bounding_box_diagonal(CCV_bad);
            std::cout<< "finished, maxerr:" << max_error/dd << std::endl; 
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