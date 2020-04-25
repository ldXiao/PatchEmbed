#include <map>
#include <vector>
#include <iostream>

#include <igl/bfs_orient.h>
#include <regex>
#include <filesystem>
#include <string>
#include <igl/read_triangle_mesh.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <igl/writeMESH.h>
#include <igl/writeOBJ.h>
#include <igl/upsample.h>
#include <igl/facet_components.h>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>
#include <igl/bounding_box_diagonal.h>
#include <Eigen/Core>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include "bcclean.h"
#include "Betti.h"
#include "edge.h"
#include "patch.h"
#include "MatchMakerTree.h"
#include "fTetwild.h"
#include "polyline_distance.h"
#include "params.h"
#include "orientation_check.h"
#include "CellularGraph.h"
#include "MatchMakerDynamic.h"
#include <cxxopts.hpp>
#include <igl/writePLY.h>
#include <nlohmann/json.hpp>
#include <unordered_map>
#include <tuple>
#include <yaml-cpp/yaml.h>
#include <chrono> 
using VFL = std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi>;

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
            ("u, upsp", "upsample stages", cxxopts::value<int>())
            ("b, btthreshold", "backtrack threshold", cxxopts::value<double>())
            ("d, data_root","data root", cxxopts::value<std::string>())
            ("t, tracing", "tracing",cxxopts::value<std::string>());

    auto args = options.parse(argc, argv);
    // Load a mesh in OBJ format
    bool re_tet = false;
    std::string data_root, tracing;
    json param_json;
    bcclean::params param;
    int stop_eng = 10;
    {

        param.data_root = args["data_root"].as<std::string>();
        param.iden = false;
        param.upsp = args["upsp"].as<int>();
        param.debug = true;
        param.guard_len_r = 0;
        re_tet = false;
        param.edge_len_r = 0.01;
        tracing = args["tracing"].as<std::string>();
        param.tracing = tracing;
        param.stop_eng = 10;
        param.merge_threshold = 0;
        param.backtrack_threshold = args["btthreshold"].as<double>();
    }
    std::string bad_mesh_file, face_label_dmat, face_label_yml;
    std::regex r(".*trimesh.*\\.obj");
    std::regex yr(".*features.*\\.yml");
    for (const auto & entry : std::filesystem::directory_iterator(param.data_root))
    {
        if(std::regex_match(entry.path().c_str(), r))
        {
            bad_mesh_file = entry.path().c_str();
            std::printf("got bad mesh: %s\n",bad_mesh_file.c_str());
        }
        if(std::regex_match(entry.path().c_str(),yr))
        {
            face_label_yml = entry.path().c_str();   
        }
    }
    
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
    {
        std::unordered_map<int, int> face_label_dict; 
        YAML::Node conf = YAML::LoadFile(face_label_yml);
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
        FL_bad = Eigen::VectorXi::Constant(count, 0);
        for(auto item: face_label_dict)
        {
            FL_bad(item.first) = item.second;
        } 
        igl::writeDMAT(param.data_root+"/feat.dmat", FL_bad);
    }

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
        Eigen::MatrixXi CCF_copy = CCF_bad;
        Eigen::VectorXi C;
        igl::bfs_orient(CCF_copy,CCF_bad,C);
        
        Eigen::MatrixXd CCV_good;
        Eigen::MatrixXi CCF_good;
        Eigen::VectorXi CCFL_good;
        std::string output_file_bad, output_file_good, output_label_good;
        output_file_bad = param.data_root + "/"+ "CC"+std::to_string(cc)+"-bad.obj";
        output_file_good = param.data_root +"/"+"CC"+std::to_string(cc)+"edg_len_r"+std::to_string(param.edge_len_r)+"-good.obj";
        output_label_good = param.data_root +"/"+"CC"+std::to_string(cc)+"edg_len_r"+std::to_string(param.edge_len_r); +"-label.dmat";
        auto file_logger = spdlog::basic_logger_mt("basic_logger", param.data_root + "/"+ "CC"+std::to_string(cc)+"/logs.txt",true);
        bool file_exists = false;
        bool finished_before =true;
        std::string out_FL  =  param.data_root + "/"+ "CC"+std::to_string(cc)+"/FL_"+param.tracing+".dmat";
        if(std::filesystem::exists("helloworld.txt")){
            Eigen::MatrixXi testFL;
            igl::readDMAT(out_FL, testFL);
            if(testFL.minCoeff()!=-1)
            {
                continue; //go to next component
            }
        }
        for (const auto & entry : std::filesystem::directory_iterator(param.data_root))
        {
            if(entry.path().string()==output_file_good)
            {
                file_exists = true;
                std::printf("good mesh exists for components %d\n", cc);
            }
        }
        
        if(!file_exists || re_tet){

            std::string out_tet_path= param.data_root+"/CC"+std::to_string(cc)+"/tet.mesh";
            bcclean::Tet::fTetwild(CCV_bad, CCF_bad, param.edge_len_r,stop_eng, CCV_good, CCF_good, out_tet_path);
            igl::writeOBJ(output_file_good, CCV_good, CCF_good);
        }
        else{
            igl::read_triangle_mesh(output_file_good, CCV_good, CCF_good); 
        }
        bcclean::flip_orientation_ifnecessary( CCV_good, CCF_good, CCV_bad, CCF_bad);
        bcclean::patch::SetStatics(CCV_bad, CCF_bad, CCFL_bad, CClabel_num);
        bcclean::CollectPatches();
        std::string CC_work_dir = param.data_root+"/CC"+std::to_string(cc);
        std::filesystem::create_directories(CC_work_dir);
        std::ofstream o3(CC_work_dir+"/result"+tracing+".json");
        json result_json;
        int betti_bad=bcclean::Betti(CCV_bad, CCF_bad);
        int betti_good=bcclean::Betti(CCV_good, CCF_good);
        if(betti_bad!= betti_good)
        {
            std::cout << "Inconsisitant topology, abort" << std::endl;
            result_json["consistant topology"] = false;
            continue;
        }
        {
            Eigen::VectorXi CCFC;
            igl::facet_components(CCF_good,CCFC);
            if(CCFC.maxCoeff()>0)
            {
                std::cout << "Inconsisitant topology, abort" << std::endl;
                result_json["consistant topology"] = false;
                continue;
            }
        }
        result_json["consistant topology"] = true;

        if(param.upsp> 0)
        {
            igl::upsample(CCV_good, CCF_good, param.upsp);
        }
        result_json["upsp"]= param.upsp;
        int label_num = bcclean::patch::FL_mod.maxCoeff()+1;
        if(label_num == 1)
        {
            result_json["succeed"] = true;
            result_json["maxerr"] = -1;
            o3 << result_json;
            continue;
        }
        bool succeed= false;
        bcclean::params param_copy = param;
        param_copy.data_root = CC_work_dir;
        int vini = CCV_good.rows();
        bcclean::CellularGraph cg;
        int remain_patch_num = -1;
        cg = bcclean::CellularGraph::GenCellularGraph(bcclean::patch::Vbase, bcclean::patch::Fbase, bcclean::patch::FL_mod);
        if(cg._nodes.size()==0)
        {
            result_json["succeed"] = true;
            result_json["maxerr"] = -1;
            o3 << result_json;
            continue;
        }
        bcclean::result_measure rm;
        // / Get starting timepoint 
        auto start = std::chrono::high_resolution_clock::now(); 
        if(tracing=="loop"){
            try{
                succeed = bcclean::MatchMaker::MatchMakerPatch(cg,CCV_good, CCF_good, CCFL_good, param_copy, file_logger, rm);
            }
            catch(...)
            {
                std::cout<< "failed" <<std::endl;
                succeed = false;
            }
        } else if (tracing == "tree")
        {
            try{
            succeed=bcclean::MatchMaker::MatchMakerTree(cg,CCV_good, CCF_good, CCFL_good, param_copy, file_logger, rm); 
            }
            catch(...)
            {
                std::cout<< "failed" <<std::endl;
                succeed = false;
            }
        }
        auto stop = std::chrono::high_resolution_clock::now();
        if(param.debug)
        {
            
            std::string path_file;
            if(param.tracing == "tree")
            {
                path_file ="/debug_paths_tree.json";
            }
            else
            {
                path_file = "/debug_paths.json";
            }
            std::ifstream i1(CC_work_dir+path_file);
            std::ifstream i2(CC_work_dir+"/debug_path_bad.json");
            
            json path_good, path_bad;
            path_good = json::parse(i1);
            path_bad  = json::parse(i2);
            double max_error= -1;
            for(auto item: path_good.items())
            {
                double err;
                try{
                    err = bcclean::Eval::hausdorff1d(bcclean::patch::Vbase, path_bad[item.key()], CCV_good, item.value());
                }
                catch(...)
                {
                    continue;
                }
                if(max_error< err)
                {
                    max_error = err;
                }
            }
            double dd = igl::bounding_box_diagonal(bcclean::patch::Vbase);
            std::cout<< "finished, maxerr:" << max_error/dd << std::endl; 
            result_json["succeed"]= succeed;
            result_json["maxerr"]= max_error/dd;
            auto duration = std::chrono::duration_cast<std::chrono::minutes>(stop - start); 
            result_json["duration"] = duration.count();
            result_json["vnum_ini"] = vini;
            result_json["vnum_fin"] = CCV_good.rows();
            result_json["remain_patch_num"] =rm.remain_patch_num;
            result_json["area_rel_diff"] = rm.area_rel_diff;
            result_json["peri_rel_diff"] = rm.peri_rel_diff;
            o3 << result_json;
        }

    }

    return 0;
}
