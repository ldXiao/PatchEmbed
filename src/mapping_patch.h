#ifndef BCCLEAN_MAPPING_PATCH_H
#define BCCLEAN_MAPPING_PATCH_H
#include <map>
#include <vector>
#include <Eigen/Core>
#include <tuple>
namespace bcclean{
    class node {
        public:
        Eigen::MatrixXd _position;
        int _total_label_num;
        int _occupied_label_num;
        std::map<int, int> _label_occupy_dict;
        bool of_same_type(const node & b);
        bool at_same_position(const Eigen::MatrixXd & position);
        bool initialize(const int total_label_num, const Eigen::MatrixXd & position, const std::vector<int> labels);
    };

    std::vector<std::vector<node>> build_label_nodes_list(
        const Eigen::MatrixXd &V, 
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXi &FL);

    void extract_label_patch_mesh(
        const Eigen::MatrixXd& V, 
        const Eigen::MatrixXi& F, 
        const Eigen::MatrixXi& FL, 
        const int lb_in, 
        Eigen::MatrixXd& V_i, 
        Eigen::MatrixXi& F_i);

    void map_vertices_to_regular_polygon(
        const Eigen::MatrixXd &V, 
        const Eigen::MatrixXi & F, 
        std::vector<node> & nodes, 
        Eigen::VectorXi & bnd,
        Eigen::MatrixXd & bnd_uv,
        std::vector<node>& ordered_nodes);

    class mapping_patch{
        public:
            int label;
            int node_num;
            int total_label_num;
            // global setting
            Eigen::MatrixXd V_raw;
            Eigen::MatrixXi F_raw;
            Eigen::VectorXi bnd_raw;
            //inputed raw mesh boundary degeneracy might exist
            Eigen::MatrixXd V_ndg;
            Eigen::MatrixXi F_ndg;
            Eigen::VectorXi bnd_ndg;
            // modified nondegenerate mesh use this mesh to project to X-Y plane
            
            Eigen::VectorXi V_bridge_raw2ndg;
            Eigen::VectorXi F_bridge_raw2ndg;
            Eigen::VectorXi V_bridge_ndg2raw;
            Eigen::VectorXi F_bridge_ndg2raw;
            // internal correspondence relation in both directions
            

            std::map<int, node> ordered_nodes_list;
            //contains the information index on bnd_ndg and the corresponding nodes

            void
    };
}
#endif //BCCLEAN_MAPPING_PATCH_H