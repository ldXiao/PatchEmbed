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
    bool cyc_flip_mapping(
        std::vector<node> & nodes, 
        std::vector<node> & target_nodes, 
        std::map<int,int> &mapping);
    class mapping_patch{
        public:
            int label;
            int node_num;
            int total_label_num;
            bool loop_direction=true;
            // global setting
            Eigen::MatrixXd _V_raw;
            Eigen::MatrixXi _F_raw;
            Eigen::VectorXi _bnd_raw;
            //inputed raw mesh boundary degeneracy might exist
            Eigen::MatrixXd _V_ndg;
            Eigen::MatrixXi _F_ndg;
            Eigen::VectorXi _bnd_ndg; // boundary vertices indices in loop order
            // modified nondegenerate mesh use this mesh to project to X-Y plane
            // temporarily, we use (_V_ndg, _F_ndg) for mapping only and don't split the corresponding tet mesh
            // TODO add tet split info
            
            Eigen::VectorXi _V_bridge_raw2ndg;
            Eigen::VectorXi _F_bridge_raw2ndg;
            Eigen::VectorXi _V_bridge_ndg2raw;
            Eigen::VectorXi _F_bridge_ndg2raw;
            // internal correspondence relation in both directions

            Eigen::MatrixXd _V_uv;
            Eigen::MatrixXi _F_uv; // bnd_uv should be the same as bnd_ndg
            Eigen::MatrixXd _bnd_uv; // bnd should be the same as bnd_ndg
            // TODO remove ndg triple and preserve only uv ones

            std::vector<int> _nails; // nail index on bnd will not change once initialized
            std::map<int, node> _nails_nodes_dict; //nail index-> node map will not change once initialized
            std::vector<int> _ccw_ordered_nails; // will totate or flip if two patch nodes list does not match
                                                 // store ccw drawing sequence of nodes, store nails       
            std::map<int, double> _edge_arc_ratio_list; //will change once initialized
            //contains the information index on bnd_ndg and the corresponding nodes
            bool build_patch(const Eigen::MatrixXd &Vi, const Eigen::MatrixXi & Fi, std::vector<node> & nodes, int lb_in);
            // void flip_loop_direction(); not yet implemented

    };
}
#endif //BCCLEAN_MAPPING_PATCH_H