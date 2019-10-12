#include <Eigen/Core>
#include <vector>
#include <map>
#include "node.h"
#include "edge.h"
namespace bcclean
{
    class patch{
        // all the patches is assumed to be planar
        public:
        static Eigen::MatrixXd Vbase;
        static Eigen::MatrixXi Fbase;
        static std::vector<node> node_list; // list of all nodes
        static std::vector<edge> edge_list; // list of all edges;
        static int total_label_num;
        static Eigen::VectorXi FL; // a copy of initial face labels
        static Eigen::VectorXi FL_mod; // Face labels after modification
        static void SetStatics(const Eigen::MatrixXd & Vbase, const Eigen::MatrixXi & Fbase, const Eigen::VectorXi & FL, size_t total_label_num);
        Eigen::MatrixXd Vraw;
        Eigen::MatrixXi Fraw;
        // base mesh for the patch extracted from the label infomation
        int label;
        bool SimplyConnected;
        Eigen::VectorXi VI;
        Eigen::VectorXi FI;
        std::vector<std::vector<int> > loop_edges; // indices of  the edges in the loop
        // finally we wan to make sure each patch contains only one loop
        // and the loop is not self-touching
   };
   void CollectPatches();
} // namespace nam bcclean