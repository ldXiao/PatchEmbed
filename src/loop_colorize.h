#ifndef BCCLEAN_LOOP_COLORIZE_H
#define BCCLEAN_LOOP_COLORIZE_H
#include <Eigen/Core>
#include <vector>
#include <igl/vertex_triangle_adjacency.h> 
#include <queue>
#include <igl/triangle_triangle_adjacency.h>
namespace bcclean
{
    void loop_colorize(
    const Eigen::MatrixXd& V, 
    const Eigen::MatrixXi & F, 
    const std::vector<std::vector<int> > & TEdges,
    const int face_seed,
    const int lb,
    Eigen::VectorXi & FL)
    {
        // for unlabeled faced FL should be -1
        //label the faces encompassed by the loop with lb in FL
        assert(FL.rows() == F.rows());
        assert(TEdges.size() == F.rows());
        assert(FL(face_seed)== -1); //the seed faces has to be uninitialized at first
        // use BFS to search over TTi starting with face_seed
        Eigen::MatrixXi TT, TTi;
        igl::triangle_triangle_adjacency(F, TT, TTi);
        // start with face 0 use TT and TTi info to get connected components
        // traverse all faces connected to face_seed and do not cross cuts label them to be lb;
        int root_face = face_seed;
        std::queue<int> search_queue;
        search_queue.push(face_seed);
        while(search_queue.size()!=0){
            int cur_face = search_queue.front();
            search_queue.pop(); // remove head
            FL(cur_face) = lb;
            Eigen::RowVector3i adjs = TT.row(cur_face);
            for(int j =0 ; j <3 ; ++j){
                int face_j = adjs(j);
                if(face_j == -1) {
                    continue;
                }
                // int vidx, vidy;
                // vidx = Fraw(cur_face, j);
                // vidy  = Fraw(cur_face, (j+1)%3);
                if(TEdges[cur_face][j] == -1){
                    // the edge is not in VCuts
                    if(FL(face_j)== -1){
                        // not visited before;
                        FL(face_j) = lb;
                        search_queue.push(face_j);
                    }
                }
            }
        }
    }
}
#endif //BCCLEAN_LOOP_COLORIZE_H