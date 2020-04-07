#include "loop_colorize.h"
namespace bcclean
{
    std::pair<int, double> loop_colorize(
    const Eigen::MatrixXd& V, 
    const Eigen::MatrixXi & F, 
    const std::vector<std::vector<int> > & TEdges,
    const int face_seed,
    const int lb,
    Eigen::VectorXi & FL)
    {
        // return a dble area of the the total colored region
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
        double DblA = 0;
        Eigen::RowVector3d a_ = V.row(F(root_face,1)) - V.row(F(root_face,0));
        Eigen::RowVector3d b_ = V.row(F(root_face,2)) - V.row(F(root_face,0));
        DblA += (a_.cross(b_)).norm();
        std::queue<int> search_queue;
        search_queue.push(face_seed);
        int count = 1;
        std::map<int, bool> visit_dict;
        for(int fidx= 0 ; fidx < F.rows(); ++fidx)
        {
            visit_dict[fidx] = false;
        }
        visit_dict[root_face] = true;
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
                    if(FL(face_j)== -1 && !(visit_dict[face_j])){
                        // not visited before;
                        search_queue.push(face_j);
                        visit_dict[face_j] = true;
                        FL(cur_face) = lb;
                        Eigen::RowVector3d a__ = V.row(F(face_j,1)) - V.row(F(face_j, 0));
                        Eigen::RowVector3d b__ = V.row(F(face_j,2)) - V.row(F(face_j,0));
                        DblA += (a__.cross(b__)).norm();
                        count += 1;
                    }
                }
            }
        }

        return std::make_pair(count,DblA);
    }
}
