#ifndef BCCLEAN_LOOP_COLORIZE_H
#define BCCLEAN_LOOP_COLORIZE_H
#include <Eige/Core>
#include <vector>
#include <igl/vertex_triangle_adjacency.h> 
namespace bcclean
{
    void loop_colorize(
    const Eigen::MatrixXd& V, 
    const Eigen::MatrixXi & F, 
    const std::vector<int>& loop, 
    const int lb,
    Eigen::VectorXi & FL)
    {
        // for unlabeled faced FL should be -1
        //label the faces encompassed by the loop with lb in FL
        assert(FL.rows() == F.rows());
        std::vector<std::vector<bool> > TCuts(F.rows());
        for(int fcount = 0; fcount  < F.rows(); ++fcount)
        {
            TCuts[fcount] = {false, false, false};
        }
        //initialize TCuts
        std::vector<std::vector<int > > VF;// #V list of lists of incident faces (adjacency list)
        std::vector<std::vector<int > > VFi;  //  #V list of lists of index of incidence within incident faces listed in VF
        igl::vertex_triangle_adjacency(V, F, VF, VFi);
        for(int loop_idx=0; loop_idx < loop.size(); ++loop_idx){
            int uidx = loop[loop_idx];
            int vidx = loop[(loop_idx+1)%loop.size()];
            std::vector<int> inter(VF[uidx].size()+ VF[vidx].size());
            auto it = std::set_intersection(VF[uidx].begin(), VF[uidx].end(), VF[vidx].begin(), VF[vidx].end(), inter.begin());
            inter.resize(it-inter.begin());
            // there should be only one comman adjacent triangle for boundary vertices
            for(auto trg: inter){
                for(int edgpos =0; edgpos < 3 ; ++edgpos){
                    int uuidx = F(trg, edgpos);
                    int vvidx = F(trg, (edgpos+1)%3);
                    if(uuidx == uidx && vvidx == vidx){
                        TCuts[trg][edgpos] = true;
                    }
                    if(uuidx == vidx && vvidx == uidx){
                        TCuts[trg][edgpos] = true;
                    }
                }
            }
        }
    }
}
#endif //BCCLEAN_LOOP_COLORIZE_H