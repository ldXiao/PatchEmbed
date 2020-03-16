#include "degenerate_clean.h"
#include <utility>
#include <map>
#include <igl/triangle_triangle_adjacency.h>
#include <vector>
#include <queue>
namespace bcclean {
    void _triangle_degeneracy(
        const Eigen::MatrixXd & V,
        const Eigen::RowVector3i & T,
        double & degeneracy,
        int & bottom_idx
    )
    {
        // bottom means the largest edge
        double maxb = -1;
        double hdb = -1;
        int e0 = 0;
        int vx = T(e0);
        int vy = T((e0+1)%3);
        int vz = T((e0+2)%3);
        Eigen::RowVector3d X, Y, Z;
        X = V.row(vx);
        Y = V.row(vy);
        Z = V.row(vz);
        double Area;
        Area = ((Y-X).cross(Z- X)).norm();
        for(auto ej : {0,1,2})
        {
            int vx = T(ej);
            int vy = T((ej+1)%3);
            int vz = T((ej+2)%3);
            Eigen::RowVector3d X, Y, Z;
            X = V.row(vx);
            Y = V.row(vy);
            Z = V.row(vz);
            double b;
            b = (Y-X).norm();
            if(b > maxb)
            {
                maxb = b;
                bottom_idx = ej;
            }
        }
        degeneracy = Area / (maxb * maxb);
    }

    std::vector<std::pair<int,int> > collect_degenerate_triangle(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        const double threshold
    )
    {
        std::vector<std::pair<int, int> > res;
        for(int fidx =0 ; fidx < F.rows(); fidx++)
        {
            int bottomidx;
            double degeneracy;
            _triangle_degeneracy(V, F.row(fidx), degeneracy, bottomidx);
            if(degeneracy < threshold)
            {
                res.push_back(std::make_pair(fidx, bottomidx));
            }
        }
        return res;
    }

    bool detect_inconsistent_degentrig(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        const Eigen::MatrixXi & TT, 
        const Eigen::VectorXi & FL,
        const std::vector<std::pair<int, int > > & degen_list,
        const std::map<int, bool> & handled_dict,
        std::pair<int,int> & merge_pair
    )
    {
        bool find = false;
        for(auto item: degen_list)
        {
            if(! handled_dict.at(item.first))
            {
                if(FL(item.first)!= FL(TT(item.first, item.second)))
                {
                    merge_pair = item;
                    return true;
                }
                else if(
                    (FL(item.first)!= FL(TT(item.first, (item.second+1)%3))) && 
                    (FL(item.first)!= FL(TT(item.first, (item.second+2)%3)))
                )
                {
                    int vx = (item.second);
                    int vy = (item.second+1)%3;
                    int vz = (item.second+2)%3;
                    Eigen::RowVectorXd X, Y, Z;
                    X = V.row(vx);
                    Y = V.row(vy);
                    Z = V.row(vz);
                    double norm1, norm2;
                    norm1 = (Y-Z).norm();
                    norm2 = (Z-X).norm();
                    if (norm1 > norm2)
                    {
                        merge_pair = std::make_pair(item.first, vy);
                        return true;
                    }
                    else 
                    {
                        merge_pair = std::make_pair(item.first, vz);
                        return true;
                    }
                }
            }
        }
        return false;
    }

    void connected_relabel(
        const Eigen::MatrixXi & TT,
        Eigen::VectorXi & FL
    )
    {
        std::map<int, bool> label_visited;
        int label_num = 0;
        for(int fidx =0; fidx < FL.rows(); fidx ++)
        {
            int lb = FL(fidx);
            if(label_visited.find(lb)== label_visited.end())
            {
                label_visited[lb]= false;
                label_num +=1;
            }
        }
        Eigen::VectorXi F_visit=Eigen::VectorXi::Constant(FL.rows(),-1);
        int root = 0;
        bool finish = false;
        int cc_count = 0;
        while(!finish)
        {
            int curlb = FL(root);
            int fill_lb = cc_count;
            cc_count +=1;
            // if(label_visited.at(curlb))
            // {
            //     fill_lb = label_num;
            //     label_num ++;
            // }
            std::queue<int> visit_queue;
            visit_queue.push(root);
            F_visit(root)= 1;
            std::vector<int> component_list;
            while(visit_queue.size()!=0)
            {
                int cur_f = visit_queue.front();
                for(int j = 0 ; j  < TT.cols(); j++)
                {
                    int fjdx = TT(cur_f, j);
                    if(FL(fjdx)==curlb && (F_visit(fjdx)==-1))
                    {
                        F_visit(fjdx)= 1;
                        visit_queue.push(fjdx);
                    }
                }
                FL(cur_f) = fill_lb;
                visit_queue.pop();
            }

            finish = true;
            for(int jj = 0 ; jj < FL.rows(); jj++)
            {
                if(F_visit(jj)==-1)
                {
                    root = jj;
                    finish = false;
                    break;
                }
            }
        }
    }

    void degenerate_clean(
        const Eigen::MatrixXd & V,
        const Eigen::MatrixXi & F,
        Eigen::VectorXi & FL,
        const double threshold
    )
    {
        std::vector<std::pair<int, int> > degen_list = collect_degenerate_triangle(V, F, threshold);
        Eigen::MatrixXi TT;
        igl::triangle_triangle_adjacency(F, TT);
        std::map<int, bool> handled_dict;
        for(auto item: degen_list)
        {
            handled_dict[item.first]= false;
        }
        std::pair<int, int> merge_pair;
        while(detect_inconsistent_degentrig(V, F, TT, FL, degen_list,handled_dict, merge_pair))
        {
            FL(merge_pair.first) = FL(TT(merge_pair.first, merge_pair.second));
            handled_dict[merge_pair.first]= true;
        }
        connected_relabel(TT, FL);
        return;
    }
}