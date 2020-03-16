#ifndef BCCLEAN_PARAMS_H
#define BCCLEAN_PARAMS_H
#include <string>
namespace bcclean{
    class params{
        public:
        std::string data_root;
        bool debug;
        bool iden;
        double edge_len_r;
        double guard_len_r;
        double stop_eng;
        double merge_threshold;
        double backtrack_threshold;
        double area_threshold;
        int upsp; 
        params(){
            data_root = ".";
            debug = true;
            iden = false;
            edge_len_r = 0.01;
            guard_len_r = 0.01;
            stop_eng = 10;
            upsp = 0;
            backtrack_threshold = 1.0;
            area_threshold = 2.0;
            merge_threshold = 0.01;
        }
    };
}
#endif // BCCLEAN_PARAM_H