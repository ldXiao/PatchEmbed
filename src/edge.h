#ifndef BCCLEAN_EDGE_H
#define BCCLEAN_EDGE_H
#include <Eigen/Core>
#include <vector>
#include <map>
#include "node.h"
namespace bcclean{
    class edge {
        public:
        node head;
        node tail;
        const int total_label_num;
    }
}

#endif // BCCLEAN_EDGE_H