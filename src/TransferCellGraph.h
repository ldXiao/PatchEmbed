#include "CellularGraph.h"
#include "TraceComplex.h"
namespace bcclean{
namespace Bijection{
    void TransferCellGraph(
        const  CellularGraph & cg,
        const  MatchMaker::TraceComplex & tc,
        CellularGraph & cgt
    );
}
}