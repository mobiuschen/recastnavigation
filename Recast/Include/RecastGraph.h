#ifndef RECASTGRAPH_H
#define RECASTGRAPH_H

#include "RecastAlloc.h"

typedef unsigned short GraphID;
typedef unsigned short Weight;

class rcContext;
struct rcPolyMesh;
struct rcGraph;

static const unsigned short MAX_POLY_NUM = 5120;
static const unsigned short RC_GRAPH_ID_NULL = 0;

struct rcGraph
{
    unsigned short  poly;
    unsigned short  nverts;
    GraphID         id;
    GraphID*        verts;      ///< In level 0, the element is poly id. In higher level, the element is graph id. [Element: index * nverts]
    Weight*         weights;    ///< The weight of vertices. [Length: #nverts]
    Weight*         edgeMatrix;
};

struct rcGraphSet
{
    rcGraph*    graphs;     ///< All the graphs in different level.
    int         ngraphs;
    GraphID*    topGraphs;  /// Top graph of diferent levels.
    int         maxedge;
};


rcGraphSet* rcAllocGraphSet(rcContext* ctx, int npoly, int level);
bool        rcBuildGraphSet(rcContext* ctx, const rcPolyMesh& pmesh, rcGraphSet& graphSet, const int level);
void        rcFreeGraphSet(rcGraphSet* pGraphSet);

//rcGraph* rcAllocGraph(rcGraph* graphPool, unsigned short& poolSize, GraphID& retGraphID);
void        rcFreeGraph(rcGraph* pGraph);
bool        rcBuildGraph(rcContext* ctx, rcGraph& graph, const GraphID id, const unsigned short nverts, const rcAllocHint hint);

#endif