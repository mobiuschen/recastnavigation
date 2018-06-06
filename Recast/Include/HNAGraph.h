#ifndef HNAGRAPH_H
#define HNAGRAPH_H

#include "RecastAlloc.h"
#include <vector>

typedef int Index;
typedef unsigned short Weight;

static const int RC_INVALID_INDEX = -1;

class rcContext;


struct rcHNAVertex
{
    Index           ipoly;  ///The index of poly in rcPolyMesh.
    Weight          vwgt;   ///The weight of vertex
    size_t          nedges; ///The size of the adjacency list of v
    Index           iedges; ///The index into Adjncy that is the beginning of the adjacency list of v
    Weight          cewgt;  ///the weight of the edges that have been contracted to create vertex (if vertexis a multinode)
    Weight          adjwgt; ///the sum of the weight of the edges adjacent to v
};

struct rcHNAGraph
{
    rcHNAVertex*    vtxs;
    Weight*         adjncy; ///the adjacency lists of the vertices
    size_t          nvt;    ///The number of vertices
};

rcHNAGraph* rcAllocHNAGraph(rcAllocHint allocHint);
bool        rcFreeGraph(rcHNAGraph* pGraph);
bool        rcBuildGraphHNA(rcContext* ctx, const rcPolyMesh& pmesh, rcHNAGraph& graph);
bool        rcBuildGraphHNA(rcContext* ctx, rcHNAGraph& graph, int nverts, rcAllocHint allocHint);

#endif