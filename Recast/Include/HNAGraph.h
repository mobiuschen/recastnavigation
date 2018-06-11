#ifndef HNAGRAPH_H
#define HNAGRAPH_H

#include "RecastAlloc.h"

typedef int Index;
typedef unsigned short Weight;
typedef Index* Map;
typedef Index* Match;

static const int RC_INVALID_INDEX = -1;

class rcContext;
struct rcPolyMesh;


struct rcHNAVertex
{
    unsigned short  ipoly;  ///The index of poly in rcPolyMesh.
    Weight          vwgt;   ///The weight of vertex
    size_t          nedges; ///The size of the adjacency list of v
    Index           iedges; ///The index into Adjncy that is the beginning of the adjacency list of v
    Weight          cewgt;  ///the weight of the edges that have been contracted to create vertex (if vertexis a multinode)
    Weight          adjwgt; ///the sum of the weight of the edges adjacent to v
};

struct rcHNAGraph
{
    rcHNAVertex*    vtxs;
    Weight*         adjncy; ///the adjacency table of the vertices
    size_t          nvt;    ///The number of vertices
};

rcHNAGraph* rcAllocHNAGraph(rcAllocHint allocHint);
bool        rcFreeGraph(rcHNAGraph* pGraph);
bool        rcBuildGraphHNA(rcContext* ctx, rcPolyMesh& pmesh, rcHNAGraph& graph);
bool        rcBuildGraphHNA(rcContext* ctx, rcHNAGraph& graph, int nverts, rcAllocHint allocHint);

bool        rcGraphSetEdge(rcContext* ctx, const rcHNAGraph& graph, Index v1, Index v2, Weight ewgt);
Weight      rcGraphGetEdge(rcContext* ctx, const rcHNAGraph& graph, Index v1, Index v2);
rcHNAGraph* rcApplyMatching(rcContext* ctx, const rcHNAGraph& srcGraph,
                            const Match match, const Map map);

bool        heavyEdgeMatch(rcContext* ctx, 
                           const rcHNAGraph& graph, 
                           const Index* queue,
                           Match const retMatch,
                           Map const retPartition,
                           size_t& retNewVertNum);
bool        heavyEdgeMatch(rcContext* ctx,
                           const rcHNAGraph& graph,
                           Match const retMatch,
                           Map const retPartition,
                           size_t& retNewVertNum);

#endif