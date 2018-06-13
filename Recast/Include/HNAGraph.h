#ifndef HNAGRAPH_H
#define HNAGRAPH_H

#include "RecastAlloc.h"

typedef int Index;
typedef int Weight;
typedef Index* Map;
typedef Index* Match;
typedef const Index* ConstMap;
typedef const Index* ConstMatch;

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
bool        rcBuildGraphHNA(rcContext* ctx, const rcPolyMesh& pmesh, rcHNAGraph& graph);
bool        rcBuildGraphHNA(rcContext* ctx, rcHNAGraph& graph, int nverts, rcAllocHint allocHint);

bool        rcGraphSetEdge(rcContext* ctx, const rcHNAGraph& graph, Index v1, Index v2, Weight ewgt);
Weight      rcGraphGetEdge(rcContext* ctx, const rcHNAGraph& graph, Index v1, Index v2);
rcHNAGraph* rcApplyMatching(rcContext* ctx, const rcHNAGraph& srcGraph,
                            ConstMatch match, ConstMap map, size_t newVertNum);

bool        heavyEdgeMatch(rcContext* ctx, 
                           const rcHNAGraph& graph, 
                           const Index* queue,
                           Match retMatch,
                           Map retPartition,
                           size_t& retNewVertNum);

///Given a partition P, the number of edges whose incident vertices belong to dierent subsets is called the **Edge-Cut** of the partition.
Weight      calcEdgeCut(rcContext* ctx, const rcHNAGraph& graph, ConstMap partition);

/// The gain g(v) of a vertex v is defined as the reduction on the edge-cut if vertex v moves from one partition to the other.
Weight      calcKLGain(rcContext* ctx, const rcHNAGraph& graph, Index v, ConstMap partition);


/// The functions below are only for testing.
bool        heavyEdgeMatch(rcContext* ctx,
                           const rcHNAGraph& graph,
                           Match const retMatch,
                           Map const retPartition,
                           size_t& retNewVertNum);

#endif