#ifndef RECASTHNA_H
#define RECASTHNA_H

#include "RecastAlloc.h"


typedef unsigned short GraphID;
typedef unsigned short Weight;

class rcContext;

struct rcPolyMesh;
struct rcGraphHNA;

static const int RC_INVALID_INDEX = -1;
static const unsigned int RC_INVALID_VERTEX = 0xFFFF;

struct rcVertex
{
    int ipoly;  ///The index of poly in rcPolyMesh.
    Weight vwgt;   ///The weight of vertex
    unsigned short nedges; ///The size of the adjacency list of v
    int iedges; ///The index into Adjncy that is the beginning of the adjacency list of v
    Weight cewgt;  ///the weight of the edges that have been contracted to create vertex (if vertexis a multinode)
    Weight adjwgt; ///the sum of the weight of the edges adjacent to v
};

struct rcGraphHNA
{
    rcVertex* vtxs;
    Weight* adjncy; ///the adjacency lists of the vertices
    int nvt;    ///The number of vertices
};


rcGraphHNA* rcAllocGraph(rcAllocHint hint);

bool rcFreeGraph();

bool rcBuildGraphHNA(rcContext* ctx, const rcPolyMesh& pmesh, rcGraphHNA& graph);

bool rcBuildGraphHNA(rcContext* ctx, rcGraphHNA& graph, int nverts, rcAllocHint allocHint);

#endif