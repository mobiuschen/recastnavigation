#ifndef RECASTGRAPH_H
#define RECASTGRAPH_H


typedef unsigned short GraphID;
typedef unsigned short Weight;

class rcContext;
struct rcPolyMesh;
struct rcGraph;
//struct rcGraphEdge;

const unsigned short MAX_POLY_NUM = 5120;

struct rcGraph
{
    unsigned short poly;
    GraphID id;
    GraphID* verts;         ///< In level 0, the element is poly id. In higher level, the element is graph id. [Element: index * nverts]
    Weight* weights;        ///< The weight of vertices. [Length: #nverts]
    Weight* edgeMatrix;
    int nverts;
};

//struct rcGraphEdge
//{
//    GraphID adjvert;
//    Weight weight;
//    rcGraphEdge* next;
//};

struct rcGraphSet
{
    rcGraph* graphs;            ///< All the graphs in different level.
    int ngraphs;
    GraphID* topGraphs;             /// Top graph of diferent levels.
    int maxedge;
};


rcGraphSet* rcAllocGraphSet(rcContext* ctx, int npoly, int level);
bool rcBuildGraphSet(rcContext* ctx, const rcPolyMesh& pmesh, rcGraphSet& graphSet, const int level);
void rcFreeGraphSet(rcGraphSet* pGraphSet);
void rcFreeGraph(rcGraph* pGraph);
#endif