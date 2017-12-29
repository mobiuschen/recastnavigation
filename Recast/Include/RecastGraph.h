#ifndef RECASTGRAPH_H
#define RECASTGRAPH_H

// G0, rcGraph.verts == nullptr;
// G1, rcGraph.verts ?[Element: index]
// G1, rcGraph.verts ?[Element: rcGraph*] 
// G1, rcGraph.verts ?[Element: rcGraph] 

typedef unsigned short GraphID;
typedef unsigned short Weight;

struct rcGraph;


struct rcGraph
{
    unsigned short level;
    GraphID id;
    GraphID* verts;         ///< In level 0, the element is poly id. In higher level, the element is graph id. [Element: index * nverts]
    Weight* weights;        ///< The weight of vertices. [Length: #nverts]
    Weight* edgeMatrix;
    int nverts;
};

struct rcGraphEdge
{
    GraphID adjvert;
    Weight weight;
    rcGraphEdge* next;
};

struct rcGraphSet
{
    rcPolyMesh* pmesh;          ///
    rcGraph* graphs;            ///< All the graphs in different level.
    int ngraphs;
    int* topGraphs;  /// Top graph of diferent levels.
    int maxedge;
};

struct rcPolyMesh;
class rcContext;

rcGraph* rcAllocGraph();
bool rcBuildGraph(rcContext* ctx, const rcPolyMesh& pmesh, rcGraph& graph);
void rcFreeGraph(rcGraph* pGraph);
#endif