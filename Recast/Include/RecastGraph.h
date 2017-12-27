#ifndef RECASTGRAPH_H
#define RECASTGRAPH_H


struct rcGraph
{
    unsigned short poly;        ///< The poly which this graph represents. It's valid only on level 0.
    rcGraph* verts;             ///< The graph vertices. [Length: #nverts]
    unsigned short* weights;    ///< The weight of vertices. [Length: #nverts]
    unsigned short* edgeMatrix; ///< The graph edges. [Length: #nedges * #nedges]. The value of every edge is the weight.
    int nverts;                 ///< The number of graph vertices.
    int nedges;                 ///< The number of edges
};


struct rcPolyMesh;
class rcContext;

rcGraph* rcAllocGraph();
bool rcBuildGraph(rcContext* ctx, const rcPolyMesh& pmesh, rcGraph& graph);
void rcFreeGraph(rcGraph* pGraph);
#endif