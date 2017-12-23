#ifndef RECASTGRAPH_H
#define RECASTGRAPH_H


struct rcGraph
{
    unsigned short poly;        ///< The poly which this graph represents. It's valid only on level 0.
    rcGraph* verts;             ///< The graph vertices. [Length: #nverts]
    unsigned short* weights;    ///< The weight of vertices. [Length: #nverts]
    unsigned short* edges;      ///< The graph edges. [Length: #nverts * 2]
    int nverts;                 ///< The number of graph vertices.
};


struct rcPolyMesh;
class rcContext;

rcGraph* rcAllocGraph();
bool rcBuildGraph(rcContext* ctx, const rcPolyMesh* pmesh, rcGraph& graph);

#endif