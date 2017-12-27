#include <string.h>
#include "Recast.h"
#include "RecastAlloc.h"
#include "RecastGraph.h"

rcGraph* rcAllocGraph()
{
    rcGraph* result = (rcGraph*)rcAlloc(sizeof(rcGraph), RC_ALLOC_PERM);
    memset(result, 0, sizeof(rcGraph));
    return result;
}

bool rcBuildGraph(rcContext* ctx, const rcPolyMesh& pmesh, rcGraph& graph)
{
    unsigned int length = 0;
    const int npolys = pmesh.npolys;
    const int nvp = pmesh.nvp;
    graph.nverts = npolys;


    length = graph.nverts;
    graph.verts = (rcGraph*)rcAlloc(sizeof(rcGraph) * length, RC_ALLOC_PERM);
    if (!graph.verts)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'graph.verts' (%d).", length);
        return false;
    }

    length = graph.nverts;
    graph.weights = (unsigned short*)rcAlloc(sizeof(unsigned short) * length, RC_ALLOC_PERM);
    if (!graph.weights)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'graph.weights' (%d).", length);
        return false;
    }
    memset(graph.weights, 1, sizeof(unsigned short) * length);

    length = npolys * npolys;
    graph.edgeMatrix = (unsigned short*)rcAlloc(sizeof(unsigned short) * length, RC_ALLOC_PERM);
    if (!graph.edgeMatrix)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'graph.edgeMatrix' (%d).", length);
        return false;
    }
    memset(graph.edgeMatrix, 0, sizeof(unsigned short) * length);

    //In Level 0, every graph represents a poly.
    for (int i = 0; i < npolys; i++)
    {
        unsigned short* p = &(pmesh.polys[i * nvp * 2]);
        rcGraph& graphVert = graph.verts[i];
        // set vertex graph
        graphVert.verts = nullptr;
        graphVert.weights = nullptr;
        graphVert.edgeMatrix = nullptr;
        graphVert.nverts = 0;
        graphVert.poly = (unsigned short)i;

        unsigned short* adjs = &(p[nvp]);
        // build edges
        for (int j = 0, m = nvp; j < m; j++)
        {
            int pAdj = adjs[j];
            if (p[j] == RC_MESH_NULL_IDX) break;
            if (pAdj & 0x8000) continue;

            graph.edgeMatrix[i * npolys + pAdj] = 1;
            graph.edgeMatrix[pAdj * npolys + i] = 1;
        }//for
    }//for

    return true;
}


void rcFreeGraph(rcGraph* pGraph)
{
    if (pGraph == nullptr)
        return;

    rcFree(pGraph->edgeMatrix);
    rcFree(pGraph->verts);
    rcFree(pGraph->weights);
    rcFree(pGraph);
}