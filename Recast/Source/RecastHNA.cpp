#include <string.h>
#include "Recast.h"
#include "RecastHNA.h"


rcGraph *rcAllocGraph()
{
    rcGraph *graph = (rcGraph *) rcAlloc(sizeof(rcGraph), RC_ALLOC_PERM);
    memset(graph, 0, sizeof(rcGraph));
    return graph;
}

bool rcBuildGraph(rcContext *ctx, const rcPolyMesh &pmesh, rcGraph &graph)
{
    bool result = false;
    const int npolys = pmesh.npolys;
    int nvts = 0;
    int allocLen = 0;

    rcVertex *vtxs = nullptr;
    Weight *adjMatrix = nullptr;

    allocLen = npolys;
    vtxs = (rcVertex *) rcAlloc(sizeof(rcVertex) * allocLen, RC_ALLOC_PERM);
    if (vtxs == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'vtxs' (%d)", allocLen);
        goto Exit0;
    }
    memset(vtxs, RC_INVALID_VERTEX, sizeof(rcVertex) * allocLen);

    allocLen = npolys * npolys;
    adjMatrix = (Weight *) rcAlloc(sizeof(Weight) * allocLen, RC_ALLOC_PERM);
    if (adjMatrix == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'adjMatrix' (%d)", allocLen);
        goto Exit0;
    }
    memset(adjMatrix, 0, sizeof(Weight) * allocLen);


    for (int i = 0; i < npolys; i++)
    {
        if (pmesh.polys[i] == RC_MESH_NULL_IDX)
            break;

        unsigned short *poly = pmesh.polys + i * pmesh.nvp * 2;
        unsigned short *padj = poly + pmesh.nvp;
        for (int j = 0; j < pmesh.nvp; j++)
        {
            unsigned short ipoly = padj[j];
            if (ipoly == RC_MESH_NULL_IDX)
                continue;

            if (adjMatrix[i * npolys + ipoly] == 0)
            {
                adjMatrix[i * npolys + ipoly] = 1;
                adjMatrix[ipoly * npolys + i] = 1;
            }
        }

        rcVertex &v = vtxs[i];
        v.ipoly = i;
        v.vwgt = 1;
        v.nedges = 0; // set in next step
        v.iedges = i * npolys;
        v.cewgt = 0;
        v.adjwgt = 0; // set in next step
        nvts++;
    }

    for (int i = 0; i < nvts; i++)
    {
        rcVertex &v = vtxs[i];
        for (int j = 0; j < nvts; j++)
        {
            if (adjMatrix[i * npolys + j] == 0)
            {
                v.adjwgt += adjMatrix[i * npolys + j];
                v.nedges++;
            }
        }//for
    }//for


    graph.vtxs = vtxs;
    graph.adjncy = adjMatrix;
    graph.nvt = nvts;

    result = true;
Exit0:
    if (!result)
    {
        if (vtxs != nullptr)
        {
            rcFree(vtxs);
            vtxs = nullptr;
        }

        if (adjMatrix != nullptr)
        {
            rcFree(adjMatrix);
            adjMatrix = nullptr;
        }
    }
    return result;
}
