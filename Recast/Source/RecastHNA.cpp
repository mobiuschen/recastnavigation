#include <string.h>
#include <stdlib.h>
#include "Recast.h"
#include "RecastHNA.h"


rcGraph* rcAllocGraph()
{
    rcGraph* graph = (rcGraph*)rcAlloc(sizeof(rcGraph), RC_ALLOC_PERM);
    memset(graph, 0, sizeof(rcGraph));
    return graph;
}

bool rcBuildGraph(rcContext* ctx, const rcPolyMesh& pmesh, rcGraph& graph)
{
    bool result = false;
    const int npolys = pmesh.npolys;
    int allocLen = 0;

    rcVertex* vtxs = nullptr;
    Weight* adjlist = nullptr;


    allocLen = npolys;
    vtxs = (rcVertex*)rcAlloc(sizeof(rcVertex) * allocLen, RC_ALLOC_PERM);
    if (vtxs == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'vtxs' (%d)", allocLen);
        goto Exit0;
    }
    memset(vtxs, 0, sizeof(rcVertex) * allocLen);

    allocLen = npolys * npolys;
    adjlist = (Weight*)rcAlloc(sizeof(Weight) * allocLen, RC_ALLOC_PERM);
    if (adjlist == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'adjlist' (%d)", allocLen);
        goto Exit0;
    }
    memset(adjlist, 0, sizeof(Weight) * allocLen);


    for (int i = 0; i < npolys; i++)
    {
        unsigned short* poly = pmesh.polys + i * pmesh.nvp * 2;
        unsigned short* padj = poly + pmesh.nvp;
        for (int j = 0; j < npolys; j++)
        {
            unsigned short d = padj[j];

            if (d == RC_MESH_NULL_IDX)
                break;


        }

        rcVertex& v = vtxs[i];
        v.ipoly = i;
        v.vwgt = 1;
        v.nedges = 0; //TODO
        v.iedges = RC_INVALID_INDEX;
        v.cewgt = 0;
        v.adjwgt = 0;
    }


    graph.vtxs = vtxs;
    graph.adjncy = adjlist;
    graph.nvt = npolys;

    result = true;
Exit0:
    if (!result)
    {
        if (vtxs != nullptr)
        {
            rcFree(vtxs);
            vtxs = nullptr;
        }

        if (adjlist != nullptr)
        {
            rcFree(adjlist);
            adjlist = nullptr;
        }
    }
    return result;
}
