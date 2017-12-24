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
    int maxEdages = pmesh.npolys == 0 ? 0 : (pmesh.npolys - 1) * 2;
    int nedges = 0;

    graph.nverts = pmesh.npolys;


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

    length = maxEdages * 2;
    graph.edges = (unsigned short*)rcAlloc(sizeof(unsigned short) * length, RC_ALLOC_PERM);
    if (!graph.edges)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'graph.edges' (%d).", length);
        return false;
    }

    length = maxEdages * 2;
    rcScopedDelete<unsigned short> edges((unsigned short*)rcAlloc(sizeof(unsigned short) * length, RC_ALLOC_TEMP));
    if (!edges)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'edges' (%d).", length);
        return false;
    }
    memset(edges, RC_MESH_NULL_IDX, sizeof(unsigned short) * length);

    const int nvp = pmesh.nvp;
    //In Level 0, every graph represents a poly.
    for (int i = 0, n = pmesh.npolys; i < n; i++)
    {
        unsigned short* p = &(pmesh.polys[i * nvp * 2]);
        rcGraph graphVert = graph.verts[i];
        // set vertex
        graphVert.verts = nullptr;
        graphVert.weights = nullptr;
        graphVert.edges = nullptr;
        graphVert.nverts = 0;
        graphVert.poly = (unsigned short)i;

        // build edges
        for (int j = 0, m = nvp; j < m; j++)
        {
            unsigned short* pAdj = &(p[nvp]);
            if (p[j] == RC_MESH_NULL_IDX) break;
            if (pAdj[j] & 0x8000) continue;

            int bigP = i > pAdj[j] ? pAdj[j] : i;
            int smallP = i > pAdj[j] ? i : pAdj[j];

            bool findEdge = false;
            for (int k = 0; k < nedges; k++)
            {
                unsigned short* edge = &(edges[k * 2]);
                if (edge[0] == RC_MESH_NULL_IDX || edge[1] == RC_MESH_NULL_IDX)
                    break;

                if (edge[0] == smallP && edge[1] == bigP)
                {
                    findEdge = true;
                    break;
                }
            }//for

            if (!findEdge)
            {
                edges[nedges] = (unsigned short)smallP;
                edges[nedges + 1] = (unsigned short)bigP;
                nedges++;
            }
        }//for
    }//for

    length = nedges * 2;
    graph.edges = (unsigned short*)rcAlloc(sizeof(unsigned short) * length, RC_ALLOC_PERM);
    graph.nedges = nedges;
    memcpy(graph.edges, edges, sizeof(unsigned short) * length);

    return true;
}


void rcFreeGraph(rcGraph* pGraph)
{
    if (pGraph == nullptr)
        return;

    rcFree(pGraph->edges);
    rcFree(pGraph->verts);
    rcFree(pGraph->weights);
    rcFree(pGraph);
}