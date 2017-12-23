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

bool rcBuildGraph(rcContext* ctx, const rcPolyMesh* pmesh, rcGraph& graph)
{
    unsigned int length = 0;
    int maxEdages = 0;
    int nedges = 0;

    if (pmesh == nullptr)
        return false;

    graph.nverts = pmesh->npolys;
    maxEdages = (graph.nverts - 1) * 2;

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

    //In Level 0, every graph represents a poly.
    for (int i = 0, n = pmesh->npolys; i < n; i++)
    {
        unsigned short* p = &(pmesh->polys[i * pmesh->nvp * 2]);
        unsigned short* polyVerts = &(pmesh->polys[i * pmesh->nvp * 2]);
        unsigned short* polyAdjs = &(p[pmesh->nvp]);

        // set vertex
        graph.verts[i].verts = nullptr;
        graph.verts[i].weights = nullptr;
        graph.verts[i].edges = nullptr;
        graph.verts[i].nverts = 0;
        graph.verts[i].poly = (unsigned short)i;

        // build edges
        for (int j = 0, m = pmesh->nvp; j < m; j++)
        {
            if (polyVerts[j] == RC_MESH_NULL_IDX || polyAdjs[j] == RC_MESH_NULL_IDX)
                break;

            unsigned short adjPoly = polyAdjs[j];
            int bigP = i > adjPoly ? adjPoly : i;
            int smallP = i > adjPoly ? i : adjPoly;

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

            if (findEdge)
                break;

            edges[nedges] = (unsigned short)smallP;
            edges[nedges + 1] = (unsigned short)bigP;
            nedges++;
        }//for
    }//for

    length = nedges * 2;
    graph.edges = (unsigned short*)rcAlloc(sizeof(unsigned short) * length, RC_ALLOC_PERM);
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