#include <string.h>
#include <stdlib.h>
#include "Recast.h"
#include "RecastHNA.h"
#include "RecastGraph.h"

typedef int* Map;
typedef int* Match;

bool partitionGraph(rcContext* ctx, const rcGraphHNA& graph);

bool coarsening(rcContext* ctx, const rcGraphHNA& graph, const int condK, const float condR);

bool initialPartition();

bool uncoarseningPhase();


bool coarsenOnce(rcContext* ctx, const rcGraphHNA& curGraph, rcGraphHNA& retCoarserGraph, Map& retMap, Match& retMatch);
bool heavyEdgeMatch(rcContext* ctx, const rcGraphHNA& graph, int* retMatch, int* retMap, int& retNewVertNum);

rcGraphHNA* rcAllocGraph(rcAllocHint allocHint)
{
    rcGraphHNA* graph = (rcGraphHNA*) rcAlloc(sizeof(rcGraphHNA), allocHint);
    memset(graph, 0, sizeof(rcGraphHNA));
    return graph;
}

bool rcFreeGraph(rcGraphHNA* graph)
{
    rcFree(graph->vtxs);
    graph->vtxs = nullptr;

    rcFree(graph->adjncy);
    graph->adjncy = nullptr;

    rcFree(graph);
    graph = nullptr;

    return true;
}


bool rcBuildGraphHNA(rcContext* ctx, rcGraphHNA& graph, int nverts, rcAllocHint allocHint)
{
    bool result = false;
    int nvts = 0;
    int allocLen = 0;
    rcVertex* vtxs = nullptr;
    Weight* adjMatrix = nullptr;

    allocLen = nverts;
    vtxs = (rcVertex*) rcAlloc(sizeof(rcVertex) * allocLen, allocHint);
    if (vtxs == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'vtxs' (%d)", allocLen);
        goto Exit0;
    }
    memset(vtxs, 0, sizeof(rcVertex) * allocLen);

    allocLen = nverts * nverts;
    adjMatrix = (Weight*) rcAlloc(sizeof(Weight) * allocLen, allocHint);
    if (adjMatrix == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'adjMatrix' (%d)", allocLen);
        goto Exit0;
    }
    memset(adjMatrix, 0, sizeof(Weight) * allocLen);

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

bool rcBuildGraphHNA(rcContext* ctx, const rcPolyMesh& pmesh, rcGraphHNA& graph)
{
    bool result = false;
    bool retCode = false;

    retCode = rcBuildGraphHNA(ctx, graph, pmesh.npolys, RC_ALLOC_TEMP);
    if (!retCode)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraphHNA: 'rcBuildGraphHNA()' fails.");
        goto Exit0;
    }

    for (int i = 0, n = graph.nvt; i < n; i++)
    {
        if (pmesh.polys[i] == RC_MESH_NULL_IDX)
            break;

        unsigned short* poly = pmesh.polys + i * pmesh.nvp * 2;
        unsigned short* padj = poly + pmesh.nvp;
        for (int j = 0; j < pmesh.nvp; j++)
        {
            unsigned short ipoly = padj[j];
            if (ipoly == RC_MESH_NULL_IDX)
                continue;

            if (graph.adjncy[i * n + ipoly] == 0)
            {
                graph.adjncy[i * n + ipoly] = 1;
                graph.adjncy[ipoly * n + i] = 1;
            }
        }

        rcVertex& v = graph.vtxs[i];
        v.ipoly = i;
        v.vwgt = 1;
        v.nedges = 0; // set in next step
        v.iedges = i * n;
        v.cewgt = 0;
        v.adjwgt = 0; // set in next step
    }

    for (int i = 0, n = graph.nvt; i < n; i++)
    {
        rcVertex& v = graph.vtxs[i];
        for (int j = 0, m = graph.nvt; j < m; j++)
        {
            if (graph.adjncy[i * n + j] > 0)
            {
                v.adjwgt += graph.adjncy[i * n + j];
                v.nedges++;
            }
        }//for
    }//for


    result = true;
Exit0:
    return result;
}


bool partitionGraph(rcContext* ctx, const rcGraphHNA& graph)
{
    bool result = false;
    bool retCode = false;

    retCode = coarsening(ctx, graph, 10, 0.8f);
    if (!retCode)
    {
        ctx->log(RC_LOG_ERROR, "partitionGraph: coarsening phase fails.");
        goto Exit0;
    }

    retCode = initialPartition();
    if (!retCode)
    {
        ctx->log(RC_LOG_ERROR, "initialPartition: coarsening phase fails.");
        goto Exit0;
    }

    retCode = uncoarseningPhase();
    if (!retCode)
    {
        ctx->log(RC_LOG_ERROR, "uncoarseningPhase: coarsening phase fails.");
        goto Exit0;
    }


    result = true;
Exit0:
    return result;
}


bool coarsenOnce(rcContext* ctx, const rcGraphHNA& curGraph, rcGraphHNA& retCoarserGraph, Map& retMap, Match& retMatch)
{
    bool result = false;
    bool retCode = false;
    int newVertNum = 0;

    retCode = heavyEdgeMatch(ctx, curGraph, retMatch, retMap, newVertNum);
    if (!retCode)
    {
        ctx->log(RC_LOG_ERROR, "coarsenOnce: 'heavyEdgeMatch()' fails");
        goto Exit0;
    }

    rcBuildGraphHNA(ctx, retCoarserGraph, newVertNum, RC_ALLOC_TEMP);
    for (int i = 0, n = curGraph.nvt; i < n; i++)
    {
        rcVertex& v = curGraph.vtxs[i];
        int iu1 = retMap[i];
        rcVertex& u = retCoarserGraph.vtxs[iu1];
        u.ipoly = RC_MESH_NULL_IDX;
        u.vwgt += v.vwgt;
        u.iedges = 0; // set in next step.

        for (int j = 0, m = curGraph.nvt; j < m; j++)
        {
            if (i == j)
                continue;

            Weight ewgt = curGraph.adjncy[v.iedges + j];
            if (ewgt == 0)
                continue;

            int iu2 = retMap[j];
            if (iu2 != iu1)
            {
                // inter edge
                retCoarserGraph.adjncy[iu1 * curGraph.nvt + iu2] += ewgt;
                retCoarserGraph.adjncy[iu2 * curGraph.nvt + iu1] += ewgt;
                u.adjwgt += ewgt;

            } else
            {
                // intra edge
                u.cewgt += ewgt;
            }
        }//for
    }//for

    for (int i = 0, n = retCoarserGraph.nvt; i < n; i++)
    {
        rcVertex& u = retCoarserGraph.vtxs[i];
        for (int j = 0, m = retCoarserGraph.nvt; j < m; j++)
        {
            if (i != j && retCoarserGraph.adjncy[u.iedges + j] > 0)
            {
                u.nedges += 1;
            }
        }
    }

    result = true;
Exit0:
    return result;
}


bool coarsening(rcContext* ctx, const rcGraphHNA& graph, const int condK, const float condR)
{
    bool result = false;
    bool retCode = false;
    const int max_level = 64;
    Match matches[max_level] = {0};
    Map maps[max_level] = {0};
    rcGraphHNA* graphs[max_level] = {nullptr};

    int level = 0;
    float radio = 1.0f;
    do
    {
        int allocLen = 0;
        const rcGraphHNA* curGraph = level == 0 ? (&graph) : graphs[level - 1];
        rcGraphHNA* coarserGraph = rcAllocGraph(RC_ALLOC_TEMP);
        if (coarserGraph == nullptr)
        {
            ctx->log(RC_LOG_ERROR, "coarsening: Out of memory 'crtGraph', level=%d", level);
            goto Exit0;
        }

        allocLen = curGraph->nvt;
        Match match = (int*) rcAlloc(sizeof(int) * allocLen, RC_ALLOC_TEMP);
        if (match == nullptr)
        {
            ctx->log(RC_LOG_ERROR, "coarsening: Out of memory 'match' (%d)", allocLen);
            goto Exit0;
        }
        memset((int*) match, RC_INVALID_VERTEX, sizeof(match));

        allocLen = curGraph->nvt;
        Map map = (int*) rcAlloc(sizeof(int) * allocLen, RC_ALLOC_TEMP);
        if (map == nullptr)
        {
            ctx->log(RC_LOG_ERROR, "coarsening: Out of memory 'map' (%d)", allocLen);
            goto Exit0;
        }
        memset((int*) map, RC_INVALID_VERTEX, sizeof(map));

        retCode = coarsenOnce(ctx, *curGraph, *coarserGraph, map, match);
        if (!retCode)
        {
            ctx->log(RC_LOG_ERROR, "'coarsenOnce' fails");
            goto Exit0;
        }

        matches[level] = match;
        maps[level] = map;
        graphs[level] = coarserGraph;
        level++;
        radio = ((float) coarserGraph->nvt) / graph.nvt;
    } while (level < max_level && radio > condR && graphs[level]->nvt > condK);


    result = true;
Exit0:
    return result;
}


bool heavyEdgeMatch(rcContext* ctx, const rcGraphHNA& graph, int* retMatch, int* retMap, int& retNewVertNum)
{
    bool result = false;
    int* shuffleVerts = nullptr;
    int allocLen = 0;

    if (retMatch == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "heavyEdgeMatch: 'retMatch' == nullptr");
        goto Exit0;
    }

    if (retMap == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "heavyEdgeMatch: 'map' == nullptr");
        goto Exit0;
    }

    allocLen = graph.nvt;
    shuffleVerts = (int*) rcAlloc(sizeof(int) * allocLen, RC_ALLOC_TEMP);
    if (shuffleVerts == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "heavyEdgeMatch: Out of memory 'shuffleVerts' (%d)", allocLen);
        goto Exit0;
    }
    //copy
    for (int i = 0, n = graph.nvt; i < n; i++)
    {
        shuffleVerts[i] = i;
    }
    /// Fisher–Yates shuffle
    for (int i = 0, n = graph.nvt - 3; i < n; i++)
    {
        int j = rand() % n;
        int temp = shuffleVerts[i];
        shuffleVerts[i] = temp;
        shuffleVerts[j] = shuffleVerts[i];
    }

    /// Heavy Edge Matching
    for (int i = 0, n = graph.nvt; i < n; i++)
    {
        int index = shuffleVerts[i];
        if (retMatch[index] == RC_INVALID_INDEX)
            continue;

        rcVertex& v = graph.vtxs[index];
        if (v.nedges == 0)
        {
            retMatch[index] = index;
            retMap[index] = retNewVertNum;
            retNewVertNum++;
            continue;
        }

        int heaviestEdge = RC_MESH_NULL_IDX;
        int maxWeight = 0;
        for (int j = 0, m = graph.nvt; j < m; j++)
        {
            Weight w = graph.adjncy[v.iedges + j];
            if (w > maxWeight)
            {
                maxWeight = w;
                heaviestEdge = j;
            }
        }//for
        retMatch[index] = heaviestEdge;
        retMatch[heaviestEdge] = index;
        retMap[index] = retNewVertNum;
        retMap[heaviestEdge] = retNewVertNum;
        retNewVertNum++;
    }//for

    result = true;
Exit0:
    return result;
}