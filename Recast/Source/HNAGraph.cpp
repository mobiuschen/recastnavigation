#include <string.h>
#include <stdlib.h>
#include "HNAGraph.h"
#include "Recast.h"
#include "RecastAssert.h"


bool shuffle(int size, const int* src, int* const dest);

rcHNAGraph* rcAllocHNAGraph(rcAllocHint allocHint)
{
    rcHNAGraph* pGraph = nullptr;        
    pGraph = (rcHNAGraph*)rcAlloc(sizeof(rcHNAGraph), allocHint);
    memset(pGraph, 0, sizeof(rcHNAGraph));
    return pGraph;
}


bool rcFreeGraph(rcHNAGraph* pGraph)
{
    if (pGraph == nullptr)
        return false;

    if (pGraph->vtxs != nullptr)
    {
        rcFree(pGraph->vtxs);
        pGraph->vtxs = nullptr;
    }

    if (pGraph->adjncy != nullptr)
    {
        rcFree(pGraph->adjncy);
        pGraph->adjncy = nullptr;
    }

    rcFree(pGraph);
    pGraph = nullptr;
    return true;
}


bool rcBuildGraphHNA(rcContext* ctx, rcHNAGraph& graph, int nverts, rcAllocHint allocHint)
{
    bool            result = false;
    int             allocLen = 0;
    rcHNAVertex*    vtxs = nullptr;
    Weight*         adjncy = nullptr;

    allocLen = nverts;
    vtxs = (rcHNAVertex*)rcAlloc(sizeof(rcHNAVertex) * allocLen, allocHint);
    if (vtxs == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'vtxs' (%d)", allocLen);
        goto Exit0;
    }
    memset(vtxs, 0, sizeof(rcHNAVertex) * allocLen);

    allocLen = nverts * nverts;
    adjncy = (Weight*)rcAlloc(sizeof(Weight) * allocLen, allocHint);
    if (adjncy == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'adjncy' (%d)", allocLen);
        goto Exit0;
    }
    memset(adjncy, 0, sizeof(Weight) * allocLen);

    graph.vtxs = vtxs;
    graph.adjncy = adjncy;
    graph.nvt = nverts;

    result = true;
Exit0:
    if (!result)
    {
        if (vtxs != nullptr)
        {
            rcFree(vtxs);
            vtxs = nullptr;
        }

        if (adjncy != nullptr)
        {
            rcFree(adjncy);
            adjncy = nullptr;
        }
    }
    return result;
}


bool rcBuildGraphHNA(rcContext* ctx, const rcPolyMesh& pmesh, rcHNAGraph& graph)
{
    bool result = false;
    bool retCode = false;

    retCode = rcBuildGraphHNA(ctx, graph, pmesh.npolys, RC_ALLOC_PERM);
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

        rcHNAVertex& v = graph.vtxs[i];
        v.ipoly = (unsigned short)i;
        v.vwgt = 1;
        v.nedges = 0; // set in next step
        v.iedges = i * n;
        v.cewgt = 0;
        v.adjwgt = 0; // set in next step
    }

    for (int i = 0, n = graph.nvt; i < n; i++)
    {
        rcHNAVertex& v = graph.vtxs[i];
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


rcHNAGraph* rcApplyMatching(rcContext* ctx, const rcHNAGraph& srcGraph, 
                            const Match match, const Map map)
{
    bool    result = false;
    size_t  retVertNum = 0;
    rcHNAGraph* pRetGraph = rcAllocHNAGraph(RC_ALLOC_PERM);

    rcIgnoreUnused(ctx);

    if (pRetGraph == nullptr)
        goto Exit0;
    
    for (int i = 0, n = srcGraph.nvt; i < n; i++)
    {
        retVertNum = (size_t)(map[i]) >= retVertNum ? (map[i] + 1) : retVertNum;
    }

    rcBuildGraphHNA(ctx, *pRetGraph, retVertNum, RC_ALLOC_PERM);

    for (int i = 0, n = pRetGraph->nvt; i < n; i++)
    {
        rcHNAVertex& v = pRetGraph->vtxs[i];
        v.ipoly = RC_MESH_NULL_IDX;
        v.vwgt = 0;     // set in next step
        v.nedges = 0;   // set in next step
        v.iedges = i * n;
        v.cewgt = 0;
        v.adjwgt = 0;   // set in next step
    }

    for (int i = 0, n = srcGraph.nvt; i < n; i++)
    {
        rcHNAVertex& v = srcGraph.vtxs[i];

        Index iu = map[i];
        rcHNAVertex& u = pRetGraph->vtxs[iu];

        u.vwgt += v.vwgt;

        Index iv2 = match[i];
        if (iv2 == i)
        {
            // solo vertex
            u.cewgt = v.cewgt;
            u.adjwgt = v.adjwgt;
        }
        else if (i < iv2)
        {
            rcHNAVertex& v2 = srcGraph.vtxs[iv2];
            Weight innerEdgeWgt = rcGraphGetEdge(ctx, srcGraph, i, iv2);
            u.cewgt = v.cewgt + v2.cewgt + innerEdgeWgt;
            u.adjwgt = v.adjwgt + v2.adjwgt - 2 * innerEdgeWgt;
        }


        for (int j = 0; j < n; j++)
        {
            Index iu2 = map[j];
            Weight wgt = 0;

            if (i <= j)
                continue;

            if (iu == iu2)
                continue;

            wgt = rcGraphGetEdge(ctx, srcGraph, i, j);
            pRetGraph->adjncy[iu * n + iu2] += wgt;
            pRetGraph->adjncy[iu2 * n + iu] += wgt;
        }
    }//for


    for (int i = 0, n = pRetGraph->nvt; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if (rcGraphGetEdge(ctx, *pRetGraph, i, j) > 0)
            {
                pRetGraph->vtxs[i].nedges++;
                pRetGraph->vtxs[j].nedges++;
            }
        }
    }//for

    result = true;
Exit0:
    if (!result)
    {
        rcFree(pRetGraph);
        pRetGraph = nullptr;
    }
    return pRetGraph;
}


bool heavyEdgeMatch(rcContext* ctx, 
                    const rcHNAGraph& graph,
                    Match const retMatch,
                    Map const retPartition,
                    size_t& retNewVertNum)
{
    bool result = false;
    size_t allocLen = graph.nvt;
    Index* vertOrders = (Index*)rcAlloc(sizeof(Index) * allocLen, RC_ALLOC_TEMP);
    
    if (vertOrders == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "heavyEdgeMatch: Out of memory 'vertOrders' (%d)", allocLen);
        goto Exit0;
    }

    for (int i = 0, n = graph.nvt; i < n; i++)
    {
        vertOrders[i] = i;
    }
    shuffle(graph.nvt, vertOrders, vertOrders);

    result = heavyEdgeMatch(ctx, graph, vertOrders, retMatch, retPartition, retNewVertNum);
Exit0:
    if (vertOrders != nullptr)
    {
        rcFree(vertOrders);
        vertOrders = nullptr;
    }
    return result;
}


bool heavyEdgeMatch(rcContext* ctx, 
                    const rcHNAGraph& graph, 
                    const Index* vertOrders,
                    Match const retMatch,
                    Map const retPartition, 
                    size_t& retNewVertNum)
{
    bool result = false;

    if (retMatch == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "heavyEdgeMatch: 'retMatch' == nullptr");
        goto Exit0;
    }

    if (retPartition == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "heavyEdgeMatch: 'retPartition' == nullptr");
        goto Exit0;
    }

    memset(retMatch, RC_INVALID_INDEX, graph.nvt * sizeof(Index));
    memset(retPartition, RC_INVALID_INDEX, graph.nvt * sizeof(Index));

    // Heavy Edge Matching
    retNewVertNum = 0;
    for (int i = 0, n = graph.nvt; i < n; i++)
    {
        int index = vertOrders[i];
        if (retMatch[index] != RC_INVALID_INDEX)
            continue;

        rcHNAVertex& v = graph.vtxs[index];
        if (v.nedges == 0)
        {
            retMatch[index] = index;
            retPartition[index] = retNewVertNum;
            retNewVertNum++;
            continue;
        }

        int heaviestEdge = index;
        int maxWeight = 0;
        for (int j = 0, m = graph.nvt; j < m; j++)
        {            
            if (index == j || retMatch[j] != RC_INVALID_INDEX)
                continue;

            Weight w = rcGraphGetEdge(ctx, graph, index, j);
            if (w > maxWeight)
            {
                maxWeight = w;
                heaviestEdge = j;
            }
        }//for
        retMatch[index] = heaviestEdge;
        retMatch[heaviestEdge] = index;
        retPartition[index] = retNewVertNum;
        retPartition[heaviestEdge] = retNewVertNum;
        retNewVertNum++;
    }//for

    result = true;
Exit0:
    return result;
}



Weight rcGraphGetEdge(rcContext* ctx, const rcHNAGraph& graph, Index v1, Index v2)
{
    rcIgnoreUnused(ctx);
    if (v1 == RC_INVALID_INDEX || v2 == RC_INVALID_INDEX)
        return 0;

    if (v1 < 0 || (size_t)v1 >= graph.nvt || v2 < 0 || (size_t)v2 >= graph.nvt)
        return 0;

    return graph.adjncy[v1 * graph.nvt + v2];
}


bool rcGraphSetEdge(rcContext* ctx, const rcHNAGraph& graph, Index v1, Index v2, Weight ewgt)
{
    rcIgnoreUnused(ctx);
    graph.adjncy[v1 * graph.nvt + v2] = ewgt;
    graph.adjncy[v2 * graph.nvt + v1] = ewgt;
    return true;
}





bool shuffle(int size, const int* src, int* const dest)
{
    /// Fisher¨CYates shuffle
    for (int i = 0, n = size; i < n; i++)
    {
        dest[i] = src[i];
    }
    for (int i = 0, n = size - 3; i < n; i++)
    {
        int j = rand() % n;
        int temp = dest[i];
        dest[i] = temp;
        dest[j] = dest[i];
    }
    return true;
}
