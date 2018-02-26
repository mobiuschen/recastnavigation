#include <string.h>
#include <stdlib.h>
#include "Recast.h"
#include "RecastAssert.h"
#include "RecastHNA.h"
#include "RecastGraph.h"

typedef int* Partition;
typedef int* Match;

static const int MAX_COARSEN_LEVEL = 64;

struct rcHNAConfig
{
    int     condK;
    float   condR;
    int     gggpTimes;
};

struct CoarsenData
{
    rcGraphHNA* levelGraphs[MAX_COARSEN_LEVEL] = {nullptr};
    Match       matches[MAX_COARSEN_LEVEL] = {0};
    Partition   partitions[MAX_COARSEN_LEVEL] = {0};
    int         nlevel;
};

struct klGainBucketCell
{
    int iv;
    Weight gain;
    klGainBucketCell* pre;
    klGainBucketCell* next;
};



bool partitionGraph(rcContext* ctx, const rcGraphHNA& graph);

bool coarsening(rcContext* ctx, const rcGraphHNA& graph, const rcHNAConfig& conf, CoarsenData& intermediateData);

bool initialPartition(rcContext* ctx, const rcGraphHNA& graph, const rcHNAConfig& conf, CoarsenData& intermediateData);

bool uncoarseningPhase();


bool greedyGraphGrowingPartition(const rcGraphHNA& graph, const Partition& retPartition);
bool coarsenOnce(rcContext* ctx, const rcGraphHNA& curGraph, rcGraphHNA& retCoarserGraph, Partition& retPartition, Match& retMatch);
bool heavyEdgeMatch(rcContext* ctx, const rcGraphHNA& graph, int* retMatch, Partition retPartition, int& retNewVertNum);
bool shuffle(const int size, const int* src, int* dest);
Weight calcKLGain(const rcGraphHNA& graph, const int iv, const Partition& p);
Weight calcEdgeCut(const rcGraphHNA& graph, const Partition& p);


bool refinePartition(const rcGraphHNA& graph, Partition& partition);
bool swapVert(const rcGraphHNA& graph, Partition& partition, klGainBucketCell* gainTbl, const int iv, const int k);
bool insertToGainBucket(klGainBucketCell* link, klGainBucketCell* item);
bool removeFromGainBucket(klGainBucketCell* item);



rcGraphHNA* rcAllocGraph(rcAllocHint allocHint)
{
    rcGraphHNA* graph = (rcGraphHNA*)rcAlloc(sizeof(rcGraphHNA), allocHint);
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
    vtxs = (rcVertex*)rcAlloc(sizeof(rcVertex) * allocLen, allocHint);
    if (vtxs == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'vtxs' (%d)", allocLen);
        goto Exit0;
    }
    memset(vtxs, 0, sizeof(rcVertex) * allocLen);

    allocLen = nverts * nverts;
    adjMatrix = (Weight*)rcAlloc(sizeof(Weight) * allocLen, allocHint);
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
    CoarsenData intermediateData;
    rcHNAConfig conf;

    conf.condK = 10;
    conf.condR = 0.8f;
    conf.gggpTimes = 4;

    retCode = coarsening(ctx, graph, conf, intermediateData);
    if (!retCode)
    {
        ctx->log(RC_LOG_ERROR, "partitionGraph: coarsening phase fails.");
        goto Exit0;
    }

    retCode = initialPartition(ctx, graph, conf, intermediateData);
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


bool initialPartition(rcContext* ctx, const rcGraphHNA& graph, const rcHNAConfig& conf,
                      CoarsenData& intermediateData)
{
    bool result = false;
    bool retCode = false;
    int allocLen = 0;
    Partition p = nullptr;
    Weight minEdgeCut = 0xffff;
    int level = intermediateData.nlevel;

    allocLen = graph.nvt;
    p = (int*)rcAlloc(sizeof(int) * allocLen, RC_ALLOC_TEMP);
    if (p == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "initialPartition: Out of memory 'p' (%d)", allocLen);
        goto Exit0;
    }
    
    for (int i = 0, n = conf.gggpTimes; i < n; i++)
    {
        rcScopedDelete<int> tPartition((Partition)rcAlloc(sizeof(int) * graph.nvt, RC_ALLOC_TEMP));
        const rcGraphHNA& g = *(intermediateData.levelGraphs[level - 1]);
        retCode = greedyGraphGrowingPartition(g, tPartition);
        if (!retCode)
        {
            ctx->log(RC_LOG_ERROR, "initialPartition: function exec fails 'greedyGraphGrowingPartition'");
            goto Exit0;
        }

        Weight w = calcEdgeCut(graph, tPartition);
        if (minEdgeCut > w)
        {
            minEdgeCut = w;
            memcpy(p, tPartition, sizeof(int) * graph.nvt);
        }
    }

    {
        rcGraphHNA* coarsestGraph = intermediateData.levelGraphs[intermediateData.nlevel - 1];
        rcAssert(coarsestGraph);
        retCode = refinePartition(*coarsestGraph, p);
        if (!retCode)
        {
            ctx->log(RC_LOG_ERROR, "initialPartition: function exec fails 'refinePartition'");
        }
    }

    result = true;
Exit0:
    if (!result && p != nullptr)
    {
        rcFree(p);
        p = nullptr;
    }
    return result;
}


bool coarsening(rcContext* ctx, const rcGraphHNA& graph, const rcHNAConfig& conf,
                CoarsenData& intermediateData)
{
    bool result = false;
    bool retCode = false;

    int level = 0;
    float radio = 1.0f;
    do
    {
        int allocLen = 0;
        const rcGraphHNA* curGraph = level == 0 ? (&graph) : intermediateData.levelGraphs[level - 1];
        rcGraphHNA* coarserGraph = rcAllocGraph(RC_ALLOC_TEMP);
        if (coarserGraph == nullptr)
        {
            ctx->log(RC_LOG_ERROR, "coarsening: Out of memory 'crtGraph', level=%d", level);
            goto Exit0;
        }

        allocLen = curGraph->nvt;
        Match match = (int*)rcAlloc(sizeof(int) * allocLen, RC_ALLOC_TEMP);
        if (match == nullptr)
        {
            ctx->log(RC_LOG_ERROR, "coarsening: Out of memory 'match' (%d)", allocLen);
            goto Exit0;
        }
        memset((int*)match, RC_INVALID_VERTEX, sizeof(match));

        allocLen = curGraph->nvt;
        Partition p = (int*)rcAlloc(sizeof(int) * allocLen, RC_ALLOC_TEMP);
        if (p == nullptr)
        {
            ctx->log(RC_LOG_ERROR, "coarsening: Out of memory 'p' (%d)", allocLen);
            goto Exit0;
        }
        memset((int*)p, RC_INVALID_VERTEX, sizeof(p));

        retCode = coarsenOnce(ctx, *curGraph, *coarserGraph, p, match);
        if (!retCode)
        {
            ctx->log(RC_LOG_ERROR, "'coarsenOnce' fails");
            goto Exit0;
        }

        intermediateData.matches[level] = match;
        intermediateData.partitions[level] = p;
        intermediateData.levelGraphs[level] = coarserGraph;
        level++;
        radio = ((float)coarserGraph->nvt) / graph.nvt;
    } while (level < MAX_COARSEN_LEVEL && radio > conf.condR && intermediateData.levelGraphs[level]->nvt > conf.condK);

    intermediateData.nlevel = level;

    result = true;
Exit0:
    return result;
}


bool coarsenOnce(rcContext* ctx, const rcGraphHNA& curGraph, rcGraphHNA& retCoarserGraph,
                 Partition& retPartition, Match& retMatch)
{
    bool result = false;
    bool retCode = false;
    int newVertNum = 0;

    retCode = heavyEdgeMatch(ctx, curGraph, retMatch, retPartition, newVertNum);
    if (!retCode)
    {
        ctx->log(RC_LOG_ERROR, "coarsenOnce: 'heavyEdgeMatch()' fails");
        goto Exit0;
    }

    rcBuildGraphHNA(ctx, retCoarserGraph, newVertNum, RC_ALLOC_TEMP);
    for (int i = 0, n = curGraph.nvt; i < n; i++)
    {
        rcVertex& v = curGraph.vtxs[i];
        int iu1 = retPartition[i];
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

            int iu2 = retPartition[j];
            if (iu2 != iu1)
            {
                // inter edge
                retCoarserGraph.adjncy[iu1 * curGraph.nvt + iu2] += ewgt;
                retCoarserGraph.adjncy[iu2 * curGraph.nvt + iu1] += ewgt;
                u.adjwgt += ewgt;

            }
            else
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


bool heavyEdgeMatch(rcContext* ctx, const rcGraphHNA& graph, int* retMatch,
                    Partition retPartition, int& retNewVertNum)
{
    bool result = false;
    int* shuffleVerts = nullptr;
    int allocLen = 0;

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

    allocLen = graph.nvt;
    shuffleVerts = (int*)rcAlloc(sizeof(int) * allocLen, RC_ALLOC_TEMP);
    if (shuffleVerts == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "heavyEdgeMatch: Out of memory 'shuffleVerts' (%d)", allocLen);
        goto Exit0;
    }
    for (int i = 0, n = graph.nvt; i < n; i++)
    {
        shuffleVerts[i] = i;
    }
    shuffle(graph.nvt, shuffleVerts, shuffleVerts);

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
            retPartition[index] = retNewVertNum;
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
        retPartition[index] = retNewVertNum;
        retPartition[heaviestEdge] = retNewVertNum;
        retNewVertNum++;
    }//for

    result = true;
Exit0:
    return result;
}


bool shuffle(const int size, const int* src, int* dest)
{
    /// Fisher–Yates shuffle
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


Partition greedyGraphGrowingPartition(const rcGraphHNA& graph, Partition& retPartition)
{
    bool result = false;
    const int nvt = graph.nvt;
    const int hnvt = nvt >> 1;
    const int ip1 = 0;
    const int ip2 = 1;

    bool visitedFlags[MAX_POLY_NUM];
    int insertedCount = 0;

    memset(retPartition, nvt, sizeof(int) * nvt);

    if (graph.nvt < 2)
    {
        goto Exit1;
    }
    else if (graph.nvt == 2)
    {
        retPartition[0] = ip2;
        goto Exit1;
    }


    //insert a random vertex in new partition
    retPartition[rand() % nvt] = ip2;
    insertedCount++;
    do
    {
        int maxGain = 0;
        int maxGainVert = RC_INVALID_INDEX;

        memset(visitedFlags, false, sizeof(bool) * nvt);
        for (int i = 0, n = nvt; i < n; i++)
        {
            if (retPartition[i] == ip1)
                continue;

            for (int j = 0, m = nvt; j < m; j++)
            {
                int newIdx = graph.adjncy[i*nvt + j];
                if (retPartition[newIdx] == retPartition[i])
                    continue;

                if (visitedFlags[newIdx])
                    continue;

                // select the max gain boundary vertex
                int gain = calcKLGain(graph, i, retPartition);
                if (gain >= maxGain)
                {
                    maxGain = gain;
                    maxGainVert = i;
                }
                visitedFlags[newIdx] = true;
            }//for
        }//for

        if (maxGainVert == RC_INVALID_INDEX)
            break;

        retPartition[maxGainVert] = ip2;
        insertedCount++;
    } while (insertedCount >= hnvt);

Exit1:
    result = true;
    //Exit0:
    if (!result && retPartition != nullptr)
    {
        rcFree(retPartition);
        retPartition = nullptr;
    }
    return retPartition;
}


Weight calcKLGain(const rcGraphHNA& graph, const int iv, const Partition& p)
{
    Weight result = 0;
    const int nvt = graph.nvt;

    for (int i = 0, n = nvt; i < n; i++)
    {
        if (i == iv)
            continue;

        Weight ewgt = graph.adjncy[iv * nvt + i];;
        if (p[iv] != p[i])
        {
            result += ewgt;
        }
        else
        {
            result -= ewgt;
        }
    }//for
    return result;
}


Weight calcEdgeCut(const rcGraphHNA& graph, const Partition& p)
{
    Weight edgeCut = 0;
    const int nvt = graph.nvt;

    for (int i = 0, n = nvt; i < n; i++)
    {
        for (int j = 0, m = nvt; j < m; j++)
        {
            if (j >= i)
                break;

            if (p[i] == p[j])
                continue;

            edgeCut += graph.adjncy[i * nvt + j];
        }
    }

    return edgeCut;
}


bool refinePartition(const rcGraphHNA& graph, Partition& partition)
{
    bool result = false;
    bool retCode = false;
    const int nvt = graph.nvt;
    klGainBucketCell* vertTable = nullptr;
    klGainBucketCell* link = nullptr;

    vertTable = (klGainBucketCell*)rcAlloc(sizeof(klGainBucketCell) * nvt, RC_ALLOC_TEMP);
    if (vertTable == nullptr)
        goto Exit0;

    for (int i = 0, n = nvt; i < n; i++)
    {
        klGainBucketCell* pCell = vertTable + i;
        pCell->iv = i;
        pCell->gain = calcKLGain(graph, i, partition);
        pCell->pre = nullptr;
        pCell->next = nullptr;
        if (link == nullptr)
        {
            link = pCell;
        }
        else
        {
            retCode = insertToGainBucket(link, pCell);
            if (!retCode)
                goto Exit0;
        }
    }

    while (link->gain > 0)
    {
        int k = partition[link->iv] == 0 ? 1 : 0;
        swapVert(graph, partition, vertTable, link->iv, k);
    }

    result = true;
Exit0:
    if (vertTable != nullptr)
    {
        rcFree(vertTable);
        vertTable = nullptr;
    }
    return result;
}


bool swapVert(const rcGraphHNA& graph, Partition& partition, klGainBucketCell* gainTbl,
              const int iv, const int k)
{
    const int nvt = graph.nvt;
    if (partition[iv] != k)
    {
        partition[iv] = k;
    }

    for (int i = 0, n = nvt; i < n; i++)
    {
        if (i == iv || graph.adjncy[iv * nvt + i] == 0)
            continue;

        // update adjacency vertex gains
        klGainBucketCell* pCell = gainTbl + i;
        klGainBucketCell* link = pCell->pre == nullptr ? pCell->next : pCell->pre;
        removeFromGainBucket(pCell);
        rcAssert(pCell != nullptr);
        rcAssert(link != nullptr);
        insertToGainBucket(link, pCell);
    }

    return true;
}


bool insertToGainBucket(klGainBucketCell* link, klGainBucketCell* item)
{
    bool result = false;
    if (link == nullptr || item == nullptr)
        goto Exit0;

    if (link->gain > item->gain)
    {
        // search back
        while (link->gain > item->gain && link->next != nullptr)
        {
            link = link->next;
        }
    }
    else if (link->gain < item->gain)
    {
        // search front
        while (link->gain < item->gain && link->pre != nullptr)
        {
            link = link->pre;
        }
    }

    if (link->gain >= item->gain)
    {
        item->pre = link;
        item->next = link->next;
        link->next = item;
        if (item->next != nullptr)
            item->next->pre = item;
    }
    else if (link->gain < item->gain)
    {
        item->pre = link->pre;
        item->next = link;
        link->pre = item;
        if (item->pre != nullptr)
            item->pre->next = item;
    }

    result = true;
Exit0:
    return result;
}


bool removeFromGainBucket(klGainBucketCell* item)
{
    bool result = false;
    if (item == nullptr)
        goto Exit0;

    if (item->pre != nullptr)
        item->pre->next = item->next;

    if (item->next != nullptr)
        item->next->pre = item->pre;

    item->next = nullptr;
    item->pre = nullptr;

    result = true;
Exit0:
    return result;
}