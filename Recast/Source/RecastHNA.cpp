#include <string.h>
#include <stdlib.h>
#include "Recast.h"
#include "RecastAssert.h"
#include "RecastHNA.h"
#include "HNAGraph.h"

static const int MAX_COARSEN_LEVEL = 64;

struct rcHNAConfig
{
    int     condK;
    float   condR;
    int     gggpTimes;
};

struct CoarsenData
{
    // coarsen phase
    rcHNAGraph* levelGraphs[MAX_COARSEN_LEVEL] = {nullptr};
    Match       matches[MAX_COARSEN_LEVEL] = {0};
    Map   partitions[MAX_COARSEN_LEVEL] = {0};
    int         nlevel;

    // init partition level
    rcHNAGraph* step2Graph = nullptr;
    Map   step2Ptt = nullptr;
};

struct klGainBucketCell
{
    int iv;
    Weight gain;
    klGainBucketCell* pre;
    klGainBucketCell* next;
};


//////////////////////////////////////////////////////////////////////////
bool partitionGraph(rcContext* ctx, const rcHNAGraph& graph);

bool coarsening(rcContext* ctx, const rcHNAGraph& graph, const rcHNAConfig& conf, CoarsenData& intermediateData);

bool initialPartition(rcContext* ctx, const rcHNAConfig& conf, CoarsenData& intermediateData);

bool uncoarseningPhase();
//////////////////////////////////////////////////////////////////////////

bool greedyGraphGrowingPartition(rcContext* ctx, const rcHNAGraph& graph,
                                 const Map& initPtt, Map& retPtt,
                                 const int iSrcIdex, const int iNewPtt);

//////////////////////////////////////////////////////////////////////////
bool coarsenOnce(rcContext* ctx, const rcHNAGraph& curGraph, rcHNAGraph& retCoarserGraph, Map& retPartition, Match& retMatch);

Weight calcKLGain(const rcHNAGraph& graph, const int iv, const int targetPartition, const Map& p);
Weight calcEdgeCut(const rcHNAGraph& graph, const Map& p);
//////////////////////////////////////////////////////////////////////////


bool refinePartition(const rcHNAGraph& graph, Map& partition, const int p1, const int p2);
bool swapVert(const rcHNAGraph& graph, Map& partition, klGainBucketCell* gainTbl, const int iv, const int k);
bool updateCellGain(klGainBucketCell* pCell, Weight gain);
bool insertToGainBucket(klGainBucketCell* link, klGainBucketCell* item);
bool removeFromGainBucket(klGainBucketCell* item);

klGainBucketCell* getFirstCell(klGainBucketCell* pCell);

bool bisectGraph(rcContext* ctx, const rcHNAConfig& conf,
                 const rcHNAGraph& graph, Map& ptt,
                 const int iTargetPtt, const int iNewPtt);

bool projectToGraph(rcContext* ctx, const Map& ptt, const rcHNAGraph& curGraph,
                    rcHNAGraph& retGraph);



bool partitionGraph(rcContext* ctx, const rcHNAGraph& graph)
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

    retCode = initialPartition(ctx, conf, intermediateData);
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
    //free intermediateData!
    return result;
}


bool initialPartition(rcContext* ctx, const rcHNAConfig& conf, CoarsenData& intermediateData)
{
    bool result = false;
    bool retCode = false;
    int allocLen = 0;
    const int nlevel = intermediateData.nlevel;

    rcAssert(intermediateData.levelGraphs);
    rcAssert(nlevel > 0);
    rcAssert(intermediateData.step2Ptt == nullptr);
    rcAssert(intermediateData.step2Graph == nullptr);

    if (intermediateData.levelGraphs == nullptr || nlevel <= 0)
        goto Exit0;


    //k-way partition
    {
        const rcHNAGraph& coarsestGraph = *(intermediateData.levelGraphs[nlevel - 1]);
        int nParts = nParts = coarsestGraph.nvt;

        allocLen = coarsestGraph.nvt;
        intermediateData.step2Ptt = (int*)rcAlloc(sizeof(int) * allocLen, RC_ALLOC_TEMP);
        if (intermediateData.step2Ptt == nullptr)
        {
            ctx->log(RC_LOG_ERROR, "initialPartition: Out of memory 'intermediateData.step2Ptt' (%d)", allocLen);
            goto Exit0;
        }

        while (nParts < conf.condK)
        {
            for (int i = 0, n = nParts; i < n; i++)
            {
                retCode = bisectGraph(ctx, conf, coarsestGraph, intermediateData.step2Ptt, i, nParts);
                if (!retCode)
                    goto Exit0;

                nParts++;
                if (nParts >= conf.condK)
                    break;
            }//for
        }//while
    }

    //build graph
    {
        rcHNAGraph* step2Graph = nullptr;
        Map& step2Patt = intermediateData.step2Ptt;
        const rcHNAGraph& coarsestGraph = *(intermediateData.levelGraphs[nlevel - 1]);
        intermediateData.step2Graph = rcAllocHNAGraph(RC_ALLOC_TEMP);
        retCode = rcBuildGraphHNA(ctx, *intermediateData.step2Graph, conf.condK, RC_ALLOC_TEMP);
        if (!retCode)
        {
            ctx->log(RC_LOG_ERROR, "initialPartition: exec function fails. 'rcBuildGraphHNA'");
            goto Exit0;
        }

        step2Graph = intermediateData.step2Graph;
        retCode = projectToGraph(ctx, step2Patt, coarsestGraph, *step2Graph);
        if (!retCode)
        {
            ctx->log(RC_LOG_ERROR, "initialPartition: exec function fails. 'projectToGraph'");
            goto Exit0;
        }
    }

    result = true;
Exit0:
    return result;
}


bool coarsening(rcContext* ctx, const rcHNAGraph& graph, const rcHNAConfig& conf,
                CoarsenData& intermediateData)
{
    bool result = false;
    bool retCode = false;

    int level = 0;
    float radio = 1.0f;
    do
    {
        int allocLen = 0;
        const rcHNAGraph* curGraph = level == 0 ? (&graph) : intermediateData.levelGraphs[level - 1];
        rcHNAGraph* coarserGraph = rcAllocHNAGraph(RC_ALLOC_TEMP);
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
        Map p = (int*)rcAlloc(sizeof(int) * allocLen, RC_ALLOC_TEMP);
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
    } while (level < MAX_COARSEN_LEVEL && radio > conf.condR && intermediateData.levelGraphs[level]->nvt > (size_t)conf.condK);

    intermediateData.nlevel = level;

    result = true;
Exit0:
    return result;
}


bool coarsenOnce(rcContext* ctx, const rcHNAGraph& curGraph, rcHNAGraph& retCoarserGraph,
                 Map& retPartition, Match& retMatch)
{
    bool result = false;
    bool retCode = false;
    size_t newVertNum = 0;

    retCode = heavyEdgeMatch(ctx, curGraph, retMatch, retPartition, newVertNum);
    if (!retCode)
    {
        ctx->log(RC_LOG_ERROR, "coarsenOnce: 'heavyEdgeMatch()' fails");
        goto Exit0;
    }

    retCode = rcBuildGraphHNA(ctx, retCoarserGraph, newVertNum, RC_ALLOC_TEMP);
    if (!retCode)
    {
        ctx->log(RC_LOG_ERROR, "coarsenOnce: exec function fails. 'rcBuildGraphHNA'");
        goto Exit0;
    }

    retCode = projectToGraph(ctx, retPartition, curGraph, retCoarserGraph);
    if (!retCode)
    {
        ctx->log(RC_LOG_ERROR, "coarsenOnce: exec function fails. 'projectToGraph'");
        goto Exit0;
    }

    result = true;
Exit0:
    return result;
}



/// Bisect a graph
bool bisectGraph(rcContext* ctx, const rcHNAConfig& conf,
                 const rcHNAGraph& graph, Map& ptt,
                 const int iTargetPtt, const int iNewPtt)
{
    bool result = false;
    bool retCode = false;
    Weight minEdgeCut = 0xFFFF;
    rcScopedDelete<int> scopedVar((int*)rcAlloc(sizeof(int) * graph.nvt, RC_ALLOC_TEMP));
    Map tempPtt = (Map)scopedVar;
    for (int j = 0, m = conf.gggpTimes; j < m; j++)
    {
        Weight edgeCut = 0;

        memset((int*)scopedVar, RC_INVALID_VERTEX, sizeof(int) * graph.nvt);
        retCode = greedyGraphGrowingPartition(ctx, graph, ptt, tempPtt, iTargetPtt, iNewPtt + 1);
        if (!retCode)
        {
            ctx->log(RC_LOG_ERROR, "bisectGraph: function exec fails 'greedyGraphGrowingPartition'");
            goto Exit0;
        }

        edgeCut = calcEdgeCut(graph, tempPtt);
        if (minEdgeCut > edgeCut)
        {
            minEdgeCut = edgeCut;
            memcpy(ptt, scopedVar, sizeof(int) * graph.nvt);
        }
    }

    retCode = refinePartition(graph, ptt, iTargetPtt, iNewPtt);
    if (!retCode)
    {
        ctx->log(RC_LOG_ERROR, "bisectGraph: function exec fails 'refinePartition'");
        goto Exit0;
    }

    result = true;
Exit0:
    return result;
}


bool greedyGraphGrowingPartition(rcContext* ctx, const rcHNAGraph& graph,
                                 const Map& initPtt, Map& retPtt,
                                 const int iSrcIdex, const int iNewPtt)
{
    bool result = false;
    int insertedCount = 0;
    int map[MAX_POLY_NUM];
    int nMap = 0;
    int hNMap = 0;
    klGainBucketCell* gainTbl = nullptr;
    klGainBucketCell* pLink = nullptr;

    memcpy(retPtt, initPtt, sizeof(int) * graph.nvt);
    for (int i = 0, n = graph.nvt; i < n; i++)
    {
        if (retPtt[i] == iSrcIdex)
        {
            map[nMap] = i;
            nMap++;
        }
    }
    hNMap = nMap >> 1;

    if (nMap < 2)
        goto Exit1;

    gainTbl = (klGainBucketCell*)rcAlloc(sizeof(klGainBucketCell) * nMap, RC_ALLOC_TEMP);
    if (gainTbl == nullptr)
        goto Exit0;

    //build gain table
    for (int i = 0, n = nMap; i < n; i++)
    {
        int iv = map[i];
        klGainBucketCell* pCell = gainTbl + i;
        pCell->iv = iv;
        pCell->gain = 0;
        pCell->pre = nullptr;
        pCell->next = nullptr;
    }

    {
        // select the first vertex randomly
        int iv = map[rand() % nMap];
        retPtt[iv] = iNewPtt;
        // update adjacency vertex gains.
        for (int i = 0, n = nMap; i < n; i++)
        {
            int iu = map[i];
            int edge = graph.adjncy[iv * graph.nvt + iu];
            if (edge > 0 && retPtt[iu] == iSrcIdex)
            {
                Weight gain = calcKLGain(graph, iu, iNewPtt, retPtt);
                klGainBucketCell* pCell = gainTbl + i;
                pCell->gain = gain;
                if (pLink == nullptr)
                {
                    pLink = pCell;
                }
                else
                {
                    insertToGainBucket(pLink, pCell);
                    pLink = getFirstCell(pLink);
                }
            }
        }
        insertedCount++;
    }

    rcAssert(pLink != nullptr);
    if (pLink == nullptr)
        goto Exit0;

    while (insertedCount < hNMap && pLink != nullptr)
    {
        int iv = pLink->iv;
        rcAssert(retPtt[iv] == iSrcIdex);
        retPtt[iv] = iNewPtt;
        insertedCount++;

        // remove inserted vertex
        if (pLink->next == nullptr)
        {
            ctx->log(RC_LOG_WARNING, "greedyGraphGrowingPartition: the graph is not continuous");
            break;
        }

        pLink = pLink->next;
        removeFromGainBucket(pLink->pre);

        // update adjacency vertex gains.
        for (int i = 0, n = nMap; i < n; i++)
        {
            int iu = map[i];
            int edge = graph.adjncy[iv * graph.nvt + iu];
            if (edge > 0 && retPtt[iu] == iSrcIdex)
            {
                Weight gain = calcKLGain(graph, iu, iNewPtt, retPtt);
                klGainBucketCell* pCell = gainTbl + i;
                updateCellGain(pCell, gain);
            }
        }

        pLink = getFirstCell(pLink);
    }

Exit1:
    result = true;
Exit0:
    if (gainTbl != nullptr)
    {
        rcFree(gainTbl);
        gainTbl = nullptr;
    }

    return result;
}


Weight calcKLGain(const rcHNAGraph& graph, const int iv,
                  const int targetPtt, const Map& partition)
{
    Weight result = 0;
    const int nvt = graph.nvt;
    const int oriPtt = partition[iv];

    rcAssert(oriPtt != targetPtt);
    if (oriPtt == targetPtt)
        goto Exit0;

    for (int i = 0, n = nvt; i < n; i++)
    {
        if (i == iv)
            continue;

        if (partition[i] == oriPtt)
        {
            result -= graph.adjncy[iv * nvt + i];
        }
        else if (partition[i] == targetPtt)
        {
            result += graph.adjncy[iv * nvt + i];
        }
    }//for

Exit0:
    return result;
}


Weight calcEdgeCut(const rcHNAGraph& graph, const Map& ptt)
{
    Weight edgeCut = 0;
    const int nvt = graph.nvt;

    for (int i = 0, n = nvt; i < n; i++)
    {
        for (int j = 0, m = nvt; j < m; j++)
        {
            if (j >= i)
                break;

            if (ptt[i] == ptt[j])
                continue;

            edgeCut += graph.adjncy[i * nvt + j];
        }
    }

    return edgeCut;
}


bool refinePartition(const rcHNAGraph& graph, Map& partition, const int ptt1, const int ptt2)
{
    bool result = false;
    bool retCode = false;
    const int nvt = graph.nvt;
    klGainBucketCell* gainTbl = nullptr;
    klGainBucketCell* pLink = nullptr;
    int map[MAX_POLY_NUM];
    int nMap = 0;

    memset(map, 0, sizeof(int) * MAX_POLY_NUM);
    for (int i = 0, n = graph.nvt; i < n; i++)
    {
        int p = partition[i];
        if (p == ptt1 || p == ptt2)
        {
            map[nMap] = i;
            nMap++;
        }
    }

    gainTbl = (klGainBucketCell*)rcAlloc(sizeof(klGainBucketCell) * nMap, RC_ALLOC_TEMP);
    if (gainTbl == nullptr)
        goto Exit0;

    for (int i = 0, n = nMap; i < n; i++)
    {
        int iv = map[i];
        int descPtt = partition[iv] == ptt1 ? ptt2 : ptt1;
        klGainBucketCell* pCell = gainTbl + i;
        pCell->iv = iv;
        pCell->gain = calcKLGain(graph, iv, descPtt, partition);
        pCell->pre = nullptr;
        pCell->next = nullptr;
        if (pLink == nullptr)
        {
            pLink = pCell;
        }
        else
        {
            retCode = insertToGainBucket(pLink, pCell);
            if (!retCode)
                goto Exit0;
        }
    }

    while (pLink->gain > 0)
    {
        int descP = partition[pLink->iv] == ptt1 ? ptt2 : ptt1;
        swapVert(graph, partition, gainTbl, pLink->iv, descP);
        pLink = getFirstCell(pLink);
    }

    result = true;
Exit0:
    if (gainTbl != nullptr)
    {
        rcFree(gainTbl);
        gainTbl = nullptr;
    }
    return result;
}


bool swapVert(const rcHNAGraph& graph, Map& ptt, klGainBucketCell* gainTbl,
              const int iv, const int destPtt)
{
    const int nvt = graph.nvt;
    const int ptt1 = ptt[iv];
    if (ptt1 != destPtt)
        ptt[iv] = destPtt;

    for (int i = 0, n = nvt; i < n; i++)
    {
        klGainBucketCell* pCell = nullptr;
        int ptt2 = 0;
        Weight gain = 0;
        if (i == iv || graph.adjncy[iv * nvt + i] == 0)
            continue;

        // update adjacency vertex gains
        pCell = gainTbl + i;
        ptt2 = ptt[i] == ptt1 ? destPtt : ptt1;
        gain = calcKLGain(graph, i, ptt2, ptt);
        updateCellGain(pCell, gain);
    }

    return true;
}


bool updateCellGain(klGainBucketCell* pCell, Weight gain)
{
    if (pCell == nullptr)
        return false;

    if (pCell->gain == gain)
        return true;

    klGainBucketCell* link = pCell->pre != nullptr ? pCell->pre : pCell->next;
    if (link == nullptr)
        return true;

    removeFromGainBucket(pCell);
    insertToGainBucket(link, pCell);
    return true;
}


bool insertToGainBucket(klGainBucketCell* pLink, klGainBucketCell* pCell)
{
    bool result = false;
    if (pLink == nullptr || pCell == nullptr)
        goto Exit0;

    if (pLink->gain > pCell->gain)
    {
        // search back
        while (pLink->gain > pCell->gain && pLink->next != nullptr)
        {
            pLink = pLink->next;
        }
    }
    else if (pLink->gain < pCell->gain)
    {
        // search front
        while (pLink->gain < pCell->gain && pLink->pre != nullptr)
        {
            pLink = pLink->pre;
        }
    }

    if (pLink->gain >= pCell->gain)
    {
        pCell->pre = pLink;
        pCell->next = pLink->next;
        pLink->next = pCell;
        if (pCell->next != nullptr)
            pCell->next->pre = pCell;
    }
    else if (pLink->gain < pCell->gain)
    {
        pCell->pre = pLink->pre;
        pCell->next = pLink;
        pLink->pre = pCell;
        if (pCell->pre != nullptr)
            pCell->pre->next = pCell;
    }

    result = true;
Exit0:
    return result;
}


bool removeFromGainBucket(klGainBucketCell* pCell)
{
    bool result = false;
    if (pCell == nullptr)
        goto Exit0;

    if (pCell->pre != nullptr)
        pCell->pre->next = pCell->next;

    if (pCell->next != nullptr)
        pCell->next->pre = pCell->pre;

    pCell->next = nullptr;
    pCell->pre = nullptr;

    result = true;
Exit0:
    return result;
}


klGainBucketCell* getFirstCell(klGainBucketCell* pCell)
{
    if (pCell == nullptr)
        return nullptr;

    while (pCell->pre != nullptr)
        pCell = pCell->pre;

    return pCell;
}


bool projectToGraph(rcContext* ctx, const Map& ptt, const rcHNAGraph& curGraph,
                    rcHNAGraph& retGraph)
{
    bool result = false;
    for (int i = 0, n = curGraph.nvt; i < n; i++)
    {
        rcHNAVertex* v = curGraph.vtxs + i;
        rcHNAVertex* u = nullptr;
        int iu1 = ptt[i];

        if ((size_t)iu1 >= retGraph.nvt)
        {
            ctx->log(RC_LOG_ERROR, "generateNewGraph: Out of length. 'retGraph'");
            goto Exit0;
        }

        u = retGraph.vtxs + iu1;
        u->ipoly = RC_MESH_NULL_IDX;
        u->vwgt += v->vwgt;
        u->iedges = 0; // set in next step.

        for (int j = 0, m = curGraph.nvt; j < m; j++)
        {
            if (i == j)
                continue;

            Weight ewgt = curGraph.adjncy[v->iedges + j];
            if (ewgt == 0)
                continue;

            int iu2 = ptt[j];
            if (iu2 != iu1)
            {
                // inter edge
                retGraph.adjncy[iu1 * curGraph.nvt + iu2] += ewgt;
                retGraph.adjncy[iu2 * curGraph.nvt + iu1] += ewgt;
                u->adjwgt += ewgt;

            }
            else
            {
                // intra edge
                u->cewgt += ewgt;
            }
        }//for
    }//for

    for (int i = 0, n = retGraph.nvt; i < n; i++)
    {
        rcHNAVertex& u = retGraph.vtxs[i];
        for (int j = 0, m = retGraph.nvt; j < m; j++)
        {
            if (i != j && retGraph.adjncy[u.iedges + j] > 0)
            {
                u.nedges += 1;
            }
        }
    }

    result = true;
Exit0:
    return result;
}