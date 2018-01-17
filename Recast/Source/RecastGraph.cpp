#include <string.h>
#include <stdlib.h>
#include "Recast.h"
#include "RecastAlloc.h"
#include "RecastGraph.h"

inline Weight isAdjVerts(const rcGraph& graph, int uidx, int vidx)
{
    return graph.edgeMatrix[uidx * graph.nverts + vidx];
}

//static int getGraphIndexInLevel(int level, const GraphID graphID, const rcGraphSet& graphSet)
//{
//    if (level < graphSet.ngraphs)
//        return -1;
//
//    rcGraph& topGraph = graphSet.graphs[level];
//    int index = -1;
//    for (int i = 0, n = topGraph.nverts; i < n; i++)
//    {
//        if (topGraph.verts[i] == graphID)
//        {
//            index = i;
//            break;
//        }
//    }
//
//    return index;
//}


static int getGraphIndex(const rcGraph& graph, const GraphID childID)
{
    int index = -1;
    for (int i = 0, n = graph.nverts; i < n; i++)
    {
        if (graph.verts[i] == childID)
        {
            index = i;
            break;
        }
    }

    return index;
}

static int getRandomGraphVert(const rcGraph& graph, bool* travelledAry)
{
    const int nverts = graph.nverts;
    int startIdx = rand() % nverts;
    int selectIdx = startIdx;
    bool findFlag = false;
    for (int i = 0, n = graph.nverts; i < n; i++)
    {
        selectIdx = (startIdx + i) % graph.nverts;
        if (!travelledAry[selectIdx])
        {
            findFlag = true;
            break;
        }
    }//for

    if (!findFlag)
    {
        return -1;
    }

    return selectIdx;
}


static Weight calcEdgeWeight(const rcGraph& g1, const rcGraph& g2, const rcGraph& lowGraph)
{
    //g1 and g2 must not have same sub graph.
    bool result = false;
    Weight totalWeight = 0;
    for (int i = 0, n = g1.nverts; i < n; i++)
    {
        GraphID childID1 = g1.verts[i];
        int childIdx1 = getGraphIndex(lowGraph, childID1);
        if (childIdx1 < 0)
            goto Exit0;

        for (int j = 0, m = g2.nverts; j < m; j++)
        {
            GraphID childID2 = g2.verts[j];
            int childIdx2 = getGraphIndex(lowGraph, childID2);
            if (childIdx2)
                goto Exit0;

            totalWeight += lowGraph.edgeMatrix[childIdx1 * lowGraph.nverts + childIdx2];
        }
    }

    result = true;
Exit0:
    if (!result)
        totalWeight = 0;
    return totalWeight;
}


//static bool mergeGraphs(rcContext* ctx, rcGraph& dest, const rcGraph& src, rcGraphSet& graphSet)
//{
//    int level = dest.level;
//    if (dest.level >= src.level || level <= 0)
//    {
//        ctx->log(RC_LOG_ERROR, "TODO");
//        return false;
//    }
//
//    rcGraph& topGraph = graphSet.graphs[level];
//    int index1;
//    int index2;
//    for (int i = 0, n = topGraph.nverts; i < n; i++)
//    {
//        if (index1 >= 0 && index2 >= 0)
//            break;
//
//        if (topGraph.verts[i] == dest.id)
//            index1 = i;
//
//        if (topGraph.verts[i] == src.id)
//            index2 = i;
//    }
//
//    if (index1 == index2)
//    {
//        ctx->log(RC_LOG_ERROR, "TODO");
//        return false;
//    }
//
//    const int nverts = topGraph.nverts;
//    if (topGraph.edgeMatrix[index1 * nverts + index2] <= 0)
//    {
//        // These two are not neighbour.
//        ctx->log(RC_LOG_ERROR, "TODO");
//        return false;
//    }
//
//    if (src.verts == 0)
//        return true;
//
//    // handle the edge cuts
//    topGraph.edgeMatrix[index1 * nverts + index2] = 0;
//    topGraph.edgeMatrix[index2 * nverts + index1] = 0;
//    for (int i = 0, n = nverts; i < n; i++)
//    {
//        Weight weight = topGraph.edgeMatrix[index2 * nverts + i];
//        if (i != index1 && weight > 0)
//        {
//            topGraph.edgeMatrix[index1 * nverts + i] += weight;
//            topGraph.edgeMatrix[i * nverts + index1] += weight;
//        }
//    }
//
//    //handle the intra edges and vertcies
//    const int newVertsNum = dest.nverts + src.nverts;
//    // TODO need to optimize
//    {
//        // Expand the size of edge matrix.
//        int length = newVertsNum * newVertsNum;
//        Weight* edgeMatrix = (Weight*)rcAlloc(sizeof(Weight) * length, RC_ALLOC_PERM);
//        for (int i = 0, n = dest.nverts; i < n; i++)
//        {
//            for (int j = 0, m = dest.nverts; j < m; j++)
//            {
//                edgeMatrix[i*length + j] = dest.edgeMatrix[i*dest.nverts + j];
//                edgeMatrix[j*length + i] = dest.edgeMatrix[j*dest.nverts + i];
//            }
//        }
//
//        rcFree(dest.edgeMatrix);
//        dest.edgeMatrix = edgeMatrix;
//    }
//
//    // Need to check the sub-verts adjacent realtion in lower level graph.
//    const rcGraph& lowerTopGraph = graphSet.graphs[level - 1];
//    for (int i = 0, n = src.nverts; i < n; i++)
//    {
//        const GraphID newVertID = src.verts[i];
//        const int newVertIdx = getGraphIndexInLevel(level - 1, newVertID, graphSet);
//
//        for (int j = 0, m = dest.nverts; j < m; j++)
//        {
//            const GraphID oriVertID = dest.verts[j];
//            const int oriVertIdx = getGraphIndexInLevel(level - 1, oriVertID, graphSet);
//            const Weight w = lowerTopGraph.edgeMatrix[oriVertIdx * lowerTopGraph.nverts + newVertIdx];
//            if (w > 0)
//            {
//                // The ori vert is adjacent
//                // new intra edge
//                dest.edgeMatrix[j * newVertsNum + dest.nverts] = w;
//                dest.edgeMatrix[dest.nverts * newVertsNum + j] = w;
//            }
//        }
//
//        // new vertex
//        dest.verts[dest.nverts] = newVertID;
//        dest.nverts++;
//    }
//
//}

static bool coarseningPhase(rcContext* ctx,
                            const rcGraph& graph,
                            int& retNverts)
{
    bool result = false;
    bool retCode = false;
    const unsigned short nverts = graph.nverts;
    const unsigned short poolCapacity = graph.nverts * 10;
    unsigned short poolSize = 0;
    rcScopedDelete<rcGraph> graphPool((rcGraph*)rcAlloc(sizeof(rcGraph) * poolCapacity, RC_ALLOC_TEMP));
    GraphID matchings[MAX_POLY_NUM];
    Weight matchingWeights[MAX_POLY_NUM];
    int nmatching = 0;


    {
        bool matchedFlags[MAX_POLY_NUM];

        for (int i = 0, n = graph.nverts; i < n; i++)
        {
            // select a random vertex u
            int uIdx = getRandomGraphVert(graph, matchedFlags);
            if (uIdx < 0)
                break;

            // select adjacent vertex v
            // Heavy Edge Matching (HEM)
            int startIdx = rand() % graph.nverts;
            int vIdx = startIdx;
            Weight maxWeight = 0;
            for (int j = 0, m = nverts; j < m; j++)
            {
                int idx = (startIdx + j) % nverts;
                //select a max weight neighbor
                Weight w = isAdjVerts(graph, uIdx, idx);
                if (!matchedFlags[idx] && w > maxWeight)
                {
                    vIdx = idx;
                    maxWeight = w;
                }
            }

            if (maxWeight > 0)
            {
                // merge to a new graph
                GraphID gid = RC_GRAPH_ID_NULL;
                rcGraph* child = rcAllocGraph(graphPool, poolSize, gid);
                if (child == nullptr)
                    goto Exit0;

                retCode = rcBuildGraph(ctx, *child, gid, 2, RC_ALLOC_TEMP);
                if (retCode)
                    goto Exit0;

                child->id = gid;
                child->verts[0] = graph.verts[uIdx];
                child->verts[1] = graph.verts[vIdx];
                child->edgeMatrix[1] = maxWeight;
                child->edgeMatrix[2] = maxWeight;
                child->weights[0] = graph.weights[uIdx];
                child->weights[1] = graph.weights[vIdx];
                matchedFlags[uIdx] = true;
                matchedFlags[vIdx] = true;
                matchings[nmatching] = child->id;
                matchingWeights[nmatching] = child->weights[0] + child->weights[1];
                nmatching++;
            }
        }//for

         // unmatch vertex would remain as solo vertex
        for (int i = 0, n = nverts; i < n; i++)
        {
            if (!matchedFlags[i])
            {
                GraphID gid = RC_GRAPH_ID_NULL;
                rcGraph* child = rcAllocGraph(graphPool, poolSize, gid);
                if (child == nullptr)
                    goto Exit0;

                retCode = rcBuildGraph(ctx, *child, gid, 1, RC_ALLOC_TEMP);
                if (retCode)
                    goto Exit0;
                                
                child->id = gid;
                child->verts[0] = graph.verts[i];
                child->weights[0] = graph.weights[i];
                matchedFlags[i] = true;
                matchings[nmatching] = child->id;
                matchingWeights[nmatching] = child->weights[0];
                nmatching++;
            }
        }//for
    }


    {
        GraphID gid = RC_GRAPH_ID_NULL;
        rcGraph* highGraph = rcAllocGraph(graphPool, poolSize, gid);
        const rcGraph* lowGraph = &graph;
        if (highGraph == nullptr || gid == RC_GRAPH_ID_NULL)
            goto Exit0;

        retCode = rcBuildGraph(ctx, *highGraph, gid, nmatching, RC_ALLOC_TEMP);
        if (!retCode)
        {
            ctx->log(RC_LOG_ERROR, "coarseningPhase: build graph fails 'highGraph'");
            goto Exit0;
        }

        for (int i = 0, n = highGraph->nverts; i < n; i++)
        {
            GraphID childID = matchings[i];
            highGraph->verts[i] = childID;
            highGraph->weights[i] = matchingWeights[i];
            highGraph->edgeMatrix[i * nmatching + i] = 0;
        }//for

         //calculate coarser graph edgeMatrix
        for (int i = 0, n = highGraph->nverts; i < n; i++)
        {
            GraphID gid1 = highGraph->verts[i];
            rcGraph& child1 = graphPool[gid1];
            for (int j = 0, m = nmatching; j < m; m++)
            {
                if (i == j)
                    continue;

                GraphID gid2 = highGraph->verts[i];
                rcGraph& child2 = graphPool[gid2];
                Weight edgeWeight = calcEdgeWeight(child1, child2, *lowGraph);
                highGraph->edgeMatrix[i * highGraph->nverts + j] = edgeWeight;
                highGraph->edgeMatrix[j * highGraph->nverts + i] = edgeWeight;
            }//for
        }
    }

    retNverts = nmatching;

    result = true;
Exit0:
    return result;
}

static bool initPartitionPhase(rcContext* ctx, const rcGraphSet& graphSet, const rcGraph& graph, int k)
{
    // multilevel recursive bisection(MLRB) algorithm
    //TODO
    return true;
}

static bool uncoarseningPhase(rcContext* ctx, const rcGraphSet& graphSet, const rcGraph& graph)
{
    //TODO
    return true;
}


static bool partitionGraph(rcContext* ctx, rcGraphSet& graphSet, rcGraph& graph, int nparts, const rcGraph& lowerGraph)
{
    //Multilevel k-way Partitioning Scheme (MLkP)
    //https://glaros.dtc.umn.edu/gkhome/fetch/papers/mlJPDC98.pdf

    //Step1. coarsening phase
    //coarseningPhase(ctx, graphSet, lowerGraph);
    //Step2. initial partitioniong phase
    initPartitionPhase(ctx, graphSet, lowerGraph, nparts);
    //Step3. uncoarsening phase
    uncoarseningPhase(ctx, graphSet, lowerGraph);
    return true;
}

rcGraphSet* rcAllocGraphSet(rcContext* ctx, int npoly, int level)
{
    int length = 0;
    int ngraphs = npoly * 2 + level + 1;
    rcGraphSet* result = (rcGraphSet*)rcAlloc(sizeof(rcGraphSet), RC_ALLOC_PERM);
    memset(result, 0, sizeof(rcGraphSet));


    length = ngraphs;
    result->graphs = (rcGraph*)rcAlloc(sizeof(rcGraph) * length, RC_ALLOC_PERM);
    if (!result->graphs)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'result->graphs' (%d).", length);
        return false;
    }
    memset(result->graphs, 0, sizeof(rcGraph) * length);

    length = level;
    result->topGraphs = (GraphID*)rcAlloc(sizeof(GraphID) * length, RC_ALLOC_PERM);
    if (!result->topGraphs)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'result->topGraphs' (%d).", length);
        return false;
    }
    memset(result->topGraphs, 0, sizeof(GraphID) * length);

    length = ngraphs;
    result->graphs = (rcGraph*)rcAlloc(sizeof(rcGraph) * length, RC_ALLOC_PERM);
    if (!result->graphs)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'result->graphs' (%d).", length);
        return false;
    }
    memset(result->graphs, 0, sizeof(rcGraph) * length);

    return result;
}


bool rcBuildGraphSet(rcContext* ctx, const rcPolyMesh& pmesh, rcGraphSet& graphSet, const int level)
{
    unsigned int length = 0;
    const int npolys = pmesh.npolys;
    const int nvp = pmesh.nvp;
    const int ngraphs = npolys + level + 1;
    graphSet.ngraphs = ngraphs;

    int initedGraphCount = level + 1;
    for (int i = 0, n = level + 1; i < n; i++)
    {
        const rcGraph* lowerGraph = i > 0 ? &(graphSet.graphs[i - 1]) : nullptr;
        rcGraph* crtLevelGraph = &(graphSet.graphs[i]);
        graphSet.topGraphs[i] = (GraphID)i;
        const float factor = 0.1;
        int subGraphNum = i == 0 ? pmesh.npolys : (int)(lowerGraph->nverts * factor);

        length = subGraphNum;
        crtLevelGraph->verts = (GraphID*)rcAlloc(sizeof(GraphID) * length, RC_ALLOC_PERM);
        if (!crtLevelGraph->verts)
        {
            ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'graph->verts' (%d).", length);
            return false;
        }
        memset(crtLevelGraph->verts, RC_MESH_NULL_IDX, sizeof(GraphID) * length);

        length = subGraphNum;
        crtLevelGraph->weights = (Weight*)rcAlloc(sizeof(Weight) * length, RC_ALLOC_PERM);
        if (!crtLevelGraph->weights)
        {
            ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'graph->weights' (%d).", length);
            return false;
        }
        memset(crtLevelGraph->weights, 1, sizeof(Weight) * length);

        length = subGraphNum * subGraphNum;
        crtLevelGraph->edgeMatrix = (Weight*)rcAlloc(sizeof(Weight) * length, RC_ALLOC_PERM);
        if (!crtLevelGraph->edgeMatrix)
        {
            ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'graph->edgeMatrix' (%d).", length);
            return false;
        }
        memset(crtLevelGraph->edgeMatrix, 0, sizeof(Weight) * length);

        if (i == 0)
        {
            //In Level 0, every graph represents a poly.
            for (int j = 0; j < npolys; j++)
            {
                unsigned short* p = &(pmesh.polys[j * nvp * 2]);
                rcGraph& vert = graphSet.graphs[initedGraphCount];
                vert.id = (GraphID)initedGraphCount;
                initedGraphCount++;

                crtLevelGraph->verts[j] = vert.id;
                crtLevelGraph->nverts++;

                // set vertex graph
                vert.verts = nullptr;
                vert.weights = nullptr;
                vert.edgeMatrix = nullptr;
                vert.nverts = 0;
                vert.poly = (unsigned short)j;

                unsigned short* adjs = &(p[nvp]);
                // build edges
                for (int k = 0, l = nvp; k < l; k++)
                {
                    int pAdj = adjs[k];
                    if (p[k] == RC_MESH_NULL_IDX) break;
                    if (pAdj & 0x8000) continue;

                    crtLevelGraph->edgeMatrix[j * npolys + pAdj] = 1;
                    crtLevelGraph->edgeMatrix[pAdj * npolys + j] = 1;
                }//for
            }//for
        }
        else
        {
            int nparts = (int)(lowerGraph->nverts * 0.1);
            partitionGraph(ctx, graphSet, *crtLevelGraph, nparts, *lowerGraph);
        }
    }//for

    return true;
}

void rcFreeGraphSet(rcGraphSet* pData)
{
    if (pData == nullptr)
        return;

    for (int i = 0, n = pData->ngraphs; i < n; i++)
    {
        rcFreeGraph((pData->graphs) + i);
    }
    rcFree(pData->graphs);
    rcFree(pData);
}


rcGraph* rcAllocGraph(rcGraph* graphPool, unsigned short& poolSize, GraphID& retGraphID)
{
    retGraphID = RC_GRAPH_ID_NULL;
    rcGraph* graph = &(graphPool[poolSize]);
    retGraphID = poolSize;
    poolSize++;
    memset(graph, 0, sizeof(rcGraph));
    return graph;
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

bool rcBuildGraph(rcContext* ctx, rcGraph& graph, const GraphID id, const unsigned short nverts, const rcAllocHint hint)
{
    bool result = false;
    int length = 0;
    graph.id = id;
    graph.nverts = nverts;
    graph.poly = 0;

    length = nverts * nverts;
    graph.edgeMatrix = (Weight*)rcAlloc(sizeof(Weight)*length, hint);
    if (graph.edgeMatrix == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'graph.edgeMatrix' (%d).", length);
        goto Exit0;
    }
    memset(graph.edgeMatrix, 0, sizeof(Weight) * length);

    length = nverts;
    graph.weights = (Weight*)rcAlloc(sizeof(Weight)*length, hint);
    if (graph.weights == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'graph.weights' (%d).", length);
        goto Exit0;
    }
    memset(graph.weights, 0, sizeof(Weight) * length);

    length = nverts;
    graph.verts = (GraphID*)rcAlloc(sizeof(GraphID)*length, hint);
    if (graph.verts == nullptr)
    {
        ctx->log(RC_LOG_ERROR, "rcBuildGraph: Out of memory 'graph.verts' (%d).", length);
        goto Exit0;
    }
    memset(graph.verts, RC_GRAPH_ID_NULL, sizeof(GraphID) * length);

    result = true;
Exit0:
    if (!result)
    {
        if (graph.edgeMatrix != nullptr)
        {
            rcFree(graph.edgeMatrix);
            graph.edgeMatrix = nullptr;
        }

        if (graph.weights != nullptr)
        {
            rcFree(graph.weights);
            graph.weights = nullptr;
        }

        if (graph.verts != nullptr)
        {
            rcFree(graph.verts);
            graph.verts = nullptr;
        }
    }
    return result;
}


