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


bool rcGraphMerge(rcContext* ctx, rcGraph& g1, rcGraph& g2, rcGraphSet& graphSet)
{
    int level = g1.level;
    if (g1.level >= g2.level || level <= 0)
    {
        ctx->log(RC_LOG_ERROR, "TODO");
        return false;
    }

    rcGraph& topGraph = graphSet.graphs[level];
    int index1;
    int index2;
    for (int i = 0, n = topGraph.nverts; i < n; i++)
    {
        if (index1 >= 0 && index2 >= 0)
            break;

        if (topGraph.verts[i] == g1.id)
            index1 = i;

        if (topGraph.verts[i] == g2.id)
            index2 = i;
    }

    if (index1 == index2)
    {
        ctx->log(RC_LOG_ERROR, "TODO");
        return false;
    }

    const int nverts = topGraph.nverts;
    int newVertsNum = g1.nverts + g2.nverts;

    if (topGraph.edgeMatrix[index1 * nverts + index2] <= 0)
    {
        // These two are not adjacent.
        ctx->log(RC_LOG_ERROR, "TODO");
        return false;
    }

    if (g2.verts == 0)
        return true;

    // handle edge cut
    topGraph.edgeMatrix[index1 * nverts + index2] = 0;
    topGraph.edgeMatrix[index2 * nverts + index1] = 0;
    for (int i = 0, n = nverts; i < n; i++)
    {
        Weight weight = topGraph.edgeMatrix[index2 * nverts + i];
        if (i != index1 && weight > 0)
        {
            topGraph.edgeMatrix[index1 * nverts + i] += weight;
            topGraph.edgeMatrix[i * nverts + index1] += weight;
        }
    }
    //handle intra edge
    // TODO need to optimize
    {
        int length = newVertsNum * newVertsNum;
        rcScopedDelete<Weight> tempEdgeMatrix((Weight*)rcAlloc(sizeof(Weight) * length, RC_ALLOC_TEMP));
        for (int i = 0, n = g1.nverts; i < n; i++)
        {
            for (int j = 0, m = g1.nverts; j < m; j++)
            {
                tempEdgeMatrix[i*length + j] = g1.edgeMatrix[i*g1.nverts + j];
                tempEdgeMatrix[j*length + i] = g1.edgeMatrix[j*g1.nverts + i];
            }
        }
        memcpy(g1.edgeMatrix, tempEdgeMatrix, sizeof(tempEdgeMatrix));
    }

    rcGraph& lowerTopGraph = graphSet.graphs[level - 1];
    for (int i = 0, n = g2.nverts; i < n; i++)
    {
        const GraphID childID2 = g1.verts[i];
        const int cindex2 = getGraphIndexInLevel(level - 1, childID2, graphSet);

        for (int j = 0, m = g1.nverts; j < m; j++)
        {
            const GraphID childID1 = g1.verts[j];
            const int cindex1 = getGraphIndexInLevel(level - 1, childID1, graphSet);
            const Weight w = lowerTopGraph.edgeMatrix[cindex1 * lowerTopGraph.nverts + cindex2];
            if (w > 0)
            {
                // new intra edge
                g1.edgeMatrix[j * newVertsNum + g1.nverts] = w;
                g1.edgeMatrix[g1.nverts * newVertsNum + j] = w;
            }
        }

        // new vertex
        g1.verts[g1.nverts] = childID2;
        g1.nverts++;
    }

}

int getGraphIndexInLevel(int level, const GraphID graphID, const rcGraphSet& graphSet)
{
    if (level < graphSet.ngraphs)
        return -1;

    rcGraph& topGraph = graphSet.graphs[level];
    int index = -1;
    for (int i = 0, n = topGraph.nverts; i < n; i++)
    {
        if (topGraph.verts[i] == graphID)
        {
            index = i;
            break;
        }
    }

    return index;
}