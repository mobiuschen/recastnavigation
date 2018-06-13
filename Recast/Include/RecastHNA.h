#ifndef RECASTHNA_H
#define RECASTHNA_H

#include "RecastAlloc.h"
#include "HNAGraph.h"
#include <vector>

typedef unsigned short GraphID;

class rcContext;

struct rcPolyMesh;

static const unsigned int RC_INVALID_VERTEX = 0xFFFF;
static const int MAX_POLY_NUM = 100000;

struct rcHNAMatch
{
    Index   match[MAX_POLY_NUM];
    size_t  size;
};


struct rcHNAMap
{
    Index   map[MAX_POLY_NUM];
    size_t  size;
};

struct klGainBucketLink
{
    std::vector<Index>  iv;
    Weight              gain;
    klGainBucketLink*   pre;
    klGainBucketLink*   next;
};


#endif