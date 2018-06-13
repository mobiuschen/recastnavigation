#include <string.h>
#include "HNAGraphGainBucket.h"
#include "Recast.h"
#include "RecastAlloc.h"


klGainBucketLink::klGainBucketLink(size_t n):
    indexList(nullptr),
    gain(0),
    pre(nullptr),
    next(nullptr)
{
    rcIgnoreUnused(n);
}




klGainBucketLink::~klGainBucketLink()
{
    rcFree(indexList);
    indexList = nullptr;
}

