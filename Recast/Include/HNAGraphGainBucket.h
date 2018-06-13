#ifndef HNAGRAPHGAINBUCKET_H
#define HNAGRAPHGAINBUCKET_H

#include "HNAGraph.h"

struct klGainBucketLink
{
    Index*              indexList;
    Weight              gain;
    klGainBucketLink*   pre;
    klGainBucketLink*   next;

public:
    klGainBucketLink(size_t n);
    ~klGainBucketLink();
};




bool    createKLGainBucketLink(rcContext* ctx, const rcHNAGraph& graph);

#endif
