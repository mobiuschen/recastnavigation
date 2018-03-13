#ifndef RECASTSAMPLEHIERARCHICALNAVMESH_H
#define RECASTSAMPLEHIERARCHICALNAVMESH_H
#include "Sample_SoloMesh.h"
#include "DetourNavMesh.h"
#include "Recast.h"


class Sample_HierarchicalNavMesh : public Sample_SoloMesh {

private:
public:
	Sample_HierarchicalNavMesh();
	virtual ~Sample_HierarchicalNavMesh();

	virtual bool handleBuild();
protected:
	void cleanup();
};
#endif