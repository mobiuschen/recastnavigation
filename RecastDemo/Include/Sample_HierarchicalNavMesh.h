#ifndef RECASTSAMPLEHIERARCHICALNAVMESH_H
#define RECASTSAMPLEHIERARCHICALNAVMESH_H
#include "Sample_SoloMesh.h"
#include "DetourNavMesh.h"
#include "Recast.h"
#include "RecastGraph.h"


class Sample_HierarchicalNavMesh : public Sample_SoloMesh {

private:
	rcGraph m_graphSet;
public:
	Sample_HierarchicalNavMesh();
	virtual ~Sample_HierarchicalNavMesh();

	virtual bool handleBuild();
protected:
	void cleanup();
};
#endif