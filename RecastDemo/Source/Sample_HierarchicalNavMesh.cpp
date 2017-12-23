
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include "InputGeom.h"
#include "Sample_HierarchicalNavMesh.h"

Sample_HierarchicalNavMesh::Sample_HierarchicalNavMesh()
{

}

Sample_HierarchicalNavMesh::~Sample_HierarchicalNavMesh()
{
	cleanup();
}


bool Sample_HierarchicalNavMesh::handleBuild()
{
	return true;
}

void Sample_HierarchicalNavMesh::cleanup()
{

}

