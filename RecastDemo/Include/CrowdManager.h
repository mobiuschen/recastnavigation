//
// Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

#ifndef CROWDMANAGER_H
#define CROWDMANAGER_H

#include "DetourNavMeshQuery.h"
#include "DetourObstacleAvoidance.h"
#include "ValueHistory.h"


class ProximityGrid
{
	int m_maxItems;
	float m_cellSize;
	float m_invCellSize;

	struct Item
	{
		unsigned short id;
		short x,y;
		unsigned short next;
	};
	Item* m_pool;
	int m_poolHead;
	int m_poolSize;
	
	unsigned short* m_buckets;
	int m_bucketsSize;
	
	int m_bounds[4];
	
public:
	ProximityGrid();
	~ProximityGrid();
	
	bool init(const int maxItems, const float cellSize);
	
	void clear();
	
	void addItem(const unsigned short id,
				 const float minx, const float miny,
				 const float maxx, const float maxy);
	
	int queryItems(const float minx, const float miny,
				   const float maxx, const float maxy,
				   unsigned short* ids, const int maxIds) const;
	
	int getItemCountAt(const int x, const int y) const;
	const int* getBounds() const { return m_bounds; }
	const float getCellSize() const { return m_cellSize; }
};



static const int AGENT_MAX_PATH = 256;
static const int AGENT_MAX_CORNERS = 4;
static const int AGENT_MAX_TRAIL = 64;
static const int AGENT_MAX_LOCALSEGS = 32;
static const int AGENT_MAX_NEIS = 8;

static const unsigned int PATHQ_INVALID = 0;

enum PathQueueRequestState
{
	PATHQ_STATE_INVALID,
	PATHQ_STATE_WORKING,
	PATHQ_STATE_READY,
};

typedef unsigned int PathQueueRef;

class PathQueue
{
	struct PathQuery
	{
		// Path find start and end location.
		float startPos[3], endPos[3];
		dtPolyRef startRef, endRef;
		// Result.
		dtPolyRef path[AGENT_MAX_PATH];
		bool ready;
		int npath;
		PathQueueRef ref;
		const dtQueryFilter* filter; // TODO: This is potentially dangerous!
		int keepalive;
	};
	
	static const int MAX_QUEUE = 8;
	PathQuery m_queue[MAX_QUEUE];
	PathQueueRef m_nextHandle;
	
	int m_delay;
	
public:
	PathQueue();
	~PathQueue();
	
	void update(dtNavMeshQuery* navquery);
	PathQueueRef request(dtPolyRef startRef, dtPolyRef endRef,
						 const float* startPos, const float* endPos, 
						 const dtQueryFilter* filter);
	int getRequestState(PathQueueRef ref);
	int getPathResult(PathQueueRef ref, dtPolyRef* path, const int maxPath);
};

// TODO: Remove dynamics stuff from mover (velocity, new pos).
// Mover should just store current location, target location and
// path corridor, corners and collision segments.

class Mover
{
	float m_pos[3];
	float m_target[3];
	float m_radius, m_height;
	
	float m_dvel[3];
	float m_nvel[3];
	float m_vel[3];
	
	float m_npos[3];
	
	float m_pathOptimizationRange;
	float m_collisionQueryRange;

	float m_localCenter[3];
	float m_localSegs[AGENT_MAX_LOCALSEGS*6];
	int m_localSegCount;

	dtPolyRef m_path[AGENT_MAX_PATH];
	int m_npath;
	
	float m_cornerVerts[AGENT_MAX_CORNERS*3];
	unsigned char m_cornerFlags[AGENT_MAX_CORNERS];
	dtPolyRef m_cornerPolys[AGENT_MAX_CORNERS];
	int m_ncorners;
	
public:
	Mover();
	~Mover();
	
	void init(dtPolyRef ref, const float* pos, const float radius, const float height,
			  const float collisionQueryRange, const float pathOptimizationRange);
	
	void updateLocalNeighbourhood(dtNavMeshQuery* navquery, const dtQueryFilter* filter);
	void updateCorners(dtNavMeshQuery* navquery, const dtQueryFilter* filter, float* opts = 0, float* opte = 0);
	void integrate(const float maxAcc, const float dt);
	void updatePosition(dtNavMeshQuery* navquery, const dtQueryFilter* filter);
	float getDistanceToGoal(const float range) const;
	void calcSmoothSteerDirection(float* dvel);
	void calcStraightSteerDirection(float* dvel);
	void appendLocalCollisionSegments(dtObstacleAvoidanceQuery* obstacleQuery);

	void setDesiredVelocity(const float* dvel);
	void setNewVelocity(const float* nvel);
	void setNewPos(const float* npos);

	inline int getCornerCount() const { return m_ncorners; }
	inline const float* getCornerPos(int i) const { return &m_cornerVerts[i*3]; }
	
	inline const float* getLocalCenter() const { return m_localCenter; }
	inline int getLocalSegmentCount() const { return m_localSegCount; }
	inline const float* getLocalSegment(int i) const { return &m_localSegs[i*6]; }
	
	inline const float* getPos() const { return m_pos; }
	inline const float* getNewPos() const { return m_npos; }
	inline const float* getVelocity() const { return m_vel; }
	inline const float* getDesiredVelocity() const { return m_dvel; }
	inline float getRadius() const { return m_radius; }
	inline float getHeight() const { return m_height; }
	inline float getCollisionQueryRange() const { return m_collisionQueryRange; }
	
	inline void setCorridor(const float* target, const dtPolyRef* path, int npath);
	inline const float* getCorridorTarget() const { return m_target; }
	inline const dtPolyRef* getCorridor() const { return m_path; }
	inline int getCorridorCount() const { return m_npath; } 	
};

struct Agent
{
	unsigned char active;
	
	Mover mover;
	
	float maxspeed;
	float t;
	float var;
	
	float opts[3], opte[3];
	float disp[3];
	
	float trail[AGENT_MAX_TRAIL*3];
	int htrail;
};


enum UpdateFlags
{
	CROWDMAN_ANTICIPATE_TURNS = 1,
	CROWDMAN_USE_VO = 2,
	CROWDMAN_DRUNK = 4,
};


class CrowdManager
{
	static const int MAX_AGENTS = 32;
	Agent m_agents[MAX_AGENTS];
	dtObstacleAvoidanceDebugData* m_vodebug[MAX_AGENTS];
	dtObstacleAvoidanceQuery* m_obstacleQuery;
	PathQueue m_pathq;
	ProximityGrid m_grid;
	
	float m_ext[3];
	dtQueryFilter m_filter;

	int m_totalTime;
	int m_rvoTime;
	int m_sampleCount;

	enum MoveRequestState
	{
		MR_TARGET_REQUESTING,
		MR_TARGET_WAITING_FOR_PATH,
		MR_TARGET_VALID,
		MR_TARGET_FAILED,
	};
	
	struct MoveRequest
	{
		int idx;
		dtPolyRef ref;
		float pos[3];
		unsigned char state;
		PathQueueRef pathqRef;
	};
	MoveRequest m_moveRequests[MAX_AGENTS];
	int m_moveRequestCount;
	
	int getNeighbours(const float* pos, const float height, const float range,
					  const Agent* skip, Agent** result, const int maxResult);

public:
	CrowdManager();
	~CrowdManager();
	
	void reset();
	const Agent* getAgent(const int idx);
	const int getAgentCount() const;
	int addAgent(const float* pos, const float radius, const float height, dtNavMeshQuery* navquery);
	void removeAgent(const int idx);
	bool requestMoveTarget(const int idx, dtPolyRef ref, const float* pos);
	
	int getActiveAgents(Agent** agents, const int maxAgents);
	
	void update(const float dt, unsigned int flags, dtNavMeshQuery* navquery);
	
	const dtQueryFilter* getFilter() const { return &m_filter; }
	const float* getQueryExtents() const { return m_ext; }
	
	const dtObstacleAvoidanceDebugData* getVODebugData(const int idx) const { return m_vodebug[idx]; }	
	inline int getTotalTime() const { return m_totalTime; }
	inline int getRVOTime() const { return m_rvoTime; }
	inline int getSampleCount() const { return m_sampleCount; }
	const ProximityGrid* getGrid() const { return &m_grid; }

};


#endif // CROWDMANAGER_H