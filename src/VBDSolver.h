#pragma once

#include "defines.h"
#include "half-edge-mesh.h"
#include <iostream>
#include <unordered_set>
#include "tet-mesh.h"


struct SolverParams {
public:
	float frameDt;
	int subSteps;
	vec3 g;
	int iterCount;

	int currMaterial;

	float k;
	float restLen;

	float u;
	float lambda;

	float m;

	float kc;
	float collisionThreshold;

	bool operator==(const SolverParams& other) const {
		return frameDt == other.frameDt && subSteps == other.subSteps && g == other.g && iterCount == other.iterCount && currMaterial == other.currMaterial && k == other.k && restLen == other.restLen && u == other.u && lambda == other.lambda && other.m == m && kc == other.kc && collisionThreshold == other.collisionThreshold;
			//&& (other.constraintGroupName == constraintGroupName || (other.constraintGroupName != nullptr && constraintGroupName != nullptr && strcmp(constraintGroupName, other.constraintGroupName) == 0));
	}
};

typedef int PhysicsMaterialID;
#define SIMPLE_SPRING 0
#define STVK_CLOTH 1
#define TET_SPRING 2
#define TET_NEOHOOK 3

struct FaceInfo {
	float restArea;
	mat2 invRestShape;
	array<int, 3> vertIDs; // vertIDs[i] = id => id's local triangle index is i
};

class VBDSolver {
public:
	VBDSolver(const vector<SolverParams>* params);

	void ResetSimulation(uPtr<HalfEdgeMesh> newStartPoseMesh = nullptr, uPtr<HalfEdgeMesh> collisionMeshSource = nullptr);
	void SimulateUpToFrame(uint frameIndex);
	inline HalfEdgeMesh* GetMesh(uint frameIndex) {
		if (frameIndex >= cachedPoses.size())
			std::cerr << "Getting Mesh From Solver doesn't have Frame!" << std::endl;
		return cachedPoses[frameIndex].get();
	}

	void VBDSolver::TruncateSimulation(int lastValidFrame);
	void VBDSolver::SetDirty();

	//
	/*int iterCount = 5;
	float frameDt = 1.0f / 24.0f;*/
	inline float stepDt() {
		return P().frameDt / (1.0f * P().subSteps);
	}
private:
	bool useTetMesh = false;

	int simulatingFrame;
	const vector<SolverParams>* cachedParams;
	inline const SolverParams& P() { return (*cachedParams)[simulatingFrame]; }

	bool dirty = false;
	vector<uPtr<HalfEdgeMesh>> cachedPoses;
	vector<uPtr<TetMesh>> cachedTetPoses;
	int lastSimulatedFrame;

	void ComputeCollisionForceAndHessian(Vertex* vert, vec3& collisionForce, mat3& collisionHessian);

	void SimulateOneFrame();
	void SimulateOneFrameTri();
	void SimulateOneFrameTet();
	vec3 PredictPosition(Vertex* vert, vec3 externalPos);
	vec3 PredictPositionCloth(const HalfEdgeMesh& mesh, Vertex* vert, vec3 externalPos);
	vec3 PredictPositionTetSpring(TetMesh&, Vertex* vert, vec3 externalPos);
	vec3 PredictPositionTetNeoHook(TetMesh&, Vertex* vert, vec3 externalPos);
	
	float planeHeight = -2.0f;
	float planeTilt = 0.0f; // angle in degrees
	void ComputePlaneCollision(vec3 planeNormal, vec3 planePoint, Vertex* vert, vec3& collisionForce, mat3& collisionHessian);
	void ComputeTriangleCollision(Vertex* vert, Vertex* a, Vertex* b, Vertex* c, vec3& collisionForce, mat3& collisionHessian);
	uPtr<HalfEdgeMesh> collisionMesh = nullptr;
	bool enableCollisionMesh = false;
	bool enableCollisionPlane = false;

	// For StVK Cloth, different ComputeHessian/ComputeForce functions can be written for different materials, but the 'element' changes too often to generalize (simple spring uses vert, cloth stvk uses triangle, later materials will use tetrahedrons)
	void ComputeClothFaceInfo();
	mat3 ComputeClothNeighborHessian(const HalfEdgeMesh& mesh, Face* face, Vertex*);
	vec3 ComputeClothNeighborForce(const HalfEdgeMesh& mesh, Face* face, Vertex* v);

	// For NeoHookean
	void ComputeNeoHookForceAndHessian(Vertex* vert, Tet* tet, vec3& force, mat3& hessian);

	//
	void ComputeInertiaForceAndHessian(vec3 vertPos, vec3 externalPos, vec3& inertiaForce, mat3& inertiaHessian);

	std::unordered_map<int, FaceInfo> facesInfo;
	std::unordered_set<int> constrainedVerts;
};