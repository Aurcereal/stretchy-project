#pragma once

#include "defines.h"
#include "half-edge-mesh.h"
#include <iostream>
#include <unordered_set>



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

	float kc = 1e6;
	float collisionThreshold = 0.1;

	bool operator==(const SolverParams& other) const {
		return frameDt == other.frameDt && subSteps == other.subSteps && g == other.g && iterCount == other.iterCount && currMaterial == other.currMaterial && k == other.k && restLen == other.restLen && u == other.u && lambda == other.lambda && other.m == m && kc == other.kc && collisionThreshold == other.collisionThreshold;
			//&& (other.constraintGroupName == constraintGroupName || (other.constraintGroupName != nullptr && constraintGroupName != nullptr && strcmp(constraintGroupName, other.constraintGroupName) == 0));
	}
};

typedef int PhysicsMaterialID;
#define SIMPLE_SPRING 0
#define STVK_CLOTH 1

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
	//uint subSteps = 1;
	//vec3 g = vec3(0.0f, -0.98f, 0.0f);
	//float m = 1.0f;

	//// For Simple Spring Only
	//float k = 150.0f;
	//float restLen = 0.3;

	//// For StVK Cloth Only
	//float u = 1.0f;
	//float lambda = 1.0f;

	//
	PhysicsMaterialID currMaterial = SIMPLE_SPRING;
private:
	int simulatingFrame;
	const vector<SolverParams>* cachedParams;
	inline const SolverParams& P() { return (*cachedParams)[simulatingFrame]; }

	bool dirty = false;
	vector<uPtr<HalfEdgeMesh>> cachedPoses;
	int lastSimulatedFrame;

	void SimulateOneFrame();
	vec3 PredictPosition(Vertex* vert, vec3 externalPos);
	vec3 PredictPositionCloth(const HalfEdgeMesh& mesh, Vertex* vert, vec3 externalPos);

	void ComputePlaneCollision(vec3 planeNormal, vec3 planePoint, Vertex* vert, vec3& collisionForce, mat3& collisionHessian);
	void ComputeTriangleCollision(Vertex* vert, Vertex* a, Vertex* b, Vertex* c, vec3& collisionForce, mat3& collisionHessian);
	uPtr<HalfEdgeMesh> collisionMesh = nullptr;
	bool enableCollisionMesh = false;

	// For StVK Cloth, different ComputeHessian/ComputeForce functions can be written for different materials, but the 'element' changes too often to generalize (simple spring uses vert, cloth stvk uses triangle, later materials will use tetrahedrons)
	void ComputeFaceInfo();
	mat3 ComputeHessian(const HalfEdgeMesh& mesh, Face* face, Vertex*);
	vec3 ComputeForce(const HalfEdgeMesh& mesh, Face* face, Vertex* v);

	std::unordered_map<int, FaceInfo> facesInfo;
	std::unordered_set<int> constrainedVerts;
};