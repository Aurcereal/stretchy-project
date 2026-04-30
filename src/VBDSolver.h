#pragma once

#include "defines.h"
#include "half-edge-mesh.h"
#include <iostream>
#include <unordered_set>

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
	VBDSolver();

	void ResetSimulation(uPtr<HalfEdgeMesh> newStartPoseMesh = nullptr);
	void SimulateUpToFrame(uint frameIndex);
	inline HalfEdgeMesh* GetCurrentMesh() {
		if (lastSimulatedFrame >= cachedPoses.size()) 
			std::cerr << "Getting Mesh From Uninitialized Solver!" << std::endl;
		return cachedPoses[lastSimulatedFrame].get();
	}

	//
	int iterCount = 5;
	float frameDt = 1.0f / 24.0f;
	inline float stepDt() {
		return frameDt / (1.0f * subSteps);
	}
	uint subSteps = 1;
	vec3 g = vec3(0.0f, -0.98f, 0.0f);
	float m = 1.0f;

	// For Simple Spring Only
	float k = 150.0f;
	float restLen = 0.3;

	// For StVK Cloth Only
	float u = 1.0f;
	float lambda = 1.0f;

	//
	PhysicsMaterialID currMaterial = SIMPLE_SPRING;
private:
	vector<uPtr<HalfEdgeMesh>> cachedPoses;
	int lastSimulatedFrame;

	void SimulateOneFrame();
	vec3 PredictPosition(Vertex* vert, vec3 externalPos);
	vec3 PredictPositionCloth(const HalfEdgeMesh& mesh, Vertex* vert, vec3 externalPos);

	// For StVK Cloth, different ComputeHessian/ComputeForce functions can be written for different materials, but the 'element' changes too often to generalize (simple spring uses vert, cloth stvk uses triangle, later materials will use tetrahedrons)
	void ComputeFaceInfo();
	mat3 ComputeHessian(const HalfEdgeMesh& mesh, Face* face, Vertex*);
	vec3 ComputeForce(const HalfEdgeMesh& mesh, Face* face, Vertex* v);

	std::unordered_map<int, FaceInfo> facesInfo;
	std::unordered_set<int> constrainedVerts;
};