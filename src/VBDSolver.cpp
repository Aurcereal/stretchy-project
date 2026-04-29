#include "VBDSolver.h"
#include <gtc/matrix_transform.hpp>
#include "helper/math.h"

VBDSolver::VBDSolver() : startPoseMesh(nullptr), lastSimulatedMesh(nullptr), lastSimulatedFrame(0) {}

void VBDSolver::ResetSimulation(uPtr<HalfEdgeMesh> newStartPoseMesh) {
	if (newStartPoseMesh == nullptr && startPoseMesh == nullptr) {
		std::cerr << "ERROR: ResetSimulation called without new start pose mesh when we don't have one!" << std::endl;
		return;
	}

	if (newStartPoseMesh != nullptr) {
		facesInfo.clear();

		startPoseMesh = mkU<HalfEdgeMesh>(*newStartPoseMesh);
		startPoseMesh->TriangulateAllFaces();
		ComputeFaceInfo();
	}

	lastSimulatedFrame = 0;
	lastSimulatedMesh = mkU<HalfEdgeMesh>(*startPoseMesh); // Copy
}
//void VBDSolver::ResetSimulation(uPtr<HalfEdgeMesh> newStartPoseMesh) {
//	if (newStartPoseMesh != nullptr) {
//		startPoseMesh = std::move(newStartPoseMesh);
//	}
//	if (startPoseMesh == nullptr) {
//		std::cerr << "ERROR: ResetSimulation called without new start pose mesh when we don't have one!" << std::endl;
//	}
//
//	lastSimulatedFrame = 0;
//	lastSimulatedMesh = mkU<HalfEdgeMesh>(*startPoseMesh); // Copy
//	lastSimulatedMesh->TriangulateAllFaces();
//
//	// This assumes we triangulated (we did on import)
//	facesInfo.clear();
//	ComputeFaceInfo();
//}

void VBDSolver::SimulateUpToFrame(uint frameIndex) {
	// Past data is obsolete
	if (lastSimulatedFrame > frameIndex) {
		ResetSimulation();
	}

	int numUpdates = frameIndex - lastSimulatedFrame;
	for (int i = 0; i < numUpdates; i++) {
		SimulateOneFrame();
	}
}

bool IsConstrained(Vertex* vert) {
	return abs(vert->pos.x) >= 0.97 && vert->pos.y > -0.01;//vert->pos.y > 0.75f;
}

vec3 VBDSolver::PredictPosition(Vertex* vert, vec3 externalPos) {
	if (IsConstrained(vert)) return vert->pos; // TODO need newest changes with map

	vec3 inertiaForce = -m / (dt * dt) * (vert->pos - externalPos);
	mat3 inertiaHessian = m / (dt * dt) * glm::identity<mat3>();

	vec3 neighborForce = vec3(0);
	mat3 neighborHessian = mat3(0);

	HalfEdge* currEdge = vert->incomingEdge;
	do {
		Vertex* neighborVert = currEdge->sym->nextVertex;

		vec3 d = vert->pos - neighborVert->pos;
		float l = length(d);
		vec3 dNormalized = normalize(d);
		mat3 dNormalizedOuterProd = glm::outerProduct(dNormalized, dNormalized);

		neighborForce += -k * (l - restLen) * dNormalized;
		neighborHessian += k * (dNormalizedOuterProd + (1.0f / l) * (l - restLen) * (glm::identity<mat3>() - dNormalizedOuterProd));

		currEdge = currEdge->next->sym;
	} while (currEdge != vert->incomingEdge);

	vec3 force = inertiaForce + neighborForce;
	mat3 hessian = inertiaHessian + neighborHessian;

	vec3 deltaX = glm::inverse(hessian) * force;
	return vert->pos + deltaX;
}

vec3 VBDSolver::PredictPositionCloth(const HalfEdgeMesh& mesh, Vertex* vert, vec3 externalPos) {
	if (IsConstrained(vert)) return vert->pos;

	vec3 inertiaForce = -m / (dt * dt) * (vert->pos - externalPos);
	mat3 inertiaHessian = m / (dt * dt) * glm::identity<mat3>();

	vec3 neighborForce = vec3(0);
	mat3 neighborHessian = mat3(0);

	HalfEdge* currEdge = vert->incomingEdge;
	do {
		Face* currFace = currEdge->face;

		neighborForce += ComputeForce(mesh, currFace, vert);
		neighborHessian += ComputeHessian(mesh, currFace, vert);

		currEdge = currEdge->next->sym;
	} while (currEdge != vert->incomingEdge);

	vec3 force = inertiaForce + neighborForce;
	mat3 hessian = inertiaHessian + neighborHessian;

	vec3 deltaX = glm::inverse(hessian) * force;
	return vert->pos + deltaX;
}

void VBDSolver::ComputeFaceInfo() {
	// Foreach face, compute restArea and Dm^-1 using basis, record which vertices are which
	for (const uPtr<Face>& f : startPoseMesh->faces) {
		facesInfo[f->id] = FaceInfo();
		FaceInfo* fi = &facesInfo[f->id];

		int i = 0;
		HalfEdge* currEdge = f->edge;
		do {
			assert(i <= 2); // Triangles ONLY

			fi->vertIDs[i] = currEdge->nextVertex->id;

			currEdge = currEdge->next;
			++i;
		} while (currEdge != f->edge);

		array<vec3, 3> vp = array<vec3, 3>();
		for (int i = 0; i < 3; i++) {
			vp[i] = startPoseMesh->vertices[fi->vertIDs[i]]->pos;
		}

		// Forming 2D orthonormal basis out of triangle
		vec3 basisX = normalize(vp[1] - vp[0]);
		vec3 basisY = normalize(cross(basisX, cross(basisX, vp[2] - vp[0])));

		// Use above basis to make MATERIAL COORDINATES, like uvs but specifically chosen to maintain angles btwn edges & avoid non-uniform scale that might happen with normal uvs (its ok if our material coordinates are arbitrarily rotated or translated or uniformly scaled, but nothing else is ok)
		array<vec2, 3> materialVP = array<vec2, 3>();
		for (int i = 0; i < 3; i++) {
			materialVP[i] = vec2(dot(vp[i], basisX), dot(vp[i], basisY));
		}

		// Inverse rest shape
		mat2 restShape = mat2x2(materialVP[1] - materialVP[0], materialVP[2] - materialVP[0]);
		fi->invRestShape = inverse(restShape);

		// Rest area
		fi->restArea = 0.5f * length(cross(vp[1] - vp[0], vp[2] - vp[0]));
	}
}

mat3 VBDSolver::ComputeHessian(const HalfEdgeMesh& mesh, Face* face, Vertex* v) {
	//
	const FaceInfo& fInfo = facesInfo[face->id];
	array<vec3, 3> vp;
	for (int i = 0; i < 3; i++) {
		vp[i] = mesh.vertices[fInfo.vertIDs[i]]->pos;
	}

	//
	int localVertexID = -1;
	for (int i = 0; i < 3; i++) {
		if (fInfo.vertIDs[i] == v->id) {
			localVertexID = i;
			break;
		}
	}
	assert(localVertexID != -1);

	mat3 dxMat = glm::identity<mat3>();
	mat3 hessian = mat3();
	for (int i = 0; i < 3; i++) {
		vec3 dx = dxMat[i];

		glm::mat2x3 dDs;
		vec2 g;
		switch (localVertexID) {
		case 0:
			g = -(fInfo.invRestShape[0] + fInfo.invRestShape[1]);
			dDs = mat2x3(-dx, -dx); break;
		case 1:
			g = fInfo.invRestShape[0];
			dDs = mat2x3(dx, vec3(0)); break;
		case 2:
			g = fInfo.invRestShape[1];
			dDs = mat2x3(vec3(0), dx); break;
		}

		// dx/dX, change in world space pos per change in rest space pos
		mat2x3 deformationMat = mat2x3(vp[1] - vp[0], vp[2] - vp[0]) * fInfo.invRestShape;
		mat2x3 dDeformationMat = dDs * fInfo.invRestShape;

		mat2x2 strain = 0.5f * (transpose(deformationMat) * deformationMat - glm::identity<mat2x2>());
		mat2x2 dStrain = 0.5f * (transpose(dDeformationMat) * deformationMat + transpose(deformationMat) * dDeformationMat);

		mat2x3 dStress = dDeformationMat * (2 * u * strain + lambda * trace(strain) * glm::identity<mat2x2>()) + deformationMat * (2 * u * dStrain + lambda * trace(dStrain) * glm::identity<mat2x2>());
		vec3 df = -fInfo.restArea * dStress * g;

		hessian[i] = df;
	}

	return hessian;
}

vec3 VBDSolver::ComputeForce(const HalfEdgeMesh& mesh, Face* face, Vertex* v) {
	//
	const FaceInfo& fInfo = facesInfo[face->id];
	array<vec3, 3> vp;
	for (int i = 0; i < 3; i++) {
		vp[i] = mesh.vertices[fInfo.vertIDs[i]]->pos;
	}

	//
	int localVertexID = -1;
	for (int i = 0; i < 3; i++) {
		if (fInfo.vertIDs[i] == v->id) {
			localVertexID = i;
			break;
		}
	}
	assert(localVertexID != -1);

	vec2 g;
	switch (localVertexID) {
	case 0:
		g = -(fInfo.invRestShape[0] + fInfo.invRestShape[1]); break;
	case 1:
		g = fInfo.invRestShape[0]; break;
	case 2:
		g = fInfo.invRestShape[1]; break;
	}

	// dx/dX, change in world space pos per change in rest space pos
	mat2x3 deformationMat = mat2x3(vp[1] - vp[0], vp[2] - vp[0]) * fInfo.invRestShape;

	mat2x2 strain = 0.5f * (transpose(deformationMat) * deformationMat - glm::identity<mat2x2>());

	mat2x3 stress = deformationMat * (2 * u * strain + lambda * trace(strain) * glm::identity<mat2x2>());
	vec3 f = -fInfo.restArea * stress * g;

	return f;
}

void VBDSolver::SimulateOneFrame() {
	if (lastSimulatedMesh == nullptr) {
		std::cerr << "ERROR: SimulateOneFrame() called with no lastSimulatedMesh" << std::endl;
		return;
	}
	lastSimulatedFrame++;

	// Predict external positions & save positions
	vector<vec3> oldPositions(lastSimulatedMesh->vertices.size());
	vector<vec3> externalPredictedPositions(lastSimulatedMesh->vertices.size());
	for (int i = 0; i < lastSimulatedMesh->vertices.size(); i++) {
		vec3 externalAcc = g;
		oldPositions[i] = lastSimulatedMesh->vertices[i]->pos;
		externalPredictedPositions[i] = lastSimulatedMesh->vertices[i]->pos + lastSimulatedMesh->vertices[i]->vel * dt + externalAcc * dt * dt;
	}

#ifdef SCRATCH
	// Scratch
	uPtr<HalfEdgeMesh> scratchMesh = mkU<HalfEdgeMesh>(*lastSimulatedMesh);

	// Gaussian Seidel Iterations
	HalfEdgeMesh* currOldMesh, * currNewMesh;
	for (int i = 0; i < iterCount; i++) {
		currOldMesh = i % 2 == 0 ? lastSimulatedMesh.get() : scratchMesh.get();
		currNewMesh = i % 2 == 0 ? scratchMesh.get() : lastSimulatedMesh.get();

		for (int i = 0; i < scratchMesh->vertices.size(); i++) {
			currNewMesh->vertices[i]->pos = PredictPosition(currOldMesh->vertices[i].get(), externalPredictedPositions[i]);
		}
	}

	// Velocity update
	for (int i = 0; i < currNewMesh->vertices.size(); i++) {
		currNewMesh->vertices[i]->vel = (1.0f / dt) * (currNewMesh->vertices[i]->pos - oldPositions[i]);
	}

	// Copy back
	lastSimulatedMesh = mkU<HalfEdgeMesh>(*currNewMesh); // no need to copy again
#else
	for (int i = 0; i < iterCount; i++) {
		for (int i = 0; i < lastSimulatedMesh->vertices.size(); i++) {
			Vertex* v = lastSimulatedMesh->vertices[i].get();

			switch (currMaterial) {
			case SIMPLE_SPRING:
				v->pos = PredictPosition(v, externalPredictedPositions[i]);
				break;

			case STVK_CLOTH:
				v->pos = PredictPositionCloth(*lastSimulatedMesh, v, externalPredictedPositions[i]);
				break;
			}
		}
	}

	for (int i = 0; i < lastSimulatedMesh->vertices.size(); i++) {
		Vertex* v = lastSimulatedMesh->vertices[i].get();
		v->vel = (1.0f / dt) * (v->pos - oldPositions[i]);
	}
#endif
}