#include "VBDSolver.h"
#include <gtc/matrix_transform.hpp>
#include "helper/math.h"

bool IsConstrained(Vertex* vert) {
	return abs(vert->pos.x) >= 0.97 && vert->pos.y > -0.01;//vert->pos.y > 0.75f;
}

VBDSolver::VBDSolver() : cachedPoses(), lastSimulatedFrame(0) {}

void VBDSolver::ResetSimulation(uPtr<HalfEdgeMesh> newStartPoseMesh) {
	if (newStartPoseMesh == nullptr && cachedPoses.size() == 0) {
		std::cerr << "ERROR: ResetSimulation called without new start pose mesh when we don't have one!" << std::endl;
		return;
	}

	cachedPoses.resize(1);

	if (newStartPoseMesh != nullptr) {
		facesInfo.clear();
		constrainedVerts.clear();

		cachedPoses[0] = mkU<HalfEdgeMesh>(*newStartPoseMesh);
		cachedPoses[0]->TriangulateAllFaces();
		ComputeFaceInfo();

		for (const uPtr<Vertex>& v : cachedPoses[0]->vertices) {
			if (IsConstrained(v.get()))
				constrainedVerts.insert(v->id);
		}
	}

	lastSimulatedFrame = 0;
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
	if (lastSimulatedFrame > frameIndex) {
		// Use the past cache
		lastSimulatedFrame = frameIndex;
		cachedPoses.resize(lastSimulatedFrame + 1);
	} else {
		// Need to simulate more
		int numUpdates = frameIndex - lastSimulatedFrame;
		for (int i = 0; i < numUpdates; i++) {
			SimulateOneFrame();
		}
	}
}

vec3 VBDSolver::PredictPosition(Vertex* vert, vec3 externalPos) {
	if (vert->constrained) return vert->pos; // TODO need newest changes with map

	vec3 inertiaForce = -m / (stepDt() * stepDt()) * (vert->pos - externalPos);
	mat3 inertiaHessian = m / (stepDt() * stepDt()) * glm::identity<mat3>();

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
	if (vert->constrained) return vert->pos;

	vec3 inertiaForce = -m / (stepDt() * stepDt()) * (vert->pos - externalPos);
	mat3 inertiaHessian = m / (stepDt() * stepDt()) * glm::identity<mat3>();

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
	for (const uPtr<Face>& f : cachedPoses[0]->faces) {
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
			vp[i] = cachedPoses[0]->vertices[fi->vertIDs[i]]->pos;
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
	if (cachedPoses.size() <= lastSimulatedFrame) {
		std::cerr << "ERROR: SimulateOneFrame() called with not enough cachedPoses" << std::endl;
		return;
	}

	uPtr<HalfEdgeMesh> simulatingMesh = mkU<HalfEdgeMesh>(*cachedPoses[lastSimulatedFrame]);

	for (int i = 0; i < subSteps; i++) {
		// Predict external positions & save positions
		vector<vec3> oldPositions(simulatingMesh->vertices.size());
		vector<vec3> externalPredictedPositions(simulatingMesh->vertices.size());
		for (int i = 0; i < simulatingMesh->vertices.size(); i++) {
			vec3 externalAcc = g;
			oldPositions[i] = simulatingMesh->vertices[i]->pos;
			externalPredictedPositions[i] = simulatingMesh->vertices[i]->pos + simulatingMesh->vertices[i]->vel * stepDt() + externalAcc * stepDt() * stepDt();
		}

		for (int i = 0; i < iterCount; i++) {
			for (int i = 0; i < simulatingMesh->vertices.size(); i++) {
				Vertex* v = simulatingMesh->vertices[i].get();

				switch (currMaterial) {
				case SIMPLE_SPRING:
					v->pos = PredictPosition(v, externalPredictedPositions[i]);
					break;

				case STVK_CLOTH:
					v->pos = PredictPositionCloth(*simulatingMesh, v, externalPredictedPositions[i]);
					break;
				}
			}
		}

		for (int i = 0; i < simulatingMesh->vertices.size(); i++) {
			Vertex* v = simulatingMesh->vertices[i].get();
			v->vel = (1.0f / stepDt()) * (v->pos - oldPositions[i]);
			v->pos = oldPositions[i] + 0.98f * v->vel * stepDt(); // DAMPING
		}
	}

	++lastSimulatedFrame;
	cachedPoses.resize(lastSimulatedFrame + 1);
	cachedPoses[lastSimulatedFrame] = std::move(simulatingMesh);

}