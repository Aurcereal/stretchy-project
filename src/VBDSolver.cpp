#include "VBDSolver.h"
#include <gtc/matrix_transform.hpp>
#include "helper/math.h"

bool IsConstrained(Vertex* vert) {
	return abs(vert->pos.x) >= 0.97 && vert->pos.y > -0.01;//vert->pos.y > 0.75f;
}

VBDSolver::VBDSolver(const vector<SolverParams>* params) : cachedParams(params), cachedPoses(), lastSimulatedFrame(0) {}

void VBDSolver::ResetSimulation(uPtr<HalfEdgeMesh> newStartPoseMesh, uPtr<HalfEdgeMesh> collisionMeshSource) {
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
		cachedPoses[0]->ComputeRestLengths();
		ComputeFaceInfo();
		

		for (const uPtr<Vertex>& v : cachedPoses[0]->vertices) {
			if (IsConstrained(v.get()))
				constrainedVerts.insert(v->id);
		}
	}

	if (collisionMeshSource != nullptr) {
		collisionMesh = mkU<HalfEdgeMesh>(*collisionMeshSource);
		collisionMesh->TriangulateAllFaces();
		enableCollisionMesh = true;
	}

	lastSimulatedFrame = 0;
}

void VBDSolver::TruncateSimulation(int lastValidFrame) {
	if (lastSimulatedFrame <= lastValidFrame) {
		return;
	}
	else {
		cachedPoses.resize(lastValidFrame + 1);
		lastSimulatedFrame = lastValidFrame;
	}
}
void VBDSolver::SetDirty() {
	dirty = true;
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


void VBDSolver::ComputePlaneCollision(vec3 planeNormal, vec3 planePoint, Vertex* vert, vec3& collisionForce, mat3& collisionHessian) {
	float d = fmax(0.0f, dot(planePoint - vert->pos, planeNormal));
	if (d > 0.0f) {

		// calculate force and hessian
		vec3 floorCollisionForce = P().kc * d * planeNormal;
		mat3 floorCollisionHessian = P().kc * glm::outerProduct(planeNormal, planeNormal);

		// output
		collisionForce += floorCollisionForce;
		collisionHessian += floorCollisionHessian;

	}
}

void VBDSolver::ComputeTriangleCollision(Vertex* vert, Vertex* a, Vertex* b, Vertex* c, vec3& collisionForce, mat3& collisionHessian) {
	if (vert == a || vert == b || vert == c) return;

	vec3 normal = normalize(cross(b->pos - a->pos, c->pos - a->pos));
	float d = dot(vert->pos - a->pos, normal);

	if (d > -P().collisionThreshold && d < 0.0f) {
		vec3 proj = vert->pos - d * normal;

		vec3 ab = b->pos - a->pos, ac = c->pos - a->pos, ap = proj - a->pos;
		float d00 = dot(ab, ab), d01 = dot(ab, ac), d11 = dot(ac, ac);
		float d20 = dot(ap, ab), d21 = dot(ap, ac);
		float denom = d00 * d11 - d01 * d01;
		float v = (d11 * d20 - d01 * d21) / denom;
		float w = (d00 * d21 - d01 * d20) / denom;
		float u = 1.0f - v - w;

		if (u < 0.f || v < 0.f || w < 0.f) return;

		collisionForce += P().kc * (-d) * normal;
		collisionHessian += P().kc * outerProduct(normal, normal);
	}
}

void VBDSolver::SimulateUpToFrame(uint frameIndex) {
	if (dirty && frameIndex == 1) {
		dirty = false;
		ResetSimulation();
	}

	if (frameIndex < lastSimulatedFrame) {
		return;
	} else {
		// Need to simulate more
		int numUpdates = frameIndex - lastSimulatedFrame;
		for (int i = 0; i < numUpdates; i++) {
			SimulateOneFrame();
		}
	}
}

vec3 VBDSolver::PredictPosition(Vertex* vert, vec3 externalPos) {
	if (constrainedVerts.count(vert->id) != 0) return vert->pos;

	vec3 inertiaForce = -P().m / (stepDt() * stepDt()) * (vert->pos - externalPos);
	mat3 inertiaHessian = P().m / (stepDt() * stepDt()) * glm::identity<mat3>();

	vec3 neighborForce = vec3(0);
	mat3 neighborHessian = mat3(0);

	vec3 collisionForce = vec3(0);
	mat3 collisionHessian = mat3(0);

	ComputePlaneCollision(vec3(0, 1, 0), vec3(0, -3, 0), vert, collisionForce, collisionHessian);

	// collision against mesh triangles
	if (enableCollisionMesh && collisionMesh != nullptr) {
		for (const uPtr<Face>& f : collisionMesh->faces) {
			// Rightnow: assuming mesh is triangle only
			Vertex* a = f->edge->nextVertex;
			Vertex* b = f->edge->next->nextVertex;
			Vertex* c = f->edge->next->next->nextVertex;
			ComputeTriangleCollision(vert, a, b, c, collisionForce, collisionHessian);
		}
	}

	HalfEdge* currEdge = vert->incomingEdge;
	do {
		Vertex* neighborVert = currEdge->sym->nextVertex;

		vec3 d = vert->pos - neighborVert->pos;
		float l = length(d);
		vec3 dNormalized = normalize(d);
		mat3 dNormalizedOuterProd = glm::outerProduct(dNormalized, dNormalized);

		float edgeRestLength = cachedPoses[0]->GetRestLength(vert, neighborVert) * P().restLen;

		neighborForce += -P().k * (l - edgeRestLength) * dNormalized;
		neighborHessian += P().k * (dNormalizedOuterProd + (1.0f / l) * (l - edgeRestLength) * (glm::identity<mat3>() - dNormalizedOuterProd));

		currEdge = currEdge->next->sym;
	} while (currEdge != vert->incomingEdge);

	//inertiaForce = vec3(0);
	//inertiaHessian = mat3(0);

	//neighborForce = vec3(0);
	//neighborHessian = mat3(0);

	//collisionForce = vec3(0);
	//collisionHessian = mat3(0);

	vec3 force = inertiaForce + neighborForce + collisionForce;
	mat3 hessian = inertiaHessian + neighborHessian + collisionHessian;

	vec3 deltaX = glm::inverse(hessian) * force;
	return vert->pos + deltaX;
}

vec3 VBDSolver::PredictPositionCloth(const HalfEdgeMesh& mesh, Vertex* vert, vec3 externalPos) {
	if (vert->constrained) return vert->pos;

	vec3 inertiaForce = -P().m / (stepDt() * stepDt()) * (vert->pos - externalPos);
	mat3 inertiaHessian = P().m / (stepDt() * stepDt()) * glm::identity<mat3>();

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

		mat2x3 dStress = dDeformationMat * (2 * P().u * strain + P().lambda * trace(strain) * glm::identity<mat2x2>()) + deformationMat * (2 * P().u * dStrain + P().lambda * trace(dStrain) * glm::identity<mat2x2>());
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

	mat2x3 stress = deformationMat * (2 * P().u * strain + P().lambda * trace(strain) * glm::identity<mat2x2>());
	vec3 f = -fInfo.restArea * stress * g;

	return f;
}

void VBDSolver::SimulateOneFrame() {
	if (cachedPoses.size() <= lastSimulatedFrame) {
		std::cerr << "ERROR: SimulateOneFrame() called with not enough cachedPoses" << std::endl;
		return;
	}

	uPtr<HalfEdgeMesh> simulatingMesh = mkU<HalfEdgeMesh>(*cachedPoses[lastSimulatedFrame]);
	simulatingFrame = lastSimulatedFrame + 1;

	for (int s = 0; s < P().subSteps; s++) {
		// Predict external positions & save positions
		vector<vec3> oldPositions(simulatingMesh->vertices.size());
		vector<vec3> externalPredictedPositions(simulatingMesh->vertices.size());
		for (int i = 0; i < simulatingMesh->vertices.size(); i++) {
			vec3 externalAcc = P().g;
			oldPositions[i] = simulatingMesh->vertices[i]->pos;
			externalPredictedPositions[i] = simulatingMesh->vertices[i]->pos + simulatingMesh->vertices[i]->vel * stepDt() + externalAcc * stepDt() * stepDt();
		}

		for (int j = 0; j < P().iterCount; j++) {
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
			v->vel *= 0.98f; // DAMPING
		}
	}

	++lastSimulatedFrame;
	cachedPoses.resize(lastSimulatedFrame + 1);
	cachedPoses[lastSimulatedFrame] = std::move(simulatingMesh);

}