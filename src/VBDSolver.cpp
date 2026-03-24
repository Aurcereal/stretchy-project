#include "VBDSolver.h"

VBDSolver::VBDSolver() : startPoseMesh(nullptr), lastSimulatedMesh(nullptr), lastSimulatedFrame(0) {}

void VBDSolver::ResetSimulation(uPtr<HalfEdgeMesh> newStartPoseMesh) {
	if (newStartPoseMesh != nullptr) {
		startPoseMesh = std::move(newStartPoseMesh);
	}
	if (startPoseMesh == nullptr) {
		std::cerr << "ERROR: ResetSimulation called without new start pose mesh when we don't have one!" << std::endl;
	}

	lastSimulatedFrame = 0;
	lastSimulatedMesh = mkU<HalfEdgeMesh>(*startPoseMesh); // Copy
}

void VBDSolver::SimulateUpToFrame(uint frameIndex) {
	// Past data is obsolete
	if (lastSimulatedFrame > frameIndex) {
		ResetSimulation();
	}

	for (int i = lastSimulatedFrame; i < frameIndex; i++) {
		SimulateOneFrame();
	}
}

constexpr vec3 g = vec3(0.0f, -0.98f, 0.0f);

void VBDSolver::SimulateOneFrame() {
	if (lastSimulatedMesh == nullptr) {
		std::cerr << "ERROR: SimulateOneFrame() called with no lastSimulatedMesh" << std::endl;
	}
	lastSimulatedFrame++;

	// Copy to new
	uPtr<HalfEdgeMesh> newMesh = mkU<HalfEdgeMesh>(*lastSimulatedMesh);

	// Simulate
	for (int i = 0; i < newMesh->vertices.size(); i++) {
		Vertex* v = newMesh->vertices[i].get();
		if (v->pos.y < -2.0f) {
			v->vel.y *= -1.0f;
			v->pos.y = -2.0f;
		}
		v->vel += dt * (g - 0.05f * v->pos);
		v->pos += dt * v->vel;
	}

	// Copy back
	lastSimulatedMesh = std::move(newMesh);
}