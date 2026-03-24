#pragma once

#include "defines.h"
#include "half-edge-mesh.h"
#include <iostream>

constexpr float dt = 1.0f / 24.0f; // TODO: make based on fps or an input

class VBDSolver {
public:
	VBDSolver();

	void ResetSimulation(uPtr<HalfEdgeMesh> newStartPoseMesh = nullptr);
	void SimulateUpToFrame(uint frameIndex);
	uPtr<HalfEdgeMesh> lastSimulatedMesh;
private:
	uPtr<HalfEdgeMesh> startPoseMesh;
	int lastSimulatedFrame;

	void SimulateOneFrame();
};