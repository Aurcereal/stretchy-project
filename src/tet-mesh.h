#pragma once
#include "defines.h"
#include "half-edge-mesh.h"
#include <vector>
#include <array>
#include "tetgen.h"

using namespace glm;

struct Tet {
    Vertex* v[4];
    mat3 DmInv; // rest pose inverse
};

class TetMesh {


public:

    vector<uPtr<Vertex>> vertices;
    vector<Tet> tets;
    vector<array<uint32_t, 3>> surfaceTris;

    unordered_map<uint64_t, float> restLengths;
    unordered_map<uint32_t, vector<Vertex*>> vertexNeighbors;
    unordered_map<uint32_t, vector<Tet*>> vertexTets;

    void PreCompute();

    TetMesh() = default;
    TetMesh(const TetMesh&);

	friend class VBDSolver;

    // Build from TetGen output
    void FromTetgenOutput(const tetgenio& out, const HalfEdgeMesh&);

    // Extract surface as flat arrays for rendering
    void ToRenderArrays(vector<vec3>* outPositions, vector<vec3>* outNormals, vector<uint32_t>* outIndices);

	// Convert half-edge mesh
    void FromHalfEdge(const HalfEdgeMesh& heMesh);
    uPtr<HalfEdgeMesh> ToHalfEdge() const;

};

int VertexPairID(Vertex* a, Vertex* b);