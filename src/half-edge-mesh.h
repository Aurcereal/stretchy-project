#pragma once

#include "defines.h"
#include <vector>
#include <glm.hpp>

#include <SOP/SOP_Node.h>

using namespace std;
using namespace glm;

struct Vertex;
struct Face;
struct HalfEdge;

struct Vertex {
    uint id;

    GA_Offset pointOffset; // id within old mesh so points can maintain arb attributes that don't need to carry into sim, debatable whether should hvae this or just reconstruct from position ourself
    vec3 pos;
    HalfEdge* incomingEdge;
};

struct HalfEdge {
    uint id;

    HalfEdge* next;
    HalfEdge* sym;
    Vertex* nextVertex;
    Face* face;
};

struct Face {
    uint id;

    HalfEdge* edge;
};

class HalfEdgeMesh {
public:
    void CreateFromGUDetail(const GU_Detail*);
private:
    uint pointCountBound;
    vector<uPtr<Vertex>> vertices;
    vector<uPtr<HalfEdge>> halfEdges;
    vector<uPtr<Face>> faces;
    unordered_map<GA_Offset, uint> pointOffsetToIndex;

    inline Face* addFace() {
        faces.push_back(mkU<Face>());
        uint id = faces.size() - 1;
        faces[id]->id = id;
        return faces[id].get();
    }
    inline HalfEdge* addEdge() {
        halfEdges.push_back(mkU<HalfEdge>());
        uint id = halfEdges.size() - 1;
        halfEdges[id]->id = id;
        return halfEdges[id].get();
    }
    inline Vertex* addVertex(GA_Offset pointOffset, vec3 pos) {
        vertices.push_back(mkU<Vertex>());

        Vertex* v = vertices[vertices.size() - 1].get();
        v->pointOffset = pointOffset;
        v->pos = pos;
        v->id = vertices.size() - 1;

        // Update map
        pointOffsetToIndex[pointOffset] = vertices.size() - 1;

        return v;
    }
    inline Vertex* getVertexAtOffset(GA_Offset offset) {
        assert(pointOffsetToIndex.count(offset) > 0);
        return vertices[pointOffsetToIndex[offset]].get();
    }

    inline uint vertexPairToID(Vertex* v1, Vertex* v2) {
        return glm::min(v1->pointOffset, v2->pointOffset) * pointCountBound + glm::max(v1->pointOffset, v2->pointOffset);
    }
};
