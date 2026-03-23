#include "half-edge-mesh.h"

#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <GU/GU_Detail.h>
#include <GU/GU_PrimPoly.h>

void HalfEdgeMesh::CreateFromGUDetail(const GU_Detail* geo) {
    pointCountBound = geo->getNumPointOffsets();
    std::vector<std::unordered_set<GA_Offset>> adjList(pointCountBound);
    std::unordered_map<uint, HalfEdge*> symMap;

    // Create vertices
    for (GA_Iterator pointIter = GA_Iterator(geo->getPointRange()); !pointIter.atEnd(); ++pointIter) {
        GA_Offset pointOffset = pointIter.getOffset();
        UT_Vector3 hpos = geo->getPos3(pointOffset);
        vec3 pos = vec3(hpos.x(), hpos.y(), hpos.z());

        addVertex(pointOffset, pos);
    }

    // Create half edge mesh
    for (GA_Iterator primIterator = GA_Iterator(geo->getPrimitiveRange()); !primIterator.atEnd(); ++primIterator) {
        Face* f = addFace();

        const GEO_Primitive* prim = geo->getGEOPrimitive(*primIterator);
        int numVertices = prim->getVertexCount();

        HalfEdge* currEdge, * firstEdge;
        HalfEdge* prevEdge = nullptr;
        Vertex* prevVertex = getVertexAtOffset(prim->getPointOffset(numVertices - 1));
        Vertex* currVertex;

        for (int i = 0; i < numVertices; i++) {
            currVertex = getVertexAtOffset(prim->getPointOffset(i));

            // Create edge & next
            currEdge = addEdge();
            currEdge->face = f;
            currEdge->nextVertex = currVertex;
            currVertex->incomingEdge = currEdge;
            if (prevEdge == nullptr) {
                // First edge
                firstEdge = currEdge;
                f->edge = firstEdge;
            }
            else {
                prevEdge->next = currEdge;
            }

            // Sym
            uint pairID = vertexPairToID(currVertex, prevVertex);
            if (symMap.count(pairID) != 0) {
                currEdge->sym = symMap[pairID];
                symMap[pairID]->sym = currEdge;
            }
            else {
                symMap[pairID] = currEdge;
            }

            // Adjacency list
            int pntOffA = prim->getPointOffset(i);
            int pntOffB = prim->getPointOffset((i + 1) % numVertices);

            adjList[pntOffA].insert(pntOffB);
            adjList[pntOffB].insert(pntOffA);

            // Update prev
            prevVertex = currVertex;
            prevEdge = currEdge;
        }

        // Link to start
        currEdge->next = firstEdge;
    }

    // Debug log the half edge mesh
    std::cerr << "Debugging Half Edge Mesh " << std::endl;
    for (const uPtr<Face>& face : faces) {
        std::cerr << "\tFace " << face->id << std::endl;

        HalfEdge* currEdge = face->edge;
        do {
            std::cerr << "Edge " << currEdge->id << " has Vertex " << currEdge->nextVertex->id << " and NextPtr " << currEdge->next->id << " and Face " << currEdge->face->id << " and Sym " << currEdge->sym->id << std::endl;
            currEdge = currEdge->next;
        } while (currEdge != face->edge);
    }
}