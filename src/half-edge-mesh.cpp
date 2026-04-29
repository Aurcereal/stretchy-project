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
    /*std::cerr << "Debugging Half Edge Mesh " << std::endl;
    for (const uPtr<Face>& face : faces) {
        std::cerr << "\tFace " << face->id << std::endl;

        HalfEdge* currEdge = face->edge;
        do {
            std::cerr << "Edge " << currEdge->id << " has Vertex " << currEdge->nextVertex->id << " and NextPtr " << currEdge->next->id << " and Face " << currEdge->face->id << " and Sym " << currEdge->sym->id << std::endl;
            currEdge = currEdge->next;
        } while (currEdge != face->edge);
    }*/
}

void HalfEdgeMesh::LoadIntoExistingTopologicallySameHoudiniMesh(GU_Detail* geo) {
    for (const uPtr<Vertex>& vert : vertices) {
        geo->setPos3(vert->pointOffset, UT_Vector3(vert->pos.x, vert->pos.y, vert->pos.z));
    }
}

HalfEdgeMesh::HalfEdgeMesh(const HalfEdgeMesh& other) {
    this->pointCountBound = other.pointCountBound;

    unordered_map<int, HalfEdge*> halfEdgeIDToPtr;
    unordered_map<int, Face*> faceIDToPtr;
    unordered_map<int, Vertex*> vertexIDToPtr;

    for (const uPtr<Vertex>& oVert : other.vertices) {
        vertices.push_back(mkU<Vertex>(*oVert));
        vertexIDToPtr[vertices[vertices.size() - 1]->id] = vertices[vertices.size() - 1].get();
    }
    for (const uPtr<HalfEdge>& oEdge : other.halfEdges) {
        halfEdges.push_back(mkU<HalfEdge>(*oEdge));
        halfEdgeIDToPtr[halfEdges[halfEdges.size() - 1]->id] = halfEdges[halfEdges.size() - 1].get();
    }
    for (const uPtr<Face>& oFace : other.faces) {
        faces.push_back(mkU<Face>(*oFace));
        faceIDToPtr[faces[faces.size() - 1]->id] = faces[faces.size() - 1].get();
    }

    for (uPtr<Vertex>& v : vertices) {
        v->incomingEdge = halfEdgeIDToPtr[v->incomingEdge->id];
    }
    for (uPtr<HalfEdge>& h : halfEdges) {
        h->next = halfEdgeIDToPtr[h->next->id];
        h->sym = h->sym == nullptr ? nullptr : halfEdgeIDToPtr[h->sym->id];
        h->nextVertex = vertexIDToPtr[h->nextVertex->id];
        h->face = faceIDToPtr[h->face->id];
    }
    for (uPtr<Face>& f : faces) {
        f->edge = halfEdgeIDToPtr[f->edge->id];
    }
}

HalfEdgeMesh::HalfEdgeMesh() : pointCountBound() {}

void HalfEdgeMesh::TriangulateAllFaces() {
    int initialFaceCount = faces.size();
    for (int i = 0; i < initialFaceCount; ++i) {
        TriangulateFace(faces[i].get());
    }
}

void HalfEdgeMesh::TriangulateFace(Face* face) {
    HalfEdge* startEdge = face->edge;
    HalfEdge* currentEdge = startEdge->next;

    HalfEdge* prevEdge = startEdge;
    while (prevEdge->next != startEdge) prevEdge = prevEdge->next;
    Vertex* v1 = prevEdge->nextVertex;

    HalfEdge* currNextForNew = startEdge;

    while (currentEdge->next->next != startEdge) {

        // Half Edge pointing to start vertex
        HalfEdge* newEdge1 = addEdge();
        newEdge1->next = currNextForNew;
        newEdge1->nextVertex = v1; v1->incomingEdge = newEdge1;

        // Half Edge pointing away from start vertex
        HalfEdge* newEdge2 = addEdge();
        newEdge2->next = currentEdge->next;
        newEdge2->nextVertex = currentEdge->nextVertex;  currentEdge->nextVertex->incomingEdge = newEdge2;

        currentEdge->next = newEdge1;

        newEdge1->sym = newEdge2;
        newEdge2->sym = newEdge1;

        if (currentEdge == startEdge->next) {
            // Connect to existing face
            currentEdge->next->face = face;
            face->edge = currentEdge->next;
        }
        else {
            // Connect to new face
            Face* newFace = addFace();

            currentEdge->face = newFace; newFace->edge = currentEdge;
            currentEdge->next->face = newFace;
            currentEdge->next->next->face = newFace;
        }

        currNextForNew = newEdge2;
        currentEdge = newEdge2->next;

    }

    Face* newFace = addFace();
    currentEdge->next->next = currNextForNew;

    currentEdge->face = newFace;
    currentEdge->next->face = newFace;
    currentEdge->next->next->face = newFace; newFace->edge = currentEdge->next->next;
}

void TriangulateConvexFace(Face* f, vector<vec3>* positions, vector<vec3>* colors, vector<vec3>* normals, vector<uint32_t>* indices) {

    uint sideCount = 0;

    vec3 faceNormal = normalize(cross(
        f->edge->next->nextVertex->pos - f->edge->nextVertex->pos,
        f->edge->next->next->nextVertex->pos - f->edge->nextVertex->pos));

    HalfEdge* startEdge = f->edge;
    HalfEdge* currEdge = startEdge;
    do {
        positions->push_back(currEdge->nextVertex->pos);
        colors->push_back(vec3(1.0f)); // TODO: Can delete
        normals->push_back(faceNormal);
        currEdge = currEdge->next;
        sideCount++;
    } while (currEdge != startEdge);

    for (int i = 0; i < sideCount - 2; i++) {
        indices->push_back(positions->size() - sideCount);
        indices->push_back(positions->size() - sideCount + i + 1);
        indices->push_back(positions->size() - sideCount + i + 2);
    }

}