#include "tet-mesh.h"

#include <iostream>


TetMesh::TetMesh(const TetMesh& other) {

    for (const uPtr<Vertex>& v : other.vertices) {
        vertices.push_back(mkU<Vertex>(*v));
    }

    for (const Tet& t : other.tets) {
        Tet newTet;
        newTet.DmInv = t.DmInv;
        for (int i = 0; i < 4; i++) {
            newTet.v[i] = vertices[t.v[i]->id].get();
        }
        tets.push_back(newTet);
    }

    surfaceTris = other.surfaceTris;

    PreCompute(); // copy constructor 
    // silly
    restLengths = other.restLengths;
}

#include <deque>

void TetMesh::FromTetgenOutput(const tetgenio& out, const HalfEdgeMesh& heMesh) {

    std::cout << "FromTetgenOutput: " << out.numberofpoints << " points, "
        << out.numberoftetrahedra << " tets, "
        << out.numberoftrifaces << " trifaces" << std::endl;

    vertices.clear();
    tets.clear();
    surfaceTris.clear();

    // Vertices
    for (int i = 0; i < out.numberofpoints; i++) {
        auto v = mkU<Vertex>();
        v->pos = vec3(
            (float)out.pointlist[i * 3 + 0],
            (float)out.pointlist[i * 3 + 1],
            (float)out.pointlist[i * 3 + 2]
        );
        v->id = i;
        vertices.push_back(std::move(v));
    }

    // Tets
    for (int i = 0; i < out.numberoftetrahedra; i++) {
        Tet t;
        for (int j = 0; j < 4; j++) {
            t.v[j] = vertices[out.tetrahedronlist[i * 4 + j]].get();
        }
        tets.push_back(t);
    }

    // Surface triangles
    for (int i = 0; i < out.numberoftrifaces; i++) {
        surfaceTris.push_back({
            (uint32_t)out.trifacelist[i * 3 + 0],
            (uint32_t)out.trifacelist[i * 3 + 1],
            (uint32_t)out.trifacelist[i * 3 + 2]
            });
    }

    // Setup BFS to propogate vel
    deque<int> visitQueue;
    vector<bool> accountedFor(vertices.size(), false);

    // set is surface flag
    for (auto& v : vertices) {
        v->isSurface = false;
    }
    for (int i = 0; i < out.numberoftrifaces; i++) {
        for (int j = 0; j < 3; j++) {
            int id = out.trifacelist[i * 3 + j];

            vertices[id]->isSurface = true;
            vertices[id]->pointOffset = heMesh.vertices[id]->pointOffset;
            vertices[id]->vel = heMesh.vertices[id]->vel;

            accountedFor[id] = true;
            visitQueue.push_back(id);
        }
    }

    PreCompute();

    // Run BFS to propogate vel
    while (visitQueue.size() > 0) {
        int vert = visitQueue[0];
        visitQueue.pop_front();

        vector<Vertex*>& neighbors = vertexNeighbors[vert];
        for (Vertex* neighbor : neighbors) {
            if (!accountedFor[neighbor->id]) {
                neighbor->vel = vertices[vert]->vel;

                accountedFor[neighbor->id] = true;
                visitQueue.push_back(neighbor->id);
            }
        }
    }
}

void TetMesh::PreCompute() {
    const int tetEdges[6][2] = { {0,1},{0,2},{0,3},{1,2},{1,3},{2,3} };
    for (Tet& tet : tets) {

        // Edges + rest lengths + neighbors
        for (auto& e : tetEdges) {
            Vertex* va = tet.v[e[0]];
            Vertex* vb = tet.v[e[1]];
            uint64_t key = VertexPairID(va, vb);
            if (restLengths.count(key) == 0) {
                restLengths[key] = glm::length(va->pos - vb->pos);
                vertexNeighbors[va->id].push_back(vb);
                vertexNeighbors[vb->id].push_back(va);
            }
        }

        // DmInv
        vec3 col0 = tet.v[1]->pos - tet.v[0]->pos;
        vec3 col1 = tet.v[2]->pos - tet.v[0]->pos;
        vec3 col2 = tet.v[3]->pos - tet.v[0]->pos;
        mat3 Dm = mat3(col0, col1, col2);
        tet.DmInv = glm::inverse(Dm);

        // Vertex tet map
        for (int i = 0; i < 4; i++) {
            vertexTets[tet.v[i]->id].push_back(&tet);
        }
    }
}

void TetMesh::ToRenderArrays(vector<vec3>* outPositions, vector<vec3>* outNormals, vector<uint32_t>* outIndices) {
    outPositions->clear();
    outNormals->clear();
    outIndices->clear();

    for (const auto& tri : surfaceTris) {
        vec3 a = vertices[tri[0]]->pos;
        vec3 b = vertices[tri[2]]->pos;
        vec3 c = vertices[tri[1]]->pos;

        vec3 normal = normalize(cross(b - a, c - a));

        uint32_t base = outPositions->size();
        outPositions->push_back(a);
        outPositions->push_back(b);
        outPositions->push_back(c);
        outNormals->push_back(normal);
        outNormals->push_back(normal);
        outNormals->push_back(normal);
        outIndices->push_back(base);
        outIndices->push_back(base + 1);
        outIndices->push_back(base + 2);
    }
}

void TetMesh::FromHalfEdge(const HalfEdgeMesh& heMesh) {
    std::cout << "Vertices: " << heMesh.vertices.size() << std::endl;
    std::cout << "Faces: " << heMesh.faces.size() << std::endl;

    tetgenio in, out;

    // Use unique vertices directly
    in.numberofpoints = heMesh.vertices.size();
    in.pointlist = new REAL[in.numberofpoints * 3];
    for (int i = 0; i < heMesh.vertices.size(); i++) {
        in.pointlist[i * 3 + 0] = heMesh.vertices[i]->pos.x;
        in.pointlist[i * 3 + 1] = heMesh.vertices[i]->pos.y;
        in.pointlist[i * 3 + 2] = heMesh.vertices[i]->pos.z;
    }

    // Build face list from half-edge faces
    in.numberoffacets = heMesh.faces.size();
    in.facetlist = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];
    for (int i = 0; i < heMesh.faces.size(); i++) {
        Face* f = heMesh.faces[i].get();
        Vertex* a = f->edge->nextVertex;
        Vertex* b = f->edge->next->nextVertex;
        Vertex* c = f->edge->next->next->nextVertex;

        tetgenio::facet& tf = in.facetlist[i];
        tf.numberofpolygons = 1;
        tf.polygonlist = new tetgenio::polygon[1];
        tf.numberofholes = 0;
        tf.holelist = nullptr;
        tetgenio::polygon& p = tf.polygonlist[0];
        p.numberofvertices = 3;
        p.vertexlist = new int[3];
        p.vertexlist[0] = a->id;
        p.vertexlist[1] = b->id;
        p.vertexlist[2] = c->id;
        in.facetmarkerlist[i] = 0;
    }

    tetgenbehavior b;
    b.parse_commandline((char*)"pq1.4Y");

    std::cout << "Calling tetrahedralize..." << std::endl;
    tetrahedralize(&b, &in, &out);
    std::cout << "Done. Points: " << out.numberofpoints << " Tets: " << out.numberoftetrahedra << std::endl;

    FromTetgenOutput(out, heMesh);
}

uPtr<HalfEdgeMesh> TetMesh::ToHalfEdge() const {
    auto heMesh = mkU<HalfEdgeMesh>();
    unordered_map<int, HalfEdge*> vertexPairToEdge;

    // Add vertices
    for (const auto& v : vertices) {
        if (!v->isSurface) continue;
        /*auto newVert = mkU<Vertex>();
        newVert->pos = v->pos;
        newVert->id = v->id;*/
        auto newVert = heMesh->addVertex(v->pointOffset, v->pos);
        newVert->id = v->id;
        newVert->pointOffset = v->pointOffset;
        newVert->vel = v->vel;
        //heMesh->vertices.push_back(std::move(newVert));
    }

    // Add faces from surface tris
    for (const auto& tri : surfaceTris) {
        Face* face = heMesh->addFace();
        HalfEdge* prevEdge = nullptr;
        HalfEdge* firstEdge = nullptr;

        int triVerts[3] = { (int)tri[0], (int)tri[1], (int)tri[2] };

        for (int i = 0; i < 3; i++) {
            int currVert = triVerts[i];
            int nextVert = triVerts[(i + 1) % 3];

            HalfEdge* edge = heMesh->addEdge();
            edge->nextVertex = heMesh->vertices[nextVert].get();
            heMesh->vertices[nextVert]->incomingEdge = edge;
            edge->face = face;
            face->edge = edge;

            // Sym linking
            int pairID = std::min(currVert, nextVert) * vertices.size() + std::max(currVert, nextVert);
            if (vertexPairToEdge.count(pairID)) {
                edge->sym = vertexPairToEdge[pairID];
                vertexPairToEdge[pairID]->sym = edge;
            }
            else {
                vertexPairToEdge[pairID] = edge;
            }

            if (prevEdge) prevEdge->next = edge;
            else firstEdge = edge;
            prevEdge = edge;
        }
        prevEdge->next = firstEdge;
        face->edge = firstEdge;
    }

    return heMesh;
}