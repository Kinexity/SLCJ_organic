#include "geometry.h"
#include <map>
#include <unordered_map>
#include "G4SystemOfUnits.hh"

const double tolerance = 1e-8;

std::vector<G4ThreeVector> getVertices(const G4TriangularFacet* facet) {
	std::vector<G4ThreeVector> vertices;
	auto numberOfVertices = facet->GetNumberOfVertices();
	for (int vertexIndex = 0; vertexIndex < numberOfVertices; vertexIndex++) {
		auto vertex = facet->GetVertex(vertexIndex);
		vertices.push_back(vertex);
	}
	return vertices;
}

// check if facets share a vertex
bool shareVertex(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2) {
	// Check each vertex of facet1
	for (auto vertex1 : getVertices(facet1)) {

		// Check each vertex of facet2
		for (auto vertex2 : getVertices(facet2)) {

			if (vertex1 == vertex2) {
				// Found a shared vertex
				return true;
			}
		}
	}

	// No shared vertex found
	return false;
}

// function checking if two G4TriangularFacet objects are in the same plane
bool areFacetsCoplanar(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2) {
	// Compute the normal vectors of the planes defined by the facets
	G4ThreeVector normal1 = facet1->GetSurfaceNormal().unit();
	G4ThreeVector normal2 = facet2->GetSurfaceNormal().unit();

	// Check if the normal vectors are parallel (within a tolerance)
	G4double dotProduct = std::abs(normal1.dot(normal2));
	G4double planeDistance1 = std::abs(facet1->GetVertex(0) * normal1); // distance of facet's plane from the point (0,0,0)
	G4double planeDistance2 = std::abs(facet2->GetVertex(0) * normal2); // distance of facet's plane from the point (0,0,0)
	return std::abs(dotProduct - 1.0) < tolerance && std::abs(planeDistance1 - planeDistance2) < tolerance; // Adjust tolerance as needed
}

// find intersection distance between line and facet's plane
G4double findIntersectionDistance(const G4ThreeVector startPoint, const G4ThreeVector direction, const G4TriangularFacet* facet) {
	G4ThreeVector normal = facet->GetSurfaceNormal().unit();
	G4ThreeVector vertex = facet->GetVertex(0);

	G4double dotProduct = normal.dot(direction);

	G4double distance = (vertex - startPoint).dot(normal) / dotProduct;
	return distance;
}

// generate vectors 
std::vector<std::pair<G4ThreeVector, G4ThreeVector>> generateVertexVectorPairs(const G4TriangularFacet* facet) {
	std::vector<std::pair<G4ThreeVector, G4ThreeVector>> vertexPairs;
	std::vector<G4ThreeVector> pointVectors;

	// get vertices
	auto vertices = getVertices(facet);
	auto numberOfVertices = vertices.size();

	//calculate middle points, middle point vectors and edge vectors
	for (int vertexIndex = 0; vertexIndex < numberOfVertices; vertexIndex++) {
		auto middlePoint = (vertices[(vertexIndex + 2) % numberOfVertices] + vertices[(vertexIndex + 1) % numberOfVertices]) / 2; // middle point
		vertexPairs.emplace_back(vertices[vertexIndex], middlePoint - vertices[vertexIndex]); // middle point vector
		vertexPairs.emplace_back(vertices[vertexIndex], vertices[(vertexIndex + 1) % numberOfVertices] - vertices[vertexIndex]); // edge vector
	}

	return vertexPairs;
}

// generate vectors 
std::vector<std::pair<G4ThreeVector, G4ThreeVector>> generateVertexEdgeVectorPairs(const G4TriangularFacet* facet) {
	std::vector<std::pair<G4ThreeVector, G4ThreeVector>> vertexPairs;
	std::vector<G4ThreeVector> pointVectors;

	// get vertices
	auto vertices = getVertices(facet);
	auto numberOfVertices = vertices.size();

	//calculate middle points, middle point vectors and edge vectors
	for (int vertexIndex = 0; vertexIndex < numberOfVertices; vertexIndex++) {
		vertexPairs.emplace_back(vertices[vertexIndex], vertices[(vertexIndex + 1) % numberOfVertices] - vertices[vertexIndex]); // edge vector
	}

	return vertexPairs;
}

// check if facets share an edge
bool containsEdge(const G4ThreeVector startPoint, const G4ThreeVector endPoint, const G4TriangularFacet* facet1) {
	// Check each pair of vertices between the facets
	for (int i = 0; i < 3; ++i) {
		G4ThreeVector vertex1 = facet1->GetVertex(i);
		G4ThreeVector vertex2 = facet1->GetVertex((i + 1) % 3);

		if ((vertex1 == startPoint && vertex2 == endPoint) ||
			(vertex1 == endPoint && vertex2 == startPoint)) {
			// Found a matching edge
			return true;
		}
	}

	// No shared edge found
	return false;
}

// check if vector is in facet's plane
bool checkIfVectorIsOnFacetPlane(const G4ThreeVector startPoint, const G4ThreeVector direction, const G4TriangularFacet* facet) {
	G4ThreeVector normal = facet->GetSurfaceNormal().unit();
	G4ThreeVector vertex = facet->GetVertex(0);

	G4double absDotProduct = std::abs(normal.dot(direction)); // 0 if vector points perpendicularly to the plane
	auto absIntersectionDistance = std::abs(normal.dot(startPoint - vertex)); // 0 if vector starts in the plane

	return absDotProduct < tolerance&& absIntersectionDistance < tolerance;
}

bool checkIfFacetHasVectorInOtherFacetPlane(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2) {
	auto vertexVectorPairs1 = generateVertexEdgeVectorPairs(facet1);
	for (auto [vertex, vector] : vertexVectorPairs1) {
		if (checkIfVectorIsOnFacetPlane(vertex, vector, facet2)) {
			return true;
		}
	}
	auto vertexVectorPairs2 = generateVertexEdgeVectorPairs(facet2);
	for (auto [vertex, vector] : vertexVectorPairs2) {
		if (checkIfVectorIsOnFacetPlane(vertex, vector, facet1)) {
			return true;
		}
	}
	return false;
};

// check if coplanar (!) point is placed inside the triangle
bool isPointInsideFacet(const G4ThreeVector& point, const G4TriangularFacet* facet) {
	std::vector<G4ThreeVector> vertices, edgeVectors, pointVectors;
	std::vector<G4double> dotProducts;
	auto numberOfVertices = facet->GetNumberOfVertices();
	for (int vertexIndex = 0; vertexIndex < numberOfVertices; vertexIndex++) {
		auto vertex = facet->GetVertex(vertexIndex);
		vertices.push_back(vertex);
		pointVectors.push_back(point - vertex);
	}

	for (int vertexIndex = 0; vertexIndex < numberOfVertices; vertexIndex++) {
		edgeVectors.push_back(vertices[(vertexIndex + 1) % numberOfVertices] - vertices[vertexIndex]);
	}

	G4ThreeVector normal = facet->GetSurfaceNormal().unit();

	for (int vertexIndex = 0; vertexIndex < facet->GetNumberOfVertices(); vertexIndex++) {
		dotProducts.push_back(normal.dot(edgeVectors[vertexIndex].cross(pointVectors[vertexIndex])));
	}

	return (dotProducts[0] >= 0 && dotProducts[1] >= 0 && dotProducts[2] >= 0) ||
		(dotProducts[0] <= 0 && dotProducts[1] <= 0 && dotProducts[2] <= 0);
}

// check if vector intersects facet's plane
bool doesVectorIntersectFacet(const G4ThreeVector startPoint, const G4ThreeVector vector, const G4TriangularFacet* facet) {
	if (containsEdge(startPoint, startPoint + vector, facet)) {
		return false;
	}
	auto vectorDirection = vector.unit();
	G4double intersectionDistance = findIntersectionDistance(startPoint, vectorDirection, facet);
	auto vectorLenght = vector.mag();

	// check if the vector actually intersects the plane or it's only the line it's on
	if (intersectionDistance < 0. || intersectionDistance > vectorLenght) {
		return false;
	}
	auto intersectionPoint = startPoint + intersectionDistance * vectorDirection;
	return isPointInsideFacet(intersectionPoint, facet);
}

// check if facets share an edge
bool shareEdge(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2) {
	// Check each pair of vertices between the facets
	for (int i = 0; i < 3; ++i) {
		G4ThreeVector vertex1 = facet1->GetVertex(i);
		G4ThreeVector vertex2 = facet1->GetVertex((i + 1) % 3);

		if (containsEdge(vertex1, vertex2, facet2)) {
			// Found a matching edge
			return true;
		}
	}

	// No shared edge found
	return false;
}

bool doFacetsIntersect(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2) {
	if (checkIfFacetHasVectorInOtherFacetPlane(facet1, facet2)) {
		return false;
	}
	auto vertexVectorPairs1 = generateVertexVectorPairs(facet1);
	for (auto [vertex, vector] : vertexVectorPairs1) {
		if (doesVectorIntersectFacet(vertex, vector, facet2)) {
			return true;
		}
	}
	auto vertexVectorPairs2 = generateVertexVectorPairs(facet2);
	for (auto [vertex, vector] : vertexVectorPairs2) {
		if (doesVectorIntersectFacet(vertex, vector, facet1)) {
			return true;
		}
	}
	return false;
}

G4double getPointDistanceFromFacetPlane(const G4ThreeVector point, const G4ThreeVector direction, const G4TriangularFacet* facet) {
	return findIntersectionDistance(point, direction, facet);
};

G4ThreeVector shiftPointToFacetPlaneXY(const G4ThreeVector point, const G4TriangularFacet* facet) {
	auto normal = facet->GetSurfaceNormal().unit();
	normal.setZ(0.); // no z axis shift
	auto distance = getPointDistanceFromFacetPlane(point, normal, facet);
	return point + distance * normal;
}

std::vector<std::vector<G4TriangularFacet*>> groupCoplanarFacets(const std::vector<G4TriangularFacet*>& facets)
{
	std::unordered_map<G4TriangularFacet*, int> facetGroups;
	std::vector<std::vector<G4TriangularFacet*>> result;

	int groupIndex = 0;
	for (const auto& facet : facets) {
		bool foundGroup = false;
		for (auto& group : result) {
			if (areFacetsCoplanar(facet, group[0])) {
				facetGroups[facet] = groupIndex;
				group.push_back(facet);
				foundGroup = true;
				break;
			}
		}
		if (!foundGroup) {
			facetGroups[facet] = groupIndex;
			result.push_back({ facet });
			groupIndex++;
		}
	}

	return result;
}

std::vector<G4TriangularFacet*> spanFacets(std::set<G4ThreeVector> vertices) {
	std::vector<G4TriangularFacet*> facets;
	std::vector<G4ThreeVector> verts;
	auto vertex = *vertices.begin();
	vertices.erase(vertex);
	auto baselineVertex = *vertices.begin();
	auto baselineVec = baselineVertex - vertex;
	std::map<double, G4ThreeVector> angles;
	for (auto vertex2 : vertices) {
		auto vec = vertex2 - vertex;
		auto cs = baselineVec.dot(vec) / (vec.mag() * baselineVec.mag());
		auto ss = baselineVec.cross(vec).mag() / (vec.mag() * baselineVec.mag());
		auto angle = std::atan2(ss, cs);
		angles[angle] = vertex2;
	}
	for (auto [angle, vrt] : angles) {
		verts.push_back(vrt);
	}
	for (int i = 0; i < verts.size() - 1; i++) {
		facets.push_back(new G4TriangularFacet(vertex, verts[i], verts[i + 1], ABSOLUTE));

	}
	return facets;
}

bool isPointFacetVertex(const G4ThreeVector& point, const G4TriangularFacet* facet) {
	auto vertices = getVertices(facet);
	return std::any_of(vertices.begin(), vertices.end(), [&](G4ThreeVector vertex) { return point == vertex; });
}

std::map<G4ThreeVector, std::array<G4ThreeVector, 3>>
getVertexFacetsMap(std::vector<G4ThreeVector>& vertices, std::vector<std::vector<G4TriangularFacet*>>& coplanarGroups)
{
	std::map<G4ThreeVector, std::array<G4ThreeVector, 3>> vertexFacetsMap;

	// Iterate over each vertex
	for (const auto& vertex : vertices) {
		std::array<G4ThreeVector, 3> facetTuple;
		int i = 0;
		// Iterate over each coplanar group
		for (const auto& group : coplanarGroups) {
			// Find the first facet in the group that contains the vertex
			auto iter = std::find_if(group.begin(), group.end(), [&](const G4TriangularFacet* facet) {
				return isPointFacetVertex(vertex, facet);
				});

			if (iter != group.end()) {
				// Store the facet in the tuple
				facetTuple[i++] = (*iter)->GetSurfaceNormal();
			}
			if (i == 3) {
				break;
			}
		}
		// Add the vertex and the corresponding facet tuple to the map
		vertexFacetsMap[vertex] = facetTuple;
	}

	return vertexFacetsMap;
}

G4ThreeVector calculateInternalPoint(std::tuple<G4ThreeVector, std::array<G4ThreeVector, 3>> vertexNormals)
{
	auto [vertex, normals] = vertexNormals;

	for (auto& normal : normals) {
		if (vertex.dot(normal) < 0) {
			normal *= -1;
		}
	}

	// Calculate the determinant of the matrix formed by the normal vectors
	double det = normals[0].dot(normals[1].cross(normals[2]));

	// Get the distance from the origin to each plane
	double d1 = normals[0].dot(vertex - 1 * mm * normals[0]);
	double d2 = normals[1].dot(vertex - 1 * mm * normals[1]);
	double d3 = normals[2].dot(vertex - 1 * mm * normals[2]);

	// Calculate the intersection point using Cramer's rule
	double x = (d1 * normals[1].cross(normals[2]).x() + d2 * normals[2].cross(normals[0]).x() + d3 * normals[0].cross(normals[1]).x()) / det;
	double y = (d1 * normals[1].cross(normals[2]).y() + d2 * normals[2].cross(normals[0]).y() + d3 * normals[0].cross(normals[1]).y()) / det;
	double z = (d1 * normals[1].cross(normals[2]).z() + d2 * normals[2].cross(normals[0]).z() + d3 * normals[0].cross(normals[1]).z()) / det;

	return G4ThreeVector(x, y, z);
}