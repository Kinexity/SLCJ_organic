#pragma once
#ifndef geometry_h
#define geometry_h
#include "G4TriangularFacet.hh"
#include <vector>
#include <set>
#include <tuple>
#include <array>

std::vector<G4ThreeVector> getVertices(const G4TriangularFacet* facet);

// check if facets share a vertex
bool shareVertex(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2);

// function checking if two G4TriangularFacet objects are in the same plane
bool areFacetsCoplanar(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2);

// find intersection distance between line and facet's plane
G4double findIntersectionDistance(const G4ThreeVector startPoint, const G4ThreeVector direction, const G4TriangularFacet* facet);

// generate vectors 
std::vector<std::pair<G4ThreeVector, G4ThreeVector>> generateVertexVectorPairs(const G4TriangularFacet* facet);

// generate vectors 
std::vector<std::pair<G4ThreeVector, G4ThreeVector>> generateVertexEdgeVectorPairs(const G4TriangularFacet* facet);

// check if facets share an edge
bool containsEdge(const G4ThreeVector startPoint, const G4ThreeVector endPoint, const G4TriangularFacet* facet1);

// check if vector is in facet's plane
bool checkIfVectorIsOnFacetPlane(const G4ThreeVector startPoint, const G4ThreeVector direction, const G4TriangularFacet* facet);

// check if coplanar (!) point is placed inside the triangle
bool isPointInsideFacet(const G4ThreeVector& point, const G4TriangularFacet* facet);

// check if vector intersects facet's plane
bool doesVectorIntersectFacet(const G4ThreeVector startPoint, const G4ThreeVector vector, const G4TriangularFacet* facet);

// check if facets share an edge
bool shareEdge(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2);

bool doFacetsIntersect(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2);

G4double getPointDistanceFromFacetPlane(const G4ThreeVector point, const G4ThreeVector direction, const G4TriangularFacet* facet);

G4ThreeVector shiftPointToFacetPlaneXY(const G4ThreeVector point, const G4TriangularFacet* facet);

std::vector<std::vector<G4TriangularFacet*>> groupCoplanarFacets(const std::vector<G4TriangularFacet*>& facets);

std::vector<G4TriangularFacet*> spanFacets(std::set<G4ThreeVector> vertices);

std::map<G4ThreeVector, std::array<G4ThreeVector, 3>>
getVertexFacetsMap(std::vector<G4ThreeVector>& vertices, std::vector<std::vector<G4TriangularFacet*>>& coplanarGroups);

G4ThreeVector calculateInternalPoint(std::tuple<G4ThreeVector, std::array<G4ThreeVector, 3>> vertexNormals);

#endif // !geometry_h
