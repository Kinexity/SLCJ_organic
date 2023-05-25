// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
/////////////////////////////////////////////////////////////////////////
//
// V.Guadilla 2019
//
////////////////////////////////////////////////////////////////////////

#include "SLCJDetectorConstruction.hh"

// GEANT4 //
#include "globals.hh"
#include "G4Transform3D.hh"
#include "G4NistManager.hh"

#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4NistManager.hh"
#include <fstream>
#include <format>
#include <map>
#include <tuple>
#include <set>

#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"
#include <G4UnionSolid.hh>
//#include "TMath.h"
#include "G4Polycone.hh"

#include "F02ElectricFieldSetup.hh"

// Global file with the events
std::ifstream geoInputFile;

inline std::string filename_string(std::string path_str) {
	return path_str.substr(path_str.rfind("\\") + 1, path_str.size() - path_str.rfind("\\") - 1);
};

#define _endl_ " (" << filename_string(__FILE__) << "; " << __LINE__ << ")" << '\n'
#define checkpoint std::cout << "checkpoint" << _endl_

const double tolerance = 1e-8;

// check if facets share a vertex
bool shareVertex(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2) {
	// Check each vertex of facet1
	for (int i = 0; i < 3; ++i) {
		G4ThreeVector vertex1 = facet1->GetVertex(i);

		// Check each vertex of facet2
		for (int j = 0; j < 3; ++j) {
			G4ThreeVector vertex2 = facet2->GetVertex(j);

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
bool areFacetsInSamePlane(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2) {
	// Compute the normal vectors of the planes defined by the facets
	G4ThreeVector normal1 = facet1->GetSurfaceNormal().unit();
	G4ThreeVector normal2 = facet2->GetSurfaceNormal().unit();

	// Check if the normal vectors are parallel (within a tolerance)
	G4double dotProduct = std::abs(normal1.dot(normal2));
	G4double planeDistance1 = std::abs(facet1->GetCircumcentre() * normal1); // distance of facet's plane from the point (0,0,0)
	G4double planeDistance2 = std::abs(facet2->GetCircumcentre() * normal2); // distance of facet's plane from the point (0,0,0)
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
	std::vector<G4ThreeVector> vertices, pointVectors;

	// get vertices
	auto numberOfVertices = facet->GetNumberOfVertices();
	for (int vertexIndex = 0; vertexIndex < numberOfVertices; vertexIndex++) {
		auto vertex = facet->GetVertex(vertexIndex);
		vertices.push_back(vertex);
	}

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
	std::vector<G4ThreeVector> vertices, pointVectors;

	// get vertices
	auto numberOfVertices = facet->GetNumberOfVertices();
	for (int vertexIndex = 0; vertexIndex < numberOfVertices; vertexIndex++) {
		auto vertex = facet->GetVertex(vertexIndex);
		vertices.push_back(vertex);
	}

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

// find intersection distance between line and facet's plane
bool checkIfVectorIsPerpendicularToFacetPlane(const G4ThreeVector startPoint, const G4ThreeVector direction, const G4TriangularFacet* facet) {
	G4ThreeVector normal = facet->GetSurfaceNormal().unit();
	G4ThreeVector vertex = facet->GetVertex(0);

	G4double dotProduct = normal.dot(direction);

	return dotProduct < tolerance && findIntersectionDistance(startPoint, direction, facet) < tolerance;
}

bool checkIfFacetHasVectorPerpendicularToOtherFacetPlane(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2) {
	auto vertexVectorPairs1 = generateVertexEdgeVectorPairs(facet1);
	for (auto [vertex, vector] : vertexVectorPairs1) {
		if (checkIfVectorIsPerpendicularToFacetPlane(vertex, vector, facet2)) {
			return true;
		}
	}
	auto vertexVectorPairs2 = generateVertexEdgeVectorPairs(facet2);
	for (auto [vertex, vector] : vertexVectorPairs2) {
		if (checkIfVectorIsPerpendicularToFacetPlane(vertex, vector, facet1)) {
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
bool doesVectorIntersectFacetPlane(const G4ThreeVector startPoint, const G4ThreeVector vector, const G4TriangularFacet* facet) {
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

		for (int j = 0; j < 3; ++j) {
			G4ThreeVector otherVertex1 = facet2->GetVertex(j);
			G4ThreeVector otherVertex2 = facet2->GetVertex((j + 1) % 3);

			if ((vertex1 == otherVertex1 && vertex2 == otherVertex2) ||
				(vertex1 == otherVertex2 && vertex2 == otherVertex1)) {
				// Found a matching edge
				return true;
			}
		}
	}

	// No shared edge found
	return false;
}

bool doFacetsIntersect(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2) {
	if (checkIfFacetHasVectorPerpendicularToOtherFacetPlane(facet1, facet2)) {
		return false;
	}
	auto vertexVectorPairs1 = generateVertexVectorPairs(facet1);
	for (auto [vertex, vector] : vertexVectorPairs1) {
		if (doesVectorIntersectFacetPlane(vertex, vector, facet2)) {
			return true;
		}
	}
	auto vertexVectorPairs2 = generateVertexVectorPairs(facet2);
	for (auto [vertex, vector] : vertexVectorPairs2) {
		if (doesVectorIntersectFacetPlane(vertex, vector, facet1)) {
			return true;
		}
	}
	return false;
}

bool doFacetsIntersect2(const G4TriangularFacet* facet1, const G4TriangularFacet* facet2) {
	auto normal = facet1->GetSurfaceNormal().unit();
	//std::cout << normal << '\n';
	G4ThreeVector vertexComp = facet1->GetVertex(0);
	G4ThreeVector vertex1 = facet2->GetVertex(0);
	//std::cout << vertex1 << '\n';
	G4ThreeVector vertex2 = facet2->GetVertex(1);
	//std::cout << vertex2 << '\n';
	G4ThreeVector vertex3 = facet2->GetVertex(2);
	//std::cout << vertex3 << '\n';

	G4ThreeVector vector1 = vertex1 - vertexComp;
	//std::cout << vector1 << '\n';
	G4ThreeVector vector2 = vertex2 - vertexComp;
	//std::cout << vector2 << '\n';
	G4ThreeVector vector3 = vertex3 - vertexComp;
	//std::cout << vector3 << '\n';

	G4double dotProduct1 = normal * vector1;
	//std::cout << dotProduct1 << '\n';
	G4double dotProduct2 = normal * vector2;
	//std::cout << dotProduct2 << '\n';
	G4double dotProduct3 = normal * vector3;
	//std::cout << dotProduct3 << '\n';

	return !((dotProduct1 >= -tolerance && dotProduct2 >= -tolerance && dotProduct3 >= -tolerance) ||
		(dotProduct1 <= tolerance && dotProduct2 <= tolerance && dotProduct3 <= tolerance));
}

G4double getPointDistanceFromFacetPlane(const G4ThreeVector point, const G4ThreeVector direction, const G4TriangularFacet* facet) {
	return findIntersectionDistance(point, direction, facet);
};

G4ThreeVector shiftPointToFacetPlaneXY(const G4ThreeVector point, const G4TriangularFacet* facet) {
	auto normal = facet->GetSurfaceNormal().unit();
	normal.setZ(0.);
	auto distance = getPointDistanceFromFacetPlane(point, normal, facet);
	return point + distance * normal;
}

SLCJDetectorConstruction::SLCJDetectorConstruction(G4double foodV) : foodVolume(foodV) {}

G4VPhysicalVolume* SLCJDetectorConstruction::Construct() {

	G4NistManager* nistManager = G4NistManager::Instance();

	//----------------------------------------------------
	// Materials definitions
	//----------------------------------------------------

	G4String name, symbol;             //a=mass of a mole;
	G4double a, z, density;            //z=mean number of protons;  
	G4int ncomponents, natoms;
	G4double fractionmass;
	G4double abundance;
	G4State state;
	G4double pressure;
	G4double temperature;

	G4int iz, n;                       //iz=number of protons  in an isotope; 
 //                                      // n=number of nucleons in an isotope;

 //
 // define Elements
 //

   // H
	G4Element* H = nistManager->FindOrBuildElement("H");

	// He
	G4Element* He = nistManager->FindOrBuildElement("He");

	// O
	G4Element* O = nistManager->FindOrBuildElement("O");

	// B
	G4Element* B = nistManager->FindOrBuildElement("B");

	// Al
	G4Element* Al = nistManager->FindOrBuildElement("Al");

	// C
	G4Element* C = nistManager->FindOrBuildElement("C");

	// Si
	G4Element* Si = nistManager->FindOrBuildElement("Si");

	// K
	G4Element* K = nistManager->FindOrBuildElement("K");

	// Cr
	G4Element* Cr = nistManager->FindOrBuildElement("Cr");

	// Fe
	G4Element* Fe = nistManager->FindOrBuildElement("Fe");

	// Ni
	G4Element* Ni = nistManager->FindOrBuildElement("Ni");

	// Sb
	G4Element* Sb = nistManager->FindOrBuildElement("Sb");

	// Ti
	G4Element* Ti = nistManager->FindOrBuildElement("Ti");

	// Ar
	G4Element* Ar = nistManager->FindOrBuildElement("Ar");

	// Cl
	G4Element* Cl = nistManager->FindOrBuildElement("Cl");

	//F
	G4Element* F = nistManager->FindOrBuildElement("F");

	// Na
	G4Element* Na = nistManager->FindOrBuildElement("Na");

	// Mg
	G4Element* Mg = nistManager->FindOrBuildElement("Mg");

	// Ca
	G4Element* Ca = nistManager->FindOrBuildElement("Ca");

	// Ce
	G4Element* Ce = nistManager->FindOrBuildElement("Ce");

	// La
	G4Element* La = nistManager->FindOrBuildElement("La");

	// Cu
	G4Element* Cu = nistManager->FindOrBuildElement("Cu");

	// Cs
	G4Element* Cs = nistManager->FindOrBuildElement("Cs");

	// N
	G4Element* N = nistManager->FindOrBuildElement("N");

	// Ge
	G4Element* Ge = nistManager->FindOrBuildElement("Ge");

	// Br
	G4Element* Br = nistManager->FindOrBuildElement("Br");

	//
	// define simple materials
	//

	   // Al
	density = 2.700 * g / cm3;
	a = 26.98 * g / mole;
	G4Material* Aluminium = new G4Material(name = "Aluminium", z = 13., a, density);

	// Cu
	density = 8.960 * g / cm3;
	G4Material* Copper = new G4Material(name = "Copper", density, ncomponents = 1);
	Copper->AddElement(Cu, natoms = 1);


	// Stainless Steel
	density = 8.00 * g / cm3;
	G4Material* SST = new G4Material(name = "SST", density, ncomponents = 3);
	SST->AddElement(Cr, fractionmass = 8.0 * perCent);
	SST->AddElement(Fe, fractionmass = 74.0 * perCent);
	SST->AddElement(Ni, fractionmass = 18.0 * perCent);

	//Plastik  

	G4Material* Plastik = new G4Material(name = "Plastik", density = 0.900 * g / cm3, ncomponents = 2);
	Plastik->AddElement(C, natoms = 3);
	Plastik->AddElement(H, natoms = 6);

	//BC404 

	G4Material* BC404 = new G4Material(name = "BC404", density = 1.032 * g / cm3, ncomponents = 2);
	BC404->AddElement(C, natoms = 10);
	BC404->AddElement(H, natoms = 11);

	// PCB                  //Bakelit definition*****
	density = 1.45 * g / cm3;
	G4Material* PCB = new G4Material(name = "PCB", density, ncomponents = 3);
	PCB->AddElement(C, natoms = 9);
	PCB->AddElement(H, natoms = 9);
	PCB->AddElement(O, natoms = 1);

	// Peek
	density = 1.26 * g / cm3;
	G4Material* Peek = new G4Material(name = "Peek", density, ncomponents = 3);
	Peek->AddElement(C, natoms = 19);
	Peek->AddElement(H, natoms = 12);
	Peek->AddElement(O, natoms = 3);

	// Borosilicated glass
	density = 2.30 * g / cm3;
	G4Material* BGL = new G4Material(name = "BGL", density, ncomponents = 5);
	BGL->AddElement(B, fractionmass = 5.3 * perCent);
	BGL->AddElement(O, fractionmass = 54.9 * perCent);
	BGL->AddElement(Al, fractionmass = 1.7 * perCent);
	BGL->AddElement(Si, fractionmass = 32.7 * perCent);
	BGL->AddElement(K, fractionmass = 5.4 * perCent);

	// Mu-metal
	density = 8.58 * g / cm3;
	G4Material* MuMetal = new G4Material(name = "MuMetal", density, ncomponents = 4);
	MuMetal->AddElement(Cr, fractionmass = 2.0 * perCent);
	MuMetal->AddElement(Fe, fractionmass = 18.0 * perCent);
	MuMetal->AddElement(Ni, fractionmass = 75.0 * perCent);
	MuMetal->AddElement(Cu, fractionmass = 5.0 * perCent);

	// Photocathode Material 
	// (Bialkali K2CsSb)
	density = 2.00 * g / cm3;
	G4Material* K2CsSb = new G4Material(name = "K2CsSb", density, ncomponents = 3);
	K2CsSb->AddElement(K, natoms = 2);
	K2CsSb->AddElement(Cs, natoms = 1);
	K2CsSb->AddElement(Sb, natoms = 1);

	// Mylar
	density = 1.39 * g / cm3;
	G4Material* Mylar = new G4Material(name = "Mylar", density, ncomponents = 3);
	Mylar->AddElement(C, natoms = 10);
	Mylar->AddElement(H, natoms = 8);
	Mylar->AddElement(O, natoms = 4);

	//Fe2O3
	density = 5.24 * g / cm3;
	G4Material* Fe2O3 = new G4Material(name = "Fe2O3", density, ncomponents = 2);
	Fe2O3->AddElement(Fe, natoms = 2);
	Fe2O3->AddElement(O, natoms = 3);

	//Magnetic media
	density = 2.93 * g / cm3;
	G4Material* MagneticMedia = new G4Material(name = "MagneticMedia", density, ncomponents = 2);
	MagneticMedia->AddMaterial(Fe2O3, fractionmass = 40.0 * perCent);
	MagneticMedia->AddMaterial(Mylar, fractionmass = 60.0 * perCent);

	//EJ-200  for the beta plastic detector 
	density = 1.023 * g / cm3;
	G4Material* EJ200 = new G4Material(name = "EJ200", density, ncomponents = 2);
	EJ200->AddElement(H, fractionmass = 0.0847);
	EJ200->AddElement(C, fractionmass = 0.9153);

	//EJ-510  for the reflective paint 
	density = 0.013 * g / cm3;
	G4Material* EJ510 = new G4Material(name = "EJ510", density, ncomponents = 4);
	EJ510->AddElement(H, fractionmass = 0.0290);
	EJ510->AddElement(C, fractionmass = 0.1719);
	EJ510->AddElement(O, fractionmass = 0.3886);
	EJ510->AddElement(Ti, fractionmass = 0.4105);

	// ABS (Acrylonitrile Butadiene Styrene)
	density = 1.04 * g / cm3;
	G4Material* ABS = new G4Material(name = "ABS", density, ncomponents = 3);
	ABS->AddElement(C, natoms = 15);
	ABS->AddElement(H, natoms = 17);
	ABS->AddElement(N, natoms = 1);

	//POLYETHYLENE
	//density = 0.94*g/cm3;
	density = 0.4 * g / cm3;
	G4Material* PE = new G4Material(name = "Polyethylene  ", density, ncomponents = 2);
	PE->AddElement(H, fractionmass = 14.3711 * perCent);
	PE->AddElement(C, fractionmass = 85.6289 * perCent);

	//POLYVINYL CHLORIDE: PVC
	density = 1.3 * g / cm3;
	G4Material* PVC = new G4Material(name = "Polyvinylchloride  ", density, ncomponents = 3);
	PVC->AddElement(H, fractionmass = 4.8380 * perCent);
	PVC->AddElement(C, fractionmass = 38.4360 * perCent);
	PVC->AddElement(Cl, fractionmass = 56.7260 * perCent);


	G4Material* Silicon =
		new G4Material("Silicon", z = 14., a = 28.09 * g / mole, density = 2.330 * g / cm3);


	// Ge
	density = 5.323 * g / cm3;
	G4Material* Germanium = new G4Material(name = "Germanium", density, ncomponents = 1);
	Germanium->AddElement(Ge, natoms = 1);

	//polyethylene source 2017Bi
	density = 1.05 * g / cm3;
	G4Material* plastic_Bi = new G4Material(name = "plastic_Bi", density, ncomponents = 2);
	plastic_Bi->AddElement(H, fractionmass = 0.498);
	plastic_Bi->AddElement(C, fractionmass = 0.502);

	////////////////////////////////////////////////////////////////////////////////////////////  
	// Nitrogen
	density = 1.25053 * mg / cm3;
	G4Material* Nitrogen = new G4Material(name = "N2", density, ncomponents = 1, kStateGas, 296.15 * kelvin, 1 * atmosphere);
	Nitrogen->AddElement(N, 2);

	// Argon 
	density = 1.7836 * mg / cm3;
	G4Material* Argon = new G4Material(name = "Argon", density, ncomponents = 1, kStateGas, 296.15 * kelvin, 1 * atmosphere);
	Argon->AddElement(Ar, 1);

	// Helium 
	density = 0.1786 * mg / cm3;
	G4Material* Helium = new G4Material(name = "Helium", density, ncomponents = 1, kStateGas, 296.15 * kelvin, 1 * atmosphere);
	Helium->AddElement(He, 1);

	// CF4 
	density = 3.72 * mg / cm3;
	G4Material* CF4 = new G4Material(name = "CF4", density, ncomponents = 2, kStateGas, 296.15 * kelvin, 1 * atmosphere);
	CF4->AddElement(C, 1);
	CF4->AddElement(F, 4);

	// Define the Polystyrene material
	density = 1.06 * g / cm3;
	G4Material* Polystyrene = new G4Material(name = "Polystyrene", density, ncomponents = 2);
	Polystyrene->AddElement(C, 8);
	Polystyrene->AddElement(H, 8);

	//////////////////////////////////////////////////////////////////
	//Vacuum
	G4Material* Vacuum = nistManager->FindOrBuildMaterial("G4_Galactic");

	//Air
	G4Material* Air = nistManager->FindOrBuildMaterial("G4_AIR");

	//Water
	G4Material* water = nistManager->FindOrBuildMaterial("G4_WATER");

	//<--------------------------------------------------------------------------------------------------------------------------------->
	//<--------------------------------------------------------------World-------------------------------------------------------------->
	//<--------------------------------------------------------------------------------------------------------------------------------->

	G4double WorldSize = 100. * cm;

	G4Box* solidWorld = new G4Box("World", WorldSize / 2, WorldSize / 2, WorldSize / 2);
	G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, Vacuum, "World");  //VACUUM
	//G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld,Air,"World");       //AIR
	G4VPhysicalVolume* physiWorld = new G4PVPlacement(0, G4ThreeVector(), "World", logicWorld, NULL, false, 0);


	G4ThreeVector SLCJ_position = G4ThreeVector(0. * cm, 0. * cm, 0. * mm);


	// Define the dimensions of the mother volume
	G4double motherVolumeX = 50 * cm;
	G4double motherVolumeY = 50 * cm;
	G4double motherVolumeZ = 50 * cm;

	//Mother volumen
	G4Box* solidmother_SLCJ = new G4Box("mother_SLCJ", motherVolumeX / 2, motherVolumeY / 2, motherVolumeZ / 2);
	G4LogicalVolume* logicmother_SLCJ = new G4LogicalVolume(solidmother_SLCJ, Vacuum, "mother_SLCJ");
	G4ThreeVector position_mother_SLCJ = G4ThreeVector(0. * mm, 0. * cm, 0. * mm);
	G4VPhysicalVolume* physimother_SLCJ = new G4PVPlacement(0, position_mother_SLCJ, "mother_SLCJ", logicmother_SLCJ, physiWorld, false, 0);

	//<--------------------------------------------------------------------------------------------------------------------------------->
	//<--------------------------------------------------------Container geometry-------------------------------------------------------->
	//<--------------------------------------------------------------------------------------------------------------------------------->

	// Define external dimensions and materials
	const G4double
		wallThickness = 1.0 * mm,
		middleDepth = 46 * mm,
		topBelowEdgeWidth = 48.8 * mm,
		height = 26 * mm,
		edgeHeight = 5 * mm,
		bottomWidth = 46.6 * mm,
		widthBelowEdgeDifference = topBelowEdgeWidth - bottomWidth,
		widthDifference = height / (height - edgeHeight) * widthBelowEdgeDifference,
		topProtrusion = widthDifference / 2,
		topWidth = topBelowEdgeWidth + widthDifference,
		topEdgeWidth = 50.1 * mm,
		topEdgeLenght = 72 * mm,
		edgeMargin = topEdgeWidth - topBelowEdgeWidth,
		frontTopDepth = topEdgeLenght - edgeMargin,
		frontBottomDepth = 64 * mm,
		frontBottomWidth = 27 * mm,
		frontTopWidth = 22 * mm,
		//frontWall_Height_ = 21 * mm,
		frontMiddleHeight = 5 * mm,//height - frontWall_Height_ * std::sqrt(1 - std::pow(topProtrusion / height, 2));
		frontMiddleWidth = 21.4 * mm;

	auto rightToLeft = [](G4ThreeVector v) {
		return G4ThreeVector(-v.x(), v.y(), v.z());
	};

	// set of external vertices
	std::set<G4ThreeVector> externalContainerVerticesSet;

	{
		//those are right external vertices
		externalContainerVerticesSet.emplace(topWidth / 2, 0, height);
		externalContainerVerticesSet.emplace(bottomWidth / 2, topProtrusion, 0);
		G4ThreeVector u(topWidth / 2, middleDepth, height);
		externalContainerVerticesSet.insert(u);
		G4ThreeVector v(bottomWidth / 2, middleDepth, 0);
		externalContainerVerticesSet.insert(v);
		G4ThreeVector w(frontTopWidth / 2, frontTopDepth, height);
		externalContainerVerticesSet.insert(w);
		{ // because of measurment imprecision there is a little shift between points' real coordinates and those in the code. Two last points have x and y offset added to account for that.
			auto tempFacet = new G4TriangularFacet(u, v, w, ABSOLUTE);
			externalContainerVerticesSet.emplace(shiftPointToFacetPlaneXY(G4ThreeVector(frontBottomWidth / 2, frontBottomDepth, 0), tempFacet));
			externalContainerVerticesSet.emplace(shiftPointToFacetPlaneXY(G4ThreeVector(frontMiddleWidth / 2, frontTopDepth - topProtrusion * (height - frontMiddleHeight) / height, frontMiddleHeight), tempFacet));
		}
	}

	//this generates left external verices from right ones
	externalContainerVerticesSet.merge([&]()->std::set<G4ThreeVector> {
		std::set<G4ThreeVector> leftVertices;
		for (auto vertex : externalContainerVerticesSet) {
			leftVertices.insert(rightToLeft(vertex));
		}
		return leftVertices;
		}());

	std::vector<G4ThreeVector> externalContainerVerticesVector;
	std::transform(externalContainerVerticesSet.begin(), externalContainerVerticesSet.end(), std::back_inserter(externalContainerVerticesVector),
		[](const auto& element) { return element; });

	// Declare main body solid before adding facets
	G4TessellatedSolid* mainBodySolid = new G4TessellatedSolid("mainBodySolid");

	// vector of body facets
	std::vector<G4TriangularFacet*> externalContainerFacetsVector;

	// vertex offset to keep (0,0,0) point inside the container
	G4ThreeVector vertexOffset(0, -0 * cm, 0);

	// Loop over all possible triangular facets and add them to the tessellated solid
	for (size_t i = 0; i < externalContainerVerticesVector.size() - 2; i++) {
		for (size_t j = i + 1; j < externalContainerVerticesVector.size() - 1; j++) {
			for (size_t k = j + 1; k < externalContainerVerticesVector.size(); k++) {
				externalContainerFacetsVector.push_back(new G4TriangularFacet(
					externalContainerVerticesVector[i] + vertexOffset,
					externalContainerVerticesVector[j] + vertexOffset,
					externalContainerVerticesVector[k] + vertexOffset,
					ABSOLUTE));
			}
		}
	}


	// Function to check if a facet intersects with any other facet while not being in the same plane
	auto isIntersectingAndNotInSamePlane = [&](G4TriangularFacet* facet, std::vector<G4TriangularFacet*>& facets) {
		for (const auto& otherFacet : facets) {
			if (facet != otherFacet && !areFacetsInSamePlane(facet, otherFacet) && doFacetsIntersect(facet, otherFacet))
				return true;
		}
		return false;
	};

	std::cout << "Le size: " << externalContainerFacetsVector.size() << '\n';

	if (false) {

		// set of external vertices
		std::vector<G4ThreeVector> testVec;

		//those are right external vertices
		testVec.emplace_back(topWidth / 2, 0, height);
		testVec.emplace_back(-frontTopWidth / 2, frontTopDepth, height);
		testVec.emplace_back(frontTopWidth / 2, frontTopDepth, height);

		auto facet = new G4TriangularFacet(testVec[0], testVec[1], testVec[2], ABSOLUTE);
		std::cout << "My facet:\n";
		for (int i = 0; i < 3; i++) {
			std::cout << facet->GetVertex(i) << '\n';
		}

		int howManyTimesTrue = 0, howManyChecked = 0, howManyFacetsTrue = 0;
		for (auto otherFacet : externalContainerFacetsVector) {
			//	const auto counterCopy = howManyTimesTrue;
			//auto otherFacet = externalContainerFacetsVector[12];
					//std::cout << howManyTimesTrue++ << '\n';
			if (facet != otherFacet && !(areFacetsInSamePlane(facet, otherFacet) || shareEdge(facet, otherFacet) || shareVertex(facet, otherFacet)) && doFacetsIntersect(facet, otherFacet)) {
				std::cout << "Intersecting facet:\n";
				for (int i = 0; i < 3; i++) {
					std::cout << otherFacet->GetVertex(i) << '\n';
				}
				exit(0);
			}
			//	howManyFacetsTrue += (howManyTimesTrue != counterCopy);
		}
		//std::cout << "Le counter: " << howManyTimesTrue << '\t' << howManyChecked << '\t' << howManyFacetsTrue << '\n';
	}

	// Remove facets from the vector that intersect with at least one other facet while not being in the same plane
	std::erase_if(externalContainerFacetsVector,
		[&](G4TriangularFacet* facet) {
			return std::any_of(externalContainerFacetsVector.begin(), externalContainerFacetsVector.end(), [&](G4TriangularFacet* otherFacet) {
				return facet != otherFacet && !(areFacetsInSamePlane(facet, otherFacet) || shareEdge(facet, otherFacet) || shareVertex(facet, otherFacet)) && doFacetsIntersect(facet, otherFacet);
				});
		});

	std::cout << "Le size: " << externalContainerFacetsVector.size() << '\n';

	decltype(externalContainerFacetsVector) externalContainerFacetsVector2;
	auto temp = std::vector(externalContainerFacetsVector.rbegin(), externalContainerFacetsVector.rend());
	for (auto facet : temp) {
		if (std::none_of(externalContainerFacetsVector2.begin(), externalContainerFacetsVector2.end(), [&](G4TriangularFacet* otherFacet) { return areFacetsInSamePlane(facet, otherFacet) && doFacetsIntersect(facet, otherFacet); })) {
			externalContainerFacetsVector2.push_back(facet);
		}
	}

	std::cout << "Le size - removed overlaps: " << externalContainerFacetsVector2.size() << '\n';

	for (auto facet : externalContainerFacetsVector2) {
		mainBodySolid->AddFacet((G4VFacet*)facet);
	}

	// Set the solid as closed to ensure proper inside/outside determination
	mainBodySolid->SetSolidClosed(true);

	// Create a logical volume for mainBodySolid
	G4LogicalVolume* mainBodyLogical = new G4LogicalVolume(mainBodySolid, Polystyrene, "mainBodyLogical");

	// Create a physical placement for mainBodyLogical
	G4RotationMatrix* rotation = new G4RotationMatrix();
	G4ThreeVector translation(0, 0, 0);  // Set the translation vector as per your requirement
	G4Transform3D transform(*rotation, translation);

	G4VPhysicalVolume* mainBodyPhysical = new G4PVPlacement(transform, mainBodyLogical, "mainBodyPhysical", logicWorld, false, 0);


	G4double
		length = 10.0 * cm,
		width = 10.0 * cm,
		cornerLength = 1.0 * cm,
		neckRadius = 0.5 * cm,
		neckHeight = 1.0 * cm;

	// Define neck dimensions
	G4double
		neckLength = 22.0 * mm,
		neckOuterRadius = 10.0 * mm,
		neckInnerRadius = neckOuterRadius - wallThickness;

	// Create neck solid
	G4Tubs* neckSolid = new G4Tubs("tubeSolid", neckInnerRadius, neckOuterRadius, neckLength / 2.0, 0.0, CLHEP::twopi);

	// Define cap dimensions
	G4double
		capOuterRadius = 12.3 * mm,
		capHeight = 14.5 * mm,
		capInnerRadius = capOuterRadius - wallThickness,
		capNarrowingInnerRadius = 10.0 * mm,
		capNarrowingHeight = 1.0 * mm;

	// Create solid for cap main body
	G4Tubs* capBodySolid = new G4Tubs("capBodySolid", 0.0, capOuterRadius, capHeight / 2.0, 0.0, CLHEP::twopi);

	// Create solid for cap narrowing
	G4Cons* capNarrowingSolid = new G4Cons("capNarrowingSolid", 0.0, capInnerRadius, 0.0, capNarrowingInnerRadius, capNarrowingHeight / 2.0, 0.0, CLHEP::twopi);

	// Create union solid for cap
	G4UnionSolid* capSolid = new G4UnionSolid("capSolid", capBodySolid, capNarrowingSolid, 0, G4ThreeVector(0.0, 0.0, capHeight / 2.0 - capNarrowingHeight / 2.0));

	//enum Location {
	//	Internal = -1,
	//	External = 0
	//};
	//enum Height {
	//	Top = 1,
	//	Bottom = -1,
	//	Cuttoff = 0,
	//};
	//enum Depth {
	//	Back = -1,
	//	Middle = 0,
	//	Front = 1
	//};
	//enum Width {
	//	Left = -1,
	//	Right = 1
	//};
	//
	//// Define the vertices of the container
	//// internal/external - if the vertex is on the inside or the outside
	//// top/bottom - top is the direction towards which bottleneck is angled
	//// back/middle/front - if the vertex is on the opposite side to bottlneck, between back and front or on the side of bottleneck
	//// left/right - determined based on looking from the back in the direction of bottleneck with it pointing upwards
	//G4ThreeVector externalTopBackRightVertex(25.0 * mm, 25.0 * mm, 15.0 * mm);
	//G4ThreeVector externalTopBackLeftVertex(-25.0 * mm, 25.0 * mm, 15.0 * mm);
	//G4ThreeVector externalTopFrontLeftVertex(-25.0 * mm, 25.0 * mm, -15.0 * mm);
	//G4ThreeVector externalTopFrontRightVertex(25.0 * mm, 25.0 * mm, -15.0 * mm);
	//G4ThreeVector internalTopBackRightVertex(23.0 * mm, 23.0 * mm, 13.0 * mm);
	//G4ThreeVector internalTopBackLeftVertex(-23.0 * mm, 23.0 * mm, 13.0 * mm);
	//G4ThreeVector internalTopFrontLeftVertex(-23.0 * mm, 23.0 * mm, -13.0 * mm);
	//G4ThreeVector internalTopFrontRightVertex(23.0 * mm, 23.0 * mm, -13.0 * mm);
	//G4ThreeVector externalBottomBackRightVertex(25.0 * mm, -25.0 * mm, 15.0 * mm);
	//G4ThreeVector externalBottomBackLeftVertex(-25.0 * mm, -25.0 * mm, 15.0 * mm);
	//G4ThreeVector externalBottomFrontLeftVertex(-25.0 * mm, -25.0 * mm, -15.0 * mm);
	//G4ThreeVector externalBottomFrontRightVertex(25.0 * mm, -25.0 * mm, -15.0 * mm);
	//G4ThreeVector internalBottomBackRightVertex(23.0 * mm, -23.0 * mm, 13.0 * mm);
	//G4ThreeVector internalBottomBackLeftVertex(-23.0 * mm, -23.0 * mm, 13.0 * mm);
	//G4ThreeVector internalBottomFrontLeftVertex(-23.0 * mm, -23.0 * mm, -13.0 * mm);
	//G4ThreeVector internalBottomFrontRightVertex(23.0 * mm, -23.0 * mm, -13.0 * mm);
	//G4ThreeVector externalBottleneckTopVertex(0.0 * mm, 35.0 * mm, 0.0 * mm);
	//G4ThreeVector externalBottleneckBottomVertex(0.0 * mm, -25.0 * mm, 0.0 * mm);
	//G4ThreeVector internalBottleneckTopVertex(0.0 * mm, 33.0 * mm, 0.0 * mm);
	//G4ThreeVector internalBottleneckBottomVertex(0.0 * mm, -23.0 * mm, 0.0 * mm);
	//
	//
	//G4ThreeVector v0(0, 0, 0);
	//G4ThreeVector v1(1, 0, 0);
	//G4ThreeVector v2(0.5, std::sqrt(3) / 2, 0);
	//G4ThreeVector v3(0.5, std::sqrt(3) / 6, std::sqrt(2. / 3.));
	//
	//// Define the triangular facets of the container
	//G4TriangularFacet* facet1 = new G4TriangularFacet(v0, v1, v2, ABSOLUTE);
	//G4TriangularFacet* facet2 = new G4TriangularFacet(v0, v3, v1, ABSOLUTE);
	//G4TriangularFacet* facet3 = new G4TriangularFacet(v1, v3, v2, ABSOLUTE);
	//G4TriangularFacet* facet4 = new G4TriangularFacet(v2, v3, v0, ABSOLUTE);
	//
	//// Add the facets to the tessellated solid
	//mainBodySolid->AddFacet(facet1);
	//mainBodySolid->AddFacet(facet2);
	//mainBodySolid->AddFacet(facet3);
	//mainBodySolid->AddFacet(facet4);
	//
	//// Create a logical volume and place it in the simulation
	//G4LogicalVolume* tetrahedronLV = new G4LogicalVolume(mainBodySolid, Air, "TetrahedronLV");
	//G4VPhysicalVolume* tetrahedronPV = new G4PVPlacement(0, G4ThreeVector(), tetrahedronLV, "TetrahedronPV", logicWorld, false, 0);
	//
	//
	//
	//bool checkOverlaps = true;
	//
	//// Create shapes
	//G4Box* cuboid = new G4Box("cuboid", length / 2, width / 2, height / 2);
	//
	//G4Tubs* bottleNeck = new G4Tubs("bottleNeck", 0, neckRadius, neckHeight / 2, 0, 360.0 * deg);
	//
	//G4SubtractionSolid* cutCuboid1 = new G4SubtractionSolid("cutCuboid1", cuboid, bottleNeck);
	//
	//G4ThreeVector translate1(length / 2 - cornerLength, width / 2 - cornerLength, 0);
	//
	//G4SubtractionSolid* cutCuboid2 = new G4SubtractionSolid("cutCuboid2", cutCuboid1, bottleNeck, 0, translate1);
	//
	//G4ThreeVector translate2(-length / 2 + cornerLength, width / 2 - cornerLength, 0);
	//
	//G4SubtractionSolid* cutCuboid3 = new G4SubtractionSolid("cutCuboid3", cutCuboid2, bottleNeck, 0, translate2);
	//
	//G4ThreeVector translate3(-length / 2 + cornerLength, -width / 2 + cornerLength, 0);
	//
	//G4SubtractionSolid* cutCuboid4 = new G4SubtractionSolid("cutCuboid4", cutCuboid3, bottleNeck, 0, translate3);
	//
	//// Create logical volumes
	//G4LogicalVolume* logicalShape = new G4LogicalVolume(cutCuboid4, Air, "logicalShape");
	//
	//G4ThreeVector translate4(0, 0, (height - neckHeight) / 2);
	//
	////new G4PVPlacement(0, translate4, "bottleNeck", logicalShape, bottleNeck, false, 0);
	//
	//G4double insideLength = length - 2 * wallThickness;
	//G4double insideWidth = width - 2 * wallThickness;
	//G4double insideHeight = height - neckHeight - 2 * wallThickness - 2 * cornerLength;
	//
	//G4Box* insideShape = new G4Box("insideShape", insideLength / 2, insideWidth / 2, insideHeight / 2);
	//
	//G4SubtractionSolid* finalShape = new G4SubtractionSolid("finalShape", logicalShape->GetSolid(), insideShape);
	//
	//G4LogicalVolume* logicalInside = new G4LogicalVolume(insideShape, Air, "logicalInside");
	//
	//G4ThreeVector translate5(0, 0, -insideHeight / 2);
	//
	////new G4PVPlacement(0, translate5, "inside", logicalInside, finalShape, false, 0, checkOverlaps);





	////<-------------------------------------------------------Old shit------------------------------------------------------------------>
	//
	//// Define the dimensions of the active volume
	//G4double activeVolumeX = 20 * cm;
	//G4double activeVolumeY = 33 * cm;
	//G4double activeVolumeZ = 21 * cm;
	//
	//// Define the dimensions of the external volume
	//G4double externalVolumeX = 22.2 * cm;
	//G4double externalVolumeY = 34.2 * cm;
	//G4double externalVolumeZ = 22 * cm;
	//
	//// Create a solid volume for the active volume
	//G4Box* activeVolumeSolid = new G4Box("activeVolumeSolid", activeVolumeX / 2, activeVolumeY / 2, activeVolumeZ / 2);
	//
	//// Create a solid volume for the external volume
	//G4Box* externalVolumeSolid = new G4Box("externalVolumeSolid", externalVolumeX / 2, externalVolumeY / 2, externalVolumeZ / 2);
	//
	//// Subtract the active volume from the external volume to create detector walls
	//G4SubtractionSolid* wallsSolid = new G4SubtractionSolid("wallsSolid", externalVolumeSolid, activeVolumeSolid /*, 0, activeVolumeShift*/);
	//
	//// Assign the Stesalit material to the walls of the detector
	//G4LogicalVolume* wallsLogical = new G4LogicalVolume(wallsSolid, Polystyrene, "wallsLogical");
	//
	//// Assign the GasSLCJ to the active volume of the detector
	//G4LogicalVolume* activeVolumeLogical = new G4LogicalVolume(activeVolumeSolid, Air, "activeVolumeLogical");
	//
	//G4VPhysicalVolume* physiSLCJ = new G4PVPlacement(0, SLCJ_position, "SLCJ_activeVolume", activeVolumeLogical, physimother_SLCJ, false, 0);
	//G4VPhysicalVolume* physiWallsSLCJ = new G4PVPlacement(0, SLCJ_position, "SLCJ_walls", wallsLogical, physimother_SLCJ, false, 0);




	// Construct the field creator - this will register the field it creates
	//F02ElectricFieldSetup* fieldSetup = new F02ElectricFieldSetup();



	//********************************************************************************************************************/


	//
	/// visualization attributes
	//


	logicmother_SLCJ->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());


	// red color for detector 
	G4VisAttributes* Att2 = new G4VisAttributes(G4Colour::Red());
	// blue for gas
	G4VisAttributes* Att3 = new G4VisAttributes(G4Colour::Blue());
	// grey for PMT
	G4VisAttributes* Att4 = new G4VisAttributes(G4Colour::Grey());
	//white
	G4VisAttributes* Att5 = new G4VisAttributes(G4Colour::White());
	//magenta
	G4VisAttributes* Att6 = new G4VisAttributes(G4Colour::Magenta());
	// black
	G4VisAttributes* Att_black = new G4VisAttributes(G4Colour::Black());
	// green
	G4VisAttributes* Att8 = new G4VisAttributes(G4Colour::Green());
	// orange
	G4VisAttributes* Att_orange = new G4VisAttributes(G4Colour(1, 0.5, 0, 0.8));
	// dark green
	G4VisAttributes* Att10 = new G4VisAttributes(G4Colour(0.5, 1.0, 0.0));


	G4VisAttributes* Att_pale_yellow = new G4VisAttributes(G4Colour(255. / 255, 204. / 255, 102. / 255));
	G4VisAttributes* Att_light_grey = new G4VisAttributes(G4Colour(224. / 255, 224. / 255, 224. / 255));
	G4VisAttributes* Att_light_blue = new G4VisAttributes(G4Colour(145. / 255, 202. / 255, 249. / 255));
	//G4VisAttributes* Att_pale_yellow = new G4VisAttributes(G4Colour(255./255, 204./255, 102./255));



	//SLCJLogicalVolume->SetVisAttributes(Att9);
	//SLCJLogicalVolume->SetVisAttributes(G4VisAttributes::GetInvisible());
	//activeVolumeLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
	//wallsLogical->SetVisAttributes(Att_pale_yellow);
	mainBodyLogical->SetVisAttributes(Att_pale_yellow);

	//
	// always return the physical World
	//

	return physiWorld;
}

const G4double SLCJDetectorConstruction::getFoodVolume() {
	return foodVolume;
}

void SLCJDetectorConstruction::saveDetails(std::filesystem::path p) {

	std::cout << p << _endl_;
	// Open file for writing physical volumes
	std::ofstream outfile;
	outfile.open(p / "volumes.csv");
	outfile << "Physical Volume,Material,Shape\n";

	checkpoint;
	// Loop over all the physical volumes
	auto& physicalVolumeStore = *G4PhysicalVolumeStore::GetInstance();
	for (const auto& physicalVolume : physicalVolumeStore) {
		const auto& material = *physicalVolume->GetLogicalVolume()->GetMaterial();
		const auto& shape = *physicalVolume->GetLogicalVolume()->GetSolid();
		outfile << std::format("{},{},{}\n",
			physicalVolume->GetName(), material.GetName(), shape.GetName());
	}
	outfile.close();

	checkpoint;
	// Open file for writing materials
	outfile.open(p / "materials.csv");
	outfile << "Material,Density (g/cm3)\n";

	// Loop over all the materials
	auto& materialTable = *G4Material::GetMaterialTable();
	for (const auto& material : materialTable) {
		outfile << std::format("{},{:.3f}\n",
			material->GetName(), material->GetDensity() / (g / cm3));
	}
	outfile.close();

	checkpoint;
	// Open file for writing elements
	outfile.open(p / "elements.csv");
	outfile << "Element,Atomic Number\n";
	std::ofstream outfile_isotopes(p / "isotopes.csv");
	outfile_isotopes << "Isotope,Atomic Number,Atomic Mass (g/mol),Natural Abundance\n";
	// Loop over all the elements
	auto& elementTable = *G4Element::GetElementTable();
	for (const auto& element : elementTable) {
		outfile << std::format("{},{}\n",
			element->GetName(), element->GetZ());
		auto& isotopeVector = *element->GetIsotopeVector();
		auto relativeAbundanceVector = element->GetRelativeAbundanceVector();
		// Loop over all the isotopes
		for (const auto& isotope : isotopeVector) {
			outfile_isotopes << std::format("{},{},{:.3f},{:.3f}\n",
				isotope->GetName(), isotope->GetZ(), isotope->GetA() / (g / mole),
				relativeAbundanceVector[isotope->GetIndex()]);
		}
	}
	outfile.close();
	outfile_isotopes.close();

}
