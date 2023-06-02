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
#include "G4RayTracer.hh"
#include <fstream>
#include <format>
#include <map>
#include <tuple>
#include <set>
#include <boost/math/tools/roots.hpp>

#include "G4IntersectionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"
#include <G4UnionSolid.hh>
//#include "TMath.h"
#include "G4Polycone.hh"
#include "geometry.h"

#include "F02ElectricFieldSetup.hh"

// Global file with the events
std::ifstream geoInputFile;

inline std::string filename_string(std::string path_str) {
	return path_str.substr(path_str.rfind("\\") + 1, path_str.size() - path_str.rfind("\\") - 1);
};

#define _endl_ " (" << filename_string(__FILE__) << "; " << __LINE__ << ")" << '\n'
#define checkpoint std::cout << "checkpoint" << _endl_

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

	// vertex offset to keep (0,0,0) point inside the container
	G4ThreeVector vertexOffset(0, -frontBottomDepth / 2, -2 * mm);

	//<--------------------------------------------------------Bottle neck solid construction-------------------------------------------->

	// Define neck dimensions
	G4double
		neckLength = 22.0 * mm,
		neckOuterRadius = 10.0 * mm,
		neckInnerRadius = neckOuterRadius - wallThickness;

	// Create neck solid
	G4Tubs* neckSolid = new G4Tubs("neckSolid", 0, neckOuterRadius, neckLength / 2.0, 0.0, CLHEP::twopi);

	// Create neck solid
	G4Tubs* neckInternalSolid = new G4Tubs("neckSolid", 0, neckInnerRadius, neckLength / 2.0, 0.0, CLHEP::twopi);

	// Define cap dimensions
	G4double
		capOuterRadius = 12.3 * mm,
		capHeight = 14.5 * mm,
		capInnerRadius = capOuterRadius - wallThickness,
		capNarrowingInnerRadius = 10.0 * mm,
		capNarrowingHeight = 1.0 * mm;

	// Create solid for cap main body
	G4Tubs* capBodySolid = new G4Tubs("capBodySolid", capInnerRadius, capOuterRadius, capHeight / 2.0, 0.0, CLHEP::twopi);

	// Create solid for cap main body
	G4Tubs* capClosingSolid = new G4Tubs("capClosingSolid", 0, capInnerRadius, wallThickness / 2.0, 0.0, CLHEP::twopi);

	// Create solid for cap narrowing
	G4Cons* capNarrowingSolid = new G4Cons("capNarrowingSolid", 0.0, capInnerRadius, 0.0, capNarrowingInnerRadius, capNarrowingHeight / 2.0, 0.0, CLHEP::twopi);

	// Create union solid for cap
	G4UnionSolid* capOpenSolid = new G4UnionSolid("capOpenSolid", capBodySolid, capNarrowingSolid, 0, G4ThreeVector(0.0, 0.0, capHeight / 2.0 - capNarrowingHeight / 2.0));

	// Create union solid for cap
	G4UnionSolid* capSolid = new G4UnionSolid("capSolid", capBodySolid, capClosingSolid, 0, G4ThreeVector(0.0, 0.0, capHeight / 2 - wallThickness / 2.0));


	// Create a logical volume for capSolid
	G4LogicalVolume* capLogical = new G4LogicalVolume(capSolid, Polystyrene, "capLogical");

	// define neck position
	const G4double
		neckEndBottomHeight = 10 * mm,
		neckHeightDifference = neckEndBottomHeight - frontMiddleHeight,
		neckAngleSine = neckHeightDifference / neckLength,
		neckAngle = std::asin(neckAngleSine);

	// Create a physical placement for neckLogical
	G4RotationMatrix* neckRotation = new G4RotationMatrix();
	neckRotation->rotateX(90 * deg + neckAngle);
	G4RotationMatrix* capRotation = new G4RotationMatrix();
	capRotation->rotateX(-90 * deg + neckAngle);
	G4ThreeVector neckTranslation(0, topEdgeLenght + neckLength / 2 - neckOuterRadius * neckAngleSine + 3 * mm, frontMiddleHeight + (height - frontMiddleHeight) / 2 + 3 * mm);
	G4ThreeVector capTranslation(0, topEdgeLenght + neckLength / 2 - neckOuterRadius * neckAngleSine - 3 * mm + (neckLength - capHeight + 2 * wallThickness) / 2 * std::cos(neckAngle), frontMiddleHeight + (height - frontMiddleHeight) / 2 + 3 * mm + (neckLength - capHeight + 2 * wallThickness) / 2 * neckAngleSine);
	G4Transform3D neckTransform(*neckRotation, neckTranslation + vertexOffset);
	G4Transform3D capTransform(*capRotation, capTranslation);

	G4VPhysicalVolume* capPhysical = new G4PVPlacement(capTransform, capLogical, "capPhysical", logicWorld, false, 0);

	//<--------------------------------------------------------External solid construction----------------------------------------------->

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
	G4TessellatedSolid* mainBodyExternalSolid = new G4TessellatedSolid("mainBodyExternalSolid");

	// vector of body facets
	std::vector<G4TriangularFacet*> externalContainerConstructionFacetsVector;

	for (auto& vertex : externalContainerVerticesVector) {
		vertex += vertexOffset;
	}

	// Loop over all possible triangular facets and add them to the tessellated solid
	for (size_t i = 0; i < externalContainerVerticesVector.size() - 2; i++) {
		for (size_t j = i + 1; j < externalContainerVerticesVector.size() - 1; j++) {
			for (size_t k = j + 1; k < externalContainerVerticesVector.size(); k++) {
				externalContainerConstructionFacetsVector.push_back(new G4TriangularFacet(
					externalContainerVerticesVector[i],
					externalContainerVerticesVector[j],
					externalContainerVerticesVector[k],
					ABSOLUTE));
			}
		}
	}

	// Remove facets from the vector that intersect with at least one other facet while not being in the same plane
	std::erase_if(externalContainerConstructionFacetsVector,
		[&](G4TriangularFacet* facet) {
			return std::any_of(externalContainerConstructionFacetsVector.begin(), externalContainerConstructionFacetsVector.end(), [&](G4TriangularFacet* otherFacet) {
				return facet != otherFacet && !(areFacetsCoplanar(facet, otherFacet) || shareEdge(facet, otherFacet) || shareVertex(facet, otherFacet)) && doFacetsIntersect(facet, otherFacet);
				});
		});

	auto groups = groupCoplanarFacets(externalContainerConstructionFacetsVector);

	decltype(externalContainerConstructionFacetsVector) externalContainerFacetsVector;

	for (auto& group : groups) {
		std::set<G4ThreeVector> vertices;
		for (auto facet : group) {
			auto verticesFacet = getVertices(facet);
			for (auto vertex : verticesFacet) {
				vertices.insert(vertex);
			}
		}
		auto facets = spanFacets(vertices);
		for (auto facet : facets) {
			externalContainerFacetsVector.push_back(facet);
			mainBodyExternalSolid->AddFacet((G4VFacet*)facet);
		}
	}

	// Set the solid as closed to ensure proper inside/outside determination
	mainBodyExternalSolid->SetSolidClosed(true);

	//<--------------------------------------------------------Internal solid construction----------------------------------------------->

	decltype(externalContainerConstructionFacetsVector) internalContainerFacetsVector;

	auto mapVertexNormals = getVertexFacetsMap(externalContainerVerticesVector, groups);

	std::map<G4ThreeVector, G4ThreeVector> pointsTranslated;

	for (auto vertexNormals : mapVertexNormals) {
		pointsTranslated[vertexNormals.first] = calculateInternalPoint(vertexNormals);
	}
	checkpoint;
	G4TessellatedSolid* mainBodyInternalSolid = new G4TessellatedSolid("mainBodyInternalSolid");

	for (auto facet : externalContainerFacetsVector) {
		auto vertices = getVertices(facet);
		mainBodyInternalSolid->AddFacet(new G4TriangularFacet(pointsTranslated[vertices[0]], pointsTranslated[vertices[1]], pointsTranslated[vertices[2]], ABSOLUTE));
	}

	// Set the solid as closed to ensure proper inside/outside determination
	mainBodyInternalSolid->SetSolidClosed(true);


	auto mainBodyExternalSolidWithNeck = new G4UnionSolid("mainBodyExternalSolidWithNeck", mainBodyExternalSolid, neckSolid, neckTransform);

	auto mainBodyInternalSolidWithNeck = new G4UnionSolid("mainBodyExternalSolidWithNeck", mainBodyInternalSolid, neckInternalSolid, neckTransform);

	//std::cout << "Volume: " << mainBodyExternalSolidWithNeck->GetCubicVolume() / cm3 << "cm3\n";

	// Create a logical volume for mainBodySolid
	G4LogicalVolume* mainBodyInternalLogical = new G4LogicalVolume(mainBodyInternalSolidWithNeck, Air, "mainBodyInternalLogical");
	// Create a logical volume for mainBodySolid
	G4LogicalVolume* mainBodyExternalLogical = new G4LogicalVolume(mainBodyExternalSolidWithNeck, Polystyrene, "mainBodyExternalLogical");

	// Create a physical placement for mainBodyLogical
	G4RotationMatrix* rotation = new G4RotationMatrix();
	G4ThreeVector translation = -vertexOffset;  // Set the translation vector as per your requirement
	G4Transform3D transform(*rotation, translation);

	G4VPhysicalVolume* mainBodyPhysical = new G4PVPlacement(rotation, -vertexOffset, mainBodyExternalLogical, "mainBodyExternalPhysical", logicWorld, false, 0);

	G4VPhysicalVolume* mainBodyInternalPhysical = new G4PVPlacement(rotation, G4ThreeVector(), mainBodyInternalLogical, "mainBodyInternalPhysical", mainBodyExternalLogical, false, 0);

	//<--------------------------------------------------------Cells and food construction----------------------------------------------->

	const G4double cellsLayerThickness = 2 * um;

	auto cellsBaseSolid = new G4Box("cellsBaseSolid", topEdgeWidth / 2, topEdgeLenght / 2, cellsLayerThickness / 2);

	auto cellsSolid = new G4UnionSolid("cellsSolid", mainBodyInternalSolidWithNeck, cellsBaseSolid, nullptr, G4ThreeVector(0, 0, -2*mm));

	// Create a logical volume for mainBodySolid
	G4LogicalVolume* cellsLogical = new G4LogicalVolume(cellsSolid, water, "cellsLogical");

	G4VPhysicalVolume* cellsPhysical = new G4PVPlacement(nullptr, -vertexOffset - G4ThreeVector(0,0,0)*mm, cellsLogical, "cellsPhysical", logicWorld, false, 0);
	// Construct the field creator - this will register the field it creates
	//F02ElectricFieldSetup* fieldSetup = new F02ElectricFieldSetup();



	//********************************************************************************************************************/


	//
	/// visualization attributes
	//


	logicmother_SLCJ->SetVisAttributes(G4VisAttributes::GetInvisible());
	logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());


	// red color for detector 
	G4VisAttributes* Att_red = new G4VisAttributes(G4Colour::Red());
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
	G4VisAttributes* Att_green = new G4VisAttributes(G4Colour::Green());
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
	mainBodyInternalLogical->SetVisAttributes(Att_light_blue);
	mainBodyExternalLogical->SetVisAttributes(Att_pale_yellow);
	capLogical->SetVisAttributes(Att_green);
	cellsLogical->SetVisAttributes(Att_red);

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
