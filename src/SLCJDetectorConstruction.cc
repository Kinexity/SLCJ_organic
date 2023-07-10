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
#include <functional>

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


	//
	// define Elements
	//

	// H
	G4Element* H = nistManager->FindOrBuildElement("H");

	// C
	G4Element* C = nistManager->FindOrBuildElement("C");

	// O
	G4Element* O = nistManager->FindOrBuildElement("O");

	// Ti
	G4Element* Ti = nistManager->FindOrBuildElement("Ti");

	//
	// define simple materials
	//

	// Define the Polystyrene material
	density = 1.06 * g / cm3;
	G4Material* Polystyrene = new G4Material(name = "Polystyrene", density, ncomponents = 2);
	Polystyrene->AddElement(C, 8);
	Polystyrene->AddElement(H, 8);


	// Define the PMMA material
	density = 1.18 * g / cm3;
	G4Material* PMMA = new G4Material(name = "PMMA", density, ncomponents = 3);
	PMMA->AddElement(C, 5);
	PMMA->AddElement(H, 8);
	PMMA->AddElement(O, 2);


	// Create the RW3 slab phantom material
	density = 1.045 * g / cm3;
	G4Material* RW3PhantomMaterial = new G4Material("RW3PhantomMaterial", density, 4);
	RW3PhantomMaterial->AddElement(H, 7.59 * perCent);
	RW3PhantomMaterial->AddElement(C, 90.41 * perCent);
	RW3PhantomMaterial->AddElement(O, 0.80 * perCent);
	RW3PhantomMaterial->AddElement(Ti, 1.2 * perCent);

	//////////////////////////////////////////////////////////////////
	//Vacuum
	G4Material* Vacuum = nistManager->FindOrBuildMaterial("G4_Galactic");

	//Air
	G4Material* Air = nistManager->FindOrBuildMaterial("G4_AIR");

	//Water
	G4Material* Water = nistManager->FindOrBuildMaterial("G4_WATER");

	//<--------------------------------------------------------------------------------------------------------------------------------->
	//<--------------------------------------------------------------World-------------------------------------------------------------->
	//<--------------------------------------------------------------------------------------------------------------------------------->

	//G4double WorldSize = 20. * cm;

	// Define the dimensions of the mother volume
	G4double WorldSizeX = 10 * cm;
	G4double WorldSizeY = 20 * cm;
	G4double WorldSizeZ = 10 * cm;

	G4Box* solidWorld = new G4Box("World", WorldSizeX / 2, WorldSizeY / 2, WorldSizeZ / 2);
	G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, Vacuum, "World");  //VACUUM
	G4VPhysicalVolume* physiWorld = new G4PVPlacement(0, G4ThreeVector(), "World", logicWorld, NULL, 0, true);


	G4ThreeVector SLCJ_position = G4ThreeVector();


	// Define the dimensions of the mother volume
	G4double motherVolumeX = 10 * cm;
	G4double motherVolumeY = 20 * cm;
	G4double motherVolumeZ = 10 * cm;

	//Mother volume
	G4Box* solidmother_SLCJ = new G4Box("mother_SLCJ", motherVolumeX / 2, motherVolumeY / 2, motherVolumeZ / 2);
	G4LogicalVolume* logicmother_SLCJ = new G4LogicalVolume(solidmother_SLCJ, Vacuum, "mother_SLCJ");
	G4ThreeVector position_mother_SLCJ = G4ThreeVector(0. * mm, 0. * cm, 0. * mm);
	G4VPhysicalVolume* physimother_SLCJ = new G4PVPlacement(0, position_mother_SLCJ, "mother_SLCJ", logicmother_SLCJ, physiWorld, 0, true);

	//<--------------------------------------------------------------------------------------------------------------------------------->
	//<--------------------------------------------------------Container geometry------------------------------------------------------->
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
		frontMiddleHeight = 5 * mm,
		frontMiddleWidth = 21.4 * mm;

	// vertex offset to keep (0,0,0) point inside the container
	G4ThreeVector vertexOffset(0, -frontBottomDepth / 2, -2 * mm);

	//<--------------------------------------------------------Bottle neck solid construction-------------------------------------------->

	// Define neck dimensions
	G4double
		neckLength = 24.0 * mm,
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
	G4ThreeVector neckTranslation(0, topEdgeLenght + neckLength / 2 - neckOuterRadius * neckAngleSine - 4 * mm, (height + frontMiddleHeight) / 2 + 2 * mm);
	G4ThreeVector capTranslation = neckTranslation + (neckLength - capHeight + 2 * wallThickness) / 2 * G4ThreeVector(0, std::cos(neckAngle), neckAngleSine);
	G4Transform3D neckTransform(*neckRotation, neckTranslation + vertexOffset);
	G4Transform3D capTransform(*capRotation, capTranslation);

	G4VPhysicalVolume* capPhysical = new G4PVPlacement(capTransform, "capPhysical", capLogical, physimother_SLCJ, false, 0);

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
			auto normal = facet->GetSurfaceNormal();
			if (normal * facet->GetVertex(0) < 0) {
				facet->SetSurfaceNormal(-normal);
			}
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

	G4TessellatedSolid* mainBodyInternalSolid = new G4TessellatedSolid("mainBodyInternalSolid");

	for (auto facet : externalContainerFacetsVector) {
		auto vertices = getVertices(facet);
		auto newFacet = new G4TriangularFacet(pointsTranslated[vertices[0]], pointsTranslated[vertices[1]], pointsTranslated[vertices[2]], ABSOLUTE);
		auto normal = newFacet->GetSurfaceNormal();
		if (normal * newFacet->GetVertex(0) < 0) {
			newFacet->SetSurfaceNormal(-normal);
		}
		mainBodyInternalSolid->AddFacet(newFacet);
	}

	// Set the solid as closed to ensure proper inside/outside determination
	mainBodyInternalSolid->SetSolidClosed(true);


	auto mainBodyExternalSolidWithNeck = new G4UnionSolid("mainBodyExternalSolidWithNeck", mainBodyExternalSolid, neckSolid, neckTransform);

	auto mainBodyInternalSolidWithNeck = new G4UnionSolid("mainBodyInternalSolidWithNeck", mainBodyInternalSolid, neckInternalSolid, neckTransform);


	// Create a logical volume for mainBodySolid
	G4LogicalVolume* mainBodyInternalLogical = new G4LogicalVolume(mainBodyInternalSolidWithNeck, Air, "mainBodyInternalLogical");
	// Create a logical volume for mainBodySolid
	G4LogicalVolume* mainBodyExternalLogical = new G4LogicalVolume(mainBodyExternalSolidWithNeck, Polystyrene, "mainBodyExternalLogical");

	// Create a physical placement for mainBodyLogical
	G4RotationMatrix* rotation = new G4RotationMatrix();
	G4ThreeVector translation = -vertexOffset;  // Set the translation vector as per your requirement
	G4Transform3D transform(*rotation, translation);

	G4VPhysicalVolume* mainBodyPhysical = new G4PVPlacement(rotation, -vertexOffset, "mainBodyExternalPhysical", mainBodyExternalLogical, physimother_SLCJ, false, 0);

	G4ThreeVector mainBodyInternalPosition(0, 0, 0);
	G4VPhysicalVolume* mainBodyInternalPhysical = new G4PVPlacement(rotation, mainBodyInternalPosition, "mainBodyInternalPhysical", mainBodyInternalLogical, mainBodyPhysical, false, 0);

	//<--------------------------------------------------------Cells and food construction----------------------------------------------->

	const G4double cellsLayerThickness = 2 * um;
	G4double foodHeight = 3 * mm; //initial value

	auto cellsBaseSolid = new G4Box("cellsBaseSolid", topEdgeWidth / 2, topEdgeLenght / 2, cellsLayerThickness / 2);

	G4ThreeVector cellsSolidPositionOffset = G4ThreeVector(0, 0, cellsLayerThickness / 2 - 1 * mm);

	auto cellsSolid = new G4IntersectionSolid("cellsSolid", mainBodyInternalSolidWithNeck, cellsBaseSolid, nullptr, cellsSolidPositionOffset);

	// food height -> food volume function
	std::function fn = [&](G4double fHeight)->G4double {
		G4ThreeVector foodSolidPositionOffset = cellsSolidPositionOffset + G4ThreeVector(0, 0, cellsLayerThickness / 2 + fHeight / 2);
		std::unique_ptr<G4Box> foodBaseSolid = std::make_unique<G4Box>("foodBaseSolid", topEdgeWidth / 2, topEdgeLenght / 2, fHeight / 2);
		std::unique_ptr<G4IntersectionSolid> foodSolid = std::make_unique<G4IntersectionSolid>("cellsSolid", mainBodyInternalSolidWithNeck, foodBaseSolid.get(), nullptr, foodSolidPositionOffset);
		return foodVolume - foodSolid->GetCubicVolume();
	};

	// approximate food height to get expected food volume
	foodHeight = boost::math::tools::bisect(fn, 2 * mm, 4 * mm, boost::math::tools::eps_tolerance<G4double>(10)).second;

	G4ThreeVector foodSolidPositionOffset = cellsSolidPositionOffset + G4ThreeVector(0, 0, cellsLayerThickness / 2 + foodHeight / 2);

	auto foodBaseSolid = new G4Box("foodBaseSolid", topEdgeWidth / 2, topEdgeLenght / 2, foodHeight / 2);

	auto foodSolid = new G4IntersectionSolid("cellsSolid", mainBodyInternalSolidWithNeck, foodBaseSolid, nullptr, foodSolidPositionOffset);

	//std::cout << "Cells surface area: " << cellsSolid->GetCubicVolume() / cellsLayerThickness / cm2 << "cm^2\n";
	std::cout << std::format("Food height: {} mm\nFood volume: {} cm^3\n", foodHeight / mm, foodSolid->GetCubicVolume() / cm3);

	// logical volume for cellsSolid
	G4LogicalVolume* cellsLogical = new G4LogicalVolume(cellsSolid, Water, "cellsLogical");

	G4LogicalVolume* foodLogical = new G4LogicalVolume(foodSolid, Water, "foodLogical");

	G4VPhysicalVolume* cellsPhysical = new G4PVPlacement(rotation, G4ThreeVector(0, 0, 0) * mm, "cellsPhysical", cellsLogical, mainBodyInternalPhysical, false, 0);

	G4VPhysicalVolume* foodPhysical = new G4PVPlacement(rotation, G4ThreeVector(0, 0, 0) * mm, "foodPhysical", foodLogical, mainBodyInternalPhysical, false, 0);

	//<--------------------------------------------------------------------------------------------------------------------------------->
	//<--------------------------------------------------Applicator and phantom geometry------------------------------------------------>
	//<--------------------------------------------------------------------------------------------------------------------------------->

	const G4double
		topPhantomThickness = 0.5 * cm,
		bottomPhantomThickness = 11.5 * cm,
		applicatorLenght = 60 * cm,
		applicatorDiameter = 10 * cm,
		phantomSide = 12 * cm;

	G4Tubs* applicatorSolid = new G4Tubs("applicatorSolid", 0.0, applicatorDiameter / 2, applicatorLenght / 2, 0.0, CLHEP::twopi);

	G4LogicalVolume* applicatorLogicalVolume = new G4LogicalVolume(applicatorSolid, PMMA, "applicatorLogicalVolume");

	G4Box* topPhantomSolid = new G4Box("BottomPhantomSolid", phantomSide / 2, phantomSide / 2, topPhantomThickness / 2);

	G4LogicalVolume* topPhantomLogical = new G4LogicalVolume(topPhantomSolid, RW3PhantomMaterial, "topPhantomLogical");

	G4Box* bottomPhantomSolid = new G4Box("BottomPhantomSolid", phantomSide / 2, phantomSide / 2, bottomPhantomThickness / 2);

	G4LogicalVolume* bottomPhantomLogical = new G4LogicalVolume(bottomPhantomSolid, RW3PhantomMaterial, "bottomPhantomLogical");

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
