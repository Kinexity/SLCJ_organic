// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN01DetectorConstruction.hh,v 1.2 1999/12/15 14:49:18 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
//////////////////////////////////////////////////////////////////////////
//
// V. Guadilla 2021
// Simple SLCJ detector geometry

#ifndef SLCJDetectorConstruction_h
#define SLCJDetectorConstruction_h 1

// STL //
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <array>

#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"

class G4VPhysicalVolume;
class F02ElectricFieldSetup;

class SLCJDetectorConstruction : public G4VUserDetectorConstruction
{
public:
	SLCJDetectorConstruction() = default;
	SLCJDetectorConstruction(G4double foodV);
	~SLCJDetectorConstruction() = default;
	G4VPhysicalVolume* Construct();
	const G4double getFoodVolume();
	// this function saves details of the simulation into a file - physical volumes, materials and elements
	void saveDetails(std::filesystem::path p);
private:
	F02ElectricFieldSetup* fEmFieldSetup;
	G4String              header1, header2, header3;
	G4double T = 0., P = 0., E = 0., d = 0.;
	
	const G4double foodVolume = 10 * cm3;
};

#endif

