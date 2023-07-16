// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02RunAction.hh,v 1.3 1999/12/15 14:49:20 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
//
/////////////////////////////////////////////////////////////////
//
//  V Guadilla 2021
//
/////////////////////////////////////////////////////////////////

#ifndef SLCJRunAction_h
#define SLCJRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <filesystem>
#include "G4Timer.hh"
#include <G4ThreeVector.hh>

class G4Run;

class G4Timer;

class SLCJRunAction : public G4UserRunAction
{
public:
	SLCJRunAction(bool iG);
	~SLCJRunAction() = default;

public:
	void BeginOfRunAction(const G4Run*);
	void EndOfRunAction(const G4Run*);

	void fillOut(std::vector<std::array<G4double, 4>>& energyDeps);
	void fillOut(std::vector<std::pair<G4String, G4ThreeVector>>& geantinoPos);

	void setEventFilePath(std::filesystem::path energyP);

	bool getIsGeantino();
private:

	std::unique_ptr<G4Timer> timer;

	std::fstream
		eventEnergyDepositFile;

	uint32_t eventIndex;

	std::filesystem::path
		eventEnergyDepositFilePath;

	std::fstream file;

	std::vector<std::pair<G4String, G4ThreeVector>> geantinoPosGlobal;

	const bool isGeantino;
};

#endif






