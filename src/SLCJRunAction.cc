// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02RunAction.cc,v 1.3 1999/12/15 14:49:22 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
////////////////////////////////////////////////////////////////
//
//  V. Guadilla 2021
///////////////////////////////////////////////////////////////

#include "SLCJRunAction.hh"

#include "G4Run.hh"

#include "G4ios.hh"
#include "fstream"
#include "iomanip"
#include "Randomize.hh"
#include "time.h"
#include <set>
#include <map>


#include "G4SystemOfUnits.hh"
// #include "g4std/fstream"
// #include "g4std/iomanip"

using namespace std;


SLCJRunAction::SLCJRunAction() {
	timer = std::make_unique<G4Timer>();
}

void SLCJRunAction::BeginOfRunAction(const G4Run*) {

	eventEnergyDepositFile.open(eventEnergyDepositFilePath, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);

	//Start CPU timer
	timer->Start();
	eventIndex = 0;
}

void SLCJRunAction::EndOfRunAction(const G4Run*)
{
	eventEnergyDepositFile.close();
	//Stop timer and get CPU time
	timer->Stop();
	G4double cputime = timer->GetRealElapsed();
	cout << " CPU time = " << cputime << " s" << endl;


}

void SLCJRunAction::fillOut(std::vector<std::array<G4double, 4>>& energyDeps) {
	eventEnergyDepositFile.write((char*)energyDeps.data(), energyDeps.size() * sizeof(std::array<G4double, 4>));
}

void SLCJRunAction::setEventFilePath(std::filesystem::path energyP) {
	eventEnergyDepositFilePath = energyP;
	std::cout << eventEnergyDepositFilePath << '\n';
}
