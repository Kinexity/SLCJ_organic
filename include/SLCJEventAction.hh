// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02EventAction.hh,v 1.3 1999/12/15 14:49:20 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
///////////////////////////////////////////////////////////////////
//
// V. Guadilla 2021
// Accumulated the energy deposited in sensitive volumes
//
//////////////////////////////////////////////////////////////////

#ifndef SLCJEventAction_h
#define SLCJEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <array>
#include <vector>

class G4Event;
class SLCJRunAction;

class SLCJEventAction : public G4UserEventAction
{
public:
	SLCJEventAction(SLCJRunAction*);
	~SLCJEventAction();

public:
	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);
	void addEdep(G4double Edep, G4double x, G4double y, G4double z);
	void add_E_i(G4int nCrystal, G4double edep);

private:
	SLCJRunAction* runAction;
	G4int Range;

	std::vector<std::array<G4double, 4>>
		EnergyDeposit;
	std::array<G4double, 20>
		TotalEnergyDepositCrystal = { 0 };
};

#endif
