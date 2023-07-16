// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02SteppingAction.cc,v 1.3 1999/12/15 14:49:22 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
/////////////////////////////////////////////////////////////////////
//
//  V.Guadilla 2019
//
//  Collect the energy deposited in the sensitive SLCJ volume
///////////////////////////////////////////////////////////////////

#include "SLCJSteppingAction.hh"
#include "SLCJEventAction.hh"
#include "G4SteppingManager.hh"

#include "G4SystemOfUnits.hh"

SLCJSteppingAction::SLCJSteppingAction(SLCJEventAction* EvAct) :eventAction(EvAct) {}

inline std::string filename_string(std::string path_str) {
	return path_str.substr(path_str.rfind("\\") + 1, path_str.size() - path_str.rfind("\\") - 1);
};

#define _endl_ " (" << filename_string(__FILE__) << "; " << __LINE__ << ")" << '\n'
#define checkpoint std::cout << "checkpoint" << _endl_


void SLCJSteppingAction::UserSteppingAction(const G4Step* aStep)
{

	G4double edep = aStep->GetTotalEnergyDeposit();

	G4StepPoint* prePoint = aStep->GetPreStepPoint();
	G4StepPoint* postPoint = aStep->GetPostStepPoint();
	G4TouchableHandle touch = prePoint->GetTouchableHandle();
	const std::string currentVolumeName = touch->GetVolume()->GetName();
	const std::string currentMaterialName = touch->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName();
	const G4VProcess* process = postPoint->GetProcessDefinedStep();
	std::string processName = process->GetProcessName();

	G4double EdepStep;
	G4int begin;

	G4String nameP = aStep->GetTrack()->GetDefinition()->GetParticleName();

	if (nameP == "geantino") {
		eventAction->addGeantinoPosition(currentVolumeName, postPoint->GetPosition());
	}
	else {
		if (edep > 0.0 && currentVolumeName == "cellsPhysical") {
			auto prePos = prePoint->GetPosition();
			auto postPos = postPoint->GetPosition();
			auto deltaPos = postPos - prePos;
			auto pos = prePos + G4UniformRand() * deltaPos; // position randomization simplified
			eventAction->addEdep(edep / keV, pos.getX() / mm, pos.getY() / mm, pos.getZ() / mm);
		}
	}

}


