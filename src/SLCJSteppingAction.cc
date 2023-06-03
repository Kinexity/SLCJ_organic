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
	//std::cout << currentVolumeName << '\t' << processName << '\n';

	G4String nameP = aStep->GetTrack()->GetDefinition()->GetParticleName();

	//if (edep>0.0 & currentVolumeName == "GasSLCJ" & nameP=="alpha"){
	//if (edep > 0.0 & currentVolumeName == "GasSLCJ") {
	//
	//	G4StepPoint* prePoint = aStep->GetPreStepPoint();
	//	G4StepPoint* postPoint = aStep->GetPostStepPoint();
	//
	//	G4double x1 = prePoint->GetPosition().x();
	//	G4double x2 = postPoint->GetPosition().x();
	//	G4double y1 = prePoint->GetPosition().y();
	//	G4double y2 = postPoint->GetPosition().y();
	//	G4double z1 = prePoint->GetPosition().z();
	//	G4double z2 = postPoint->GetPosition().z();
	//
	//	G4double x = x1 + G4UniformRand() * (x2 - x1);
	//	G4double y = y1 + G4UniformRand() * (y2 - y1);
	//	G4double z = z1 + G4UniformRand() * (z2 - z1);
	//
	//	eventAction->addEdep(edep / keV, x / mm, y / mm, z / mm);
	//	//G4cout<<edep/keV<<"    "<<x/mm<<"    "<<y/mm<<"    "<<z/mm<<G4endl;
	//}

	eventAction->addGeantinoPosition(currentVolumeName, postPoint->GetPosition());

	if (edep > 0.0 && currentVolumeName == "organicMaterialLogical") {

		G4int nCrystal = touch->GetCopyNumber(1); //N will be the number of levels up, we have to check it to pickup the index of CeBr3 crystal
		//eventAction->add_E_i(nCrystal, edep / keV);
	}

}


