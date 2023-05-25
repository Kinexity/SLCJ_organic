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


#include "G4SystemOfUnits.hh"
// #include "g4std/fstream"
// #include "g4std/iomanip"

using namespace std;


SLCJRunAction::SLCJRunAction()
{
	timer = std::make_unique<G4Timer>();


	///////////////////////////////////////////////////////////////////////////////////	

}

void SLCJRunAction::BeginOfRunAction(const G4Run*)
{
	//eventTotalDepositFile.open(eventTotalDepositFilePath.string() + ".csv", std::ios_base::out | std::ios_base::trunc);
	//eventStepsDepositFile.open(eventStepsDepositFilePath.string() + ".csv", std::ios_base::out | std::ios_base::trunc);
	eventTotalDepositFileBinary.open(eventTotalDepositFilePath.string() + ".bin", std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
	//eventStepsDepositFileBinary.open(eventStepsDepositFilePath.string() + ".bin", std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
	//Start CPU timer
	timer->Start();
	eventIndex = 0;
}

void SLCJRunAction::EndOfRunAction(const G4Run*)
{

	//eventTotalDepositFile.close();
	//eventStepsDepositFile.close();
	eventTotalDepositFileBinary.close();
	//eventStepsDepositFileBinary.close();
	//Stop timer and get CPU time
	timer->Stop();
	G4double cputime = timer->GetRealElapsed();
	cout << " CPU time = " << cputime << " s" << endl;


}

void SLCJRunAction::fillOut(std::vector<std::array<G4double, 4>>& EnergyDeposit, std::array<G4double, 20>& EnergyGammaCrystals)
{
	//if (EnergyDeposit.size() > 0) {
	//	eventStepsDepositFile << eventIndex << '\t' << EnergyDeposit.size() << '\n';
	//	for (auto& EnergyDeposit_i : EnergyDeposit) {
	//		eventStepsDepositFile << EnergyDeposit_i[0] << ',' << EnergyDeposit_i[1] << ',' << EnergyDeposit_i[2] << ',' << EnergyDeposit_i[3] << '\n';
	//	}
	//}

	eventTotalDepositFileBinary.write((char*)EnergyGammaCrystals.data(), EnergyGammaCrystals.size() * sizeof(G4double));

	//for (auto& EnergyDepositOneCrystal : EnergyGammaCrystals) {
	//	eventTotalDepositFile << EnergyDepositOneCrystal << '\t';
	//}
	//eventTotalDepositFile << "\n";

	eventIndex++;
	if (eventIndex % 10000 == 0) {
		std::cout << eventIndex << '\n';
	}
}

void SLCJRunAction::setEventFilePath(std::filesystem::path totalP, std::filesystem::path stepsP) {
	eventTotalDepositFilePath = totalP;
	eventStepsDepositFilePath = stepsP;
	std::cout << eventTotalDepositFilePath << '\n' << eventStepsDepositFilePath << '\n';
}
