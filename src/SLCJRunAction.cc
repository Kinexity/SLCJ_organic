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


SLCJRunAction::SLCJRunAction(bool iG) : isGeantino(iG) {
	timer = std::make_unique<G4Timer>();
}

void SLCJRunAction::BeginOfRunAction(const G4Run*) {

	if (isGeantino) {
		file.open("geantino.txt", std::ios_base::out | std::ios_base::trunc);
	}
	else {
		eventEnergyDepositFile.open(eventEnergyDepositFilePath, std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
	}

	//Start CPU timer
	timer->Start();
	eventIndex = 0;
}

void SLCJRunAction::EndOfRunAction(const G4Run*)
{
	if (isGeantino) {
		std::set<G4String> vols;
		for (auto [vol, pos] : geantinoPosGlobal) {
			vols.insert(vol);
		}
		std::map<G4String, std::array<int, 3>> vol_colours;
		int i = 0;
		//std::cout << vols.size() << '\n';
		for (auto vol : vols) {
			//std::cout << vol << '\n';
			vol_colours[vol] = { 255 * (i % 2), 255 * (i / 2 % 2), 255 * (i / 4 % 2) };
			std::cout << vol << '\t' << G4ThreeVector(255 * (i % 2), 255 * (i / 2 % 2), 255 * (i / 4 % 2)) << '\n';
			i++;
		}
		for (auto [vol, pos] : geantinoPosGlobal) {
			auto colour = vol_colours[vol];
			file << std::format("{}\t{}\t{}\t{}\t{}\t{}\n", pos.getX(), pos.getY(), pos.getZ(), colour[0], colour[1], colour[2]);
		}
		file.close();
	}
	else {
		eventEnergyDepositFile.close();
	}
	//Stop timer and get CPU time
	timer->Stop();
	G4double cputime = timer->GetRealElapsed();
	cout << " CPU time = " << cputime << " s" << endl;


}

void SLCJRunAction::fillOut(std::vector<std::array<G4double, 4>>& energyDeps) {
	eventEnergyDepositFile.write((char*)energyDeps.data(), energyDeps.size() * sizeof(std::array<G4double, 4>));
}

void SLCJRunAction::fillOut(std::vector<std::pair<G4String, G4ThreeVector>>& geantinoPos) {
	for (auto temp : geantinoPos) {
		geantinoPosGlobal.push_back(temp);
	}
}

void SLCJRunAction::setEventFilePath(std::filesystem::path energyP) {
	eventEnergyDepositFilePath = energyP;
	std::cout << eventEnergyDepositFilePath << '\n';
}

bool SLCJRunAction::getIsGeantino() {
	return isGeantino;
}
