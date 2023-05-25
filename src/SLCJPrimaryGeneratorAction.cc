/////////////////////////////////////////////////////////////////////////
//
// V. Guadilla 2021
/////////////////////////////////////////////////////////////////////////

#include "SLCJPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GeneralParticleSource.hh" 
#include "G4IonTable.hh"
#include "G4RandomDirection.hh"

#include "G4ThreeVector.hh"
#include "globals.hh"
#include "Randomize.hh"

#include "G4ios.hh"
#include "fstream"
#include "iomanip"
#include <cstdlib>

#include "G4SystemOfUnits.hh"

#include "SLCJRunAction.hh"

using namespace std;

inline std::string filename_string(std::string path_str) {
	return path_str.substr(path_str.rfind("\\") + 1, path_str.size() - path_str.rfind("\\") - 1);
};

#define _endl_ " (" << filename_string(__FILE__) << "; " << __LINE__ << ")" << '\n'
#define checkpoint std::cout << "checkpoint" << _endl_

SLCJPrimaryGeneratorAction::SLCJPrimaryGeneratorAction(SLCJRunAction* RunAct)
	:runAction(RunAct)
{

	//choose the Random engine 
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
	//set random seed with system time (it must change from one run to another!!! otherwise we simulate all the time the same in the SLCJ case!!)
	G4long seed = time(NULL);
	CLHEP::HepRandom::setTheSeed(seed);
	//////////////////////////////////////////////////////////////////////////////

	//////////Reading the input data for primary generator///////////

	ifstream evenInputInformation;
	evenInputInformation.open("../../../particles3.data");
	if (!evenInputInformation.is_open()) {
		std::cout << "\n\nNO EVENT INPUT INFORMATION FILE FOUND!!!" << _endl_;
		exit(1);
	}
	std::string header1, header2;

	// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*
	//Geantino for geometry checking
	G4int n_particle = 1;
	particleGun = new G4ParticleGun(n_particle);

	// particle type
	G4ParticleTable particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;

	// particle type
	particleGun->SetParticleDefinition(particleTable->FindParticle(particleName = "geantino"));
	particleGun->SetParticleEnergy(1.0  GeV);
	particleGun->SetParticleMomentumDirection(G4RandomDirection());
	*/

	/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Simulation of physical particles

	// particles per event
	G4int n_particle = 1;
	particleGun = std::make_unique<G4ParticleGun>(n_particle);

	// particle type
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;

	/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	electron = particleTable->FindParticle(particleName = "electron");

}



void SLCJPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

	/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*
	//Geantino for geometry checking

	// particle definitions
	G4ParticleTable particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;

	particleGun->SetParticleDefinition(particleTable->FindParticle(particleName = "geantino"));
	particleGun->SetParticleEnergy(1.0  GeV);
	particleGun->SetParticleMomentumDirection(G4RandomDirection());

	particleGun->GeneratePrimaryVertex(anEvent);
	*/
	/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Simulation of physical particles

	x = 0.0 * mm;
	y = 0.0 * mm;
	z = 0.0 * mm;

	std::array tpl = { E / keV, x / mm, y / mm, z / mm };
	metaFile.write((char*)tpl.data(), sizeof(tpl));


	G4ThreeVector momentumDirection(0,0,-1);
	particleGun->SetParticleDefinition(electron);

	particleGun->SetParticlePosition(G4ThreeVector(x, y, z));
	particleGun->SetParticleMomentumDirection(momentumDirection);
	particleGun->SetParticleEnergy(E);
	particleGun->GeneratePrimaryVertex(anEvent);
}

void SLCJPrimaryGeneratorAction::setRunPath(std::filesystem::path runPath) {
	metaFile.open(runPath / "metadata.bin", std::ios_base::binary | std::ios_base::out);
}

G4double& SLCJPrimaryGeneratorAction::getEnergy() {
	return E;
}


