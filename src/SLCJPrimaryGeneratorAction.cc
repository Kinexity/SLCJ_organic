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
#include "G4Electron.hh"

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

	electron = G4Electron::Definition();
	geantino = particleTable->FindParticle(particleName = "geantino");
}

// random position with symmetry and uniformity - check before use, probably wrong
G4ThreeVector GetRandomPositionInCircle2(double radius) {
	double phi = G4UniformRand() * CLHEP::twopi;  // Random angle between 0 and 2*pi
	double r = sqrt(G4UniformRand()) * radius;  // Random radius between 0 and radius

	// Apply <3% symmetry to the angle
	phi += G4RandGauss::shoot(0, 0.03 * CLHEP::twopi);
	if (phi > CLHEP::twopi)
		phi -= CLHEP::twopi;
	else if (phi < 0)
		phi += CLHEP::twopi;

	// Apply <5% uniformity to the radius
	double deltaR = radius * 0.05;
	r += G4RandFlat::shoot(-deltaR, deltaR);
	if (r > radius)
		r = radius;
	else if (r < 0)
		r = 0;

	G4ThreeVector position;
	position.setRhoPhiZ(r, phi, 0);
	return position;
}

G4ThreeVector GetRandomPositionInCircle(double radius) {
	double phi = G4UniformRand() * CLHEP::twopi;  // Random angle between 0 and 2*pi
	double r = sqrt(G4UniformRand()) * radius;  // Random radius between 0 and radius
	G4ThreeVector position;
	position.setRhoPhiZ(r, phi, 0);
	return position;
}

void SLCJPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

	/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/*
	//Geantino for geometry checking

	// particle definitions
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;

	particleGun->SetParticleDefinition(geantino);
	particleGun->SetParticlePosition(G4ThreeVector(0, 4, 1) * cm);
	particleGun->SetParticleEnergy(1.0 * GeV);
	particleGun->SetParticleMomentumDirection(G4RandomDirection());
	//particleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1));

	particleGun->GeneratePrimaryVertex(anEvent);
	*/
	/// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Simulation of physical particles
	///*

	// offset to put center of electron beam over the center of cell sample
	G4ThreeVector initialPositionOffset(0, 3.5 * cm, 4 * cm);

	auto particlePosition = GetRandomPositionInCircle(fieldDiameter / 2) + initialPositionOffset;

	std::array tpl = { E / keV, x / mm, y / mm, z / mm };
	metaFile.write((char*)tpl.data(), sizeof(tpl));

	// -Z electron direction
	G4ThreeVector momentumDirection(0, 0, -1);
	particleGun->SetParticleDefinition(electron);

	particleGun->SetParticlePosition(G4ThreeVector(x, y, z));
	particleGun->SetParticleMomentumDirection(momentumDirection);
	particleGun->SetParticlePosition(particlePosition);
	particleGun->SetParticleEnergy(E);
	particleGun->GeneratePrimaryVertex(anEvent);
	//*/
}

void SLCJPrimaryGeneratorAction::setRunPath(std::filesystem::path runPath) {
	metaFile.open(runPath / "metadata.bin", std::ios_base::binary | std::ios_base::out);
}

G4double& SLCJPrimaryGeneratorAction::getEnergy() {
	return E;
}


