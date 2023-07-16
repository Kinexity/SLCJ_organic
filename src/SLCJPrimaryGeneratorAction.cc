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

G4ThreeVector GetRandomPositionInCircle(double radius) {
	double phi = G4UniformRand() * CLHEP::twopi;  // Random angle between 0 and 2*pi
	double r = sqrt(G4UniformRand()) * radius;  // Random radius between 0 and radius
	G4ThreeVector position;
	position.setRhoPhiZ(r, phi, 0);
	return position;
}

void SLCJPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

	if (runAction->getIsGeantino()) {
		particleGun->SetParticleDefinition(geantino);
		particleGun->SetParticlePosition(G4ThreeVector());
		particleGun->SetParticleEnergy(1.0 * GeV);
		particleGun->SetParticleMomentumDirection(G4RandomDirection());
	}
	else {
		// offset to put center of electron beam over the center of cell sample on top of PMMA tube
		G4ThreeVector initialPositionOffset(0, 0, 74.6 / 2 * cm);

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
	}
	particleGun->GeneratePrimaryVertex(anEvent);
	//*/
}

void SLCJPrimaryGeneratorAction::setRunPath(std::filesystem::path runPath) {
	metaFile.open(runPath / "metadata.bin", std::ios_base::binary | std::ios_base::out);
}

G4double& SLCJPrimaryGeneratorAction::getEnergy() {
	return E;
}


