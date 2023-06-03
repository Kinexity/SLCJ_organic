/////////////////////////////////////////////////////////////////////////
//
//     04-2002  A. Algora--> SLCJDetectorConstruction.cc
//
//     Test of Geometries                                              //
/////////////////////////////////////////////////////////////////////////

#ifndef SLCJPrimaryGeneratorAction_h
#define SLCJPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include <array>
#include <filesystem>
#include <fstream>

class G4Event;
class SLCJRunAction;

extern std::ifstream eventInputFile;

class SLCJPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	SLCJPrimaryGeneratorAction(SLCJRunAction*);
	~SLCJPrimaryGeneratorAction() = default;

public:
	void GeneratePrimaries(G4Event* anEvent);
	void setRunPath(std::filesystem::path runPath);
	G4ParticleGun* GetParticleGun() { return particleGun.get(); };
	G4double& getEnergy();
private:
	SLCJRunAction* runAction;
	std::unique_ptr<G4ParticleGun> particleGun;

	std::fstream metaFile;
	G4double E, x, y, z;
	G4double fieldDiameter = 10 * cm;

	G4ParticleDefinition
		*electron,
		*geantino;

	std::array<G4ParticleDefinition*, 4> particleDefinitions;
};

#endif
