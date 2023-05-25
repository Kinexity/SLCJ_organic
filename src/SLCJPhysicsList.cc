//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file electromagnetic/TestEm11/src/SLCJPhysicsList.cc
/// \brief Implementation of the SLCJPhysicsList class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SLCJPhysicsList.hh"
#include "SLCJPhysicsListMessenger.hh"

#include "PhysListEmStandard.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"

#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"


G4ThreadLocal StepMax* SLCJPhysicsList::fStepMaxProcess = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SLCJPhysicsList::SLCJPhysicsList(std::string PLname) : G4VModularPhysicsList(), physicsListName(PLname)
{
	fMessenger = std::make_unique<SLCJPhysicsListMessenger>(this);

	SetVerboseLevel(1);

	// EM physics
	if (physicsListName == "local") {
		fEmSLCJPhysicsList = std::make_unique<PhysListEmStandard>("local");
	}
	else if (physicsListName == "emlivermore") {
		fEmSLCJPhysicsList = std::make_unique<G4EmLivermorePhysics>();
	}
	else if (physicsListName == "emstandard_opt3") {
		fEmSLCJPhysicsList = std::make_unique<G4EmStandardPhysics_option3>();
	}
	else if (physicsListName == "empenelope") {
		fEmSLCJPhysicsList = std::make_unique<G4EmPenelopePhysics>();
	}
	G4LossTableManager::Instance();
	SetDefaultCutValue(1. * mm);
	//SetDefaultCutValue(0.1*mm);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SLCJPhysicsList::ConstructParticle()
{
	G4BosonConstructor  pBosonConstructor;
	pBosonConstructor.ConstructParticle();

	G4LeptonConstructor pLeptonConstructor;
	pLeptonConstructor.ConstructParticle();

	G4MesonConstructor pMesonConstructor;
	pMesonConstructor.ConstructParticle();

	G4BaryonConstructor pBaryonConstructor;
	pBaryonConstructor.ConstructParticle();

	G4IonConstructor pIonConstructor;
	pIonConstructor.ConstructParticle();

	G4ShortLivedConstructor pShortLivedConstructor;
	pShortLivedConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ProcessManager.hh"

void SLCJPhysicsList::ConstructProcess()
{
	// transportation
	//
	AddTransportation();

	// electromagnetic physics list
	//
	fEmSLCJPhysicsList->ConstructProcess();

	// Em options
	//
	G4EmParameters* param = G4EmParameters::Instance();
	param->SetBuildCSDARange(true);

	// decay process
	//
	AddDecay();

	// radioactive decay Process
	//
	AddRadioactiveDecay();

	// step limitation (as a full process)
	//  
	AddStepMax();

	// Ion Gas models
	//AddIonGasModels();

	// Add Physics List
	//
	AddSLCJPhysicsList();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SLCJPhysicsList::AddSLCJPhysicsList()
{

	if (fEmSLCJPhysicsList != nullptr) return;

	if (physicsListName == "local") {
		fEmSLCJPhysicsList = std::make_unique<PhysListEmStandard>(physicsListName);
	}
	else if (physicsListName == "emstandard_opt0") {
		fEmSLCJPhysicsList = std::make_unique<G4EmStandardPhysics>();
	}
	else if (physicsListName == "emstandard_opt1") {
		fEmSLCJPhysicsList = std::make_unique<G4EmStandardPhysics_option1>();
	}
	else if (physicsListName == "emstandard_opt2") {
		fEmSLCJPhysicsList = std::make_unique<G4EmStandardPhysics_option2>();
	}
	else if (physicsListName == "emstandard_opt3") {
		fEmSLCJPhysicsList = std::make_unique<G4EmStandardPhysics_option3>();
	}
	else if (physicsListName == "emstandard_opt4") {
		fEmSLCJPhysicsList = std::make_unique<G4EmStandardPhysics_option4>();
	}
	else if (physicsListName == "emstandardSS") {
		fEmSLCJPhysicsList = std::make_unique<G4EmStandardPhysicsSS>();
	}
	else if (physicsListName == "emstandardGS") {
		fEmSLCJPhysicsList = std::make_unique<G4EmStandardPhysicsGS>();
	}
	else if (physicsListName == "emstandardWVI") {
		fEmSLCJPhysicsList = std::make_unique<G4EmStandardPhysicsWVI>();
	}
	else if (physicsListName == "emlowenergy") {
		fEmSLCJPhysicsList = std::make_unique<G4EmLowEPPhysics>();
	}
	else if (physicsListName == "emlivermore") {
		fEmSLCJPhysicsList = std::make_unique<G4EmLivermorePhysics>();
	}
	else if (physicsListName == "empenelope") {
		fEmSLCJPhysicsList = std::make_unique<G4EmPenelopePhysics>();
	}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"

void SLCJPhysicsList::AddDecay()
{
	// decay process
	//
	G4Decay* fDecayProcess = new G4Decay();

	auto particleIterator = GetParticleIterator();
	particleIterator->reset();
	while ((*particleIterator)()) {
		G4ParticleDefinition* particle = particleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();

		if (fDecayProcess->IsApplicable(*particle) && !particle->IsShortLived()) {

			pmanager->AddProcess(fDecayProcess);

			// set ordering for PostStepDoIt and AtRestDoIt
			pmanager->SetProcessOrdering(fDecayProcess, idxPostStep);
			pmanager->SetProcessOrdering(fDecayProcess, idxAtRest);

		}
	}
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4PhysicsListHelper.hh"
#include "G4RadioactiveDecay.hh"
#include "G4GenericIon.hh"
#include "G4NuclideTable.hh"

void SLCJPhysicsList::AddRadioactiveDecay()
{
	G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();

	radioactiveDecay->SetARM(true);                //Atomic Rearangement

	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
	ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());

	// mandatory for G4NuclideTable
	//
	G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1 * picosecond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StepMax.hh"

void SLCJPhysicsList::AddStepMax()
{
	// Step limitation seen as a process
	fStepMaxProcess = new StepMax();

	G4double stepMax = 0.1;

	//G4double stepMax=0.01; //as in 2008 IEEE Nuclear Science Symposium Conference Record 
	fStepMaxProcess->SetMaxStep2(stepMax);

	auto particleIterator = GetParticleIterator();
	particleIterator->reset();
	while ((*particleIterator)()) {
		G4ParticleDefinition* particle = particleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();

		if (fStepMaxProcess->IsApplicable(*particle))
		{
			pmanager->AddDiscreteProcess(fStepMaxProcess);
		}
	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//////////////////////////////////////////////////////////////////
/// Adds the ion gas model


#include "G4IonFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4BraggIonGasModel.hh"
#include "G4BetheBlochIonGasModel.hh"

#include "G4PhysicalConstants.hh"

#include "G4EmConfigurator.hh"

void SLCJPhysicsList::AddIonGasModels() {
	G4EmConfigurator* em_config =
		G4LossTableManager::Instance()->EmConfigurator();
	auto particleIterator = GetParticleIterator();
	particleIterator->reset();
	while ((*particleIterator)())
	{
		G4ParticleDefinition* particle = particleIterator->value();
		G4String partname = particle->GetParticleName();
		if (partname == "alpha" || partname == "He3" || partname == "GenericIon" || partname == "proton") {
			G4BraggIonGasModel* mod1 = new G4BraggIonGasModel();
			G4BetheBlochIonGasModel* mod2 = new G4BetheBlochIonGasModel();
			G4double eth = 2. * MeV * particle->GetPDGMass() / proton_mass_c2;
			em_config->SetExtraEmModel(partname, "ionIoni", mod1, "", 0.0, eth,
				new G4IonFluctuations());
			em_config->SetExtraEmModel(partname, "ionIoni", mod2, "", eth, 100 * TeV,
				new G4UniversalFluctuation());

		}
	}
}

std::string SLCJPhysicsList::getPhysicsListName() {
	return fEmSLCJPhysicsList->GetPhysicsName();
}
