// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: exampleN01.cc,v 1.2 1999/12/15 14:49:18 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
//
// --------------------------------------------------------------
//      GEANT 4 - exampleN01
//
//      For information related to this code contact:
//      CERN, IT Division, ASD Group
// --------------------------------------------------------------
//
/////////////////////////////////////////////////////////////////
//
//  Modified: 16-01-2013 A. Algora to run with the newer versions
//
//  To run simple SLCJ simulations
//
////////////////////////////////////////////////////////////////

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef G4UI_USE_TCSH
#include "G4UItcsh.hh"
#endif

#ifdef G4VIS_USE
// #include "SLCJVisManager.hh"
#include "G4VisExecutive.hh"
#endif

#include "SLCJDetectorConstruction.hh"
#include "SLCJPhysicsList.hh"
#include "SLCJPrimaryGeneratorAction.hh"
#include "SLCJRunAction.hh"
#include "SLCJEventAction.hh"
#include "SLCJSteppingAction.hh"

#include "Randomize.hh"
#include "globals.hh"


#include "G4ios.hh"
#include <fstream>
#include <iomanip>
#include <filesystem>
#include <memory>
#include <chrono>
#include <format>
#include <algorithm>
#include <numeric>

#include <boost/program_options.hpp>

inline std::string filename_string(std::string path_str) {
	return path_str.substr(path_str.rfind("\\") + 1, path_str.size() - path_str.rfind("\\") - 1);
};

#define _endl_ " (" << filename_string(__FILE__) << "; " << __LINE__ << ")" << '\n'
#define checkpoint std::cout << "checkpoint" << _endl_

namespace po = boost::program_options;


int main(int argc, char** argv) {

	bool
		skipIfDataExists = false,
		dataOverwrite = false,
		isGeantino = true;
	G4double
		foodVolume = 10 * cm3,
		cutValue = 0.01 * mm,
		energy = 9 * MeV;
	std::string
		physicsListName = "emlivermore";
	G4int
		numberOfEvent = 1000000;

	auto start = std::chrono::high_resolution_clock::now();

	// Define options and their descriptions
	po::options_description desc("Allowed options");
	desc.add_options()
		("skipIfDataExists", po::bool_switch(&skipIfDataExists), "Skip if data exists")
		("dataOverwrite", po::bool_switch(&dataOverwrite), "Data overwrite")
		("isGeantino", po::bool_switch(&isGeantino), "Is geantino")
		("foodVolume", po::value<double>(&foodVolume), "Food volume")
		("cutValue", po::value<double>(&cutValue), "Cut value")
		("energy", po::value<double>(&energy), "Energy")
		("physicsListName", po::value<std::string>(&physicsListName), "Physics list name")
		("numberOfEvent", po::value<int>(&numberOfEvent), "Number of events");

	// Parse command line arguments
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	// apply units if necessery
	if (vm.count("energy")) {
		energy *= MeV;
	}
	if (vm.count("foodVolume")) {
		foodVolume *= cm3;
	}
	if (vm.count("cutValue")) {
		cutValue *= mm;
	}

	// Print the values
	std::cout << "skipIfDataExists: " << std::boolalpha << skipIfDataExists << std::endl;
	std::cout << "dataOverwrite: " << std::boolalpha << dataOverwrite << std::endl;
	std::cout << "isGeantino: " << std::boolalpha << isGeantino << std::endl;
	std::cout << "foodVolume: " << foodVolume / cm3 << std::endl;
	std::cout << "cutValue: " << cutValue / mm << std::endl;
	std::cout << "energy: " << energy / MeV << std::endl;
	std::cout << "physicsListName: " << physicsListName << std::endl;
	std::cout << "numberOfEvent: " << numberOfEvent << std::endl;


	std::filesystem::path directory = std::filesystem::current_path();
	if (!std::filesystem::exists(directory)) {
		std::cout << "Set correct home directory!" << _endl_;
		exit(1);
	}
	auto resultsDirectoryPath = directory / "results_SLCJ";
	if (!std::filesystem::exists(resultsDirectoryPath)) {
		std::filesystem::create_directory(resultsDirectoryPath);
	}


	// Choose the random engine and initialize
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanluxEngine);

	G4long Seed = 2193585;
	G4int Lux = 3;
	CLHEP::HepRandom::setTheSeed(Seed, Lux);

	checkpoint;
#ifdef  WIN32
	if (false) {
		_putenv_s("G4ENSDFSTATEDATA", "C:\\Program Files (x86)\\Geant4 10.5\\data\\G4ENSDFSTATE2.2");
		system("set G4ENSDFSTATEDATA");
		_putenv_s("G4LEVELGAMMADATA", "C:\\Program Files(x86)\\Geant4 10.5\\data\\PhotonEvaporation5.3");
		system("set G4LEVELGAMMADATA");
		_putenv_s("G4RADIOACTIVEDATA", "C:\\Program Files (x86)\\Geant4 10.5\\data\\RadioactiveDecay5.3");
		system("set G4RADIOACTIVEDATA");
		_putenv_s("G4LEDATA", "C:\\Program Files (x86)\\Geant4 10.5\\data\\G4EMLOW7.7");
		system("set G4LEDATA");
	}
#endif //  WIN32
	// construct the default run manager
	std::unique_ptr<G4RunManager> runManager = std::make_unique<G4RunManager>();

	checkpoint;
	// set mandatory initialization classes
	auto SLCJdetector = new SLCJDetectorConstruction(foodVolume);
	runManager->SetUserInitialization(SLCJdetector);
	auto SLCJphysList = new SLCJPhysicsList(physicsListName);
	runManager->SetUserInitialization(SLCJphysList);

	checkpoint;
	// set aditional user action classes
	SLCJRunAction* SLCJrun = new SLCJRunAction(isGeantino);
	runManager->SetUserAction(SLCJrun);
	SLCJEventAction* SLCJevent = new SLCJEventAction(SLCJrun);
	runManager->SetUserAction(SLCJevent);
	SLCJSteppingAction* SLCJstep = new SLCJSteppingAction(SLCJevent);
	runManager->SetUserAction(SLCJstep);

	checkpoint;
	// set mandatory user action class
	SLCJPrimaryGeneratorAction*
		SLCJgun = new SLCJPrimaryGeneratorAction(SLCJrun);
	runManager->SetUserAction(SLCJgun);

	checkpoint;
	// initialize G4 kernel
	runManager->Initialize();

	// Visualization, if you choose to have it!
#ifdef G4VIS_USE
//   G4VisManager* visManager = new SLCJVisManager;
	std::unique_ptr<G4VisManager> visManager = std::make_unique<G4VisExecutive>();
	visManager->Initialize();

	G4VVisManager::GetConcreteInstance()->SetDrawOverlapsFlag(true);
#endif
	SLCJphysList->SetDefaultCutValue(cutValue);

	G4String pname = "e-";
	std::string runDirectoryName;
	std::filesystem::path runDirectoryPath;
	std::string paramString = std::format("{}_{}mm_{}cm3_{}MeV",
		SLCJphysList->getPhysicsListName(),
		SLCJphysList->GetCutValue(pname) / mm,
		SLCJdetector->getFoodVolume() / cm3,
		energy / MeV);
	// find first free index to create data dump directory
	for (int i = 0;; i++) {
		runDirectoryName = std::format("event_{}_{}",
			paramString,
			i);
		runDirectoryPath = resultsDirectoryPath / runDirectoryName;
		if (std::filesystem::exists(runDirectoryPath) && skipIfDataExists) {
			std::cout << "Data exists. Aborting simulation." << _endl_;
			return 0;
		}
		if (!std::filesystem::exists(runDirectoryPath) || dataOverwrite) {
			break;
		}
	}
	std::filesystem::create_directory(runDirectoryPath);
	SLCJdetector->saveDetails(runDirectoryPath);
	SLCJgun->setRunPath(runDirectoryPath);

	// change verbosity here if needed
	G4UImanager* UI = G4UImanager::GetUIpointer();
	UI->ApplyCommand("/run/verbose 0");      // Run level
	UI->ApplyCommand("/event/verbose 0");    // Event generation level
	UI->ApplyCommand("/tracking/verbose 0"); // Tracking level
	UI->ApplyCommand("/process/verbose 0");  // Physics processes level
	UI->ApplyCommand("/geometry/verbose 0"); // Geometry level

	// create paths to simulation data files
	auto partialFileName = std::format("event_{}_",
		paramString);
	auto eventEnergyDepositFileName = partialFileName + "energyDeposit.bin";
	auto eventEnergyDepositFilePath = runDirectoryPath / eventEnergyDepositFileName;
	SLCJrun->setEventFilePath(eventEnergyDepositFilePath);
	// start a run
	checkpoint;
	SLCJgun->getEnergy() = energy; //set energy for each run
	runManager->BeamOn(numberOfEvent);

	auto stop = std::chrono::high_resolution_clock::now();
	std::cout << double((stop - start).count()) / 1e9 << '\n';
	std::cout << runDirectoryPath << '\n';
	// job termination
	// uncomment the following line if you want VRML visualisation (it won't work properly because of complex geometry)
	//UI->ApplyCommand("/control/execute vis.mac");
	return 0;
}
