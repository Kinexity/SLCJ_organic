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

#include "F02ElectricFieldSetup.hh"

#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4DormandPrince745.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "fstream"

using namespace std;

// Global file with the events
ifstream geoInputFile2;

//  Constructor:

F02ElectricFieldSetup::F02ElectricFieldSetup() {

	//////////Reading the input data for primary generator///////////

	ifstream geoInputFile2("../geo.data");

	geoInputFile2 >> header1;
	geoInputFile2 >> gas[0] >> gas[1] >> gas[2];
	geoInputFile2 >> header2;
	geoInputFile2 >> fgas[0] >> fgas[1] >> fgas[2];
	geoInputFile2 >> header3;
	geoInputFile2 >> T >> P >> E;
	geoInputFile2.close();
	G4cout << E << " V/cm" << G4endl;
	//////////////////////////////////////////////////////////////////


	fEMfield = new G4UniformElectricField(G4ThreeVector(0.0, 0.0, 1e-3 * E * kilovolt / cm)); //ex 400V/cm

	// Create an equation of motion for this field
	fEquation = new G4EqMagElectricField(fEMfield);

	G4int nvar = 8;

	// Create the Runge-Kutta 'stepper' (different methods):
	//fStepper = new G4ClassicalRK4( fEquation, nvar ); //the default one
	//fStepper = new G4DormandPrince745( fEquation, nvar );
	fStepper = new G4CashKarpRKF45(fEquation, nvar);
	//fStepper = new G4ExplicitEuler( fEquation, nvar );
	//fStepper = new G4ImplicitEuler( fEquation, nvar );
	//fStepper = new G4SimpleRunge( fEquation, nvar );
	//fStepper = new G4SimpleHeum( fEquation, nvar );

	// Get the global field manager
	fFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();

	// Set this field to the global field manager
	fFieldManager->SetDetectorField(fEMfield);

	fMinStep = 0.010 * mm; // minimal step of 10 microns

	// The driver will ensure that integration is control to give acceptable integration error
	fIntgrDriver = new G4MagInt_Driver(fMinStep, fStepper, fStepper->GetNumberOfVariables());
	fChordFinder = new G4ChordFinder(fIntgrDriver);
	fFieldManager->SetChordFinder(fChordFinder);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//Destructor:

F02ElectricFieldSetup::~F02ElectricFieldSetup()
{
	delete fChordFinder;  fChordFinder = nullptr;
	delete fStepper;      fStepper = nullptr;
	delete fEquation;     fEquation = nullptr;
	delete fEMfield;      fEMfield = nullptr;
}

