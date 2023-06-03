// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02RunAction.hh,v 1.3 1999/12/15 14:49:20 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
//
/////////////////////////////////////////////////////////////////
//
//  V Guadilla 2021
//
/////////////////////////////////////////////////////////////////

#ifndef SLCJRunAction_h
#define SLCJRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <array>
#include <filesystem>
#include "G4Timer.hh"
#include <G4ThreeVector.hh>

class G4Run;

class G4Timer;

class SLCJRunAction : public G4UserRunAction
{
  public:
    SLCJRunAction();
    ~SLCJRunAction() = default;

  public:
    void BeginOfRunAction(const G4Run*);
    void EndOfRunAction(const G4Run*);


    void fillOut(std::vector<std::array<G4double, 4>>& EnergyDeposit, std::array<G4double, 20>& EnergyGammaCrystals);
    void fillOut(std::vector<std::pair<G4String, G4ThreeVector>>& geantinoPos);
    
    
    void setEventFilePath(std::filesystem::path totalP, std::filesystem::path stepsP);
private:
   
    std::unique_ptr<G4Timer> timer;
    
    std::fstream 
        eventTotalDepositFile,
        eventStepsDepositFile,
        eventTotalDepositFileBinary,
        eventStepsDepositFileBinary;
   
    uint32_t eventIndex;

    std::filesystem::path 
        eventTotalDepositFilePath,
        eventStepsDepositFilePath;

    std::fstream file;

    std::vector<std::pair<G4String, G4ThreeVector>> geantinoPosGlobal;
};

#endif






