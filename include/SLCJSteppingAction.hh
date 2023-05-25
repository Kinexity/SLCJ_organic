// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02SteppingAction.hh,v 1.3 1999/12/15 14:49:20 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//
//////////////////////////////////////////////////////////////////////
//
//  V. Guadilla 2021 
//  Get energy deposited in sensitive volume
//


#ifndef SLCJSteppingAction_h
#define SLCJSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include <string>

class SLCJEventAction;

class SLCJSteppingAction : public G4UserSteppingAction
{
public:
	SLCJSteppingAction(SLCJEventAction*);
	~SLCJSteppingAction() = default;

	void UserSteppingAction(const G4Step*);

private:
	SLCJEventAction* eventAction;
};

#endif
