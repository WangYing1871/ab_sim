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
/// \file analysis/shared/include/SteppingAction.hh
/// \brief Definition of the SteppingAction class
//
//
// $Id: SteppingAction.hh 68015 2013-03-13 13:27:27Z gcosmo $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "EventAction.hh"
#include "globals.hh"
#include "Constant.h"

class DetectorConstruction;
class TrackingAction;
class G4LogicalVolume;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(DetectorConstruction*, EventAction*, TrackingAction*);
  virtual ~SteppingAction();
  
  static SteppingAction* Instance();

  virtual void UserSteppingAction(const G4Step*);
    
  void Reset();

private:
  static SteppingAction* fgInstance;
  DetectorConstruction* fDetector;
  EventAction*             fEventAction;
  TrackingAction*          fTrackingAction;  
  G4LogicalVolume*         fScoringVolume;
  G4double                 fEdep;
  G4double                 fStepLen;
  //added on 2023.10.10, function to get the induced charge and its arrival time
  void DriftOneElectron(G4ThreeVector steppos, G4double edep);
  void FillOneChnSignal(int xpos, int ypos, double edep, double drift_distance);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
