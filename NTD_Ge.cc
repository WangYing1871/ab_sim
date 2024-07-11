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
/// \file DBDecay.cc
/////
// $Id: DBDecay.cc  2015-09-10 19:07:37  $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// #ifdef G4MULTITHREADED
// #include "G4MTRunManager.hh"
// #else
#include "G4RunManager.hh"
// #endif

#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
//#include "ActionInitialization.hh"
#include "SteppingVerbose.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "TrackingAction.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

int main(int argc,char** argv) {
 
  //choose the Random engine
  G4cout<<"<<------------choose the Random engine-------------------->>"<<G4endl;
  G4cout<<"<<------------Number of threads: "<< G4Threading::G4GetNumberOfCores() <<G4endl;
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
   G4long seed = time(NULL);
  CLHEP::HepRandom::setTheSeed(seed);


  // Construct the default run manager
// #ifdef G4MULTITHREADED
//   G4MTRunManager* runManager = new G4MTRunManager;
//   runManager->SetNumberOfThreads(G4Threading::G4GetNumberOfCores());
// #else
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
  G4RunManager* runManager = new G4RunManager;
// #endif  


  // set mandatory initialization classes
  //
  G4cout<<"<<------------set mandatory initialization classes-------------------->>"<<G4endl;
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);
  PhysicsList* physics = new PhysicsList;
  //runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserInitialization(physics);
  // set user action classes
  // 
  /*
  ActionInitialization* actionInitialization = new ActionInitialization(detector);
  runManager->SetUserInitialization(actionInitialization);  
  */
  
  PrimaryGeneratorAction* primary = new PrimaryGeneratorAction(detector);
  runManager->SetUserAction(primary);
  
  HistoManager* histo = new HistoManager();
  //SteppingVerbose* stepV = new SteppingVerbose();
  
  RunAction* runAction = new RunAction(primary,histo);
  runManager->SetUserAction(runAction);
  
  
  EventAction* eventAction = new EventAction(histo);
  runManager->SetUserAction(eventAction);
  
  TrackingAction* trackingAction = new TrackingAction(runAction,eventAction,histo);
  runManager->SetUserAction(trackingAction);
  
  SteppingAction* steppingAction = new SteppingAction(detector,eventAction,trackingAction);
  runManager->SetUserAction(steppingAction);
  
  //Initialize G4 kernel
  G4cout<<"<<------------Initialize G4 kernel-------------------->>"<<G4endl;
  runManager->Initialize();
    
  // get the pointer to the User Interface manager 
  G4cout<<"<<------------get the pointer to the User Interface manager-------------------->>"<<G4endl;

  G4VisManager* visManager = new G4VisExecutive(argc, argv);
  visManager->Initialize();

  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc!=1)   // batch mode  
    { 
      G4cout<<"<<------------batch mode-------------------->>"<<G4endl;
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);  
    }
    
  else           // define visualization and UI terminal for interactive mode 
    { 
      G4cout<<"<<------------define visualization and UI terminal for interactive mode-------------------->>"<<G4endl;
     G4UIExecutive * ui = new G4UIExecutive(argc,argv);      
     UI->ApplyCommand("/control/execute vis.mac");          
     ui->SessionStart();
     delete ui;
    }

  delete visManager;

  // job termination
  //
  G4cout<<"<<------------job termination-------------------->>"<<G4endl;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo..... 

