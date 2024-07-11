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
/// \file DBDecay/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 68030 2015-03-13 13:51:27Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
//#include "EventMessenger.hh"

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include <iomanip>
#include <HistoManager.hh>
#include <TTree.h>
#include "RunAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(HistoManager* histo)
    :G4UserEventAction(),
    fPrintModulo(1000000),fDecayChain(),fHistoManager_Event(histo)
{
    G4cout<<"<<------------EventAction::EventAction(HistoManager* histo)-------------------->>"<<G4endl;
    //fHistoManager_Event = new HistoManager();
    //fEventMessenger = new EventMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
    G4cout<<"<<------------EventAction::~EventAction()-------------------->>"<<G4endl;
    //delete fHistoManager_Event;
    //delete fEventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* evt)
{
    // G4cout<<"<<------------EventAction::BeginOfEventAction(const G4Event*)-------------------->>"<<G4endl;
    fDecayChain = " ";
    fHistoManager_Event->fParticleInfo.reset();
    fTrackEdepInSV = 0.;
    fCount = true;
    ntracks = 0;
    TrackParentID.clear();
    fEventStartPosition = evt->GetPrimaryVertex()->GetPosition();
    // G4cout<<"begin of event"<<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
    // G4cout<<"<<------------EventAction::EndOfEventAction(const G4Event* evt)-------------------->>"<<G4endl;
    G4int evtNb = evt->GetEventID(); 
    //printing survey
    //
    if (evtNb%fPrintModulo == 0) 
        G4cout << "\n end of event " << std::setw(6) << evtNb 
            << " :" + fDecayChain << G4endl;
    // G4cout<<"end of event "<<fHistoManager_Event->fParticleInfo.nTrack<<" "<<fHistoManager_Event->fParticleInfo.fTrackTime[0]<<G4endl;
    //Fill this event, only fill when there is at least one track!!
    if(fHistoManager_Event->fParticleInfo.nTrack>0 && fHistoManager_Event->AllFilesOutput) fHistoManager_Event->fNtuple->Fill();

    auto analysis = G4RootAnalysisManager::Instance();
    if(fCount){
        analysis->FillH1(9,fTrackEdepInSV);
    }

    // retrieve trajectory info
    G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if (trajectoryContainer) {
        n_trajectories = trajectoryContainer->entries();
        // G4cout << " There are " << n_trajectories << " trajectories stored in this event." << G4endl;
        //for each kept track, search for its parent particle's name
        if(TrackParentID.size()!=0){
            // G4cout << " Among them are " << TrackParentID.size() << " kept  trajectories." << G4endl;
            for(size_t i=0; i<TrackParentID.size(); i++) {
                G4int ParentParticleID;
                for(int j=0;j<n_trajectories;j++){
                    ParentParticleID = (*trajectoryContainer)[j]->GetTrackID();
                    if(ParentParticleID == TrackParentID[i]) {
                        // find the parent track
                        G4String ParentParticleName = (*trajectoryContainer)[j]->GetParticleName();
                        // G4cout << "The parent particle is: " << ParentParticleName << ", its ID is: " << ParentParticleID 
                            // << ", it's position in the vector is : " << j << G4endl;
                        if(strcmp(ParentParticleName,"gamma") == 0) G4cout << "This is a gamma induced track! "<< G4endl;
                        break;
                    }
                }
                // G4String ParentParticleName = (*trajectoryContainer)[TrackParentID[i]]->GetParticleName();
                // G4int ParentParticleID = (*trajectoryContainer)[TrackParentID[i]]->GetTrackID();
                // G4cout << "The parent particle is: " << ParentParticleName << ", its ID is: " << ParentParticleID 
                //         << ", which should equal to: " << TrackParentID[i] << G4endl;
                // if(strcmp(ParentParticleName,"gamma") == 0) G4cout << "This is a gamma induced track! "<< G4endl;
            }
        }
    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


