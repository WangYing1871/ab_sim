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
/// \file DBDecay/include/EventAction.hh
/// \brief Definition of the EventAction class
//
// $Id: EventAction.hh 68017 2013-03-13 13:29:53Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "HistoManager.hh"
class EventMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
public:
  EventAction(HistoManager*);
  ~EventAction();
  
public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void   EndOfEventAction(const G4Event*);
  
  void SetPrintModulo(G4int val)   {fPrintModulo = val;};
  void AddDecayChain(G4String val) {fDecayChain += val;};
  void UpdateParticleInfo(ParticleInfo* fParticleInfo_Event){
    fHistoManager_Event->fParticleInfo.nTrack++;
    
    // fHistoManager_Event->fParticleInfo.fParticleZ.push_back(fParticleInfo_Event->fParticleZ.back());
    // fHistoManager_Event->fParticleInfo.fParticleA.push_back(fParticleInfo_Event->fParticleA.back());
    // fHistoManager_Event->fParticleInfo.fParticleMass.push_back(fParticleInfo_Event->fParticleMass.back());
    // fHistoManager_Event->fParticleInfo.fParticleCode.push_back(fParticleInfo_Event->fParticleCode.back());
    // fHistoManager_Event->fParticleInfo.fTrackKineticEnergy.push_back(fParticleInfo_Event->fTrackKineticEnergy.back());
    // fHistoManager_Event->fParticleInfo.fTrackVertexPosX.push_back(fParticleInfo_Event->fTrackVertexPosX.back());
    // fHistoManager_Event->fParticleInfo.fTrackVertexPosY.push_back(fParticleInfo_Event->fTrackVertexPosY.back());
    // fHistoManager_Event->fParticleInfo.fTrackVertexPosZ.push_back(fParticleInfo_Event->fTrackVertexPosZ.back());
    // fHistoManager_Event->fParticleInfo.fTrackVertexDirX.push_back(fParticleInfo_Event->fTrackVertexDirX.back());
    // fHistoManager_Event->fParticleInfo.fTrackVertexDirY.push_back(fParticleInfo_Event->fTrackVertexDirY.back());
    // fHistoManager_Event->fParticleInfo.fTrackVertexDirZ.push_back(fParticleInfo_Event->fTrackVertexDirZ.back());
    // fHistoManager_Event->fParticleInfo.fTrackID.push_back(fParticleInfo_Event->fTrackID.back());
    // fHistoManager_Event->fParticleInfo.fParentID.push_back(fParticleInfo_Event->fParentID.back());
    

    // //fHistoManager_Event->fParticleInfo.fTrackEdepInCuCollimation.push_back(fParticleInfo_Event->fTrackEdepInCuCollimation.back());
    // // fHistoManager_Event->fParticleInfo.fTrackEdepInSource.push_back(fParticleInfo_Event->fTrackEdepInSource.back());
    // //fHistoManager_Event->fParticleInfo.fTrackEdepInHole1.push_back(fParticleInfo_Event->fTrackEdepInHole1.back());
    // // fHistoManager_Event->fParticleInfo.fTrackEdepInPET.push_back(fParticleInfo_Event->fTrackEdepInPET.back());
    // // fHistoManager_Event->fParticleInfo.fTrackEdepInGas.push_back(fParticleInfo_Event->fTrackEdepInGas.back());
    // fHistoManager_Event->fParticleInfo.fTrackTime.push_back(fParticleInfo_Event->fTrackTime.back());

    //new branches
    fHistoManager_Event->fParticleInfo.fStartPosX.push_back(fParticleInfo_Event->fStartPosX.back());
    fHistoManager_Event->fParticleInfo.fStartPosY.push_back(fParticleInfo_Event->fStartPosY.back());
    fHistoManager_Event->fParticleInfo.fTrackLength.push_back(fParticleInfo_Event->fTrackLength.back());
    fHistoManager_Event->fParticleInfo.fTrackLength_MM.push_back(fParticleInfo_Event->fTrackLength_MM.back());
    fHistoManager_Event->fParticleInfo.fDriftDistance.push_back(fParticleInfo_Event->fDriftDistance.back());
    fHistoManager_Event->fParticleInfo.fXTrackLength.push_back(fParticleInfo_Event->fXTrackLength.back());
    fHistoManager_Event->fParticleInfo.fYTrackLength.push_back(fParticleInfo_Event->fYTrackLength.back());
    fHistoManager_Event->fParticleInfo.fCotTheta.push_back(fParticleInfo_Event->fCotTheta.back());
    fHistoManager_Event->fParticleInfo.fTrackEnergy.push_back(fParticleInfo_Event->fTrackEnergy.back());
    fHistoManager_Event->fParticleInfo.fTrackEnergy_MM.push_back(fParticleInfo_Event->fTrackEnergy_MM.back());
    fHistoManager_Event->fParticleInfo.fMaxEdepPosition.push_back(fParticleInfo_Event->fMaxEdepPosition.back());
    fHistoManager_Event->fParticleInfo.fMaxEdepPositionZ.push_back(fParticleInfo_Event->fMaxEdepPositionZ.back());

    fHistoManager_Event->fParticleInfo.fTrackStartPos.push_back(fParticleInfo_Event->fTrackStartPos.back());
    fHistoManager_Event->fParticleInfo.fEventStartPos.push_back(fParticleInfo_Event->fEventStartPos.back());
    fHistoManager_Event->fParticleInfo.fTrackVertexPosX_per_step.push_back(fParticleInfo_Event->fTrackVertexPosX_per_step.back());
    fHistoManager_Event->fParticleInfo.fTrackVertexPosY_per_step.push_back(fParticleInfo_Event->fTrackVertexPosY_per_step.back());
    fHistoManager_Event->fParticleInfo.fTrackVertexPosZ_per_step.push_back(fParticleInfo_Event->fTrackVertexPosZ_per_step.back());
    fHistoManager_Event->fParticleInfo.fTrackdE_dx_per_step.push_back(fParticleInfo_Event->fTrackdE_dx_per_step.back());
    fHistoManager_Event->fParticleInfo.fTrackLength_per_step.push_back(fParticleInfo_Event->fTrackLength_per_step.back());

  }
  void AddEdep_ScoringVolume(G4double edep){
    fTrackEdepInSV +=edep;
    return;
  }
  void SetCount(G4bool countflag){fCount=countflag;}
  G4bool GetCount(){return fCount;}
  void SetTime0(G4double t){fTime=t;}
  G4double GetTime0(){return fTime;}
  void ScoringVolumeTrackCounts(){ntracks++;}
  G4int Getntracks(){return ntracks;}
  void AddTrackParentID(G4int id){TrackParentID.push_back(id);}
  G4ThreeVector   fEventStartPosition;

  
  private:
    G4int           fPrintModulo;
    G4String        fDecayChain;  
    G4double        fTrackEdepInSV; 
    G4bool          fCount;          
    G4double        fTime;
    G4int           ntracks;      
    std::vector<G4int>     TrackParentID;        //this vector is used to fill the kept tracks per event, 10 is set that the track number is not larger than 10
  //EventMessenger* fEventMessenger;
  HistoManager* fHistoManager_Event;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
