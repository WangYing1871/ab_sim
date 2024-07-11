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
//
// $Id: HistoManager.hh 68017 2013-03-13 13:29:53Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include <vector>
#include "G4RootAnalysisManager.hh"
#include <G4ThreeVector.hh>
#include "Constant.h"


using namespace std;
using namespace TPCsystem;

class TTree;
class TFile;
//const G4int kMAXTrack=5000;//should be checked!!!
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ParticleInfo
{
public:
  G4int nTrack;
  std::vector<G4double> fParticleZ;
  std::vector<G4double> fParticleA;
  std::vector<G4double> fParticleMass;
  std::vector<G4double> fParticleCode;
  std::vector<G4double> fTrackKineticEnergy;
  std::vector<G4double> fTrackVertexPosX;
  std::vector<G4double> fTrackVertexPosY;
  std::vector<G4double> fTrackVertexPosZ;
//  std::vector<G4ThreeVector> fTrackVertexPos;
  std::vector<G4double> fTrackVertexDirX;
  std::vector<G4double> fTrackVertexDirY;
  std::vector<G4double> fTrackVertexDirZ;
//  std::vector<G4ThreeVector> fTrackVertexDir;
  std::vector<G4double> fTrackID;
  std::vector<G4double> fParentID;
  // std::vector<G4double> fTrackEdepInSource;
  //std::vector<G4double> fTrackEdepInCuCollimation;
  //std::vector<G4double> fTrackEdepInHole1;
  // std::vector<G4double> fTrackEdepInPET;
  // std::vector<G4double> fTrackEdepInGas;
  std::vector<G4double> fTrackTime;

  //new variables
  std::vector<G4double> fStartPosX;
  std::vector<G4double> fStartPosY;
  std::vector<G4double> fTrackLength;
  std::vector<G4double> fTrackLength_MM;
  std::vector<G4double> fDriftDistance;
  std::vector<G4double> fXTrackLength;            //the track length projected onto x axis
  std::vector<G4double> fYTrackLength;            //the track length projected onto y axis
  std::vector<G4double> fCotTheta;                //cotangent of the incident angle
  std::vector<G4double> fTrackEnergy;
  std::vector<G4double> fTrackEnergy_MM;
  std::vector<G4double> fMaxEdepPosition;
  std::vector<G4double> fMaxEdepPositionZ;
  std::vector<std::vector<G4double> > fTrackStartPos;
  std::vector<std::vector<G4double> > fEventStartPos;
  std::vector<std::vector<G4double> > fTrackVertexPosX_per_step;
  std::vector<std::vector<G4double> > fTrackVertexPosY_per_step;
  std::vector<std::vector<G4double> > fTrackVertexPosZ_per_step;
  std::vector<std::vector<G4double> > fTrackdE_dx_per_step;
  std::vector<std::vector<G4double> > fTrackLength_per_step;
  
  void reset()
  {
    nTrack=0;
    fParticleZ.clear();
    fParticleA.clear();
    fParticleMass.clear();
    fParticleCode.clear();
    fTrackKineticEnergy.clear();
    fTrackVertexPosX.clear();
    fTrackVertexPosY.clear();
    fTrackVertexPosZ.clear();
    fTrackVertexDirX.clear();
    fTrackVertexDirY.clear();
    fTrackVertexDirZ.clear();
    fTrackID.clear();
    fParentID.clear();
    //fTrackEdepInCuCollimation.clear();
    //fTrackEdepInHole1.clear();
    // fTrackEdepInPET.clear();
    // fTrackEdepInSource.clear();
    // fTrackEdepInGas.clear();
    fTrackTime.clear();

    fStartPosX.clear();
    fStartPosY.clear();
    fTrackLength.clear();
    fTrackLength_MM.clear();
    fDriftDistance.clear();
    fXTrackLength.clear();
    fYTrackLength.clear();
    fCotTheta.clear();
    fTrackEnergy.clear();
    fTrackEnergy_MM.clear();
    fMaxEdepPosition.clear();
    fMaxEdepPositionZ.clear();

    fTrackStartPos.clear();
    fEventStartPos.clear();
    fTrackVertexPosX_per_step.clear();
    fTrackVertexPosY_per_step.clear();
    fTrackVertexPosZ_per_step.clear();
    fTrackdE_dx_per_step.clear();
    fTrackLength_per_step.clear();

  };

  

  ParticleInfo():
    nTrack(0){
    fParticleZ.clear();
    fParticleA.clear();
    fParticleMass.clear();
    fParticleCode.clear();
    fTrackKineticEnergy.clear();
    fTrackVertexPosX.clear();
    fTrackVertexPosY.clear();
    fTrackVertexPosZ.clear();
    fTrackVertexDirX.clear();
    fTrackVertexDirY.clear();
    fTrackVertexDirZ.clear();
    fTrackID.clear();
    fParentID.clear();
    //fTrackEdepInCuCollimation.clear();
    //fTrackEdepInHole1.clear();
    // fTrackEdepInPET.clear();
    // fTrackEdepInSource.clear();
    // fTrackEdepInGas.clear();
    fTrackTime.clear();

    fStartPosX.clear();
    fStartPosY.clear();
    fTrackLength.clear();
    fTrackLength_MM.clear();
    fDriftDistance.clear();
    fXTrackLength.clear();
    fYTrackLength.clear();
    fCotTheta.clear();
    fTrackEnergy.clear();
    fTrackEnergy_MM.clear();
    fMaxEdepPosition.clear();
    fMaxEdepPositionZ.clear();

    fTrackStartPos.clear();
    fEventStartPos.clear();
    fTrackVertexPosX_per_step.clear();
    fTrackVertexPosY_per_step.clear();
    fTrackVertexPosZ_per_step.clear();
    fTrackdE_dx_per_step.clear();
    fTrackLength_per_step.clear();

    }
};

class TrackInfo
{
public:
  std::vector<G4double> fStepVertexPosX;
  std::vector<G4double> fStepVertexPosY;
  std::vector<G4double> fStepVertexPosZ;
  std::vector<G4double> fStepdE_dx;
  std::vector<G4double> fStepTrackLen;
  
  void reset()
  {
    fStepVertexPosX.clear();
    fStepVertexPosY.clear();
    fStepVertexPosZ.clear();
    fStepdE_dx.clear();
    fStepTrackLen.clear();
  };

  TrackInfo(){
    fStepVertexPosX.clear();
    fStepVertexPosY.clear();
    fStepVertexPosZ.clear();
    fStepdE_dx.clear();
    fStepTrackLen.clear();
  }
};

//New class: added on 2023.10.09, variables for digitization output
class RawRootData
{
public:
  int event=0;
  int nHits=0;
  int Fec[Tch];
  int Chip[Tch];
  int Chn[Tch];
  int ADC[Tch][Nsp];
  float sumADC[Tch];
  float maxADC[Tch];
  float maxPoint[Tch];
  float summaxADC;
  int pixelX[Tch];
  int pixelY[Tch];
  float startpos_x;
  float startpos_y;
  float startpos_z;
  float kinE_start;
  
  void reset()
  {
    // event=0;             //event number is not reset in the entire run, so that the output event is like: 0, 1, 2, 3, ...
    nHits=0;
    memset(Fec,0,sizeof(Fec));
    memset(Chip,0,sizeof(Chip));
    memset(Chn,0,sizeof(Chn));
    memset(ADC,0,sizeof(ADC));
    memset(sumADC,0,sizeof(sumADC));
    memset(maxADC,0,sizeof(maxADC));
    memset(maxPoint,0,sizeof(maxPoint));
    summaxADC = 0;
    memset(pixelX,0,sizeof(pixelX));
    memset(pixelY,0,sizeof(pixelY));
    startpos_x = 0;
    startpos_y = 0;
    startpos_z = 0;
    kinE_start = 0;
  };

  RawRootData(){
    event=0;
    nHits=0;
    memset(Fec,0,sizeof(Fec));
    memset(Chip,0,sizeof(Chip));
    memset(Chn,0,sizeof(Chn));
    memset(ADC,0,sizeof(ADC));
    memset(sumADC,0,sizeof(sumADC));
    memset(maxADC,0,sizeof(maxADC));
    memset(maxPoint,0,sizeof(maxPoint));
    summaxADC = 0;
    memset(pixelX,0,sizeof(pixelX));
    memset(pixelY,0,sizeof(pixelY));
    startpos_x = 0;
    startpos_y = 0;
    startpos_z = 0;
    kinE_start = 0;
  }
};

class HistoManager
{
public:
   HistoManager();
  ~HistoManager();
  void save();
  void book();
  void FillTrackGraph(TrackInfo*, G4int, G4int, G4int);
  void SaveRawRootData(int waveform_X[Tch][Nsp], int waveform_Y[Tch][Nsp],G4ThreeVector position,G4double energy);
  ParticleInfo fParticleInfo;
  RawRootData fRawRootData;
private:
  
  void Book();
  G4String fFileName;
public:  
  TFile* fRootFile;
  TFile* fGraphRootFile;
  TFile* fRawRootFile;
  TTree* fNtuple;
  TTree* fNtuple2;
  G4bool AllFilesOutput;      //choose whether to generate output MC root file and track file
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

