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
/// \file DBDecay/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorAction.cc 78307 2013-12-11 10:55:57Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//==========================
//ordinary setup



#include "PrimaryGeneratorAction.hh"

#if simulation_type != 2

#include "PrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "G4HEPEvtInterface.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "G4ios.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
 : G4VUserPrimaryGeneratorAction(),
   fGParticleSource(),
   fDetector(det)
{
  G4cout<<"<<------------PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)-------------------->>"<<G4endl;
  const char* filename = "pythia_event.data";
  HEPEvt = new G4HEPEvtInterface(filename);
  fGParticleSource  = new G4GeneralParticleSource();

  messenger = new PrimaryGeneratorMessenger(this);
  useHEPEvt = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  G4cout<<"<<------------PrimaryGeneratorAction::~PrimaryGeneratorAction()-------------------->>"<<G4endl;
  delete HEPEvt;
  delete fGParticleSource;
  delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // G4cout<<"<<------------PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)-------------------->>"<<G4endl;
  // set random particle  position
  //  
  G4double x0 = 1.*cm,   y0 = 1.*cm,   z0 = 1.*cm;
  G4double dx0= 1.*cm, dy0= 1.*cm, dz0= 0.1*cm;
  x0 = dx0*(G4UniformRand()-0.5);
  y0 = dy0*(G4UniformRand()-0.5);
  z0 = dz0*(G4UniformRand()-0.5);
  HEPEvt->SetParticlePosition(G4ThreeVector(x0,y0,z0));
//  HEPEvt->SetParticlePosition(G4ThreeVector(0,0,0));

  // create vertex
  //
  if(useHEPEvt)
    { HEPEvt->GeneratePrimaryVertex(anEvent);}
  else
   { fGParticleSource->GeneratePrimaryVertex(anEvent);}
//     G4cout<<fGParticleSource->GetParticlePosition()<<G4endl;}

   
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#else


//==================================================
// for muon generation

#include "PrimaryGeneratorMessenger.hh"
#include "G4Event.hh"
#include "G4HEPEvtInterface.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "EcoMug.h"

#include <stdio.h>

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0),
   mu_plus(0),
   mu_minus(0),
   fDetector(det)
{
  G4cout<<"<<------------PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)-------------------->>"<<G4endl;
  // const char* filename = "pythia_event.data";
  // HEPEvt = new G4HEPEvtInterface(filename);
  fMuonGen.SetUseSky();
  fMuonGen.SetSkySize({{300.*mm, 300.*mm}});
  fMuonGen.SetSkyCenterPosition({{0., 0., 15.}});
  fMuonGen.CalulateFlux();

  fParticleGun  = new G4ParticleGun(1);
  mu_plus = G4ParticleTable::GetParticleTable()->FindParticle("mu+");
  mu_minus = G4ParticleTable::GetParticleTable()->FindParticle("mu-");

  messenger = new PrimaryGeneratorMessenger(this);
  // useHEPEvt = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  G4cout<<"<<------------PrimaryGeneratorAction::~PrimaryGeneratorAction()-------------------->>"<<G4endl;
  delete HEPEvt;
  delete fParticleGun;
  delete messenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // G4cout<<"<<------------PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)-------------------->>"<<G4endl;
  // set random particle  position
  //  
  // G4double x0 = 1.*cm,   y0 = 1.*cm,   z0 = 1.*cm;
  // G4double dx0= 1.*cm, dy0= 1.*cm, dz0= 0.1*cm;
  // x0 = dx0*(G4UniformRand()-0.5);
  // y0 = dy0*(G4UniformRand()-0.5);
  // z0 = dz0*(G4UniformRand()-0.5);
  // HEPEvt->SetParticlePosition(G4ThreeVector(x0,y0,z0));
//  HEPEvt->SetParticlePosition(G4ThreeVector(0,0,0));

  // create vertex
  //
  // if(useHEPEvt)
  //   { HEPEvt->GeneratePrimaryVertex(anEvent);}
  // else
  //  { fGParticleSource->GeneratePrimaryVertex(anEvent);}
//     G4cout<<fGParticleSource->GetParticlePosition()<<G4endl;}

  fMuonGen.Generate();
  array<double,3>  muon_pos = fMuonGen.GetGenerationPosition();
  double muon_ptot = fMuonGen.GetGenerationMomentum();
  double muon_theta = fMuonGen.GetGenerationTheta();
  double muon_phi = fMuonGen.GetGenerationPhi();
  double muon_charge = fMuonGen.GetCharge();

  fParticleGun->SetParticlePosition(G4ThreeVector(
                                              muon_pos[0]*mm,
                                              muon_pos[1]*mm,
                                              muon_pos[2]*mm
                                              ));

  fParticleGun->SetParticleMomentum(G4ParticleMomentum(
                                        muon_ptot*sin(muon_theta)*cos(muon_phi)*GeV,
                                        muon_ptot*sin(muon_theta)*sin(muon_phi)*GeV,
                                        muon_ptot*cos(muon_theta)*GeV
                                        ));

  if (muon_charge > 0){
      fParticleGun->SetParticleDefinition(mu_plus);
  }
  else
  {
      fParticleGun->SetParticleDefinition(mu_minus);
  }

  fParticleGun->GeneratePrimaryVertex(anEvent);
   
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
