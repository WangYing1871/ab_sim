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
/// \file analysis/shared/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
// $Id: SteppingAction.cc 68015 2013-03-13 13:27:27Z gcosmo $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "TrackingAction.hh"

#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
// #include "TRandom.h"
#include <random>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//
SteppingAction *SteppingAction::fgInstance = 0;
SteppingAction *SteppingAction::Instance()
{
  // G4cout<<"<<------------SteppingAction::Instance()-------------------->>"<<G4endl;
  return fgInstance;
}
//
SteppingAction::SteppingAction(DetectorConstruction *det, EventAction *event, TrackingAction *tracking)
    : G4UserSteppingAction(),
      fDetector(det),
      fEventAction(event),
      fTrackingAction(tracking)
{
  // G4cout<<"<<------------SteppingAction::SteppingAction(DetectorConstruction* det,TrackingAction* tracking) -------------------->>"<<G4endl;
  fgInstance = this;
  fEdep = 0;
  fStepLen = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{
  //  G4cout<<"<<------------SteppingAction::~SteppingAction()-------------------->>"<<G4endl;
  fgInstance = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step *aStep)
{
  // G4cout<<"<<------------SteppingAction::UserSteppingAction(const G4Step* aStep)-------------------->>"<<G4endl;
  // get volume of the current step
  G4String preVolumeName = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  // if(preVolumeName != "collimation" && preVolumeName != "source" && preVolumeName != "Hole1" && preVolumeName != "PET" && preVolumeName != "Gas")
  //   return;

  // collect energy and track length step by step
  G4double edep = aStep->GetTotalEnergyDeposit()-aStep->GetNonIonizingEnergyDeposit();
  G4double steplen = aStep->GetStepLength();

  /*
  if(edep>0 && preVolumeName == "collimation"){
    fTrackingAction->AddEdep_CuCollimation(edep);
  }


  if(edep>0 && preVolumeName == "Hole1"){
    fTrackingAction->AddEdep_Hole1(edep);
  }
*/

  // if(edep>0 && preVolumeName == "PET"){
  //   fTrackingAction->AddEdep_PET(edep);
  // }
  // if(edep>0 && preVolumeName == "source"){
  //   fTrackingAction->AddEdep_Source(edep);
  // }

  // if(edep>0 && preVolumeName == "Gas"){
  //   fTrackingAction->AddEdep_Gas(edep);
  // }

  // add this steps' energy lost if it is in the Effective Gas volume
  if (preVolumeName == "Gas" || preVolumeName == "GasEff")
  {
    fEventAction->AddEdep_ScoringVolume(edep);
  }

  // NEW PART: if this step does NOT begin in the scoring volume, drops it
  if (fTrackingAction->GetSelectTrackFlag() && fTrackingAction->GetNotFilteredFlag())
  {
    // if the previous tracks are ALL NOT in the periphery volume and the anticoincident MM region
    if (/* preVolumeName == "Gas" ||  */preVolumeName == "GasEff2")
      fTrackingAction->SetNotFilteredFlag(false); // this means the track get into the frame
    if (preVolumeName == "Gas" || preVolumeName == "GasEff")
    {
      fTrackingAction->SetHitSVFlag(true); // this means this track get into the gas volume

      // //Add secondary particles' energy
      // G4TrackVector* trackvector = astep->GetSecondary();
      // for(int i=0;i<trackvector.size();i++){
      //   fTrackingAction->AddEdep_ScoringVolume(trackvector[i]->GetKineticEnergy());
      // }

      // record the step points of this track
      G4ThreeVector stepVertexPos = aStep->GetPreStepPoint()->GetPosition();
      // calculate the total energy transfer in this step, including the energy transferred to secondaries
      G4double TotalEneTransfer = aStep->GetPreStepPoint()->GetKineticEnergy() - aStep->GetPostStepPoint()->GetKineticEnergy();

      fTrackingAction->AddTracklen_ScoringVolume(steplen);
      fTrackingAction->AddTrackEdep_SV(TotalEneTransfer);

      // this part sum the edep until the length is long enough, then fill
      fEdep += TotalEneTransfer;
      fStepLen += steplen;
      if (fStepLen > 0.5 * mm || aStep->GetPostStepPoint()->GetKineticEnergy() == 0)
      {
        G4double dE_dx = fEdep / fStepLen;
        
        // Get the total track length until this step
        //  G4double TotalTrackLen = aStep->GetTrack()->GetTrackLength();

        if (stepVertexPos.x() > fTrackingAction->MaxPosition[0])
          fTrackingAction->MaxPosition[0] = stepVertexPos.x();
        if (stepVertexPos.x() < fTrackingAction->MinPosition[0])
          fTrackingAction->MinPosition[0] = stepVertexPos.x();
        if (stepVertexPos.y() > fTrackingAction->MaxPosition[1])
          fTrackingAction->MaxPosition[1] = stepVertexPos.y();
        if (stepVertexPos.y() < fTrackingAction->MinPosition[1])
          fTrackingAction->MinPosition[1] = stepVertexPos.y();
        if (stepVertexPos.z() > fTrackingAction->MaxPosition[2])
          fTrackingAction->MaxPosition[2] = stepVertexPos.z();
        if (stepVertexPos.z() < fTrackingAction->MinPosition[2])
          fTrackingAction->MinPosition[2] = stepVertexPos.z();

        if (dE_dx > (fTrackingAction->MaxEdep))
        {
          fTrackingAction->MaxEdep = dE_dx;
          fTrackingAction->MaxEdepPos = fTrackingAction->GetTrackLenInSV();
          fTrackingAction->MaxEdepPosZ = stepVertexPos.z();
        }

        // record step informations: the step points of this track, the tracklen & dE/dx, for plotting
        fTrackingAction->AddStepInfo(stepVertexPos, dE_dx, fTrackingAction->GetTrackLenInSV());

        // // added part on 2023.10.09, to implement the process of digitization
        // // using the XY plane projected position of the step and its energy deposition to get the hit strip and recored it.
        // if (fTrackingAction->fDigitization)
        // {
        //   DriftOneElectron(stepVertexPos, fEdep);
        // }

        fEdep = 0;    // reset
        fStepLen = 0; // reset
      }

      // added part on 2023.10.09, to implement the process of digitization
      // using the XY plane projected position of the step and its energy deposition to get the hit strip and recored it.
      if (fTrackingAction->fDigitization)
      {
        DriftOneElectron(stepVertexPos, TotalEneTransfer);
      }
    }
    else if (preVolumeName == "GasEff2")
    {
      // this means the step is in the anti-coincidence MM
      G4double TotalEneTransfer = aStep->GetPreStepPoint()->GetKineticEnergy() - aStep->GetPostStepPoint()->GetKineticEnergy();

      fTrackingAction->AddTracklen_MMVolume(steplen);
      fTrackingAction->AddTrackEdep_MM(TotalEneTransfer);
    }
  }
}

void SteppingAction::Reset()
{
}

void SteppingAction::DriftOneElectron(G4ThreeVector steppos, G4double edep)
{
  // G4cout << "===========>> starting DriftOneElectron()" << G4endl;
  // the diffusion effect during drifting can be added later
  double drift_distance = -steppos.z();
  if(drift_distance < 0) return;
  double sigma_l = dl*sqrt(drift_distance*0.1)*10;         //in mm
  double sigma_t = dt*sqrt(drift_distance*0.1)*10;         //in mm
  if(drift_distance>70.) {
    G4cout << "This maybe an error: drift distance is out of range!" << G4endl;
    G4cout << "Drift distance (before diffusion): " << drift_distance << G4endl;
    G4cout << "Step position: " << steppos.x() << ", " << steppos.y() << ", " << steppos.z() <<G4endl;
  }

  // TRandom *randomGenerator = new TRandom();
  std::default_random_engine generator;
  std::normal_distribution<double> distributionL(0.0, sigma_l);
  std::normal_distribution<double> distributionT(0.0, sigma_t);

  // double sigma_buf;
  double sigma1, sigma2;
    
  // sigma1 = randomGenerator->Gaus(0, sigma_l);
  // sigma1 = G4RandGauss::shoot(0, sigma_l);
  sigma1 = distributionL(generator);

  drift_distance += sigma1;
  
  // sigma1 = randomGenerator->Gaus(0, sigma_t);
  // sigma2 = randomGenerator->Gaus(0, sigma_t);
  // sigma1 = G4RandGauss::shoot(0, sigma_t);
  // sigma2 = G4RandGauss::shoot(0, sigma_t);
  sigma1 = distributionT(generator);
  sigma2 = distributionT(generator);

  double projected_x = steppos.x()+sigma1;
  double projected_y = steppos.y()+sigma2;

  G4RootAnalysisManager* analysis = G4RootAnalysisManager::Instance();
  analysis->FillH2(0,projected_x,projected_y);

  // double velocity = v_drift;

  // double x_inplane = projected_x + nch/2*chnwidth;
  // double y_inplane = projected_y + nch/2*chnwidth;

  // change the coordinate system for simplier mapping
  double x_rot = sqrt(2) / 2 * (projected_x - projected_y);
  double y_rot = sqrt(2) / 2 * (projected_x + projected_y);

  int x_pos = (int)((x_rot + 0.5 * sqrt(2) / 2 * chnwidth) / (sqrt(2) / 2 * chnwidth));
  if ((x_rot + 0.5 * sqrt(2) / 2 * chnwidth) < 0)
    x_pos = x_pos - 1;
  int y_pos = (int)((y_rot + 0.5 * sqrt(2) / 2 * chnwidth) / (sqrt(2) / 2 * chnwidth));
  if ((y_rot + 0.5 * sqrt(2) / 2 * chnwidth) < 0)
    y_pos = y_pos - 1;

  // int chn;
  
  if(!on_plane_spread){
    FillOneChnSignal(x_pos, y_pos, edep, drift_distance);
    // if ((x_pos + y_pos) % 2 == 0)
    // {
    //   // this means it's on a Y strip
    //   chn = (y_pos - x_pos) / 2 + nch / 2;
    //   if (chn >= 0 && chn < nch)
    //   {
    //     // G4cout << "X position: " << steppos.x() << ", Y position: " << steppos.y() << ", Y strip number:" << chn << G4endl;
    //     fTrackingAction->charge_Y[chn].push_back(edep*1e6/E_ion);
    //     fTrackingAction->time_Y[chn].push_back(drift_distance / velocity * 100); // mm/(cm/us)*100--> ns
    //     fTrackingAction->IsEmpty = false;
    //   }
    // }
    // else
    // {
    //   // this means it's on a X strip
    //   chn = (y_pos + x_pos - 1) / 2 + nch / 2;
    //   if (chn >= 0 && chn < nch)
    //   {
    //     // G4cout << "Edep: " << edep << ", ionized charges: " << edep*1e6/E_ion << G4endl;
    //     // G4cout << "X position: " << steppos.x() << ", Y position: " << steppos.y() << ", X strip number:" << chn << G4endl;
    //     fTrackingAction->charge_X[chn].push_back(edep*1e6/E_ion);
    //     fTrackingAction->time_X[chn].push_back(drift_distance / velocity * 100); // mm/(cm/us)*100--> ns
    //     fTrackingAction->IsEmpty = false;
    //   }
    // }
  }
  else{
    double x_center_pos = x_pos*chnwidth/sqrt(2);
    double y_center_pos = y_pos*chnwidth/sqrt(2);

    double dx_plus, dx_minus, dy_plus, dy_minus;
    double r_this, r_dx_plus, r_dx_minus, r_dy_plus, r_dy_minus, r_sum;

    // G4cout << "distance to center: " << fabs(x_rot-x_center_pos) << " and " << fabs(y_rot-y_center_pos) << G4endl;
    // if(fabs(x_rot-x_center_pos)>sqrt(2)/4*chnwidth || fabs(y_rot-y_center_pos)>sqrt(2)/4*chnwidth) {
    //   G4cout << "Error! The position is out of the reconstructed strip!" << G4endl;
    // }
    // try to dimulate the spreading of charge due to resistive layer
    if ((x_pos + y_pos) % 2 == 0)
    {
      // this means it's on a Y strip
      // this stip's central line function is: y-ypos*chnwidth/sqrt(2) = k(x-x_pos*chnwidth/sqrt(2))
      // chn = (y_pos - x_pos) / 2 + nch / 2;
      dx_plus = 0.5*chnwidth-(x_rot+y_rot-x_center_pos-y_center_pos)/sqrt(2);
      dx_minus = 0.5*chnwidth+(x_rot+y_rot-x_center_pos-y_center_pos)/sqrt(2);
      dy_plus = chnwidth+(x_rot-y_rot-x_center_pos+y_center_pos)/sqrt(2);
      dy_minus = chnwidth-(x_rot-y_rot-x_center_pos+y_center_pos)/sqrt(2);

      //decide charge sharing ratio of different strips, using a gaussian function
      r_this = 1;
      r_dx_plus = exp(-(pow(dx_plus, 2) / (2*pow(sigma_spread, 2))));
      r_dx_minus = exp(-(pow(dx_minus, 2) / (2*pow(sigma_spread, 2))));
      r_dy_plus = exp(-(pow(dy_plus, 2) / (2*pow(sigma_spread, 2))));
      r_dy_minus = exp(-(pow(dy_minus, 2) / (2*pow(sigma_spread, 2))));
      r_sum = r_this+r_dx_plus+r_dx_minus+r_dy_plus+r_dy_minus;

      //fill the signal for all five channels according to their ratios
      FillOneChnSignal(x_pos, y_pos, edep*r_this/r_sum, drift_distance);
      FillOneChnSignal(x_pos+1, y_pos, edep*r_dx_plus/r_sum, drift_distance);
      FillOneChnSignal(x_pos-1, y_pos, edep*r_dx_minus/r_sum, drift_distance);
      FillOneChnSignal(x_pos-1, y_pos+1, edep*r_dy_plus/r_sum, drift_distance);
      FillOneChnSignal(x_pos+1, y_pos-1, edep*r_dy_minus/r_sum, drift_distance);
    }
    else
    {
      // this means it's on a X strip
      // this stip's central line function is: y-ypos*chnwidth/sqrt(2) = k(x-x_pos*chnwidth/sqrt(2))
      // chn = (y_pos - x_pos) / 2 + nch / 2;
      dx_plus = chnwidth-(x_rot+y_rot-x_center_pos-y_center_pos)/sqrt(2);
      dx_minus = chnwidth+(x_rot+y_rot-x_center_pos-y_center_pos)/sqrt(2);
      dy_plus = 0.5*chnwidth+(x_rot-y_rot-x_center_pos+y_center_pos)/sqrt(2);
      dy_minus = 0.5*chnwidth-(x_rot-y_rot-x_center_pos+y_center_pos)/sqrt(2);

      //decide charge sharing ratio of different strips, using a gaussian function
      r_this = 1;
      r_dx_plus = exp(-(pow(dx_plus, 2) / (2*pow(sigma_spread, 2))));
      r_dx_minus = exp(-(pow(dx_minus, 2) / (2*pow(sigma_spread, 2))));
      r_dy_plus = exp(-(pow(dy_plus, 2) / (2*pow(sigma_spread, 2))));
      r_dy_minus = exp(-(pow(dy_minus, 2) / (2*pow(sigma_spread, 2))));
      r_sum = r_this+r_dx_plus+r_dx_minus+r_dy_plus+r_dy_minus;

      //fill the signal for all five channels according to their ratios
      FillOneChnSignal(x_pos, y_pos, edep*r_this/r_sum, drift_distance);
      FillOneChnSignal(x_pos+1, y_pos+1, edep*r_dx_plus/r_sum, drift_distance);
      FillOneChnSignal(x_pos-1, y_pos-1, edep*r_dx_minus/r_sum, drift_distance);
      FillOneChnSignal(x_pos, y_pos+1, edep*r_dy_plus/r_sum, drift_distance);
      FillOneChnSignal(x_pos, y_pos-1, edep*r_dy_minus/r_sum, drift_distance);
    }
  }
}

void SteppingAction::FillOneChnSignal(int x_pos, int y_pos, double edep, double drift_distance)
{
  int chn;
  double velocity = v_drift;
  if ((x_pos + y_pos) % 2 == 0)
    {
      // this means it's on a Y strip
      chn = (y_pos - x_pos) / 2 + nch / 2;
      if (chn >= 0 && chn < nch)
      {
        // G4cout << "This is a Y channel hit, channel number is: " << chn << ", energy deposition: " << edep << G4endl;
        fTrackingAction->charge_Y[chn].push_back(edep*1e6/E_ion);
        fTrackingAction->time_Y[chn].push_back(drift_distance / velocity * 100); // mm/(cm/us)*100--> ns
        fTrackingAction->IsEmpty = false;
      }
    }
    else
    {
      // this means it's on a X strip
      chn = (y_pos + x_pos - 1) / 2 + nch / 2;
      if (chn >= 0 && chn < nch)
      {
        // G4cout << "This is a X channel hit, channel number is: " << chn << ", energy deposition: " << edep << G4endl;
        fTrackingAction->charge_X[chn].push_back(edep*1e6/E_ion);
        fTrackingAction->time_X[chn].push_back(drift_distance / velocity * 100); // mm/(cm/us)*100--> ns
        fTrackingAction->IsEmpty = false;
      }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
