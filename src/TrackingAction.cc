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
/// \file DBDecay/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//
// $Id: TrackingAction.cc 78307 2013-12-11 10:55:57Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include <HistoManager.hh>
//#include "Run.hh"
#include "EventAction.hh"
#include "TrackingMessenger.hh"
#include "G4SystemOfUnits.hh"

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "SteppingAction.hh"
#include "RunAction.hh"
#include "G4RunManager.hh"
#include <iostream>
#include <iomanip>
#include "math.h"
#include <random>
#include "Randomize.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace TPCsystem;

TrackingAction::TrackingAction(RunAction*runAct,EventAction* EA,HistoManager* histo)
    :G4UserTrackingAction(),
    fRun(runAct),fEvent(EA),fHistoManager_Track(histo),fTrackMessenger(0),
    fFullChain(true)
{
    G4cout<<"<<------------TrackingAction::TrackingAction(RunAction*runAct,EventAction* EA)-------------------->>"<<G4endl;
    fTrackMessenger = new TrackingMessenger(this);
    nCounts = 0;

    // charge_X.resize(nch);
    // time_X.resize(nch);
    // charge_Y.resize(nch);
    // time_Y.resize(nch);
    memset(waveform_X,0,sizeof(waveform_X));
    memset(waveform_Y,0,sizeof(waveform_Y));
    memset(waveform_X_int,0,sizeof(waveform_X_int));
    memset(waveform_Y_int,0,sizeof(waveform_Y_int));
    IsEmpty = true;
    //fSteppingVerbose_Tracking = new SteppingVerbose();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::~TrackingAction()
{
    G4cout<<"<<------------TrackingAction::~TrackingAction()-------------------->>"<<G4endl;
    delete fTrackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
    // G4cout<<"<<------------TrackingAction::PreUserTrackingAction(const G4Track* track)-------------------->>"<<G4endl;
    /*
       Run* run 
       = static_cast<Run*>(G4RunManager::GetRunManager()->GetCurrentRun());
       */   
    TrackStartPos.clear();          //reset the track start position
    EventStartPos.clear();

    //added on 2023.10.09, for digitization
    fDigitization = true;

    G4ParticleDefinition* particle = track->GetDefinition();
    G4String name   = particle->GetParticleName();
    fCharge = particle->GetPDGCharge();
    fBaryonNo = particle->GetBaryonNumber();
    fMass   = particle->GetPDGMass();  
    fParticleEnCode= particle->GetPDGEncoding();
    G4double Ekin = track->GetKineticEnergy();
    G4ThreeVector trackVertexPos = track->GetVertexPosition();
    G4ThreeVector trackVertexDir = track->GetVertexMomentumDirection();
    G4int ID      = track->GetTrackID();
    G4int parentID = track->GetParentID();
    volume_name = track->GetLogicalVolumeAtVertex()->GetName();
    if(ID!=1) {
        creator_process = track->GetCreatorProcess()->GetProcessName();
    }
    else {
        creator_process = "Primary";
    }


    //record the track start position
    TrackStartPos.push_back(trackVertexPos.x());
    TrackStartPos.push_back(trackVertexPos.y());
    TrackStartPos.push_back(trackVertexPos.z());

    //record the event start position
    EventStartPos.push_back(fEvent->fEventStartPosition.x());
    EventStartPos.push_back(fEvent->fEventStartPosition.y());
    EventStartPos.push_back(fEvent->fEventStartPosition.z());

    

    //fTrackEdepInCuCollimation=0.;
    //fTrackEdepInHole1=0.;
    // fTrackEdepInPET=0.;
    // fTrackEdepInSource=0.;
    // fTrackEdepInGas=0.;
    fTracklenInSV = 0.;
    fTracklenInMM = 0.;
    fTrackEnergyInSV = 0.;
    fTrackEnergyInMM = 0.;
    MinPosition[0]= 100*cm;
    MaxPosition[0] = -100*cm;
    MinPosition[1]= 100*cm;
    MaxPosition[1] = -100*cm;
    MinPosition[2]= 100*cm;
    MaxPosition[2] = -100*cm;
    MaxEdep = 0.;
    MaxEdepPos = 0.;
    MaxEdepPosZ = 0.;
    fSelectTrack = true;          //at the beginning of a track, set the default value to be true
    fNotFiltered=true;      //at the beginning of a track, set the default value to be true
    fHitSV=false;      //at the beginning of a track, set the default value to be false
    fTrackInfo_Stepping.reset();    //reset the recorded step points at each new track
    if(!IsEmpty) ClearChannelBuffer();           //clear the tracking digitization info of the former track

    //
    //-----------------update------------------------
    if(ID==1)
        fParticleInfo_Tracking.reset();
    // fParticleInfo_Tracking.nTrack=1;                 //this line is useless
    // fParticleInfo_Tracking.fParticleZ.push_back(fCharge);
    // fParticleInfo_Tracking.fParticleA.push_back(fBaryonNo);
    // fParticleInfo_Tracking.fParticleMass.push_back(fMass);
    // fParticleInfo_Tracking.fParticleCode.push_back(fParticleEnCode);
    // fParticleInfo_Tracking.fTrackKineticEnergy.push_back(Ekin);
    // fParticleInfo_Tracking.fTrackVertexPosX.push_back(trackVertexPos.x());
    // fParticleInfo_Tracking.fTrackVertexPosY.push_back(trackVertexPos.y());
    // fParticleInfo_Tracking.fTrackVertexPosZ.push_back(trackVertexPos.z());

    // fParticleInfo_Tracking.fTrackVertexDirX.push_back(trackVertexDir.x());
    // fParticleInfo_Tracking.fTrackVertexDirY.push_back(trackVertexDir.y());
    // fParticleInfo_Tracking.fTrackVertexDirZ.push_back(trackVertexDir.z());
    // fParticleInfo_Tracking.fTrackID.push_back(ID);
    // fParticleInfo_Tracking.fParentID.push_back(parentID);
    //------------------------------------------------
    //

    G4bool condition = false;

    //count particles
    //
    //if (ID>1) 
    fRun->ParticleCount(name, Ekin);

    //energy spectrum
    //
    G4int ih = 0;
    if (particle == G4Electron::Electron()||
            particle == G4Positron::Positron())  ih = 1;
    else if (particle == G4NeutrinoE::NeutrinoE()||
            particle == G4AntiNeutrinoE::AntiNeutrinoE()) ih = 2;
    else if (particle == G4Gamma::Gamma()) ih = 3;
    else if (particle == G4Alpha::Alpha()) ih = 4;
    else if (fCharge > 2.) ih = 5;
    if (ih) G4RootAnalysisManager::Instance()->FillH1(ih, Ekin);

    //fill the primary particle spectrum (used for direct particle source)
    if  (ID==1) {
        G4RootAnalysisManager::Instance()->FillH1(11, Ekin);
        E_primary = Ekin;
    }            
    //fill the primary beta spectrum (used for radioactive decay source(Co-60/Sr-90/...))
    // if  (ID<5 && parentID==1 && particle == G4Alpha::Alpha()) G4RootAnalysisManager::Instance()->FillH1(12, Ekin);        
    if  (parentID==2 && particle == G4Electron::Electron()) G4RootAnalysisManager::Instance()->FillH1(12, Ekin); 

    //Ion
    //
    if (fCharge > 2.) {
        //build decay chain
        if (ID == 1) fEvent->AddDecayChain(name);
        else       fEvent->AddDecayChain(" ---> " + name);
        // 
        //full chain: put at rest; if not: kill secondary      
        G4Track* tr = (G4Track*) track;
        if (fFullChain)  tr->SetTrackStatus(fStopButAlive);
        else if (ID>1) tr->SetTrackStatus(fStopAndKill);
    }

    //example of saving random number seed of this fEvent, under condition
    //
    ////condition = (ih == 3);
    if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent();

    //if the track start with a very low energy or this is not a electron track, drop this track
    KinEnergy_start = track->GetKineticEnergy();

    //Important: the following part should be switched for alpha/beta track information extracting!!!
    // if(KinEnergy_start < 10*keV || name != "alpha") fSelectTrack = false;
#if simulation_type == 0
    if(/* KinEnergy_start < 10*keV ||  */name != "e-") fSelectTrack = false;
    if(creator_process == "eIoni" || creator_process == "ionIoni"/*  || name != "gamma" */) fSelectTrack = false;                 //for Fe55 source simulation
#elif simulation_type ==1
    if(KinEnergy_start < 10*keV || name != "alpha") fSelectTrack = false;
#else
    if(name != "mu+" && name != "mu-") fSelectTrack = false;         // for cosmic ray simulation
#endif
    // if(KinEnergy_start < 10*keV || (name != "mu+" && name != "mu-")) fSelectTrack = false;         // for cosmic ray simulation
    // if(ID!=1) fSelectTrack = false;                //Test: only save primary beta tracks (for run12_2 only)!
    // if(parentID!=1 && parentID!=2) fSelectTrack = false;           //Test2: only save Sr90/Y90 decay primary betas (for Sr90 beta source run only)!

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
    // G4cout<<"<<------------TrackingAction::PostUserTrackingAction(const G4Track* track)-------------------->>"<<G4endl;
    //keep only ions
    //
    //if (fCharge < 3. ) return;
    /*
       Run* run 
       = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
       */
    G4RootAnalysisManager* analysis = G4RootAnalysisManager::Instance();

    //get time
    //   
    G4double time = track->GetGlobalTime();
    G4int ID = track->GetTrackID();
    if (ID == 1) fRun->PrimaryTiming(time);        //time of life of primary ion  

    //energy and momentum balance (from secondaries)
    //
    const std::vector<const G4Track*>* secondaries 
        = track->GetStep()->GetSecondaryInCurrentStep();
    size_t nbtrk = (*secondaries).size();
    if (nbtrk) {
        //there are secondaries --> it is a decay
        //
        //balance    
        G4double EkinTot = 0.;
        G4ThreeVector Pbalance = - track->GetMomentum();
        for (size_t itr=0; itr<nbtrk; itr++) {
            const G4Track* trk = (*secondaries)[itr];
            EkinTot += trk->GetKineticEnergy();
            //exclude gamma desexcitation from momentum balance
            if (trk->GetDefinition() != G4Gamma::Gamma())         
                Pbalance += trk->GetMomentum();                 
        }
        G4double Pbal = Pbalance.mag(); 
        fRun->Balance(EkinTot,Pbal);  
        analysis->FillH1(6,EkinTot);
        analysis->FillH1(7,Pbal);
    }

    //no secondaries --> end of chain    
    //  
    if (!nbtrk) {
        fRun->EventTiming(time);                     //total time of life
        analysis->FillH1(8,time);
    }

    //
    //------------------------
    // fParticleInfo_Tracking.fTrackEdepInSource.push_back(fTrackEdepInSource);
    //fParticleInfo_Tracking.fTrackEdepInCuCollimation.push_back(fTrackEdepInCuCollimation);
    //fParticleInfo_Tracking.fTrackEdepInHole1.push_back(fTrackEdepInHole1);
    // fParticleInfo_Tracking.fTrackEdepInPET.push_back(fTrackEdepInPET);
    // fParticleInfo_Tracking.fTrackEdepInGas.push_back(fTrackEdepInGas);
    // fParticleInfo_Tracking.fTrackTime.push_back(time);
    // G4cout<<time<<G4endl;
    // fEvent->UpdateParticleInfo(&fParticleInfo_Tracking);
    // G4cout<<"---GlobalTime="<<std::setprecision(18)<<time<<G4endl;

    //whether reject this track based on the former tracks in the event
    G4bool fReject = false;
    if(!fEvent->GetCount() && (std::abs(time-fEvent->GetTime0()) < 500)) fReject = true;


    G4String LastVolumeName = track->GetVolume()->GetName();
    if(LastVolumeName == "GasEff2") fNotFiltered = false;
    
    // generate random point to introduce the effect of efficiency of anticoincidence
    double inefficiency = G4UniformRand();

    if (fSelectTrack && (fNotFiltered || ( inefficiency < 0. )) && fHitSV/*  && LastVolumeName != "Gas"  */ && !fReject)
    {   //count conditions:
        //  this e- track does NOT cross the Al frame
        //  this track has at least part of it in the gas volume

        // if(!fEvent->GetCount()) {
        //     //this event is counted
        //     fRun->ScoringVolumeEventCounts();
        //     fEvent->SetCount(true);
        // }
        fEvent->ScoringVolumeTrackCounts();

        G4int nTracks = fEvent->Getntracks();
        G4int nEvents = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
        // nCounts++;
        // fHistoManager_Track->FillTrackGraph(&fTrackInfo_Stepping, nEvents, nTracks);

        // G4double EnergyInGas = fParticleInfo_Tracking.fTrackKineticEnergy.back();       //the beginning KinE of this track which equals to its energy lost in Gas
        // analysis->FillH1(9,fTrackEdepInSV);

        // G4double tracklength = track->GetTrackLength();
        analysis->FillH1(10,fTracklenInSV);

        //fill the x, y position of the start point
        G4double x_pos = fTrackInfo_Stepping.fStepVertexPosX[0];
        G4double y_pos = fTrackInfo_Stepping.fStepVertexPosY[0];
        

        //fill the tree branches

        if(MaxPosition[2]-MinPosition[2]>0.){
            analysis->FillH1(13, x_pos);
            analysis->FillH1(14, y_pos);
            analysis->FillH1(15, KinEnergy_start);

            fParticleInfo_Tracking.fStartPosX.push_back(x_pos);
            fParticleInfo_Tracking.fStartPosY.push_back(y_pos);
            fParticleInfo_Tracking.fTrackLength.push_back(fTracklenInSV);
            fParticleInfo_Tracking.fTrackLength_MM.push_back(fTracklenInMM);
            fParticleInfo_Tracking.fXTrackLength.push_back(MaxPosition[0]-MinPosition[0]);
            fParticleInfo_Tracking.fYTrackLength.push_back(MaxPosition[1]-MinPosition[1]);
            fParticleInfo_Tracking.fDriftDistance.push_back(MaxPosition[2]-MinPosition[2]);
            fParticleInfo_Tracking.fCotTheta.push_back((MaxPosition[2]-MinPosition[2])/(sqrt((MaxPosition[0]-MinPosition[0])*(MaxPosition[0]-MinPosition[0])+(MaxPosition[1]-MinPosition[1])*(MaxPosition[1]-MinPosition[1]))));
            fParticleInfo_Tracking.fTrackEnergy.push_back(fTrackEnergyInSV);
            fParticleInfo_Tracking.fTrackEnergy_MM.push_back(fTrackEnergyInMM);
            fParticleInfo_Tracking.fMaxEdepPosition.push_back(MaxEdepPos/fTracklenInSV);
            fParticleInfo_Tracking.fMaxEdepPositionZ.push_back((MaxPosition[2]-MaxEdepPosZ)/(MaxPosition[2]-MinPosition[2]));
            //Fill the track points information from the "fTrackInfo_Stepping" class
            fParticleInfo_Tracking.fTrackStartPos.push_back(TrackStartPos);
            fParticleInfo_Tracking.fEventStartPos.push_back(EventStartPos);
            fParticleInfo_Tracking.fTrackVertexPosX_per_step.push_back(fTrackInfo_Stepping.fStepVertexPosX);
            fParticleInfo_Tracking.fTrackVertexPosY_per_step.push_back(fTrackInfo_Stepping.fStepVertexPosY);
            fParticleInfo_Tracking.fTrackVertexPosZ_per_step.push_back(fTrackInfo_Stepping.fStepVertexPosZ);
            fParticleInfo_Tracking.fTrackdE_dx_per_step.push_back(fTrackInfo_Stepping.fStepdE_dx);
            fParticleInfo_Tracking.fTrackLength_per_step.push_back(fTrackInfo_Stepping.fStepTrackLen);

            if((MaxPosition[2]-MaxEdepPosZ)/(MaxPosition[2]-MinPosition[2])>=0){
                //Fill the first 100 tracks that fits the following condition
                fEvent->UpdateParticleInfo(&fParticleInfo_Tracking);
                if(nCounts<100)
                {
                    fHistoManager_Track->FillTrackGraph(&fTrackInfo_Stepping, nEvents, nTracks, nCounts);
                    nCounts++;
                }

                //added on 2023.10.10, to fill the digitized waveform and data
                if(fDigitization) {
                    // G4cout << "========== digitization begins =============" << G4endl;
                    FillChannelWaveforms();
                }
            }

            // if(TrackStartPos[2]>-10*mm){
            //     G4cout << "Start Position: " << TrackStartPos[0] << ", " << TrackStartPos[1] << ", " << TrackStartPos[2] << G4endl;
            //     G4cout << "The parent ID of this track: " << track->GetParentID() << " previous track ID: " << ParentID << G4endl;
            //     if(track->GetParentID() == ParentID) G4cout << "Previous track particle name: " << ParentTrackParticleName << G4endl;
            // }

            // G4cout  << "~~~~~~~ Volume name of this track's start position: " << volume_name << "\n" 
            //         << "------- creator process name of this track: " << creator_process << "\n" 
            //         << "------- initial kinetic energy of this track: " << KinEnergy_start << " MeV\n" 
            //         << "======= start position of this track: " << TrackStartPos[0] << ", " << TrackStartPos[1] << ", " << TrackStartPos[2] << "\n" << G4endl;
            
            
            if(volume_name == "PCB") {analysis->FillH1(16, KinEnergy_start);}
            else if(volume_name == "GasEff") {analysis->FillH1(17, KinEnergy_start);}
            else if(volume_name == "PET" || volume_name == "Poly" || volume_name == "PET2") {analysis->FillH1(18, KinEnergy_start);}
            else if(volume_name == "Cu_board") {analysis->FillH1(19, KinEnergy_start);}
            else {analysis->FillH1(20, KinEnergy_start);}

            //divided by different creator process
            if(creator_process == "compt") {analysis->FillH1(21, KinEnergy_start);}
            else if(creator_process == "phot") {analysis->FillH1(22, KinEnergy_start);}
            else if(creator_process == "eIoni") {analysis->FillH1(23, KinEnergy_start);}
            else {analysis->FillH1(24, KinEnergy_start);}
            

            //NEW PART: add the track parent id into the vector
            fEvent->AddTrackParentID(track->GetParentID());

            // fEvent->UpdateParticleInfo(&fParticleInfo_Tracking);

            // //added on 2023.10.10, to fill the digitized waveform and data
            // if(fDigitization) {
            //     // G4cout << "========== digitization begins =============" << G4endl;
            //     FillChannelWaveforms();
            // }
        }

    }
    else if(track->GetDefinition()->GetParticleName() == "e-") {
        //this means there's at least 1 previous electron track in this event that doesn't pass the selection
        if(track->GetGlobalTime()>fEvent->GetTime0()) {fEvent->SetTime0(track->GetGlobalTime());}
        fEvent->SetCount(false);
    }
    
    // // renew the previous track particle name 
    // ParentTrackParticleName = track->GetDefinition()->GetParticleName();
    // ParentID = track->GetTrackID();

}

void TrackingAction::FillChannelWaveforms()
{
    // G4cout << "============FillChannelWaveforms() begin ============= " << G4endl;
    memset(waveform_X,0,sizeof(waveform_X));
    memset(waveform_Y,0,sizeof(waveform_Y));
    memset(waveform_X_int,0,sizeof(waveform_X_int));
    memset(waveform_Y_int,0,sizeof(waveform_Y_int));
    // G4cout << "first value: " << time_X[0][0] << G4endl;
    // G4cout << "sizeof the vector each row: " << time_X[0].size() << G4endl;
    double minT_x = -1; // 

    // 遍历二维向量
    for (unsigned int row = 0; row < nch; ++row) {
        for (unsigned int col = 0; col < time_X[row].size(); ++col) {
            if(minT_x == -1) minT_x = time_X[row][col];
            if (time_X[row][col] < minT_x) {
                minT_x = time_X[row][col];
            }
        }
    }

    double minT_y = -1; // 

    // 遍历二维向量
    for (unsigned int row = 0; row < nch; ++row) {
        for (unsigned int col = 0; col < time_Y[row].size(); ++col) {
            if(minT_y == -1) minT_y = time_Y[row][col];
            if (time_Y[row][col] < minT_y) {
                minT_y = time_Y[row][col];
            }
        }
    }

    double minT = (minT_x<minT_y?minT_x:minT_y);

    for(int i=0;i<nch;i++){
        if(charge_X[i].size()!=time_X[i].size() || charge_Y[i].size()!=time_Y[i].size()){
            std::cout<<"Error: The size of charge and time vector is not equal!"<<std::endl;
            return;
        }
        if(charge_X[i].size()>0){
            bool isefficent = false;
            for(unsigned int j=0;j<charge_X[i].size();j++){
                
                int delta_T = (int)((time_X[i][j]-minT)/40);
                //this "100" term is to ensure the first arrived signal to sit near 100 time point
                if(delta_T>=(Nsp-100)) continue;
                int time_bin_start = -delta_T-100+peak_pos;
                //fill the waveform of this channel
                for(int ADC=0;ADC<Nsp;ADC++){
                    if(time_bin_start<0 || time_bin_start>(int)(sizeof(response_func)/sizeof(response_func[0])) || charge_X[i][j]==0) waveform_X[i][ADC]+=0;
                    else {
                        waveform_X[i][ADC] += 1.60218e-19/timebinwidth*gain*factor*1/CSAgain*charge_X[i][j]*response_func[time_bin_start];
                        isefficent = true;
                    }
                    time_bin_start++;
                }
            }
            //this process is to ensure the first point of a valid waveform is always non zero,
            // which is convinient for later filling data
            if(isefficent && waveform_X[i][0]==0) waveform_X[i][0]=-1;
            //save the waveform in int format also
            for(int k=0;k<Nsp;k++){
                waveform_X_int[i][k] = static_cast<int>(waveform_X[i][k]);
            }
            // if(isefficent) {
            //     int maxpoint = 0;
            //     int maxvalue_int = 0;
            //     double maxvalue_double = 0;
            //     for(int k=0;k<Nsp;k++) {
            //         if(waveform_X[i][k] > maxvalue_double) {
            //             maxpoint = k;
            //             maxvalue_double = waveform_X[i][k];
            //             maxvalue_int = waveform_X_int[i][k];
            //         }
            //     }
            //     // if(maxvalue_int ==0) {
            //     //     G4cout << "This waveform's maximum value is 0, it's maximum originally is: " << maxvalue_double << " at " << maxpoint << G4endl;
            //     // }
            // }
        }

        if(charge_Y[i].size()>0){
            bool isefficent = false;
            for(unsigned int j=0;j<charge_Y[i].size();j++){
                
                int delta_T = (int)((time_Y[i][j]-minT)/40);
                //this "100" term is to ensure the first arrived signal to sit near 100 time point
                if(delta_T>=(Nsp-100)) continue;
                int time_bin_start = -delta_T-100+peak_pos;
                //fill the waveform of this channel
                for(int ADC=0;ADC<Nsp;ADC++){
                    if(time_bin_start<0 || time_bin_start>(int)(sizeof(response_func)/sizeof(response_func[0])) || charge_Y[i][j]==0) waveform_Y[i][ADC]+=0;
                    else {
                        waveform_Y[i][ADC] += 1.60218e-19/timebinwidth*gain*factor*1/CSAgain*charge_Y[i][j]*response_func[time_bin_start];
                        isefficent = true;
                    }
                    time_bin_start++;
                }
            }
            //this process is to ensure the first point of a valid waveform is always non zero,
            // which is convinient for later filling data
            if(isefficent && waveform_Y[i][0]==0) waveform_Y[i][0]=-1;
            //save the waveform in int format also
            for(int k=0;k<Nsp;k++){
                waveform_Y_int[i][k] = static_cast<int>(waveform_Y[i][k]);
            }
            
            // if(isefficent) {
            //     int maxpoint = 0;
            //     int maxvalue_int = 0;
            //     double maxvalue_double = 0;
            //     for(int k=0;k<Nsp;k++) {
            //         if(waveform_Y[i][k] > maxvalue_double) {
            //             maxpoint = k;
            //             maxvalue_double = waveform_Y[i][k];
            //             maxvalue_int = waveform_Y_int[i][k];
            //         }
            //     }
            //     // if(maxvalue_int ==0) {
            //     //     G4cout << "This waveform's maximum value is 0, it's maximum originally is: " << maxvalue_double << " at " << maxpoint << G4endl;
            //     // }
            // }
        }
    }
    // //using float form of the waveform
    // fHistoManager_Track->SaveRawRootData(waveform_X,waveform_Y);

    G4ThreeVector position = G4ThreeVector(fTrackInfo_Stepping.fStepVertexPosX[0],fTrackInfo_Stepping.fStepVertexPosY[0],fTrackInfo_Stepping.fStepVertexPosZ[0]);
    //using int form of the waveform
    fHistoManager_Track->SaveRawRootData(waveform_X_int,waveform_Y_int,position,E_primary);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

