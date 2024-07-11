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
/// \file DBDecay/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
//
// $Id: PhysicsList.cc 78307 2013-12-11 10:55:57Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleTypes.hh"
#include "G4IonConstructor.hh"
#include "G4RadioactiveDecay.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4StepLimiter.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicsListHelper.hh"

#include "G4KleinNishinaModel.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"

#include "G4MuMultipleScattering.hh"

#include "G4hMultipleScattering.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

// particles

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// gamma 
#include "G4PhotoElectricEffect.hh" 
#include "G4LivermorePhotoElectricModel.hh" 
#include "G4ComptonScattering.hh" 
#include "G4LivermoreComptonModel.hh" 
#include "G4GammaConversion.hh" 
#include "G4LivermoreGammaConversionModel.hh" 
#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh" 
// e- 
#include "G4eIonisation.hh" 
#include "G4LivermoreIonisationModel.hh" 
#include "G4UniversalFluctuation.hh" 
#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh" 
// e+ 
#include "G4eplusAnnihilation.hh"

// mu

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

// hadrons, ions

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

	PhysicsList::PhysicsList()
: G4VUserPhysicsList()
{
	G4cout<<"<<------------PhysicsList::PhysicsList()-------------------->>"<<G4endl;
	//add new units for radioActive decays
	// 
	const G4double minute = 60*second;
	const G4double hour   = 60*minute;
	const G4double day    = 24*hour;
	const G4double year   = 365*day;
	new G4UnitDefinition("minute", "min", "Time", minute);
	new G4UnitDefinition("hour",   "h",   "Time", hour);
	new G4UnitDefinition("day",    "d",   "Time", day);
	new G4UnitDefinition("year",   "y",   "Time", year);        
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{ G4cout<<"<<------------PhysicsList::~PhysicsList()-------------------->>"<<G4endl;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PhysicsList::ConstructParticle()
{
	G4cout<<"<<------------PhysicsList::ConstructParticle()-------------------->>"<<G4endl;
	G4BosonConstructor  pBosonConstructor;
	pBosonConstructor.ConstructParticle();

	G4LeptonConstructor pLeptonConstructor;
	pLeptonConstructor.ConstructParticle();

	G4MesonConstructor pMesonConstructor;
	pMesonConstructor.ConstructParticle();

	G4BaryonConstructor pBaryonConstructor;
	pBaryonConstructor.ConstructParticle();

	G4IonConstructor pIonConstructor;
	pIonConstructor.ConstructParticle();

	G4ShortLivedConstructor pShortLivedConstructor;
	pShortLivedConstructor.ConstructParticle(); 
	/*
	// pseudo-particles
	G4Geantino::GeantinoDefinition();

	// gamma
	G4Gamma::GammaDefinition();

	// leptons
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();

	G4NeutrinoE::NeutrinoEDefinition();
	G4AntiNeutrinoE::AntiNeutrinoEDefinition();

	// baryons
	G4Proton::ProtonDefinition();
	G4Neutron::NeutronDefinition();  

	// ions
	G4IonConstructor iConstructor;
	iConstructor.ConstructParticle();  
	*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
	G4cout<<"<<------------PhysicsList::ConstructProcess()-------------------->>"<<G4endl;
	AddTransportation();

	G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();

	radioactiveDecay->SetARM(true);               //Atomic Rearangement

	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();  
	ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());

	// Deexcitation (in case of Atomic Rearangement)
	//
	G4UAtomicDeexcitation* de = new G4UAtomicDeexcitation();
	de->SetFluo(true);
	de->SetAuger(true);   
	de->SetPIXE(true);  
	G4LossTableManager::Instance()->SetAtomDeexcitation(de);  

	ConstructEMProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
	G4cout<<"<<------------PhysicsList::SetCuts()-------------------->>"<<G4endl;
	//SetCutsWithDefault();
	G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250*eV, 1*GeV);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructEMProcess()
{

	//revised according to examples/extended/medical/electronScattering/src/PhysicsListEmStandardGS.cc
	G4cout<<"<<------------PhysicsList::ConstructEMProcess()-------------------->>"<<G4endl;
	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

	// Add standard EM Processes
	//
	auto theParticleIterator=GetParticleIterator();
	theParticleIterator->reset();
	while( (*theParticleIterator)() ){
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();

		// Step limit applied to all particles:
    	pmanager->AddDiscreteProcess(new G4StepLimiter);

		//Applicability range for Livermore models
		//for higher energies, the Standard models are used   
		G4double highEnergyLimit = 1*GeV;
		if (particleName == "gamma") {
			// gamma         

			G4PhotoElectricEffect* phot = new G4PhotoElectricEffect();
			G4LivermorePhotoElectricModel* 
				photModel = new G4LivermorePhotoElectricModel();
			photModel->SetHighEnergyLimit(highEnergyLimit);
			phot->AddEmModel(0, photModel);
			pmanager->AddDiscreteProcess(phot);

			G4ComptonScattering* compt = new G4ComptonScattering();
			G4LivermoreComptonModel* 
				comptModel = new G4LivermoreComptonModel();
			comptModel->SetHighEnergyLimit(highEnergyLimit);
			compt->AddEmModel(0, comptModel);
			pmanager->AddDiscreteProcess(compt);

			G4GammaConversion* conv = new G4GammaConversion();
			G4LivermoreGammaConversionModel* 
				convModel = new G4LivermoreGammaConversionModel();
			convModel->SetHighEnergyLimit(highEnergyLimit);
			conv->AddEmModel(0, convModel);
			pmanager->AddDiscreteProcess(conv);

			G4RayleighScattering* rayl = new G4RayleighScattering();
			G4LivermoreRayleighModel* 
				raylModel = new G4LivermoreRayleighModel();
			raylModel->SetHighEnergyLimit(highEnergyLimit);
			rayl->AddEmModel(0, raylModel);
			pmanager->AddDiscreteProcess(rayl);

		} else if (particleName == "e-") {
			//electron
			ph->RegisterProcess(new G4eMultipleScattering(), particle);
			pmanager->AddProcess(new G4eMultipleScattering, -1, 1,1);

			G4eIonisation* eIoni = new G4eIonisation();
			G4LivermoreIonisationModel* 
				eIoniModel = new G4LivermoreIonisationModel();
			eIoniModel->SetHighEnergyLimit(highEnergyLimit); 
			eIoni->AddEmModel(0, eIoniModel, new G4UniversalFluctuation() );
			pmanager->AddProcess(eIoni,                   -1, 2, 2);

			G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
			G4LivermoreBremsstrahlungModel* 
				eBremModel = new G4LivermoreBremsstrahlungModel();
			eBremModel->SetHighEnergyLimit(highEnergyLimit);
			eBrem->AddEmModel(0, eBremModel);
			pmanager->AddProcess(eBrem,                   -1, 3, 3);

		} else if (particleName == "e+") {
			//positron
			ph->RegisterProcess(new G4eMultipleScattering(), particle);
			pmanager->AddProcess(new G4eMultipleScattering, -1, 1,1);
			pmanager->AddProcess(new G4eIonisation,       -1, 2, 2);
			pmanager->AddProcess(new G4eBremsstrahlung,   -1, 3, 3);
			pmanager->AddProcess(new G4eplusAnnihilation,  0,-1, 4);

		} else if( particleName == "mu+" || 
				particleName == "mu-"    ) {
			//muon  
			ph->RegisterProcess(new G4MuMultipleScattering(), particle);
			pmanager->AddProcess(new G4MuMultipleScattering, -1, 1,1); 
			pmanager->AddProcess(new G4MuIonisation,      -1, 2, 2);
			pmanager->AddProcess(new G4MuBremsstrahlung,  -1, 3, 3);
			pmanager->AddProcess(new G4MuPairProduction,  -1, 4, 4);
      

		} else if( particleName == "alpha" || particleName == "GenericIon" ) { 
			ph->RegisterProcess(new G4hMultipleScattering(), particle);
			pmanager->AddProcess(new G4hMultipleScattering,  -1, 1,1);
			pmanager->AddProcess(new G4ionIonisation,     -1, 2, 2);
			pmanager->AddProcess(new G4NuclearStopping,		-1, 3, -1);				//Add nuclear stopping process for alpha

		} else if ((!particle->IsShortLived()) &&
				(particle->GetPDGCharge() != 0.0) && 
				(particle->GetParticleName() != "chargedgeantino")) {
			//all others charged particles except geantino
			pmanager->AddProcess(new G4hMultipleScattering, -1,1,1);
			pmanager->AddProcess(new G4hIonisation,       -1, 2, 2);
			ph->RegisterProcess(new G4hMultipleScattering(), particle);
		}
	}
}
