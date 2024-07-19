#include "G4ParticleTypes.hh"
#include "G4IonConstructor.hh"
#include "G4RadioactiveDecay.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4UserLimits.hh"
#include "G4StepLimiter.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicsListHelper.hh"

#include "G4KleinNishinaModel.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"

#include "G4hMultipleScattering.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"
//particles
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
//gamma
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh" 
#include "G4GammaConversion.hh" 
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"
//e-
#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4UniversalFluctuation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
//e+
#include "G4eplusAnnihilation.hh"
//mu
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuMultipleScattering.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4SystemOfUnits.hh"

#include "physics_list.hh"

physics_list::physics_list(): G4VUserPhysicsList(){
	//const G4double minute = 60*second;
	//const G4double hour   = 60*minute;
	//const G4double day    = 24*hour;
	//const G4double year   = 365*day;
	//new G4UnitDefinition("minute", "min", "Time", minute);
	//new G4UnitDefinition("hour",   "h",   "Time", hour);
	//new G4UnitDefinition("day",    "d",   "Time", day);
	//new G4UnitDefinition("year",   "y",   "Time", year);
}

/*override*/
void physics_list::ConstructParticle(){
  G4BosonConstructor  BosonConstructor; 
  G4LeptonConstructor LeptonConstructor; 
  G4MesonConstructor  MesonConstructor;
  G4BaryonConstructor BaryonConstrcutor;
  G4IonConstructor    IonConstructor;
  G4ShortLivedConstructor ShortLivedConstructor;
  BosonConstructor.ConstructParticle();
  LeptonConstructor.ConstructParticle();
  MesonConstructor.ConstructParticle();
  BaryonConstrcutor.ConstructParticle();
  IonConstructor.ConstructParticle();
  ShortLivedConstructor.ConstructParticle(); }

/*override*/
void physics_list::ConstructProcess(){
  AddTransportation();
  auto* radioactive_decay = new G4RadioactiveDecay();
  radioactive_decay->SetARM(true);
  auto* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  ph->RegisterProcess(radioactive_decay,G4GenericIon::GenericIon());

  auto* atomic_deexcitation = new G4UAtomicDeexcitation();
  atomic_deexcitation->SetFluo(true);
  atomic_deexcitation->SetAuger(true);
  atomic_deexcitation->SetPIXE(true);
  G4LossTableManager::Instance()->SetAtomDeexcitation(atomic_deexcitation);

  construct_em_process();
}

void physics_list::construct_em_process(){
  auto* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  auto piter = GetParticleIterator();
  piter->reset();
  while((*piter)()){
    auto* p = piter->value();
    auto* pmanager = p->GetProcessManager();
    auto  pname = p->GetParticleName();
    pmanager->AddDiscreteProcess(new G4StepLimiter);
    double highenergy_limit = 1.*GeV;
    if (pname=="gamma"){
      auto* pe_eff = new G4PhotoElectricEffect();
      auto* pe_model = new G4LivermorePhotoElectricModel();
      pe_model->SetHighEnergyLimit(highenergy_limit);
      pe_eff->AddEmModel(0,pe_model);
      pmanager->AddDiscreteProcess(pe_eff);

      auto* compt = new G4ComptonScattering();
      auto* compt_model = new G4LivermoreComptonModel();
      compt_model->SetHighEnergyLimit(highenergy_limit);
      compt->AddEmModel(0,compt_model);
      pmanager->AddDiscreteProcess(compt);

      auto* conv = new G4GammaConversion();
      auto* conv_model = new G4LivermoreGammaConversionModel();
      conv_model->SetHighEnergyLimit(highenergy_limit);
      conv->AddEmModel(0,conv_model);
      pmanager->AddDiscreteProcess(conv);

      auto* rayl = new G4RayleighScattering();
      auto* rayl_model = new G4LivermoreRayleighModel();
      rayl_model->SetHighEnergyLimit(highenergy_limit);
      rayl->AddEmModel(0,rayl_model);
      pmanager->AddDiscreteProcess(rayl);}
    else if(pname=="e-"){
      ph->RegisterProcess(new G4eMultipleScattering(),p);
      pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);

      auto* e_ioni = new G4eIonisation();
      auto* e_ioni_model = new G4LivermoreIonisationModel();
      e_ioni_model->SetHighEnergyLimit(highenergy_limit);
      e_ioni->AddEmModel(0,e_ioni_model,new G4UniversalFluctuation());
      pmanager->AddProcess(e_ioni,-1,2,2);

      auto* e_brem = new G4eBremsstrahlung();
      auto* e_brem_model = new G4LivermoreBremsstrahlungModel();
      e_brem_model->SetHighEnergyLimit(highenergy_limit);
      e_brem->AddEmModel(0,e_brem_model);
      pmanager->AddProcess(e_brem,-1,3,3); }
    else if(pname=="e+"){
      ph->RegisterProcess(new G4eMultipleScattering(),p);
      pmanager->AddProcess(new G4eMultipleScattering,-1,1,1);
      pmanager->AddProcess(new G4eIonisation,-1,2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,-1,3,3);
      pmanager->AddProcess(new G4eplusAnnihilation,0,-1,4); }
    else if(pname=="mu+" || pname=="mu-"){
      ph->RegisterProcess(new G4MuMultipleScattering(),p);
      pmanager->AddProcess(new G4MuMultipleScattering,-1,1,1);
      pmanager->AddProcess(new G4MuIonisation,-1,2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung,-1,3,3);
      pmanager->AddProcess(new G4MuPairProduction,-1,4,4); }
    else if(pname=="alpha" || pname=="GenericIon"){
      ph->RegisterProcess(new G4hMultipleScattering(),p);
      pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);
      pmanager->AddProcess(new G4ionIonisation,-1,2,2);
      pmanager->AddProcess(new G4NuclearStopping,-1,3,-1); }
    else if((!p->IsShortLived()) &&
        (p->GetPDGCharge()!=0.0) &&
        (p->GetParticleName()!="chargedgeantino")){
      pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);
      pmanager->AddProcess(new G4hIonisation,-1,2,2);
      ph->RegisterProcess(new G4hMultipleScattering(),p); }
  }
}
