#include <iostream>

#include "G4Event.hh"
#include "G4HEPEvtInterface.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "detector_construction.hh"
#include "primary_generator_action.hh"
#include "primary_generator_messager.hh"
#include "util.h"

primary_generator_action::primary_generator_action(detector_construction* det):
  G4VUserPrimaryGeneratorAction()
  ,m_detector(det){
  m_HEP_evt = new G4HEPEvtInterface(util::config::HEPEvt_Interface_filename);
  m_gparticle_source = new G4GeneralParticleSource();
  m_message = new primary_generator_messager(this);
  _is_use_HEP_evt = true;
    std::cout<<"wangying\n";
}

primary_generator_action::~primary_generator_action(){
  delete m_HEP_evt;
  delete m_gparticle_source;
  delete m_message;
}

/* override */
void primary_generator_action::GeneratePrimaries(G4Event* evt){
  double x0 = 1.*cm*(G4UniformRand()-.5);
  double y0 = 1.*cm*(G4UniformRand()-.5);
  double z0 = .5*cm*(G4UniformRand()-.5);
  m_HEP_evt->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  if (_is_use_HEP_evt) m_HEP_evt->GeneratePrimaryVertex(evt);
  else m_gparticle_source->GeneratePrimaryVertex(evt);
}
//....ooOO0OOoo........ooOO0OOoo=========ooOO0OOoo........ooOO0OOoo....

#include "EcoMug.h"
#include "G4Geantino.hh"
#include "G4IonTable.hh"

primary_generator_muon::primary_generator_muon(detector_construction* det):
  G4VUserPrimaryGeneratorAction()
  ,m_detector(det){
  m_muon_gen.SetUseSky();
  m_muon_gen.SetSkySize({300.*mm,300.*mm});
  m_muon_gen.SetSkyCenterPosition({0.,0.,15.});
  m_muon_gen.CalulateFlux();

  m_particle_gun = new G4ParticleGun(1);
  m_mu_plus = G4ParticleTable::GetParticleTable()->FindParticle("mu+");
  m_mu_minus = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
}

primary_generator_muon::~primary_generator_muon(){ delete m_HEP_evt; delete m_particle_gun; }

/* override */
void primary_generator_muon::GeneratePrimaries(G4Event* evt){
  m_muon_gen.Generate();
  auto muon_pos = m_muon_gen.GetGenerationPosition();
  double muon_ptot = m_muon_gen.GetGenerationMomentum();
  double muon_theta = m_muon_gen.GetGenerationTheta();
  double muon_phi = m_muon_gen.GetGenerationPhi();
  double muon_charge = m_muon_gen.GetCharge();

  m_particle_gun->SetParticlePosition(G4ThreeVector(
        muon_pos[0]*mm ,muon_pos[1]*mm ,muon_pos[1]*mm));
  m_particle_gun->SetParticleMomentum(G4ParticleMomentum(
        muon_ptot*sin(muon_theta)*cos(muon_phi)*GeV
        ,muon_ptot*sin(muon_theta)*sin(muon_phi)*GeV
        ,muon_ptot*cos(muon_theta)*GeV
        ));
  m_particle_gun->SetParticleDefinition(
      G4ParticleTable::GetParticleTable()->FindParticle(muon_charge>0 ? "mu+" : "mu-"));
  m_particle_gun->GeneratePrimaryVertex(evt); }
