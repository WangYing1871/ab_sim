#ifndef primary_generator_action_HH
#define primary_generator_action_HH 1 

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4VPrimaryGenerator;

class detector_construction;
class primary_generator_messager;

class primary_generator_action : public G4VUserPrimaryGeneratorAction{
public:
  typedef primary_generator_action self_t;
  typedef G4VUserPrimaryGeneratorAction base_t;
  primary_generator_action(detector_construction*);
  ~primary_generator_action();

public:
  virtual void GeneratePrimaries(G4Event*) override;
  G4GeneralParticleSource* GetParticleGun() const{return m_gparticle_source;}

private:
  G4VPrimaryGenerator* m_HEP_evt;
  primary_generator_messager* m_message;
  bool _is_use_HEP_evt;
  G4GeneralParticleSource* m_gparticle_source;
  detector_construction* m_detector;

public:
  inline void use_HEP_evt_generator(bool f) {_is_use_HEP_evt = f;}
  inline bool is_HEP_evt_generator() {return _is_use_HEP_evt;}
  inline void set_detector(detector_construction* v) {m_detector = v;}
  
};
//---------------------------------------------------------------------
#include "EcoMug.h"
class primary_generator_muon : public G4VUserPrimaryGeneratorAction{ 
public:
  typedef primary_generator_muon self_t;
  typedef G4VUserPrimaryGeneratorAction base_t;

public:
  primary_generator_muon(detector_construction*);
  ~primary_generator_muon();

  virtual void GeneratePrimaries(G4Event*) override;
  G4ParticleGun const* GetParticleGun() {return m_particle_gun;}

private:
  G4VPrimaryGenerator* m_HEP_evt;
  bool _is_use_HEP_evt;
  detector_construction* m_detector;

  G4ParticleGun* m_particle_gun;
  G4ParticleDefinition* m_mu_plus,* m_mu_minus;
  EcoMug m_muon_gen;

public:
  inline void use_HEP_evt_generator(bool f) {_is_use_HEP_evt = f;}
  inline bool is_HEP_evt_generator() {return _is_use_HEP_evt;}
};



#endif 
