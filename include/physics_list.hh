#ifndef physics_list_HH
#define physics_list_HH 1
#include "globals.hh"
#include "G4VUserPhysicsList.hh"
#include "G4SystemOfUnits.hh"

class physics_list : public G4VUserPhysicsList{
  typedef physics_list self_t;
  typedef G4VUserPhysicsList base_t;

public:
  physics_list();
  virtual ~physics_list() = default;

protected:

  virtual void ConstructParticle() override;
  virtual void ConstructProcess() override;
  virtual void SetCuts() override{
	  G4ProductionCutsTable::GetProductionCutsTable()
      ->SetEnergyRange(250*eV, 1*GeV); }

public:
  void construct_em_process();

};

#endif
