#ifndef primary_generator_messager_HH
#define primary_generator_messager_HH 1 
#include "G4UImessenger.hh"
#include "globals.hh"
class G4UIdirectory;
class G4UIcmdWithAString;

class primary_generator_action;

class primary_generator_messager : public G4UImessenger{
public:
  primary_generator_messager(primary_generator_action*);
  ~primary_generator_messager();

public:
  typedef primary_generator_messager self_t;
  typedef G4UImessenger base_t;

public:
  virtual void SetNewValue(G4UIcommand*, G4String) override;
  virtual G4String GetCurrentValue(G4UIcommand*);

private:
  primary_generator_action* m_primary_generator;
  G4UIdirectory* m_directory;
  G4UIcmdWithAString* m_gen_cmd;
};
#endif
