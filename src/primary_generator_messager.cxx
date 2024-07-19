#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"

#include "primary_generator_action.hh"
#include "primary_generator_messager.hh"

primary_generator_messager::primary_generator_messager(primary_generator_action* pga):
  m_primary_generator(pga){
  m_directory = new G4UIdirectory("/det/");
  m_directory->SetGuidance("NTD_Ge detector control commands");

  m_gen_cmd = new G4UIcmdWithAString("/det/generator",this);
  m_gen_cmd->SetGuidance("Select primary generator");
  m_gen_cmd->SetParameterName("generator", true);
  m_gen_cmd->SetDefaultValue("PYTHIA");
  m_gen_cmd->SetCandidates("PYTHIA fGParticleSource");
}

primary_generator_messager::~primary_generator_messager(){
  delete m_gen_cmd; delete m_directory; }

/* override */
void primary_generator_messager::SetNewValue(G4UIcommand* command, G4String v){
  if (command==m_gen_cmd)
    m_primary_generator->use_HEP_evt_generator(v=="PYTHIA");
}

G4String primary_generator_messager::GetCurrentValue(G4UIcommand* command){
  if (command==m_gen_cmd)
    return m_primary_generator->is_HEP_evt_generator() ? "PYTHIA" : "User";
  return {};
}
