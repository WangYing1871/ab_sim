#include "G4UIcmdWithABool.hh"

#include "tracking_messenger.hh"
#include "tracking_action.hh"

tracking_messenger::tracking_messenger(tracking_action* v):
  G4UImessenger()
  ,m_track_action(v){
  m_track_cmd = new G4UIcmdWithABool("/ab_sim/fullChain",this);
  m_track_cmd->SetGuidance("allow full decay chain");
  m_track_cmd->SetParameterName("flag",true);
  m_track_cmd->SetDefaultValue(true);
}

tracking_messenger::~tracking_messenger(){
  delete m_track_cmd; }

/* override */
void tracking_messenger::SetNewValue(G4UIcommand* command
    ,G4String v){
  //if (command==m_track_cmd)
  //  m_track_action->set_full_chain()
}
