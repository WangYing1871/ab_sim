#ifndef tracking_messenger_HH
#define tracking_messenger_HH 1 
#include "G4UImessenger.hh"
class G4UIcmdWithABool;
class tracking_action;

class tracking_messenger : public G4UImessenger{
public:
  typedef tracking_messenger self_t;
  typedef G4UImessenger base_t;

  tracking_messenger(tracking_action*);
  ~tracking_messenger();

  virtual void SetNewValue(G4UIcommand*,G4String) override;

private:
  tracking_action* m_track_action;
  G4UIcmdWithABool* m_track_cmd;
};
#endif
