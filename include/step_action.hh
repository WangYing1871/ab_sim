#ifndef step_action_HH
#define step_action_HH 1 
#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
class G4LogicalVolume;
class detector_construction;
class tracking_action;
class event_action;

class step_action : public G4UserSteppingAction{
public:
  typedef step_action self_t;
  typedef G4UserSteppingAction base_t;

  step_action();
  ~step_action();

  virtual void UserSteppingAction(G4Step const*) override;

  inline void set_detector(detector_construction* v) {m_detector = v;}
  inline void set_event_action(event_action* v) {m_event_action = v;}
  inline void set_tracking_action(tracking_action* v) {m_track_action = v;}

private:
  detector_construction* m_detector;
  event_action* m_event_action;
  tracking_action* m_track_action;
  G4LogicalVolume* m_scoring_volume_LV;
  double m_edep=0., m_step_len=0.;

  void reset();
  void drift_one_electron(G4ThreeVector const&, double edep);
  void fill_one_channel_signal(int,int,double,double);

public:
  static self_t* instance();

private:
  static self_t* m_s_step_action;

};


#endif
