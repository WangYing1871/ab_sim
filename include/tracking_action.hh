#ifndef tracking_action_HH
#define tracking_action_HH 1 
#include "G4UserTrackingAction.hh"
#include "globals.hh"
#include "util.h"

class step_action;
class event_action;
class tracking_messenger;

class tracking_action : public G4UserTrackingAction{

  virtual void PreUserTrackingAction(G4Track const*);
  virtual void PostUserTrackingAction(G4Track const*);

public:
  typedef tracking_action self_t;
  typedef G4UserTrackingAction base_t;

  tracking_action();
  ~tracking_action();

public:
  double max_position[3];
  double min_position[3];

  std::vector<double> m_charge_x[TPCsystem::nch];
  std::vector<double> m_time_x[TPCsystem::nch];
  std::vector<double> m_charge_y[TPCsystem::nch];
  std::vector<double> m_time_y[TPCsystem::nch];

  double m_waveform_x[TPCsystem::Tch][TPCsystem::Nsp];
  double m_waveform_y[TPCsystem::Tch][TPCsystem::Nsp];

  bool _is_not_filter;
  bool _is_select_track = true;
  bool _is_hit_sv;
  bool _is_empty = true;
  bool _is_digitalize = false;

  double m_max_dedx;
  double m_max_dedx_posz;

  int m_count;

  std::string m_create_process_name;

  bool is_not_filtered() const {return _is_not_filter;}
  bool is_select_track() const {return _is_select_track;}
  bool is_hit_sv() const {return _is_hit_sv;}
  void hit_sv(bool v) {_is_hit_sv = v;}
  void not_filtered(bool v) {_is_not_filter = v;}
 
private:
  typedef G4ThreeVector pos_t;


  event_action* m_event;
  tracking_messenger* m_message;

  util::track_info m_track_info;
  //util::particle_info m_particle_info;
  util::data_buffer* m_data_buffer;

  double m_track_len_in_sv;
  double m_track_len_in_MM;
  double m_track_edep_in_sv;
  double m_track_edep_in_MM;

  double m_kinE_start;

  pos_t m_vertex_position;

public:
  inline void set_event_action(event_action* v) {m_event=v;}

public:
  void digitalize();

  void add_step_info(G4ThreeVector const& pos, double dedx, double len){
    m_track_info.m_step_vertex_px.emplace_back(pos.x());
    m_track_info.m_step_vertex_py.emplace_back(pos.y());
    m_track_info.m_step_vertex_pz.emplace_back(pos.z());
    m_track_info.m_step_dedx.emplace_back(dedx);
    m_track_info.m_step_track_len.emplace_back(len);
  }

  double track_len_in_sv() {return m_track_len_in_sv;}
  void add_track_len_in_sv(double v) {m_track_len_in_sv += v;}
  double track_len_in_MM() {return m_track_len_in_MM;}
  void add_track_len_in_MM(double v) {m_track_len_in_MM += v;}
  double track_edep_in_sv() {return m_track_edep_in_sv;}
  void add_track_edep_in_sv(double v) {m_track_edep_in_sv += v;}
  double track_edep_in_MM() {return m_track_edep_in_MM;}
  void add_track_edep_in_MM(double v) {m_track_edep_in_MM += v;}

  void is_digitalize(bool v) {_is_digitalize = v;}

  void reset();
  void clear_digi();
  



}; 
#endif
