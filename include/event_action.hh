#ifndef event_action_HH
#define event_action_HH 1 
#include <string>
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "util.h"
class run_action;

class event_action : public G4UserEventAction{
  typedef event_action self_t;
  typedef G4UserEventAction base_t;
  typedef util::particle_info particle_info_t;

public:
  event_action();
  ~event_action() = default;

public:
  virtual void BeginOfEventAction(G4Event const*) override;
  virtual void EndOfEventAction(G4Event const*) override;
  
public:
  void print_modulo(int v) {m_print_modulo = v;}
  void add_decay_chain(std::string const& v) {m_decay_chain += v;}


  void add_edep_sv(double v) {m_track_edep_sv += v;}

  void count(bool v) {_is_count = v;}
  bool count() const {return _is_count;}
  void time(double t) {m_time = t;}
  double time() const {return m_time;}

  void add_sv_track() {m_tracks++;}
  int tracks() const {return m_tracks;}
  void add_parent_track_id(int v) {m_track_parent_id.emplace_back(v);}

private:
  bool _is_count;

  int m_print_modulo;
  int m_tracks;
  double m_edep_sv;
  double m_time;
  double m_track_edep_sv;
  std::string m_decay_chain;
  G4ThreeVector m_start_position;

  std::vector<int> m_track_parent_id;
  run_action* m_run;

public:
  void set_run_action(run_action* v) {m_run = v;}
  run_action* get_run_action() const {return m_run;}
  void update(util::particle_info_data const&);
};
#endif
