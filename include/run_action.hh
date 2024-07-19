#ifndef run_action_HH
#define run_action_HH 1 

#include <unordered_map>
#include <unordered_set>
#include <array>

#include "G4UserRunAction.hh"
#include "globals.hh"

#include "util.h"

class G4Run;
//class primary_generator_action;
class G4VUserPrimaryGeneratorAction;


class run_action : public G4UserRunAction{
  typedef run_action self_t;
  typedef G4UserRunAction base_t;
  typedef std::array<double,3> sum_min_max_t;

public:
  run_action();
  ~run_action()=default;

  inline void primary_timing(double v) {m_primary_time+=v;}
  void event_timing(double v);

public:
  virtual void BeginOfRunAction(G4Run const*) override;
  virtual void EndOfRunAction(G4Run const*) override;

  inline void record_particle(std::string const& name, double ekin){
    m_particles_recorder[name].insert(ekin); }
  void balance(double ekin, double pbal);

  void set_primary_generator(G4VUserPrimaryGeneratorAction* v) {m_primary = v;}
  void set_data_buf(util::data_buffer* v) {m_data_buf = v;}
  util::data_buffer* data_buf() const {return m_data_buf;}

private:
  //primary_generator_action* m_primary;
  G4VUserPrimaryGeneratorAction* m_primary;
  util::data_buffer* m_data_buf;

  int m_decay_count, m_time_count; int m_event_count;
  double m_primary_time;
  std::unordered_map<std::string,std::multiset<double>> m_particles_recorder;
  sum_min_max_t m_ekin_tot, m_pbalance, m_event_time;
};

#endif
