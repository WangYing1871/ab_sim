#include <iomanip>

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "run_action.hh"
#include "primary_generator_action.hh"

void run_action::balance(double ekin, double pbal){
  m_decay_count++;
}

void run_action::event_timing(double v){
  m_event_count++;
  m_event_time[0] += v;
  if (m_event_count==1) m_event_time[1]=m_event_time[2]=v;
  if (v<m_event_time[1]) m_event_time[1] = v;
  if (v>m_event_time[2]) m_event_time[2] = v;
}

run_action::run_action(){
}

/* override */
void run_action::BeginOfRunAction(G4Run const*){
  m_decay_count = m_time_count = 0;
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  m_data_buf->open();
  m_data_buf->book();
  
}

/* override */
void run_action::EndOfRunAction(G4Run const* run){
  auto evt_no = run->GetNumberOfEvent();
  if (evt_no==0) return;

  for(auto&& [x,y] : m_particles_recorder){
    //...

  }
  m_particles_recorder.clear();
  m_data_buf->save();
  

  
}
