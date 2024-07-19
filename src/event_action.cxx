#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

#include "TFile.h"
#include "TTree.h"

#include "run_action.hh"
#include "event_action.hh"


event_action::event_action(){
}


void event_action::update(util::particle_info_data const& v){
  m_run->data_buf()->m_particles_info.n_track++;
  //m_run->data_buf()->m_particles_info.append(v);
}

void event_action::BeginOfEventAction(G4Event const* evt){
  m_decay_chain = " ";
  m_track_edep_sv = 0.;
  _is_count = true;
  m_tracks = 0;
  m_track_parent_id.clear();
  m_start_position = evt->GetPrimaryVertex()->GetPosition();

  m_run->data_buf()->clear();
}

void event_action::EndOfEventAction(G4Event const* ev){
  //int env_id = ev->GetEventID();
  //...
  //if (m_run->data_buf()->m_particles_info.n_track>0)
  //  m_run->data_buf()->m_ntuple->Fill();




  //...
  auto* trajectories = ev->GetTrajectoryContainer();
  if (trajectories){
    for (size_t i=0; i<m_track_parent_id.size(); ++i){
      for (std::size_t j=0; j<trajectories->entries(); ++j){
        if (m_track_parent_id[i]==trajectories->operator[](j)->GetTrackID()){
          auto pp_name = trajectories->operator[](j)->GetParticleName();
          if (pp_name=="gamma"){

          }
          break;
        }
      }
    }
  }

  info_out(m_run->data_buf()->m_particles_info.n_track);

  m_run->data_buf()->fill();


}
