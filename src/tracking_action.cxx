#include <set>

#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "run_action.hh"
#include "event_action.hh"
#include "tracking_action.hh"
#include "tracking_messenger.hh"

void tracking_action::clear_digi(){
  for(auto&& x : m_charge_x) x.clear();
  for(auto&& x : m_charge_y) x.clear();
  for(auto&& x : m_time_x) x.clear();
  for(auto&& x : m_time_y) x.clear();
  for(auto&& x : m_waveform_x) for(auto&& y : x)  y = 0.;
  for(auto&& x : m_waveform_y) for(auto&& y : x)  y = 0.;
  _is_empty = true;
}

void tracking_action::reset(){
  m_track_len_in_sv = m_track_len_in_MM
  = m_track_edep_in_sv = m_track_edep_in_MM
  = m_max_dedx = m_max_dedx_posz = 0.;
  for (auto&& x : max_position) x = -99.*cm;
  for (auto&& x : min_position) x = 99.*cm;
  m_track_info.reset();


}


tracking_action::tracking_action():
  G4UserTrackingAction(){
  m_message = new tracking_messenger(this);

  
}

tracking_action::~tracking_action(){
  delete m_message;
}


/* override */
void tracking_action::PreUserTrackingAction(G4Track const* track){
  auto* run_action = m_event->get_run_action();
  _is_not_filter = true;
  _is_select_track = true;
  _is_hit_sv = false;

  auto particle = track->GetDefinition();
  auto trk_vertex_pos = track->GetVertexPosition();
  auto trk_vertex_dir = track->GetVertexMomentumDirection();
  auto trk_id = track->GetTrackID();
  auto trk_parent_id = track->GetParentID();

  //create_process_name = trk_id==1 ? "primary" : track->GetCreatorProcess()->GetProcessName();
  reset();
  if (!_is_empty) clear_digi();
  //if (trk_id==1) m_particle_info.reset();

  std::string pname = particle->GetParticleName();
  m_kinE_start = track->GetKineticEnergy();
  run_action->record_particle(pname,m_kinE_start);
  double charge = particle->GetPDGCharge();

  m_vertex_position = track->GetVertexPosition();

  if (charge>2.){
    m_event->add_decay_chain(trk_id==1 ? pname : ("--->"+pname));
    auto* tk_ptr = const_cast<G4Track*>(track);
    tk_ptr->SetTrackStatus(trk_id>1 ? fStopAndKill : fStopButAlive); }

   auto* data_record = util::data_viewer::instance();
   if (particle==G4Electron::Electron() || particle==G4Positron::Positron()) 
     data_record->fill_h1("e+/e-",m_kinE_start);
   if (particle==G4NeutrinoE::NeutrinoE() || G4AntiNeutrinoE::AntiNeutrinoE())
     data_record->fill_h1("neutrino/anti-neutrino(e)",m_kinE_start);
   if (particle==G4Gamma::Gamma()) data_record->fill_h1("gamma",m_kinE_start);
   if (particle==G4Alpha::Alpha()) data_record->fill_h1("alpha",m_kinE_start);
   if (charge>2.) data_record->fill_h1("ion ?",m_kinE_start);
   if (trk_parent_id==2 && particle==G4Electron::Electron())
     data_record->fill_h1("primary beta",m_kinE_start);

  if (m_kinE_start<10*keV || pname=="alpha") _is_select_track = false;
}

/* override */
void tracking_action::PostUserTrackingAction(G4Track const* track){
  auto* run_action = m_event->get_run_action();
  auto* data_record = util::data_viewer::instance();
  auto glb_time = track->GetGlobalTime();
  auto trk_id = track->GetTrackID();
  if (trk_id==1) run_action->primary_timing(glb_time);

  auto secondaries = track->GetStep()->GetSecondaryInCurrentStep();
  size_t sectrk_no = (*secondaries).size();
  if(sectrk_no){
    double ekin_tot = 0.; G4ThreeVector p_banlance = -track->GetMomentum();
    for (size_t i=0; i<sectrk_no; ++i){
      G4Track const* trk = (*secondaries)[i];
      ekin_tot += trk->GetKineticEnergy();
      if (trk->GetDefinition() != G4Gamma::Gamma()) p_banlance += trk->GetMomentum();
    }
    double p = p_banlance.mag();
    run_action->balance(ekin_tot,p);
    //...
  }else{
    run_action->event_timing(glb_time);
    //...
  }

  //...
  bool reject = !m_event->count() && std::abs(glb_time-m_event->time())<500;
  auto vol_name = track->GetVolume()->GetName();
  if (vol_name=="MM_gas_pv") _is_not_filter = false;
  if (_is_select_track && (_is_not_filter || G4UniformRand()<0) 
      && _is_hit_sv && !reject){
    m_event->add_sv_track();
    int evt_trks = m_event->tracks();
    int evt_id = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    //...
    double x_pos = m_track_info.m_step_vertex_px[0];
    double y_pos = m_track_info.m_step_vertex_py[0];

    if (max_position[2]-min_position[2]>0.){
      util::particle_info_data buf;
      buf.x_pos = x_pos; buf.y_pos = y_pos;
      buf.tl_in_sv = m_track_len_in_sv; buf.tl_in_MM = m_track_len_in_MM;
      buf.edep_in_sv = m_track_edep_in_sv; buf.edep_in_MM = m_track_edep_in_MM;
      buf.tl_x = max_position[0]-min_position[0];
      buf.tl_y = max_position[1]-min_position[1];
      buf.tl_z = max_position[2]-min_position[2];
      buf.x = m_vertex_position.x();
      buf.y = m_vertex_position.y();
      buf.z = m_vertex_position.z();
      G4ThreeVector p0(max_position[0],max_position[1],max_position[2]);
      G4ThreeVector p1(min_position[0],min_position[1],min_position[2]);
      buf.cot = util::cot(p0,p1);
      buf.max_dedx = m_max_dedx; buf.max_dedx_z = m_max_dedx_posz;
      buf.steps = m_track_info;
      if ((max_position[2]-m_max_dedx_posz)/(max_position[2]-min_position[2])>=0){
        m_event->update(buf);
        if (m_count<100){
          //m_data_buffer->fill_track_graph(&m_track_info,evt_id,evt_trks,m_count); TODO
          m_count++;
        }
        if (_is_digitalize) digitalize();
      }
      //...

      m_event->add_parent_track_id(track->GetParentID());
    }
  }else if(track->GetDefinition()->GetParticleName()=="e-"){
    if(track->GetGlobalTime()>m_event->time()) m_event->time(track->GetGlobalTime()); }
}
//---------------------------------------------------------------------
#include "constant.h"
void tracking_action::digitalize(){
  typedef std::vector<double>(aarray_t)[TPCsystem::nch];
  auto const& find_min = [](aarray_t const& data)->double{
    std::multiset<double> min_set;
    for (auto&& x : data) min_set.insert(*std::min_element(x.begin(),x.end()));
    return *min_set.begin(); };
  double min_tx = find_min(m_time_x);
  double min_ty = find_min(m_time_y);
  double min_t = std::min(min_tx,min_ty);

  using TPCsystem::timebinwidth;
  using TPCsystem::gain;
  using TPCsystem::factor;
  using TPCsystem::CSAgain;
  using TPCsystem::response_func;
  const int response_func_sz = 
    (int)(sizeof(TPCsystem::response_func)/sizeof(TPCsystem::response_func[0]));

  typedef decltype(m_time_x) time_tt;
  typedef decltype(m_charge_x) charge_t;
  typedef decltype(m_waveform_x) waveform_t;

  auto const& conv = [&](time_tt const& tms, charge_t const& charges, waveform_t& wave){
    for (int i=0; i<TPCsystem::nch; ++i){
      bool is_eff = false;
      for (size_t j=0; j<charges[i].size(); ++j){
        int delta_t = (int)((tms[i][j]-min_t)/40);
        if (delta_t>=(TPCsystem::Nsp-100)) continue;
        int tm_binstart = -delta_t-100+TPCsystem::peak_pos;
        for (int k=0; k<TPCsystem::Nsp; ++k){
          if (tm_binstart>0 && tm_binstart<response_func_sz && charges[i][j]!=0){
            wave[i][k] += 1.60218e-19/timebinwidth*gain*factor/CSAgain*charges[i][j]*
              response_func[tm_binstart];
            is_eff = true; }
          tm_binstart++; }
      }
      if (is_eff && wave[i][0]==0) wave[i][0] = -1; }
  };
  conv(m_time_x,m_charge_x,m_waveform_x);
  conv(m_time_y,m_charge_y,m_waveform_y);
  //...
}
