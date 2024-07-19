#include <random>

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Step.hh"

#include "step_action.hh"
#include "detector_construction.hh"
#include "tracking_action.hh"
#include "event_action.hh"

step_action* step_action::m_s_step_action = nullptr;
step_action* step_action::instance(){ return m_s_step_action; }

step_action::step_action():
  G4UserSteppingAction(){
  m_s_step_action = this;
}

step_action::~step_action(){ m_s_step_action = nullptr; }

/* override */
void step_action::UserSteppingAction(G4Step const* step){
  G4String pre_vm = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  double edep = step->GetTotalEnergyDeposit()-step->GetNonIonizingEnergyDeposit();
  double step_len = step->GetStepLength();
  if (pre_vm=="gas_pv" || pre_vm=="gas_eff_pv") m_event_action->add_edep_sv(edep);
  if (m_track_action->is_select_track() && m_track_action->is_not_filtered()){
    if (pre_vm=="MM_gas_pv") m_track_action->not_filtered(false);
    if (pre_vm=="gas_pv" || pre_vm=="gas_eff_pv"){
      m_track_action->hit_sv(true);
      auto vertex_pos = step->GetPreStepPoint()->GetPosition();
      double ene_transfer = step->GetPreStepPoint()->GetKineticEnergy()-
        step->GetPostStepPoint()->GetKineticEnergy();
      m_track_action->add_track_len_in_sv(step_len);
      m_track_action->add_track_edep_in_sv(edep);
      m_edep += ene_transfer; m_step_len += step_len;
      if (m_step_len>0.5*mm || step->GetPostStepPoint()->GetKineticEnergy()==0){
        double dEdx = m_edep/m_step_len;
        if (vertex_pos.x()>m_track_action->max_position[0]) 
          m_track_action->max_position[0] = vertex_pos.x();
        if (vertex_pos.y()>m_track_action->max_position[1]) 
          m_track_action->max_position[1] = vertex_pos.y();
        if (vertex_pos.z()>m_track_action->max_position[2]) 
          m_track_action->max_position[2] = vertex_pos.z();
        if (vertex_pos.x()<m_track_action->min_position[0]) 
          m_track_action->min_position[0] = vertex_pos.x();
        if (vertex_pos.y()<m_track_action->min_position[1]) 
          m_track_action->min_position[1] = vertex_pos.y();
        if (vertex_pos.z()<m_track_action->min_position[2]) 
          m_track_action->min_position[2] = vertex_pos.z();
        if (dEdx > m_track_action->m_max_dedx){
          m_track_action->m_max_dedx = dEdx;
          //...
          m_track_action->m_max_dedx_posz = vertex_pos.z();
        }
        m_track_action->add_step_info(vertex_pos,dEdx,m_track_action->track_len_in_sv());
        m_edep = 0.; m_step_len = 0.;
      }
      if (m_track_action->_is_digitalize) drift_one_electron(vertex_pos,ene_transfer);
    }else if (pre_vm=="MM_gas_pv"){
      double ene_transfer = step->GetPreStepPoint()->GetKineticEnergy()-
        step->GetPostStepPoint()->GetKineticEnergy();
      m_track_action->add_track_len_in_MM(step_len);
      m_track_action->add_track_edep_in_MM(ene_transfer);
    }
  }
}

void step_action::reset(){
}

void step_action::drift_one_electron(G4ThreeVector const& pos, double edep){
  double drift_distance = -pos.z();
  if (drift_distance<0) return;
  double sigma_l = TPCsystem::dl*std::sqrt(drift_distance*0.1)*10;
  double sigma_t = TPCsystem::dt*std::sqrt(drift_distance*0.1)*10;
  if (drift_distance>70.){
  }
  std::default_random_engine generator;
  std::normal_distribution<double> distributionL(0.0,sigma_l);
  std::normal_distribution<double> distributionT(0.0,sigma_t);
  
  drift_distance += distributionL(generator);

  double projected_x = pos.x()+distributionT(generator);
  double projected_y = pos.y()+distributionT(generator);
  double x_rot = std::sqrt(2)/2*(projected_x-projected_y);
  double y_rot = std::sqrt(2)/2*(projected_x+projected_y);

  //FIXME
  double sqrt2 = std::sqrt(2);
  int x_pos = (int)((x_rot+.5*sqrt2/2*TPCsystem::chnwidth)/(sqrt2/2*TPCsystem::chnwidth));
  int y_pos = (int)((y_rot+.5*sqrt2/2*TPCsystem::chnwidth)/(sqrt2/2*TPCsystem::chnwidth));
  if ((x_pos+0.5*sqrt2/2*TPCsystem::chnwidth)<0) x_pos -= 1;
  if ((y_pos+0.5*sqrt2/2*TPCsystem::chnwidth)<0) x_pos -= 1;
  
  if (!TPCsystem::on_plane_spread){
    fill_one_channel_signal(x_pos,y_pos,edep,drift_distance);
  }else{
    double x_center = x_pos*TPCsystem::chnwidth/sqrt2;
    double y_center = y_pos*TPCsystem::chnwidth/sqrt2;
    if ((x_pos+y_pos)%2==0){
      double dx_plus = 0.5*TPCsystem::chnwidth-(x_rot+y_rot-x_center-y_center)/sqrt2;
      double dx_minus = 0.5*TPCsystem::chnwidth+(x_rot+y_rot-x_center-y_center)/sqrt2;
      double dy_plus = TPCsystem::chnwidth+(x_rot-y_rot-x_center+y_center)/sqrt2;
      double dy_minus = TPCsystem::chnwidth-(x_rot-y_rot-x_center+y_center)/sqrt2;
      using namespace std;
      using TPCsystem::sigma_spread;
      double r_dx_plus = exp(-(pow(dx_plus,2)/(2*pow(sigma_spread,2))));
      double r_dx_minus = exp(-(pow(dx_minus,2)/(2*pow(sigma_spread,2))));
      double r_dy_plus = exp(-(pow(dy_plus,2)/(2*pow(sigma_spread,2))));
      double r_dy_minus = exp(-(pow(dy_minus,2)/(2*pow(sigma_spread,2))));
      double r_this = 1.;
      double r_sum = r_this+r_dx_plus+r_dx_minus+r_dy_plus+r_dy_minus;
      fill_one_channel_signal(x_pos,y_pos,edep*r_this/r_sum,drift_distance);
      fill_one_channel_signal(x_pos+1,y_pos,edep*r_dx_plus/r_sum,drift_distance);
      fill_one_channel_signal(x_pos-1,y_pos,edep*r_dx_minus/r_sum,drift_distance);
      fill_one_channel_signal(x_pos-1,y_pos+1,edep*r_dx_plus/r_sum,drift_distance);
      fill_one_channel_signal(x_pos+1,y_pos-1,edep*r_dx_minus/r_sum,drift_distance);
    }else{
      double dx_plus = TPCsystem::chnwidth-(x_rot+y_rot-x_center-y_center)/sqrt2;
      double dx_minus = TPCsystem::chnwidth+(x_rot+y_rot-x_center-y_center)/sqrt2;
      double dy_plus = 0.5*TPCsystem::chnwidth+(x_rot-y_rot-x_center+y_center)/sqrt2;
      double dy_minus = 0.5*TPCsystem::chnwidth-(x_rot-y_rot-x_center+y_center)/sqrt2;
      using namespace std;
      using TPCsystem::sigma_spread;
      double r_dx_plus = exp(-(pow(dx_plus,2)/(2*pow(sigma_spread,2))));
      double r_dx_minus = exp(-(pow(dx_minus,2)/(2*pow(sigma_spread,2))));
      double r_dy_plus = exp(-(pow(dy_plus,2)/(2*pow(sigma_spread,2))));
      double r_dy_minus = exp(-(pow(dy_minus,2)/(2*pow(sigma_spread,2))));
      double r_this = 1.;
      double r_sum = r_this+r_dx_plus+r_dx_minus+r_dy_plus+r_dy_minus;
      fill_one_channel_signal(x_pos,y_pos,edep*r_this/r_sum,drift_distance);
      fill_one_channel_signal(x_pos+1,y_pos+1,edep*r_dx_plus/r_sum,drift_distance);
      fill_one_channel_signal(x_pos-1,y_pos-1,edep*r_dx_minus/r_sum,drift_distance);
      fill_one_channel_signal(x_pos,y_pos+1,edep*r_dx_plus/r_sum,drift_distance);
      fill_one_channel_signal(x_pos,y_pos-1,edep*r_dx_minus/r_sum,drift_distance);
    }
  }
}

void step_action::fill_one_channel_signal(int x_pos, int y_pos, double edep, double drift_distance){
  using TPCsystem::nch;
  double drift_speed = TPCsystem::v_drift;
  if ((x_pos+y_pos)%2==0){
    int chn = (y_pos-x_pos)/2 + nch/2;
    if (chn>=0 && chn<nch){
      m_track_action->m_charge_y[chn].emplace_back(edep*1e6/TPCsystem::E_ion);
      m_track_action->m_time_y[chn].emplace_back(drift_distance/drift_speed*100);
      m_track_action->_is_empty = false; }
  }else{
    int chn = (y_pos+x_pos-1)/2 + nch/2;
    if (chn>=0 && chn<nch){
      m_track_action->m_charge_x[chn].emplace_back(edep*1e6/TPCsystem::E_ion);
      m_track_action->m_time_x[chn].emplace_back(drift_distance/drift_speed*100);
      m_track_action->_is_empty = false; }
  }
}

