#ifndef util_H
#define util_H 1 

#include <iostream>
#define info_out(X) std::cout<<"==> "<<__LINE__<<" "<<#X<<" |"<<(X)<<"|\n"
#include <vector>
#include <string>
#include <cstring>
#include <array>
#include <functional>
#include <type_traits>
#include <unordered_map>

#include "G4ThreeVector.hh"

#include "constant.h"
class TFile;
class TTree;
namespace util{

namespace gemo{
double cot(G4ThreeVector const&, G4ThreeVector const&);


}

using namespace gemo;

namespace config{
static char const* const HEPEvt_Interface_filename = "pythia_event.data";
}

struct data_viewer{
  typedef data_viewer self_t;

public:
  static self_t* instance(){
    if (m_s_data_viewer) return m_s_data_viewer;
    static data_viewer tmp;
    m_s_data_viewer = std::addressof(tmp);
    return m_s_data_viewer; }

  self_t* fill_h1(std::string const&, double);
  self_t* fill_h2(std::string const&, double,double);
  self_t* fill_h3(std::string const&, double,double,double);

  template <size_t v>
  using array_t = typename std::conditional<v==1,std::vector<double>
     ,typename std::conditional<v==2,std::vector<std::array<double,2>>,
        std::vector<std::array<double,3>>>::type>::type;
  template <class _tp,size_t _v=1>
  _tp get(std::string const& name, std::function<_tp(array_t<_v>&)> f){
    if constexpr(_v==1){
      if (auto iter = m_h1s.find(name); iter != m_h1s.end())
        return std::invoke(f,*iter);
      throw std::invalid_argument("... 1"); }
    if constexpr(_v==2){
      if (auto iter = m_h2s.find(name); iter != m_h2s.end())
        return std::invoke(f,*iter);
      throw std::invalid_argument("... 2"); }
    if constexpr(_v==3){
      if (auto iter = m_h3s.find(name); iter != m_h3s.end())
        return std::invoke(f,*iter);
      throw std::invalid_argument("... 3"); }
  }


private:
  static data_viewer* m_s_data_viewer;
  data_viewer() = default;
  ~data_viewer() noexcept = default;

  std::unordered_map<std::string,std::vector<double>> m_h1s;
  std::unordered_map<std::string,std::vector<std::array<double,2>>> m_h2s;
  std::unordered_map<std::string,std::vector<std::array<double,3>>> m_h3s;
};

class track_info{
  typedef std::vector<double> dyn_array_t;
public:
  dyn_array_t m_step_vertex_px;
  dyn_array_t m_step_vertex_py;
  dyn_array_t m_step_vertex_pz;
  dyn_array_t m_step_dedx;
  dyn_array_t m_step_track_len;

  void reset(){
    m_step_vertex_px.clear();
    m_step_vertex_py.clear();
    m_step_vertex_pz.clear();
    m_step_dedx.clear();
    m_step_track_len.clear(); }
  track_info() = default;
};
struct particle_info_data{
  double x_pos;
  double y_pos;
  double tl_in_sv, tl_in_MM, edep_in_sv, edep_in_MM;
  double tl_x, tl_y, tl_z;
  double max_dedx, max_dedx_z;
  double cot;
  double x,y,z;
  track_info steps;
};
class particle_info{
public:
  int n_track=0;
  typedef std::vector<double> dyn_array_t;
  typedef std::vector<std::vector<double>> dyn_aarray_t;
  dyn_array_t m_particles_Z;
  dyn_array_t m_particles_A;
  dyn_array_t m_particles_mass;
  dyn_array_t m_particles_code;

  dyn_array_t m_track_kinetic_energy;
  dyn_array_t m_track_vertex_x;
  dyn_array_t m_track_vertex_y;
  dyn_array_t m_track_vertex_z;

  dyn_array_t m_track_id;
  dyn_array_t m_parent_id;

  dyn_array_t m_track_time;

  dyn_array_t m_start_px;
  dyn_array_t m_start_py;
  dyn_array_t m_track_length, m_track_length_MM;
  dyn_array_t m_drift_distance;
  dyn_array_t m_x_track_length;
  dyn_array_t m_y_track_length;
  dyn_array_t m_cot_theta;
  dyn_array_t m_track_energy, m_track_energy_MM;
  dyn_array_t m_maxedep_pos;
  dyn_array_t m_maxedep_posZ;

  dyn_aarray_t m_track_start_pos;
  dyn_aarray_t m_event_start_pos;
  dyn_aarray_t m_track_vertexposx_perstep;
  dyn_aarray_t m_track_vertexposy_perstep;
  dyn_aarray_t m_track_vertexposz_perstep;
  dyn_aarray_t m_track_Ededx_perstep;
  dyn_aarray_t m_track_length_perstep;

  void append(particle_info_data const& v){
    //m_start_px.emplace_back(v.x_pos);
    //m_start_py.emplace_back(v.y_pos);
    //m_track_length.emplace_back(v.tl_in_sv);
    //m_track_length_MM.emplace_back(v.tl_in_MM);
    //m_track_energy.emplace_back(v.edep_in_sv);
    //m_track_energy_MM.emplace_back(v.edep_in_MM);
    //m_cot_theta.emplace_back(v.cot);
    //m_maxedep_pos.emplace_back(v.max_dedx);
    //m_maxedep_posZ.emplace_back(v.max_dedx_z);
    //m_drift_distance.emplace_back(v.tl_z);
    //m_x_track_length.emplace_back(v.tl_x);
    //m_y_track_length.emplace_back(v.tl_y);

    //m_track_start_pos.emplace_back(v.x);
    //m_track_start_pos.emplace_back(v.y);
    //m_track_start_pos.emplace_back(v.z);

    //m_track_vertexposx_perstep.emplace_back(v.steps.m_step_vertex_px);
    //m_track_vertexposy_perstep.emplace_back(v.steps.m_step_vertex_py);
    //m_track_vertexposz_perstep.emplace_back(v.steps.m_step_vertex_pz);
    //m_track_Ededx_perstep.emplace_back(v.steps.m_step_dedx);
    //m_track_length_perstep.emplace_back(v.steps.m_step_track_len);
  }

  void reset(){
    n_track = 0;
    m_particles_Z.clear();
    m_particles_A.clear();
    m_particles_mass.clear();
    m_particles_code.clear();
    m_track_kinetic_energy.clear();
    m_track_vertex_x.clear();
    m_track_vertex_y.clear();
    m_track_vertex_z.clear();
    m_track_id.clear();
    m_parent_id.clear();
    m_track_time.clear();
    m_start_px.clear();
    m_start_py.clear();
    m_track_length.clear(), m_track_length_MM.clear();
    m_drift_distance.clear();
    m_x_track_length.clear();
    m_y_track_length.clear();
    m_cot_theta.clear();
    m_track_energy.clear(), m_track_energy_MM.clear();
    m_maxedep_pos.clear();
    m_maxedep_posZ.clear();

    m_track_start_pos.clear();
    m_event_start_pos.clear();
    m_track_vertexposx_perstep.clear();
    m_track_vertexposy_perstep.clear();
    m_track_vertexposz_perstep.clear();
    m_track_Ededx_perstep.clear();
    m_track_length_perstep.clear(); }

  particle_info()=default;

  friend class particle_info_data;
};


namespace root{
using namespace TPCsystem;
struct raw_data{

  int event_id=0;
  int n_hits=0;
  int fec[Tch];
  int chip[Tch];
  int chn[Tch];
  int adc[Tch][Nsp];
  float sum_adc[Tch];
  float max_adc[Tch];
  float max_point[Tch];
  float summax_adc;
  int pixelX[Tch];
  int pixelY[Tch];
  float startpos_x;
  float startpos_y;
  float startpos_z;
  float kinE_start;

  void reset(){ memset(this,0,sizeof(raw_data)); }
  raw_data() = default;
};

class data_buffer{
public:
  data_buffer();
  ~data_buffer();

  void save();
  void book();
  void open();
  void fill();
  void clear();

  void fill_track_graph(track_info*,int,int,int);
  void save_rawdata(int waveform_X[Tch][Nsp], int waveform_Y[Tch][Nsp]
      ,G4ThreeVector position, double energy);

  void book_1();
  std::string m_file_name;

public:
  TFile* m_root_file,* m_graph_root_file,* m_raw_root_file;
  std::unordered_map<std::string,TTree*> m_trees;
  bool _is_all_files_output = false;
 
public:
  particle_info m_particles_info;
  raw_data raw_data_frame;

public:
  void file_name(std::string const& v) {m_file_name = v;}

public:
  struct data_playload{
    int a;
    double b;
    std::vector<double> c;
  };
  data_playload m_data_playload;


};

}

using namespace root;
}
#endif
