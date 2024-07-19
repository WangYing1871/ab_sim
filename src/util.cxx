#include <cmath>
#include "util.h"

#include "TFile.h"
#include "TTree.h"
#include "TPolyLine3D.h"
#include "TGraph.h"

namespace util::gemo{
double cot(G4ThreeVector const& p0, G4ThreeVector const& p1){
  double x = std::sqrt(std::pow(p0.x()-p1.x(),2)+std::pow(p0.y()-p1.y(),2));
  return x==0. ? std::nan("") : (p0.z()-p1.z())/x;
  
}

}

namespace util{
data_viewer* data_viewer::m_s_data_viewer = nullptr;

data_viewer* data_viewer::fill_h1(std::string const& name, double v){
  m_h1s[name].emplace_back(v); return this;}
data_viewer* data_viewer::fill_h2(std::string const& name, double v0,double v1){
  m_h2s[name].emplace_back(std::array<double,2>{v0,v1}); return this;}
data_viewer* data_viewer::fill_h3(std::string const& name, double v0,double v1,double v2){
  m_h3s[name].emplace_back(std::array<double,3>{v0,v1,v2}); return this;}

}
namespace util::root{
data_buffer::data_buffer(){
  book_1();
}

void data_buffer::book_1(){
  
}

void data_buffer::fill(){
  for (auto&& [x,y] : m_trees) if (y) y->Fill();
}

void data_buffer::open(){
  if (m_file_name.size()==0) m_file_name = "ab_sim.root";
  m_root_file = new TFile(m_file_name.c_str(),"recreate");
  m_root_file->cd();
}
void data_buffer::book(){
  TTree* data_tree = new TTree("data","data");
  
  //TTree* data_tree 
  //data_tree->Branch("a",&m_data_playload.a);
  //data_tree->Branch("b",&m_data_playload.b);
  //data_tree->Branch("c",&m_data_playload.c);
  //m_trees.emplace("data_tree",data_tree);


  //std::string datafile, trackfile, rawfile;
  //if (auto iter_pos = m_file_name.rfind("/"); iter_pos != std::string::npos){
  //  //datafile = m_file_name.substr(0,iter_pos)+
  //}
  
}

void data_buffer::save(){
  for (auto&& [x,y] : m_trees) if (y) y->Write();
  m_root_file->Write(); m_root_file->Close();



}
void data_buffer::clear(){
  m_particles_info.reset();
}

void data_buffer::fill_track_graph(track_info* trk_info
  ,int events, int tracks, int counts){
  if (!_is_all_files_output) return;

}

void data_buffer::save_rawdata(int waveform_x[Tch][Nsp]
    ,int waveform_y[Tch][Nsp]
    ,G4ThreeVector position
    ,double energy){
}

}
