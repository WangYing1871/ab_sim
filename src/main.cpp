//--------------------------------Stamp-------------------------------
//^-^ Author: Zhi Heng            Email: 2835516101@qq.com     
//^-^ Time: 2024-07-12 15:03:29   Posi: Hefei
//^-^ File: main.cpp
//--------------------------------------------------------------------
#include <algorithm>
#include <numeric>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <ctime>

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

#include "detector_construction.hh"
#include "primary_generator_action.hh"
#include "run_action.hh"
#include "event_action.hh"
#include "tracking_action.hh"
#include "step_action.hh"
#include "physics_list.hh"

#include "util.h"


int main(int argc, char* argv[]){
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheSeed(std::time(0));

  auto* run_manager = new G4RunManager;
  auto* detector = new detector_construction;
  run_manager->SetUserInitialization(detector);
  auto* physics = new physics_list;
  run_manager->SetUserInitialization(physics);

  auto* primary_generator = new primary_generator_action(detector);
  //primary_generator->SetVerboseLevel(1);
  run_manager->SetUserAction(primary_generator);

  auto* run = new run_action();
  auto* data_buffer = new util::data_buffer();
  std::string out_fname = "ab_sim.root";
  data_buffer->file_name(out_fname);
  run->set_data_buf(data_buffer);

  run->set_primary_generator(primary_generator);
  run_manager->SetUserAction(run);

  auto* event = new event_action();
  run_manager->SetUserAction(event);
  event->set_run_action(run);

  auto* track = new tracking_action();
  //track->is_digitalize(true);
  track->set_event_action(event);
  run_manager->SetUserAction(track);

  auto* step = new step_action();
  step->set_detector(detector);
  step->set_tracking_action(track);
  step->set_event_action(event);
  run_manager->SetUserAction(step);

  run_manager->Initialize();
  auto* vis_manager = new G4VisExecutive(argc,argv);
  vis_manager->Initialize();
  auto* ui = G4UImanager::GetUIpointer();
  
  if (argc!=1){
    G4String mac_name = argv[1];
    ui->ApplyCommand("/control/execute "+mac_name);
  }else{
    auto* ui_exec = new G4UIExecutive(argc,argv);
    ui->ApplyCommand("/control/execute vis.mac");
    ui_exec->SessionStart();
    delete ui_exec;
  }

  delete vis_manager;
  delete run_manager;
  return 0;
}
