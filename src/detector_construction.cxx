#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4GlobalMagFieldMessenger.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4AutoDelete.hh"

#include "detector_construction.hh"

G4ThreadLocal G4GlobalMagFieldMessenger* detector_construction::m_s_magfield_messenger = 0;


void detector_construction::init_materials(){
  auto* g4nm = G4NistManager::Instance();
  m_materials.emplace("Vacuum",g4nm->FindOrBuildMaterial("G4_Galactic"));
  m_materials.emplace("Al",g4nm->FindOrBuildMaterial("G4_Al"));
  m_materials.emplace("Cu",g4nm->FindOrBuildMaterial("G4_Cu"));
  m_materials.emplace("Pb",g4nm->FindOrBuildMaterial("G4_Pb"));
  m_materials.emplace("Ge",g4nm->FindOrBuildMaterial("G4_Ge"));
  m_materials.emplace("Air",g4nm->FindOrBuildMaterial("G4_AIR"));
  m_materials.emplace("CO2",g4nm->FindOrBuildMaterial("G4_CARBON_DIOXIDE"));
  m_materials.emplace("Ar",g4nm->FindOrBuildMaterial("G4_Ar"));
  m_materials.emplace("Ne",g4nm->FindOrBuildMaterial("G4_Ne"));
 // m_materials.emplace("Xe",g4nm->FindOrBuildMaterial("G4_Xe"));

	auto* C  = new G4Element("Carbon",    "C",  6., 12.00*g/mole);
	auto* H  = new G4Element("Hydrogen",  "H",  1., 1.00*g/mole);
	auto* O  = new G4Element("Oxygen",    "O" , 8. ,16.00*g/mole);
	auto* N  = new G4Element("Nitrogen",  "N" , 7. ,14.00*g/mole);
	auto* Si = new G4Element("Si",        "Si", 14.,28.00*g/mole);
	auto* Fe = new G4Element("Iron",      "Fe", 26.,55.85*g/mole);
	auto* Mn = new G4Element("Manganese", "Mn", 25.,54.94*g/mole);
	auto* Cr = new G4Element("Chromium",  "Cr", 24.,52.00*g/mole);
	auto* Ni = new G4Element("Nickel",    "Ni", 28.,58.7*g/mole);

  auto* PET = new G4Material("PET",1.38*g/cm3,3);
  PET->AddElement(C,10); PET->AddElement(H,8); PET->AddElement(O,4);

  auto* epoxy = new G4Material("epoxy",1.2*g/cm3,2);
  epoxy->AddElement(H,2); epoxy->AddElement(C,2);
  auto* SiO2 = new G4Material("SiO2",2.200*g/cm3,2);
  SiO2->AddElement(Si,1); SiO2->AddElement(O,2);
  auto* FR4 = new G4Material("FR4",1.86*g/cm3,2);
  FR4->AddMaterial(epoxy,47.2*perCent); FR4->AddMaterial(epoxy,52.8*perCent);

  auto* ployacrylate = new G4Material("ployacrylate",1.39*g/cm3,3);
  ployacrylate->AddElement(C,3); ployacrylate->AddElement(H,4); ployacrylate->AddElement(O,2);

  auto* atlasgas = new G4Material("atlasgas",1.79e-3*g/cm3,2);
  atlasgas->AddMaterial(m_materials["CO2"],7.0*perCent); 
  atlasgas->AddMaterial(m_materials["Ar"],93.0*perCent);

  auto* C4H10 = new G4Material("C4H10",2.48*mg/cm3,2);
  C4H10->AddElement(C,4); C4H10->AddElement(H,10);

  auto* Ar_C4H10 = new G4Material("Ar+C4H10",(0.035*2.48+96.5*1.78)*mg/cm3,2);
  Ar_C4H10->AddMaterial(C4H10,3.5*perCent);
  Ar_C4H10->AddMaterial(m_materials["Ar"],96.5*perCent);

  m_materials.emplace("Ar+C4H10",Ar_C4H10);

  auto* kapton = new G4Material("kapton",1.43*g/cm3,4);
  kapton->AddElement(C,22); kapton->AddElement(H,10);
  kapton->AddElement(N,2); kapton->AddElement(O,5);

  auto* wood = new G4Material("wood",0.9*g/cm3,3);
  wood->AddElement(H,4); wood->AddElement(O,1); wood->AddElement(C,2);

	auto* mat304steel = new G4Material("Stainless steel 304",7.999*g/cm3,6);
	mat304steel->AddElement(Mn, 0.02); mat304steel->AddElement(Si, 0.01);
	mat304steel->AddElement(Cr, 0.19); mat304steel->AddElement(Ni, 0.10);
	mat304steel->AddElement(Fe, 0.6792); mat304steel->AddElement(C, 0.0008);

  m_materials.emplace("PET",PET);
  m_materials.emplace("FR4",FR4);
  m_materials.emplace("ployacrylate",ployacrylate);
  m_materials.emplace("atlasgas",atlasgas);
  m_materials.emplace("kapton",kapton);
  m_materials.emplace("wood",wood);
  m_materials.emplace("mat304steel",mat304steel);
}

/* override */
G4VPhysicalVolume* detector_construction::Construct(){
  init_materials();

  auto world_rmax = 100*cm;
  auto* world_sv = new G4Orb("world_sv",world_rmax);
  auto* world_lv = new G4LogicalVolume(world_sv,m_materials["Air"],"world_lv");
  auto* world_pv = new G4PVPlacement(0,G4ThreeVector(),world_lv,"world_pv",0,false,0);

  double frame_sx = 26.9*cm;
  double frame_sy = 26.9*cm;

  double TPC_sz = 6.5*cm;

  auto* TPC_sv = new G4Box("TPC_sv",.5*frame_sx,.5*frame_sy,.5*TPC_sz);
  auto* TPC_lv = new G4LogicalVolume(TPC_sv,m_materials["Ar+C4H10"],"TPC_lv");
  auto TPC_positon = G4ThreeVector(0,0,-.5*TPC_sz);
  new G4PVPlacement(0,TPC_positon,TPC_lv,"TPC_Pv",world_lv,false,0);

  double frame1_height = 18.3*mm;
  double frame1_thickness = 41.*mm;
  double frame2_height = 34.7*mm;
  double frame2_thickness = 12.*mm;
  double frame3_height = 12.*mm;
  double frame3_inner_radius = 58.*mm;

  G4ThreeVector frame1_position(0,0,-.5*frame1_height);
  auto* frame1out_sv = new G4Box("frame1out_sv"
      ,.5*frame_sx,.5*frame_sy,.5*frame1_height);
  auto* frame1in_sv = new G4Box("frame1in_sv"
      ,.5*frame_sx-frame1_thickness,.5*frame_sy-frame1_thickness,.5*frame1_height);
  auto* frame1_sv = new G4SubtractionSolid("frame1_sv",frame1out_sv,frame1in_sv);
  auto* frame1_lv = new G4LogicalVolume(frame1_sv,m_materials["Al"],"frame1_lv");
  new G4PVPlacement(0,frame1_position-TPC_positon,frame1_lv,"frame1_pv",TPC_lv,false,0);

  G4ThreeVector frame2_position(0,0,-frame1_height-0.5*frame2_height);
  auto* frame2out_sv = new G4Box("frame2out_sv"
      ,.5*frame_sx,.5*frame_sy,.5*frame2_height);
  auto* frame2in_sv = new G4Box("frame2in_sv" 
      ,.5*frame_sx-frame2_thickness,.5*frame_sy-frame2_thickness,.5*frame2_height);
  auto* frame2_sv= new G4SubtractionSolid("frame2_sv",frame2out_sv,frame2in_sv);
  auto* frame2_lv = new G4LogicalVolume(frame2_sv,m_materials["Al"],"frame2_lv");
  new G4PVPlacement(0,frame2_position-TPC_positon,frame2_lv ,"frame2_pv",TPC_lv,false,0);

  G4ThreeVector frame3_position(0.,0.,-frame1_height-frame2_height-0.5*frame3_height);
  auto* frame3out_sv = new G4Box("frame3out_sv",.5*frame_sx,.5*frame_sy,.5*frame3_height);
  auto* frame3in_sv = new G4Tubs("frame3in_sv",0,frame3_inner_radius
      ,.5*frame3_height,0.,360.);
  auto* frame3_sv = new G4SubtractionSolid("frame3_sv",frame3out_sv,frame3in_sv);
  auto* frame3_lv = new G4LogicalVolume(frame3_sv,m_materials["Al"],"frame3_lv");
  new G4PVPlacement(0,frame3_position-TPC_positon,frame3_lv,"frame3_pv",TPC_lv,false,0);

  double window_thickness = 5*mm;
  double window_radius = 35*mm;
  G4ThreeVector window_position(0,0,-TPC_sz+.5*window_thickness);
  auto* window_sv = new G4Tubs("window_sv",window_radius
      ,frame3_inner_radius,0.5*window_thickness,0.,360.);
  auto* window_lv = new G4LogicalVolume(window_sv,m_materials["Al"],"window_lv");
  new G4PVPlacement(0,window_position-TPC_positon,window_lv,"window_pv",TPC_lv,false,0);

  double pet_sz = 3.e-3*mm;
  G4ThreeVector film_position(0,0,-TPC_sz+window_thickness+.5*pet_sz);
  auto* film_sv = new G4Tubs("film_sv",0,frame3_inner_radius,.5*pet_sz,0.,360.);
  auto* film_lv = new G4LogicalVolume(film_sv,m_materials["PET"],"film_lv");
  new G4PVPlacement(0,film_position-TPC_positon,film_lv,"film_pv",TPC_lv,false,0);

  double ring_thickness = 3*mm;
  G4ThreeVector ring_position(0,0,-TPC_sz+window_thickness+pet_sz+.5*ring_thickness);
  auto* ring_sv = new G4Tubs("ring_sv",window_radius,frame3_inner_radius
      ,.5*ring_thickness,0,360);
  auto* ring_lv = new G4LogicalVolume(ring_sv,m_materials["Al"],"ring_lv");
  new G4PVPlacement(0,ring_position-TPC_positon,ring_lv,"ring_pv",TPC_lv,false,0);

  G4ThreeVector step_position(0,0,-TPC_sz+0.5*window_thickness);
  auto* step_sv = new G4Tubs("step_sv",0,window_radius,.5*window_thickness,0,360);
  auto* step_lv = new G4LogicalVolume(step_sv,m_materials["Air"],"step_lv");
  new G4PVPlacement(0,step_position-TPC_positon,step_lv,"step_pv",TPC_lv,false,0);

  //TODO 
  //if (){ 
  //}

	double field_cage_size = 16.6*cm;
	double field_cage_thinkness = 0.16*cm;
  double field_cage_height = 5.*cm;
  G4ThreeVector field_cage_position(0,0,-0.5*field_cage_height);
  auto* field_cageout_sv = new G4Box("field_cageout_sv"
      ,.5*field_cage_size+field_cage_thinkness
      ,.5*field_cage_size+field_cage_thinkness
      ,.5*field_cage_height);
  auto* field_cagein_sv = new G4Box("field_cagein_sv"
      ,.5*field_cage_size
      ,.5*field_cage_size
      ,.5*field_cage_height);
  auto* field_cage_sv = new G4SubtractionSolid("field_cage_sv",field_cageout_sv,field_cagein_sv);
  auto* field_cage_lv = new G4LogicalVolume(field_cage_sv,m_materials["FR4"],"field_cage_lv");
  new G4PVPlacement(0,field_cage_position-TPC_positon,field_cage_lv,"field_cage_pv",TPC_lv,false,0);

  double TPC_gas_thickness = TPC_sz-frame3_height;
  double detector_sx = 15*cm;
  double detector_sy = 15*cm;
  G4ThreeVector gas_position(0,0,-.5*TPC_gas_thickness);
  auto* gas_region0_sv = new G4Box("gas_region0_sv",0.5*detector_sx,0.5*detector_sy
      ,0.5*TPC_gas_thickness);
  auto* gas_region1_sv = new G4Tubs("gas_region1_sv",0,frame3_inner_radius
      ,0.5*(frame3_height-window_thickness-pet_sz-ring_thickness)
      ,0.
      ,360.);
  auto* gas_region2_sv = new G4Tubs("gas_region2_sv",0,window_radius,.5*ring_thickness,0.,360.);
  auto* empty_rm = new G4RotationMatrix;
  G4ThreeVector shift1(0,0
      ,-.5*(TPC_gas_thickness+frame3_height-window_thickness-pet_sz-ring_thickness));
  G4ThreeVector shift2 (0,0
      ,-.5*(TPC_gas_thickness+2*frame3_height-2*window_thickness-2*pet_sz-ring_thickness));
  auto* gas_sv = new G4UnionSolid("gas_sv"
      ,new G4UnionSolid("tmp00_sv",gas_region0_sv,gas_region1_sv,empty_rm,shift1)
      ,gas_region2_sv,empty_rm,shift2);
  auto* gas_lv = new G4LogicalVolume(gas_sv,m_materials["Ar+C4H10"],"gas_lv");
  new G4PVPlacement(0,gas_position-TPC_positon,gas_lv,"gas_pv",TPC_lv,false,0);

  double gas_eff_sx=12.*cm;
  double gas_eff_sy=12.*cm;
  G4ThreeVector gas_eff_position = G4ThreeVector(0,0,-.5*TPC_gas_thickness);
  auto* gas_eff_region0_sv = new G4Box("gas_eff_region0_sv",
      .5*gas_eff_sx,.5*gas_eff_sy,.5*TPC_gas_thickness);
  auto* gas_eff_sv = new G4UnionSolid("gas_eff_sv"
      ,new G4UnionSolid("gas_eff_region0+1_sv",gas_eff_region0_sv,gas_region1_sv,empty_rm,shift1)
      ,gas_region2_sv,empty_rm,shift2);
  auto* gas_eff_lv = new G4LogicalVolume(gas_eff_sv,m_materials["Ar+C4H10"],"gas_eff_lv");
  new G4PVPlacement(0
      ,G4ThreeVector(),gas_eff_lv,"gas_eff_pv",gas_lv,false,0);

  double PCB_sx = 23.3*cm;
  double PCB_sy = 23.3*cm;
  double PCB_thickness = 0.23*cm;
  G4ThreeVector PCB_frame_position(0,0,0.5*PCB_thickness);
  auto* PCB_frame_sv = new G4Box("PCB_frame_sv",.5*frame_sx,.5*frame_sy,.5*PCB_thickness);
  auto* PCB_frame_lv = new G4LogicalVolume(PCB_frame_sv,m_materials["Al"],"PCB_frame_lv");
  new G4PVPlacement(0,PCB_frame_position,PCB_frame_lv,"PCB_frame_pv",world_lv,false,0);
  auto* PCB_sv = new G4Box("PCB_sv",.5*PCB_sx,.5*PCB_sy,.5*PCB_thickness);
  auto* PCB_lv = new G4LogicalVolume(PCB_sv,m_materials["FR4"],"PCB_lv");
  new G4PVPlacement(0,G4ThreeVector(),PCB_lv,"PCB_pv",PCB_frame_lv,false,0);

  double MM_sz = 0.8*cm;
	double MM_gas_thickness = 0.6*cm;
	double MM_effective_sx = 15.*cm;
	//double MM_effective_sy = 15.*cm; /*unused*/
  double MM_frame_sx = 26.9*cm;
  double MM_frame_sy = 20.6*cm;
  double MM_frame_thickness_x = 1.5*cm;
  double MM_frame_thickness_y = 1.2*cm;
  double inner_frame_sy = 18.2*cm;
  G4ThreeVector MM_position(0,0,PCB_thickness+.5*MM_sz);
  auto* MM_sv = new G4Box("MM_sv",.5*MM_frame_sx,.5*MM_frame_sy,.5*MM_sz);
  auto* MM_lv = new G4LogicalVolume(MM_sv,m_materials["Al"],"MM_lv");
  new G4PVPlacement(0,MM_position,MM_lv,"MM_pv",world_lv,false,0);

  G4ThreeVector MM_gas0_position(0,0,PCB_thickness+.5*MM_gas_thickness);
  auto* MM_gas0_sv = new G4Box("MM_gas0_sv"
      ,.5*MM_frame_sx-MM_frame_thickness_x
      ,.5*inner_frame_sy-MM_frame_thickness_y
      ,.5*MM_gas_thickness);
  auto* MM_gas0_lv = new G4LogicalVolume(MM_gas0_sv,m_materials["Ar+C4H10"],"MM_gas0_lv");
  new G4PVPlacement(0,MM_gas0_position-MM_position,MM_gas0_lv,"MM_gas0_pv",MM_lv,false,0);

  auto* MM_gas_sv = new G4Box("MM_gas_sv",.5*MM_effective_sx,.5*MM_effective_sx,.5*MM_gas_thickness);
  auto* MM_gas_lv = new G4LogicalVolume(MM_gas_sv,m_materials["Ar+C4H10"],"MM_gas_lv");
  new G4PVPlacement(0,G4ThreeVector(),MM_gas_lv,"MM_gas_pv",MM_gas0_lv,false,0);

  double max_step = 1.0*mm;
  m_step_limits = new G4UserLimits(max_step);
  gas_eff_lv->SetUserLimits(m_step_limits);
  gas_lv->SetUserLimits(m_step_limits);

  auto* vis_attributes = new G4VisAttributes(G4Colour(0.9,0,0));
  vis_attributes->SetVisibility(false);
  //window_lv->SetVisAttributes(vis_attributes);
  world_lv->SetVisAttributes(vis_attributes);

  gas_lv->SetVisAttributes(new G4VisAttributes(G4Colour(1,0,1)));
  gas_eff_lv->SetVisAttributes(new G4VisAttributes(G4Colour(0,1,0)));  
  MM_gas_lv->SetVisAttributes(new G4VisAttributes(G4Color(1,0,1)));
  film_lv->SetVisAttributes(new G4VisAttributes(G4Colour(0.5,0,1)));
  return world_pv; }
 
/* override */
void detector_construction::ConstructSDandField(){
  m_s_magfield_messenger = new G4GlobalMagFieldMessenger(G4ThreeVector());
  m_s_magfield_messenger->SetVerboseLevel(1);
  G4AutoDelete::Register(m_s_magfield_messenger); }
