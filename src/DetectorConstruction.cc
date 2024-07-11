//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file radioactivedecay/rdecay01/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
// $Id: DetectorConstruction.cc 78307 2013-12-11 10:55:57Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "SteppingAction.hh"

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
#include "G4GlobalMagFieldMessenger.hh"
#include "G4UserLimits.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4AutoDelete.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0;

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{G4cout<<"<<------------DetectorConstruction::DetectorConstruction()-------------------->>"<<G4endl;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{G4cout<<"<<------------DetectorConstruction::~DetectorConstruction()-------------------->>"<<G4endl;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	G4cout<<"<<------------DetectorConstruction::Construct()-------------------->>"<<G4endl;
	//
	// World volume
	//   

	// Material ---> Vacuum  
	G4Material* Vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4Material* Al = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
	G4Material* Cu = G4NistManager::Instance()->FindOrBuildMaterial("G4_Cu"); 
	G4Material* Pb = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");  
	G4Material* Ge = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ge"); 
	G4Material* Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
	//Material --->CO2
	G4Material* CO2 = G4NistManager::Instance()->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
	//Material --->Ar
	G4Material* Ar = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ar");
	//Material --->Ne
	G4Material* Ne = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ne");
	G4Material* Xe = G4NistManager::Instance()->FindOrBuildMaterial("G4_Xe");

	//Materials defined by self
	G4String symbol;
	G4double a;                                              // atomic mass
	G4double z;                                              // atomic number
	G4double density;
	G4int ncomponents, natoms;
	G4Element* C  = new G4Element("Carbon",     symbol= "C",  z= 6.,  a= 12.00*g/mole);
	G4Element* H  = new G4Element("Hydrogen",   symbol= "H",  z= 1.,  a= 1.00*g/mole);
	G4Element* O  = new G4Element("Oxygen",     symbol= "O" , z= 8. , a= 16.00*g/mole);
	G4Element* N  = new G4Element("Nitrogen",     symbol= "N" , z= 7. , a= 14.00*g/mole);
	G4Element* Si = new G4Element("Si",			symbol= "Si", z= 14., a= 28.00*g/mole);
	G4Element* Fe = new G4Element("Iron", "Fe", z= 26., a= 55.85*g/mole);
	G4Element* Mn = new G4Element("Manganese", "Mn", z= 25., a= 54.94*g/mole);
	G4Element* Cr = new G4Element("Chromium", "Cr", z= 24., a= 52.00*g/mole);
	G4Element* Ni = new G4Element("Nickel", "Ni", z= 28., a= 58.7*g/mole);
	//Material--->PET C10H8O4
	G4Material* PET = new G4Material("PET", density=1.38*g/cm3, ncomponents=3); //PET 200um
	PET->AddElement(C, natoms=10);
	PET->AddElement(H, natoms=8);
	PET->AddElement(O, natoms=4);

	//Define FR4, from https://agenda.infn.it/event/14179/contributions/23405/attachments/16712/18966/Geant_FR4_Raffaella_v1.pdf
	//epoxy
	G4Material* Epoxy = new G4Material("Epoxy", density = 1.2*g/cm3, ncomponents = 2);
	Epoxy->AddElement(H, natoms = 2);
	Epoxy->AddElement(C, natoms = 2);
	//SiO2
	G4Material* SiO2 = new G4Material("SiO2", density = 2.200*g/cm3, ncomponents = 2);
	SiO2->AddElement(Si, natoms = 1);
	SiO2->AddElement(O, natoms = 2);
	//Material--->FR4
	G4Material* FR4 = new G4Material("FR4", density = 1.86*g/cm3, ncomponents = 2);
	FR4->AddMaterial(Epoxy, 47.2*perCent);
	FR4->AddMaterial(SiO2, 52.8*perCent);

	G4Material* Polyacrylate = new G4Material("Polyacrylate", density=1.39*g/cm3, ncomponents=3); //PET 200um
	Polyacrylate->AddElement(C, natoms=3);
	Polyacrylate->AddElement(H, natoms=4);
	Polyacrylate->AddElement(O, natoms=2);
	//Material ---> Atlasgas
	G4Material* Atlasgas = new G4Material("Atlasgas", density=1.79e-3*g/cm3, ncomponents=2);   //atlas gas 10mm
	Atlasgas->AddMaterial(CO2, 7.0*perCent);
	Atlasgas->AddMaterial(Ar,  93.0*perCent);

	//Material ---> C4H10
	G4Material* iC4H10 = new G4Material("iC4H10",  density=2.487e-3*g/cm3, ncomponents=2);  
	//density should be checked!!!
	iC4H10->AddElement(C, natoms=4);
	iC4H10->AddElement(H, natoms=10);

	//Material ---> Ne+C4H10
	G4Material* NeiC4H10 = new G4Material("NeiC4H10", density=0.98e-3*g/cm3, ncomponents=2);   
	NeiC4H10->AddMaterial(iC4H10,  5.0*perCent);
	NeiC4H10->AddMaterial(Ne,     95.0*perCent);

	//Material ---> Ar+iC4H10				//density is calculated by volume fraction 96.5/3.5
	G4Material* AriC4H10 = new G4Material("AriC4H10", density = (0.001782 * 0.965 + 0.00251 * 0.035)*g/cm3, ncomponents = 2);
	AriC4H10->AddMaterial(iC4H10, 	3.5*perCent);
	AriC4H10->AddMaterial(Ar,	  	96.5*perCent);

	//Material ---> Kapton
	G4Material* Kapton = new G4Material("Kapton", density = 1.43*g/cm3, ncomponents = 4);
	Kapton->AddElement(C, natoms = 22);
	Kapton->AddElement(H, natoms = 10);
	Kapton->AddElement(N, natoms = 2);
	Kapton->AddElement(O, natoms = 5);

	// wood
	G4Material* Wood = new G4Material("Wood", density=0.9*g/cm3, ncomponents=3);
	Wood->AddElement(H , 4);
	Wood->AddElement(O , 1);
	Wood->AddElement(C , 2);

	G4Material* mat304steel = new G4Material("Stainless steel 304", 7.999*g/cm3, 6);
	mat304steel->AddElement(Mn, 0.02);
	mat304steel->AddElement(Si, 0.01);
	mat304steel->AddElement(Cr, 0.19);
	mat304steel->AddElement(Ni, 0.10);
	mat304steel->AddElement(Fe, 0.6792);
	mat304steel->AddElement(C, 0.0008);

	//---------------construct detector---------------
	// Full sphere shape
	G4double solidWorld_rmax = 100*cm;
	G4Orb* solidWorld = new G4Orb(
			"World",                   // its name
			solidWorld_rmax);                 // its size 

	G4LogicalVolume* logicWorld = new G4LogicalVolume(
			solidWorld,             // its solid
			Air,                 // its material
			"World");               // its name
	G4VPhysicalVolume* physicalWorld = new G4PVPlacement(
			0,                        // no rotation
			G4ThreeVector(),          // at (0,0,0)
			logicWorld,               // its logical volume
			"World",                  // its name
			0,                        // its mother  volume
			false,                    // no boolean operation
			0);                       // copy number


	//common parameters of the aluminum frame
	G4double FrameSizeX = 26.9*cm;
	G4double FrameSizeY = 26.9*cm;
	G4double FrameSizeX_MM = 26.9*cm;
	G4double FrameSizeY_MM = 20.6*cm;
	G4double FrameThicknessX_MM = 1.5*cm;
	G4double FrameThicknessY_MM = 1.2*cm;
	G4double InnerFrameSizeX = 18.2*cm;
	G4double InnerFrameSizeY = 18.2*cm;
	G4double DetectorSizeX = 15*cm;
	G4double DetectorSizeY = 15*cm;
	G4double FrameStepHeight = 0.05*cm;
	G4double WindowSizeX = 16.6*cm;
	G4double WindowSizeY = 16.6*cm;

	// parameters of PET
	G4double petsizeX = InnerFrameSizeX;
	G4double petsizeY = InnerFrameSizeY;
	G4double petsizeZ = (3e-4)*cm;

	// Double-sided adhesive
	// parameters of Polyacrylate
	G4double polysizeX = InnerFrameSizeX;
	G4double polysizeY = InnerFrameSizeY;
	G4double polysizeZ = (1e-4)*cm;

	// parameters of TPC
	G4double frame1_height = 18.3*mm;
	G4double frame1_thickness = 41*mm;
	G4double frame2_height = 34.7*mm;
	G4double frame2_thickness = 12*mm;
	G4double frame3_height = 12*mm;
	G4double frame3_inner_radius = 58*mm;
	G4double window_radius = 35*mm;
	G4double window_thickness = 5*mm;
	G4double ring_thickness = 3*mm;

	G4double TPCSizeZ = 6.5*cm;	
	G4double TPCGasThickness = TPCSizeZ-frame3_height;
	G4double EffectiveSizeX = 12.*cm;
	G4double EffectiveSizeY = 12.*cm;

	// parameters of the micromegas
	G4double MMSizeZ = 0.8*cm;
	G4double MMGasThickness = 0.6*cm;
	G4double MMEffectiveSizeX = 15.*cm;
	G4double MMEffectiveSizeY = 15.*cm;

	// parameters of the PCB board
	G4double PCBsizeX = 23.3*cm;
	G4double PCBsizeY = 23.3*cm;
	G4double PCBdeltaY = 0*cm;
	G4double PCBthickness = 0.23*cm;

	// parameters of Cu board
	G4double CusizeX = PCBsizeX;
	G4double CusizeY = PCBsizeY;
	G4double CudeltaY = PCBdeltaY;
	// G4double Cuthickness = 10e-4*cm;
	G4double Cuthickness = 10e-4*cm;
	
	// parameter of the Cu board below TPC
	G4double gap = 0*cm;
	G4double Cuthickness2 = 1*cm;

	// parameter of the field cage
	G4double FieldCageSize = 16.6*cm;
	G4double FieldCageThickness = 0.16*cm;
	G4double FieldCageHeight = 5*cm;

	G4double source_thickness = 2*mm;
	G4double source_radius = 1*cm;

	//shielding shell parameters
	G4double ShellSize = 30*cm;
	G4double ShellThickness = 2*cm;
	G4double ShellHeight = 15*cm;

	// G4double Cuthickness3 = 2*cm;

	//whether to put the source under the TPC
	G4bool IfSource = false;

	// distance of the source bottom from the film
	G4double distance = 58*mm;


	
/*
*/


//=======Copper board below TPC (if exists)=====================

	// G4ThreeVector positionCubrd2 = G4ThreeVector(0., CudeltaY, -TPCSizeZ-gap-0.5*Cuthickness2);	

	// G4Box* solidCubrd2 = new G4Box("Cu_board2",                                    // its name
	// 		0.5*FrameSizeX, 0.5*FrameSizeY, 0.5*Cuthickness2);                      // its size

	// G4LogicalVolume* logicCubrd2 = new G4LogicalVolume(
	// 		solidCubrd2,                                    // its solid
	// 		Pb,                                    // its material
	// 		"Cu_board2");                                      // its name

	// new G4PVPlacement(
	// 		0,                                           // no rotation
	// 		positionCubrd2,                                 // at (0,0,0)
	// 		logicCubrd2,                                    // its logical volume
	// 		"Cu_board2",                                       // its name
	// 		logicWorld,                                  // its mother  volume
	// 		false,                                       // no boolean operation
	// 		0);                                          // copy number

//=============================================

//=======Copper board above anti-coincident detector (if exists)=====================

	// G4ThreeVector positionCubrd3 = G4ThreeVector(0., CudeltaY, PCBthickness+MMSizeZ+gap+0.5*Cuthickness2);	

	// G4Box* solidCubrd3 = new G4Box("Cu_board3",                                    // its name
	// 		0.5*FrameSizeX, 0.5*FrameSizeY, 0.5*Cuthickness2);                      // its size

	// G4LogicalVolume* logicCubrd3 = new G4LogicalVolume(
	// 		solidCubrd3,                                    // its solid
	// 		Pb,                                    // its material
	// 		"Cu_board3");                                      // its name

	// new G4PVPlacement(
	// 		0,                                           // no rotation
	// 		positionCubrd3,                                 // at (0,0,0)
	// 		logicCubrd3,                                    // its logical volume
	// 		"Cu_board3",                                       // its name
	// 		logicWorld,                                  // its mother  volume
	// 		false,                                       // no boolean operation
	// 		0);                                          // copy number

//=============================================

//=============Added part: a copper shell around the TPC ==================

	// G4ThreeVector positionShell = G4ThreeVector(0., 0., -2*cm);		//(0,0,-3.35)cm center, top at z=0

	// G4Box* solidShellOut = new G4Box("ShellOut",                                    // its name
	// 		0.5*ShellSize+ShellThickness, 0.5*ShellSize+ShellThickness, 0.5*ShellHeight+ShellThickness);                      // its size
	
	// G4Box* solidShellIn = new G4Box("ShellIn",                                    // its name
	// 		0.5*ShellSize, 0.5*ShellSize, 0.5*ShellHeight);                      // its size

	// G4SubtractionSolid* solidShell = new G4SubtractionSolid("Shell",solidShellOut,solidShellIn);
	
	// G4LogicalVolume* logicShell = new G4LogicalVolume(
	// 		solidShell,                                    // its solid
	// 		Pb,                                    // its material
	// 		"Shell");                                      // its name

	// new G4PVPlacement(
	// 		0,                                           // no rotation
	// 		positionShell,                                 // at (0,0,0)
	// 		logicShell,                                    // its logical volume
	// 		"Shell",                                       // its name
	// 		logicWorld,                                  // its mother  volume
	// 		false,                                       // no boolean operation
	// 		0);                                          // copy number

	// G4double Althickness = 4*mm;
	
	// G4Box* solidAlShellIn = new G4Box("AlShellIn",                                    // its name
	// 		0.5*ShellSize-Althickness, 0.5*ShellSize-Althickness, 0.5*ShellHeight-Althickness);                      // its size

	// G4SubtractionSolid* solidAlShell = new G4SubtractionSolid("AlShell",solidShellIn,solidAlShellIn);
	
	// G4LogicalVolume* logicAlShell = new G4LogicalVolume(
	// 		solidAlShell,                                    // its solid
	// 		Al,                                    // its material
	// 		"AlShell");                                      // its name

	// new G4PVPlacement(
	// 		0,                                           // no rotation
	// 		positionShell,                                 // at (0,0,0)
	// 		logicAlShell,                                    // its logical volume
	// 		"AlShell",                                       // its name
	// 		logicWorld,                                  // its mother  volume
	// 		false,                                       // no boolean operation
	// 		0);                                          // copy number

//==============================================================

//============Add wooden desk under the TPC ===============

	// G4double distance_desk = 7.6*cm;
	// G4double thickness_desk = 5*cm;
	// G4double size_desk = 50*cm;
	
	// G4ThreeVector positionDesk = G4ThreeVector(0., 0., -TPCSizeZ-distance_desk-0.5*thickness_desk);	

	// G4Box* solidDesk = new G4Box("Desk",                                    // its name
	// 		0.5*size_desk, 0.5*size_desk, 0.5*thickness_desk);                      // its size

	// G4LogicalVolume* logicDesk = new G4LogicalVolume(
	// 		solidDesk,                                    // its solid
	// 		Wood,                                    // its material
	// 		"Desk");                                      // its name

	// new G4PVPlacement(
	// 		0,                                           // no rotation
	// 		positionDesk,                                 // at (0,0,0)
	// 		logicDesk,                                    // its logical volume
	// 		"Desk",                                       // its name
	// 		logicWorld,                                  // its mother  volume
	// 		false,                                       // no boolean operation
	// 		0);                                          // copy number

//===================================================

//=============PART1: The TPC detector=========

	// Gas chamber
	G4ThreeVector positionTPC = G4ThreeVector(0., 0., -0.5*TPCSizeZ);		//(0,0,-3.35)cm center, top at z=0

	G4Box* solidTPC = new G4Box("TPCVolume",                                    // its name
			0.5*FrameSizeX, 0.5*FrameSizeY, 0.5*TPCSizeZ);                      // its size

	G4LogicalVolume* logicTPC = new G4LogicalVolume(
			solidTPC,                                    // its solid
			AriC4H10,                                    // its material
			"TPCVolume");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionTPC,                                 // at (0,0,0)
			logicTPC,                                    // its logical volume
			"TPCVolume",                                       // its name
			logicWorld,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	// ===== several parts of the Aluminum frame for the TPC side ====

	

	G4ThreeVector position_frame1 = G4ThreeVector(0., 0., -0.5*frame1_height);
	G4Box* solidFrame1Out = new G4Box("Frame1Out",                                    // its name
			0.5*FrameSizeX, 0.5*FrameSizeY, 0.5*frame1_height);                      // its size
	
	G4Box* solidFrame1In = new G4Box("Frame1In",                                    // its name
			0.5*FrameSizeX-frame1_thickness, 0.5*FrameSizeY-frame1_thickness, 0.5*frame1_height);                      // its size

	G4SubtractionSolid* solidFrame1 = new G4SubtractionSolid("Frame1",solidFrame1Out,solidFrame1In);
	
	G4LogicalVolume* logicFrame1 = new G4LogicalVolume(
			solidFrame1,                                    // its solid
			Al,                                    // its material
			"Frame1");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			position_frame1-positionTPC,                                 // at (0,0,0)
			logicFrame1,                                    // its logical volume
			"Frame1",                                       // its name
			logicTPC,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	G4ThreeVector position_frame2 = G4ThreeVector(0., 0., -frame1_height-0.5*frame2_height);
	G4Box* solidFrame2Out = new G4Box("Frame2Out",                                    // its name
			0.5*FrameSizeX, 0.5*FrameSizeY, 0.5*frame2_height);                      // its size
	
	G4Box* solidFrame2In = new G4Box("Frame2In",                                    // its name
			0.5*FrameSizeX-frame2_thickness, 0.5*FrameSizeY-frame2_thickness, 0.5*frame2_height);                      // its size

	G4SubtractionSolid* solidFrame2 = new G4SubtractionSolid("Frame2",solidFrame2Out,solidFrame2In);
	
	G4LogicalVolume* logicFrame2 = new G4LogicalVolume(
			solidFrame2,                                    // its solid
			Al,                                    // its material
			"Frame2");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			position_frame2-positionTPC,                                 // at (0,0,0)
			logicFrame2,                                    // its logical volume
			"Frame2",                                       // its name
			logicTPC,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	G4ThreeVector position_frame3 = G4ThreeVector(0., 0., -frame1_height-frame2_height-0.5*frame3_height);
	G4Box* solidFrame3Out = new G4Box("Frame3Out",                                    // its name
			0.5*FrameSizeX, 0.5*FrameSizeY, 0.5*frame3_height);                      // its size
	
	G4Tubs* solidFrame3In = new G4Tubs(
			"Frame3In",                   // its name
			0,                 // r_min 
			frame3_inner_radius,                 // r_max
			0.5*frame3_height,				//height
			0.,
			360.);

	G4SubtractionSolid* solidFrame3 = new G4SubtractionSolid("Frame3",solidFrame3Out,solidFrame3In);
	
	G4LogicalVolume* logicFrame3 = new G4LogicalVolume(
			solidFrame3,                                    // its solid
			Al,                                    // its material
			"Frame3");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			position_frame3-positionTPC,                                 // at (0,0,0)
			logicFrame3,                                    // its logical volume
			"Frame3",                                       // its name
			logicTPC,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	//the central round window

	G4ThreeVector positionwindow = G4ThreeVector(0., 0., -TPCSizeZ+0.5*window_thickness);

	G4Tubs* solidwindow = new G4Tubs(
			"Window",                   // its name
			window_radius,                 // r_min 
			frame3_inner_radius,                 // r_max
			0.5*window_thickness,				//height
			0.,
			360.);

	G4LogicalVolume* logicwindow = new G4LogicalVolume(
			solidwindow,                                    // its solid
			Al,                                    // its material
			"Window");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionwindow-positionTPC,                                 // at (0,0,0)
			logicwindow,                                    // its logical volume
			"Window",                                       // its name
			logicTPC,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	// film under the window
	G4ThreeVector positionfilm = G4ThreeVector(0., 0., -TPCSizeZ+window_thickness+0.5*petsizeZ);

	G4Tubs* solidfilm = new G4Tubs(
			"Film",                   // its name
			0,                 // r_min 
			frame3_inner_radius,                 // r_max
			0.5*petsizeZ,				//height
			0.,
			360.);

	G4LogicalVolume* logicfilm = new G4LogicalVolume(
			solidfilm,                                    // its solid
			PET,                                    // its material
			"Film");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionfilm-positionTPC,                                 // at (0,0,0)
			logicfilm,                                    // its logical volume
			"Film",                                       // its name
			logicTPC,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	// ring part below the window, for installation
	G4ThreeVector positionring = G4ThreeVector(0., 0., -TPCSizeZ+window_thickness+petsizeZ+0.5*ring_thickness);

	G4Tubs* solidring = new G4Tubs(
			"Ring",                   // its name
			window_radius,                 // r_min 
			frame3_inner_radius,                 // r_max
			0.5*ring_thickness,				//height
			0.,
			360.);

	G4LogicalVolume* logicring = new G4LogicalVolume(
			solidring,                                    // its solid
			Al,                                    // its material
			"Ring");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionring-positionTPC,                                 // at (0,0,0)
			logicring,                                    // its logical volume
			"Ring",                                       // its name
			logicTPC,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	// air outside of the film region
	G4ThreeVector positionStep = G4ThreeVector(0., 0., -TPCSizeZ+0.5*window_thickness);

	G4Tubs* solidStep = new G4Tubs(
			"Step",                   // its name
			0,                 // r_min 
			window_radius,                 // r_max
			0.5*window_thickness,				//height
			0.,
			360.);

	G4LogicalVolume* logicStep = new G4LogicalVolume(
			solidStep,                                    // its solid
			Air,                                    // its material
			"Step");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionStep-positionTPC,                                 // at (0,0,0)
			logicStep,                                    // its logical volume
			"Step",                                       // its name
			logicTPC,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number


	//=================Added part: a Aluminum plate to simulate the beta source ==========================

	if(IfSource){
		if(distance<5*mm) {
			// this means the source is inside of the TPC volume
			G4ThreeVector positionsource = G4ThreeVector(0., 0., 0.5*window_thickness-distance-0.5*source_thickness);

			G4Tubs* solidsource = new G4Tubs(
					"beta_source",                   // its name
					0,                 // r_min 
					source_radius,                 // r_max
					0.5*source_thickness,				//height
					0.,
					360.);

			G4LogicalVolume* logicsource = new G4LogicalVolume(
					solidsource,                                    // its solid
					Al,                                    // its material
					"beta_source");                                      // its name

			new G4PVPlacement(
					0,                                           // no rotation
					positionsource,                                 // at (0,0,0)
					logicsource,                                    // its logical volume
					"beta_source",                                       // its name
					logicStep,                                  // its mother  volume
					false,                                       // no boolean operation
					0);                                          // copy number
		}
		else {
			// this means the source is outside of the TPC volume
			G4ThreeVector positionsource = G4ThreeVector(0., 0., -TPCSizeZ+window_thickness-distance-0.5*source_thickness);

			G4Tubs* solidsource = new G4Tubs(
					"beta_source",                   // its name
					0,                 // r_min 
					source_radius,                 // r_max
					0.5*source_thickness,				//height
					0.,
					360.);

			G4LogicalVolume* logicsource = new G4LogicalVolume(
					solidsource,                                    // its solid
					Al,                                    // its material
					"beta_source");                                      // its name

			new G4PVPlacement(
					0,                                           // no rotation
					positionsource,                                 // at (0,0,0)
					logicsource,                                    // its logical volume
					"beta_source",                                       // its name
					logicWorld,                                  // its mother  volume
					false,                                       // no boolean operation
					0);                                          // copy number
		}
	}
	//=============================================



	// ========= TPC part frame constructed! ===========


	// field cage
	G4ThreeVector positionFieldCage = G4ThreeVector(0., 0., -0.5*FieldCageHeight);

	G4Box* solidFieldCageOut = new G4Box("FieldCageOut",                                    // its name
			0.5*FieldCageSize+FieldCageThickness, 0.5*FieldCageSize+FieldCageThickness, 0.5*FieldCageHeight);                      // its size
	
	G4Box* solidFieldCageIn = new G4Box("FieldCageIn",                                    // its name
			0.5*FieldCageSize, 0.5*FieldCageSize, 0.5*FieldCageHeight);                      // its size

	G4SubtractionSolid* solidFieldCage = new G4SubtractionSolid("FieldCage",solidFieldCageOut,solidFieldCageIn);
	
	G4LogicalVolume* logicFieldCage = new G4LogicalVolume(
			solidFieldCage,                                    // its solid
			FR4,                                    // its material
			"FieldCage");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionFieldCage-positionTPC,                                 // at (0,0,0)
			logicFieldCage,                                    // its logical volume
			"FieldCage",                                       // its name
			logicTPC,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number


	// detector sensitive volume
	G4ThreeVector positionGas = G4ThreeVector(0., 0., -0.5*TPCGasThickness);

	G4Box* solidGasBox = new G4Box("GasBox",                                    // its name
			0.5*DetectorSizeX, 0.5*DetectorSizeY, 0.5*TPCGasThickness);                      // its size

	G4Tubs* solidGasTub1 = new G4Tubs(
				"GasTub1",                   // its name
				0,                 // r_min 
				frame3_inner_radius,                 // r_max
				0.5*(frame3_height-window_thickness-petsizeZ-ring_thickness),				//height
				0.,
				360.);
	
	G4Tubs* solidGasTub2 = new G4Tubs(
				"GasTub2",                   // its name
				0,                 // r_min 
				window_radius,                 // r_max
				0.5*ring_thickness,				//height
				0.,
				360.);

	G4RotationMatrix *rm = new G4RotationMatrix;
	G4ThreeVector shift1 = G4ThreeVector(0., 0., -0.5*(TPCGasThickness+frame3_height-window_thickness-petsizeZ-ring_thickness));
	G4ThreeVector shift2 = G4ThreeVector(0., 0., -0.5*(TPCGasThickness+2*frame3_height-2*window_thickness-2*petsizeZ-1*ring_thickness));

	G4UnionSolid* solidGasAdd1 = new G4UnionSolid("GasAdd1",solidGasBox,solidGasTub1,rm,shift1);
	G4UnionSolid* solidGas = new G4UnionSolid("Gas",solidGasAdd1,solidGasTub2,rm,shift2);

	G4LogicalVolume* logicGas = new G4LogicalVolume(
			solidGas,                                    // its solid
			AriC4H10,                                    // its material
			"Gas");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionGas-positionTPC,                                 // at (0,0,0)
			logicGas,                                    // its logical volume
			"Gas",                                       // its name
			logicTPC,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number


	//selected TPC sensitive volume
	G4ThreeVector positionGasEff = G4ThreeVector(0., 0., -0.5*TPCGasThickness);
	G4Box* solidGasEffBox = new G4Box("GasEffBox",                                    // its name
			0.5*EffectiveSizeX, 0.5*EffectiveSizeY, 0.5*TPCGasThickness);                      // its size

	G4UnionSolid* solidGasEffAdd1 = new G4UnionSolid("GasEffAdd1",solidGasEffBox,solidGasTub1,rm,shift1);
	G4UnionSolid* solidGasEff = new G4UnionSolid("GasEff",solidGasEffAdd1,solidGasTub2,rm,shift2);

	G4LogicalVolume* logicGasEff = new G4LogicalVolume(
			solidGasEff,                                    // its solid
			AriC4H10,                                    // its material
			"GasEff");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			G4ThreeVector(),                                 // at (0,0,0)
			logicGasEff,                                    // its logical volume
			"GasEff",                                       // its name
			logicGas,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	
	

	//	PCB board

	G4ThreeVector positionPCBFrame = G4ThreeVector(0., 0., 0.5*PCBthickness);				//at (0,0,-6.8365)		
	G4Box* solidPCBFrame = new G4Box("PCBFrame",                                    // its name
			0.5*FrameSizeX, 0.5*FrameSizeY, 0.5*PCBthickness);                      // its size

	G4LogicalVolume* logicPCBFrame = new G4LogicalVolume(
			solidPCBFrame,                                    // its solid
			Al,                                    // its material
			"PCBFrame");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionPCBFrame,                                 // at (0,0,0)
			logicPCBFrame,                                    // its logical volume
			"PCBFrame",                                       // its name
			logicWorld,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	G4Box* solidPCB = new G4Box("PCB",                                    // its name
			0.5*PCBsizeX, 0.5*PCBsizeY, 0.5*PCBthickness);                      // its size

	G4LogicalVolume* logicPCB = new G4LogicalVolume(
			solidPCB,                                    // its solid
			FR4,                                    // its material
			"PCB");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			G4ThreeVector(),                                 // at (0,0,0)
			logicPCB,                                    // its logical volume
			"PCB",                                       // its name
			logicPCBFrame,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	//	Copper board

	// G4ThreeVector positionCubrd = G4ThreeVector(0., CudeltaY, PCBthickness+0.5*Cuthickness);	

	// G4Box* solidCubrd = new G4Box("Cu_board",                                    // its name
	// 		0.5*CusizeX, 0.5*CusizeY, 0.5*Cuthickness);                      // its size

	// G4LogicalVolume* logicCubrd = new G4LogicalVolume(
	// 		solidCubrd,                                    // its solid
	// 		Cu,                                    // its material
	// 		"Cu_board");                                      // its name

	// new G4PVPlacement(
	// 		0,                                           // no rotation
	// 		positionCubrd,                                 // at (0,0,0)
	// 		logicCubrd,                                    // its logical volume
	// 		"Cu_board",                                       // its name
	// 		logicWorld,                                  // its mother  volume
	// 		false,                                       // no boolean operation
	// 		0);                                          // copy number
			
	//	PCB board

	// G4ThreeVector positionPCB2 = G4ThreeVector(0., PCBdeltaY, PCBthickness+Cuthickness+0.5*PCBthickness);				//at (0,0,-6.8365)		
	// G4Box* solidPCB2 = new G4Box("PCB2",                                    // its name
	// 		0.5*PCBsizeX, 0.5*PCBsizeY, 0.5*PCBthickness);                      // its size

	// G4LogicalVolume* logicPCB2 = new G4LogicalVolume(
	// 		solidPCB2,                                    // its solid
	// 		Kapton,                                    // its material
	// 		"PCB2");                                      // its name

	// new G4PVPlacement(
	// 		0,                                           // no rotation
	// 		positionPCB2,                                 // at (0,0,0)
	// 		logicPCB2,                                    // its logical volume
	// 		"PCB2",                                       // its name
	// 		logicWorld,                                  // its mother  volume
	// 		false,                                       // no boolean operation
	// 		0);                                          // copy number

	
	// anticoincident MM

	G4ThreeVector positionMM = G4ThreeVector(0., 0., PCBthickness+0.5*MMSizeZ);

	G4Box* solidMM = new G4Box("MM",                                    // its name
			0.5*FrameSizeX_MM, 0.5*FrameSizeY_MM, 0.5*MMSizeZ);                      // its size

	G4LogicalVolume* logicMM = new G4LogicalVolume(
			solidMM,                                    // its solid
			Al,                                    // its material
			"MM");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionMM,                                 // at (0,0,0)
			logicMM,                                    // its logical volume
			"MM",                                       // its name
			logicWorld,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	// full gas volume of MM

	G4ThreeVector positionMMGas0 = G4ThreeVector(0., 0., PCBthickness+0.5*MMGasThickness);		
	
	G4Box* solidMMGas0 = new G4Box("MMGas0",
			0.5*FrameSizeX_MM-FrameThicknessX_MM, 0.5*InnerFrameSizeY-FrameThicknessY_MM, 0.5*MMGasThickness);

	G4LogicalVolume* logicMMGas0 = new G4LogicalVolume(
			solidMMGas0,                                    // its solid
			AriC4H10,                                    // its material
			"MMGas0");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			positionMMGas0-positionMM,                                 // at (0,0,0)
			logicMMGas0,                                    // its logical volume
			"MMGas0",                                       // its name
			logicMM,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number

	// sensitive volume of the MM		
	
	G4Box* solidMMGas = new G4Box("GasEff2",
			0.5*MMEffectiveSizeX, 0.5*MMEffectiveSizeX, 0.5*MMGasThickness);

	G4LogicalVolume* logicMMGas = new G4LogicalVolume(
			solidMMGas,                                    // its solid
			AriC4H10,                                    // its material
			"GasEff2");                                      // its name

	new G4PVPlacement(
			0,                                           // no rotation
			G4ThreeVector(),                                 // at (0,0,0)
			logicMMGas,                                    // its logical volume
			"GasEff2",                                       // its name
			logicMMGas0,                                  // its mother  volume
			false,                                       // no boolean operation
			0);                                          // copy number





// ==========================================================



// ===============================================================


	//-----Set the step limits in the Gas volume-------------
	G4double maxStep = 1.0*mm;
	fStepLimits = new G4UserLimits(maxStep);
	logicGasEff->SetUserLimits(fStepLimits);
	logicGas->SetUserLimits(fStepLimits);
	//-------------------------------



	G4VisAttributes* visAttributes = new G4VisAttributes(G4Colour(0.9,0.0,0.0));
	visAttributes->SetVisibility(false);
	logicWorld->SetVisAttributes(visAttributes);

	// visAttributes = new G4VisAttributes(G4Colour(0.0,0.0,1.0)); // red
	// logicCollimation->SetVisAttributes(visAttributes);

	visAttributes = new G4VisAttributes(G4Colour(1.0,0.0,1.0)); 
	logicGas->SetVisAttributes(visAttributes);
	visAttributes = new G4VisAttributes(G4Colour(0.0,1.0,0.0)); 
	logicGasEff->SetVisAttributes(visAttributes);
	visAttributes = new G4VisAttributes(G4Colour(1.0,0.0,1.0)); 
	logicMMGas->SetVisAttributes(visAttributes);
	visAttributes = new G4VisAttributes(G4Colour(0.5,0.0,1.0)); 
	logicfilm->SetVisAttributes(visAttributes);


	//
	//always return the physical World
	//  
	return physicalWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField(){
	// Create global magnetic field messenger.
  	// Uniform magnetic field is then created automatically if
  	// the field value is not zero.
  	G4ThreeVector fieldValue = G4ThreeVector();
  	fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  	fMagFieldMessenger->SetVerboseLevel(1);
	
  	// Register the field messenger for deleting
  	G4AutoDelete::Register(fMagFieldMessenger);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
