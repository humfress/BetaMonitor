// BM_Detector.cc - Implementation of BM_Detector class: BetaMon detector and geometry

// Field handling
#include "G4Mag_UsualEqRhs.hh"
#include "G4EqMagElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"

// Integrators
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4PropagatorInField.hh"

// Geometry
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVParameterised.hh"
#include "G4TransportationManager.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4GDMLParser.hh"

#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "BM_Detector.hh"
#include "BM_SD.hh"
#include "BM_SteppingAction.hh"

/*
  NOTE: Nov 2025 - many block comments & other old geometries were removed to
  make the code more readable.  IF you need info on an old geometry: refer to ../test/BM_Detector.cc.

  running now:
  "SiPMBlake": SmallScin, Al, cone flange

  older geometries:
  OrigMon : top round, bottom round, square, Cu, normal flange
  SiPMWinston :  LargeScin, Cu, normal flange
  TeleOld ("telescope") : thinTele, Cu, normal flange
  TeleTAMU : TAMUdE, ?, normal flange?
  Beamline : adds the flange and decay volume geometry
  TabSource : adds table source geometry

  TODO: implement options to turn on/off certain volumes.
  Need a separate function for detector, bore, and source

  TODO: rework such that port is at 0,0,0?
*/

// MagField* BM_Detector::fMagneticField = 0;
G4FieldManager *BM_Detector::fFieldMgr = 0;
G4UniformMagField *BM_Detector::fMagneticField = 0;

BM_Detector::BM_Detector() : G4VUserDetectorConstruction(), vacuumLV(nullptr), 
vacuumWindowLV(nullptr), aScintillatorLV(nullptr), bScintillatorLV(nullptr), logicWorld(nullptr) {}

BM_Detector::~BM_Detector() {}

G4VPhysicalVolume *BM_Detector::Construct()
{
  G4cout << "Running sim with Blake's SiPM geometry: SmallScin, Al, cone flange" << G4endl;

  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // NIST material manager
  G4NistManager *nist = G4NistManager::Instance();
  G4Material *Fe = nist->FindOrBuildMaterial("G4_Fe");
  G4Material *PVT = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  // G4Material *Tape = nist->FindOrBuildMaterial("G4_POLYVINYL_CHLORIDE");
  G4Material *Mylar = nist->FindOrBuildMaterial("G4_MYLAR");
  G4Material *Al = nist->FindOrBuildMaterial("G4_Al");
  G4Material *Cr = nist->FindOrBuildMaterial("G4_Cr");
  G4Material *Ni = nist->FindOrBuildMaterial("G4_Ni");
  G4Material *Cu = nist->FindOrBuildMaterial("G4_Cu");
  G4Material *Mn = nist->FindOrBuildMaterial("G4_Mn");
  // G4Material *Si = nist->FindOrBuildMaterial("G4_Si");
  // G4Material *Kap = nist->FindOrBuildMaterial("G4_KAPTON");
  G4double density = 8030. * mg / cm3;
  G4Material *Stainless_Steel = new G4Material("Stainless_Steel", density, 4);

  // auto pcb = new G4Element("Vetronite", 2. * g / cm3, 2); // Fiber glass
  auto el_PCB_Si = new G4Element("Silicon", "Si", 14., 28.0855 * g / mole);
  auto el_PCB_O = new G4Element("Oxygen", "O", 8., 15.9994 * g / mole);
  auto mat_PCB = new G4Material("Vetronite", 2. * g / cm3, 2); // Fiber glass
  mat_PCB->AddElement(el_PCB_Si, 1);
  mat_PCB->AddElement(el_PCB_O, 2);
  Stainless_Steel->AddMaterial(Fe, 72 * perCent);
  Stainless_Steel->AddMaterial(Cr, 18 * perCent);
  Stainless_Steel->AddMaterial(Ni, 8 * perCent);
  Stainless_Steel->AddMaterial(Mn, 2 * perCent);

  // Vacuum
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008 * g / mole;
  G4double density2 = 1.e-25 * g / cm3;
  G4double temperature = 2.73 * kelvin;
  G4double pressure = 1.3332e-8 * pascal;
  G4Material *Vacuum = new G4Material("interGalactic", atomicNumber, massOfMole, density2, kStateGas, temperature, pressure);

  // World
  G4double world_sizeXY = 40 * cm;
  G4double world_sizeZ = 45 * cm;
  G4Material *world_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Box *solidWorld = new G4Box("World", 0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);
  logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World"); // (G4 volume instance, material, name)

  // set maximum step - SetUserLimits(Limits) should be called for every G4LogicalVolume we create
  G4UserLimits *Limits = new G4UserLimits(10 * um); 
  logicWorld->SetUserLimits(Limits);

  // for reference - give the parameter definitions for G4PVPlacement once
  G4VPhysicalVolume *physWorld =
      new G4PVPlacement(0,               // no rotation
                        G4ThreeVector(), // origin (default at (0,0,0))
                        logicWorld,      // its logical volume
                        "World",         // its name
                        0,               // its mother volume
                        false,           // no boolean operation
                        0,               // copy number
                        checkOverlaps);  // overlaps checking

  // ===========================================================================
  // === Declare Geometry components (each is assoc. with a G4LogicalVolume) ===
  // NOTE: To temporarily remove something from the geometry, comment out the G4LV, GVPVPLacement, and SetUserLimits lines.
  
  // === T Pipe Inner Vacuum (a.k.a. Decay Volume) == 
  G4double cyl_hdv = 6.75 * 2.54 * cm;    // height of decay volume 16.51 standard
  G4double Tdv_r1i = 1.7399 * cm;     // T Pipe inner radius
  G4double Tdv_h1 = 6.75 * 2.54 * cm;     // T Pipe Major axis length
  G4double Tdv_h2 = 6.75 * 2.54 / 2 * cm; // T Pipe Minor axis length
  G4double flange_width = 1.27 * cm; // Flange width
  G4double tol  = 1e-6 * mm; // generic tolerance used throughout
  G4double scint_z = 3 * mm; // scintillator thickness
  G4double scint_xy =  20 * mm; // scintillator size 
  G4double aScint2window_z = 8.4074 * mm; // distance from scintillatorA to window 
  G4double bScint2source_z = 2.642 * cm; // distacne from scintillatorB to calibration source
  G4double Tdv_r1o = 3.81 / 2 * cm;       // T Pipe outer radius
  G4double cyl_r2dvo = 6.9088 / 2 * cm; // Flange Outer Radius
  G4double cyl_r1dvi = 3.4798 / 2 * cm; // Inner Radius of T Pipe Decay Volume 
  G4double cyl_hc = 0.00762 * cm;     // thickness of aluminum window
  G4double mylar_t = 0.076 * mm; // scint mylar thickness
  G4double aluminum_t = 0.0001 * mm; // scint aluminum thickness 
  G4double sipmLid_ri = 10.16 / 2 * mm; //siomLid inner radius 
  G4double sipmLid_ro = 120.65 / 2 * mm; // sipmlid outer radius 
  G4double sipmLid_z = 8.89 * mm; // sipmLid thickness
  G4double sipmSourceLip = 0.12064 * cm;
  G4double source_ro = 12.65 / 2 * mm; // source assembly outer radius 
  G4double scintHolder_xy = 22.2 * mm; // scint holder xy
  G4double scintHolderWall_t = 0.85 * mm; // thickness of scint holder walls
  G4double scintHolderLip_z = 0.5 * mm; // lip thickness that scint sits on
  G4double scintHolderLip_xy = 2.1 * mm; // lip xy - how much the lip extends in the x and y direction beyond the scint (so that the scint can rest on it)
  G4double scintHolderBase_ro = 45.72 / 2 * mm; 
  G4double scintHolderBase_z = 1.575 * mm;
  G4double scintHolder_z = bScint2source_z + 4*(scint_z/2 + mylar_t + aluminum_t * 2) - sipmSourceLip - scintHolderBase_z + scintHolderLip_z; // scint holder z based on 
  G4double scintHolderBase_ri = 20 / 2 * mm;
  G4double scintSleeve_xy = scintHolder_xy - 2 * scintHolderWall_t;
  G4double scintSleeve_z = 10 * mm; 
  G4double scintSleeve_t = 3.5 * mm;
  G4double source_r = 0.939 / 2 * cm; // radius of mylar disk that the source is on
  G4double sourceMy_z = 0.00032 * 2 * cm; // thickness of mylar disk that the source is on
  G4double sourceAlo_z = 0.31749 * cm; 
  G4double sourceAl_z = 0.635 * cm * 2;
  G4double sourceAl2_r = 0.31623 * cm;
  G4double sourceAl2_cutoutdepth = 0.635 * cm;
  G4double sourceSr_z = 1e-3*mm; // used to establish volume in which we can generate particles using GPS.


  G4ThreeVector TransT(- Tdv_h2 / 2 - flange_width / 2, 0, 0);
  G4RotationMatrix *yRotT = new G4RotationMatrix;
  yRotT->rotateY(90.0 * deg); // Rotates 90 degrees
  G4Tubs *solidShapeT1i = new G4Tubs("Pipe1i",      // 
                                      0,            // rmin
                                      Tdv_r1i,      // rmax
                                      Tdv_h1 / 2 + flange_width,  // delta-z
                                      0,            // start-phi
                                      360 * deg);   // delta-phi
  G4Tubs *solidShapeT2i = new G4Tubs("Pipe2i", 0, Tdv_r1i, Tdv_h2 / 2 + flange_width / 2, 0, 360 * deg);
  G4UnionSolid *Pipei = new G4UnionSolid("InnerTPipe", solidShapeT1i, solidShapeT2i, yRotT, TransT);

 // create logical volume for decay volume and place it in world
  vacuumLV = new G4LogicalVolume(Pipei, Vacuum, "Vacuum"); 
  new G4PVPlacement(yRotT, G4ThreeVector(0 * cm, 0 * cm, Tdv_h2 + flange_width),
                    vacuumLV, "Vacuum", logicWorld, false, 1, checkOverlaps);
  vacuumLV->SetUserLimits(Limits);


  // === T Pipe Outer (steel) ===
  // includes a 1NM tolerance between the envelope defined above and the decay volume t-pipe. 
  G4ThreeVector TransT2(-Tdv_h2 / 2, 0, 0);
  G4Tubs *solidShapeT1o = new G4Tubs("Pipe1o", 0, Tdv_r1o, Tdv_h1 / 2., 0, 360 * deg);
  G4Tubs *solidShapeT1i_tol = new G4Tubs("Pipe1i_tol", 0, Tdv_r1i + tol, Tdv_h1 / 2 + flange_width / 2, 0, 360 * deg);
  G4Tubs *solidShapeT2o = new G4Tubs("Pipe2o", 0, Tdv_r1o, Tdv_h2 / 2, 0, 360 * deg);
  G4Tubs *solidShapeT2i_tol = new G4Tubs("Pipe2i_tol", 0, Tdv_r1i + tol, Tdv_h2 / 2 + flange_width / 2, 0, 360 * deg);
  G4UnionSolid *Pipeo = new G4UnionSolid("OuterTPipe", solidShapeT1o, solidShapeT2o, yRotT, TransT2);
  G4SubtractionSolid *TPipe1 = new G4SubtractionSolid("TDecayVolume", Pipeo, solidShapeT1i_tol);
  G4SubtractionSolid *TPipe = new G4SubtractionSolid("TDecayVolume", TPipe1, solidShapeT2i_tol, yRotT, TransT2);
  
  G4LogicalVolume *steelTPipeLV = new G4LogicalVolume(TPipe, Stainless_Steel, "SteelTPipe");
  new G4PVPlacement(yRotT, G4ThreeVector(0, 0, Tdv_h2 + flange_width),
                    steelTPipeLV, "SteelTPipe", logicWorld, false, 0, checkOverlaps);
  steelTPipeLV->SetUserLimits(Limits);


  // === T Pipe Flanges (steel) === 

  G4Tubs *Flange = new G4Tubs("Flange", cyl_r1dvi, cyl_r2dvo, flange_width / 2, 0, 360 * deg);
  
  // leftmost flange (nearest beta monitor)
  G4LogicalVolume *tpipeFlangeScintLV = new G4LogicalVolume(Flange, Stainless_Steel, "tPipeFlangeScint");
  new G4PVPlacement(0, G4ThreeVector(0 *  cm, 0 * cm, -flange_width / 2), 
                   tpipeFlangeScintLV, "tPipeFlangeScint", logicWorld,
                    false, 0, checkOverlaps);
  tpipeFlangeScintLV->SetUserLimits(Limits);

  // second to leftmost flange (sandwiches the "window")
  G4LogicalVolume *tPipeFlange1LV = new G4LogicalVolume(Flange, Stainless_Steel, "tPipeFlange1");
  new G4PVPlacement(0, G4ThreeVector(0 * cm, 0 * cm, flange_width / 2), 
                    tPipeFlange1LV, "tPipeFlange1", logicWorld, false, 0, checkOverlaps);
  tPipeFlange1LV->SetUserLimits(Limits);

  // right flange
  G4LogicalVolume *tPipeFlange2LV = new G4LogicalVolume(Flange, Stainless_Steel, "tPipeFlange2");
  new G4PVPlacement(yRotT, G4ThreeVector(Tdv_h1/2 + flange_width / 2, 0, cyl_hdv / 2 + flange_width), tPipeFlange2LV, "tPipeFlange2",
                    logicWorld, false, 0, checkOverlaps);
  tPipeFlange2LV->SetUserLimits(Limits);

  // middle flange
  G4LogicalVolume *tPipeFlange3LV = new G4LogicalVolume(Flange, Stainless_Steel, "tPipeFlange3");
  new G4PVPlacement(yRotT, G4ThreeVector(-Tdv_h1/2 - flange_width / 2, 0, cyl_hdv / 2 + flange_width),
                    tPipeFlange3LV, "tPipeFlange3", logicWorld, false, 0, checkOverlaps);
  tPipeFlange3LV->SetUserLimits(Limits);


  // === Vacuum Window (copper) - between two left flanges ===
  // G4double cyl_r2c = 3.556 / 2 * cm; // copper seal radius //Tdv_r1i = 3.4798 / 2 * cm;
  // G4double cyl_hkap = 0.0012 * cm;   // thickness of kapton (trials - 1:0.0012, 2: 0.00075, 3: 0.006, 4: 0.0127)
  G4Tubs *vacuumWindowDisk = new G4Tubs("VacuumWindowDisk", 0. * cm, Tdv_r1i, cyl_hc / 2., 0, 360 * deg);
  vacuumWindowLV = new G4LogicalVolume(vacuumWindowDisk, Al, "VacuumWindow");
  new G4PVPlacement(0, G4ThreeVector(0 * cm, 0 * cm, -cyl_hc / 2), vacuumWindowLV, "VacuumWindow", logicWorld, false, 2, checkOverlaps);
  vacuumWindowLV->SetUserLimits(Limits);
  //Wedge shape at window
  G4double wedge_ri = 0.78 * 2.54 /2 * cm; // inner radius of wedge (matches inner radius of T pipe)
  const G4int nZ = 3;
  G4double zPlane[nZ] = { 0*mm, 1.0*mm, 7.0*mm };
  G4double rInner[nZ] = { wedge_ri, wedge_ri, Tdv_r1i-tol };
  G4double rOuter[nZ] = { Tdv_r1i, Tdv_r1i, Tdv_r1i };

  // full ring
  auto solidRing = new G4Polycone(
      "WedgeRevolvedRing",
      0.0*deg,      // start phi
      360.0*deg,    // delta phi
      nZ,
      zPlane,
      rInner,
      rOuter
  );

  auto ringLV = new G4LogicalVolume(solidRing, Stainless_Steel, "WedgeRevolvedRingLV");
  
  new G4PVPlacement(nullptr, G4ThreeVector(), ringLV, "WedgeRevolvedRingPV",
                  logicWorld, false, 0, checkOverlaps);

  // === Scintillators (A & B) ===

  // bare scintillator
  G4Box *SmallScin = new G4Box("SmallScin", scint_xy / 2, scint_xy / 2, scint_z / 2);
  
  // the scintillators are wrapped in aluminized mylar (a sandwich of Al - My - Al)


  G4Box *SmallScinM = new G4Box("SmallScinM", // mylar
                                scint_xy / 2 + mylar_t + aluminum_t, 
                                scint_xy / 2 + mylar_t + aluminum_t,
                                scint_z / 2 + mylar_t + aluminum_t);
  G4Box *SmallScinMA = new G4Box("SmallScinMA", // aluminum 
                                  scint_xy / 2 + aluminum_t, 
                                  scint_xy / 2 + aluminum_t,
                                  scint_z / 2 + aluminum_t);                                  
  G4Box *SmallScinMA2 = new G4Box("SmallScinMAOuter", // aluminum
                                  scint_xy / 2 + mylar_t + aluminum_t * 2,
                                  scint_xy / 2 + mylar_t + aluminum_t * 2,
                                  scint_z / 2  + mylar_t + aluminum_t * 2);                                                             
  
  // inner scintillator "A" (inner = closer to window)
  // 8.4074 mm is the distnce between the window and the surface of the scint. 
  G4ThreeVector aScintillatorPos = G4ThreeVector(0 * cm, 0 * cm, - aScint2window_z - (scint_z / 2  + mylar_t + aluminum_t * 2));
  aScintillatorLV = new G4LogicalVolume(SmallScin, PVT, "AScintillator");
  new G4PVPlacement(0, aScintillatorPos, aScintillatorLV, "AScintillator", logicWorld, false, 3, checkOverlaps);
  aScintillatorLV->SetUserLimits(Limits);
  
  // inner aluminum layer for scint A
  G4SubtractionSolid *solidAlMyI = new G4SubtractionSolid("InnerAl", SmallScinMA, SmallScin);
  G4LogicalVolume *logicAlMylarAli1 = new G4LogicalVolume(solidAlMyI, Al, "Ali_sq1");
  new G4PVPlacement(0, aScintillatorPos, logicAlMylarAli1, "Ali_sq1", logicWorld, false, 0, checkOverlaps);
  logicAlMylarAli1->SetUserLimits(Limits);

  // mylar layer for scint A
  G4SubtractionSolid *solidAlMy = new G4SubtractionSolid("Mylar", SmallScinM, SmallScinMA);
  G4LogicalVolume *logicAlMylarMy1 = new G4LogicalVolume(solidAlMy, Mylar, "Mylar_sq1");
  new G4PVPlacement(0, aScintillatorPos, logicAlMylarMy1, "Mylar_sq1", logicWorld, false, 0, checkOverlaps);
  logicAlMylarMy1->SetUserLimits(Limits);

  // outer aluminum layer for scint A -- colors Scint A "blue"
  G4SubtractionSolid *solidAlMyO = new G4SubtractionSolid("OuterAl", SmallScinMA2, SmallScinM);
  G4LogicalVolume *logicAlMylarAlo1 = new G4LogicalVolume(solidAlMyO, Al, "Alo_sq1");
  new G4PVPlacement(0, aScintillatorPos, logicAlMylarAlo1, "Alo_sq1", logicWorld, false, 0, checkOverlaps);
  logicAlMylarAlo1->SetUserLimits(Limits);
  
  // outer scintillator "B" (outer = farther from window)
  // G4ThreeVector bScintillatorPos = G4ThreeVector(0 * cm, 0 * cm, -(3.317 * mm + (3.0 + 0.0762) * 3 / 2 * mm));

  // G4ThreeVector bScintillatorPos = G4ThreeVector(0 * cm, 0 * cm,  - aScint2window_z - 4 / 2 * (scint_z + 2*mylar_t + 4*aluminum_t));
  G4ThreeVector bScintillatorPos = G4ThreeVector(0 * cm, 0 * cm, - aScint2window_z - 3 * (scint_z / 2  + mylar_t + aluminum_t * 2));
  bScintillatorLV = new G4LogicalVolume(SmallScin, PVT, "BScintillatorLV");
  new G4PVPlacement(0, bScintillatorPos, bScintillatorLV, "BScintillator", logicWorld, false, 4, checkOverlaps);
  bScintillatorLV->SetUserLimits(Limits);

  // aluminum layer for scint B
  G4SubtractionSolid *solidAlMyI2 = new G4SubtractionSolid("InnerAl2", SmallScinMA, SmallScin);
  G4LogicalVolume *logicAlMylarAli2 = new G4LogicalVolume(solidAlMyI2, Al, "Ali_sq2");
  new G4PVPlacement(0, bScintillatorPos, logicAlMylarAli2, "Ali_sq2", logicWorld, false, 0, checkOverlaps);
  logicAlMylarAli2->SetUserLimits(Limits);          
  
  // mylar layer for scint B
  G4SubtractionSolid *solidAlMy2 = new G4SubtractionSolid("Mylar2", SmallScinM, SmallScinMA);
  G4LogicalVolume *logicAlMylarMy2 = new G4LogicalVolume(solidAlMy2, Mylar, "Mylar_sq2");
  new G4PVPlacement(0, bScintillatorPos, logicAlMylarMy2, "Mylar_sq2", logicWorld, false, 0, checkOverlaps);
  logicAlMylarMy2->SetUserLimits(Limits);

  // outer aluminum layer for scint B. -- colors scint B red
  G4SubtractionSolid *solidAlMyO2 = new G4SubtractionSolid("OuterAl2", SmallScinMA2, SmallScinM);
  G4LogicalVolume *logicAlMylarAlo2 = new G4LogicalVolume(solidAlMyO2, Al, "Alo_sq2");
  new G4PVPlacement(0, bScintillatorPos, logicAlMylarAlo2, "Alo_sq2", logicWorld, false, 0, checkOverlaps);
  logicAlMylarAlo2->SetUserLimits(Limits);

  
  // === Copper Back Plate / "SiPM Lid" ===
  G4RotationMatrix *Rot0 = new G4RotationMatrix;
  Rot0->rotateX(0 * rad); // seems unnecessary

  G4Tubs *sipmLidBase = new G4Tubs("sipmLidBase", 
                                    sipmLid_ri, //rmin
                                    sipmLid_ro, //rmax
                                    sipmLid_z / 2, //delta-z
                                    0, //start phi
                                    360 * deg); //end phi
  G4Tubs *sipmLidSource = new G4Tubs("sipmLidSource", 
                                      0, //rmin
                                      source_ro + tol, //rmax
                                      sipmLid_z / 2, //delta-z
                                      0, //start phi
                                      360 * deg); //end phi
  // source cutout is moved in the negative z dir by the half the width of the 
  // -32.079 (position of front of sipmLidBase) - (-34.873 [pos of source mylar] + 3.1749/2 [half width of sourceAlo ring]) = 1.2065 cm

  G4SubtractionSolid *sipmLid = new G4SubtractionSolid("sipmLid", sipmLidBase, sipmLidSource, Rot0, G4ThreeVector(0, 0, -sipmSourceLip));
  G4LogicalVolume *logicSiPMLid = new G4LogicalVolume(sipmLid, Cu, "sipm_Lid");
  // TODO pick up geometry fixing here. 
  
  G4ThreeVector sipmLidPos = G4ThreeVector(0, 0, -aScint2window_z - bScint2source_z - 4*(scint_z/2 + mylar_t + aluminum_t * 2) + sipmSourceLip - sipmLid_z / 2);
  new G4PVPlacement(0, sipmLidPos, logicSiPMLid, "sipm_Lid",
                    logicWorld, false, 0, checkOverlaps);
  logicSiPMLid->SetUserLimits(Limits);


  // === Scintillator Holder - (purple block with base, green sleeve, cutout) ===

  G4Box *scintHolderblock = new G4Box("scintHolderbase",
                                      scintHolder_xy / 2, scintHolder_xy / 2,
                                      scintHolder_z / 2);
  G4Box *scintHoldercutout1 = new G4Box("scintHoldercutout1",
                                        scintHolder_xy / 2 - scintHolderWall_t,
                                        scintHolder_xy / 2 - scintHolderWall_t,
                                        scintHolder_z / 2 - scintHolderLip_z);
  G4Box *scintHoldercutout2 = new G4Box("scintHoldercutout2",
                                        scintHolder_xy / 2 - scintHolderWall_t - scintHolderLip_xy,
                                        scintHolder_xy / 2 - scintHolderWall_t - scintHolderLip_xy,
                                        scintHolder_z / 2 + 10);
  G4SubtractionSolid *scintHolder1 = new G4SubtractionSolid("scintHolder1", scintHolderblock, scintHoldercutout2);
  G4SubtractionSolid *scintHolder = new G4SubtractionSolid("scintHolder", scintHolder1, scintHoldercutout1);
  G4LogicalVolume *logicScintHolder = new G4LogicalVolume(scintHolder, mat_PCB, "scint_holder");
  G4ThreeVector scintHolderPos = G4ThreeVector(0, 0, -aScint2window_z + scintHolderLip_z - scintHolder_z / 2);
  new G4PVPlacement(0, scintHolderPos, logicScintHolder, "scint_holder",
                    logicWorld, false, 0, checkOverlaps);
  logicScintHolder->SetUserLimits(Limits);

  G4Tubs *scintHolderBase = new G4Tubs("scintHolderbase", scintHolderBase_ri, scintHolderBase_ro, scintHolderBase_z / 2, 0, 360 * deg);
  G4LogicalVolume *logicScintHolderBase = new G4LogicalVolume(scintHolderBase, mat_PCB, "scint_holderBase");
  G4ThreeVector scintHolderBasePos = G4ThreeVector(0, 0, -aScint2window_z - bScint2source_z - 4*(scint_z/2 + mylar_t + aluminum_t * 2) + sipmSourceLip + scintHolderBase_z / 2);
  new G4PVPlacement(0, scintHolderBasePos, logicScintHolderBase, "scint_holderBase",
                    logicWorld, false, 0, checkOverlaps);
  logicScintHolderBase->SetUserLimits(Limits);
  

  G4Box *scintSleevebase = new G4Box("scintSleevebase",  scintSleeve_xy / 2,  scintSleeve_xy / 2, scintSleeve_z / 2);
  G4Box *scintSleevecutout = new G4Box("scintHoldercutout2",  scintSleeve_xy / 2 - scintSleeve_t, scintSleeve_xy / 2 - scintSleeve_t, scintSleeve_z);
  G4SubtractionSolid *scintSleeve = new G4SubtractionSolid("scintSleeve", scintSleevebase, scintSleevecutout);
  G4LogicalVolume *logicScintSleeve = new G4LogicalVolume(scintSleeve, Al, "scint_Sleeve");

  G4ThreeVector scintSleevePos = G4ThreeVector(0 * cm, 0 * cm, - aScint2window_z - 4*(scint_z / 2  + mylar_t + aluminum_t * 2) - scintSleeve_z / 2);
  new G4PVPlacement(0, scintSleevePos, logicScintSleeve, "scint_Sleeve",
                    logicWorld, false, 0, checkOverlaps);
  logicScintSleeve->SetUserLimits(Limits);


  // === Calibration source (typically Sr90) on "back" copper plane ===


  // source_ro for outer radius of source cutout in the SiPM lid is defined above as 12.65/2 mm.

   // source cylinder (aluminum)
  G4Tubs *SourceAl1 = new G4Tubs("SourceAl1", 
                                  0 * cm,         // rmin
                                  source_ro,  // rmax
                                  sourceAl_z / 2,     // delta-z
                                  0, 360 * deg);  // start-phi, delta-phi
  G4Tubs *SourceAl2 = new G4Tubs("SourceAl2",  
                                  0 * cm, 
                                  sourceAl2_r, 
                                  sourceAlo_z, // just a coincidence that this is same as sourceAlo_z. All that matter is that it matches ztrans below.  
                                  0, 360 * deg);  
  G4Tubs *SourceAlo = new G4Tubs("SourceAlo", 
                                  source_r - tol, 
                                  source_ro, 
                                  sourceAlo_z / 2, 
                                  0, 360 * deg);

  G4LogicalVolume *logicSourceAlo = new G4LogicalVolume(SourceAlo, Al, "SourceAlo");
  G4ThreeVector sourceAloPos = G4ThreeVector(0, 0, -aScint2window_z - bScint2source_z - 4*(scint_z/2 + mylar_t + aluminum_t * 2) - sourceAlo_z / 2 - 1*tol);
  new G4PVPlacement(0, sourceAloPos, logicSourceAlo, "SourceAlo", logicWorld, false, 0, checkOverlaps);
  logicSourceAlo->SetUserLimits(Limits);

  G4RotationMatrix *yRot = new G4RotationMatrix; 
  yRot->rotateY(0 * rad);
  G4ThreeVector zTrans(0, 0, -sourceAl2_cutoutdepth + sourceAlo_z / 2);
  G4SubtractionSolid *SourceAl = new G4SubtractionSolid("SourceAl", SourceAl1, SourceAl2, yRot, zTrans);
  // G4UnionSolid *SourceAlUnion = new G4UnionSolid("SourceAlUnion", SourceAl, SourceAlo, 0, G4ThreeVector(0 * cm, 0 * cm, sourceAl_z / 2 + sourceAlo_z / 2 - tol));
  G4LogicalVolume *logicSourceAl = new G4LogicalVolume(SourceAl, Al, "SourceAl");
  G4ThreeVector sourceAlPos = G4ThreeVector(0, 0, -aScint2window_z - bScint2source_z - 4*(scint_z/2 + mylar_t + aluminum_t * 2) - sourceAlo_z - sourceAl_z / 2 - 2*tol);
  new G4PVPlacement(0, sourceAlPos, logicSourceAl, "SourceAl", logicWorld, false, 0, checkOverlaps);
  logicSourceAl->SetUserLimits(Limits);


  // mylar window
  G4Tubs *SourceMy = new G4Tubs("SourceMy", 0 * cm, source_r-tol, sourceMy_z / 2, 0, 360 * deg);
  G4LogicalVolume *logicSourceMy = new G4LogicalVolume(SourceMy, Mylar, "SourceMy");
  G4ThreeVector sourceMyPos = G4ThreeVector(0 * cm, 0 * cm, -aScint2window_z - bScint2source_z - 4*(scint_z/2 + mylar_t + aluminum_t * 2) - sourceAlo_z / 2);
  new G4PVPlacement(0, sourceMyPos, logicSourceMy, "SourceMy", logicWorld, false, 0, checkOverlaps);
  logicSourceMy->SetUserLimits(Limits);


  G4Tubs *SourceCal = new G4Tubs("SourceCal", 0 * cm, source_r-tol, sourceSr_z / 2, 0, 360 * deg);
  G4LogicalVolume *logicSourceCal= new G4LogicalVolume(SourceCal, Al, "SourceCal");
  // Placed in logicWorld (not logicSiPMLid) so GPS /gps/pos/confine can find it.
  // The subtraction solid parent prevents the navigator from descending into it.
  // Absolute position = sipmLid centre + (0,0,0) offset = -(32.079+8.89/2)*mm + 7.6477*mm
  // G4ThreeVector sourceCalPos = G4ThreeVector(0, 0, -aScint2window_z - bScint2source_z - 4*(scint_z/2 + mylar_t + aluminum_t * 2) - sourceAlo_z/2 + sourceSr_z / 2 + tol);
  G4ThreeVector sourceCalPos = G4ThreeVector(0, 0, -aScint2window_z - bScint2source_z - 4*(scint_z/2 + mylar_t + aluminum_t * 2) - sourceAlo_z + sourceSr_z / 2);
  new G4PVPlacement(0, sourceCalPos,
                    logicSourceCal, "SourceCal", logicWorld, false, 1, checkOverlaps);
  logicSourceCal->SetUserLimits(Limits); 

  

  // Remove geometry_export.gdml if it exists to avoid G4GDMLParser exception
  std::remove("../output/geometry_export.gdml");
  G4GDMLParser parser;
  parser.Write("../output/geometry_export.gdml", physWorld, true); // true = store auxiliary info

  return physWorld;
}

void BM_Detector::ConstructSDandField()
{
  /*
  NOTE: although this function is not explicitly called in BetaMon cc/hh files, it is still called.
    ConstructSDandField() is invoked in G4RunManager::InitializeGeometry() alongside Construct().
    If run with multithreading mode, it is invoked for each thread additionally from G4WorkerRunManager::InitializeGeometry().
  https://geant4-forum.web.cern.ch/t/constructsdandfield-in-multi-threaded-mode/2986
  */
  SDMan = G4SDManager::GetSDMpointer();

  G4VSensitiveDetector *bScintillatorSD = new BM_SD("BScintillatorSD", "BScintillatorHC"); // trigger
  G4VSensitiveDetector *aScintillatorSD = new BM_SD("AScintillatorSD", "AScintillatorHC");
  G4VSensitiveDetector *windowFoilSD = new BM_SD("WindowFoilSD", "WindowFoilHC"); // window foil
  G4VSensitiveDetector *vacuumSD = new BM_SD("VacuumSD", "VacuumHC");

  // Add the silicon detectors to the Sens.Det.Management
  SDMan->AddNewDetector(aScintillatorSD);
  SDMan->AddNewDetector(bScintillatorSD);
  SDMan->AddNewDetector(vacuumSD);
  SDMan->AddNewDetector(windowFoilSD);


  // Turn on the sensitive detectors (1 - window, 6 - vacuum, 3 - scint, 5 - trig(?) )
  vacuumWindowLV->SetSensitiveDetector(windowFoilSD);
  vacuumLV->SetSensitiveDetector(vacuumSD);
  aScintillatorLV->SetSensitiveDetector(aScintillatorSD);
  bScintillatorLV->SetSensitiveDetector(bScintillatorSD);

  // Magnetic field
  G4double amplitude = 0. * gauss;
  G4double theta = 60. * degree;
  G4MagneticField *MagneticField = new G4UniformMagField(G4ThreeVector(amplitude * sin(theta * 3.141592653 / 180), 
                                                          0., amplitude * cos(theta * 3.141592653 / 180)));
  G4FieldManager *globalFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  globalFieldMgr->SetDetectorField(MagneticField);
  globalFieldMgr->CreateChordFinder(MagneticField);
  vacuumLV->SetFieldManager(globalFieldMgr, false);

  fScoringVolume = logicWorld;
}