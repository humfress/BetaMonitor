// BM_Telescope.cc - Implementation of BM_Telescope class: BetaMon telescope and geometry

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

#include "BM_Telescope.hh"
#include "BM_SD.hh"
#include "BM_SteppingAction.hh"


// MagField* BM_Telescope::fMagneticField = 0;
G4FieldManager *BM_Telescope::fFieldMgr = 0;
G4UniformMagField *BM_Telescope::fMagneticField = 0;

BM_Telescope::BM_Telescope() : G4VUserDetectorConstruction(), vacuumLV(nullptr), 
vacuumWindowLV(nullptr), aScintillatorLV(nullptr), bScintillatorLV(nullptr), logicWorld(nullptr) {}

BM_Telescope::~BM_Telescope() {}

G4VPhysicalVolume *BM_Telescope::Construct()
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
  G4Material *Glass = nist->FindOrBuildMaterial("G4_GLASS_PLATE");
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
  G4double aScint_z = 2 * mm; // scintillator thickness
  G4double aScint_d = 50 * mm; // scintillator A diameter
  G4double bScint_z = 50 * mm; // scintillator thickness
  G4double bScint_d = 50 * mm; // scintillator B diameter
  G4double aScint2window_z = 8.4074 * mm; // distance from scintillatorA to window 
  G4double bScint2window_z = aScint2window_z + aScint_z + 3.0*mm; // distance from scintillatorB to window
  G4double bScint2source_z = 2.642 * cm; // distacne from scintillatorB to calibration source
  G4double Tdv_r1o = 3.81 / 2 * cm;       // T Pipe outer radius
  G4double cyl_r2dvo = 6.9088 / 2 * cm; // Flange Outer Radius
  G4double cyl_r1dvi = 3.4798 / 2 * cm; // Inner Radius of T Pipe Decay Volume 
  G4double cyl_hc = 0.2 * cm;     // thickness of glass window
  G4double window_d = 20 *mm; // diameter of window
  G4double mylar_t = 0.076 * mm; // scint mylar thickness
  G4double aluminum_t = 0.0001 * mm; // scint aluminum thickness 
  G4double tpipe_dz = flange_width;// used to slide whole t-pipe and vacuum assembly along the z axis. 
 


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
  new G4PVPlacement(yRotT, G4ThreeVector(0 * cm, 0 * cm, tpipe_dz + Tdv_h2 + flange_width),
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
  new G4PVPlacement(yRotT, G4ThreeVector(0, 0, tpipe_dz + Tdv_h2 + flange_width),
                    steelTPipeLV, "SteelTPipe", logicWorld, false, 0, checkOverlaps);
  steelTPipeLV->SetUserLimits(Limits);


  // === T Pipe Flanges (steel) === 

  G4Tubs *Flange = new G4Tubs("Flange", cyl_r1dvi, cyl_r2dvo, flange_width / 2, 0, 360 * deg);
  G4Tubs *FlangeScint = new G4Tubs("FlangeScint", window_d/2, cyl_r2dvo, flange_width / 2, 0, 360 * deg);
  // leftmost flange (nearest beta monitor)
  G4LogicalVolume *tpipeFlangeScintLV = new G4LogicalVolume(FlangeScint, Stainless_Steel, "tPipeFlangeScint");
  new G4PVPlacement(0, G4ThreeVector(0 *  cm, 0 * cm, flange_width / 2), 
                   tpipeFlangeScintLV, "tPipeFlangeScint", logicWorld,
                    false, 0, checkOverlaps);
  tpipeFlangeScintLV->SetUserLimits(Limits);

  // second to leftmost flange (sandwiches the "window")
  G4LogicalVolume *tPipeFlange1LV = new G4LogicalVolume(Flange, Stainless_Steel, "tPipeFlange1");
  new G4PVPlacement(0, G4ThreeVector(0 * cm, 0 * cm, tpipe_dz + flange_width / 2), 
                    tPipeFlange1LV, "tPipeFlange1", logicWorld, false, 0, checkOverlaps);
  tPipeFlange1LV->SetUserLimits(Limits);

  // right flange
  G4LogicalVolume *tPipeFlange2LV = new G4LogicalVolume(Flange, Stainless_Steel, "tPipeFlange2");
  new G4PVPlacement(yRotT, G4ThreeVector(Tdv_h1/2 + flange_width / 2, 0, tpipe_dz + cyl_hdv / 2 + flange_width), tPipeFlange2LV, "tPipeFlange2",
                    logicWorld, false, 0, checkOverlaps);
  tPipeFlange2LV->SetUserLimits(Limits);

  // middle flange
  G4LogicalVolume *tPipeFlange3LV = new G4LogicalVolume(Flange, Stainless_Steel, "tPipeFlange3");
  new G4PVPlacement(yRotT, G4ThreeVector(-Tdv_h1/2 - flange_width / 2, 0, tpipe_dz + cyl_hdv / 2 + flange_width),
                    tPipeFlange3LV, "tPipeFlange3", logicWorld, false, 0, checkOverlaps);
  tPipeFlange3LV->SetUserLimits(Limits);


  // === Vacuum Window (copper) - between two left flanges ===
  // G4double cyl_r2c = 3.556 / 2 * cm; // copper seal radius //Tdv_r1i = 3.4798 / 2 * cm;
  // G4double cyl_hkap = 0.0012 * cm;   // thickness of kapton (trials - 1:0.0012, 2: 0.00075, 3: 0.006, 4: 0.0127)
  G4Tubs *vacuumWindowDisk = new G4Tubs("VacuumWindowDisk", 0. * cm, window_d/2, cyl_hc / 2., 0, 360 * deg);
  vacuumWindowLV = new G4LogicalVolume(vacuumWindowDisk, Glass, "VacuumWindow");
  new G4PVPlacement(0, G4ThreeVector(0 * cm, 0 * cm, cyl_hc / 2), vacuumWindowLV, "VacuumWindow", logicWorld, false, 2, checkOverlaps);
  vacuumWindowLV->SetUserLimits(Limits);

  // === Scintillators (A & B) ===

  // // bare scintillator
  // G4Box *SmallScin = new G4Box("SmallScin", scint_xy / 2, scint_xy / 2, scint_z / 2);
  // thin scintillator A, tubular
  G4Tubs *ThinScinA = new G4Tubs("ThinScinA", 0, aScint_d / 2, aScint_z / 2, 0, 360 * deg);
  G4Tubs *ThinScinA_Al = new G4Tubs("ThinScinA_Al", 0, aScint_d / 2 + aluminum_t, aScint_z / 2 + aluminum_t, 0, 360 * deg);
  
  // inner scintillator "A" (inner = closer to window)
  // 8.4074 mm is the distnce between the window and the surface of the scint. 
  G4ThreeVector aScintillatorPos = G4ThreeVector(0 * cm, 0 * cm, - aScint2window_z - (aScint_z / 2  + aluminum_t));
  aScintillatorLV = new G4LogicalVolume(ThinScinA, PVT, "AScintillator");
  new G4PVPlacement(0, aScintillatorPos, aScintillatorLV, "AScintillator", logicWorld, false, 3, checkOverlaps);
  aScintillatorLV->SetUserLimits(Limits);
  
  // inner aluminum layer for scint A
  G4SubtractionSolid *solidAlMyI = new G4SubtractionSolid("InnerAl", ThinScinA_Al, ThinScinA);
  G4LogicalVolume *logicAlMylarAli1 = new G4LogicalVolume(solidAlMyI, Al, "Ali_sq1");
  new G4PVPlacement(0, aScintillatorPos, logicAlMylarAli1, "Ali_sq1", logicWorld, false, 0, checkOverlaps);
  logicAlMylarAli1->SetUserLimits(Limits);

  G4Tubs *ThicScinB = new G4Tubs("ThicScinB", 0, bScint_d / 2, bScint_z / 2, 0, 360 * deg);
  G4Tubs *ThicScinB_Al = new G4Tubs("ThicScinB_Al", 0, bScint_d / 2 + aluminum_t, bScint_z / 2 + aluminum_t, 0, 360 * deg);
  // the scintillators are wrapped in aluminized mylar (a sandwich of Al - My - Al)
                                                            


  
  // outer scintillator "B" (outer = farther from window)
  // G4ThreeVector bScintillatorPos = G4ThreeVector(0 * cm, 0 * cm, -(3.317 * mm + (3.0 + 0.0762) * 3 / 2 * mm));

  // G4ThreeVector bScintillatorPos = G4ThreeVector(0 * cm, 0 * cm,  - aScint2window_z - 4 / 2 * (scint_z + 2*mylar_t + 4*aluminum_t));
  G4ThreeVector bScintillatorPos = G4ThreeVector(0 * cm, 0 * cm, - bScint2window_z - (bScint_z / 2  + aluminum_t * 2));
  bScintillatorLV = new G4LogicalVolume(ThicScinB, PVT, "BScintillatorLV");
  new G4PVPlacement(0, bScintillatorPos, bScintillatorLV, "BScintillator", logicWorld, false, 4, checkOverlaps);
  bScintillatorLV->SetUserLimits(Limits);

  // aluminum layer for scint B
  G4SubtractionSolid *solidAlMyI2 = new G4SubtractionSolid("InnerAl2", ThicScinB_Al, ThicScinB);
  G4LogicalVolume *logicAlMylarAli2 = new G4LogicalVolume(solidAlMyI2, Al, "Ali_sq2");
  new G4PVPlacement(0, bScintillatorPos, logicAlMylarAli2, "Ali_sq2", logicWorld, false, 0, checkOverlaps);
  logicAlMylarAli2->SetUserLimits(Limits);          


  // Remove geometry_telescope_export.gdml if it exists to avoid G4GDMLParser exception
  std::remove("../output/geometry_telescope_export.gdml");
  G4GDMLParser parser;
  parser.Write("../output/geometry_telescope_export.gdml", physWorld, true); // true = store auxiliary info

  return physWorld;
}

void BM_Telescope::ConstructSDandField()
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