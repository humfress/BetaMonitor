// BM_Detector.cc - GDML-driven detector construction for BetaMon

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4SDManager.hh"
#include "G4TransportationManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4GDMLParser.hh"
#include "G4GenericMessenger.hh"

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "BM_Detector.hh"
#include "BM_SD.hh"
#include "BM_SteppingAction.hh"

#include <cstdlib>

G4FieldManager *BM_Detector::fFieldMgr = 0;
G4UniformMagField *BM_Detector::fMagneticField = 0;

BM_Detector::BM_Detector() : G4VUserDetectorConstruction(), vacuumLV(nullptr), 
vacuumWindowLV(nullptr), aScintillatorLV(nullptr), bScintillatorLV(nullptr), logicWorld(nullptr),
fGdmlFile("../output/geometry_export.gdml")
{
  fMessenger = std::make_unique<G4GenericMessenger>(this, "/bm/det/", "Detector controls");
  fMessenger->DeclareProperty("gdmlFile", fGdmlFile,
                              "Set GDML geometry path (set before /run/initialize).");
}

BM_Detector::~BM_Detector() {}

G4VPhysicalVolume *BM_Detector::Construct()
{
  G4cout << "Loading detector geometry from GDML." << G4endl;

  static const G4String kDefaultGdmlPath = "../output/geometry_export.gdml";
  const char *gdmlPathEnv = std::getenv("BM_GDML_FILE");
  G4String gdmlPath = fGdmlFile.empty() ? kDefaultGdmlPath : fGdmlFile;
  if (gdmlPath == kDefaultGdmlPath && gdmlPathEnv && gdmlPathEnv[0] != '\0')
  {
    gdmlPath = G4String(gdmlPathEnv);
  }

  G4cout << "GDML input file: " << gdmlPath << G4endl;

  G4GDMLParser parser;
  parser.SetOverlapCheck(true);
  parser.Read(gdmlPath, false);

  G4VPhysicalVolume *physWorld = parser.GetWorldVolume();
  if (!physWorld)
  {
    G4Exception("BM_Detector::Construct", "BMDET001", FatalException,
                "Failed to read GDML world volume.");
  }

  logicWorld = physWorld->GetLogicalVolume();

  auto lvStore = G4LogicalVolumeStore::GetInstance();
  vacuumLV = lvStore->GetVolume("Vacuum", false);
  vacuumWindowLV = lvStore->GetVolume("VacuumWindow", false);
  aScintillatorLV = lvStore->GetVolume("AScintillator", false);
  bScintillatorLV = lvStore->GetVolume("BScintillatorLV", false);
  if (!bScintillatorLV)
  {
    bScintillatorLV = lvStore->GetVolume("BScintillator", false);
  }

  if (!vacuumLV || !vacuumWindowLV || !aScintillatorLV || !bScintillatorLV)
  {
    G4Exception("BM_Detector::Construct", "BMDET002", FatalException,
                "Missing one or more expected logical volumes in GDML (Vacuum, VacuumWindow, AScintillator, BScintillatorLV). ");
  }

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