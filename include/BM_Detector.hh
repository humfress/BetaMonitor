// BM_Detector.hh - header for BM_Detector class: BetaMon detector and geometry

#ifndef BM_DETECTOR_H
#define BM_DETECTOR_H

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"

#include "G4UniformMagField.hh"
#include "G4SDManager.hh"

#include <memory>

class G4VPhysicalVolume;
class G4LogicalVolume;
class MagField;
class G4GenericMessenger;

/// Detector construction class to define materials and geometry.

class BM_Detector : public G4VUserDetectorConstruction
{
public:
  BM_Detector();
  virtual ~BM_Detector();

  virtual G4VPhysicalVolume *Construct();
  virtual void ConstructSDandField();

  G4LogicalVolume *GetScoringVolume() const { return fScoringVolume; }
  G4SDManager *SDMan;

  G4LogicalVolume *vacuumLV;
  G4LogicalVolume *vacuumWindowLV;
  G4LogicalVolume *aScintillatorLV;
  G4LogicalVolume *bScintillatorLV;
  G4LogicalVolume *logicWorld;

  static G4UniformMagField *fMagneticField;
  static G4FieldManager *fFieldMgr;

protected:
  G4LogicalVolume *fScoringVolume;
  G4String fGdmlFile;
  std::unique_ptr<G4GenericMessenger> fMessenger;
};

#endif
