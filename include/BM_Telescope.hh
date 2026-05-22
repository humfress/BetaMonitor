// BM_Telescope.hh - header for BM_Telescope class: BetaMon telescope and geometry

#ifndef BM_TELESCOPE_H
#define BM_TELESCOPE_H

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"

#include "G4UniformMagField.hh"
#include "G4SDManager.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class MagField;

/// Detector construction class to define materials and geometry.

class BM_Telescope : public G4VUserDetectorConstruction
{
public:
  BM_Telescope();
  virtual ~BM_Telescope();

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
};

#endif
