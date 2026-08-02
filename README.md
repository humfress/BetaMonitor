# BetaMonitor Usage
This GEANT4 simulation is designed to populate a vacuum cell with decay products from an input source, and measure the energy deposited within a detector located near the cell's window. The goal is to measure the rate of detection from the source intensity.

## Geant4 Details
This implementation of BetaMonitor uses Geant4 11.4.0. 

Required G4 installation used:

- GEANT4_INSTALL_DATA: ON
    - required for use of G4RadioactiveDecayModule
- GEANT4_USE_GDML: ON
    - Required for use of GDML export function of geometry.
- If using G4 GUI:
    - GEANT4_USE_QT: ON
    - GEANT4_USE_RAYTRACER_X11: ON
- Options used for local implementation (should not be critical for BetaMonitor):
    - GEANT4_USE_SYSTEM_EXPAT: ON
    - GEANT4_USE_SYSTEM_ZLIB: ON
    - GEANT4_BUILD_MULTITHREADED: ON
        - Current implementation is single thread, so this is not required. 


## Set up & Compilation
To run this simulation create a ```build/``` directory and ```cmake``` pointing toward the BetaMonitor directory. Following this utilize the ```make``` command within the build directory whenever changes are made to the files in the ```include/``` and ```src/``` directories. There will then be an executable created named BetaMon within the build directory. 

In the build directory, running the following command will launch the BetaMonitor program in the G4 GUI. 

```./BetaMon```

## Running Simulations via Macro. 
The details of the simulation are controlled via macro. In /BetaMonitor/Scripts/Macros there are several pre-configured macros for simple simulations. An example of the macro file for 6He is repeated below:


/run/initialize 

This time threshold is important for capturing decay products that are appropriate for the timescale of the experiement For a given species, look at the decay chain and confirm if this value needs to be adjusted to something shorter (ex: for 241Am, a value of ~1e4 years will prevent an unstable daughter nuclear with T_1/2 of 1e5 years from dominating the decay spectrum)

/process/had/rdm/thresholdForVeryLongDecayTime 1.0e+60 year

The General Particle Source (GPS) module is used to initialize single ions with whatever parameters are desired. In this case, they are created with 0 initial energy, so that all subsequent energy comes purely from the decay processed. 

/gps/particle ion
/gps/ene/mono 0 eV
/gps/ion 2 6 0 0  # Z=2 (He), A=6, Q=0, E=0
/gps/number 1

GPS is also used to confine the particles to a particle volume. First they are randomly generated in a parallelepiped with halfx, halfy, and halfz values that encompass the volume of interest. Then only particles in the volume name ("Vacuum" in the example below) satisfy the conditions of the simulation and are counted as a run. Making the generation volume too much larger than the confinement volume can slow down the particle generation process. 

/gps/pos/type Volume
/gps/pos/shape Para
/gps/pos/halfx 20 cm
/gps/pos/halfy 2 cm
/gps/pos/halfz 10 cm
/gps/pos/centre 0 0 8 cm
/gps/pos/confine Vacuum
/gps/ang/type iso

Finally, this sets the number of primaries that are generated in total. 

/run/beamOn 1000000

To run a particle macro from the build directory, a command like the one below will direct the output to a .root file with the specified name: 

./BetaMon ../scripts/macros/run_source_6He.mac ../output/6He_1e6Primaries_output.root

## Merging files
Included in the ```BetaMonitorAnalysis.ipynb``` is the ```load2root``` function which utilizes the uproot library to load in the series of files associated with one simulation and condenses them, after making the cut that the main detector (archaically called SQ) must have some energy deposited by the event. The new root file will contain all events that deposit energy in the detector in the ```[Tree]``` branch, and the total number of simulated events in the ```[Size]``` branch. These can be read in individually via uproot. After condensing these files, one can delete the original separated files to free up disk space.

# Changing Sources
To change the source, GPS settings should be adjusted in a macro file. Two key attributes must be addressed: the source ion itself, and the geometry/volume of where it is generated:

For the source ion, the user must set the Z and A values of the desired isotope using \gps\ion Z A Q E (ex: for 6He, \gps\ion 2 6 0 0)

For the geometry, the user must set a bounding volume based on the overall BetaMonitor geometry (relative to the world origin), and then pick a confinement volume. For primaries generated in the T-Pipe, use the "Vacuum" volume. For calibration sources slotted into the back window, there is a very thin volume called SourceCal. See the included macro scripts for examples. 

Because the primary particle generation process is using the radioactive decay module (RDM), there is no need for the user to specify the particle type of the decay products (e.g., electron or positron). All decay products will be handled and simulated by RDM. 

# Detector & Decay Volume Geometry
All physical objects are defined within ```src/BM_Detector.cc```. Adding new detectors will require the creation of a new ```logicDetector``` within ```include/BM_Detector.hh```. The detector will then be assigned to a sensitive detector. To properly associate the events with the detector, adjust the HitsCollection within ```include/BM_EventAction.hh```, and ```src/BM_EventAction.cc```. If a new designation is required beyond the provided ```trigger, square, window, vac```, a new case statement will need to be added to ```src/BM_EventAction.cc```.

## Generating GDML Geometry

From the repository root, activate the Conda environment containing `pyg4ometry` and run the generator:

```bash
conda activate g4-env-py311
python scripts/python/gdml_create_betamon.py
```

This writes `dat/geometry_export_new.gdml`. Pass a path to write the geometry elsewhere:

```bash
python scripts/python/gdml_create_betamon.py dat/geometry_export.gdml
```

# BetaMonitor Implementation Details
## Data Flow
```text
BetaMon.cc (main)

  |
  +-- /run/initialize -> BM_Detector::Construct()   [build geometry, register SDs]
  |
  +-- /run/beamOn N  -> for each event:
    BM_PrimaryGenerator::GeneratePrimaries()  [fire particle]
    |
    +-- particle steps through geometry
      BM_SteppingAction::UserSteppingAction()  [intercept each step]
      BM_SD::ProcessHits()                     [record hits per volume]
      |
      +-- BM_Hit (data object)
    |
    BM_EventAction::BeginOfEventAction()   [look up HC IDs]
    BM_EventAction::EndOfEventAction()     [harvest all HCs -> BM_Output]
      |
      +-- BM_Output::Fill() / FillPrimaryInput()  [write ROOT TTrees]
    |
  BM_RunAction::BeginOfRunAction()    [OpenFile]
  BM_RunAction::EndOfRunAction()      [CloseFile]
```

  ## Files

  ### BetaMon.cc
  Sets up program and runtime. Configures CLI arguments. Created a G4RunManager. 

  Registers physics lists. Current implementation uses QBBC. Adds radioactive decay physics. 

  Registers all user actions (BM_PrimaryGenerator, BM_EventAction, BM_SteppingAction, BM_RunAction).

  ### BM_Detector.cc (G4VUserDetectorConstruction)
  Builds physical geometry and registers sensitive detectors. 

  Sensitive Detectors are numbered in accordance with their likely order of interaction for a primary particle generated in the T-Pipe:
  - Vacuum (1)
  - Window (2)
  - Scintillator A (3)
  - Scintillator B (4)

  ### BM_PrimaryGenerator.cc (G4VUserPrimaryGeneratorAction)
  Sets up the primary particle generation using GPS. All particle details are handled in the macro, not here. 

  ### BM_SD.cc (G4VSensitiveDetector)
  For every step that happens in a sensitive volume, BM_SD will create a "hit" and collect hits into a "Hits Collection." For each primary event, each volume will have its own BM_HitsCollection.

  ### BM_Hit.cc (G4VHit)
  Defines variables and data extraction methods for "hit" class used in BM_SD.

  ### BM_EventAction.cc (G4UserEventAction)
  After a primary event, BM_EventAction gathers "track-wise" information about energy deposition in the sensitive volumes and writes this data to the ROOT file. 

  ### BM_SteppingAction.cc (G4UserSteppingAction)
  Intercepts all steps in all volumes. Main use is to log Primary Particle information for any particles created by RDM. There is also some special handling of sulfur and chlorine ions. 

  ### BM_Output.cc
Configures the output ROOT file. The ROOT file has two trees:
  - simData: one row per track per sensitive volume with the following columns:
    - pid: particle ID
    - eventid: global eventid, unique for each generated primary ion event. 
    - trackid: for each event, tracks are numbers for individual particles that are tracked through the simulation. 
    - parentid: for secondary particles, the trackID of the parent particle. 
    - volumeid: volume ID of the 4 sensitive volumes. 
    - primaryenergy: energy of the particle at the track vertex. 
    - inenergy: energy of the tracked particle when it enters the sensitive volume.
    - kineticenergy: energy of the tracked volume at the last tracked step in the sensitive volume.
    - depenergy: deposited energy in the sensitive volume of the track.
    - exited: boolean flag for if the particle exited the sensitive volume. 
    - x, y, z: global position of the particle's entry point into the sensitive volume. 
    - px, py, pz: particle entry momentum into the sensitive volume.
    - pathlen: total path length of the track of the particle in the sensitive volume. 
- primaryInput: one row per radioactive decay product as captured in BM_SteppingAction.cc
    - pid: particle ID
    - eventid: global eventid, unique for each generated primary ion event. Can be used to coincide events between the two trees. 
    - x, y, z: global position of where the particles were generated. 
    - px, py, pz: Initial momentum of all primaries. 
    - primaryenergy: initial kinetic energy of all primaries. 

  ### BM_RunAction.cc (G4UserRunAction)
  Controls run. Opens (and later closes) ROOT file. 