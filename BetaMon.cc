#include <fstream>
#include <string>
#include <chrono>

// #define SEG_DBG 
#ifdef SEG_DBG
  #include <stdio.h>
  #include <execinfo.h>
  #include <signal.h>
  #include <stdlib.h>
  #include <unistd.h>
#endif

#include "G4RootAnalysisManager.hh"
#include "G4Types.hh"
#include "G4UImanager.hh"

#define G4VIZ
#ifdef G4VIZ
  #include "G4VisExecutive.hh"
  #include "G4UIExecutive.hh"
#endif
#include "G4RunManager.hh"

#include "BM_RunAction.hh"
#include "BM_EventAction.hh"
#include "BM_SteppingAction.hh"
//#include "BM_PhysicsList.hh" // removed
// #include "QBBC.hh"  // a standard phys list
#include "G4PhysListFactory.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4DecayPhysics.hh"
#include "BM_Telescope.hh"
#include "BM_PrimaryGenerator.hh"
#include "BM_Output.hh"

#ifdef SEG_DBG
void handler(int sig) // for segfault
{
   void* array[10];
   size_t size;
   size = backtrace(array, 10);
   fprintf(stderr, "Error: signal-%d:\n", sig);
   backtrace_symbols_fd(array, size, STDERR_FILENO);
   exit(1);
}
#endif

int main(int argc, char** argv)
{ 
  // NOTE: this version of the BetaMon sim was changed to single-threaded
  // operation.  (change to main branch to find David McClain's multi-threaded version).

  auto now = std::chrono::system_clock::now();
  auto ts_start = std::chrono::system_clock::to_time_t(now);
  G4cout << "Sim started: " << std::ctime(&ts_start) << G4endl;
  int unix_start = static_cast<int>(ts_start);

  // review input arguments
  G4cout << "Running with " << argc << " arguments:";
  for (int i = 0; i < argc; ++i) {
    G4cout << i << ": " << argv[i] << "  ";
  }
  G4cout << G4endl;
    
  #ifdef SEG_DBG
    std::cout << "Using our segfault debugger...\n";
    signal(SIGSEGV, handler); // using the segfault handler
  #endif

  G4RunManager* runManager = new G4RunManager;
  runManager->SetVerboseLevel(1);

  // initialize detector
  runManager->SetUserInitialization(new BM_Telescope());
  
  // initialize physics list - use a standard one for now
  G4PhysListFactory factory;
  // G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("QGSP_BIC_HP");
  G4VModularPhysicsList* physicsList = factory.GetReferencePhysList("QBBC"); 
  physicsList->RegisterPhysics(new G4RadioactiveDecayPhysics());


  physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);
  
  // set user actions
  runManager->SetUserAction(new BM_PrimaryGenerator());
  runManager->SetUserAction(new BM_EventAction());
  runManager->SetUserAction(new BM_SteppingAction());
  runManager->SetUserAction(new BM_RunAction());
  
  // initialize ROOT output (example: ./BetaMon run1.mac ./output/outfile.root)
  G4String outname = "./output/test.root";
  if (argc == 3) {
    outname = argv[2];
  }
  G4cout << "Output file: " << outname << G4endl;
  BM_Output::Instance()->SetFilename(outname);
  
  // initialize visualization
  #ifdef G4VIZ
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
  #endif

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  if (argc != 1)
  {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }
  else
  {
    // interactive mode : define UI session
    G4UIExecutive* ui = 0;
    ui = new G4UIExecutive(argc, argv);
    UImanager->ApplyCommand("/control/macroPath ../scripts/macros/");
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
    // #endif
  }

  // Job termination
  // Free the store. user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
  #ifdef G4VIZ
    delete visManager;
  #endif
  
  delete runManager;

  now = std::chrono::system_clock::now();
  auto ts_stop = std::chrono::system_clock::to_time_t(now);
  G4cout << "Sim complete: " << std::ctime(&ts_stop) << G4endl;
  int unix_stop = static_cast<int>(ts_stop);

  G4cout << "Time elapsed: " << unix_stop - unix_start << " seconds.\n";

  return 0;
}
