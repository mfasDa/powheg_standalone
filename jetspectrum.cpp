#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/PowhegHooks.h"

#include "TH1.h"
using namespace Pythia8

int main(int argc, const char **argv){
  // Generator
  Pythia pythia;

  // Load configuration file
  pythia.readFile("main31.cmnd");

  // Read in main settings
  int nEvent      = pythia.settings.mode("Main:numberOfEvents");
  int nError      = pythia.settings.mode("Main:timesAllowErrors");
  // Read in key POWHEG merging settings
  int vetoMode    = pythia.settings.mode("POWHEG:veto");
  int MPIvetoMode = pythia.settings.mode("POWHEG:MPIveto");
  bool loadHooks  = (vetoMode > 0 || MPIvetoMode > 0);

  // Add in user hooks for shower vetoing
  PowhegHooks *powhegHooks = NULL;
  if (loadHooks) {

    // Set ISR and FSR to start at the kinematical limit
    if (vetoMode > 0) {
      pythia.readString("SpaceShower:pTmaxMatch = 2");
      pythia.readString("TimeShower:pTmaxMatch = 2");
    }

    // Set MPI to start at the kinematical limit
    if (MPIvetoMode > 0) {
      pythia.readString("MultipartonInteractions:pTmaxMatch = 2");
    }

    powhegHooks = new PowhegHooks();
    pythia.setUserHooksPtr((UserHooks *) powhegHooks);
  }

  // Initialise and list settings
  pythia.init();
  while (true) {

    // Generate the next event
    if (!pythia.next()) {

      // If failure because reached end of file then exit event loop
      if (pythia.info.atEndOfFile()) break;

      // Otherwise count event failure and continue/exit as necessary
      cout << "Warning: event " << iEvent << " failed" << endl;
      if (++iError == nError) {
        cout << "Error: too many event failures.. exiting" << endl;
        break;
      }

      continue;
    }

    // logics here

  }
}