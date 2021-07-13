// simulation software for the Panda Barrel DIRC prototype
// contact: r.dzhygadlo@gsi.de

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "TROOT.h"

#include "PrtPhysicsList.h"
#include "PrtDetectorConstruction.h"

#include "PrtActionInitialization.h"
#include "time.h"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "TApplication.h"

#include "PrtRun.h"
#include "PrtManager.h"
#include "PrtLutReco.h"
#include "../../prttools/PrtTools.h"

namespace {
void PrintUsage() {
  G4cerr << " Usage: " << G4endl;
  G4cerr << " Prt [-m macro ] [-u UIsession] [-t nThreads] [-r seed] " << G4endl;
  G4cerr << "   note: -t option is available only for multi-threaded mode." << G4endl;
}
} // namespace

int main(int argc, char **argv) {
  for (G4int i = 1; i < argc; i = i + 2) std::cout << argv[i] << "  " << argv[i + 1] << std::endl;

  // Evaluate arguments
  if (argc > 100) {
    PrintUsage();
    return 1;
  }

#ifdef G4MULTITHREADED
  G4int nThreads = 1;
#endif
  TApplication theApp("App", 0, 0);

  G4String macro, events, geometry, radiator, physlist, session, geomTheta, geomPhi, batchmode,
    lensId, particle = "mix_pip", momentum, testVal1, testVal2, testVal3, prismStepX, prismStepY,
            beamZ, beamX, timeSigma, beamDimension, mcpLayout;
  TString infile = "", lutfile = "", pdffile = "", outfile = "";
  G4int firstevent(0), runtype(0), study(0), fid(0), verbose(0);

  G4long myseed = 0;
  for (G4int i = 1; i < argc; i = i + 2) {
    if (G4String(argv[i]) == "-m") macro = argv[i + 1];
    //    else if ( G4String(argv[i]) == "-u" ) session   = argv[i+1];
    else if (G4String(argv[i]) == "-seed") myseed = atol(argv[i + 1]);
    else if (G4String(argv[i]) == "-o") outfile = argv[i + 1];
    else if (G4String(argv[i]) == "-i") infile = argv[i + 1];
    else if (G4String(argv[i]) == "-u") lutfile = argv[i + 1];
    else if (G4String(argv[i]) == "-pdf") pdffile = argv[i + 1];
    else if (G4String(argv[i]) == "-g") geometry = argv[i + 1];
    else if (G4String(argv[i]) == "-h") radiator = argv[i + 1];
    else if (G4String(argv[i]) == "-a") geomTheta = argv[i + 1];
    else if (G4String(argv[i]) == "-phi") geomPhi = argv[i + 1];
    else if (G4String(argv[i]) == "-b") batchmode = argv[i + 1];
    else if (G4String(argv[i]) == "-f") firstevent = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-e") events = argv[i + 1];
    else if (G4String(argv[i]) == "-l") lensId = argv[i + 1];
    else if (G4String(argv[i]) == "-x") particle = argv[i + 1];
    else if (G4String(argv[i]) == "-p") momentum = argv[i + 1];
    else if (G4String(argv[i]) == "-w") physlist = argv[i + 1];
    else if (G4String(argv[i]) == "-r") runtype = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-study") study = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-fid") fid = atoi(argv[i + 1]);
    else if (G4String(argv[i]) == "-z") beamDimension = argv[i + 1];
    else if (G4String(argv[i]) == "-c") mcpLayout = argv[i + 1];
    else if (G4String(argv[i]) == "-t1") testVal1 = argv[i + 1];
    else if (G4String(argv[i]) == "-t2") testVal2 = argv[i + 1];
    else if (G4String(argv[i]) == "-t3") testVal3 = argv[i + 1];
    else if (G4String(argv[i]) == "-gsx") prismStepX = argv[i + 1];
    else if (G4String(argv[i]) == "-gsy") prismStepY = argv[i + 1];
    else if (G4String(argv[i]) == "-gz") beamZ = argv[i + 1];
    else if (G4String(argv[i]) == "-gx") beamX = argv[i + 1];
    else if (G4String(argv[i]) == "-tr") timeSigma = argv[i + 1];
    else if (G4String(argv[i]) == "-v") verbose = atoi(argv[i + 1]);
#ifdef G4MULTITHREADED
    else if (G4String(argv[i]) == "-t") {
      nThreads = G4UIcommand::ConvertToInt(argv[i + 1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }

  if (runtype == 1) {
    particle = "opticalphoton";
    momentum = " 3.18e-09";
    geomTheta = "180";
    geomPhi = "0";
  }

  if (outfile == "" && runtype == 0 || runtype == 5) outfile = "hits.root"; // simulation
  if (outfile == "" && (runtype == 1 || runtype == 7 || runtype == 11))
    outfile = "../data/lut.root";                                 // lookup table generation
  if (outfile == "" && runtype == 6) outfile = "focalplane.root"; // focal plane simulation

  if (batchmode.size()) gROOT->SetBatch(kTRUE);
  if (!events.size()) events = "1";

  PrtTools t;
  PrtRun *run = t.find_run(study, fid);

  if (runtype == 2 || runtype == 3 || runtype == 4) {
    if (infile == "") {
      infile = t.get_inpath();
      std::cout << "--- infile  " << infile << std::endl;
    }
    run = t.get_run(infile);
    if (outfile == "") {
      outfile = t.get_outpath();
      if (run->getStudy() == 0) outfile = "reco.root";
      std::cout << "--- outfile  " << outfile << std::endl;
    }
  }

  run->setRunType(runtype);

  if (momentum.size()) run->setMomentum(atof(momentum));
  if (physlist.size()) run->setPhysList(atoi(physlist));
  if (geometry.size()) run->setGeometry(atoi(geometry));
  if (radiator.size()) run->setRadiator(atoi(radiator));
  if (lensId.size()) run->setLens(atoi(lensId));
  if (mcpLayout.size()) run->setPmtLayout(atoi(mcpLayout));
  if (beamDimension.size()) run->setBeamSize(atof(beamDimension));
  if (testVal1.size()) run->setTest1(atof(testVal1));
  if (testVal2.size()) run->setTest2(atof(testVal2));
  if (testVal3.size()) run->setTest3(atof(testVal3));
  if (geomTheta.size()) run->setTheta(atof(geomTheta));
  if (geomPhi.size()) run->setPhi(atof(geomPhi));
  if (prismStepX.size()) run->setPrismStepX(atof(prismStepX));
  if (prismStepY.size()) run->setPrismStepY(atof(prismStepY));
  if (beamX.size()) run->setBeamX(atof(beamX));
  if (beamZ.size()) run->setBeamZ(atof(beamZ));
  if (timeSigma.size()) run->setTimeSigma(atof(timeSigma));

  PrtManager::Instance(outfile, run);

  std::cout << "=== Run info:  " << std::endl << run->getInfo() << std::endl;

  if (runtype == 2 || runtype == 3 || runtype == 4) {
    PrtLutReco *reco = new PrtLutReco(infile, lutfile, pdffile, verbose);
    reco->Run(firstevent, atoi(events));
    return 0;
  }

  G4Random::setTheEngine(new CLHEP::RanecuEngine);

#ifdef G4MULTITHREADED
  G4MTRunManager *runManager = new G4MTRunManager;
  if (nThreads > 0) runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager *runManager = new G4RunManager;
#endif

  G4Random::setTheSeed(myseed);

  runManager->SetUserInitialization(new PrtDetectorConstruction());
  runManager->SetUserInitialization(new PrtPhysicsList());
  runManager->SetUserInitialization(new PrtActionInitialization());
  runManager->Initialize();

  G4VisManager *visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager *UImanager = G4UImanager::GetUIpointer();

  if (macro.size()) {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + macro);
  } else {
    //  UImanager->ApplyCommand("/control/execute ../prt.mac");

    UImanager->ApplyCommand("/gun/direction 0 0 1");
    UImanager->ApplyCommand("/gun/position 0 0 -100 cm");
    UImanager->ApplyCommand("/Prt/geom/prtRotation 90 deg");
    UImanager->ApplyCommand("/gun/direction 0 0 1");
  }

  if (particle.size()) {
    int pdgid = 0;
    if (particle == "mix_pip") PrtManager::Instance()->getRun()->setPid(10001);
    else if (particle == "mix_pik") PrtManager::Instance()->getRun()->setPid(10002);
    else {
      G4String command = "/gun/particle ";
      UImanager->ApplyCommand(command + particle);
      if (particle == "proton") pdgid = 2212;
      if (particle == "pi+") pdgid = 211;
      if (particle == "pi0") pdgid = 111;
      if (particle == "kaon+") pdgid = 321;
      if (particle == "kaon-") pdgid = -321;
      if (particle == "mu-") pdgid = 13;
      if (particle == "e-") pdgid = 11;
      if (particle == "opticalphoton") pdgid = 0;
      PrtManager::Instance()->getRun()->setPid(pdgid);
    }
  }

  // if(momentum.size()) UImanager->ApplyCommand( "/gun/momentumAmp "+momentum);

  if (batchmode.size()) { // batch mode
    UImanager->ApplyCommand("/run/beamOn " + events);
  } else { // UI session for interactive mode

    G4UIExecutive *ui = new G4UIExecutive(argc, argv, "Qt");
    UImanager->ApplyCommand("/control/execute ../vis.mac");
    if (ui->IsGUI()) UImanager->ApplyCommand("/control/execute gui.mac");
    UImanager->ApplyCommand("/run/beamOn " + events);
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;

  return 0;
}
