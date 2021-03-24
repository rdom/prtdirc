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

//#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
//#endif

//#ifdef G4UI_USE
#include "G4UIExecutive.hh"
//#endif

#include "TApplication.h"

#include "PrtManager.h"
#include "PrtLutReco.h"

namespace {
  void PrintUsage() {
    G4cerr<<" Usage: "<<G4endl;
    G4cerr<<" Prt [-m macro ] [-u UIsession] [-t nThreads] [-r seed] "<<G4endl;
    G4cerr<<"   note: -t option is available only for multi-threaded mode."<<G4endl;
  }
}

int main(int argc,char** argv)
{
  for ( G4int i=1; i<argc; i=i+2 ) std::cout<< argv[i] << "  "<<argv[i+1] <<std::endl;
  
  // Evaluate arguments
  if ( argc > 100 ) {
    PrintUsage();
    return 1;
  }
 
#ifdef G4MULTITHREADED
  G4int nThreads = 1;
#endif
  TApplication theApp("App", 0, 0);

  G4String macro, events, geometry, radiator, physlist, outfile, 
    session,geomAng,geomPhi,batchmode,lensId,particle,momentum,testVal1,testVal2,testVal3,
    prismStepX,prismStepY,beamZ,beamX,timeRes,
    beamDimension, mcpLayout, infile = "hits.root", lutfile = "../data/lut.root";
  G4int firstevent(0), runtype(0), study(0), verbose(0);

  G4long myseed = 345354;
  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro     = argv[i+1];
    //    else if ( G4String(argv[i]) == "-u" ) session   = argv[i+1];
    else if ( G4String(argv[i]) == "-r" ) myseed    = atol(argv[i+1]);
    else if ( G4String(argv[i]) == "-o" ) outfile   = argv[i+1];
    else if ( G4String(argv[i]) == "-i" ) infile    = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) lutfile   = argv[i+1];
    else if ( G4String(argv[i]) == "-g" ) geometry  = argv[i+1];
    else if ( G4String(argv[i]) == "-h" ) radiator  = argv[i+1];
    else if ( G4String(argv[i]) == "-a" ) geomAng   = argv[i+1];
    else if ( G4String(argv[i]) == "-phi" ) geomPhi   = argv[i+1];
    else if ( G4String(argv[i]) == "-b" ) batchmode = argv[i+1];
    else if ( G4String(argv[i]) == "-f" ) firstevent= atoi(argv[i+1]);
    else if ( G4String(argv[i]) == "-e" ) events    = argv[i+1];
    else if ( G4String(argv[i]) == "-l" ) lensId    = argv[i+1];
    else if ( G4String(argv[i]) == "-x" ) particle  = argv[i+1];
    else if ( G4String(argv[i]) == "-p" ) momentum  = argv[i+1];
    else if ( G4String(argv[i]) == "-w" ) physlist  = argv[i+1];
    else if ( G4String(argv[i]) == "-s" ) runtype   = atoi(argv[i+1]);
    else if ( G4String(argv[i]) == "-study" ) study   = atoi(argv[i+1]);
    else if ( G4String(argv[i]) == "-z" ) beamDimension  = argv[i+1];
    else if ( G4String(argv[i]) == "-c" ) mcpLayout = argv[i+1];
    else if ( G4String(argv[i]) == "-t1" ) testVal1   = argv[i+1];
    else if ( G4String(argv[i]) == "-t2" ) testVal2   = argv[i+1];
    else if ( G4String(argv[i]) == "-t3" ) testVal3   = argv[i+1];
    else if ( G4String(argv[i]) == "-gsx" ) prismStepX   = argv[i+1];
    else if ( G4String(argv[i]) == "-gsy" ) prismStepY   = argv[i+1];
    else if ( G4String(argv[i]) == "-gz" ) beamZ    = argv[i+1];
    else if ( G4String(argv[i]) == "-gx" ) beamX    = argv[i+1];
    else if ( G4String(argv[i]) == "-tr" ) timeRes  = argv[i+1];
    else if ( G4String(argv[i]) == "-v" ) verbose   = atoi(argv[i+1]);
#ifdef G4MULTITHREADED
    else if ( G4String(argv[i]) == "-t" ) {
      nThreads = G4UIcommand::ConvertToInt(argv[i+1]);
    }
#endif
    else {
      PrintUsage();
      return 1;
    }
  }

  if(outfile=="" && runtype == 0) outfile = "hits.root"; // simulation
  if(outfile=="" && (runtype == 1 || runtype == 5 || runtype == 11)) outfile = "../data/lut.root";  // lookup table generation
  if(outfile=="" && runtype == 6) outfile = "focalplane.root";  // focal plane simulation
  if(outfile=="" && (runtype == 2 || runtype == 3)) outfile = "reco.root"; // reconstruction

  if(batchmode.size()) gROOT->SetBatch(kTRUE);
  if(!events.size()) events = "1";
  PrtManager::Instance(outfile,runtype);


  PrtManager::Instance()->SetStudy(study);
  if(physlist.size()) PrtManager::Instance()->SetPhysList(atoi(physlist));
  if(geometry.size()) PrtManager::Instance()->SetGeometry(atoi(geometry));
  if(radiator.size()) PrtManager::Instance()->SetRadiator(atoi(radiator));
  if(lensId.size())   PrtManager::Instance()->SetLens(atoi(lensId));
  if(mcpLayout.size())PrtManager::Instance()->SetMcpLayout(atoi(mcpLayout));
  if(beamDimension.size())   PrtManager::Instance()->SetBeamDimension(atof(beamDimension));
  if(testVal1.size())   PrtManager::Instance()->SetShift(atof(testVal1));
  if(testVal1.size())   PrtManager::Instance()->SetTest1(atof(testVal1));
  if(testVal2.size())   PrtManager::Instance()->SetTest2(atof(testVal2));
  if(testVal3.size())   PrtManager::Instance()->SetTest3(atof(testVal3));
  if(geomAng.size())   PrtManager::Instance()->SetAngle(atof(geomAng));
  if(geomPhi.size())   PrtManager::Instance()->SetPhi(atof(geomPhi));
  if(prismStepX.size())   PrtManager::Instance()->SetPrismStepX(atof(prismStepX));
  if(prismStepY.size())   PrtManager::Instance()->SetPrismStepY(atof(prismStepY));
  if(beamX.size())   PrtManager::Instance()->SetBeamX(atof(beamX));
  if(beamZ.size())   PrtManager::Instance()->SetBeamZ(atof(beamZ));
  if(timeRes.size())   PrtManager::Instance()->SetTimeRes(atof(timeRes));

  //if(beamDimension.size())PrtManager::Instance()->SetTest1(atoi(beamDimension));

  if(runtype == 2 || runtype == 3 || runtype == 4){
    PrtLutReco * reco = new PrtLutReco(infile.c_str(),lutfile.c_str(),verbose); 
    reco->Run(firstevent, atoi(events));
    return 0;
  }

  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
  if ( nThreads > 0 ) runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager * runManager = new G4RunManager;
#endif
  
  std::cout<<"SEED "<<myseed <<std::endl;
  G4Random::setTheSeed(myseed);
  // Seed the random number generator manually
  // if(myseed==345354) myseed = time(NULL);

  
  // Detector construction
  runManager->SetUserInitialization(new PrtDetectorConstruction());
  // Physics list
  runManager->SetUserInitialization(new PrtPhysicsList());
  // User action initialization
  runManager->SetUserInitialization(new PrtActionInitialization(outfile));
  // Initialize G4 kernel
  runManager->Initialize();

  //#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
  //#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( macro.size() ) {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }else {
    //  UImanager->ApplyCommand("/control/execute ../prt.mac");
    
    UImanager->ApplyCommand("/gun/direction 0 0 1");
    UImanager->ApplyCommand("/gun/position 0 0 -100 cm");
    UImanager->ApplyCommand("/Prt/geom/prtRotation 90 deg");
    UImanager->ApplyCommand("/gun/direction 0 0 1");
  }
   
  if ( particle.size() ) {
    int pdgid = 0;
    if(particle=="mix_pip") PrtManager::Instance()->SetMix(1);
    else if(particle=="mix_pik") PrtManager::Instance()->SetMix(2);
    else{
      G4String command = "/gun/particle ";
      UImanager->ApplyCommand(command+particle);
      if(particle=="proton") pdgid = 2212;
      if(particle=="pi+") pdgid = 211;
      if(particle=="pi0") pdgid = 111;
      if(particle=="kaon+") pdgid = 321;
      if(particle=="kaon-") pdgid = -321;
      if(particle=="mu-") pdgid = 13;
      if(particle=="e-") pdgid = 11;
    }    
    PrtManager::Instance()->SetParticle(pdgid);
  }

  if(momentum.size()) UImanager->ApplyCommand( "/gun/momentumAmp "+momentum);

  if(batchmode.size()){ // batch mode
    UImanager->ApplyCommand("/run/beamOn "+events);
  }else{  // UI session for interactive mode

    G4UIExecutive * ui = new G4UIExecutive(argc,argv,"Qt");
    UImanager->ApplyCommand("/control/execute ../vis.mac");
    if (ui->IsGUI()) UImanager->ApplyCommand("/control/execute gui.mac");
    UImanager->ApplyCommand("/run/beamOn "+events);
    //UImanager->ApplyCommand("/vis/ogl/printEPS");
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;

  return 0;
}

