#include "TInterpreter.h"
#include "G4SystemOfUnits.hh"

#include "PrtManager.h"
#include "PrtHit.h"
#include "PrtLutNode.h"

PrtManager *PrtManager::fInstance = NULL;

PrtManager::PrtManager(TString filename, PrtRun *run) {

  fOutName = filename;
  fRun = run;
  fRunType = fRun->getRunType();
  fEvent = new PrtEvent();

  fRootFile = new TFile(filename, "RECREATE");
  fRun->setMc(true);
  
  fRunTree = new TTree("header", "run info");
  fRunTree->Branch("PrtRun", "PrtRun", &fRun, 64000, 2);
  
  switch (fRunType) {
  
    case 0:
    case 5:
    case 6:
    case 20:   // NEW uniform sampling mode
      fTree = new TTree("data", "Prototype hits tree");
      fTree->Branch("PrtEvent", "PrtEvent", &fEvent, 64000, 2);
      break;
  
    case 1:
    case 7:
    case 11:
      fLut = new TClonesArray("PrtLutNode");
      fTree = new TTree("prtlut", "Look-up table for the geometrical reconstruction.");
      fTree->Branch("LUT", &fLut, 256000, 2);
  
      for (Long64_t i = 0; i < fRun->getNpmt()*fRun->getNpix(); i++) {
        new ((*fLut)[i]) PrtLutNode(i);
      }
      break;
  
    default:
      std::cerr << "ERROR: Unsupported fRunType = "
                << fRunType << std::endl;
      std::exit(EXIT_FAILURE);
  }

  fnX1 = TVector3(1, 0, 0);
  fnY1 = TVector3(0, 1, 0);
  fCriticalAngle = asin(1.00028 / 1.47125);

  std::cout << "PrtManager has been successfully initialized. " << std::endl;
}

PrtManager *PrtManager::Instance(TString outfile, PrtRun *run) {
  if (!fInstance) {
    std::cout << "Info in (PrtManager::Instance): Making a new instance. " << std::endl;
    fInstance = new PrtManager(outfile, run);
  }
  return fInstance;
}

void PrtManager::addEvent(PrtEvent event) {

  switch (static_cast<PrtRunType>(fRunType)) {

    case PrtRunType::Data:
    case PrtRunType::Data5:
    case PrtRunType::Data6:
    case PrtRunType::UniformGun:
      fEvent = new PrtEvent(event);
      break;

    default:
      break;  // LUT modes do nothing
  }
}

void PrtManager::addHit(PrtHit hit,
                        TVector3 localpos,
                        TVector3 digipos,
                        TVector3 vertex) {

  switch (static_cast<PrtRunType>(fRunType)) {

    case PrtRunType::Data:
    case PrtRunType::Data5:
    case PrtRunType::Data6:
    case PrtRunType::UniformGun:
      fEvent->addHit(hit);
      break;

    case PrtRunType::Lut:
    case PrtRunType::Lut7:
    case PrtRunType::Lut11: {

      if (fMomentum.Angle(fnX1) > fCriticalAngle &&
          fMomentum.Angle(fnY1) > fCriticalAngle) {

        int ch = hit.getChannel();
        double time = hit.getLeadTime();

        if (fRunType == 11)
          time -= fEvent->getTime();

        ((PrtLutNode*)(fLut->At(ch)))
          ->AddEntry(ch, fMomentum,
                     hit.getPathInPrizm(),
                     0, time,
                     localpos, digipos,
                     0, vertex);
      }
      break;
    }

    default:
      std::cerr << "ERROR: Unsupported runType in addHit()" << std::endl;
      break;
  }
}

void PrtManager::addHit(PrtHit hit) {

  switch (fRunType) {

    case 0:
    case 5:
    case 6:
    case 20:
      fEvent->addHit(hit);
      break;

    default:
      break;  // LUT modes ignore this
  }
}

void PrtManager::fill() {

  switch (static_cast<PrtRunType>(fRunType)) {

    case PrtRunType::Data:
    case PrtRunType::Data5:
    case PrtRunType::Data6:
    case PrtRunType::UniformGun:
      fTree->Fill();
      fEvent->Clear();
      break;

    default:
      break;
  }
}

void PrtManager::save(){
  if (fRootFile) {
    fRunTree->Fill();
    fRootFile->Write();
  }
}

void PrtManager::fillLut() {

  switch (static_cast<PrtRunType>(fRunType)) {

    case PrtRunType::Lut:
    case PrtRunType::Lut7:
    case PrtRunType::Lut11:
      fTree->Fill();
      break;

    default:
      break;
  }
}
