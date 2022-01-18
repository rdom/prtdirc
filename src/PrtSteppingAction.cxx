#include "PrtSteppingAction.h"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"

#include "G4TrackingManager.hh"
#include "PrtManager.h"

PrtSteppingAction::PrtSteppingAction() : G4UserSteppingAction() {
  fScintillationCounter = 0;
  fCerenkovCounter = 0;
  fEventNumber = -1;
  fRunType = PrtManager::Instance()->getRun()->getRunType();
}

PrtSteppingAction::~PrtSteppingAction() {}

void PrtSteppingAction::UserSteppingAction(const G4Step *step) {

  G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  if (eventNumber != fEventNumber) {
    fEventNumber = eventNumber;
    fScintillationCounter = 0;
    fCerenkovCounter = 0;
  }

  G4Track *track = step->GetTrack();

  // int parentId = track->GetParentID();  
  // std::cout<<"parentId   "<<parentId <<std::endl;
  // std::cout<<"length  "<<track->GetTrackLength() << " step num "<< track->GetCurrentStepNumber()
  // <<std::endl;

  if (track->GetCurrentStepNumber() > 50000 || track->GetTrackLength() > 10000) {
    track->SetTrackStatus(fStopAndKill);
  }
  
  if (fRunType == 11 || fRunType == 1) {
    G4String aname = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
    G4String bname = step->GetPostStepPoint()->GetPhysicalVolume()->GetName();

    if (fRunType == 11 && aname == "wBar" && bname == "wOpticalGrease") {
      G4ThreeVector dir = step->GetPreStepPoint()->GetMomentum();
      TVector3 v(dir.x(), dir.y(), dir.z());
      v.RotateY(-(TMath::Pi() - 20 * TMath::Pi() / 180.));
      PrtManager::Instance()->getEvent()->setMomentum(v);
      PrtManager::Instance()->getEvent()->setTime(step->GetPreStepPoint()->GetLocalTime());
    }

    // stop photons at the edge of the lens for LUT
    if (fRunType == 1 && aname == "wLens3" && bname == "wDirc") {
      track->SetTrackStatus(fStopAndKill);
    }
  }

  // if(aname=="Bar" && bname=="ExpHall" ) track->SetTrackStatus(fStopAndKill);
  // if(step->GetPreStepPoint()->GetPosition().z()>1200 ) track->SetTrackStatus(fStopAndKill);
  // G4String ParticleName =
  // track->GetDynamicParticle()->GetParticleDefinition()->GetParticleName(); if (ParticleName ==
  // "opticalphoton") return; const std::vector<const G4Track *> *secondaries =
  // step->GetSecondaryInCurrentStep(); if (secondaries->size() > 0) {
  //   for (unsigned int i = 0; i < secondaries->size(); ++i) {
  //     if (secondaries->at(i)->GetParentID() > 0) {
  //       if (secondaries->at(i)->GetDynamicParticle()->GetParticleDefinition() ==
  //       G4OpticalPhoton::OpticalPhotonDefinition()) {
  //         if (secondaries->at(i)->GetCreatorProcess()->GetProcessName() == "Scintillation")
  //         fScintillationCounter++; if (secondaries->at(i)->GetCreatorProcess()->GetProcessName()
  //         == "Cerenkov") fCerenkovCounter++;
  //       }
  //     }
  //   }
  // }
}
