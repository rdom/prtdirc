#include "PrtOpBoundaryProcess.h"

#include "PrtManager.h"

PrtOpBoundaryProcess::PrtOpBoundaryProcess() : G4OpBoundaryProcess() {
  fLensId = PrtManager::Instance()->getRun()->getLens();
  fRunType = PrtManager::Instance()->getRun()->getRunType();
}

G4VParticleChange *PrtOpBoundaryProcess::PostStepDoIt(const G4Track &aTrack, const G4Step &aStep) {
  G4StepPoint *pPreStepPoint = aStep.GetPreStepPoint();
  G4StepPoint *pPostStepPoint = aStep.GetPostStepPoint();
  G4VParticleChange *particleChange = G4OpBoundaryProcess::PostStepDoIt(aTrack, aStep);

  G4String aname = pPreStepPoint->GetPhysicalVolume()->GetName();
  G4String bname = pPostStepPoint->GetPhysicalVolume()->GetName();

  // int parentId = aTrack.GetParentID();
  // if(parentId==1) particleChange->ProposeTrackStatus(fStopAndKill);

  // ideal focusing
  if (fLensId == 10) {
    G4String ParticleName = aTrack.GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
    if (ParticleName == "opticalphoton") {
      double endofbar = 1250 / 2.;
      G4ThreeVector theGlobalPoint1 = pPostStepPoint->GetPosition();
      G4TouchableHistory *touchable = (G4TouchableHistory *)(pPostStepPoint->GetTouchable());
      G4ThreeVector lpoint = touchable->GetHistory()->GetTransform(1).TransformPoint(theGlobalPoint1);

      if (lpoint.getZ() < endofbar + 0.0001 && lpoint.getZ() > endofbar - 0.0001) {
        G4ThreeVector ww = pPreStepPoint->GetTouchableHandle()->GetHistory()->GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0, 0, endofbar));
        if (aname != "wBar")
          particleChange->ProposeTrackStatus(fStopAndKill);
        else {
          G4Navigator *theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
          theNavigator->LocateGlobalPointWithinVolume(ww);
          aParticleChange.ProposePosition(ww.getX(), ww.getY(), ww.getZ());
        }
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      }
    }
  }

  // bar surface scattering
  // if(1){
  //   G4String ParticleName = aTrack.GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
  //   if (ParticleName == "opticalphoton" && aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wBar"){
  //     G4double z = 0.5*PrtManager::Instance()->getRun()->getRadiatorL();
  //     G4double w = 0.5*PrtManager::Instance()->getRun()->getRadiatorW();
  //     G4double h = 0.5*PrtManager::Instance()->getRun()->getRadiatorH();

  //     G4ThreeVector theGlobalPoint1 = pPostStepPoint->GetPosition();
  //     G4TouchableHistory* touchable = (G4TouchableHistory*)(pPostStepPoint->GetTouchable());
  //     G4ThreeVector lpoint =  touchable->GetHistory()->GetTransform(1).TransformPoint(theGlobalPoint1);
  // 	{
  // 	  std::cout<<w<<" lpoint "<<lpoint<<std::endl;

  // 	  if(lpoint.getY() > w-0.01)
  // 	  {
  // 	  std::cout<<h<<" lpoint.getX() "<<lpoint.getX()<<" "<<lpoint.getY()<<std::endl;

  // 	  G4ThreeVector ww  = pPreStepPoint->GetTouchableHandle()->GetHistory()->
  // 	    GetTopTransform().Inverse().TransformPoint(G4ThreeVector(lpoint.getX(),lpoint.getY(),lpoint.getZ()));

  // 	  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  // 	  theNavigator->LocateGlobalPointWithinVolume(ww);
  // 	  aParticleChange.ProposePosition(ww.getX(), ww.getY(),ww.getZ());
  // 	  // aParticleChange.ProposeMomentumDirection(G4double Px, G4double Py, G4double Pz)

  // 	  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  // 	}
  //     }
  //   }
  // }

  if (fRunType == 1 && pPostStepPoint->GetPosition().z() < pPreStepPoint->GetPosition().z()) {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  // kill reflections from FP
  if (fRunType == 0 && aname == "wPrizm") {
    auto touchable = (G4TouchableHistory *)(pPreStepPoint->GetTouchable());
    auto pos1 = touchable->GetHistory()->GetTopTransform().TransformPoint(pPreStepPoint->GetPosition());
    auto pos2 = touchable->GetHistory()->GetTopTransform().TransformPoint(pPostStepPoint->GetPosition());
    if (pos2.y() > pos1.y()) particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if (fRunType == 5 && aname == "wDirc" && bname == "wPrizm" && GetStatus() == FresnelRefraction) {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  // kill photons outside bar and prizm
  if (GetStatus() == FresnelRefraction && bname == "wDirc") {
    // rd for air gap
    if (fLensId != 4) particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if ((aname == "wLens1" || aname == "wLens2" || aname == "wLens3") && bname == "wDirc") {
    // if(fLensId!=4) particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if ((aname == "wOpticalGreased" || aname == "wOpticalGrease") && bname == "wDirc") {
    if (fLensId == 2) particleChange->ProposeTrackStatus(fStopAndKill);
  }

  return particleChange;
}
