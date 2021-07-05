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
  
  if (fRunType == 5) {
    if (GetStatus() == TotalInternalReflection) {

      double endofbar = 1200.06 / 2.;
      G4TouchableHistory *touchable = (G4TouchableHistory *)(pPreStepPoint->GetTouchable());

      G4ThreeVector gpos = pPreStepPoint->GetPosition();
      G4ThreeVector gdir = pPreStepPoint->GetMomentumDirection();
      G4ThreeVector lpos = touchable->GetHistory()->GetTransform(1).TransformPoint(gpos);
      G4ThreeVector ldir = touchable->GetHistory()->GetTransform(1).TransformAxis(gdir);
      if (lpos.z() < endofbar - 0.01) {
        double h = 17.15;
        double w = 34.93;
        double lenz = endofbar - lpos.z();
        if (ldir.z() < 0) lenz = 4 * endofbar - lenz;
        double len = fabs(lenz / cos(ldir.theta()));
        double theta = ldir.theta();

        int signx = 1, signy = 1, signz = 1;
        double dx, newx, dy, newy;
        double x = len * ldir.x();
        double y = len * ldir.y();

	dx = x + 0.5 * h + lpos.x();
	dy = y + 0.5 * w + lpos.y();			       

        if (x > 0) {
          newx = fmod(dx, h) - 0.5 * h;
          if (dx < h) newx = x + lpos.x();
        } else {
	  dx -= h;
          newx = fmod(dx, h) + 0.5 * h;
	  if (dx > -h) newx = x + lpos.x();
        }
	
        if (y > 0) {
          newy = fmod(dy, w) - 0.5 * w;
	  if (dy < w) newy = y + lpos.y();
        } else {
	  dy -= w;
          newy = fmod(dy, w) + 0.5 * w;
	  if (dy > -w) newy = y + lpos.y();
        }

	if (ldir.z() < 0) signz = -1;
	if (int(floor(fabs(dx) / h)) % 2 == 1) signx = -1;
	if (int(floor(fabs(dy) / w)) % 2 == 1) signy = -1;

	newx *= signx;
	newy *= signy;	
	
        if (0) { // transport_efficiency
          double pi(4 * atan(1));
          double roughness(0.5); // nm
          double angleX = ldir.angle(G4ThreeVector(1, 0, 0));
          double angleY = ldir.angle(G4ThreeVector(0, 1, 0));
          if (angleX > 0.5 * pi) angleX = pi - angleX;
          if (angleY > 0.5 * pi) angleY = pi - angleY;

          double lengthx = fabs(len * ldir.x());
          double lengthy = fabs(len * ldir.y());

          int nBouncesX = (int)(lengthx) / h; // 17 bar height
          int nBouncesY = (int)(lengthy) / w; // 36 bar width

          double wavelength = 1.2398 / (aTrack.GetMomentum().mag() * 1E6) * 1000;
          double ll = wavelength * wavelength;
          double n_quartz =
            sqrt(1. + (0.696 * ll / (ll - pow(0.068, 2))) + (0.407 * ll / (ll - pow(0.116, 2))) +
                 0.897 * ll / (ll - pow(9.896, 2)));
          double bounce_probX =
            1 - pow(4 * pi * cos(angleX) * roughness * n_quartz / wavelength, 2);
          double bounce_probY =
            1 - pow(4 * pi * cos(angleY) * roughness * n_quartz / wavelength, 2);

          double totalProb = pow(bounce_probX, nBouncesX) * pow(bounce_probY, nBouncesY);

          if (G4UniformRand() > totalProb) aParticleChange.ProposeTrackStatus(fStopAndKill);
        }

        double tl = aParticleChange.GetTrueStepLength();
        aParticleChange.ProposeTrueStepLength(len);

        G4ThreeVector ww = pPreStepPoint->GetTouchableHandle()
                             ->GetHistory()
                             ->GetTopTransform()
                             .Inverse()
                             .TransformPoint(G4ThreeVector(newx, newy, endofbar - 0.01));
        G4Navigator *theNavigator =
          G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
        theNavigator->LocateGlobalPointWithinVolume(ww);
        aParticleChange.ProposePosition(ww.x(), ww.y(), ww.z());

        G4ThreeVector rdir =
          pPreStepPoint->GetTouchableHandle()
            ->GetHistory()
            ->GetTopTransform()
            .Inverse()
            .TransformAxis(G4ThreeVector(signx * ldir.x(), signy * ldir.y(), signz * ldir.z()));
        aParticleChange.ProposeMomentumDirection(rdir);

        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      }
    }
  }

  G4String aname = pPreStepPoint->GetPhysicalVolume()->GetName();
  G4String bname = pPostStepPoint->GetPhysicalVolume()->GetName();

  // int parentId = aTrack.GetParentID();
  // if(parentId==1) particleChange->ProposeTrackStatus(fStopAndKill);

  // bar surface scattering
  // if(1){
  //   G4String ParticleName =
  //   aTrack.GetDynamicParticle()->GetParticleDefinition()->GetParticleName(); if (ParticleName
  //   == "opticalphoton" && aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wBar"){
  //     G4double z = 0.5*PrtManager::Instance()->getRun()->getRadiatorL();
  //     G4double w = 0.5*PrtManager::Instance()->getRun()->getRadiatorW();
  //     G4double h = 0.5*PrtManager::Instance()->getRun()->getRadiatorH();

  //     G4ThreeVector theGlobalPoint1 = pPostStepPoint->GetPosition();
  //     G4TouchableHistory* touchable = (G4TouchableHistory*)(pPostStepPoint->GetTouchable());
  //     G4ThreeVector lpoint =
  //     touchable->GetHistory()->GetTransform(1).TransformPoint(theGlobalPoint1);
  // 	{
  // 	  std::cout<<w<<" lpoint "<<lpoint<<std::endl;

  // 	  if(lpoint.getY() > w-0.01)
  // 	  {
  // 	  std::cout<<h<<" lpoint.getX() "<<lpoint.getX()<<" "<<lpoint.getY()<<std::endl;

  // 	  G4ThreeVector ww  = pPreStepPoint->GetTouchableHandle()->GetHistory()->
  // 	    GetTopTransform().Inverse().TransformPoint(G4ThreeVector(lpoint.getX(),lpoint.getY(),lpoint.getZ()));

  // 	  G4Navigator* theNavigator =
  // G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  // 	  theNavigator->LocateGlobalPointWithinVolume(ww);
  // 	  aParticleChange.ProposePosition(ww.getX(), ww.getY(),ww.getZ());
  // 	  // aParticleChange.ProposeMomentumDirection(G4double Px, G4double Py, G4double Pz)

  // 	  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  // 	}
  //     }
  //   }
  // }

  // kill photons outside bar and prizm
  if (GetStatus() == FresnelRefraction && bname == "wDirc") {
    // rd for air gap
    if (fLensId != 4) particleChange->ProposeTrackStatus(fStopAndKill);
  }

  // kill reflections from FP
  if ((fRunType == 0 || fRunType == 5) && aname == "wPrizm") {
    auto touchable = (G4TouchableHistory *)(pPreStepPoint->GetTouchable());
    auto pos1 =
      touchable->GetHistory()->GetTopTransform().TransformPoint(pPreStepPoint->GetPosition());
    auto pos2 =
      touchable->GetHistory()->GetTopTransform().TransformPoint(pPostStepPoint->GetPosition());
    if (pos2.y() > pos1.y()) particleChange->ProposeTrackStatus(fStopAndKill);
  }

  // if ((aname == "wLens1" || aname == "wLens2" || aname == "wLens3") && bname == "wDirc") {
  //   // if(fLensId!=4) particleChange->ProposeTrackStatus(fStopAndKill);
  // }

  // ideal focusing
  if (fLensId == 10) {
    G4String ParticleName = aTrack.GetDynamicParticle()->GetParticleDefinition()->GetParticleName();
    if (ParticleName == "opticalphoton") {
      double endofbar = 1200.06 / 2.;
      G4ThreeVector theGlobalPoint1 = pPostStepPoint->GetPosition();
      G4TouchableHistory *touchable = (G4TouchableHistory *)(pPostStepPoint->GetTouchable());
      G4ThreeVector lpoint =
        touchable->GetHistory()->GetTransform(1).TransformPoint(theGlobalPoint1);

      if (lpoint.getZ() < endofbar + 0.0001 && lpoint.getZ() > endofbar - 0.0001) {
        G4ThreeVector ww = pPreStepPoint->GetTouchableHandle()
                             ->GetHistory()
                             ->GetTopTransform()
                             .Inverse()
                             .TransformPoint(G4ThreeVector(0, 0, endofbar));
        if (aname != "wBar") particleChange->ProposeTrackStatus(fStopAndKill);
        else {
          G4Navigator *theNavigator =
            G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
          theNavigator->LocateGlobalPointWithinVolume(ww);
          aParticleChange.ProposePosition(ww.getX(), ww.getY(), ww.getZ());
        }
        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      }
    }
  }

  if (fRunType == 1 && pPostStepPoint->GetPosition().z() < pPreStepPoint->GetPosition().z()) {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if (fRunType == 7 && aname == "wDirc" && bname == "wPrizm" && GetStatus() == FresnelRefraction) {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if (fLensId == 2 &&
      ((aname == "wOpticalGreased" || aname == "wOpticalGrease") && bname == "wDirc")) {
    particleChange->ProposeTrackStatus(fStopAndKill);
  }

  return particleChange;
}
