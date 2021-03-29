#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "PrtPrimaryGeneratorAction.h"
#include "PrtPrimaryGeneratorMessenger.h"
#include "PrtManager.h"

PrtPrimaryGeneratorAction::PrtPrimaryGeneratorAction() : G4VUserPrimaryGeneratorAction(), fParticleGun(0) {
  int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  // create a messenger for this class
  fGunMessenger = new PrtPrimaryGeneratorMessenger(this);

  // default kinematic
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  fParticleP = particleTable->FindParticle("proton");
  fParticlePi = particleTable->FindParticle("pi+");
  fParticleKaon = particleTable->FindParticle("kaon+");

  fParticleGun->SetParticleDefinition(fParticleP);
  fParticleGun->SetParticleTime(0.0 * ns);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0 * cm, 0.0 * cm, 0.0 * cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
  fParticleGun->SetParticleEnergy(7 * MeV);

  // int mid=-1, pid=-1;
  // G4ThreeVector vdirc,vmcp[12],vpix[64];
  // auto store = G4PhysicalVolumeStore::GetInstance();
  // for (size_t i=0;i<store->size();i++){
  //   if((*store)[i]->GetName()=="wDirc")  vdirc = (*store)[i]->GetTranslation();
  //   if((*store)[i]->GetName()=="wMcp")   vmcp[++mid] = (*store)[i]->GetTranslation();
  //   if((*store)[i]->GetName()=="wPixel") vpix[++pid] = (*store)[i]->GetTranslation();
  // }

  // for(auto m=0; m<mid; m++){
  //   for(auto p=0; p<pid; p++){
  //     gpix[m][p] = vdirc+(vmcp[m]+vpix[p]).rotateY(PrtManager::Instance()->GetAngle()*deg-180*deg);
  //   }
  // }
  iter = 0;
}

PrtPrimaryGeneratorAction::~PrtPrimaryGeneratorAction() {
  delete fParticleGun;
  delete fGunMessenger;
}

void PrtPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {

  double x, y, z;
  PrtRun *run = PrtManager::Instance()->getRun();
  double radiatorL = run->getRadiatorL();
  double radiatorW = run->getRadiatorW();
  double radiatorH = run->getRadiatorH();
  int pid = run->getPid();

  PrtManager::Instance()->addEvent(PrtEvent());

  if (pid > 10000) {
    if (pid == 211 || pid == 0) {
      if (pid == 10001) {
        fParticleGun->SetParticleDefinition(fParticleP);
        PrtManager::Instance()->getEvent()->setPid(2212);
      } else if (pid == 10002) {
        fParticleGun->SetParticleDefinition(fParticleKaon);
        PrtManager::Instance()->getEvent()->setPid(321);
      }
    } else {
      fParticleGun->SetParticleDefinition(fParticlePi);
      PrtManager::Instance()->getEvent()->setPid(211);
    }
  }

  if (run->getBeamSize() == -1) { // random momentum
    fParticleGun->SetParticleMomentum(G4ThreeVector(0, 0, 4.0 * GeV * G4UniformRand()));
  }
  if (run->getBeamSize() > 0) { // smearing and divergence
    double sigma = run->getBeamSize() * mm;
    z = fParticleGun->GetParticlePosition().z();

    // // gaussian smearing
    // x = G4RandGauss::shoot(0,sigma);
    // y = G4RandGauss::shoot(0,sigma);

    // box smearing
    x = (0.5 - G4UniformRand()) * sigma;
    y = (0.5 - G4UniformRand()) * sigma;

    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
    PrtManager::Instance()->getEvent()->setPosition(TVector3(x, y, z));
    double angle = -G4UniformRand() * M_PI;
    G4ThreeVector vec(0, 0, 1);
    vec.setTheta(G4RandGauss::shoot(0, 0.0015)); // beam divergence
    vec.setPhi(2 * M_PI * G4UniformRand());

    fParticleGun->SetParticleMomentumDirection(vec);
  }
  if (run->getRunType() == 1) { // LUT generation
    // fParticleGun->SetParticlePosition(G4ThreeVector(radiatorH*(0.5-G4UniformRand()),radiatorW*(0.5-G4UniformRand()),radiatorL/2.-0.1));
    fParticleGun->SetParticlePosition(G4ThreeVector(run->getPrismStepY(), //+5-10*G4UniformRand(),
                                                    run->getPrismStepX(), //+10-20*G4UniformRand(),
                                                    radiatorL / 2. - 0.1));

    // if(iter>3) iter=0;
    // double posX[]={-3,3,3,-3};
    // double posY[]={-6,6,-6,6};
    // fParticleGun->SetParticlePosition(G4ThreeVector(run->getPrismStepY()+posX[iter],//+5-10*G4UniformRand(),
    // 						    run->getPrismStepX()+posY[iter],//+10-20*G4UniformRand(),
    // 						    radiatorL/2.-0.1));
    // iter++;

    // fParticleGun->SetParticlePosition(G4ThreeVector(run->getPrismStepY()+8-16*G4UniformRand(),
    // 						    run->getPrismStepX()+16-32*G4UniformRand(),
    // 						    radiatorL/2.-0.1));

    double angle = -G4UniformRand() * M_PI;
    G4ThreeVector vec(0, 0, 1);
    vec.setTheta(acos(G4UniformRand()));
    vec.setPhi(2 * M_PI * G4UniformRand());

    // vec.rotateY(-M_PI/2.);
    fParticleGun->SetParticleMomentumDirection(vec);
  }
  if (run->getRunType() == 5) { // calibration light
    double shift = run->getTest3();

    fParticleGun->SetParticlePosition(G4ThreeVector(-radiatorL / 2. + 0.1 - shift, 0, 5 + tan(45 * M_PI / 180.) * shift + 25));
    double angle = -G4UniformRand() * M_PI;
    G4ThreeVector vec(0, 0, 1);
    vec.setTheta(acos(G4UniformRand()));
    vec.setPhi(2 * M_PI * G4UniformRand());

    vec.rotateY(-M_PI / 2.);
    fParticleGun->SetParticleMomentumDirection(vec);
  }
  if (run->getRunType() == 6) { // for determining focal plane of the lens
    double shiftx = -radiatorL / 2. + 0.1;
    double shifty = radiatorW / 2. - G4UniformRand() * radiatorW;
    double shiftz = radiatorH / 2. - G4UniformRand() * radiatorH;

    double angle = 0.7 * (M_PI / 2. - G4UniformRand() * M_PI);
    G4ThreeVector vec(0, 0, 1);
    vec.setTheta(angle);
    /// vec.setTheta(acos(G4UniformRand()));
    // vec.setPhi(2*M_PI*G4UniformRand());
    // std::cout<<"angle "<<angle*180/M_PI <<std::endl;

    double lensThickness = 15;
    double separation = run->getBeamSize();
    if (separation < 0.001) separation = 10;
    double rotShiftX = 0.5 * separation * std::cos(angle) + (0.5 * lensThickness + 0.1) * std::tan(angle);
    double rotShiftY = -0.5 * radiatorL + 0.1;

    // fParticleGun->SetParticlePosition(G4ThreeVector(shiftx,shifty,shiftz));
    // fParticleGun->SetParticlePosition(G4ThreeVector(rotShiftY, 0,-rotShiftX));
    fParticleGun->SetParticlePosition(G4ThreeVector(rotShiftY, 0, 0.5 * separation));

    vec.rotateY(-M_PI / 2.);
    fParticleGun->SetParticleMomentumDirection(vec);

    fParticleGun->GeneratePrimaryVertex(anEvent);
    rotShiftX = -0.5 * separation * std::cos(angle) + (0.5 * lensThickness + 0.1) * std::tan(angle);
    shiftx = -radiatorL / 2. + 0.1;
    shifty = radiatorW / 2. - G4UniformRand() * radiatorW;
    shiftz = radiatorH / 2. - G4UniformRand() * radiatorH;
    // fParticleGun->SetParticlePosition(G4ThreeVector(shiftx,shifty,shiftz));
    // fParticleGun->SetParticlePosition(G4ThreeVector(rotShiftY,0,-rotShiftX));
    fParticleGun->SetParticlePosition(G4ThreeVector(rotShiftY, 0, -0.5 * separation));
  }

  fParticleGun->GeneratePrimaryVertex(anEvent);

  G4ThreeVector dir = fParticleGun->GetParticleMomentumDirection();
  dir *= fParticleGun->GetParticleMomentum();
  PrtManager::Instance()->getEvent()->setMomentum(TVector3(dir.x(), dir.y(), dir.z()));  
}

void PrtPrimaryGeneratorAction::SetOptPhotonPolar() {
  double angle = G4UniformRand() * 360.0 * deg;
  SetOptPhotonPolar(angle);
}

void PrtPrimaryGeneratorAction::SetOptPhotonPolar(double angle) {
  if (fParticleGun->GetParticleDefinition()->GetParticleName() != "opticalphoton") {
    G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
              "the particleGun is not an opticalphoton "
           << fParticleGun->GetParticleDefinition()->GetParticleName() << G4endl;
    return;
  }

  G4ThreeVector normal(1., 0., 0.);
  G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
  G4ThreeVector product = normal.cross(kphoton);
  double modul2 = product * product;

  G4ThreeVector e_perpend(0., 0., 1.);
  if (modul2 > 0.) e_perpend = (1. / std::sqrt(modul2)) * product;
  G4ThreeVector e_paralle = e_perpend.cross(kphoton);

  G4ThreeVector polar = std::cos(angle) * e_paralle + std::sin(angle) * e_perpend;
  fParticleGun->SetParticlePolarization(polar);
}
