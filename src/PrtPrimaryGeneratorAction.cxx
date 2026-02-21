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

PrtPrimaryGeneratorAction::PrtPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(), fParticleGun(0) {
  int n_particle = 1;
  fRun = PrtManager::Instance()->getRun();
  double mom = fRun->getMomentum();
  fRadiatorL = fRun->getRadiatorL();
  fRadiatorW = fRun->getRadiatorW();
  fRadiatorH = fRun->getRadiatorH();

  fParticleGun = new G4ParticleGun(n_particle);

  // create a messenger for this class
  fGunMessenger = new PrtPrimaryGeneratorMessenger(this);

  // default kinematic
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  fParticle[4] = particleTable->FindParticle("proton");
  fParticle[3] = particleTable->FindParticle("kaon+");
  fParticle[2] = particleTable->FindParticle("pi+");
  fParticle[1] = particleTable->FindParticle("mu+");
  fParticle[0] = particleTable->FindParticle("e-");

  fParticleOP = particleTable->FindParticle("opticalphoton");

  fParticleGun->SetParticleDefinition(fParticleOP);
  fParticleGun->SetParticleTime(0.0 * ns);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0 * cm, 0.0 * cm, 0.0 * cm));
  fParticleGun->SetParticleMomentum(G4ThreeVector(0, 0, mom * GeV));

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
  //     gpix[m][p] =
  //     vdirc+(vmcp[m]+vpix[p]).rotateY(PrtManager::Instance()->GetAngle()*deg-180*deg);
  //   }
  // }

  iter = 0;
  fPid = 4;
}

PrtPrimaryGeneratorAction::~PrtPrimaryGeneratorAction() {
  delete fParticleGun;
  delete fGunMessenger;
}

void PrtPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {

  double x, y, z;

  PrtManager::Instance()->addEvent(PrtEvent());
  int pdg = fRun->getPid();

  if (pdg > 0) {
    if (pdg == 2212) fPid = 4;
    else if (pdg == 321) fPid = 3;
    else if (pdg == 211) fPid = 2;
    else if (pdg == 10001 && fPid > 2) fPid = 2;
    else if (pdg == 10001 && fPid == 2) fPid = 4;
    else if (pdg == 10002 && fPid > 2) fPid = 2;
    else if (pdg == 10002 && fPid == 2) fPid = 3;

    PrtManager::Instance()->getEvent()->setPid(fPid);
    PrtManager::Instance()->getEvent()->setTof(fRun->getTheta()); // save thata into tof variable, as tof is not used in sims
    fParticleGun->SetParticleDefinition(fParticle[fPid]);
  } else {
    fParticleGun->SetParticleDefinition(fParticleOP);
  }
  // if (fRun->getRunType() == 11)  fParticleGun->SetParticleMomentum(G4ThreeVector(0, 0, (7 - 4 * G4UniformRand()) * GeV));    
  
  if (fRun->getBeamSize() > 0) { // smearing and divergence
    double sigma = fRun->getBeamSize() * mm;
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
    double divergence = fRun->getTest2(); //0.0015
    vec.setTheta(G4RandGauss::shoot(0, divergence)); // beam divergence
    vec.setPhi(2 * M_PI * G4UniformRand());

    fParticleGun->SetParticleMomentumDirection(vec);
  }

if (fRun->getRunType() == 20) {

  G4double momentum = 0.5 + G4UniformRand()*(3.5-0.5);
  fParticleGun->SetParticleMomentum(momentum*GeV);
  
  G4ThreeVector b( std::sin(22*deg), 0, std::cos(22*deg) );

  b = b.unit();                   // bar coordinate frame
  G4ThreeVector u(0,1,0);         // bar coordinate frame
  G4ThreeVector w = b.cross(u);   // bar coordinate frame
  w = w.unit();

  // ---- sample angle ----
  double cosMin = std::cos(140*deg);
  double cosMax = std::cos(22*deg);
  double cosTheta = cosMin + G4UniformRand()*(cosMax - cosMin);
  double theta = std::acos(cosTheta);

  double phi = 0.025*M_PI*G4UniformRand();

  // ---- build direction ----
  G4ThreeVector dir =
      cosTheta * b
    + std::sin(theta)*std::cos(phi) * u
    + std::sin(theta)*std::sin(phi) * w;

  dir = dir.unit();

  G4ThreeVector entry(0, 0, 50*cm);

  G4double L = 50*cm;   // same as your upstream distance
  G4ThreeVector start = entry - L * dir;

  fParticleGun->SetParticlePosition(start);
  fParticleGun->SetParticleMomentumDirection(dir);
}


  if (fRun->getRunType() == 1) { // LUT generation
    // fParticleGun->SetParticlePosition(G4ThreeVector(fRadiatorH*(0.5-G4UniformRand()),fRadiatorW*(0.5-G4UniformRand()),fRadiatorL/2.-0.1));
    fParticleGun->SetParticlePosition(G4ThreeVector(fRun->getPrismStepY(), //+5-10*G4UniformRand(),
                                                    fRun->getPrismStepX(), //+10-20*G4UniformRand(),
                                                    fRadiatorL / 2. - 0.1));

    // if(iter>3) iter=0;
    // double posX[]={-3,3,3,-3};
    // double posY[]={-6,6,-6,6};
    // fParticleGun->SetParticlePosition(G4ThreeVector(fRun->getPrismStepY()+posX[iter],//+5-10*G4UniformRand(),
    // 						    fRun->getPrismStepX()+posY[iter],//+10-20*G4UniformRand(),
    // 						    fRadiatorL/2.-0.1));
    // iter++;

    // fParticleGun->SetParticlePosition(G4ThreeVector(fRun->getPrismStepY()+8-16*G4UniformRand(),
    // 						    fRun->getPrismStepX()+16-32*G4UniformRand(),
    // 						    fRadiatorL/2.-0.1));

    double angle = -G4UniformRand() * M_PI;
    G4ThreeVector v(0, 0, 1);
    v.setTheta(acos(G4UniformRand()));
    v.setPhi(2 * M_PI * G4UniformRand());

    fParticleGun->SetParticleMomentumDirection(v);
    PrtManager::Instance()->setMomentum(TVector3(v.x(), v.y(), v.z()));

    // for (int i = 0; i < 10; i++) {
    //   for (int j = 0; j < 20; j++) {
    //     fParticleGun->SetParticlePosition(G4ThreeVector(
    //       fRun->getPrismStepY() + 5 - i, fRun->getPrismStepX() + 10 - j, fRadiatorL / 2. - 0.1));
    //     fParticleGun->SetParticleMomentumDirection(v);
    //     fParticleGun->GeneratePrimaryVertex(anEvent);
    //   }
    // }
  }
  if (fRun->getRunType() == 7) { // calibration light
    double shift = fRun->getTest3();

    fParticleGun->SetParticlePosition(
      G4ThreeVector(-fRadiatorL / 2. + 170 - shift, 0, tan(33 * deg) * shift - 30));
    double angle = -G4UniformRand() * M_PI;
    G4ThreeVector vec(0, 0, 1);
    vec.setTheta(acos(G4UniformRand()));
    vec.setPhi(2 * M_PI * G4UniformRand());

    vec.rotateY(-M_PI / 2.);
    fParticleGun->SetParticleMomentumDirection(vec);
  }
  if (fRun->getRunType() == 6) { // for determining focal plane of the lens
    double shiftx = -fRadiatorL / 2. + 0.1;
    double shifty = fRadiatorW / 2. - G4UniformRand() * fRadiatorW;
    double shiftz = fRadiatorH / 2. - G4UniformRand() * fRadiatorH;

    double angle = 0.7 * (M_PI / 2. - G4UniformRand() * M_PI);
    G4ThreeVector vec(0, 0, 1);
    vec.setTheta(angle);
    /// vec.setTheta(acos(G4UniformRand()));
    // vec.setPhi(2*M_PI*G4UniformRand());
    // std::cout<<"angle "<<angle*180/M_PI <<std::endl;

    double lensThickness = 15;
    double separation = fRun->getBeamSize();
    if (separation < 0.001) separation = 10;
    double rotShiftX =
      0.5 * separation * std::cos(angle) + (0.5 * lensThickness + 0.1) * std::tan(angle);
    double rotShiftY = -0.5 * fRadiatorL + 0.1;

    // fParticleGun->SetParticlePosition(G4ThreeVector(shiftx,shifty,shiftz));
    // fParticleGun->SetParticlePosition(G4ThreeVector(rotShiftY, 0,-rotShiftX));
    fParticleGun->SetParticlePosition(G4ThreeVector(rotShiftY, 0, 0.5 * separation));

    vec.rotateY(-M_PI / 2.);
    fParticleGun->SetParticleMomentumDirection(vec);

    fParticleGun->GeneratePrimaryVertex(anEvent);
    rotShiftX = -0.5 * separation * std::cos(angle) + (0.5 * lensThickness + 0.1) * std::tan(angle);
    shiftx = -fRadiatorL / 2. + 0.1;
    shifty = fRadiatorW / 2. - G4UniformRand() * fRadiatorW;
    shiftz = fRadiatorH / 2. - G4UniformRand() * fRadiatorH;
    // fParticleGun->SetParticlePosition(G4ThreeVector(shiftx,shifty,shiftz));
    // fParticleGun->SetParticlePosition(G4ThreeVector(rotShiftY,0,-rotShiftX));
    fParticleGun->SetParticlePosition(G4ThreeVector(rotShiftY, 0, -0.5 * separation));
  }

  fParticleGun->GeneratePrimaryVertex(anEvent);

  G4ThreeVector dir = fParticleGun->GetParticleMomentumDirection();
  dir *= fParticleGun->GetParticleMomentum() * 0.001; // GeV/c

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
