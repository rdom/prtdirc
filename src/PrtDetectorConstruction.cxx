
#include "PrtDetectorConstruction.h"

#include "G4Material.hh"
#include "G4SDManager.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Trap.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4SubtractionSolid.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "PrtManager.h"
#include "PrtTriggerSD.h"
#include "PrtBarSD.h"
#include "PrtPrizmSD.h"
#include "PrtPixelSD.h"

PrtDetectorConstruction::PrtDetectorConstruction() : G4VUserDetectorConstruction() {

  fRun = PrtManager::Instance()->getRun();
  fRunType = fRun->getRunType();
  fGeomId = fRun->getGeometry();
  fMcpLayout = fRun->getPmtLayout();
  fLensId = fRun->getLens();
  fRadiatorId = fRun->getRadiator();
  fPrismStepX = fRun->getPrismStepX();
  fTest1 = fRun->getTest1();
  fTest2 = fRun->getTest2();
  fTest3 = fRun->getTest3();

  fNRow = 3;
  fNCol = 5;

  fHall[0] = 1500;
  fHall[1] = 500;
  fHall[2] = 6000;

  // fBar[0] = 17; fBar[1] = 32; fBar[2] =1250;
  fBar[0] = 17.15;
  fBar[1] = 34.93;
  fBar[2] = 1200.06; // InSync #3
  // fBar[0] = 17.1; fBar[1] = 35.0; fBar[2] =1224.0; // InSync #0
  // fBar[0] = 17.9; fBar[1] = 35.9; fBar[2] =1200.7; // Zygo //2017

  fMirror[0] = 20;
  fMirror[1] = 40;
  fMirror[2] = 1;
  fPrizm[0] = 170;
  fPrizm[1] = 300;
  fPrizm[3] = 50;
  fPrizm[2] = fPrizm[3] + fPrizm[1] * tan(45 * deg);
  fMcpTotal[0] = fMcpTotal[1] = 53 + 6;
  fMcpTotal[2] = 1;
  fMcpActive[0] = fMcpActive[1] = 53;
  fMcpActive[2] = 1;
  fLens[0] = fLens[1] = 40;
  fLens[2] = 10;
  if (fLensId == 2) {
    fLens[0] = 50;
    fLens[1] = 175;
    fLens[2] = 14.4;
  } else if (fLensId == 3) {
    fLens[0] = 60;
    fLens[1] = 60;
    fLens[2] = 15;
  } else if (fLensId == 4 || fLensId == 5) {
    fLens[0] = 50;
    fLens[1] = 50;
    fLens[2] = 5.7;
  } else if (fLensId == 6 || fLensId == 7 || fLensId == 8) {
    fLens[0] = 50;
    fLens[1] = 175;
    fLens[2] = 12;
  }

  if (fMcpLayout == 2012) {
    fNCol = 3;
    fPrizm[2] = 30 + fPrizm[1] * tan(30 * deg);
    fPrizm[3] = 30;
  }

  if (fMcpLayout == 2015) {
  }

  if (fMcpLayout == 2016) {
    fNCol = 3;
    fPrizm[2] = 30 + fPrizm[1] * tan(30 * deg);
    fPrizm[3] = 30;
  }

  if (fMcpLayout == 2017 || fMcpLayout == 2030) {
    fNCol = 4;
  }

  if (fMcpLayout == 2018) {
    fNRow = 2;
    fNCol = 4;
  }

  if (fMcpLayout == 20171) {
    fNRow = 2;
    fNCol = 4;
  }

  if (fMcpLayout == 2021) { // Barrel DIRC layout
    fNCol = 4;
    fPrizm[2] = 40 + fPrizm[1] * tan(33 * deg);
    fPrizm[3] = 40;
  }

  if (fRadiatorId == 2) {
    fBar[0] = 17.1;
    fBar[1] = 174.8;
    fBar[2] = 1224.9;
    fMirror[1] = 180;
  }

  if (fRadiatorId == 11) {
    fBar[0] = 17.1;
    fBar[1] = 35;
    fBar[2] = 1224.9;
    fMirror[1] = 180;
  }

  if (fRadiatorId == 20) {
    fBar[0] = 17.1;
    fBar[1] = fTest1;
    fBar[2] = 1224.9;
    fMirror[1] = 180;
  }

  if (fRadiatorId == 5) {
    fBar[0] = 17.1;
    fBar[1] = 70;
    fBar[2] = 1224.9;
    fMirror[1] = 180;
  }

  // X configuration
  double radiatorStepY = fRun->getPrismStepY();
  if (radiatorStepY != 0) radiatorStepY += fBar[0] / 2. - fPrizm[3] / 2.;
  double radiatorStepX = -(fBar[1] - fPrizm[0]) / 2. - fPrismStepX;

  fCenterShift = G4ThreeVector(0., 0., 0.);

  if (fGeomId == 2015) {
    // fPrismRadiatorStep = fPrizm[3]/2.-fBar[0]/2.;
    // fPrismRadiatorStep = fBar[0]/2.-fPrizm[3]/2.   +30.6;
    // fCenterShift =  G4ThreeVector(fBar[2]/2.,0,-132);

    fCenterShift = G4ThreeVector(0.5 * fBar[2] - 96, -0.5 * fPrizm[0] + fRun->getBeamX(), -122);
  }

  if (fGeomId == 2016) {
    fPrizm[0] = 170;
    fPrizm[1] = 300;
    fPrizm[3] = 30;
    fPrizm[2] = fPrizm[3] + fPrizm[1] * tan(30 * deg);
    fCenterShift = G4ThreeVector(0.5 * fBar[2] - 96, -0.5 * fPrizm[0] + fRun->getBeamX(), -(279 - 187.5 - fBar[0]));
  }

  fOffset = 0;
  if (fGeomId == 2017) {
    fOffset = 146;
    fPrizm[0] = 175;
    fPrizm[1] = 300;
    fPrizm[3] = 50;
    fPrizm[2] = fPrizm[3] + fPrizm[1] * tan(33 * deg);
    fCenterShift = G4ThreeVector(0.5 * fBar[2] - fOffset, -0.5 * fPrizm[0] + fRun->getBeamX(), -100);
  }

  if (fGeomId == 2018) {
    fOffset = 146;
    fPrizm[0] = 175;
    fPrizm[1] = 300;
    fPrizm[3] = 50;
    fPrizm[2] = fPrizm[3] + fPrizm[1] * tan(33 * deg);
    fCenterShift = G4ThreeVector(0.5 * fBar[2] - fOffset, -0.5 * fPrizm[0] + fRun->getBeamX(), -100);
  }

  if (fGeomId == 2021) {
    fCenterShift = G4ThreeVector(0.5 * fBar[2] - 96, -0.5 * fPrizm[0] + fRun->getBeamX(), -(279 - 187.5 - fBar[0]));
  }

  if (fRunType == 1) {
    fCenterShift = G4ThreeVector(0., 0., 0.);
  }

  if (fRunType == 6) { // focal plane scan
    fPrizm[1] = 60;
    fPrizm[3] = 500;
    fPrizm[2] = 500;
    radiatorStepX = 0;
  }

  if (fGeomId > 2014 && (fRadiatorId > 1)) radiatorStepX = -(fBar[1] - fPrizm[0]) / 2.;
  if (fRadiatorId > 2) radiatorStepX = 0;
  
  fRun->setRadiatorL(fBar[2]);
  fRun->setRadiatorW(fBar[1]);
  fRun->setRadiatorH(fBar[0]);
  fRun->setPrismStepY(radiatorStepY);
  fRun->setPrismStepX(radiatorStepX);

  fPrtRot = new G4RotationMatrix();
  // create a messenger for this class
  fGeomMessenger = new PrtDetectorConstructionMessenger(this);
  fRotAngle = 0;
}

PrtDetectorConstruction::~PrtDetectorConstruction() {}

G4VPhysicalVolume *PrtDetectorConstruction::Construct() {
  DefineMaterials();

  // ------------- Volumes --------------

  // The experimental Hall
  G4Box *gExpHall = new G4Box("gExpHall", fHall[0], fHall[1], fHall[2]);
  lExpHall = new G4LogicalVolume(gExpHall, defaultMaterial, "lExpHall", 0, 0, 0);
  double zshift = (fRun->getBeamZ() == -1) ? 0 : fRun->getBeamZ() - fOffset;
  double radiatorStepY = fRun->getPrismStepY();
  double radiatorStepX = fRun->getPrismStepX();

  G4VPhysicalVolume *wExpHall = new G4PVPlacement(0, G4ThreeVector(), lExpHall, "gExpHall", 0, false, 0);

  // The Trigger and The front material
  G4Box *gFront = new G4Box("gFront", 200., 200., 5);
  lFront = new G4LogicalVolume(gFront, frontMaterial, "lFront", 0, 0, 0);
  G4Box *gTrigger = new G4Box("gTrigger", 20., 20., 5);
  lTrigger = new G4LogicalVolume(gTrigger, frontMaterial, "lTrigger", 0, 0, 0);
  G4Box *gEdd = new G4Box("gEdd", 200, 200, 23);
  lEdd = new G4LogicalVolume(gEdd, BarMaterial, "lEdd", 0, 0, 0);

  if (fGeomId == 3 || fGeomId >= 2015) {
    new G4PVPlacement(0, G4ThreeVector(0, 0, -4500), lFront, "wFront", lExpHall, false, 0);
    // new G4PVPlacement(0,G4ThreeVector(0,0,-5500),lEdd,"wEdd",lExpHall,false,0);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 1500), lTrigger, "wTrigger", lExpHall, false, 0);
  }

  // The DIRC
  G4Box *gDirc = new G4Box("gDirc", 400., 200., fBar[2] / 2. + fPrizm[1] + 50);
  lDirc = new G4LogicalVolume(gDirc, defaultMaterial, "lDirc", 0, 0, 0);

  G4ThreeVector dircpos = G4ThreeVector(0., 0., 0.);
  if (fCenterShift.mag() != 0) {
    dircpos = fCenterShift; // G4ThreeVector(fCenterShift, 0., 0.);
    dircpos.rotateY((fRun->getTheta() - 90) * deg);
  }
  if (fGeomId == 0 || fRunType == 1) zshift = 0;
  // tilt scan
  fPrtRot->rotateY((90 - fRun->getTheta()) * deg);
  fPrtRot->rotateX(fRun->getPhi() * deg);
  wDirc = new G4PVPlacement(fPrtRot, dircpos + G4ThreeVector(-zshift, 0, 0), lDirc, "wDirc", lExpHall, false, 0);

  // The DIRC cover box
  G4Box *gCover = new G4Box("gCover", 5, 150, fBar[2] / 2.);
  lCover = new G4LogicalVolume(gCover, MirrorMaterial, "lCover", 0, 0, 0);

  if (fGeomId == 3 || fGeomId >= 2015) {
    new G4PVPlacement(0, G4ThreeVector(-100, 0, 0), lCover, "wCover", lDirc, false, 0);
  }

  // The Bar
  G4Box *gBar = new G4Box("gBar", fBar[0] / 2., fBar[1] / 2., fBar[2] / 2.);

  if (fRunType == 6) lBar = new G4LogicalVolume(gBar, OilMaterial, "lBar", 0, 0, 0);
  else lBar = new G4LogicalVolume(gBar, BarMaterial, "lBar", 0, 0, 0);

  wBar = new G4PVPlacement(0, G4ThreeVector(radiatorStepY, radiatorStepX, 0), lBar, "wBar", lDirc, false, 0);

  // radiator covered with grease
  double greased = 0 * mm;
  if (fGeomId < 2017 && fRunType != 6) {
    greased = 1.5 * mm;
    if (fLensId == 0) greased = 0.5 * mm;
    G4Box *gOpticalGreased = new G4Box("gOpticalgreased", 0.5 * fBar[0], 0.5 * fBar[1], 0.5 * greased);
    lOpticalGreased = new G4LogicalVolume(gOpticalGreased, BarMaterial, "lOpticalGreased", 0, 0, 0);
    new G4PVPlacement(0, G4ThreeVector(radiatorStepY, radiatorStepX, 0.5 * fBar[2] + 0.5 * greased), lOpticalGreased, "wOpticalGreased", lDirc, false, 0);
  }

  // Optical grease
  if (fRunType != 6) {
    double greasew = 0.1 * mm;
    if (fLensId == 0) greasew = 0.1 * mm;
    G4Box *gOpticalGrease = new G4Box("gOpticalgrease", 0.5 * fBar[0], 0.5 * fBar[1], 0.5 * greasew);
    lOpticalGrease = new G4LogicalVolume(gOpticalGrease, opticalGreaseMaterial, "lOpticalGrease", 0, 0, 0);
    new G4PVPlacement(0, G4ThreeVector(radiatorStepY, radiatorStepX, 0.5 * fBar[2] + greased + 0.5 * greasew), lOpticalGrease, "wOpticalGrease", lDirc, false, 0);
    greased += greasew;
  }

  if (fRun->getStudy() == 430) { // add cookies
    double cookiew = 2 * mm;
    G4Box *gCookie1 = new G4Box("gCookie1", 0.9 * fBar[0], 0.9 * fBar[1], 0.5 * cookiew);
    lCookie1 = new G4LogicalVolume(gCookie1, opticalCookieMaterial, "lCookie1", 0, 0, 0);
    new G4PVPlacement(0, G4ThreeVector(radiatorStepY, radiatorStepX, 0.5 * fBar[2] + greased + 0.5 * cookiew), lCookie1, "wCookie1", lDirc, false, 0);
    greased += cookiew;
  }

  // // The Mirror gap
  // double mirrorgap=0.1*mm;
  // G4Box* gMirrorGap = new G4Box("gMirrorGap",fMirror[0]/2.,fMirror[1]/2.,0.5*mirrorgap);
  // lMirrorGap = new G4LogicalVolume(gMirrorGap,defaultMaterial,"lMirrorGap",0,0,0);
  // wMirrorGap =new G4PVPlacement(0,G4ThreeVector(radiatorStepY,radiatorStepX,-fBar[2]/2.-0.5*mirrorgap),lMirrorGap,"wMirrorGap", lDirc,false,0);

  // The Mirror
  G4Box *gMirror = new G4Box("gMirror", fMirror[0] / 2., fMirror[1] / 2., fMirror[2] / 2.);
  lMirror = new G4LogicalVolume(gMirror, MirrorMaterial, "lMirror", 0, 0, 0);
  wMirror = new G4PVPlacement(0, G4ThreeVector(radiatorStepY, radiatorStepX, -fBar[2] / 2. - fMirror[2] / 2.), lMirror, "wMirror", lDirc, false, 0); //-mirrorgap

  // The Lens
  G4Box *gfbox = new G4Box("Fbox", fLens[0] / 2., fLens[1] / 2., fLens[2] / 2.);

  if (fLensId == 1) { // 2-layer spherical lens
    double r1 = 0;    // fTest1
    double lensrad1 = (r1 == 0) ? 73.58 : r1;
    double lensMinThikness = 2;

    G4ThreeVector zTrans1(0, 0, -lensrad1 + fLens[2] / 2. - lensMinThikness);
    G4Tubs *gftub = new G4Tubs("Ftub", 0, fLens[0] / 2., fLens[2] / 2., 0. * deg, 360. * deg);
    G4Sphere *gsphere = new G4Sphere("Sphere", 0, lensrad1, 0, 360 * deg, 0, 360 * deg);
    G4IntersectionSolid *gLens1 = new G4IntersectionSolid("Ftub*Sphere", gftub, gsphere, new G4RotationMatrix(), zTrans1);
    G4SubtractionSolid *gLens2 = new G4SubtractionSolid("Ftub-Sphere", gftub, gsphere, new G4RotationMatrix(), zTrans1);
    lLens1 = new G4LogicalVolume(gLens1, Nlak33aMaterial, "lLens1", 0, 0, 0); // Nlak33aMaterial
    lLens2 = new G4LogicalVolume(gLens2, BarMaterial, "lLens2", 0, 0, 0);
  }

  if (fLensId == 2) { // 2-layer cylindrical lens
    double lensrad = 73.58;
    double lensMinThikness = 8;
    G4ThreeVector zTrans(0, 0, -lensrad + fLens[2] / 2. - lensMinThikness);
    G4Tubs *gcylinder = new G4Tubs("Cylinder", 0, lensrad, fLens[1] / 2. + 1.1, 0, 360 * deg);
    G4RotationMatrix *xRot = new G4RotationMatrix();
    xRot->rotateX(M_PI / 2. * rad);
    G4IntersectionSolid *gLens1 = new G4IntersectionSolid("Fbox*Cylinder", gfbox, gcylinder, xRot, zTrans);
    G4SubtractionSolid *gLens2 = new G4SubtractionSolid("Fbox-Cylinder", gfbox, gcylinder, xRot, zTrans);
    lLens1 = new G4LogicalVolume(gLens1, Nlak33aMaterial, "lLens1", 0, 0, 0); // Nlak33aMaterial
    lLens2 = new G4LogicalVolume(gLens2, BarMaterial, "lLens2", 0, 0, 0);
  }

  if (fLensId == 334) { // 3-component spherical lens
    double lensMinThikness = 2;

    double r1 = 0; // fTest1;
    double r2 = 0; // fTest2;

    double lensrad1 = (r1 == 0) ? 47.8 : r1;
    double lensrad2 = (r2 == 0) ? 29.1 : r2;

    G4ThreeVector zTrans1(0, 0, -lensrad1 - fLens[2] / 2. + lensrad1 - sqrt(lensrad1 * lensrad1 - fLens[0] / 2. * fLens[0] / 2.) + lensMinThikness);
    G4ThreeVector zTrans2(0, 0, -lensrad2 + fLens[2] / 2. - lensMinThikness);

    G4Sphere *gsphere1 = new G4Sphere("Sphere1", 0, lensrad1, 0, 360 * deg, 0, 360 * deg);
    G4Sphere *gsphere2 = new G4Sphere("Sphere2", 0, lensrad2, 0, 360 * deg, 0, 360 * deg);

    G4IntersectionSolid *gLens1 = new G4IntersectionSolid("Fbox*Sphere1", gfbox, gsphere1, new G4RotationMatrix(), zTrans1);
    G4SubtractionSolid *gLenst = new G4SubtractionSolid("Fbox-Sphere1", gfbox, gsphere1, new G4RotationMatrix(), zTrans1);

    G4IntersectionSolid *gLens2 = new G4IntersectionSolid("gLenst*Sphere2", gLenst, gsphere2, new G4RotationMatrix(), zTrans2);
    G4SubtractionSolid *gLens3 = new G4SubtractionSolid("gLenst-Sphere2", gLenst, gsphere2, new G4RotationMatrix(), zTrans2);

    lLens1 = new G4LogicalVolume(gLens1, BarMaterial, "lLens1", 0, 0, 0);
    lLens2 = new G4LogicalVolume(gLens2, Nlak33aMaterial, "lLens2", 0, 0, 0);
    lLens3 = new G4LogicalVolume(gLens3, BarMaterial, "lLens3", 0, 0, 0);
  }

  if (fLensId == 3) { // 3-component spherical lens
    double lensMinThikness = 2.0;

    double r1 = fTest1;
    double r2 = fTest2;

    if (fRunType == 6) { // focal plane scan
      r1 = fTest1;
      r2 = fTest2;
    }

    r1 = (r1 == 0) ? 47.80 : r1;
    r2 = (r2 == 0) ? 29.12 : r2;
    double shight = 40;
    double bwidth = fLens[2] - lensMinThikness * 2;

    G4ThreeVector zTrans1(0, 0, -r1 - fLens[2] / 2. + r1 - sqrt(r1 * r1 - shight / 2. * shight / 2.) + lensMinThikness);
    G4ThreeVector zTrans2(0, 0, -r2 - fLens[2] / 2. + r2 - sqrt(r2 * r2 - shight / 2. * shight / 2.) + lensMinThikness * 2);

    G4Box *gfbox0 = new G4Box("Fbox0", 0.5 * fLens[0] + 1, 0.5 * fLens[1] + 1, 0.5 * fLens[2]);
    G4Tubs *gftub = new G4Tubs("Ftub", 0, 0.5 * fLens[0], 0.5 * fLens[2], 0. * deg, 360. * deg);
    G4Box *gfsbox = new G4Box("Fsbox", 0.5 * shight, 0.5 * fLens[1], 0.5 * fLens[2]);
    G4Tubs *gfstube = new G4Tubs("ftube", 0, 0.5 * shight, 0.5 * fLens[2], 0. * deg, 360. * deg);

    G4Sphere *gsphere1 = new G4Sphere("Sphere1", 0, r1, 0, 360 * deg, 0, 360 * deg);
    G4Sphere *gsphere2 = new G4Sphere("Sphere2", 0, r2, 0, 360 * deg, 0, 360 * deg);

    G4IntersectionSolid *gbbox = new G4IntersectionSolid("bbox", gftub, gfbox0, new G4RotationMatrix(), G4ThreeVector(0, 0, lensMinThikness * 2));
    G4IntersectionSolid *gsbox = new G4IntersectionSolid("sbox", gfstube, gfbox0, new G4RotationMatrix(), G4ThreeVector(0, 0, -lensMinThikness * 2));

    G4UnionSolid *gubox = new G4UnionSolid("unionbox", gbbox, gsbox, new G4RotationMatrix(), G4ThreeVector(0, 0, 0));

    G4IntersectionSolid *gLens1 = new G4IntersectionSolid("Lens1", gubox, gsphere1, new G4RotationMatrix(), zTrans1);
    G4SubtractionSolid *gLenst = new G4SubtractionSolid("temp", gubox, gsphere1, new G4RotationMatrix(), zTrans1);

    G4IntersectionSolid *gLens2 = new G4IntersectionSolid("Lens2", gLenst, gsphere2, new G4RotationMatrix(), zTrans2);
    G4SubtractionSolid *gLens3 = new G4SubtractionSolid("Lens3", gLenst, gsphere2, new G4RotationMatrix(), zTrans2);

    lLens1 = new G4LogicalVolume(gLens1, BarMaterial, "lLens1", 0, 0, 0);
    if (fTest3 > 1)
      lLens2 = new G4LogicalVolume(gLens2, PbF2Material, "lLens2", 0, 0, 0);
    else
      lLens2 = new G4LogicalVolume(gLens2, Nlak33aMaterial, "lLens2", 0, 0, 0);
    lLens3 = new G4LogicalVolume(gLens3, BarMaterial, "lLens3", 0, 0, 0);
  }

  if (fLensId == 4) { // Spherical lens with air gap // f =250 , d = , w = 5.7
    double r1 = 0;    // fTest1;
    double lensrad1 = (r1 == 0) ? 250 : r1;
    double lensMinThikness = 2;

    G4ThreeVector zTrans1(0, 0, -lensrad1 + fLens[2] / 2. - lensMinThikness);
    G4Tubs *gftub = new G4Tubs("Ftub", 0, fLens[0] / 2., fLens[2] / 2., 0. * deg, 360. * deg);
    G4Sphere *gsphere = new G4Sphere("Sphere", 0, lensrad1, 0, 360 * deg, 0, 360 * deg);
    G4IntersectionSolid *gLens1 = new G4IntersectionSolid("Ftub*Sphere", gftub, gsphere, new G4RotationMatrix(), zTrans1);
    lLens1 = new G4LogicalVolume(gLens1, BarMaterial, "lLens1", 0, 0, 0); // Nlak33aMaterial
  }

  if (fLensId == 5) { // Spherical lens with air gap // f =250 , d = , w = 5.7 // black edges
    double r1 = 0;    // fTest1
    double lensrad1 = (r1 == 0) ? 250 : r1;
    double lensMinThikness = 2;

    G4ThreeVector zTrans1(0, 0, -lensrad1 + fLens[2] / 2. - lensMinThikness);
    G4Tubs *gftub = new G4Tubs("Ftub", 0, fLens[0] / 2., fLens[2] / 2., 0. * deg, 360. * deg);
    G4Sphere *gsphere = new G4Sphere("Sphere", 0, lensrad1, 0, 360 * deg, 0, 360 * deg);
    G4IntersectionSolid *gLens1 = new G4IntersectionSolid("Ftub*Sphere", gftub, gsphere, new G4RotationMatrix(), zTrans1);
    G4SubtractionSolid *gLens2 = new G4SubtractionSolid("Ftub-Sphere", gftub, gsphere, new G4RotationMatrix(), zTrans1);
    lLens1 = new G4LogicalVolume(gLens1, BarMaterial, "lLens1", 0, 0, 0); // Nlak33aMaterial
    lLens2 = new G4LogicalVolume(gLens2, defaultMaterial, "lLens2", 0, 0, 0);
  }

  if (fLensId == 6) { // 3-component cylindrical lens
    double lensMinThikness = 2.0;

    double r1 = 0; // fTest1
    double r2 = 0; // fTest2

    if (fRunType == 6) { // focal plane scan
      r1 = fTest1;
      r2 = fTest2;
    }

    // RMI lens
    fLens[2] = 13.12;
    lensMinThikness = 2.51;
    double layer12 = lensMinThikness + 3.525; // lensMinThikness*2;

    // r1 = (r1==0)? 27.45: r1;
    // r2 = (r2==0)? 20.02: r2;

    r1 = (r1 == 0) ? 33 : r1;
    r2 = (r2 == 0) ? 24 : r2;
    double shight = 20;

    G4ThreeVector zTrans1(0, 0, -r1 - fLens[2] / 2. + r1 - sqrt(r1 * r1 - shight / 2. * shight / 2.) + lensMinThikness);
    G4ThreeVector zTrans2(0, 0, -r2 - fLens[2] / 2. + r2 - sqrt(r2 * r2 - shight / 2. * shight / 2.) + layer12);

    G4Box *gftub = new G4Box("ftub", 0.5 * fLens[0], 0.5 * fLens[1], 0.5 * fLens[2]);
    G4Box *gcbox = new G4Box("cbox", 0.5 * fLens[0], 0.5 * fLens[1] + 1, 0.5 * fLens[2]);
    G4ThreeVector tTrans1(0.5 * (fLens[0] + shight), 0, -fLens[2] + layer12);
    G4ThreeVector tTrans0(-0.5 * (fLens[0] + shight), 0, -fLens[2] + layer12);
    G4SubtractionSolid *tubox = new G4SubtractionSolid("tubox", gftub, gcbox, new G4RotationMatrix(), tTrans1);
    G4SubtractionSolid *gubox = new G4SubtractionSolid("gubox", tubox, gcbox, new G4RotationMatrix(), tTrans0);

    G4Tubs *gcylinder1 = new G4Tubs("Cylinder1", 0, r1, 0.5 * fLens[1], 0 * deg, 360 * deg);
    G4Tubs *gcylinder2 = new G4Tubs("Cylinder2", 0, r2, 0.5 * fLens[1] - 0.5, 0 * deg, 360 * deg);
    G4Tubs *gcylinder1c = new G4Tubs("Cylinder1c", 0, r1, 0.5 * fLens[1] + 0.5, 0 * deg, 360 * deg);
    G4Tubs *gcylinder2c = new G4Tubs("Cylinder2c", 0, r2, 0.5 * fLens[1] + 0.5, 0 * deg, 360 * deg);
    G4RotationMatrix *xRot = new G4RotationMatrix();
    xRot->rotateX(M_PI / 2. * rad);

    G4IntersectionSolid *gLens1 = new G4IntersectionSolid("Lens1", gubox, gcylinder1, xRot, zTrans1);
    G4SubtractionSolid *gLenst = new G4SubtractionSolid("temp", gubox, gcylinder1c, xRot, zTrans1);

    G4IntersectionSolid *gLens2 = new G4IntersectionSolid("Lens2", gLenst, gcylinder2, xRot, zTrans2);
    G4SubtractionSolid *gLens3 = new G4SubtractionSolid("Lens3", gLenst, gcylinder2c, xRot, zTrans2);

    lLens1 = new G4LogicalVolume(gLens1, BarMaterial, "lLens1", 0, 0, 0);
    lLens2 = new G4LogicalVolume(gLens2, Nlak33aMaterial, "lLens2", 0, 0, 0);
    lLens3 = new G4LogicalVolume(gLens3, BarMaterial, "lLens3", 0, 0, 0);
  }

  if (fLensId == 7) { // 3-component cylindrical lens
    double lensMinThikness = 2.0;

    double r1 = 0; // fTest1
    double r2 = 0; // fTest2

    if (fRunType == 6) { // focal plane scan
      r1 = fTest1;
      r2 = fTest2;
    }

    // r1 = (r1==0)? 27.45: r1;
    // r2 = (r2==0)? 20.02: r2;

    r1 = (r1 == 0) ? 33 : r1;
    r2 = (r2 == 0) ? 25 : r2;
    double shight = 20;

    G4ThreeVector zTrans1(0, 0, -r1 - fLens[2] / 2. + r1 - sqrt(r1 * r1 - shight / 2. * shight / 2.) + lensMinThikness);
    G4ThreeVector zTrans2(0, 0, -r2 - fLens[2] / 2. + r2 - sqrt(r2 * r2 - shight / 2. * shight / 2.) + lensMinThikness * 2);

    G4Box *gfboxx = new G4Box("fboxx", 0.5 * fLens[0], 0.5 * fLens[1], 0.5 * fLens[2]);
    G4Box *gcbox = new G4Box("cbox", 0.5 * fLens[0], 0.5 * fLens[1] + 1, 0.5 * fLens[2]);
    G4ThreeVector tTrans1(0.5 * (fLens[0] + shight), 0, -fLens[2] + lensMinThikness);
    G4ThreeVector tTrans0(-0.5 * (fLens[0] + shight), 0, -fLens[2] + lensMinThikness);
    G4SubtractionSolid *tubox = new G4SubtractionSolid("tubox", gfboxx, gcbox, new G4RotationMatrix(), tTrans1);
    G4SubtractionSolid *gubox = new G4SubtractionSolid("gubox", tubox, gcbox, new G4RotationMatrix(), tTrans0);

    G4Tubs *gcylinder1 = new G4Tubs("Cylinder1", 0, r1, 0.5 * fLens[1], 0 * deg, 360 * deg);
    G4Tubs *gcylinder2 = new G4Tubs("Cylinder2", 0, r2, 0.5 * fLens[1] - 0.5, 0 * deg, 360 * deg);
    G4Tubs *gcylinder1c = new G4Tubs("Cylinder1c", 0, r1, 0.5 * fLens[1] + 0.5, 0 * deg, 360 * deg);
    G4Tubs *gcylinder2c = new G4Tubs("Cylinder2c", 0, r2, 0.5 * fLens[1] + 0.5, 0 * deg, 360 * deg);
    G4RotationMatrix *xRot = new G4RotationMatrix();
    xRot->rotateX(M_PI / 2. * rad);

    G4IntersectionSolid *gLens1 = new G4IntersectionSolid("Lens1", gubox, gcylinder1, xRot, zTrans1);
    G4SubtractionSolid *gLenst = new G4SubtractionSolid("temp", gubox, gcylinder1c, xRot, zTrans1);

    G4IntersectionSolid *gLens2 = new G4IntersectionSolid("Lens2", gLenst, gcylinder2, xRot, zTrans2);
    G4SubtractionSolid *gLens3 = new G4SubtractionSolid("Lens3", gLenst, gcylinder2c, xRot, zTrans2);

    lLens1 = new G4LogicalVolume(gLens1, BarMaterial, "lLens1", 0, 0, 0);
    lLens2 = new G4LogicalVolume(gLens2, Nlak33aMaterial, "lLens2", 0, 0, 0);
    lLens3 = new G4LogicalVolume(gLens3, BarMaterial, "lLens3", 0, 0, 0);
  }

  if (fLensId == 8) { // 3-component cylindrical lens
    double lensMinThikness = 2.0;
    double r1 = 0; // fTest1
    double r2 = 0; // fTest2

    if (fRunType == 6) { // focal plane scan
      r1 = fTest1;
      r2 = fTest2;
    }

    // // thickness scan
    // double d = fTest1
    // d = (d==0)? 3: d;

    // r1 = (r1==0)? 27.45: r1;
    // r2 = (r2==0)? 20.02: r2;

    r1 = (r1 == 0) ? 33 : r1;
    r2 = (r2 == 0) ? 25 : r2;
    double shight = 0; // 19

    G4ThreeVector zTrans1(0, 0, -r1 - fLens[2] / 2. + r1 - sqrt(r1 * r1 - shight / 2. * shight / 2.) + 3.0);     // 1.5
    G4ThreeVector zTrans2(0, 0, -r2 - fLens[2] / 2. + r2 - sqrt(r2 * r2 - shight / 2. * shight / 2.) + 3.0 + 5); // 3.5

    G4Box *gfboxx = new G4Box("fboxx", 0.5 * fLens[0], 0.5 * fLens[1], 0.5 * fLens[2]);
    G4Box *gcbox = new G4Box("cbox", 0.5 * fLens[0], 0.5 * fLens[1] + 1, 0.5 * fLens[2]);
    G4ThreeVector tTrans1(0.5 * (fLens[0] + shight), 0, -fLens[2]);
    G4ThreeVector tTrans0(-0.5 * (fLens[0] + shight), 0, -fLens[2]);
    G4SubtractionSolid *tubox = new G4SubtractionSolid("tubox", gfboxx, gcbox, new G4RotationMatrix(), tTrans1);
    G4SubtractionSolid *gubox = new G4SubtractionSolid("gubox", tubox, gcbox, new G4RotationMatrix(), tTrans0);

    G4Tubs *gcylinder1 = new G4Tubs("Cylinder1", 0, r1, 0.5 * fLens[1], 0 * deg, 360 * deg);
    G4Tubs *gcylinder2 = new G4Tubs("Cylinder2", 0, r2, 0.5 * fLens[1] - 0.5, 0 * deg, 360 * deg);
    G4Tubs *gcylinder1c = new G4Tubs("Cylinder1c", 0, r1, 0.5 * fLens[1] + 0.5, 0 * deg, 360 * deg);
    G4Tubs *gcylinder2c = new G4Tubs("Cylinder2c", 0, r2, 0.5 * fLens[1] + 0.5, 0 * deg, 360 * deg);
    G4RotationMatrix *xRot = new G4RotationMatrix();
    xRot->rotateX(M_PI / 2. * rad);

    G4IntersectionSolid *gLens1 = new G4IntersectionSolid("Lens1", gubox, gcylinder1, xRot, zTrans1);
    G4SubtractionSolid *gLenst = new G4SubtractionSolid("temp", gubox, gcylinder1c, xRot, zTrans1);

    G4IntersectionSolid *gLens2 = new G4IntersectionSolid("Lens2", gLenst, gcylinder2, xRot, zTrans2);
    G4SubtractionSolid *gLens3 = new G4SubtractionSolid("Lens3", gLenst, gcylinder2c, xRot, zTrans2);

    lLens1 = new G4LogicalVolume(gLens1, BarMaterial, "lLens1", 0, 0, 0);
    lLens2 = new G4LogicalVolume(gLens2, Nlak33aMaterial, "lLens2", 0, 0, 0);
    lLens3 = new G4LogicalVolume(gLens3, BarMaterial, "lLens3", 0, 0, 0);
  }

  if (fLensId != 0 && fLensId != 10) {
    new G4PVPlacement(0, G4ThreeVector(radiatorStepY, radiatorStepX, 0.5 * fBar[2] + greased + 0.5 * fLens[2]), lLens1, "wLens1", lDirc, false, 0);
    if (fLensId != 4) new G4PVPlacement(0, G4ThreeVector(radiatorStepY, radiatorStepX, 0.5 * fBar[2] + greased + 0.5 * fLens[2]), lLens2, "wLens2", lDirc, false, 0);
    if (fLensId == 3 || fLensId == 6 || fLensId == 7 || fLensId == 8)
      new G4PVPlacement(0, G4ThreeVector(radiatorStepY, radiatorStepX, 0.5 * fBar[2] + greased + 0.5 * fLens[2]), lLens3, "wLens3", lDirc, false, 0);

  } else {
    fLens[2] = 0;
  }

  if (fRun->getStudy() == 430) {
    double cookiew = 2 * mm;
    G4Box *gCookie2 = new G4Box("gCookie2", 1.5 * fBar[0], 0.9 * fBar[1], 0.5 * cookiew);
    lCookie2 = new G4LogicalVolume(gCookie2, opticalCookieMaterial, "lCookie2", 0, 0, 0);
    new G4PVPlacement(0, G4ThreeVector(radiatorStepY, radiatorStepX, 0.5 * fBar[2] + fLens[2] + greased + 0.5 * cookiew), lCookie2, "wCookie1", lDirc, false, 0);
    greased += cookiew;
  }

  // The Prizm
  G4Trap *gPrizm = new G4Trap("gPrizm", fPrizm[0], fPrizm[1], fPrizm[2], fPrizm[3]);
  if (fRunType == 6)
    lPrizm = new G4LogicalVolume(gPrizm, OilMaterial, "lPrizm", 0, 0, 0);
  else
    lPrizm = new G4LogicalVolume(gPrizm, BarMaterial, "lPrizm", 0, 0, 0);

  G4RotationMatrix *xRot = new G4RotationMatrix();
  xRot->rotateX(M_PI / 2. * rad);
  fPrismShift = G4ThreeVector((fPrizm[2] + fPrizm[3]) / 4. - fPrizm[3] / 2., 0, 0.5 * fBar[2] + greased + 0.5 * fPrizm[1] + fLens[2]);
  fPrismShift = G4ThreeVector((fPrizm[2] + fPrizm[3]) / 4. - fPrizm[3] / 2., 0, 0.5 * fBar[2] + greased + 0.5 * fPrizm[1] + fLens[2]);
  new G4PVPlacement(xRot, fPrismShift, lPrizm, "wPrizm", lDirc, false, 0);

  if (fRunType == 7) { // calibration
    double height = 75;
    double cheight = height/tan(33 * deg);
    double cwidth = 75;
    double shift = fTest3;
    G4Trap *gPrizmC = new G4Trap("gPrizmC", cwidth, cheight, height, 0.00000000000001);
    lPrizmC = new G4LogicalVolume(gPrizmC, BarMaterial, "lPrizmC", 0, 0, 0);
    G4RotationMatrix *xRotC = new G4RotationMatrix();
    xRotC->rotateX(M_PI / 2.);
    xRotC->rotateZ(M_PI);
    G4ThreeVector fPrismShiftC = G4ThreeVector(height - height / 4. + fPrizm[3] / 2. + tan(33 * deg) * shift, 0, fBar[2] / 2. + cheight / 2. + fLens[2] + shift + 0.1);
    new G4PVPlacement(xRotC, fPrismShiftC, lPrizmC, "wPrizmC", lDirc, false, 0);
  }

  // Scaning plain
  if (false) {
    double shift = fTest3;
    if (shift < 935 - fPrizm[1]) {
      G4Box *gScan = new G4Box("gScan", 350, 100, 0.01);
      G4LogicalVolume *lScan = new G4LogicalVolume(gScan, BarMaterial, "lScan", 0, 0, 0);
      new G4PVPlacement(0, G4ThreeVector(0, 0, shift), lScan, "wScan", lBar, false, 0);
    } else {
      G4Box *gScan = new G4Box("gScan", 350, 0.01, 100);
      G4LogicalVolume *lScan = new G4LogicalVolume(gScan, BarMaterial, "lScan", 0, 0, 0);
      new G4PVPlacement(0, G4ThreeVector(0, 935 - shift - 150, 0), lScan, "wScan", lPrizm, false, 0);
    }
  }

  G4Box *gMcp;
  G4Box *gPixel;

  if (fMcpLayout > 1) {
    // The MCP
    gMcp = new G4Box("gMcp", fMcpTotal[0] / 2., fMcpTotal[1] / 2., fMcpTotal[2] / 2.);
    lMcp = new G4LogicalVolume(gMcp, BarMaterial, "lMcp", 0, 0, 0); // BarMaterial

    // The MCP Pixel
    int mcpDimx = 8;
    int mcpDimy = 8;
    if (fMcpLayout == 2030) {
      mcpDimx = 16;
      mcpDimy = 16;
    }

    if (fGeomId > 101 && fGeomId < 2000) {
      mcpDimx = fGeomId / 100;
      mcpDimy = fGeomId % 100;
    }
    gPixel = new G4Box("gPixel", fMcpActive[0] / (2 * (double)mcpDimx), fMcpActive[1] / (2 * (double)mcpDimy), fMcpActive[2] / 20.);
    lPixel = new G4LogicalVolume(gPixel, BarMaterial, "lPixel", 0, 0, 0);

    for (int i = 0; i < mcpDimx; i++) {
      for (int j = mcpDimy - 1; j >= 0; j--) {
        // if(i!=4 || j !=4 )continue; // for lut visualization
        double shiftx = i * (fMcpActive[0] / (double)mcpDimx) - fMcpActive[0] / 2. + fMcpActive[0] / (2 * (double)mcpDimx);
        double shifty = j * (fMcpActive[0] / (double)mcpDimy) - fMcpActive[0] / 2. + fMcpActive[0] / (2 * (double)mcpDimy);
        new G4PVPlacement(0, G4ThreeVector(shiftx, shifty, 0), lPixel, "wPixel", lMcp, false, mcpDimx * j + i);
      }
    }

    for (int j = 0; j < fNRow; j++) {
      for (int i = 0; i < fNCol; i++) {
        // if(i!=2 || j !=1 )continue; // for lut visualization
        double shiftx = i * (fMcpTotal[0] + 14) - fPrizm[3] / 2 + fMcpActive[0] / 2. + 2.5; // +2.5 adjustment to the prt2014
        if (j != 1) shiftx -= 14;                                                           //(3/2.)*fMcpActive[0]/8.;
        double shifty = (fMcpTotal[0] + 3) * (j - 1);

        // // cad version
        // double shiftx = i*(fMcpTotal[0]+14)-fPrizm[3]/2+fMcpActive[0]/2.;
        // if(j!=1) shiftx += 14;
        // shiftx -= 14;
        // double shifty = j*fMcpTotal[0]-fMcpTotal[0]+3*(j-1);

        if (fMcpLayout == 2012) {
          shiftx = i * (fMcpActive[0] + 2 + 6) - fBar[0] / 2. + fMcpActive[0] / 2. - 3 - 1;
          shifty = j * (fMcpActive[0] + 9 + 6) - fPrizm[0] / 2. + fMcpActive[0] / 2.; //-fPrizm[2]/2.-fPrizm[3]/2.
        }

        if (fMcpLayout == 2015) {
          double msh = 9; //(fPrizm[2]-5*fMcpTotal[0])/4.;
          shiftx = i * (fMcpTotal[0] + msh) - fPrizm[3] / 2 + fMcpActive[0] / 2.;
          // if(j==1) shiftx -= 3*fMcpActive[0]/8.;
          if (j == 1) shiftx -= (3 / 2.) * fMcpActive[0] / 8.; // i*(fMcpTotal[0]+3)-fPrizm[3]/2+fMcpActive[0]/2.+2*fMcpActive[0]/8.;
          shifty = (fMcpTotal[0] + 3) * (j - 1);
        }

        if (fMcpLayout == 2016) { // 9MCP
          double msh = 13;        //(fPrizm[2]-5*fMcpTotal[0])/4.;
          shiftx = 3 + i * (fMcpTotal[0] + msh) - fPrizm[3] / 2 + fMcpActive[0] / 2.;
          // if(j==1) shiftx -= 3*fMcpActive[0]/8.;
          if (j == 1) shiftx += (1 / 2.) * fMcpActive[0] / 8.; // i*(fMcpTotal[0]+3)-fPrizm[3]/2+fMcpActive[0]/2.+2*fMcpActive[0]/8.;
          shifty = (fMcpTotal[0] + 3) * (j - 1);
        }

        if (fMcpLayout == 201612) { // 12MCP
          double msh = 3;           //(fPrizm[2]-5*fMcpTotal[0])/4.;
          shiftx = 0 + i * (fMcpTotal[0] + msh) - fPrizm[3] / 2 + fMcpActive[0] / 2.;
          // if(j==1) shiftx -= 3*fMcpActive[0]/8.;
          if (j == 1) shiftx += (1 / 2.) * fMcpActive[0] / 8.;
          shifty = (fMcpTotal[0] + 3) * (j - 1);
        }

        if (fMcpLayout == 20171) {
          double msh = 11.2; //(fPrizm[2]-5*fMcpTotal[0])/4.;
          shiftx = i * (fMcpTotal[0] + msh) - fPrizm[3] / 2 + fMcpActive[0] / 2. + 3;
          if (j == 1) shiftx += (3 / 2.) * fMcpActive[0] / 8.;
          shifty = (fMcpTotal[0] + 3) * (j - 1);
        }

        if (fMcpLayout == 2017 || fMcpLayout == 2030) {
          double msh = 3;
          shiftx = i * (fMcpTotal[0] + msh) - fPrizm[3] / 2 + fMcpActive[0] / 2. + 3;
          shifty = (fMcpTotal[0] + 3) * (j - 1);
        }

        if (fMcpLayout == 2018) {
          double msh = 3;
          shiftx = i * (fMcpTotal[0] + msh) - fPrizm[3] / 2 + fMcpActive[0] / 2. + 3;
          // if(j==1) shiftx += (1/2.)*fMcpActive[0]/8.;
          shifty = (fMcpTotal[0] + 3) * (j - 1) + 0.5 * fMcpTotal[0] + 1.5;
        }

        if (fMcpLayout == 20171) {
          double msh = 3;
          shiftx = i * (fMcpTotal[0] + msh) - fPrizm[3] / 2 + fMcpActive[0] / 2. + 3;
          //	  if(j==1) shiftx += (1/2.)*fMcpActive[0]/8.;
          shifty = (fMcpTotal[0] + 3) * (j - 1) + fMcpActive[0] / 2. + 7;
        }

        if (fMcpLayout == 2021) {
          // double msh = 0; // true 2021
          double msh = 3; // 2021 layout for 2017 prototupe
          shiftx = i * (fMcpTotal[0] + msh) - fPrizm[3] / 2 + fMcpActive[0] / 2. + 3;

          // shifty = (fMcpTotal[0]+1)*(j-1);
          shifty = (fMcpTotal[0] + 3) * (j - 1); // 2021 layout for 2017 prototupe
          if (i == 0) {
            if (j == 0) continue;
            // shifty = (fMcpTotal[0]+1)*(j-1)-0.5*fMcpTotal[0]-0.5;
            shifty = (fMcpTotal[0] + 3) * (j - 1) - 0.5 * fMcpTotal[0] - 0.5; // 2021 layout for 2017 prototupe
          }
          new G4PVPlacement(0, G4ThreeVector(shiftx, shifty, 0.5 * fBar[2] + fPrizm[1] + 0.5 * fMcpActive[2] + greased + fLens[2]), lMcp, "wMcp", lDirc, false,
                            fNRow * i + (j - 1));
          continue;
        }

        new G4PVPlacement(0, G4ThreeVector(shiftx, shifty, 0.5 * fBar[2] + fPrizm[1] + 0.5 * fMcpActive[2] + greased + fLens[2]), lMcp, "wMcp", lDirc, false, fNRow * i + j);
        // new G4PVPlacement(0,G4ThreeVector(shiftx,shifty,0.5*fBar[2]+fPrizm[1]+0.5*fMcpActive[2]+greased+fLens[2]+0.1),lMcp,"wMcp", lDirc,false,fNRow*i+j);  //rd air gap
      }
    }
  } else {
    // for layout optimization
    // The MCP
    gMcp = new G4Box("gMcp", fPrizm[2] / 2., fPrizm[0] / 2., fMcpTotal[2] / 2.);
    lMcp = new G4LogicalVolume(gMcp, BarMaterial, "lMcp", 0, 0, 0);

    // The MCP Pixel
    if (fMcpLayout == 0) { // one prism-size mcp with one pixel
      gPixel = new G4Box("gPixel", fPrizm[2] / 2., fPrizm[0] / 2., fMcpActive[2] / 16.);
      lPixel = new G4LogicalVolume(gPixel, BarMaterial, "lPixel", 0, 0, 0);
      new G4PVPlacement(0, G4ThreeVector(0, 0, 0), lPixel, "wPixel", lMcp, false, 1);
    }

    if (fMcpLayout == 1) { // one prism-size mcp with many pixels
      int pixelId = 0;
      int mcpDimx = 80;
      int mcpDimy = 40;

      gPixel = new G4Box("gPixel", fPrizm[2] / (2 * (double)mcpDimx), fPrizm[0] / (2 * (double)mcpDimy), fMcpActive[2] / 20.);
      lPixel = new G4LogicalVolume(gPixel, BarMaterial, "lPixel", 0, 0, 0);
      for (int i = 0; i < mcpDimx; i++) {
        for (int j = 0; j < mcpDimy; j++) {
          double shiftx = i * (fPrizm[2] / (double)mcpDimx) - fPrizm[2] / 2. + fPrizm[2] / (2 * (double)mcpDimx);
          double shifty = j * (fPrizm[0] / (double)mcpDimy) - fPrizm[0] / 2. + fPrizm[0] / (2 * (double)mcpDimy);
          new G4PVPlacement(0, G4ThreeVector(shiftx, shifty, 0), lPixel, "wPixel", lMcp, false, pixelId++);
        }
      }
    }

    new G4PVPlacement(0, G4ThreeVector(fPrizm[2] / 2. - fPrizm[3] / 2., 0, fBar[2] / 2. + fPrizm[1] + fMcpActive[2] / 2. + fLens[2]), lMcp, "wMcp", lDirc, false, 1);
  }

  const int num = 36;
  double WaveLength[num];
  double PhotonEnergy[num];
  double PMTReflectivity[num];
  double EfficiencyMirrors[num];
  const double LambdaE = 2.0 * 3.14159265358979323846 * 1.973269602e-16 * m * GeV;
  for (int i = 0; i < num; i++) {
    WaveLength[i] = (300 + i * 10) * nanometer;
    PhotonEnergy[num - (i + 1)] = LambdaE / WaveLength[i];
    PMTReflectivity[i] = 0.;
    EfficiencyMirrors[i] = 0;
  }

  /***************** QUANTUM EFFICIENCY OF BURLE AND HAMAMTSU PMT'S *****/

  // ideal pmt quantum efficiency
  double QuantumEfficiencyIdial[num] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  // Burle PMT's
  double QuantumEfficiencyB[num] = {0.,   0.001, 0.002, 0.005, 0.01,  0.015, 0.02,  0.03, 0.04, 0.05, 0.06, 0.07, 0.09, 0.1,  0.13, 0.15, 0.17, 0.2,
                                    0.24, 0.26,  0.28,  0.282, 0.284, 0.286, 0.288, 0.29, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.15, 0.13, 0.12, 0.10};

  // hamamatsu pmt quantum efficiency
  double QuantumEfficiencyPMT[num] = {0.001, 0.002, 0.004, 0.007, 0.011, 0.015, 0.020, 0.026, 0.033, 0.040, 0.045, 0.056, 0.067, 0.085, 0.109, 0.129, 0.138, 0.147,
                                      0.158, 0.170, 0.181, 0.188, 0.196, 0.203, 0.206, 0.212, 0.218, 0.219, 0.225, 0.230, 0.228, 0.222, 0.217, 0.210, 0.199, 0.177};

  // these quantum efficiencies have to be multiplied by geometry
  //   efficiency of given PMT's
  //   for Hamamatsu by factor 0.7
  //   for Burle by factor 0.45
  for (int k = 0; k < 36; k++) {
    QuantumEfficiencyB[k] = QuantumEfficiencyB[k] * 0.45;
    QuantumEfficiencyPMT[k] = QuantumEfficiencyPMT[k] * .7;
  }

  // double QuantumEfficiency[num]=
  //    { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
  //      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
  //      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};

  //  double QuantumEfficiencyPMT[num] =
  //    { 1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
  //      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,
  //      1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};

  /* define quantum efficiency for burle PMT's - the same efficiency is
     assigned to pads and also to slots !!!! */

  // burle pmt - bigger slots => logicPad
  G4MaterialPropertiesTable *PhotocatodBurleMPT = new G4MaterialPropertiesTable();
  PhotocatodBurleMPT->AddProperty("EFFICIENCY", PhotonEnergy, QuantumEfficiencyB, num);
  PhotocatodBurleMPT->AddProperty("REFLECTIVITY", PhotonEnergy, PMTReflectivity, num);

  G4OpticalSurface *BurlePMTOpSurface = new G4OpticalSurface("BurlePMTOpSurface", glisur, polished, dielectric_metal);
  BurlePMTOpSurface->SetMaterialPropertiesTable(PhotocatodBurleMPT);

  // // assignment for pad
  // if(burle)
  //   new G4LogicalSkinSurface("BurlePMTSurface",logicBurPad,BurlePMTOpSurface);

  // if(burle1)
  //   new G4LogicalSkinSurface("Burle1PMTSurface",logicBur1Pad,BurlePMTOpSurface);

  /* hamamatsu pmt's - smaller slots => quantum efficiency again
     assign to slot and pad */

  fQuantumEfficiency = QuantumEfficiencyIdial; // QuantumEfficiencyPMT;//QuantumEfficiencyIdial;
  G4MaterialPropertiesTable *PhotocatodHamamatsuMPT = new G4MaterialPropertiesTable();
  PhotocatodHamamatsuMPT->AddProperty("EFFICIENCY", PhotonEnergy, fQuantumEfficiency, num);
  PhotocatodHamamatsuMPT->AddProperty("REFLECTIVITY", PhotonEnergy, PMTReflectivity, num);

  G4OpticalSurface *HamamatsuPMTOpSurface = new G4OpticalSurface("HamamatsuPMTOpSurface", glisur, polished, dielectric_metal);
  HamamatsuPMTOpSurface->SetMaterialPropertiesTable(PhotocatodHamamatsuMPT);

  // // assignment to pad
  // if(hamamatsu8500)
  new G4LogicalSkinSurface("HamamatsuPMTSurface", lPixel, HamamatsuPMTOpSurface);

  // Mirror
  G4OpticalSurface *MirrorOpSurface = new G4OpticalSurface("MirrorOpSurface", glisur, polished, dielectric_metal);

  // Mirror
  const int numPr = 44;
  const int numUv = 46;

  double mirrEnPr[] = {3.627, 3.602, 3.58,  3.557, 3.536, 3.514, 3.493, 3.471, 3.448, 3.423, 3.397, 3.37, 3.339, 3.301, 3.254, 3.195, 3.118, 3.032, 2.95,  2.873, 2.8,   2.73,
                       2.664, 2.601, 2.541, 2.483, 2.429, 2.376, 2.326, 2.278, 2.231, 2.187, 2.144, 2.11, 2.064, 2.019, 1.982, 1.947, 1.913, 1.88,  1.849, 1.818, 1.789, 1.772};

  double mirrEnUv[] = {5.571, 5.458, 5.349, 5.223, 5.07,  4.876, 4.669, 4.478, 4.303, 4.142, 3.991, 3.851, 3.72,  3.598, 3.485, 3.378,
                       3.277, 3.182, 3.092, 3.006, 2.927, 2.852, 2.78,  2.711, 2.646, 2.584, 2.524, 2.468, 2.413, 2.361, 2.312, 2.264,
                       2.219, 2.174, 2.132, 2.092, 2.053, 2.015, 1.979, 1.944, 1.91,  1.878, 1.846, 1.816, 1.785, 1.771};

  // protected aluminium
  double mirrReflPr[] = {0.5082, 0.5340, 0.5592, 0.5850, 0.6098, 0.6339, 0.6572, 0.6800, 0.7041, 0.7286, 0.7512, 0.7730, 0.7957, 0.8177, 0.8397,
                         0.8607, 0.8807, 0.8959, 0.9055, 0.9114, 0.9153, 0.9175, 0.9190, 0.9193, 0.9188, 0.9176, 0.9158, 0.9134, 0.9105, 0.9072,
                         0.9020, 0.8964, 0.8918, 0.8875, 0.8809, 0.8736, 0.8669, 0.8595, 0.8518, 0.8436, 0.8349, 0.8256, 0.8158, 0.8102};

  // UV enhanced aluminium mirror
  double mirrReflUv[] = {0.7736, 0.7975, 0.8190, 0.8401, 0.8621, 0.8836, 0.9011, 0.9116, 0.9210, 0.9261, 0.9275, 0.9282, 0.9283, 0.9270, 0.9254, 0.9235,
                         0.9210, 0.9178, 0.9144, 0.9108, 0.9068, 0.9027, 0.8984, 0.8940, 0.8896, 0.8850, 0.8803, 0.8755, 0.8709, 0.8664, 0.8621, 0.8578,
                         0.8536, 0.8498, 0.8462, 0.8426, 0.8393, 0.8362, 0.8329, 0.8299, 0.8271, 0.8244, 0.8216, 0.8188, 0.8158, 0.8146};

  // old mirror
  double ReflectivityMirrorBar[num] = {0.8755, 0.882, 0.889,  0.895,  0.9,    0.9025, 0.91, 0.913,  0.9165, 0.92,  0.923, 0.9245, 0.9285, 0.932, 0.934, 0.935, 0.936, 0.9385,
                                       0.9395, 0.94,  0.9405, 0.9405, 0.9405, 0.9405, 0.94, 0.9385, 0.936,  0.934, 0.931, 0.9295, 0.928,  0.928, 0.921, 0.92,  0.927, 0.9215};

  for (int i = 0; i < num; i++)
    ReflectivityMirrorBar[i] -= ReflectivityMirrorBar[i] * 0.08;

  G4MaterialPropertiesTable *MirrorMPT = new G4MaterialPropertiesTable();

  // MirrorMPT->AddProperty("REFLECTIVITY", mirrEnPr, mirrReflPr, numPr);
  // MirrorMPT->AddProperty("REFLECTIVITY", mirrEnUv, mirrReflUv, numUv);
  MirrorMPT->AddProperty("REFLECTIVITY", PhotonEnergy, ReflectivityMirrorBar, num);

  // MirrorMPT->AddProperty("EFFICIENCY", PhotonEnergy, EfficiencyMirrors,   num);

  MirrorOpSurface->SetMaterialPropertiesTable(MirrorMPT);
  new G4LogicalSkinSurface("MirrorSurface", lMirror, MirrorOpSurface);

  SetVisualization();

  return wExpHall;
}

void PrtDetectorConstruction::DefineMaterials() {
  G4String symbol;      // a=mass of a mole;
  double a, z, density; // z=mean number of protons;

  int ncomponents, natoms;
  double fractionmass;

  // define Elements
  G4Element *H = new G4Element("Hydrogen", symbol = "H", z = 1., a = 1.01 * g / mole);
  G4Element *C = new G4Element("Carbon", symbol = "C", z = 6., a = 12.01 * g / mole);
  G4Element *N = new G4Element("Nitrogen", symbol = "N", z = 7., a = 14.01 * g / mole);
  G4Element *O = new G4Element("Oxygen", symbol = "O", z = 8., a = 16.00 * g / mole);
  G4Element *Si = new G4Element("Silicon", symbol = "Si", z = 14., a = 28.09 * g / mole);

  G4Element *Al = new G4Element("Aluminum", symbol = "Al", z = 13., a = 26.98 * g / mole);

  // quartz material = SiO2
  G4Material *SiO2 = new G4Material("quartz", density = 2.200 * g / cm3, ncomponents = 2);
  SiO2->AddElement(Si, natoms = 1);
  SiO2->AddElement(O, natoms = 2);

  Nlak33aMaterial = new G4Material("Nlak33a", density = 4.220 * g / cm3, ncomponents = 2);
  Nlak33aMaterial->AddElement(Si, natoms = 1);
  Nlak33aMaterial->AddElement(O, natoms = 2);

  PbF2Material = new G4Material("PbF2", density = 4.220 * g / cm3, ncomponents = 2);
  PbF2Material->AddElement(Si, natoms = 1);
  PbF2Material->AddElement(O, natoms = 2);

  G4Material *Vacuum = new G4Material("interGalactic", 1., 1.008 * g / mole, 1.e-25 * g / cm3, kStateGas, 2.73 * kelvin, 3.e-18 * pascal);
  G4Material *Air = new G4Material("Air", density = 1.290 * mg / cm3, ncomponents = 2);
  Air->AddElement(N, fractionmass = 0.7);
  Air->AddElement(O, fractionmass = 0.3);

  G4Material *Aluminum = new G4Material("Aluminum", density = 2.7 * g / cm3, ncomponents = 1);
  Aluminum->AddElement(Al, fractionmass = 1.0);

  G4Material *KamLandOil = new G4Material("KamLandOil", density = 0.914 * g / cm3, ncomponents = 2);
  KamLandOil->AddElement(C, natoms = 12);
  KamLandOil->AddElement(H, natoms = 26);

  G4Material *CarbonFiber = new G4Material("CarbonFiber", density = 0.145 * g / cm3, ncomponents = 1);
  CarbonFiber->AddElement(C, fractionmass = 1.0);

  /* as I don't know the exact material composition,
     I will use Epoxyd material composition and add
     the optical property of Epotek to this material */

  G4Material *Epotek = new G4Material("Epotek", density = 1.2 * g / cm3, ncomponents = 3);

  Epotek->AddElement(C, natoms = 3);
  Epotek->AddElement(H, natoms = 5);
  Epotek->AddElement(O, natoms = 2);

  // assign main materials
  if (fGeomId == 1)
    defaultMaterial = Vacuum;
  else
    defaultMaterial = Air; // Vacuum // material of world
  frontMaterial = CarbonFiber;
  BarMaterial = SiO2;        // material of all Bars, Quartz and Window
  OilMaterial = KamLandOil;  // material of volume 1,2,3,4
  MirrorMaterial = Aluminum; // mirror material
  epotekMaterial = Epotek;   // Epotek material - glue between bars

  // ------------ Generate & Add Material Properties Table ------------

  static const double LambdaE = 2.0 * 3.14159265358979323846 * 1.973269602e-16 * m * GeV;
  const int num = 36;
  double WaveLength[num];
  double AirAbsorption[num];      // absorption value for air
  double AirRefractiveIndex[num]; // air refractive index
  double PhotonEnergy[num];       // energy of photons which correspond to the given
  // refractive or absoprtion values

  double PhotonEnergyNlak33a[76] = {1,       1.2511,  1.26386, 1.27687, 1.29016, 1.30372, 1.31758, 1.33173, 1.34619, 1.36097, 1.37607, 1.39152, 1.40731, 1.42347, 1.44,    1.45692,
                                    1.47425, 1.49199, 1.51016, 1.52878, 1.54787, 1.56744, 1.58751, 1.6081,  1.62923, 1.65092, 1.6732,  1.69609, 1.71961, 1.7438,  1.76868, 1.79427,
                                    1.82062, 1.84775, 1.87571, 1.90452, 1.93423, 1.96488, 1.99652, 2.0292,  2.06296, 2.09787, 2.13398, 2.17135, 2.21006, 2.25017, 2.29176, 2.33492,
                                    2.37973, 2.42631, 2.47473, 2.52514, 2.57763, 2.63236, 2.68946, 2.7491,  2.81143, 2.87666, 2.94499, 3.01665, 3.09187, 3.17095, 3.25418, 3.34189,
                                    3.43446, 3.53231, 3.6359,  3.74575, 3.86244, 3.98663, 4.11908, 4.26062, 4.41225, 4.57506, 4.75035, 4.93961};

  /*************************** ABSORPTION COEFFICIENTS *****************************/

  // absorption of KamLandOil per 50 cm - from jjv
  double KamLandOilAbsorption[num] = {0.97469022,  0.976603956, 0.978511548, 0.980400538, 0.982258449, 0.984072792, 0.985831062, 0.987520743, 0.989129303,
                                      0.990644203, 0.992052894, 0.993342822, 0.994501428, 0.995516151, 0.996374433, 0.997063719, 0.997571464, 0.997885132,
                                      0.997992205, 0.997880183, 0.997536591, 0.99,        0.98,        0.97,        0.96,        0.94,        0.93,
                                      0.924507,    0.89982,     0.883299,    0.85657,     0.842637,    0.77020213,  0.65727,     0.324022,    0.019192};

  // absorption of quartz per 1 m - from jjv
  double QuartzAbsorption[num] = {0.999572036, 0.999544661, 0.999515062, 0.999483019, 0.999448285, 0.999410586, 0.999369611, 0.999325013, 0.999276402,
                                  0.999223336, 0.999165317, 0.999101778, 0.999032079, 0.998955488, 0.998871172, 0.998778177, 0.99867541,  0.998561611,
                                  0.998435332, 0.998294892, 0.998138345, 0.997963425, 0.997767484, 0.997547418, 0.99729958,  0.99701966,  0.99670255,
                                  0.996342167, 0.995931242, 0.995461041, 0.994921022, 0.994298396, 0.993577567, 0.992739402, 0.991760297, 0.990610945};

  // absorption of epotek per one layer - thicknes 0.001'' - from jjv
  double EpotekAbsorption[num] = {0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999,
                                  0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999, 0.99999999,
                                  0.99999999, 0.9999,     0.9998,     0.9995,     0.999,      0.998,      0.997,      0.996,      0.9955,     0.993,      0.9871,     0.9745};

  // N-Lak 33a
  double Nlak33aAbsorption[76] = {371813,  352095,  331021,  310814,  291458,  272937,   255238,   238342,    222234,    206897,    192313,    178463,  165331,
                                  152896,  141140,  130043,  119585,  109747,  100507,   91846.3,  83743.1,   76176.7,   69126.1,   62570.2,   56488,   50858.3,
                                  45660.1, 40872.4, 36474.6, 32445.8, 28765.9, 25414.6,  22372.2,  19619.3,   17136.9,   14906.5,   12910.2,   11130.3, 9550.13,
                                  8153.3,  6924.25, 5848.04, 4910.46, 4098.04, 3398.06,  2798.54,  2288.32,   1856.99,   1494.92,   1193.28,   943.973, 739.657,
                                  573.715, 440.228, 333.94,  250.229, 185.064, 134.967,  96.9664,  68.5529,   47.6343,   32.4882,   21.7174,   14.2056, 9.07612,
                                  5.65267, 3.4241,  2.01226, 1.14403, 0.62722, 0.330414, 0.166558, 0.0799649, 0.0363677, 0.0155708, 0.00623089};

  // NLak33b from refractiveindex for 1 cm
  double Nlak33bEn[25] = {0.4959, 0.5332, 0.6293, 0.8103, 1.1696, 1.7712, 1.8785, 1.9997, 2.1376, 2.2707, 2.4796, 2.6953, 2.8436,
                          2.9520, 3.0613, 3.0996, 3.1790, 3.2627, 3.3509, 3.3968, 3.5424, 3.7121, 3.8745, 3.9994, 4.1328};
  double Nlak33bAb[25] = {0.398114, 0.679068, 0.937060, 0.985032, 0.997996, 0.997996, 0.997595, 0.997194, 0.997595, 0.997996, 0.997194, 0.994376, 0.991546,
                          0.988297, 0.982161, 0.979691, 0.971388, 0.954455, 0.928177, 0.910019, 0.820600, 0.657099, 0.455454, 0.245954, 0.158490};

  int n_PbF2 = 56;
  double en_PbF2[] = {1.55,  1.569, 1.59,  1.61,  1.631, 1.653, 1.675, 1.698, 1.722, 1.746, 1.771, 1.797, 1.823, 1.851, 1.879, 1.907, 1.937, 1.968, 2,
                      2.033, 2.066, 2.101, 2.138, 2.175, 2.214, 2.254, 2.296, 2.339, 2.384, 2.431, 2.48,  2.53,  2.583, 2.638, 2.695, 2.755, 2.818, 2.883,
                      2.952, 3.024, 3.1,   3.179, 3.263, 3.351, 3.444, 3.542, 3.647, 3.757, 3.875, 3.999, 4.133, 4.275, 4.428, 4.592, 4.769, 4.959};

  double ab_PbF2[] = {407,   403.3, 379.1, 406.3, 409.7, 408.9, 406.7, 404.7, 391.7, 397.7, 409.6, 403.7, 403.8, 409.7, 404.9, 404.2, 407.1, 411.1, 403.1,
                      406.1, 415.4, 399.1, 405.8, 408.2, 385.7, 405.6, 405.2, 401.6, 402.6, 407.1, 417.7, 401.1, 389.9, 411.9, 400.9, 398.3, 402.1, 408.7,
                      384.8, 415.8, 413.1, 385.7, 353.7, 319.1, 293.6, 261.9, 233.6, 204.4, 178.3, 147.6, 118.2, 78.7,  51.6,  41.5,  24.3,  8.8};
  double ref_PbF2[] = {1.749, 1.749, 1.75,  1.75,  1.751, 1.752, 1.752, 1.753, 1.754, 1.754, 1.755, 1.756, 1.757, 1.757, 1.758, 1.759, 1.76,  1.761, 1.762,
                       1.764, 1.765, 1.766, 1.768, 1.769, 1.771, 1.772, 1.774, 1.776, 1.778, 1.78,  1.782, 1.785, 1.787, 1.79,  1.793, 1.796, 1.8,   1.804,
                       1.808, 1.813, 1.818, 1.824, 1.83,  1.837, 1.845, 1.854, 1.865, 1.877, 1.892, 1.91,  1.937, 1.991, 1.38,  1.915, 1.971, 2.019};

  for (int i = 0; i < num; i++) {
    WaveLength[i] = (300 + i * 10) * nanometer;
    //    AirAbsorption[i] = 4*cm; // if photon in the air -> kill it immediately
    AirAbsorption[i] = 4 * cm; // if photon in the air -> kill it immediately
    // AirAbsorption[i] = 400000000*cm; //rd for air gap
    AirRefractiveIndex[i] = 1;
    PhotonEnergy[num - (i + 1)] = LambdaE / WaveLength[i];

    /* as the absorption is given per length and G4 needs
       mean free path length, calculate it here
       mean free path length - taken as probability equal 1/e
       that the photon will be absorbed */

    EpotekAbsorption[i] = (-1) / log(EpotekAbsorption[i]) * 0.001 * 2.54 * cm;
    QuartzAbsorption[i] = (-1) / log(QuartzAbsorption[i]) * 100 * cm;
    KamLandOilAbsorption[i] = (-1) / log(KamLandOilAbsorption[i]) * 50 * cm;
  }
  for (int i = 0; i < 25; i++) {
    Nlak33bEn[i] *= eV;
    Nlak33bAb[i] = (-0.8) / log(Nlak33bAb[i]) * 1 * cm; // account for glue in lens
    // std::cout << "en " << Nlak33bEn[i] <<" "<<  Nlak33bAb[i] << std::endl;    
  }

  /**************************** REFRACTIVE INDEXES ****************************/

  // only phase refractive indexes are necessary -> g4 calculates group itself !!

  double QuartzRefractiveIndex[num] = {1.456535, 1.456812, 1.4571,   1.457399, 1.457712, 1.458038, 1.458378, 1.458735, 1.459108, 1.4595,   1.459911, 1.460344,
                                       1.460799, 1.46128,  1.461789, 1.462326, 1.462897, 1.463502, 1.464146, 1.464833, 1.465566, 1.46635,  1.46719,  1.468094,
                                       1.469066, 1.470116, 1.471252, 1.472485, 1.473826, 1.475289, 1.476891, 1.478651, 1.480592, 1.482739, 1.485127, 1.487793};

  double EpotekRefractiveIndex[num] = {1.554034, 1.555575, 1.55698,  1.558266, 1.559454, 1.56056,  1.561604, 1.562604, 1.563579, 1.564547, 1.565526, 1.566536,
                                       1.567595, 1.568721, 1.569933, 1.57125,  1.57269,  1.574271, 1.576012, 1.577932, 1.580049, 1.582381, 1.584948, 1.587768,
                                       1.590859, 1.59424,  1.597929, 1.601946, 1.606307, 1.611033, 1.616141, 1.621651, 1.62758,  1.633947, 1.640771, 1.64807};

  double KamLandOilRefractiveIndex[num] = {1.433055, 1.433369, 1.433698, 1.434045, 1.434409, 1.434793, 1.435198, 1.435626, 1.436077, 1.436555, 1.4371,   1.4376,
                                           1.4382,   1.4388,   1.4395,   1.4402,   1.4409,   1.4415,   1.4425,   1.4434,   1.4444,   1.4455,   1.4464,   1.4479,
                                           1.4501,   1.450428, 1.451976, 1.453666, 1.455513, 1.45754,  1.45977,  1.462231, 1.464958, 1.467991, 1.471377, 1.475174};

  double Nlak33aRefractiveIndex[76] = {
    1.73816, 1.73836, 1.73858, 1.73881, 1.73904, 1.73928, 1.73952, 1.73976, 1.74001, 1.74026, 1.74052, 1.74078, 1.74105, 1.74132, 1.7416,  1.74189, 1.74218, 1.74249, 1.74279,
    1.74311, 1.74344, 1.74378, 1.74412, 1.74448, 1.74485, 1.74522, 1.74562, 1.74602, 1.74644, 1.74687, 1.74732, 1.74779, 1.74827, 1.74878, 1.7493,  1.74985, 1.75042, 1.75101,
    1.75163, 1.75228, 1.75296, 1.75368, 1.75443, 1.75521, 1.75604, 1.75692, 1.75784, 1.75882, 1.75985, 1.76095, 1.76211, 1.76335, 1.76467, 1.76608, 1.76758, 1.7692,  1.77093,
    1.77279, 1.7748,  1.77698, 1.77934, 1.7819,  1.7847,  1.78775, 1.79111, 1.79481, 1.79889, 1.80343, 1.8085,  1.81419, 1.82061, 1.8279,  1.83625, 1.84589, 1.85713, 1.87039};

  /* ASSIGNING REFRACTIVE AND ABSORPTION PROPERTIES TO THE GIVEN MATERIALS */

  // Quartz material => Si02
  G4MaterialPropertiesTable *QuartzMPT = new G4MaterialPropertiesTable();
  QuartzMPT->AddProperty("RINDEX", PhotonEnergy, QuartzRefractiveIndex, num);
  QuartzMPT->AddProperty("ABSLENGTH", PhotonEnergy, QuartzAbsorption, num);
  BarMaterial->SetMaterialPropertiesTable(QuartzMPT);

  // Air
  G4MaterialPropertiesTable *AirMPT = new G4MaterialPropertiesTable();
  AirMPT->AddProperty("RINDEX", PhotonEnergy, AirRefractiveIndex, num);
  AirMPT->AddProperty("ABSLENGTH", PhotonEnergy, AirAbsorption, num);
  //  assign this parameter table to the air
  defaultMaterial->SetMaterialPropertiesTable(AirMPT);

  // KamLandOil
  G4MaterialPropertiesTable *KamLandOilMPT = new G4MaterialPropertiesTable();
  KamLandOilMPT->AddProperty("RINDEX", PhotonEnergy, KamLandOilRefractiveIndex, num);
  KamLandOilMPT->AddProperty("ABSLENGTH", PhotonEnergy, KamLandOilAbsorption, num);
  OilMaterial->SetMaterialPropertiesTable(KamLandOilMPT);

  // N-Lak 33a
  for (int i = 0; i < 76; i++) {
    PhotonEnergyNlak33a[i] *= eV;
    Nlak33aAbsorption[i] *= cm; // cm to mm
  }
  G4MaterialPropertiesTable *Nlak33aMPT = new G4MaterialPropertiesTable();
  Nlak33aMPT->AddProperty("RINDEX", PhotonEnergyNlak33a, Nlak33aRefractiveIndex, 76);
  // Nlak33aMPT->AddProperty("ABSLENGTH",PhotonEnergyNlak33a, Nlak33aAbsorption, 76);
  Nlak33aMPT->AddProperty("ABSLENGTH", Nlak33bEn, Nlak33bAb, 25);
  Nlak33aMaterial->SetMaterialPropertiesTable(Nlak33aMPT);

  // PbF2
  G4MaterialPropertiesTable *PbF2MPT = new G4MaterialPropertiesTable();
  PbF2MPT->AddProperty("RINDEX", en_PbF2, ref_PbF2, n_PbF2);
  PbF2MPT->AddProperty("ABSLENGTH", en_PbF2, ab_PbF2, n_PbF2);
  PbF2Material->SetMaterialPropertiesTable(PbF2MPT);

  // Optical grease
  G4MaterialPropertiesTable *opticalGreaseMPT = new G4MaterialPropertiesTable();
  double og_en[6] = {1 * eV, 2 * eV, 3 * eV, 4 * eV, 4.2 * eV, 10 * eV};
  double og_ab[6] = {0.66 * cm, 0.66 * cm, 0.47 * cm, 0.14 * cm, 0.06 * cm, 0.02 * cm};
  double og_re[6] = {1.55, 1.56, 1.59, 1.64, 1.64, 1.64};
  opticalGreaseMPT->AddProperty("RINDEX", og_en, og_re, 6);
  opticalGreaseMPT->AddProperty("ABSLENGTH", og_en, og_ab, 6);

  opticalGreaseMaterial = new G4Material("opticalGreaseMaterial", density = 2.200 * g / cm3, ncomponents = 2);
  opticalGreaseMaterial->AddElement(Si, natoms = 1);
  opticalGreaseMaterial->AddElement(O, natoms = 2);
  opticalGreaseMaterial->SetMaterialPropertiesTable(opticalGreaseMPT);

  // Optical cookie (RTV615)
  G4MaterialPropertiesTable *opticalCookieMPT = new G4MaterialPropertiesTable();
  double oc_en[9] = {1.50 * eV, 2.00 * eV, 2.50 * eV, 3.00 * eV, 3.50 * eV, 4.00 * eV, 4.10 * eV, 4.50 * eV, 5.00 * eV};
  double oc_ab[9] = {14.2 * cm, 14.2 * cm, 14.2 * cm, 11.54 * cm, 5.29 * cm, 2.98 * cm, 2.43 * cm, 2.43 * cm, 2.43 * cm};
  double oc_re[9] = {1.406, 1.406, 1.406, 1.406, 1.406, 1.406, 1.406, 1.406, 1.406};
  opticalCookieMPT->AddProperty("RINDEX", oc_en, oc_re, 9);
  opticalCookieMPT->AddProperty("ABSLENGTH", oc_en, oc_ab, 9);

  opticalCookieMaterial = new G4Material("opticalCookieMaterial", density = 2.200 * g / cm3, ncomponents = 2);
  opticalCookieMaterial->AddElement(Si, natoms = 1);
  opticalCookieMaterial->AddElement(O, natoms = 2);
  opticalCookieMaterial->SetMaterialPropertiesTable(opticalCookieMPT);

  // Epotek Glue
  G4MaterialPropertiesTable *EpotekMPT = new G4MaterialPropertiesTable();
  EpotekMPT->AddProperty("RINDEX", PhotonEnergy, EpotekRefractiveIndex, num);
  EpotekMPT->AddProperty("ABSLENGTH", PhotonEnergy, EpotekAbsorption, num);
  // assign this parameter table to the epotek
  epotekMaterial->SetMaterialPropertiesTable(EpotekMPT);
}

void PrtDetectorConstruction::SetVisualization() {

  G4Colour blue = G4Colour(0.0, 0.0, 1.0);
  G4Colour green = G4Colour(0.0, 1.0, .0);
  G4Colour red = G4Colour(1.0, 0.0, .0);
  G4Colour DircColour = G4Colour(1., 1.0, 0.);
  G4Colour Dark = G4Colour(0., 0.05, 0.05, 0.15);
  // G4Colour Dark = G4Colour(0.,0.85,0.85,0.15);

  G4VisAttributes *waExpHall = new G4VisAttributes(DircColour);
  waExpHall->SetVisibility(false);
  lExpHall->SetVisAttributes(waExpHall);

  G4VisAttributes *waDirc = new G4VisAttributes(DircColour);
  waDirc->SetVisibility(false);
  lDirc->SetVisAttributes(waDirc);
  lCover->SetVisAttributes(waDirc);
  lFront->SetVisAttributes(waDirc);
  lEdd->SetVisAttributes(waDirc);

  G4VisAttributes *waBar = new G4VisAttributes(G4Colour(0., 1., 0.9, 0.2));
  // G4VisAttributes *waBar = new G4VisAttributes(Dark);
  waBar->SetVisibility(true);
  lBar->SetVisAttributes(waBar);

  G4VisAttributes *waMirror = new G4VisAttributes(G4Colour(1., 1., 0.9));
  waMirror->SetVisibility(true);
  lMirror->SetVisAttributes(waMirror);

  double transp = 0.15;
  if (fLensId != 0 && fLensId != 10) {
    G4VisAttributes *vaLens = new G4VisAttributes(G4Colour(0., 1., 1., transp));
    // vaLens->SetForceWireframe(true);
    // vaLens->SetForceAuxEdgeVisible(true);
    lLens1->SetVisAttributes(vaLens);
    G4VisAttributes *vaLens1 = new G4VisAttributes(G4Colour(0., 0.5, 1.0, transp));
    if (fLensId != 4) lLens2->SetVisAttributes(vaLens1);
    if (fLensId == 2) {
      lLens1->SetVisAttributes(vaLens1);
      lLens2->SetVisAttributes(vaLens);
    }
    if (fLensId == 3 || fLensId == 6 || fLensId == 7 || fLensId == 8) lLens3->SetVisAttributes(vaLens);
  }

  G4VisAttributes *waPrizm = new G4VisAttributes(G4Colour(0., 0.9, 0.9, 0.2));
  // G4VisAttributes *waPrizm = new G4VisAttributes(Dark);
  waPrizm->SetVisibility(true);
  // waPrizm->SetForceAuxEdgeVisible(true);
  // waPrizm->SetForceSolid(true);
  lPrizm->SetVisAttributes(waPrizm);

  G4VisAttributes *waMcp = new G4VisAttributes(green);
  waMcp->SetForceWireframe(true);
  lMcp->SetVisAttributes(waMcp);

  G4VisAttributes *waPixel = new G4VisAttributes(red);
  waPixel->SetForceWireframe(true);
  lPixel->SetVisAttributes(waPixel);
}

void PrtDetectorConstruction::ConstructSDandField() {
  // Sensitive detectors
  PrtPixelSD *pixelSD = new PrtPixelSD("PixelSD", "PixelHitsCollection", 0);
  G4SDManager::GetSDMpointer()->AddNewDetector(pixelSD);
  SetSensitiveDetector("lPixel", pixelSD);
  // SetSensitiveDetector("lScan",pixelSD);

  PrtPrizmSD *prizmSD = new PrtPrizmSD("PrizmSD", "PrizmHitsCollection", 0);
  G4SDManager::GetSDMpointer()->AddNewDetector(prizmSD);
  SetSensitiveDetector("lPrizm", prizmSD);

  // PrtTriggerSD* triggerSD = new PrtTriggerSD("TriggerSD", "TriggerHitsCollection", 0);
  // G4SDManager::GetSDMpointer()->AddNewDetector(triggerSD);
  // SetSensitiveDetector("lTrigger",triggerSD);

  PrtBarSD *barSD = new PrtBarSD("BarSD", "BarHitsCollection", 0);
  G4SDManager::GetSDMpointer()->AddNewDetector(barSD);
  SetSensitiveDetector("lBar", barSD);

  // Magnetic field
}

void PrtDetectorConstruction::SetRotation(double angle) {

  fPrtRot->rotateY(-fRotAngle);
  fPrtRot->rotateY(angle);
  fRotAngle = angle;

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void PrtDetectorConstruction::SetLens(int id) {

  // if(id==0){
  //   fPrismShift.setZ(fPrismShift.z()-fLens[2]);
  // }
  // std::cout<<"id  "<<id <<std::endl;

  // G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void PrtDetectorConstruction::SetQuantumEfficiency(int id) {
  const int num = 36;
  // ideal pmt quantum efficiency
  double QuantumEfficiencyIdial[num] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  // Burle PMT's
  double QuantumEfficiencyB[num] = {0.,   0.001, 0.002, 0.005, 0.01,  0.015, 0.02,  0.03, 0.04, 0.05, 0.06, 0.07, 0.09, 0.1,  0.13, 0.15, 0.17, 0.2,
                                    0.24, 0.26,  0.28,  0.282, 0.284, 0.286, 0.288, 0.29, 0.28, 0.26, 0.24, 0.22, 0.20, 0.18, 0.15, 0.13, 0.12, 0.10};

  // hamamatsu pmt quantum efficiency
  double QuantumEfficiencyPMT[num] = {0.001, 0.002, 0.004, 0.007, 0.011, 0.015, 0.020, 0.026, 0.033, 0.040, 0.045, 0.056, 0.067, 0.085, 0.109, 0.129, 0.138, 0.147,
                                      0.158, 0.170, 0.181, 0.188, 0.196, 0.203, 0.206, 0.212, 0.218, 0.219, 0.225, 0.230, 0.228, 0.222, 0.217, 0.210, 0.199, 0.177};

  if (id == 0) fQuantumEfficiency = QuantumEfficiencyIdial;
  if (id == 1) fQuantumEfficiency = QuantumEfficiencyPMT;

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}
