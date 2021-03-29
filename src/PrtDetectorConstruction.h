// -----------------------------------------
// PrtDetectorConstruction class
//
// author  : r.dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtDetectorConstruction_h
#define PrtDetectorConstruction_h 1

#include "globals.hh"
#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"

#include "PrtRun.h"
#include "PrtDetectorConstructionMessenger.h"

class PrtDetectorConstructionMessenger;

class PrtDetectorConstruction : public G4VUserDetectorConstruction {
 public:
  PrtDetectorConstruction();
  virtual ~PrtDetectorConstruction();

 public:
  virtual G4VPhysicalVolume *Construct();
  virtual void ConstructSDandField();
  void DefineMaterials();
  void SetVisualization();
  void SetRotation(double angle);
  void SetLens(int id);
  void SetQuantumEfficiency(int id);

 private:
  PrtRun *fRun;
  G4LogicalVolume *lExpHall;
  G4LogicalVolume *lFront;
  G4LogicalVolume *lEdd;
  G4LogicalVolume *lTrigger;
  G4LogicalVolume *lDirc;
  G4LogicalVolume *lBar;
  G4LogicalVolume *lOpticalGrease;
  G4LogicalVolume *lOpticalGreased;
  G4LogicalVolume *lCookie1, *lCookie2;
  G4LogicalVolume *lMirrorGap;
  G4LogicalVolume *lMirror;
  G4LogicalVolume *lLens1;
  G4LogicalVolume *lLens2;
  G4LogicalVolume *lLens3;
  G4LogicalVolume *lPrizm;
  G4LogicalVolume *lMcp;
  G4LogicalVolume *lPixel;
  G4LogicalVolume *lPrizmC;
  G4LogicalVolume *lCover;

  G4VPhysicalVolume *wBar;
  G4VPhysicalVolume *wMirrorGap;
  G4VPhysicalVolume *wMirror;
  G4VPhysicalVolume *wDirc;

  G4Material *defaultMaterial; // material for bars
  G4Material *BarMaterial;     // material for bars
  G4Material *OilMaterial;
  G4Material *MirrorMaterial;        // material of mirror
  G4Material *opticalGreaseMaterial; // material of mirror
  G4Material *opticalCookieMaterial;
  G4Material *epotekMaterial;
  G4Material *Nlak33aMaterial;
  G4Material *PbF2Material;
  G4Material *frontMaterial;

  int fNRow;
  int fNCol;
  int fRunType;
  int fGeomId;
  int fLensId;
  int fRadiatorId;
  double fPrismStepX;
  int fMcpLayout;
  int fTest1, fTest2, fTest3;
  double fHall[3];
  double fBar[3];
  double fMirror[3];
  double fPrizm[4];
  double fLens[4];
  double fMcpTotal[3];
  double fMcpActive[3];
  G4ThreeVector fPrismShift;
  G4ThreeVector fCenterShift;
  double fOffset;

  double fRotAngle;
  G4RotationMatrix *fPrtRot;
  PrtDetectorConstructionMessenger *fGeomMessenger;
  double *fQuantumEfficiency;
};

#endif
