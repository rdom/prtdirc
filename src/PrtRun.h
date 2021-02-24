// -----------------------------------------
// PrtRun.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtRun_h
#define PrtRun_h 1

#include <iostream>
#include "TObject.h"
#include "TVector3.h"

class PrtRun: public TObject  {

public:

  PrtRun(); 	//the default constructor
  ~PrtRun(){};
  void PrintInfo();

  // Accessors 
  Int_t runId()         const { return frunId; }
  Int_t studyId()       const { return fstudyId; }
  Bool_t mc()           const { return fmc; }
  Double_t theta()      const { return ftheta; }
  Double_t phi()        const { return fphi; }
  Int_t physList()      const { return fphysList; }
  Int_t pid()           const { return fpid; }
  TVector3 momentum()   const { return fmomentum; }
  TVector3 position()   const { return fposition; }
  Int_t geometry()      const { return fgeometry; }
  Int_t lens()          const { return flens; }
  Int_t trigger()       const { return ftrigger; } 
  Double_t prismStepX() const { return fprismStepX; }
  Double_t prismStepY() const { return fprismStepY; }
  Double_t beamX()      const { return fbeamX; }
  Double_t beamZ()      const { return fbeamZ; }
  Double_t timeSigma()  const { return ftimeSigma; }
  Double_t tof()        const { return ftof; }
  Double_t tofPi()      const { return ftofPi; }
  Double_t tofP()       const { return ftofP; }  
  Double_t test1()      const { return ftest1; }
  Double_t test2()      const { return ftest2; }
  
  // Mutators
  void runId(Int_t v)      { frunId = v; }
  void studyId(Int_t v)    { fstudyId = v; }
  void mc(Bool_t v)    { fmc = v; }
  void theta(Double_t v)   { ftheta = v; }
  void phi(Double_t v)     { fphi = v; }
  void physList(Int_t v) { fphysList = v; }
  void pid(Int_t v) { fpid = v; }
  void momentum(TVector3 v) { fmomentum = v; }
  void position(TVector3 v) { fposition = v; }
  void geometry(Int_t v) { fgeometry = v; }
  void lens(Int_t v) { flens = v; }
  void trigger(Int_t v) { ftrigger = v; }
  void prismStepX(Double_t v){ fprismStepX = v; }
  void prismStepY(Double_t v){ fprismStepY = v; }
  void beamX(Double_t v){ fbeamX = v; }
  void beamZ(Double_t v){ fbeamZ = v; }
  void timeSigma(Double_t v){ ftimeSigma = v; }
  void tof(Double_t v){ ftof = v; }
  void tofPi(Double_t v){ ftofPi = v; }
  void tofP(Double_t v){ ftofP = v; }  
  void test1(Double_t v) { ftest1 = v; }
  void test2(Double_t v) { ftest2 = v; }

private: 
  Int_t frunId;
  Int_t fstudyId;
  Bool_t fmc;
  Int_t fphysList;
  Int_t fpid;
  Int_t fgeometry;
  Int_t flens;
  Int_t fradiator;
  Int_t ftrigger;
  
  Double_t ftheta;
  Double_t fphi;
  TVector3 fmomentum;
  TVector3 fposition;

  Double_t fprismStepX;
  Double_t fprismStepY;
  Double_t fbeamX;
  Double_t fbeamZ;
  Double_t fbeamSize;
  Double_t ftimeSigma;

  Double_t ftof;
  Double_t ftofPi;
  Double_t ftofP;

  Double_t ftest1;
  Double_t ftest2;
  
  ClassDef(PrtRun, 0);
};
#endif
