// -----------------------------------------
// PrtManager.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtManager_h
#define PrtManager_h

#include "globals.hh"

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "TClonesArray.h"

#include "PrtEvent.h"
#include "PrtHit.h"

class PrtManager
{
  static PrtManager* fInstance;
  TFile *fRootFile;
  TTree *fTree;
  PrtEvent *fEvent;
  PrtHit *fHit;
  TH1F *fHist;

public:
  PrtManager(G4String outfile, G4int runtype);
  ~PrtManager(){};
  static PrtManager* Instance(G4String outfile="hits.root", G4int runtype=0);
  void Save()             { if(fRootFile) fRootFile->Write(); }
  void Fill();
  void FillLut();
  void AddEvent(PrtEvent event);
  void AddHit(PrtHit hit, TVector3 localpos, TVector3 digipos, TVector3 position);
  PrtEvent* Event(){ return fEvent; }
  
  // Mutators
  void SetRunType(int val){ fRunType = val; }
  void SetStudy(int val){ fStudy = val; }
  void SetPhysList(int val){ fPhysList = val; }
  void SetGeometry(int val){ fGeometry = val; }
  void SetBeamDimension(double val){ fBeamDimension = val; }
  void SetRadiator(int val){ fRadiator = val; }
  void SetLens(int val){ fLens = val; }
  void SetMcpLayout(int val){ fMcpLayout = val; }
  void SetAngle(double val){ fAngle = val; }
  void SetPhi(double val){ fPhi = val; }
  void SetRadiatorL(double val){ fRadiatorL = val; }
  void SetRadiatorW(double val){ fRadiatorW = val; }
  void SetRadiatorH(double val){ fRadiatorH = val; }
  void SetParticle(int val){ fParticle = val; }
  void SetMomentum(TVector3 val){ fMomentum = val; if(fRunType==0) fEvent->SetMomentum(fMomentum);}
  void SetCurrentCherenkov(double val){ fCurrentCherenkov = val; }
  void SetShift(double val){ fShift = val; }
  void SetTest1(double val){ fTest1 = val; }
  void SetTest2(double val){ fTest2 = val; }
  void SetTest3(double val){ fTest3 = val; }
  void SetPrismStepX(double val){ fPrismStepX = val; }
  void SetPrismStepY(double val){ fPrismStepY = val; }
  void SetRStepX(double val){ fRStepX = val; }
  void SetRStepY(double val){ fRStepY = val; }
  void SetLStepX(double val){ fLStepX = val; }
  void SetLStepY(double val){ fLStepY = val; }
  void SetBeamX(double val){ fBeamX = val; }
  void SetBeamZ(double val){ fBeamZ = val; }
  void SetTime(double val){ fTime = val; }
  void SetTimeRes(double val){ fTimeRes = val; }
  void SetInfo(TString val){ fInfo = val; }
  void AddInfo(TString val){ fInfo += val + "\n"; }
  void SetMixPiP(bool val){fMixPiP = val;}
  
  // Accessors
  int GetRunType(){ return fRunType; }
  int GetStudy(){ return fStudy; }
  int GetPhysList(){ return fPhysList; }
  int GetGeometry(){ return fGeometry; }
  double GetBeamDinsion(){ return fBeamDimension; }
  int GetRadiator(){ return fRadiator; }
  int GetLens(){ return fLens; }
  int GetMcpLayout(){ return fMcpLayout; }
  double GetAngle(){ return fAngle; }
  double GetPhi(){ return fPhi; }
  double GetRadiatorL(){ return fRadiatorL; }
  double GetRadiatorW(){ return fRadiatorW; }
  double GetRadiatorH(){ return fRadiatorH; }
  int GetParticle(){ return fParticle; }
  TVector3 GetMomentum(){ return fMomentum; }
  double GetCurrentCherenkov(){ return fCurrentCherenkov; }
  double GetShift(){ return fShift; }
  double GetTest1(){ return fTest1; }
  double GetTest2(){ return fTest2; }
  double GetTest3(){ return fTest3; }
  double GetPrismStepX(){ return fPrismStepX; }
  double GetPrismStepY(){ return fPrismStepY; }
  double GetRStepX(){ return fRStepX; }
  double GetRStepY(){ return fRStepY; }
  double GetLStepX(){ return fLStepX; }
  double GetLStepY(){ return fLStepY; }
  double GetBeamX(){ return fBeamX; }
  double GetBeamZ(){ return fBeamZ; }
  double GetTime(){ return fTime; }
  double GetTimeRes(){ return fTimeRes; }
  TString GetOutName(){return fOutName; }
  TString GetInfo() { return fInfo; }
  bool GetMixPiP(){ return fMixPiP; }

private:
  int fRunType;
  int fStudy;
  int fPhysList;
  int fGeometry;
  int fRadiator;
  int fLens;
  int fMcpLayout;
  double fAngle;
  double fPhi;
  double fRadiatorL;
  double fRadiatorW;
  double fRadiatorH;
  int fParticle;
  double fBeamDimension;
  TVector3 fMomentum;
  TClonesArray *fLut;
  TClonesArray *fTrackInfoArray;
  double fCurrentCherenkov;
  double fShift;
  double fTest1;
  double fTest2;
  double fTest3;
  double fPrismStepX;
  double fPrismStepY;
  double fRStepX;
  double fRStepY;
  double fLStepX;
  double fLStepY;
  double fBeamX;
  double fBeamZ;
  double fTime;
  double fTimeRes;
  TString fOutName;
  TString fInfo;
  bool fMixPiP;

  TVector3 fnX1;
  TVector3 fnY1;
  Double_t fCriticalAngle;
};

#endif
