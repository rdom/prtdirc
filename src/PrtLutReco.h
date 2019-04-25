// -----------------------------------------
// PrtLutReco.h
//
// Created on: 13.07.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------
// Class for reconstruction in DIRC using look-up table method
 
#ifndef PrtLutReco_h
#define PrtLutReco_h 1

#include "PrtEvent.h"
#include "PrtHit.h"

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TRandom.h"
#include "TVirtualFitter.h"
#include "TArc.h"
#include "CLHEP/Units/SystemOfUnits.h"

class PrtLutReco{

public:

  // Standard constructors
  PrtLutReco(TString infile, TString lutfile, Int_t verbose=0);

  // Destructor
  ~PrtLutReco();
  void Run(Int_t start=0, Int_t end=0);
  void drawTheoryLines();
  
private:
  Bool_t FindPeak(Double_t& cherenkovreco, Double_t& spr,Double_t a, Int_t tofpid=0);
  Int_t FindPdg(Double_t mom, Double_t cangle);
  Int_t GetEdge(Int_t mcpid, Int_t pixid);
  void SearchClusters();
  void FitRing(Double_t& x0, Double_t& y0, Double_t& theta);
  Double_t fillLnDiffPPi(Double_t cangle, Int_t tofPid, Double_t mom);
  Double_t fillLnDiffPPi2(Double_t cangle, Int_t tofPid, Double_t mom);
  void ResetHists();
  Int_t fDetectorID;  
  Double_t fBboxNum,fPipehAngle,fDphi,fBarPhi;
  TRandom fRand;
  Int_t fMethod;
  Int_t fRadiator;

  TClonesArray *fLut;
  TClonesArray *fTrackInfoArray;

  TFile *fFile; 
  TTree *fTree;
  TChain *fChain;

  PrtEvent* fEvent;
  PrtHit fHit;
  
  // Verbosity level
  Int_t fVerbose;
  Int_t nevents;
  TString fInputFile;
  TH1F *fHist;
  TH1F *fHistPi;
  TH1F *fHisti;
  TF1 *fFit;
  TSpectrum *fSpect;
  Double_t fAngleP;
  Double_t fAnglePi;
  Double_t fTest;
};

#endif
