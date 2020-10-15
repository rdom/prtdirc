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
  PrtLutReco(TString infile, TString lutfile, int verbose=0);

  // Destructor
  ~PrtLutReco();
  void Run(int start=0, int end=0);
  void drawTheoryLines();
  
private:
  bool FindPeak(double& cangle, double& spr,double& cangle_pi, double& spr_pi,double a, int tofpid=0);
  int FindPdg(double mom, double cangle);
  int GetEdge(int mcpid, int pixid);
  void SearchClusters();
  void FitRing(double& x0, double& y0, double& theta);
  double fillLnDiffPPi(double cangle, int tofPid, double mom);
  double fillLnDiffPPi2(double cangle, int tofPid, double mom);
  void ResetHists();
  int fDetectorID;  
  double fBboxNum,fPipehAngle,fDphi,fBarPhi;
  TRandom fRand;
  int fMethod;
  int fRadiator;
  int fStudyId;

  TClonesArray *fLut;
  TClonesArray *fTrackInfoArray;

  TFile *fFile; 
  TTree *fTree;
  TChain *fChain;

  PrtEvent* fEvent;
  PrtHit fHit;
  
  // Verbosity level
  int fVerbose;
  int nevents;
  TString fInputFile;
  TH1F *fHist;
  TH1F *fHistPi;
  TH1F *fHisti;
  TF1 *fFit;
  TF1 *fFunc[5];
  TSpectrum *fSpect;
  double fAngle[5];
  double fTest;
  double fCorr[8];
  TString fCorrFile;
};

#endif
