// -----------------------------------------
// PrtLutReco.h
//
// Created on: 13.07.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------
// Class for reconstruction in DIRC using look-up table method

#ifndef PrtLutReco_h
#define PrtLutReco_h 1

#include "PrtRun.h"
#include "PrtEvent.h"
#include "PrtHit.h"
#include "../../prttools/PrtTools.h"
#include "PrtManager.h"
#include "PrtLutNode.h"

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
#include "TGraph.h"
#include "CLHEP/Units/SystemOfUnits.h"

class PrtLutReco {

 public:
  // Standard constructors
  PrtLutReco(TString infile, TString lutfile,TString pdffile, int verbose = 0);

  // Destructor
  ~PrtLutReco();
  void Run(int start = 0, int end = 0);
  void drawTheoryLines();
  
 private:
  bool FindPeak(double &cangle, double &spr, double &cangle_pi, double &spr_pi, double a,
                int tofpid = 0);
  int FindPdg(double mom, double cangle);
  int GetEdge(int mcpid, int pixid);
  void SearchClusters();
  void FitRing(double &x0, double &y0, double &theta);
  double fillLnDiffPPi(double cangle, int tofPid, double mom);
  double fillLnDiffPPi2(double cangle, int tofPid, double mom);
  void ResetHists();
  int  getneighbours(int m, int p);
  void getclusters();
 
  
  PrtRun *frun;
  PrtTools ft;
  int fmaxch,fnpmt,fnpix;
  int fDetectorID;
  double fBboxNum, fPipehAngle, fDphi, fBarPhi;
  int fMethod;
  int fRadiator;
  int fLens;
  int fStudyId;
  bool fTimeImaging;
  
  TClonesArray *fLut;
  TClonesArray *fTrackInfoArray;

  TFile *fFile;
  TTree *fTree;
  TChain *fChain;

  PrtEvent *fEvent;
  PrtHit fHit;

  // Verbosity level
  int fVerbose;
  int nevents;
  int fgg_i = 0, fgg_ind = 0;
  TString fInputFile;
  TH1F *fHist;
  TH1F *fHistPi;
  TH1F *fHisti;
  TF1 *fFit;
  TF1 *fFunc[5];
  double fAngle[5];
  int fPk;
  int fCor_level;
  double fCor_angle[12] = {0}, fCor_time[12] = {0}, fCor_time_refl[2] = {0};
  double fCorrSpr;
  TString fCorrPath;
  TString fPdfPath;
  TGraph *fPdf2[3072], *fPdf4[3072];
  TH1F *fTime2[3072], *fTime4[3072];
  TH1F *hTof[5], *hTofc[5];
  TH1F *hNph[5], *hNph_ti[5];

  TH1F *hLnDiffGr[5], *hLnDiffTi[5];
  TH1F *fHist0, *fHist0d, *fHist0r, *fHist0s, *fHist0i;

  PrtLutNode *fLutNode[5000];
  TH1F *fHistMcp[15], *fHistTime[15], *fHistTimeRefl[2];
  TH1F *fHistCh[5000];

  // cluster search
  int mcpdata[15][65];
  int cluster[15][65];
  int lneighbours[65];
  int lsize = 0;

};

#endif
