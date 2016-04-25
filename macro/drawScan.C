#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

void drawScan(TString infile="hits.root"){
  // infile="/SAT/hera/had1/dervish/data/prt/study/151/beam*Sx*.root";
  // infile="/SAT/hera/had1/dirc/testbeam/2015/proc/151/beam*C.root";

  // infile="/data.local/data/jun15/beam_15177135523C.root";
  infile="$HOME/proc/152/beam_15183022858S.root";
  fSavePath = "scan3";
  PrtInit(infile,1); //digi
  
  Int_t itest(0);
  PrtHit fHit;
  for (Int_t ievent=0; ievent< fCh->GetEntries(); ievent++){
    PrtNextEvent(ievent,1000);
    for(Int_t h=0; h<prt_event->GetHitSize(); h++){
      fHit = prt_event->GetHit(h);
      Int_t mcpid = fHit.GetMcpId();
      Int_t pixid = fHit.GetPixelId()-1;
      
      Double_t time = fHit.GetLeadTime();
      if(prt_event->GetParticle() ==211) fhDigi[mcpid]->Fill(pixid%8, pixid/8);
    }
  }
  itest = fTest1+50;
  //  drawDigi("m,p,v\n",2,-2,-2);
  drawDigi("m,p,v\n",3,-2,-2);
  cDigi->cd();

  (new TPaletteAxis(0.90,0.1,0.94,0.90,((TH1 *)(fhDigi[0])->Clone())))->Draw();
   
  cDigi->SetName(Form("sc_%d_%d",fAngle,fMomentum/1000));
  canvasAdd(cDigi);  
  canvasSave(1,0);
  
}
