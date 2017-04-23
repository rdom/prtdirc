#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"
#include "TRandom.h"

void drawHP(TString infile="../build/hits.root"){

  fSavePath = "data/hit_pattern_plate16";
  CreateMap();
  PrtInit(infile,1);

  PrtHit fHit;
  Int_t nEvents=fCh->GetEntries();
  for (Int_t ievent=0; ievent< nEvents; ievent++){
    PrtNextEvent(ievent,1000);
    for(Int_t h=0; h<prt_event->GetHitSize(); h++){
      fHit = prt_event->GetHit(h);
      Int_t mcpid = fHit.GetMcpId();
      Int_t pixid = fHit.GetPixelId()-1;
      Double_t time = fHit.GetLeadTime();
      if(pixid<0) continue;
      Int_t ch = map_mpc[mcpid][pixid];
 
      if(prt_pid==4)
	fhDigi[mcpid]->Fill(pixid%8, pixid/8);
    }
  }

  drawDigi("m,p,v\n",2017,0,0);
  cDigi->cd();
  cDigi->SetName(Form("hp_sim_pi_%d_s217",(Int_t)fAngle));
  fhDigi[5]->GetZaxis()->SetLabelSize(0.06);
  (new TPaletteAxis(0.82,0.1,0.86,0.90,((TH1 *)(fhDigi[5])->Clone())))->Draw();

  canvasAdd(cDigi);
  canvasSave(0,0);
  
}
