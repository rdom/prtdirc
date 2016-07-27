#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"
#include "TRandom.h"

void drawHP(TString infile="../build/hits.root"){

  fSavePath = "hit_pattern";
  CreateMap();
  PrtInit(infile,1); //digi
  TGaxis::SetMaxDigits(3);

  Int_t total(0), itest(0);
  PrtHit fHit;
  Int_t nEvents=fCh->GetEntries(); //60k
  for (Int_t ievent=0; ievent< nEvents; ievent++){
    PrtNextEvent(ievent,1000);
    if(prt_event->GetParticle() ==2212) total++;
    for(Int_t h=0; h<prt_event->GetHitSize(); h++){
      fHit = prt_event->GetHit(h);
      Int_t mcpid = fHit.GetMcpId();
      Int_t pixid = fHit.GetPixelId()-1;
      Int_t ch = map_mpc[mcpid][pixid];

      fhDigi[mcpid]->Fill(pixid%8, pixid/8);
    }
  }
    
  drawDigi("m,p,v\n",7,0,0);
  cDigi->cd();
  fhDigi[0]->GetZaxis()->SetLabelSize(0.06);
  (new TPaletteAxis(0.82,0.1,0.86,0.90,((TH1 *)(fhDigi[0])->Clone())))->Draw();

  canvasAdd(cDigi);
  
  canvasSave(1,0);
  
}
