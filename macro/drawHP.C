#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

void drawHP(TString infile="../build/hits.root"){
  
  if(!prt_init(infile,1,"data/drawHP_2017")) return;
 
  PrtHit hit;
  for (Int_t ievent=0; ievent< prt_entries; ievent++){
    prt_nextEvent(ievent,1000);
    for(Int_t h=0; h<prt_event->GetHitSize(); h++){
      hit = prt_event->GetHit(h);
      Int_t mcpid = hit.GetMcpId();
      Int_t pixid = hit.GetPixelId()-1;
      Double_t time = hit.GetLeadTime();
      Int_t ch = map_mpc[mcpid][pixid];

      // if(mcpid%3==0 && pixid<32) continue;
      // if(mcpid%3==2 && pixid>=32) continue; 
      
      // if(prt_pid==4)
	prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
    }
  }
 
  prt_drawDigi("m,p,v\n",prt_geometry,0,0);
  prt_cdigi->SetName(Form("hp_sim_%d_%d",(Int_t)prt_theta,(Int_t)prt_test1));
  prt_canvasAdd(prt_cdigi);
  prt_cdigi_palette->Draw();
  prt_canvasSave(1,0);
}
