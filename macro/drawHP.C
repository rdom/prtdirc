#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

void drawHP(TString infile="../build/hits.root"){
  
  if(!prt_init(infile,1,"data/drawHP")) return;
 
  PrtHit hit;
  for (Int_t ievent=0; ievent< prt_entries && ievent< 50000; ievent++){
    prt_nextEvent(ievent,1000);
    for(Int_t h=0; h<prt_event->GetHitSize(); h++){
      hit = prt_event->GetHit(h);
      Int_t mcpid = hit.GetMcpId();
      Int_t pixid = hit.GetPixelId()-1;
      Double_t time = hit.GetLeadTime();
      Int_t ch = map_mpc[mcpid][pixid];
      if(prt_isBadChannel(ch)) continue;
      
      // if(mcpid%3==0 && pixid<32) continue;
      // if(mcpid%3==2 && pixid>=32) continue; 
      
      if(prt_pid==2)
	prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
    }
  }


  //i%3*4+i/3
 
  prt_drawDigi("m,p,v\n",2017,0,0);
  prt_cdigi->SetName(Form("hp_dataS317_%d_%2.1f",(Int_t)prt_theta,prt_phi));
  prt_canvasAdd(prt_cdigi);
  //prt_cdigi_palette->Draw();
  prt_canvasSave(1,0);
}
