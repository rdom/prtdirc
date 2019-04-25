#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

void drawHP(TString infile="../build/hits.root"){
  
  if(!prt_init(infile,1,"data/drawHP")) return;
  PrtHit hit;
  for (Int_t ievent=0; ievent< prt_entries; ievent++){
    prt_nextEvent(ievent,1000);
    for(Int_t h=0; h<prt_event->GetHitSize(); h++){
      hit = prt_event->GetHit(h);
      Int_t mcp = hit.GetMcpId();
      Int_t pix = hit.GetPixelId()-1;
      Double_t time = hit.GetLeadTime();
      
      Int_t ch = map_mpc[mcp][pix];
      if(prt_isBadChannel(ch)) continue;
      if(mcp>7) continue;
      if(prt_pid==4 && time<30)
	prt_hdigi[mcp]->Fill(pix%8, pix/8);
    }
  }
  
  //i%3*4+i/3

  prt_drawDigi("m,p,v\n",2018,0,0);
  //prt_cdigi->SetName(Form("hp_dataProtonS332_%d_%2.1f",(Int_t)prt_theta,prt_phi));
  prt_canvasAdd(prt_cdigi);
  prt_canvasSave(0);
}

