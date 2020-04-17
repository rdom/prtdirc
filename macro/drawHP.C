#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

void drawHP(TString infile="../build/hits.root"){
  
  if(!prt_init(infile,1)) return;

  for (int ievent=0; ievent< prt_entries; ievent++){
    prt_nextEvent(ievent,1000);
    for(auto hit : prt_event->GetHits()){
      int mcp = hit.GetMcpId();
      int pix = hit.GetPixelId()-1;
      double time = hit.GetLeadTime();      

      int ch = map_mpc[mcp][pix];      
      if(prt_isBadChannel(ch)) continue;
      if(mcp>7) continue;
      if(prt_pid==4 && time<30) prt_hdigi[mcp]->Fill(pix%8, pix/8);
    }
  }

  auto cdigi = prt_drawDigi(2018,0,0);
  prt_canvasAdd(cdigi);
  prt_canvasSave("data/drawHP",0);
}

