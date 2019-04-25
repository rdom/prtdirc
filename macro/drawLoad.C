#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

const Double_t prismangle = 33;
Double_t prismSize[] = {50+300*tan(prismangle*TMath::Pi()/180.),175};
Double_t prismShift = (prismSize[0])/2. -50/2.;
void drawPrism(Double_t x, Double_t y){
  gPad->cd();
  TBox *box2 = new TBox(x-prismSize[0]/2.,y-prismSize[1]/2.,x+prismSize[0]/2.,y+prismSize[1]/2.);
  box2->SetFillStyle(0);
  box2->SetLineColor(4);
  box2->SetLineWidth(2);
  box2->Draw("same");
}

void drawLoad(TString infile="../build/hits.root"){
  if(!prt_init(infile,1,"data/load")) return;
  
  TH2F* hHits = new TH2F("hHits",";x [mm];y [mm]",500,-40,230,500,-100,100);
  TH2F* hTime = new TH2F("hTime",";x [mm];time [ns]",500,-40,230,250,0,25);
  int angle(0), step(0);

  for (int ievent=0; ievent<prt_entries; ievent++){
    prt_nextEvent(ievent,1000);
    for(auto hit : prt_event->GetHits()){
      int mcp = hit.GetMcpId();
      int pix = hit.GetPixelId()-1;
      TVector3 pos = hit.GetPosition();      
      Double_t time = hit.GetLeadTime();

      if(fabs(pos.Y()-20)>3) continue;
      hHits->Fill(pos.X(),pos.Y());
      hTime->Fill(pos.X(),time);      
    }
  }

  gStyle->SetOptStat(0);

  prt_canvasAdd(Form("time_%d",angle),800,500);
  hTime->Draw("colz");
  
  prt_canvasAdd(Form("load_%d",angle),800,500);
  hHits->Draw("colz");
  
  drawPrism(prismShift,0);
  prt_canvasSave(1,0);
}
