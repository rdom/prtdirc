#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"
#include "TRandom.h"

void drawScan(TString infile="../build/hits.root"){

  fSavePath = "scan3";
  CreateMap();
  PrtInit(infile,1); //digi
  TGaxis::SetMaxDigits(3);

  TH1F *hTime1 = new TH1F("hTimeP",";LE time [ns]; entries [#]", 500,0,50);
  TH1F *hTime2 = new TH1F("hTimePi",";LE time [ns]; entries [#]", 500,0,50);

  TH1F *hTime1p = new TH1F("hTimepP",";LE time [ns]; entries [#]",300,0,30);
  TH1F *hTime2p = new TH1F("hTimepPi",";LE time [ns]; entries [#]", 300,0,30);


  TRandom rand;
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
      if(badcannel(ch)) continue;
      if(fHit.GetChannel()>960 || ch==0) continue;
 
      //if(prt_event->GetParticle() ==211)
	fhDigi[mcpid]->Fill(pixid%8, pixid/8);

      Double_t time = fHit.GetLeadTime();//+ rand.Gaus(0,0.4);
      if(time<10)continue;
      if(prt_event->GetParticle() ==2212){
	hTime1->Fill(time);
	if(ch==291) hTime1p->Fill(time);
      }else{
	hTime2->Fill(time);
	if(ch==291) hTime2p->Fill(time);
      }
    }
  }
    
  std::cout<<"total  "<<total <<std::endl;
  
  itest = fTest1+50;
  //  drawDigi("m,p,v\n",2,-2,-2);
  // drawDigi("m,p,v\n",3,-2,-2);
  drawDigi("m,p,v\n",3,0,0);
  cDigi->cd();
  fhDigi[0]->GetZaxis()->SetLabelSize(0.06);
  (new TPaletteAxis(0.90,0.1,0.94,0.90,((TH1 *)(fhDigi[0])->Clone())))->Draw();

  TString name = Form("152_7GeV_211_data_%d_%d",fAngle,fMomentum/1000);
  cDigi->SetName("sc_"+name);
  canvasAdd(cDigi);

  canvasAdd("cTime"+name);
  hTime1->Scale(hTime2->GetMaximum()/hTime1->GetMaximum());
  
  hTime1->SetLineColor(2);
  hTime1->Draw();
  hTime2->SetLineColor(4);
  hTime2->Draw("same");




  gStyle->SetOptStat(0);
  canvasAdd("time_152_ch291"+name,800,500);
  hTime1p->Scale(hTime2p->GetMaximum()/hTime1p->GetMaximum());
  
  hTime1p->SetLineColor(2);
  hTime1p->Draw();
  hTime2p->SetLineColor(4);
  hTime2p->Draw("same");

  TLegend *leg = new TLegend(0.65,0.7,0.8,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hTime1p,"protons ","lp");
  leg->AddEntry(hTime2p,"pions","lp");
  leg->Draw();
  
  canvasSave(0,0);
  
}
