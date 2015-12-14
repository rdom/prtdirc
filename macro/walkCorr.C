#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TMath.h"
#include "TChain.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TString.h"
#include "TArrayD.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TPaveStats.h"

#include "../src/PrtEvent.h"
#include "../src/PrtHit.h"

#define prt__beam
#include "prttools.C"

const int nmcp = 15, npix = 64, pmax = 100000;
TH1F * hPTime[nmcp][npix];
TH2F * hWcorr[nmcp][npix];

TH1F * hPTime0[nmcp][npix];
TH1F * hPTime1[nmcp][npix];
TH1F * hPTime2[nmcp][npix];

TH2F * hWcorr0[nmcp][npix];
TH2F * hWcorr1[nmcp][npix];
TH2F * hWcorr2[nmcp][npix];

// // pions
// Int_t trigger= 1920,dmcpid=0;
// Double_t xbin1= -86, xbin2= -72,xbin3= -250,xbin4= -200,xbin5= -165,xbin6= -145,ybin1= -1, ybin2= 4, ybin3= 103.1, ybin4= 108.1;

// pilas
Int_t trigger = 1952,dmcpid=0;
Double_t xbin1 = 86, xbin2 = 102, xbin3 = 20, xbin4 = 60, xbin5 = -70, xbin6 = 40,
  ybin1 = -4, ybin2 = 8,  ybin3 = 10, ybin4 = 10.3;

Double_t fbin1 = xbin1,fbin2=xbin2;
Double_t fbin3 = -1.2, fbin4=1.2;

TH1F *thist=new TH1F("Corrected","Corrected",   300,fbin1,fbin2);
TH2F *thist2=new TH2F("thist2","Uncorrected",400,fbin1,fbin2,200,fbin3,fbin4);
TH2F *thist3=new TH2F("thist3","Corrected",400,fbin1,fbin2,200,fbin3,fbin4);
TH2F *thist4=new TH2F("thist4","thist4",400,fbin1,fbin2,200,-0.5,0.5);

bool no2d = true;
Int_t cshift = 0;
TCanvas *cTime = new TCanvas("cTime","cTime",0,0,800,400);
TCanvas *cLeadingEdge = new TCanvas("cLeadingEdge","cLeadingEdge",0,0,800,400);  
double tpar, par[4] = {-0.5,0,-0.5,0};

Int_t uniqmcp = -8;
Float_t dataLeg[15][npix][pmax];
Float_t dataTo2[15][npix][pmax];
Float_t dataTo1[15][npix][pmax];

Double_t gfitted0 = 0, gfitted1 = 0;
Int_t gm, gp, idata[nmcp][npix];
TVector3 offsetLeg, offsetTot1, offsetTot2;
TSpectrum *spect;

Double_t dist(TVector3 v1, TVector3 v2){
  return sqrt((v1.X()-v2.X())*(v1.X()-v2.X())+(v1.Y()-v2.Y())*(v1.Y()-v2.Y())+(v1.Z()-v2.Z())*(v1.Z()-v2.Z()));
}

TVector3 fit(TH1F *h){
  int binmax = h->GetMaximumBin();
  double xmax = h->GetXaxis()->GetBinCenter(binmax);
  TF1 *gaust = new TF1("gaust","gaus(0)+gaus(3)",xmax-3,xmax+3);//+gaus(3)
  Double_t integral = h->Integral(h->GetXaxis()->FindBin(xmax-0.6),h->GetXaxis()->FindBin(xmax+0.6));
  Double_t xxmin, xxmax, sigma1, mean1, sigma2, mean2;
  
  if(integral>50){ 
    std::cout<<"integral  "<< integral<<std::endl;
    Int_t nfound = spect->Search(h,2,"",0.2);
    Float_t *xpeaks = spect->GetPositionX();
    if(nfound==1){
      xxmax = xmax;
      xxmin = xxmax;
      gaust =new TF1("gaust","gaus(0)",xmax-3,xmax+3);
      gaust->SetParameter(1,xpeaks[0]);
    }
    if(nfound==2){
      if(xpeaks[0]>xpeaks[1]) {
	xxmax = xpeaks[0];
	xxmin = xpeaks[1];
      }else {
	xxmax = xpeaks[1];
	xxmin = xpeaks[0];
      }
      gaust =new TF1("gaust","gaus(0)+gaus(3)",xmax-3,xmax+3);
      gaust->SetParameter(1,xxmin);
      gaust->SetParameter(4,xxmax);
    }
    h->GetXaxis()->SetRangeUser(xxmin-3, xxmax+3);

    gaust->SetParameter(2,0.3);
    gaust->SetParameter(5,0.3);

    h->Fit("gaust","","MQ",xxmin-0.7, xxmax+0.4);
    mean1 = gaust->GetParameter(1);
    sigma1 = gaust->GetParameter(2);
    //shiftHist(h,xmax);
    h->SetTitle(Form("%s , Offset %f ns",h->GetName(), xmax));
    mean2 = (nfound==1) ? gaust->GetParameter(1) : gaust->GetParameter(4);
    sigma2 = (nfound==1) ? gaust->GetParameter(2) : gaust->GetParameter(5);
  }
  return TVector3(mean1,mean2,sigma1);
}

TVector3 findPeak(TH2F* hh){
  TH1F *h = (TH1F*) hh->ProjectionY();
  TVector3 res = fit(h);
  h->Draw("same");
  return res;
}
TPaveStats *s;
void drawFitted(Double_t a, Double_t b){
  cTime->cd(4-cshift);
  gfitted1 = fit(thist).Z();
  //s->SetTextSize(0.5); 

  thist->Draw();
  thist->GetXaxis()->UnZoom();
  gPad->Update();
  s =  (TPaveStats*) thist->FindObject( "stats" );
  s->SetY1NDC(0.45);
  s->SetX1NDC(0.55);
  
  //std::cout<<"par0 = "<<a<<"   par1 = "<<b<<"   mmin = "<<mmin <<std::endl;
  
  if(!no2d){
    cTime->cd(1);
    thist2->Draw("colz");
    thist2->GetXaxis()->SetRangeUser(b-2,b+2);

    Double_t x[2], y[2];
    y[0] = -1; y[1] = 1;
    x[0] = y[0]*a+b; x[1] = y[1]*a+b;
    TGraph *gr = new TGraph(2,x,y);
    gr->SetLineColor(2);
    gr->Draw("same L");
    
    cTime->cd(2);
    thist3->Draw("colz");
    thist3->GetXaxis()->SetRangeUser(b-2,b+2);
    //cLeadingEdge->cd(2);
    //thist4->Draw("colz");
   
    cTime->Update();
    // wait(cTime);
  }
}

double getSigma(const double *xx ){
  thist->Reset();
  thist2->Reset();
  thist3->Reset();
  thist4->Reset();
  gm = (uniqmcp<0) ? gm : 0; 
  
  for(Int_t i=0; i<idata[gm][gp]; i++){
    Double_t leg = dataLeg[gm][gp][i];//-offsetLeg.X();
    Double_t ctot = dataTo1[gm][gp][i];
    //combine tot peaks
    if(abs(ctot-offsetTot1.X()) > abs(ctot-offsetTot1.Y())){
      ctot = ctot - offsetTot1.Y() +  offsetTot1.X();
      leg = (offsetTot1.X()>offsetTot1.Y()) ? leg-0.18: leg+0.18;
    }
    Double_t tot = ctot-offsetTot1.X();
    if(abs(leg) > 300 || abs(tot) > 1 ) continue; 

    Double_t corrleg = leg - tot * xx[0];
    thist->Fill(corrleg );
    thist2->Fill(leg,     tot);
    thist3->Fill(corrleg, tot);
    thist4->Fill(corrleg, dataTo2[gm][gp][i] - offsetTot2.X());
  }
  Double_t mmin = 0, xmax = thist->GetXaxis()->GetBinCenter(thist->GetMaximumBin());
  Double_t integral = thist->Integral(thist->GetXaxis()->FindBin(xmax-0.6),thist->GetXaxis()->FindBin(xmax+0.6));
  if(integral>50){
    thist->Fit("gaus","","MQ",xmax-0.5, xmax+0.5);
    TF1 *fit1 = (TF1*)thist->GetFunction("gaus");
    mmin = fit1->GetParameter(2);//abs(offsetLeg.X()-fit1->GetParameter(1));
    //mmin = fit(thist).Z();
  }
  //drawFitted(xx[0],xx[1]);
  return mmin;
}

void minimize(){
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");
  // set tolerance , etc...
  min->SetMaxFunctionCalls(10000); // for Minuit/Minuit2 
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);
  ROOT::Math::Functor f(&getSigma,4); 
  double step[4] = {0.1,0.1,0.1,0.1};
  min->SetFunction(f);
   
  min->SetVariable(0,"x1",par[0], step[0]);
  min->SetVariable(1,"y1",offsetLeg.X(), step[1]);
  min->SetVariable(2,"x2",par[2], step[2]);
  min->SetVariable(3,"y2",par[3], step[3]);

  // min->SetVariableLimits(1,offsetLeg.X()-1,offsetLeg.X()+1);
  min->FixVariable(1);
  min->FixVariable(2);
  min->FixVariable(3);
  min->SetPrintLevel(-1);

  min->Minimize(); 
  const double *xs = min->X();
  std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << min->MinValue()  << std::endl;
  drawFitted(xs[0],xs[1]);
}

void docorrections(Int_t mcpid, Int_t pixid){
  if(mcpid!=uniqmcp && uniqmcp > 0) return;
  cLeadingEdge->cd();
  hWcorr[mcpid][pixid]->GetXaxis()->SetTitle("LE time, [ns]");
  hWcorr[mcpid][pixid]->GetXaxis()->SetTitleSize(0.05);
  hWcorr[mcpid][pixid]->GetXaxis()->SetTitleOffset(0.75);
  hWcorr[mcpid][pixid]->GetYaxis()->SetTitle("signal TOT, [ns]");
  hWcorr[mcpid][pixid]->GetYaxis()->SetTitleSize(0.05);
  hWcorr[mcpid][pixid]->GetYaxis()->SetTitleOffset(0.65);
  hWcorr[mcpid][pixid]->Draw("colz");
  // cLeadingEdge->cd(2);
  // hWcorr0[mcpid][pixid]->Draw("colz");
  // cLeadingEdge->cd(3);
  // hWcorr1[mcpid][pixid]->Draw("colz");
  // cLeadingEdge->cd(4);
  // hWcorr2[mcpid][pixid]->Draw("colz");
  cLeadingEdge->Update();

  cTime->cd(3-cshift);
  offsetLeg = fit(hPTime[mcpid][pixid]);
  gfitted0 = offsetLeg.Z();
  gm = mcpid;
  gp = pixid;
  hPTime[mcpid][pixid]->Draw();
  hPTime[mcpid][pixid]->GetXaxis()->UnZoom();
  cLeadingEdge->cd(4);
  offsetTot1 = findPeak(hWcorr1[mcpid][pixid]);
  offsetTot2 = findPeak(hWcorr2[mcpid][pixid]);
  minimize();
  cTime->Update();
}

void walkCorr(TString inFile="hits_c.root", TString dirid="", TString fileid=""){
  //inFile="/SAT/hera/had1/dervish/data/dirc/pilas/pilas_cc14195184300fast.root";
  //inFile="/SAT/hera/had1/dervish/data/dirc/pilas/pilas_cc*fast.root";
  //  inFile="../../../data/cc14191172313High.root";
  //inFile="/SAT/hera/had1/dervish/data/dirc/scanangle1/1160*fast.root";
  //inFile="/SAT/hera/had1/dervish/data/dirc/scanangle2/324*fast.root";
  inFile="/SAT/hera/had1/dervish/data/dirc/pilas/pilas_cc14195184051fast.root";

  fSaveFlag = 2;
  PrtInit(inFile,1);


  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  // gStyle->SetStatW(0.35);
  // gStyle->SetStatH(0.4);
  // gROOT->ForceStyle();
  thist->SetStats(1);
  

 
  thist2->SetStats(0);
  thist2->GetXaxis()->SetTitle("LE time, [ns]");
  thist2->GetXaxis()->SetTitleSize(0.05);
  thist2->GetXaxis()->SetTitleOffset(0.95);
  thist2->GetYaxis()->SetTitle("signal TOT - mean signal TOT, [ns]");
  thist2->GetYaxis()->SetTitleSize(0.05);
  thist2->GetYaxis()->SetTitleOffset(0.95);

  thist3->SetStats(0);
  thist3->GetXaxis()->SetTitle("LE time (corrected), [ns]");
  thist3->GetXaxis()->SetTitleSize(0.05);
  thist3->GetXaxis()->SetTitleOffset(0.95);
  thist3->GetYaxis()->SetTitle("signal TOT - mean signal TOT, [ns]");
  thist3->GetYaxis()->SetTitleSize(0.05);
  thist3->GetYaxis()->SetTitleOffset(0.95);
  axisTime800x500(thist);
  thist->GetXaxis()->SetTitle("LE time (corrected), [ns]");


  const Int_t NumOfChans=2816;
  cLeadingEdge->Divide(1,4);
  if(no2d) {
    cshift =2;
    cTime->Divide(2,1);
  }else{
    cTime->Divide(2,2);
  }

  Double_t        fLeadingEdge[NumOfChans];
  Double_t        fTot[NumOfChans];
  UInt_t          nTdcId[NumOfChans];
  UInt_t          nTdcChan[NumOfChans];
  UInt_t          nMcpId[NumOfChans];
  UInt_t          nPixel[NumOfChans];
  UInt_t          nRow[NumOfChans];
  UInt_t          nCol[NumOfChans];
  Bool_t          bValid[NumOfChans];

  TChain * chain = new TChain("K");
  chain->Add(inFile);
  chain->SetBranchAddress("fLeadingEdge", fLeadingEdge);
  chain->SetBranchAddress("fTot",         fTot);
  chain->SetBranchAddress("nTdcId",       nTdcId);
  chain->SetBranchAddress("nTdcChan",     nTdcChan);
  chain->SetBranchAddress("nMcpId",       nMcpId);
  chain->SetBranchAddress("nPixel",       nPixel);
  chain->SetBranchAddress("nRow",         nRow);
  chain->SetBranchAddress("nCol",         nCol);
  chain->SetBranchAddress("bValid",       bValid);

  for(Int_t m=0; m<nmcp; m++){
    for(Int_t p=0; p<npix; p++){
      hPTime[m][p]   = new TH1F(Form("mcp%d,pixel%d",m,p),Form("mcp %d, pixel %d",m, p),  300,xbin1,xbin2);
      hPTime0[m][p]  = new TH1F(Form("mcp0%d,pixel%d",m,p),Form("mcp %d, pixel %d",m, p), 300,xbin1,xbin2);
      hPTime1[m][p]  = new TH1F(Form("mcp1m%dp%d",m,p),Form("mcp %d, pixel %d",m, p),     300,xbin3,xbin4);
      hPTime2[m][p]  = new TH1F(Form("mcp2m%dp%d",m,p),Form("mcp %d, pixel %d",m, p),     300,xbin5,xbin6);
   
      hWcorr[m][p]  = new TH2F(Form("hWcorrm%dp%d" ,m,p), Form("mcp %d, pixel %d",m, p), 200,xbin1,xbin2,200,ybin1,ybin2); 
      hWcorr0[m][p] = new TH2F(Form("hWcorr0m%dp%d",m,p), Form("mcp %d, pixel %d",m, p), 200,xbin1,xbin2,200,ybin3,ybin4); 
      hWcorr1[m][p] = new TH2F(Form("hWcorr1m%dp%d",m,p), Form("mcp %d, pixel %d",m, p), 200,xbin3,xbin4,200,ybin1,ybin2); 
      hWcorr2[m][p] = new TH2F(Form("hWcorr2m%dp%d",m,p), Form("mcp %d, pixel %d",m, p), 200,xbin5,xbin6,200,ybin3,ybin4); 
      axisTime800x500(hPTime[m][p]);
      axisTime800x500(hPTime0[m][p]);
      hPTime[m][p]->SetStats(1);
      hPTime0[m][p]->SetStats(1);
      idata[m][p]=0;
    }
  }

  Int_t entries = chain->GetEntries();
  std::cout<<"# of Entries  "<< entries <<std::endl;
  
  for (Int_t itrig=0; itrig<entries; itrig++) {
    if(itrig%10000==0) std::cout<<"Event # "<< itrig <<std::endl;
    
    chain->GetEntry(itrig);
    for(Int_t chan=0; chan<NumOfChans; chan++){
      Int_t mcpid = nMcpId[chan];
      Double_t t = fLeadingEdge[chan];
      Int_t row = nRow[chan]-1;
      Int_t col = nCol[chan]-1;
      Int_t pixid = 8*row + col;

      //      t = m*fTot[chan] + q;
	      
      // bad pixels
      if(mcpid==2  && pixid==55) continue;
      if(mcpid==2  && pixid==62) continue;
      if(mcpid==10  && pixid==49) continue;
      if(mcpid==10  && pixid==53) continue;
      if(mcpid==10  && pixid==58) continue;

      if(mcpid==10  && pixid==57) continue;
      if(mcpid==14  && pixid==27) continue;
      if(mcpid==14  && pixid==35) continue;
      

      if(bValid[chan] && mcpid<15 && (uniqmcp == mcpid || uniqmcp<0)){
	dmcpid = (uniqmcp<0) ? mcpid : 0; 
	hWcorr[mcpid][pixid]->SetTitle(Form("ch %d",chan));
	hWcorr[mcpid][pixid]->Fill(t-fLeadingEdge[trigger], fTot[chan]);
	hWcorr0[mcpid][pixid]->Fill(t-fLeadingEdge[trigger], fTot[trigger]);
	hWcorr1[mcpid][pixid]->Fill(t, fTot[chan]);
	hWcorr2[mcpid][pixid]->Fill(fLeadingEdge[trigger], fTot[trigger]);
	dataLeg[dmcpid][pixid][idata[dmcpid][pixid]]=t-fLeadingEdge[trigger];
	dataTo1[dmcpid][pixid][idata[dmcpid][pixid]]=fTot[chan];
	dataTo2[dmcpid][pixid][idata[dmcpid][pixid]]=fTot[trigger];
	idata[dmcpid][pixid]++;
	//fhDigi[mcpid]->Fill(7-pixid/8, pixid%8);
	fhDigi[mcpid]->Fill(col, row);
      	hPTime[mcpid][pixid]->Fill(t-fLeadingEdge[trigger]);
	hPTime1[mcpid][pixid]->Fill(t);
	hPTime2[mcpid][pixid]->Fill(fLeadingEdge[trigger]);
      	//hTime[m][nRow[chan]-1+8*(8-nCol[chan])]->Fill(t-fLeadingEdge[nChannelDiff]);
      }
    }
  }
  spect = new TSpectrum(2);
  TString fitres = "";
  TString path = createDir("rdata/ttt/"+dirid, "beam_data "+fileid, fSaveFlag); 
  writeInfo("digi.csv", drawDigi("m,p,v\n",1), fSaveFlag);
  Int_t angle = 0;
  for(Int_t m=0; m<nmcp; m++){
    for(Int_t p=0; p<npix; p++){
      if(hPTime[m][p]->GetEntries()<1 && (uniqmcp != m && uniqmcp>0)) continue;
      if(fSaveFlag>0) {
  	docorrections(m,p);
  	if(abs(gfitted0)<0.08 || abs(gfitted1)<0.08 || gfitted0>5 || abs(gfitted1)>5) fitres += Form("<tr><td>%d</td><td>%d</td><td>-</td><td>-</td><td>-</td></tr>\n",m,p);
  	else fitres += Form("<tr><td>%d</td><td>%d</td><td>%f</td><td>%f</td><td>%f</td></tr>\n",m,p,abs(gfitted0), abs(gfitted1), (1-abs(gfitted1)/abs(gfitted0))*100);
	
      }
      //hPTime[m][p]->SetTitle(Form("Theta %d, mcp %d, pixel %d", angle, m, p));
      //fit(hPTime[m][p]);
      cTime->cd(1);
      hPTime[m][p]->Draw();
      // cTime->Modified(); cTime->Update(); cTime->WaitPrimitive();
      save(cTime,path,Form("time_ang%dmcp%dpix%d",angle,m,p),fInfo,fSaveFlag,1); //cLeadingEdge
    }
  }
  std::cout<<"fitres  "<<fitres <<std::endl;
  

  cDigi->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
  		 "exec3event(Int_t,Int_t,Int_t,TObject*)");
}

Bool_t lock = false;
void exec3event(Int_t event, Int_t gx, Int_t gy, TObject *selected){
  TCanvas *c = (TCanvas *) gTQSender;
  TPad *pad = (TPad *) c->GetSelectedPad();
  if (!pad) return;
  Float_t x = pad->AbsPixeltoX(gx);
  Float_t y = pad->AbsPixeltoY(gy);
  x = pad->PadtoX(x);
  y = pad->PadtoY(y);
  if(event ==1 && lock) lock = false;
  else if(event ==1) lock = true;
  if(lock) return;

  if (selected->InheritsFrom(TH2::Class())){
    TH2F *hDigi = (TH2F *) selected;
    Int_t binx = hDigi->GetXaxis()->FindBin(x);
    Int_t biny = hDigi->GetYaxis()->FindBin(y);
    TString nmcp = selected->GetName();
    nmcp = nmcp(3,nmcp.Sizeof());
    Int_t mcpid = nmcp.Atoi();
    Int_t pixid = 8*(biny-1)+binx-1;
    printf("Canvas %s: event=%d, x=%d, y=%d, p=%d, selected=%d\n", nmcp.Data(), event, binx, biny, pixid,nmcp.Atoi());
    docorrections(mcpid,pixid);
  }
}

