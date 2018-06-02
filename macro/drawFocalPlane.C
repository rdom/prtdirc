#define prt__sim

#include "../src/PrtEvent.h"
#include "../src/PrtHit.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

double findVertex(TVector3 v1,TVector3 m1, TVector3 v2,TVector3 m2, TVector3* newvertex){
  TVector3 head0 = v1;
  TVector3 tail0 = v1+100*m1;

  TVector3 head1 = v2;
  TVector3 tail1 = v2+100*m2;

  TVector3 dir0 = head0 - tail0;
  TVector3 dir1 = head1 - tail1;

  TVector3 diff = tail0 - tail1;

  double a = dir0.Dot(dir0);//always >= 0
  double b = dir0.Dot(dir1);
  double c = dir1.Dot(dir1);//always >= 0
  double d = dir0.Dot(diff);
  double e = dir1.Dot(diff);
  double f = a * c - b * b;//always >= 0

  double sc;
  double tc;

  if(f<0.000001){//The lines are almost parallel
    sc = 0.0f;
    tc = b > c ? d / b : e / c;//Use the largest denominator
  }else{
    sc = (b * e - c * d) / f;
    tc = (a * e - b * d) / f;
  }
  
  TVector3 point0 = tail0 + dir0 * sc;
  TVector3 point1 = tail1 + dir1 * tc;
  newvertex->SetXYZ((point0.X()+point1.X())/2.,(point0.Y()+point1.Y())/2.,(point0.Z()+point1.Z())/2.);
  return (point1-point0).Mag();
}

//void drawFocalPlane(TString infile="../build/focalplane.root", Double_t r1 = 47.8, Double_t r2 = 29.1, Int_t it1=0, Int_t it2=0, Double_t energy=-1){
void drawFocalPlane(TString infile="../build/focalplane.root", double sep=100, Double_t r1 = 47.8, Double_t r2 = 29.1, Int_t it1=0, Int_t it2=0, Double_t energy=-1){

  if(!prt_init(infile,1, "data/fp_tmp")) return;

  Double_t lensThickness(13.12); //RMI lens
  Double_t radiatorL(1200.7); //bar
  //Double_t radiatorL(1224.9); //plate
  TVector3 res;
  TH2F *hFp1 = new TH2F("hFp1",Form("r_{1}=%2.2f    r_{2}=%2.2f;x [cm];y [cm]",r1,r2),500,0,50,600,-30,30 );
  TH2F *hFp2 = new TH2F("hFp2",Form("r_{1}=%2.2f    r_{2}=%2.2f;z [cm];y [cm]",r1,r2),500,-30,30,500,-30,50 );
  
  PrtHit hit[2];
  int count(0);
  for (int ievent=0; ievent<prt_entries; ievent++){
    prt_nextEvent(ievent,1000);
    int nhits = prt_event->GetHitSize();
    if(nhits!=2) continue;
    for(int h=0; h<nhits; h++) hit[h] = prt_event->GetHit(h);
    TVector3 sourcePos =  hit[0].GetLocalPos();
    
    double d = findVertex(hit[0].GetGlobalPos(),hit[0].GetMomentum().Unit(),hit[1].GetGlobalPos(),hit[1].GetMomentum().Unit(), &res);
    if(d<0.01){
      double x = -(res.X()+0.5*radiatorL+lensThickness)/10.;
      if(x>5){
	hFp1->Fill(x,res.Z()/10.);
	hFp2->Fill(res.Y()/10.,res.Z()/10.);
      }
      // x = -(sourcePos.X()+0.5*radiatorL+lensThickness)/10.;
      // hFp1->Fill(x,sourcePos.Z()/10.);
            
      count++;
    }
  }

  
  double eff = 100*count/(double)prt_entries;
  TString senergy = "";
  if(energy>-1){
    double m = 1, cm = 0.01, nanometer = 0.000000001, GeV = 1000000000;
    double LambdaE = 2.0 * 3.14159265358979323846 * 1.973269602e-16 * m * GeV;
    double lam = LambdaE/(energy*nanometer);

    senergy = Form("   E=%2.2f [eV] / #lambda=%2.0f [nm]",energy,lam);
  }
  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  hFp1->SetMinimum(2);
  hFp1->SetTitle(Form("d = %2.1f mm",sep));
  //hFp1->SetTitle(Form("r_{1}=%2.2f    r_{2}=%2.2f   #varepsilon=%2.0f",r1,r2,eff)+senergy);
  //prt_canvasAdd(Form("fp_%d_%d",it1,it2),800,500);
  prt_canvasAdd(Form("fp_%2.1f",sep),800,500);
  hFp1->Draw("col");
  // prt_canvasAdd(Form("fp2_%d_%d",it1,it2),600,800);
  // hFp2->Draw("colz");
  //prt_canvasSave(0,"drawFocalPlane.C",1,"fp");
  prt_canvasSave(0);
  
}
