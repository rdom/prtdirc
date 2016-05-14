#include "../../prttools/prttools.C"
void da_cangle(TString inFile = "../build/reco_spr.root"){
  fSavePath = "data/cangle";
  TChain ch("dirc"); ch.Add(inFile);
  Double_t cangle,spr,trr,nph,par1,par2,par3,par4,par5,par6,test1,test2,theta,phi; 
  Int_t tofPid;
  
  ch.SetBranchAddress("spr",&spr);
  ch.SetBranchAddress("tofPid",&tofPid);
  ch.SetBranchAddress("trr",&trr);
  ch.SetBranchAddress("nph",&nph);
  ch.SetBranchAddress("cangle",&cangle);
  ch.SetBranchAddress("par4",&par4);
  ch.SetBranchAddress("par5",&par5);
  ch.SetBranchAddress("par6",&par6);
  ch.SetBranchAddress("test1",&test1);
  ch.SetBranchAddress("test2",&test2);
  ch.SetBranchAddress("theta",&theta);
  ch.SetBranchAddress("phi",&phi);
  
  TH1F *hCangle2 = new TH1F("hCangle2",";#theta_{c} [degree];entries [#]",100,0.8,0.84);
  TH1F *hCangle3 = new TH1F("hCangle3",";#theta_{c} [degree];entries [#]",100,0.8,0.84);

  Int_t nent = ch.GetEntries();
  std::cout<<"# entries  "<< nent <<std::endl;
  for (Int_t i = 0; i < nent; i++) {
    ch.GetEvent(i);
    if(tofPid==211) hCangle2->Fill(cangle);
    if(tofPid==2212) hCangle3->Fill(cangle);
    if(tofPid==211)std::cout<<"cangle "<< cangle<<std::endl;
    
  }

  prt_normalize(hCangle2,hCangle3);
  gStyle->SetOptFit(1);
  canvasAdd("hCangle");


  hCangle2->Fit("gaus");
  hCangle2->SetLineColor(4);
  hCangle3->Fit("gaus");
  hCangle3->SetLineColor(1);

  hCangle2->Draw();
  hCangle3->Draw("same");
  
  canvasSave(1,0);
}
