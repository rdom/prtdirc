// -----------------------------------------
// PrtLutReco.cpp
//
// Created on: 13.07.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------

#include "PrtLutReco.h"
#include "PrtManager.h"
#include "PrtLutNode.h"

#define prt__sim
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

using std::cout;
using std::endl;

TH1F*  fHist0 = new TH1F("timediff",";t_{measured}-t_{calculated} [ns];entries [#]", 500,-10,10);
TH1F*  fHist0d = new TH1F("timediffd",";t_{measured}-t_{calculated} [ns];entries [#]", 500,-10,10);
TH1F*  fHist0r = new TH1F("timediffr",";t_{measured}-t_{calculated} [ns];entries [#]", 500,-10,10);
TH1F*  fHist0s = new TH1F("timediffs",";t_{measured}-t_{calculated} [ns];entries [#]", 500,-10,10);
TH1F*  fHist0i = new TH1F("timediffi",";t_{measured}-t_{calculated} [ns];entries [#]", 500,-10,10);

TH1F*  fhNph_pi = new TH1F("fhNph_pi",";detected photons [#];entries [#]", 150,0,150);
TH1F*  fhNph_p = new TH1F("fhNph_p",";detected photons [#];entries [#]", 150,0,150);
TH1F*  fHist1 = new TH1F("time1",";measured time [ns];entries [#]",   1000,0,50);
TH1F*  fHist2 = new TH1F("time2",";calculated time [ns];entries [#]", 1000,0,50);
TH1F*  fHist6 = new TH1F("time6",";measured time [ns];entries [#]", 1000,0,50);

TH1F*  hBounce = new TH1F("bounce",";number of bounces [#];photons per event [#]", 150,0,150);

TH2F*  fDiff = new TH2F("diff",";measured time [ns];t_{measured}-t_{calculated} [ns]", 300,0,30,150,-5,5);
TH2F*  fHist3 = new TH2F("time3",";calculated time [ns];measured time [ns]", 500,0,80, 500,0,40);
TH2F*  fHist4 = new TH2F("time4",";#theta_{c}sin(#varphi_{c});#theta_{c}cos(#varphi_{c})", 100,-1,1, 100,-1,1);
TH2F*  fHist5 = new TH2F("time5",";#theta_{c}sin(#varphi_{c});#theta_{c}cos(#varphi_{c})", 100,-1,1, 100,-1,1);

TH1F *hLnDiffGr4 = new TH1F("hLnDiffGr4",";ln L(p) - ln L(#pi);entries [#]",120,-60,60);
TH1F *hLnDiffGr2 = new TH1F("hLnDiffGr2",";ln L(p) - ln L(#pi);entries [#]",120,-60,60);
TH1F *hLnDiffTi4 = new TH1F("hLnDiffTi4",";ln L(p) - ln L(#pi);entries [#]",120,-60,60);
TH1F *hLnDiffTi2 = new TH1F("hLnDiffTi2",";ln L(p) - ln L(#pi);entries [#]",120,-60,60);
TH2F *hLnMap = new TH2F("hLnMap",";GR     ln L(p) - ln L(#pi);TI     ln L(p) - ln L(#pi); ",120,-60,60,120,-60,60);

TH2F *hChrom = new TH2F("chrom",";t_{measured}-t_{calculated} [ns];#Delta#theta_{C} [mrad]", 100,-1.5,1.5, 100,-30,30);
TH2F *hChromL = new TH2F("chroml",";(t_{measured}-t_{calculated})/t_{measured};#Delta#theta_{C} [mrad]", 100,-0.15,0.15, 100,-30,30);

const int nphi=80, ntheta=40;
TH2F *hLutCorrD = new TH2F("hLutCorrD",";#theta_{l}sin(#varphi_{l});#theta_{l}cos(#varphi_{l})",200,-1,1,200,-1,1);
TH2F *hLutCorrC = new TH2F("hLutCorrC",";#theta_{l}sin(#varphi_{l});#theta_{l}cos(#varphi_{l})",100,-1,1,100,-1,1);

int gg_i(0), gg_ind(0);
TGraph gg_gr;
PrtLutNode *fLutNode[prt_maxdircch];

TH1F*  fHistMcp[15];
TH1F*  fHistCh[prt_maxdircch];

//cluster search
int mcpdata[15][65];
int cluster[15][65];
int lneighbours[65];
int lsize(0);

// -----   Default constructor   -------------------------------------------
PrtLutReco::PrtLutReco(TString infile, TString lutfile, int verbose){
  fVerbose = verbose;
  fChain = new TChain("data");
  fChain->Add(infile);
  fEvent = new PrtEvent();
  fChain->SetBranchAddress("PrtEvent", &fEvent);

  fChain->SetBranchStatus("fHitArray.fParentParticleId", 0);
  fChain->SetBranchStatus("fHitArray.fNreflectionsInPrizm", 0);
  fChain->SetBranchStatus("fHitArray.fCherenkovMC", 0);
 
  fFile = new TFile(lutfile);
  fTree=(TTree *) fFile->Get("prtlut") ;
  fLut = new TClonesArray("PrtLutNode");
  fTree->SetBranchAddress("LUT",&fLut); 
  fTree->GetEntry(0);

  fHist = new TH1F("cherenkov_angle_hist",  "cherenkov angle;#theta_{C} [rad];entries [#]", 150,0.6,1); //150
  fHistPi = new TH1F("cherenkov_angle_hist_Pi",  "cherenkov angle pi;#theta_{C} [rad];entries [#]", 150,0.6,1); //150
  fHisti = new TH1F("cherenkov_angle_histi","cherenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150
  fFit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
  fFit->SetNpx(300);
  fSpect = new TSpectrum(10);
  fRadiator=1;
  fStudyId = 420;
  fMethod = PrtManager::Instance()->GetRunType();
  
  if(infile.Contains("beam_")){
    TString fileid(infile);
    fileid.Remove(0,fileid.Last('/')+1);
    fileid.Remove(fileid.Last('.')-1);
    prt_data_info = getDataInfo(fileid);
    fRadiator =  prt_data_info.getRadiatorId();
    fStudyId = prt_data_info.getStudyId();
    
    TString opath(infile);
    opath.Remove(opath.Last('/'));
    if(infile.Contains("C.root")){
      prt_savepath = opath+Form("/%dr/%d",prt_data_info.getStudyId(),prt_data_info.getFileId());
    }else{
      prt_savepath = opath+Form("/%ds/%d",prt_data_info.getStudyId(),prt_data_info.getFileId());
    }
    
  }else prt_savepath="data/sim";
  std::cout<<"prt_savePath  "<< prt_savepath <<std::endl;    
 
  for(int i=0; i<prt_maxdircch; i++){
    fLutNode[i] = (PrtLutNode*) fLut->At(i);
  }
  cout << "-I- PrtLutReco: Intialization successfull" << endl;

  for(int i=0; i<prt_nmcp; i++){
    fHistMcp[i] = new TH1F(Form("fHistMcp_%d",i),Form("fHistMcp_%d;#theta_{C} [rad];entries [#]",i), 40,-0.05,0.05); //150
  }

  for(int i=0; i<prt_maxdircch; i++){
    fHistCh[i] = new TH1F(Form("fHistCh_%d",i),Form("fHistCh_%d;#theta_{C} [rad];entries [#]",i), 150,0.6,1); //150
  }
  
  // read corrections
  fCorrPath = PrtManager::Instance()->GetOutName()+"_corr.root";    
  for(int i=0; i<prt_nmcp; i++) fCorr[i] = 0;
  if(!gSystem->AccessPathName(fCorrPath)){
    std::cout<<"------- reading  "<<fCorrPath <<std::endl;
    int pmt;
    double corr;
    TChain ch("corr"); ch.Add(fCorrPath);
    ch.SetBranchAddress("pmt",&pmt);
    ch.SetBranchAddress("corr",&corr);
    for(int i=0; i<ch.GetEntries(); i++){
      ch.GetEvent(i);
      fCorr[pmt] = (fabs(corr)<0.011)? corr: 0.00001;
      std::cout<<"pmt "<<pmt<<"  "<<corr<<std::endl;    
    }
  }else{
    std::cout<<"------- corr file not found  "<<fCorrPath <<std::endl;
  }

  fTimeImaging = (fMethod == 4)? true : false;
  fPdfPath = infile.ReplaceAll(".root",".pdf1.root");
  if(fMethod == 2) { // read pdf
    if(!gSystem->AccessPathName(fPdfPath)){
      std::cout<<"------- reading  "<<fPdfPath <<std::endl;
      TFile pdfFile(fPdfPath);
      double sigma = 250; // ps
      int binfactor=(int)(sigma/50.+0.1);
      for(int i=0; i<prt_maxdircch; i++){  
	auto hpdf2 = (TH1F*)pdfFile.Get(Form("hs_%d",i));
	auto hpdf4 = (TH1F*)pdfFile.Get(Form("hf_%d",i));
	if(sigma > 0) hpdf2->Rebin(binfactor);
	if(sigma > 0) hpdf4->Rebin(binfactor);

	fPdf2[i] = new TGraph(hpdf2);
	fPdf4[i] = new TGraph(hpdf4);
	fTimeImaging = true;
      }
    }    
  }

  for(int i=0; i<prt_maxdircch; i++){
    fTime2[i] = new TH1F(Form("hs_%d",i),"pdf;LE time [ns]; entries [#]",1000,0,50);
    fTime4[i] = new TH1F(Form("hf_%d",i),"pdf;LE time [ns]; entries [#]",1000,0,50);
  }
  std::cout<<"initialization done "<<std::endl;
  
}

// -----   Destructor   ----------------------------------------------------
PrtLutReco::~PrtLutReco(){

}

//-------------- Loop over tracks ------------------------------------------
void PrtLutReco::Run(int start, int end){
  TVector3 dird, dir, momInBar,momInBar0,posInBar,cz;
  double mom,tangle,likelihood(0),boxPhi,weight,evtime,bartime,lenz,posz,dirz,luttheta,
    barHitTime, hitTime,angdiv,dtheta,dtphi,prtangle;
  int  tofPid(0),distPid(0),likePid(0),pdgcode, evpointcount(0),total2(0),total4(0);
  int events[5]={0};
  bool reflected = kFALSE;
  gStyle->SetOptFit(111);

  TVector3 fnX1 = TVector3 (1,0,0);   
  TVector3 fnY1 = TVector3( 0,1,0);
  bool testTrRes = false;

  TString outFile = PrtManager::Instance()->GetOutName()+".root";
  double theta(0),phi(0),cangle(0),spr(0), trr(0),nph(0),nph_err(0),
    cangle_pi(0),spr_pi(0),trr_pi(0),nph_pi(0),nph_pi_err(0),
    par1(0), par2(0), par3(0), par4(0), par5(0), par6(0), test1(0), test2(0), test3(0),
    sep_gr(0),sep_gr_err(0),sep_ti(0),sep_ti_err(0),
    beamx(0),beamz(0),nnratio_p(0),nnratio_pi(0),timeRes(0);
  double minChangle(0);
  double maxChangle(1);  
  double criticalAngle = asin(1.00028/1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
  double radiatorL = (fRadiator==2)? 1224.9 : 1200; //plate : bar
  double radiatorW = (fRadiator==2)? 174.8 : 34.9;  //plate : bar
  double radiatorH = (fRadiator==2)? 1224.9 : 17.1;  //plate : bar

  prt_setRootPalette(1);
  prt_createMap();
  prt_initDigi();

  outFile.ReplaceAll("reco_",Form("reco_%d_",prt_data_info.getFileId()));
  TFile file(outFile,"recreate");
  TTree tree("reco","PrtLutReco");
  tree.Branch("mom", &mom,"mom/D");
  tree.Branch("theta",&theta,"theta/D");
  tree.Branch("phi",&phi,"phi/D");
  tree.Branch("nph",&nph,"nph/D");
  tree.Branch("nph_err",&nph_err,"nph_err/D");  
  tree.Branch("sep_gr",&sep_gr,"sep_gr/D");
  tree.Branch("sep_gr_err",&sep_gr_err,"sep_gr_err/D");
  tree.Branch("sep_ti",&sep_ti,"sep_ti/D");
  tree.Branch("sep_ti_err",&sep_ti_err,"sep_ti_err/D");
  tree.Branch("time_res",&timeRes,"timeRes/D");

  tree.Branch("cangle",&cangle,"cangle/D");
  tree.Branch("spr", &spr,"spr/D");
  tree.Branch("trr", &trr,"trr/D");  

  tree.Branch("nph_pi",&nph_pi,"nph_pi/D");
  tree.Branch("nph_pi_err",&nph_pi_err,"nph_pi_err/D");
  tree.Branch("cangle_pi",&cangle_pi,"cangle_pi/D");
  tree.Branch("spr_pi", &spr_pi,"spr_pi/D");
  tree.Branch("trr_pi", &trr_pi,"trr_pi/D");  
  
  tree.Branch("pid_tof", &tofPid,"tofPid/I");
  tree.Branch("pid_dist", &distPid,"distPid/I");
  tree.Branch("pid_lh", &likePid,"likePid/I");
  tree.Branch("likelihood",&likelihood,"likelihood/D");
  
  tree.Branch("test1",&test1,"test1/D");
  tree.Branch("test2",&test2,"test2/D");
  tree.Branch("test3",&test3,"test3/D");
  tree.Branch("nnratio_p",&nnratio_p,"nnratio_p/D");
  tree.Branch("nnratio_pi",&nnratio_pi,"nnratio_pi/D");
  tree.Branch("beamx",&beamx,"beamx/D");
  tree.Branch("beamz",&beamz,"beamz/D");
  tree.Branch("par5",&par5,"par5/D");
  tree.Branch("par6",&par6,"par6/D");  
  
  test1 = PrtManager::Instance()->GetTest1();
  test2 = PrtManager::Instance()->GetTest2();
  test3 = PrtManager::Instance()->GetTest3();
  par5 = PrtManager::Instance()->GetPrismStepX();
  par6 = PrtManager::Instance()->GetPrismStepY();
  timeRes = PrtManager::Instance()->GetTimeRes();
  fMethod = PrtManager::Instance()->GetRunType();

  int nEvents = fChain->GetEntries();
  if(end==0) end = nEvents;
  
  int pdfend = 50000;
  if(fPdfPath.Contains("S.pdf1.root")) pdfend = 5000;
  
  if(fMethod == 4) {
    start = pdfend;
    pdfend = nEvents;
    end = nEvents;
  }
  
  std::cout<<"Run started for ["<<start<<","<<end <<"]"<<std::endl;
  int nsHits(0),nsEvents(0),nHits(0), ninfit(1);

  if(start<0) {
    ninfit=abs(start);
    start=0;
  }
  
  for (int ievent = start; ievent < nEvents && (events[2] < end || events[4] < end) && ievent < pdfend; ievent++){
    int nhhits(0);
    fChain->GetEntry(ievent);
    nHits = fEvent->GetHitSize();
    bool bsim = fEvent->GetType();
    
    if(ievent%1000==0) std::cout<<"Event # "<< ievent << " has "<< nHits <<" hits "<< events[2]<<" "<<events[4]<<std::endl;	

    if(bsim) {
      posz = 0.5*radiatorL-fEvent->GetPosition().Z() + fRand.Uniform(-5,5);
    }
    else posz = fEvent->GetPosition().Z();

    double momentum = fEvent->GetMomentum().Mag();
    if(bsim) momentum /= 1000;
    tofPid = fEvent->GetParticle();
    int pid = prt_get_pid(tofPid);
    if(events[pid]>=end) continue;

    double speed = 197.0; // mm/ns
    double sigma[]={0,0,0.0081,0,0.0081};
    double angle1(0), angle2(0),sum1(0),sum2(0),sumti(0),sumti2(0),sumti4(0),range(5*sigma[2]),noise(0.2); //0.0082
        
    if(ievent-start==0){
      tree.SetTitle(fEvent->PrintInfo());
      prtangle = fEvent->GetAngle() + test1*TMath::RadToDeg(); // prt_data_info.getAngle();
      phi = fEvent->GetPhi() + test2*TMath::RadToDeg(); //prt_data_info.getPhi();
      mom = fEvent->GetMomentum().Mag();
      beamx = fEvent->GetPosition().X();
      beamz = fEvent->GetPosition().Z();
      if(bsim) beamz = 0.5*radiatorL-beamz;
      
      for(int i: {2,4}){
	fAngle[i] = acos(sqrt(momentum*momentum+ prt_mass[i]*prt_mass[i])/momentum/1.4725); //1.4738 = 370 = 3.35 // 1.4725 = 380 = 3.26
	fFunc[i] = new TF1(Form("gaus_%d",i),"[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
	fFunc[i]->SetParameter(0,1);
	fFunc[i]->SetParameter(1,fAngle[i]);
	fFunc[i]->SetParameter(2,sigma[i]);
      }
 
      std::cout<<fStudyId<<"  prtangle++  "<<prtangle<< " phi "<<phi<<" t1 "<<test1<<" t2 "<<test2<<std::endl;

      momInBar = TVector3(0,0,1);
      momInBar.RotateY(TMath::Pi() - prtangle*TMath::DegToRad());
      momInBar.RotateZ(phi*TMath::DegToRad());
      momInBar0 = momInBar;
    }

    // //smear track
    // momInBar = momInBar0;
    // double smearangle = 0.005; // prt_rand.Gaus(0,test1);
    // momInBar.RotateY(smearangle);
    // momInBar.Rotate(prt_rand.Uniform(0,TMath::TwoPi()),momInBar0);

    if(fVerbose==3){
      cz = momInBar.Unit();
      cz = TVector3(-cz.X(),cz.Y(),cz.Z());
    }
    
    // if(fMethod==2 && tofPid!=2212) continue;
	
    if(!bsim){
      int gch, ndirc(0), t2(0), t3h(0), t3v(0),
	str1(0),stl1(0),str2(0),stl2(0);
      int hodo1(0), hodo2(0);
      if(fabs(fEvent->GetMomentum().Mag()-7)<0.1){
	double tof = fEvent->GetTest1();
        if(fStudyId==403 && fMethod != 4){
	  if( pid == 4 && tof < 34.2 ) continue;
	  if( pid == 2 && tof > 33.3 ) continue;
	}
	if(fStudyId == 420 && fMethod != 4){
	  if( pid == 4 && tof < 36.9 ) continue;
	  if( pid == 2 && tof > 35.7 ) continue;
	}
	if(fStudyId == 420 && fMethod == 4){
	  if( pid == 4 && tof < 36.6 ) continue;
	  if( pid == 2 && tof > 36.0 ) continue;
	}
      }
      for(int h=0; h<nHits; h++) {
      	gch = fEvent->GetHit(h).GetChannel();
	if(gch<prt_maxdircch) ndirc++;

	if(gch==513) t2++;
	if(gch==514) t3h++;
	if(gch==515) t3v++;
	if(gch>=1089 && gch<=1106) hodo1++;
	//if(gch>=1094 && gch<=1101) hodo1++;
	//if(gch>=1097 && gch<=1098) hodo1++;
	if(gch==1140) str1++;
	if(gch==1142) stl1++;
	if(gch==1144) str2++;
	if(gch==1146) stl2++;
	
	//if(gch>=1115 && gch<=1120)
	hodo2++;      
      }
      if(ndirc<5) continue;
      if(!(t2 && t3h && t3v && hodo1 && hodo2)) continue;
      if(!(str1 && stl1 && str2 && stl2)) continue;
      //if(!(t3h && t3v)) continue;
    }
    
    // SearchClusters();
    
    //event t0 smearing
    double t0smear = fRand.Gaus(0,0.15);

    if(!bsim && fStudyId==403){      
      posz -= 15;
    }
    
    for(int h=0; h<nHits; h++) {
      fHit = fEvent->GetHit(h);
      hitTime = fHit.GetLeadTime();
      if(bsim) hitTime += fRand.Gaus(0,0.25) + t0smear; // time resol. in case it was not simulated
      else{
	if(fStudyId == 420) hitTime += 0.62;
	if(fStudyId == 403){
	  double o = 0.1;
	  if(fabs(prtangle-20)<1) o = -0.2;
	  if(fabs(prtangle-25)<1) o = -0.05;
	  if(fabs(prtangle-30)<1) o =  0.08;
	  if(fabs(prtangle-35)<1) o =  0.2;
	  if(fabs(prtangle-40)<1) o =  0.2;
	  if(fabs(prtangle-45)<1) o =  0.2;
	  if(fabs(prtangle-50)<1) o =  0.1;
	  if(fabs(prtangle-55)<1) o =  0.1;
	  if(fabs(prtangle-60)<1) o =  -0.1;
	  if(fabs(prtangle-65)<1) o =  -0.2;
	  if(fabs(prtangle-70)<1) o = -0.3;
	  if(fabs(prtangle-75)<1) o = -0.3;
	  if(fabs(prtangle-80)<1) o = -0.4;
	  if(fabs(prtangle-85)<1) o = -0.4;
	  if(fabs(prtangle-90)<1) o = -0.4;
	  if(fabs(prtangle-95)<1) o = -0.4;
	  if(fabs(prtangle-100)<1) o = -0.4;
	  if(fabs(prtangle-105)<1) o = -0.35;
	  if(fabs(prtangle-110)<1) o = -0.3;
	  if(fabs(prtangle-115)<1) o = -0.1;
	  if(fabs(prtangle-120)<1) o =  0.0;
	  if(fabs(prtangle-120)<1) o =  0.0;
	  if(fabs(prtangle-130)<1) o = +0.1;
	  if(fabs(prtangle-135)<1) o = +0.3;
	  if(fabs(prtangle-140)<1) o = +0.4;
	  hitTime += o;
	}
      }
      
      //======================================== dynamic cuts
      {
	{ //time cuts
	  if(prtangle<=75){
	    if(hitTime<7 || hitTime>45) continue;
	    reflected = kTRUE;
	  }else if(prtangle>105){
	    if(hitTime<3 || hitTime>40) continue;
	    reflected = kFALSE;
	  }else{
	    if(hitTime<11)  reflected = kFALSE; //13.5
	    else reflected = kTRUE;
	  }
	}
      }
      //==================================================
      
      if(fVerbose==3){
	// TVector3 cd = fHit.GetMomentum();
	// fHist5->Fill(cd.Theta()*TMath::Sin(cd.Phi()),cd.Theta()*TMath::Cos(cd.Phi()));
      }

      // TVector3 vv = fHit.GetMomentum();
      // vv.RotateY(prtangle*deg);
      // dirz = vv.Z();
      // if(dirz<0) reflected = kTRUE;
      // else reflected = kFALSE;

      int pixid=fHit.GetPixelId()-1;
      int mcpid=fHit.GetMcpId();
      int ch = map_mpc[mcpid][pixid];

      // if(reflected) continue;
      if(reflected) lenz = 2*radiatorL - posz;
      else lenz = posz;
      
      if(prt_isBadChannel(ch)) continue;
      int nedge=GetEdge(mcpid, pixid);
      // if(cluster[mcpid][pixid]>4) continue;
      
      bool isGoodHit_gr(0), isGoodHit_ti(0);
      int bestbounce = 0;      
      double besttangle = 0, besttdiff = 100;
      int size = fLutNode[ch]->Entries();
      
      for(int i=0; i<size; i++){
	weight = 1;//fLutNode[ch]->GetWeight(i);
	dird   = fLutNode[ch]->GetEntryCs(i,nedge); // nedge=0
        //dird  = fLutNode[ch]->GetEntry(i);
	evtime = fLutNode[ch]->GetTime(i);
	
	int pathid = fLutNode[ch]->GetPathId(i);
	bool samepath(false);
	if(bsim && pathid==fHit.GetPathInPrizm()) samepath=true;
	//if(fLutNode[ch]->GetNRefl(i)!=1 ) continue;
	//if(!samepath) continue;
 
	double lphi = dird.Phi();
	double ltheta = dird.Theta();
	if(lphi<0) lphi = TMath::TwoPi()+lphi;
	if(ltheta>TMath::PiOver2()) ltheta = TMath::Pi()-ltheta;	  
	int iphi = nphi*(lphi)/TMath::TwoPi();
	int itheta = ntheta*(ltheta)/TMath::PiOver2();
	
	for(int u=0; u<4; u++){
	  if(u == 0) dir = dird;
	  if(u == 1) dir.SetXYZ( -dird.X(), dird.Y(), dird.Z());
	  if(u == 2) dir.SetXYZ( dird.X(), -dird.Y(), dird.Z()); //no need when no divergence in vertical plane
	  if(u == 3) dir.SetXYZ( -dird.X(),-dird.Y(), dird.Z()); //no need when no divergence in vertical plane
	  if(reflected) dir.SetXYZ( dir.X(), dir.Y(), -dir.Z());
	  if(dir.Angle(fnX1) < criticalAngle || dir.Angle(fnY1) < criticalAngle) continue;

	  luttheta = dir.Theta();  
	  if(luttheta > TMath::PiOver2()) luttheta = TMath::Pi()-luttheta;

	  double len = lenz/cos(luttheta);
	  
	  // if(nx%2==0 && u==1 && 2*lenz/nx>60) continue;
	  // if(nx%2==0 && u==3 && 2*lenz/nx>60) continue;	  
	  // if(ny%2==0 && u==2 && 2*lenz/ny>150) continue;
	  // if(ny%2==0 && u==3 && 2*lenz/ny>150) continue; 
	  
	  bartime = fabs(len/speed);	  
	  double luttime = bartime+evtime;
	  double tdiff = hitTime-luttime;
	  fHist0->Fill(tdiff);
	  if(reflected) fHist0r->Fill(tdiff);
	  else fHist0d->Fill(tdiff);
	  if(samepath) fHist0i->Fill(tdiff);
	  
	  tangle = momInBar.Angle(dir)+fCorr[mcpid];
	  if(fabs(prtangle-90)<16){
	    if(reflected) if(fabs(tdiff)<1.5)  tangle -= 0.0075*tdiff; // chromatic correction
	    if(!reflected) if(fabs(tdiff)<1.5) tangle -= 0.006*tdiff; // chromatic correction	  
	  }else{
	    if(fabs(tdiff/hitTime)<0.15) tangle -= 0.035*tdiff/hitTime;
	  }

	  if(fabs(tdiff)<1.5)  tangle -= 0.002*tdiff; // chromatic correction

	  hChrom->Fill(tdiff,(tangle-fAngle[pid])*1000);
	  hChromL->Fill(tdiff/hitTime,(tangle-fAngle[pid])*1000);
	  
	  if(fabs(tdiff)<1.5+luttime*0.04 && fabs(tangle-0.815)<0.05) isGoodHit_ti = true;
	  if(fabs(tdiff)>timeRes+luttime*0.04) continue;

	  fDiff->Fill(hitTime,tdiff);
	  fHist3->Fill(fabs(luttime),hitTime);	  	  
	   
	  double w = 1; //2000/len; //w=0.5*(2-fabs(tdiff));
	  if(tangle > minChangle && tangle < maxChangle && tangle < 1.85){
	    if(tofPid==211 && fMethod==2) fHistPi->Fill(tangle ,weight);
	    else fHist->Fill(tangle, weight);

	    fHistMcp[mcpid]->Fill(tangle-fAngle[pid] ,weight);
	    fHistCh[ch]->Fill(tangle ,weight);

	    if(fabs(tangle-0.815)<0.05){
	      isGoodHit_gr = true;

	      sum1 += w*TMath::Log(fFunc[4]->Eval(tangle)+noise);
	      sum2 += w*TMath::Log(fFunc[2]->Eval(tangle)+noise);	   

	      fHist0s->Fill(tdiff);
	      double lenx = len*dir.X();
	      double leny = len*dir.Y();
	      int nx = round(lenx/radiatorH);
	      int ny = round(leny/radiatorW);
	      if(fabs(tdiff)<fabs(besttdiff)) {
		besttdiff = tdiff;
		bestbounce = nx+ny;
		besttangle = tangle;
	      }

	      fHist2->Fill(luttime);
	      hLutCorrD->Fill(ltheta*TMath::Sin(lphi),ltheta*TMath::Cos(lphi));
	    }
	    if(fVerbose==3){
	      TVector3 rdir = TVector3(-dir.X(),dir.Y(),dir.Z());
	      rdir.RotateUz(cz);	      
	      double llphi = rdir.Phi();
	      double tt =  rdir.Theta();
	      fHist4->Fill(tt*TMath::Sin(llphi),tt*TMath::Cos(llphi));

	      //for cherenckov circle fit
	      gg_gr.SetPoint(gg_i,tt*TMath::Sin(llphi),tt*TMath::Cos(llphi));
	      gg_i++;
	    }
	  }
	}
      }

      if(fTimeImaging && isGoodHit_ti){

	if(fMethod == 2){
	  double t = hitTime;
	  if(fabs(besttdiff) < 0.3) t -= besttdiff;
	  double noiseti = 1e-5;
	  double aminti4 = fPdf4[ch]->Eval(t);
	  double aminti2 = fPdf2[ch]->Eval(t);
	  sumti4 += TMath::Log((aminti4+noiseti));
	  sumti2 += TMath::Log((aminti2+noiseti));
	}

	if(fMethod == 4){ // create pdf 		  
	  double t = hitTime;
	  if(fabs(besttdiff) < 0.3) t -= besttdiff;
	  if(pid == 4){
	    total4++;
	    fTime4[ch]->Fill(t);
	  }
	  if(pid == 2){
	    total2++;
	    fTime2[ch]->Fill(t);
	  }	  
	}
      }

      if(isGoodHit_gr){
	hBounce->Fill(bestbounce);
	fHist1->Fill(hitTime);
	nhhits++;
	nsHits++;
	if(tofPid == 211) prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);	
      }
    }

    if(fMethod == 2 && fTimeImaging){ // time imaging

      sumti = 2*(sumti4-sumti2);
      if(sumti != 0){
	if(tofPid == 2212) hLnDiffTi4->Fill(sumti);       
	if(tofPid == 211) hLnDiffTi2->Fill(sumti);
      }      
    }    
    
    double sum = sum1-sum2;
    if(sum != 0){
      if(tofPid == 2212){
	fhNph_p->Fill(nhhits);
	hLnDiffGr4->Fill(sum);
      }
      if(tofPid == 211){
	fhNph_pi->Fill(nhhits);
	hLnDiffGr2->Fill(sum);
      }
      likelihood = sum;
      events[pid]++;
    }

    if(pid == 4 && sumti != 0) hLnMap->Fill(sum,sumti);
    
    // if(fVerbose==1){
    //   prt_canvasAdd("ff",800,400);
    //   fFunc[4]->Draw();
    //   fFunc[2]->SetLineColor(4);
    //   fFunc[2]->Draw("same");
      
    //   prt_waitPrimitive("ff");
    //   prt_canvasDel("ff");
    //   //prt_canvasSave(1,0);
    //   //prt_canvasDel(Form("lh_%d",gg_ind));
    // }
	
    if(fVerbose>0 &&  fMethod==3 && nsEvents%ninfit==0){
      if(nsHits>5){
        // if(tofPid==2212 && sum > 0){
	//   std::cout<<"p  "<<sum1 << "   pi "<<sum2 << "  s "<< sum<<std::endl;
	//   if(fVerbose>0)  if(!FindPeak(cangle,spr, prtangle, tofPid)) continue;
	// }

	FindPeak(cangle,spr,cangle_pi,spr_pi, prtangle, tofPid);	
	distPid = FindPdg(momentum,cangle);
	nph = nsHits/(double)ninfit;
	spr = spr*1000;
	trr = spr/sqrt(nph);
	theta = fEvent->GetAngle();
	par3 = fEvent->GetTest1();
	tree.Fill();
      }
      ResetHists();
      nsHits=0;
    }

  }

  if(fMethod == 4){ // create pdf
    std::cout<<"saving pdfs into "<<fPdfPath<<std::endl;
 
    TFile efile(fPdfPath,"RECREATE");
    for(int i=0; i<prt_maxdircch; i++){
      fTime2[i]->Scale(1/(Double_t)total2);
      fTime4[i]->Scale(1/(Double_t)total4);
      fTime2[i]->Write();
      fTime4[i]->Write();
    }
    efile.Write();
    efile.Close();
    std::cout<<"total2 "<<total2 <<" total4 "<<total4<<std::endl;
  }

  if(fMethod==2){
    TF1 *ff;
    gROOT->SetBatch(1);
    if(fhNph_p->GetEntries()>20){
      TFitResultPtr r = fhNph_p->Fit("gaus","SQN","",0,120);
      nph = r->Parameter(1);
      nph_err = r->ParError(1);
    }
    if(fhNph_pi->GetEntries()>20){
      TFitResultPtr r = fhNph_pi->Fit("gaus","SQN","",0,120);
      nph_pi = r->Parameter(1);
      nph_pi_err = r->ParError(1);
    }
    //nph = prt_fit(fhNph_pi,40,10,50,1).X();
    gROOT->SetBatch(0);
    FindPeak(cangle,spr,cangle_pi,spr_pi, prtangle);
    //nph = nsHits/(double)nsEvents;
    spr = spr*1000;
    trr = spr/sqrt(nph);
    spr_pi = spr_pi*1000;
    trr_pi = spr_pi/sqrt(nph_pi);
    theta = fEvent->GetAngle();
    par3 = fEvent->GetTest1();

    double m1,m2,s1,s2,dm1,dm2,ds1,ds2;; 
    if(hLnDiffGr4->GetEntries()>10){
      hLnDiffGr4->Fit("gaus","Q");
      ff = hLnDiffGr4->GetFunction("gaus");
      m1=ff->GetParameter(1);
      s1=ff->GetParameter(2);
      dm1=ff->GetParError(1);
      ds1=ff->GetParError(2);
    }
    if(hLnDiffGr2->GetEntries()>10){
      hLnDiffGr2->Fit("gaus","Q");
      ff = hLnDiffGr2->GetFunction("gaus");
      m2=ff->GetParameter(1);
      s2=ff->GetParameter(2);
      dm2=ff->GetParError(1);
      ds2=ff->GetParError(2);
    }
    sep_gr = (fabs(m2-m1))/(0.5*(s1+s2));

    double e1,e2,e3,e4;
    e1=2/(s1+s2)*dm1;
    e2=2/(s1+s2)*dm2;
    e3=-((2*(m1 + m2))/((s1 + s2)*(s1 + s2)))*ds1;
    e4=-((2*(m1 + m2))/((s1 + s2)*(s1 + s2)))*ds2;
    sep_gr_err=sqrt(e1*e1+e2*e2+e3*e3+e4*e4);

    if(fTimeImaging){
      if(hLnDiffTi4->GetEntries()>10){
	hLnDiffTi4->Fit("gaus","Q");
	ff = hLnDiffTi4->GetFunction("gaus");
	m1=ff->GetParameter(1);
	s1=ff->GetParameter(2);
	dm1=ff->GetParError(1);
	ds1=ff->GetParError(2);
      }
      if(hLnDiffTi2->GetEntries()>10){
	hLnDiffTi2->Fit("gaus","Q");
	ff = hLnDiffTi2->GetFunction("gaus");
	m2=ff->GetParameter(1);
	s2=ff->GetParameter(2);
	dm2=ff->GetParError(1);
	ds2=ff->GetParError(2);
      }
      sep_ti = (fabs(m2-m1))/(0.5*(s1+s2));
    
      e1=2/(s1+s2)*dm1;
      e2=2/(s1+s2)*dm2;
      e3=-((2*(m1 + m2))/((s1 + s2)*(s1 + s2)))*ds1;
      e4=-((2*(m1 + m2))/((s1 + s2)*(s1 + s2)))*ds2;
      sep_ti_err=sqrt(e1*e1+e2*e2+e3*e3+e4*e4);
    }
    

    nnratio_pi = fhNph_pi->GetEntries()/(double)end;
    nnratio_p = fhNph_p->GetEntries()/(double)end;

    if(fVerbose) {
      std::cout<<"nnratio_pi "<<nnratio_pi<<" "<<end <<"  "<< fhNph_pi->GetEntries()<<std::endl;	
      std::cout<<Form("p  SPR = %2.2F N = %2.2f +/- %2.2f",spr,nph,nph_err)<<std::endl;
      std::cout<<Form("pi SPR = %2.2F N = %2.2f +/- %2.2f",spr_pi,nph_pi,nph_pi_err)<<std::endl;
      std::cout<<"sep gr "<< sep_gr <<" +/- "<<sep_gr_err <<std::endl;
      std::cout<<"sep ti "<< sep_ti <<" +/- "<<sep_ti_err <<std::endl;

      double a = prtangle;
      TString nid = "";//Form("_%2.0f",a);

      // plots      
      { // cherenkov angle
	prt_canvasAdd("tangle"+nid,800,400);
      
	fHist->SetTitle(Form("theta %3.1f , TOF PID = %d", a, tofPid));
	fHist->SetMinimum(0);
	//fHist->Scale(1/fHist->GetMaximum());

	prt_normalize(fHist,fHistPi);
	fHistPi->SetLineColor(4);
	fHist->SetLineColor(1);
	fHist->Draw();
	fHistPi->Draw("same");
	// fFunc[4]->Draw("same");
	// fFunc[2]->Draw("same");
	fHisti->SetLineColor(kRed+2);
	if(fHisti->GetEntries()>5) fHisti->Draw("same");

	TLegend *leg = new TLegend(0.13,0.68,0.39,0.82);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->AddEntry((TObject*)0,Form("#sigma_{p} = %2.2f mrad",spr),"");
	leg->AddEntry((TObject*)0,Form("#sigma_{#pi} = %2.2f mrad",spr_pi),"");
	leg->Draw();
	  
	drawTheoryLines();
      }
      
      { // time
	prt_canvasAdd("time",800,400);
	prt_normalize(fHist1,fHist2);
	fHist1->SetTitle(Form("theta %3.1f", a));
	fHist2->SetLineColor(2);
	fHist1->Draw();
	fHist2->Draw("same");

	// fHist1->Scale(1/4000.);
	// fHist1->Rebin(10);
	// fHist1->Draw("p hist");
	  
	prt_canvasAdd("diff_time",800,400);
	fDiff->Draw("colz");

	prt_canvasAdd("time_diff"+nid,800,400);
	fHist0->SetTitle(Form("theta %3.1f", a));
	fHist0->SetLineColor(kBlack);		
	fHist0->Draw();
	fHist0d->SetLineColor(kGreen+1);
	fHist0d->Draw("same");
	fHist0r->SetLineColor(kBlue+1);
	fHist0r->Draw("same");
	fHist0s->SetLineColor(kOrange+6);
	fHist0s->Draw("same");
	fHist0i->SetLineColor(kRed+1);
	if(fHist0i->GetEntries()>5)  fHist0i->Draw("same"); 

	// prt_canvasAdd("cm"+nid,800,400);
	// fHist3->SetTitle(Form("theta %3.1f", a));
	// fHist3->Draw("colz");
      }

      { // likelihood
	prt_canvasAdd("lhood_gr",800,400);
	prt_normalize(hLnDiffGr4,hLnDiffGr2);
	hLnDiffGr4->SetLineColor(2);
	hLnDiffGr4->SetName(Form("s_%2.2f",sep_gr));
	hLnDiffGr4->Draw();
	hLnDiffGr2->SetLineColor(4);
	hLnDiffGr2->Draw("same");

	if(fTimeImaging){
	  prt_canvasAdd("lhood_ti",800,400);
	  prt_normalize(hLnDiffTi4,hLnDiffTi2);
	  hLnDiffTi4->SetLineColor(2);
	  hLnDiffTi4->SetName(Form("s_%2.2f",sep_ti));
	  hLnDiffTi4->Draw();
	  hLnDiffTi2->SetLineColor(4);
	  hLnDiffTi2->Draw("same");
	}
	
	prt_canvasAdd("lhood_map",800,800);
	hLnMap->SetStats(0);
	hLnMap->Draw("colz");
      }
      
      { // bounce
	prt_canvasAdd("bounce",800,400);
	hBounce->Scale(1/4000.);
	hBounce->Draw("hist");
      }
	
      { // lut corrections
	prt_canvasAdd("lutcorrd"+nid,600,600);
	TGaxis::SetMaxDigits(3);
	hLutCorrD->SetStats(0);
	hLutCorrD->Draw("colz");
      }

      { // chromatic corrections
	prt_canvasAdd("chrom"+nid,800,400);
	hChrom->SetStats(0);
	hChrom->Draw("colz");
	prt_canvasAdd("chroml"+nid,800,400);
	hChromL->SetStats(0);
	hChromL->Draw("colz");
      }

      { // nph
	prt_canvasAdd("nph"+nid,800,400);
	fhNph_p->SetLineColor(kRed+1);
	fhNph_p->Draw();
	fhNph_pi->SetLineColor(kBlue+1);
	fhNph_pi->Draw("same");
      }
	
      { // hp
	auto cdigi = prt_drawDigi(2018);
	cdigi->SetName("hp"+nid);
	prt_canvasAdd(cdigi);
      }
	
      if(false){
	int tmax, max=0;
	for(int m=0; m<prt_nmcp;m++){
	  prt_hdigi[m]->Rebin2D(8,8);
	  prt_hdigi[m]->GetXaxis()->SetNdivisions(0);
	  prt_hdigi[m]->GetYaxis()->SetNdivisions(0);
	  prt_hdigi[m]->GetXaxis()->SetTickLength(0);
	  prt_hdigi[m]->GetYaxis()->SetTickLength(0);
	  prt_hdigi[m]->GetXaxis()->SetAxisColor(1);
	  prt_hdigi[m]->GetYaxis()->SetAxisColor(1);
	  prt_hdigi[m]->SetMarkerSize(10);
	  tmax = prt_hdigi[m]->GetMaximum();
	  if(max<tmax) max = tmax;
	}
	for(int m=0; m<prt_nmcp;m++){	  
	  prt_hdigi[m]->Scale(1/(double)max);
	}
      }
        
      if(fVerbose==3){
	TCanvas* cCher = new TCanvas("cCher","cCher",0,0,800,400);
	cCher->Divide(2,1);
	cCher->cd(1);
     
	fHist4->SetStats(0);
	fHist4->SetTitle(Form("Calculated from LUT, #theta = %3.1f#circ", a));
	fHist4->Draw("colz");
	double x0(0), y0(0);
	FitRing(x0,y0,cangle);
	TVector3 corr(x0,y0,1-TMath::Sqrt(x0*x0+y0*y0));
	std::cout<<"Tcorr "<< corr.Theta()*1000<< "  Pcorr "<< corr.Phi() <<std::endl;

	TLegend *leg = new TLegend(0.5,0.7,0.85,0.87);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->AddEntry((TObject*)0,Form("Entries %0.0f",fHist4->GetEntries()),"");
	leg->AddEntry((TObject*)0,Form("#Delta#theta_{c} %f [mrad]",corr.Theta()*1000),"");
	leg->AddEntry((TObject*)0,Form("#Delta#varphi_{c} %f [mrad]",corr.Phi()),"");
	leg->Draw();

	TArc *arc = new TArc(x0,y0,cangle);
	arc->SetLineColor(kRed);
	arc->SetLineWidth(1);
	arc->SetFillStyle(0);
	arc->Draw();
	gg_i=0;
	gg_gr.Set(0);

	cCher->cd(2);
	gStyle->SetOptStat(1110); 
	fHist5->SetTitle(Form("True from MC, #theta = %3.1f#circ", a));
	fHist5->Draw("colz");
	
	prt_canvasAdd(cCher);
      }
      gROOT->SetBatch(0);
    }
  }

  tree.Fill();
  file.Write();
  
  prt_set_style();
  if(fVerbose>1) prt_waitPrimitive("time","");
  if(fVerbose) {
    prt_canvasSave(1,0,true);
    ResetHists(); 
  }  
}

int g_num =0;
bool PrtLutReco::FindPeak(double& cangle, double& spr,double& cangle_pi, double& spr_pi, double a, int tofpdg){
  cangle=0;
  spr=0;
  gROOT->SetBatch(1);

  if(fHist->GetEntries()>20 || fHistPi->GetEntries()>20){
    // int nfound = fSpect->Search(fHist,1,"",0.9); //0.6
    // if(nfound>0) cangle = fSpect->GetPositionX()[0];
    // else cangle =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());
    
    cangle = fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());
    if(cangle>0.84) cangle=0.82;
    fFit->SetParameters(100,cangle,0.008);
    fFit->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");      
    fFit->SetParLimits(0,10,1E6);
    fFit->SetParLimits(1,cangle-0.03,cangle+0.03); 
    fFit->SetParLimits(2,0.005,0.016);
    
    if(fMethod==3){
      fFit->SetParameter(2,0.007);
      fFit->SetParLimits(2,0.005,0.008);
      gROOT->SetBatch(1);
      fHist->Fit("fgaus","MQ","",0.6,1);
      cangle = fFit->GetParameter(1);
      spr = fFit->GetParameter(2);
      gROOT->SetBatch(!fVerbose);
      
      if(fVerbose>2){
	gROOT->SetBatch(0);
	TString nid = Form("casphi_%2.3f_mrad",PrtManager::Instance()->GetTest1());
	prt_canvasAdd(nid,800,400);
	fHist->Draw();
	drawTheoryLines();
	prt_waitPrimitive(nid);
	// prt_canvasSave(1,0);
      }
    }
    
    if(fMethod==2){
      gROOT->SetBatch(1);
        
      fHist->Fit("fgaus","MQ","",cangle-0.06,cangle+0.06);
      if(fVerbose>1) gROOT->SetBatch(0);
      cangle = fFit->GetParameter(1);
      spr = fFit->GetParameter(2);
      
      fHistPi->Fit("fgaus","MQ","",cangle-0.06,cangle+0.06);        
      cangle_pi = fFit->GetParameter(1);
      spr_pi = fFit->GetParameter(2);      
      	
      // for(int i=0; i<prt_maxdircch; i++){
      // 	prt_canvasAdd(Form("r_tangle_ch_%d",i),800,400);
      // 	fHistCh[i]->Fit("fgaus","lq","",fAngle[2]-0.03,fAngle[2]+0.03);
      // 	std::cout<<"if(ch=="<< i<<") tangle += "<<fAngle[2]-fFit->GetParameter(1)<<";" <<std::endl;	
      // 	fHistCh[i]->Draw();
      // } 
      
      { // corrections
	if(fabs(fCorr[0])<0.00000001 && fabs(fCorr[7])<0.00000001){
	  std::cout<<"Writing "<<fCorrPath<<std::endl;
	  
	  TFile fc(fCorrPath,"recreate");
	  TTree *tc = new TTree("corr","corr");
	  int pmt;
	  double corr;
	  tc->Branch("pmt",&pmt,"pmt/I");
	  tc->Branch("corr",&corr,"corr/D");
	
	  fFit->SetParameter(1,0);    // mean
	  fFit->SetParLimits(1,-0.012,0.012);
	  fFit->SetParLimits(2,0.006,0.010); // width		
	  for(int i=0; i<prt_nmcp; i++){
	    if(fHistMcp[i]->GetEntries()<10000) continue;
	    if(fVerbose>2) prt_canvasAdd(Form("r_tangle_%d",i),800,400);
	    fHistMcp[i]->Fit("fgaus","MQ","",-0.03,0.03);
	    pmt = i;
	    corr= -fFit->GetParameter(1);
	    tc->Fill();
	    std::cout<<"if(mcpid=="<< i<<") tangle += "<<corr<<";" <<std::endl;	  
	    if(fVerbose>2){
	      fHistMcp[i]->Draw();
	      drawTheoryLines();
	    }
	  }
	  
	  tc->Write();
	  fc.Write();
	  fc.Close();

	  fFit->ReleaseParameter(1);
	  fFit->ReleaseParameter(2);
	}
      }
            
    }
  }
  
  return (cangle>0 && cangle<1);
}

void PrtLutReco::ResetHists(){
  fHist->Reset();
  fHisti->Reset();
  fHist0->Reset();
  fHist0i->Reset();
  fHist1->Reset();
  fHist2->Reset();
  fHist3->Reset();
  fHist4->Reset();
  for(int m=0; m<prt_nmcp;m++) prt_hdigi[m]->Reset();
}

TF1 *lFit = new TF1("lgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.6,0.9);
TF1 *lFitPi = new TF1("lgausPi","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.6,0.9);

double PrtLutReco::fillLnDiffPPi(double cangle, int tofPid, double mom){
  if(fHist->GetEntries()>20 ){
    double angle1(0), angle2(0), sigma(0.006),range(0.015);

    // //fHist->Scale(1/fHist->GetMaximum());

    // double d1,d2, sum1(0),sum2(0);
    // int sbin = fHist->FindBin(fAngle[4]-range);
    // int ebin = fHist->FindBin(fAngle[4]+range); 
    // // fHist->GetXaxis()->GetNbins()
    // for(int i=sbin; i< ebin; i++){
    //   if(fHist->GetBinContent(i) < 0.01 ) continue;
    //   d1 = fFunc[4]->Eval(fHist->GetBinCenter(i))- fHist->GetBinContent(i);
    //   d2 = fFunc[2]->Eval(fHist->GetBinCenter(i))- fHist->GetBinContent(i);

    //   std::cout<<"f1 "<< fFunc[4]->Eval(fHist->GetBinCenter(i)) << "   f2 "<<fFunc[2]->Eval(fHist->GetBinCenter(i)) << "    v "<< fHist->GetBinContent(i) <<std::endl;
    
    //   // if(d1>0) sum1+=TMath::Log(d1);
    //   // if(d2>0) sum2+=TMath::Log(d2);
    //   sum1+=TMath::Log(fabs(d1));
    //   sum2+=TMath::Log(fabs(d2));

    // }
  
    lFit->SetRange(fAngle[4]-range,fAngle[2]+range);
    lFit->FixParameter(0,fHist->GetMaximum()-0.5);
    lFit->FixParameter(1,fAngle[4]);
    lFit->FixParameter(2,0.01);
    lFit->FixParameter(3,0);
    lFit->FixParameter(4,0.5);
 

    fHist->Fit("lgaus","Q","",fAngle[4]-range,fAngle[2]+range);
    double amin,amin2,edm,errdef;
    int nvpar,nparx;
    TVirtualFitter *fitter = TVirtualFitter::Fitter(fHist);
    fitter->GetStats(amin,edm,errdef,nvpar,nparx);
  
    // lFitPi->SetRange(fAnglePi-range,fAnglePi+range);
    // lFitPi->SetLineColor(4);
    // lFitPi->FixParameter(0,fFit->GetParameter(0));
    // lFitPi->FixParameter(1,fAnglePi);
    // lFitPi->FixParameter(2,sigma);
    // lFitPi->FixParameter(3,fFit->GetParameter(3));
    // lFitPi->FixParameter(4,fFit->GetParameter(4));

    lFitPi->SetRange(fAngle[4]-range,fAngle[2]+range);
    lFitPi->SetLineColor(4);
    lFitPi->FixParameter(0,fHist->GetMaximum()-0.5);
    lFitPi->FixParameter(1,fAngle[2]);
    lFitPi->FixParameter(2,0.01);
    lFitPi->FixParameter(3,0);
    lFitPi->FixParameter(4,0.5);

    fHist->Fit("lgausPi","Q","",fAngle[4]-range,fAngle[2]+range);
    fitter = TVirtualFitter::Fitter(fHist);
    fitter->GetStats(amin2,edm,errdef,nvpar,nparx);
  
    if(fVerbose) printf("tofPid %04d | %1.4f (%1.4f/%1.4f) likelihood is %1.2f/%1.2f \n",tofPid,cangle,fAngle[4],fAngle[2], amin, amin2);
    gg_ind++;

    if(fVerbose==1){
      prt_canvasAdd("ff",800,400);
      //prt_canvasAdd(Form("lh_%d",gg_ind),800,400);
      fHist->SetTitle(Form("%d",tofPid));
      fHist->Draw();
      lFit->SetLineColor(2);
      lFit->Draw("same");
      // fFunc[4]->Draw("same");
      // fFunc[2]->SetLineColor(4);
      // fFunc[2]->Draw("same");
     
      //if(fabs(amin-amin2)<5)
      prt_waitPrimitive("ff");
      prt_canvasDel("ff");
      //prt_canvasSave(1,0);
      //prt_canvasDel(Form("lh_%d",gg_ind));
    }
  
    return amin-amin2;
  }
  return 1000;
}

double PrtLutReco::fillLnDiffPPi2(double cangle, int tofPid, double mom){
  if(fHist->GetEntries()>20 ){
    double angle1(0), angle2(0), sigma(0.006),range(0.03);

    double d1,d2, sum1(0),sum2(0);
    int sbin = fHist->FindBin(fAngle[2]-range);
    int ebin = fHist->FindBin(fAngle[4]+range); 
    for(int i=sbin; i< ebin; i++){
      if(fHist->GetBinContent(i)<1 ) continue;
      d1 = 10*fabs(fHist->GetBinContent(i) *(fAngle[4]  - fHist->GetBinCenter(i)));
      d2 = 10*fabs(fHist->GetBinContent(i) *(fAngle[2] - fHist->GetBinCenter(i)));
      if(d1>0 && d2>0){
	std::cout<<"d1  "<<d1 << "   d2    "<< d2 <<std::endl;
	sum1+=TMath::Log(d1);
	sum2+=TMath::Log(d2);
      }
    }
  
    if(fVerbose) printf("tofPid %04d | %1.4f (%1.4f/%1.4f) likelihood is %1.2f/%1.2f \n",tofPid,cangle,fAngle[4],fAngle[2], sum1, sum2);
    gg_ind++;

    if(fVerbose==1){
      prt_canvasAdd("ff",800,400);
      //prt_canvasAdd(Form("lh_%d",gg_ind),800,400);
      fHist->SetTitle(Form("%d",tofPid));
      fHist->Draw();
      lFit->SetLineColor(2);
      lFit->Draw("same");
      // gFp->Draw("same");
      // gFpi->SetLineColor(4);
      // gFpi->Draw("same");
     
      //if(fabs(amin-amin2)<5)
      prt_waitPrimitive("ff");
      prt_canvasDel("ff");
      //prt_canvasSave(1,0);
      //prt_canvasDel(Form("lh_%d",gg_ind));
    }
  
    return sum1-sum2;
  }
  return 1000;
}

void circleFcn(int &, double *, double &f, double *par, int) {
  int np = gg_gr.GetN();
  f = 0;
  double *x = gg_gr.GetX();
  double *y = gg_gr.GetY();
  for (int i=0;i<np;i++) {
    double u = x[i] + par[0];
    double v = y[i] + par[1];
    double dr = par[2] - TMath::Sqrt(u*u+v*v);
    f += dr*dr;  
  }
  std::cout<<"fcn  "<< f<<std::endl;
  
}

void circleFcn2(int &, double *, double &f, double *par, int) {
  int np = gg_gr.GetN();
  f = 0;
  double *x = gg_gr.GetX();
  double *y = gg_gr.GetY();
  for (int i=0;i<np;i++){
    double u = x[i] + par[0];
    double v = y[i] + par[1];
    double dr = par[2] - TMath::Sqrt(u*u+v*v);
    if(dr>0.07) f += dr*dr; 
    else f += fabs(dr);
  }
}

void PrtLutReco::FitRing(double& x0, double& y0, double& theta){
  TGraph ff_gr;
  int ff_i(0);
  int np = gg_gr.GetN();
  double *x = gg_gr.GetX();
  double *y = gg_gr.GetY();
  for (int i=0;i<np;i++) {
    if( fabs(theta - TMath::Sqrt(x[i]*x[i]+y[i]*y[i]))<0.05) {
      ff_gr.SetPoint(ff_i,x[i],y[i]);
      ff_i++;
    }
  }
  gg_gr = ff_gr;

  //Fit a circle to the graph points
  TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
  TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
  fitter->SetPrecision(0.00000001);
  fitter->SetMaxIterations(1000);

  fitter->SetFCN(circleFcn);
  fitter->SetParameter(0, "x0",   0.03, 0.01, -0.05,0.05);
  fitter->SetParameter(1, "y0",   0, 0.01, -0.05,0.05);
  fitter->SetParameter(2, "R",    theta, 0.01, theta-0.05,theta+0.05);

  //fitter->FixParameter(0);
  //fitter->FixParameter(1);
  fitter->FixParameter(2);
  double arglist[1] = {0};
  fitter->ExecuteCommand("MINIMIZE", arglist, 0);

  // fitter->SetFCN(circleFcn2);
  // fitter->ExecuteCommand("MINIMIZE", arglist, 0);

  x0 = fitter->GetParameter(0);
  y0 = fitter->GetParameter(1);
  theta = fitter->GetParameter(2);
}

int PrtLutReco::FindPdg(double mom, double cangle){
  cangle =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());
  
  int pdg[]={211,2212};
  double mass[] = {0.139570,0.9382723};
  double tdiff, diff=100;
  int minid=0;
  for(int i=0; i<2; i++){
    tdiff = fabs(cangle - acos(sqrt(mom*mom + mass[i]*mass[i])/mom/1.46907)); //1.46907 - fused silica
    if(tdiff<diff){
      diff = tdiff;
      minid = i;
    }
  }
  return pdg[minid]; 
}

int PrtLutReco::GetEdge(int mcpid, int pixid){
  // int x(0),y(0), piid(pixid) , nedge(0); //new
  // for(int h=0; h<fEvent->GetHitSize(); h++) {
  // 	int pid=fEvent->GetHit(h).GetPixelId();
  // 	int mid=fEvent->GetHit(h).GetMcpId();
  // 	double tdif=fabs(hitTime-fEvent->GetHit(h).GetLeadTime());
  // 	if(mid!=mcpid || pid==piid || tdif>0.3) continue;
  // 	if(pid==piid-1 && piid%8!=0) y-=1;
  // 	if(pid==piid+1 && piid%8!=7) y+=1;

  // 	if(pid==piid+8 && piid<57) x-=1;
  // 	if(pid==piid-8 && piid>8)  x+=1;
  // }

  int x(0),y(0),piid(pixid+1),nedge(0); //old
  for(int h=0; h<fEvent->GetHitSize(); h++) {
    int pid=fEvent->GetHit(h).GetPixelId();
    int mid=fEvent->GetHit(h).GetMcpId();
    if(mid!=mcpid || pid==piid) continue;
    if(pid==piid-1 && piid%8!=1) x-=1;
    if(pid==piid+1 && piid%8!=0) x+=1;

    if(pid==piid+8 && piid<57) y+=1;
    if(pid==piid-8 && piid>8)  y-=1;
  }
      
  if(x== 0 && y== 0) nedge=0;
  if(x==-1 && y== 0) nedge=1;
  if(x==-1 && y== 1) nedge=2;
  if(x== 0 && y== 1) nedge=3;
  if(x== 1 && y== 1) nedge=4;
  if(x== 1 && y== 0) nedge=5;
  if(x== 1 && y==-1) nedge=6;
  if(x== 0 && y==-1) nedge=7;
  if(x==-1 && y==-1) nedge=8;
  
  return nedge;
}

int getneighbours(int m, int p){
  for(int i=0; i<65; i++) if(p==lneighbours[i]) return -1;
  lneighbours[lsize]=p;
  lsize++;
  for(int t=0; t<65; t++){
    if(mcpdata[m][t]){
      for(int i=0; i<65; i++) if(t==lneighbours[i]) continue;
      if((t==p-1 && p%8!=0) || (t==p+1 && p%8!=7) ||
	 (t==p+8 && p<57) || (t==p-8 && p>8)) getneighbours(m,t);
    }
  }
  return lsize;
}

void getclusters(){
  for(int m=0; m<prt_nmcp; m++){
    for(int p=0; p<65; p++){
      if(mcpdata[m][p])  cluster[m][p] = getneighbours(m,p);
      lsize=0;
      for(int i=0; i<65; i++) lneighbours[i]=0;
    }
  }
}
  
void PrtLutReco::SearchClusters(){

  for(int j=0; j<prt_nmcp; j++){
    for(int i=0; i<65; i++){
      mcpdata[j][i]=0;
      cluster[j][i]=0;
    }
  }
  
  for(int h=0; h<fEvent->GetHitSize(); h++) {
    int mid=fEvent->GetHit(h).GetMcpId();
    int pid=fEvent->GetHit(h).GetPixelId()-1;
    mcpdata[mid][pid]=1;
  }
  getclusters();
}

void PrtLutReco::drawTheoryLines(){
  gPad->Update();
  TLine *line = new TLine(0,0,0,1000);
  line->SetX1(fAngle[4]);
  line->SetX2(fAngle[4]);
  line->SetY1(gPad->GetUymin());
  line->SetY2(gPad->GetUymax());
  line->SetLineColor(kRed);
  line->Draw();

  TLine *line1 = new TLine(0,0,0,1000);
  line1->SetX1(fAngle[2]);
  line1->SetX2(fAngle[2]);
  line1->SetY1(gPad->GetUymin());
  line1->SetY2(gPad->GetUymax());
  line1->SetLineColor(kBlue);
  line1->Draw();
}
