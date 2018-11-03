// -----------------------------------------
// PrtLutReco.cpp
//
// Created on: 13.07.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------

#include "PrtLutReco.h"
#include "PrtManager.h"
#include "PrtLutNode.h"
#include "PrtTrackInfo.h"
#include "PrtPhotonInfo.h"
#include "PrtAmbiguityInfo.h"

#include "TVirtualFitter.h"
#include "TArc.h"
#include "CLHEP/Units/SystemOfUnits.h"

#define prt__sim
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

using std::cout;
using std::endl;

TH1F*  fHist0 = new TH1F("timediff",";t_{calc}-t_{measured} [ns];entries [#]", 500,-10,10);
TH1F*  fHist0i = new TH1F("timediffi",";t_{calc}-t_{measured} [ns];entries [#]", 500,-10,10);
TH1F*  fhNph = new TH1F("fhNph",";detected photons [#];entries [#]", 150,0,150);
TH1F*  fhNph_pi = new TH1F("fhNph_pi",";detected photons [#];entries [#]", 150,0,150);
TH1F*  fhNph_p = new TH1F("fhNph_p",";detected photons [#];entries [#]", 150,0,150);
TH1F*  fHist1 = new TH1F("time1",";measured time [ns];entries [#]",   1000,0,100);
TH1F*  fHist2 = new TH1F("time2",";calculated time [ns];entries [#]", 1000,0,100);
TH1F*  fHist6 = new TH1F("time6",";measured time [ns];entries [#]", 1000,0,100);

TH2F*  fHist3 = new TH2F("time3",";calculated time [ns];measured time [ns]", 500,0,80, 500,0,40);
TH2F*  fHist4 = new TH2F("time4",";#theta_{c}sin(#varphi_{c});#theta_{c}cos(#varphi_{c}", 100,-1,1, 100,-1,1);
TH2F*  fHist5 = new TH2F("time5",";#theta_{c}sin(#varphi_{c});#theta_{c}cos(#varphi_{c}", 100,-1,1, 100,-1,1);

TH1F *hLnDiffP = new TH1F("hLnDiffP",  ";ln L(p) - ln L(#pi);entries [#]",90,-30,30);
TH1F *hLnDiffPi = new TH1F("hLnDiffPi",";ln L(p) - ln L(#pi);entries [#]",90,-30,30);

TF1 *gF1 = new TF1("gaus0","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
TF1 *gF2= new TF1("gaus0","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);

Int_t gg_i(0), gg_ind(0);
TGraph gg_gr;
PrtLutNode *fLutNode[5000];

TH1F*  fHistMcp[15];
TH1F*  fHistCh[960];

//cluster search
Int_t mcpdata[15][65];
Int_t cluster[15][65];
Int_t lneighbours[65];
Int_t lsize(0);

// -----   Default constructor   -------------------------------------------
PrtLutReco::PrtLutReco(TString infile, TString lutfile, Int_t verbose){
  fVerbose = verbose;
  fChain = new TChain("data");
  fChain->Add(infile);
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
  fSpect = new TSpectrum(10);
  fRadiator=1;
  
  if(infile.Contains("beam_")){
    TString fileid(infile);
    fileid.Remove(0,fileid.Last('/')+1);
    fileid.Remove(fileid.Last('.')-1);
    prt_data_info = getDataInfo(fileid);
    fRadiator =  prt_data_info.getRadiatorId();
    
    TString opath(infile);
    opath.Remove(opath.Last('/'));
    if(infile.Contains("C.root")){
      prt_savepath = opath+Form("/%dr/%d",prt_data_info.getStudyId(),prt_data_info.getFileId());
    }else{
      prt_savepath = opath+Form("/%ds/%d",prt_data_info.getStudyId(),prt_data_info.getFileId());
    }
    
  }else prt_savepath="data/sim";
  std::cout<<"prt_savePath  "<< prt_savepath <<std::endl;    

  for(Int_t i=0; i<5000; i++){
    fLutNode[i] = (PrtLutNode*) fLut->At(i);
  }
  cout << "-I- PrtLutReco: Intialization successfull" << endl;

  for(Int_t i=0; i<prt_nmcp; i++){
    fHistMcp[i] = new TH1F(Form("fHistMcp_%d",i),Form("fHistMcp_%d;#theta_{C} [rad];entries [#]",i), 160,0.6,1); //150
  }

  for(Int_t i=0; i<960; i++){
    fHistCh[i] = new TH1F(Form("fHistCh_%d",i),Form("fHistCh_%d;#theta_{C} [rad];entries [#]",i), 100,0.6,1); //150
  }

}

// -----   Destructor   ----------------------------------------------------
PrtLutReco::~PrtLutReco(){

}

//-------------- Loop over tracks ------------------------------------------
void PrtLutReco::Run(Int_t start, Int_t end){
  TVector3 dird, dir, momInBar(0,0,1),posInBar,cz;
  Double_t mom, cangle,spr,tangle,likelihood(0),boxPhi,weight,evtime,bartime,lenz,posz,dirz,luttheta, barHitTime, hitTime;
  Int_t  tofPid(0),distPid(0),likePid(0),pdgcode, evpointcount=0;
  int events[5]={0};
  Bool_t reflected = kFALSE;
  gStyle->SetOptFit(111);

  TVector3 fnX1 = TVector3 (1,0,0);   
  TVector3 fnY1 = TVector3( 0,1,0);
  bool testTrRes = false;
  Double_t angdiv,dtheta,dtphi,prtangle;

  TString outFile = PrtManager::Instance()->GetOutName()+".root";
  Double_t theta(0),phi(0), trr(0),  nph(0),nph_err(0),
    par1(0), par2(0), par3(0), par4(0), par5(0), par6(0), test1(0), test2(0), test3(0),
    separation(0),beamx(0),beamz(0),nnratio(0),nnratio_p(0),nnratio_pi(0),timeRes(0);
  Double_t minChangle(0);
  Double_t maxChangle(1);  
  Double_t criticalAngle = asin(1.00028/1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
  Double_t radiatorL = (fRadiator==2)? 1224.9 : 1200; //plate : bar

  prt_setRootPalette(1);
  prt_createMap();
  prt_initDigi();

  outFile.ReplaceAll("reco_",Form("reco_%d_",prt_data_info.getFileId()));
  TFile file(outFile,"recreate");
  TTree tree("dirc","SPR");
  tree.Branch("mom", &mom,"mom/D");
  tree.Branch("tofPid", &tofPid,"tofPid/I");
  tree.Branch("distPid", &distPid,"distPid/I");
  tree.Branch("likePid", &likePid,"likePid/I");
  tree.Branch("spr", &spr,"spr/D");
  tree.Branch("trr", &trr,"trr/D");
  tree.Branch("nph",&nph,"nph/D");
  tree.Branch("nph_err",&nph_err,"nph_err/D");
  tree.Branch("cangle",&cangle,"cangle/D");
  tree.Branch("likelihood",&likelihood,"par3/D");
  tree.Branch("separation",&separation,"separation/D");
  tree.Branch("par5",&par5,"par5/D");
  tree.Branch("par6",&par6,"par6/D");
  tree.Branch("test1",&test1,"test1/D");
  tree.Branch("test2",&test2,"test2/D");
  tree.Branch("test3",&test3,"test3/D");
  tree.Branch("nnratio",&nnratio,"nnratio/D");
  tree.Branch("nnratio_p",&nnratio_p,"nnratio_p/D");
  tree.Branch("nnratio_pi",&nnratio_pi,"nnratio_pi/D");
  tree.Branch("theta",&theta,"theta/D");
  tree.Branch("beamx",&beamx,"beamx/D");
  tree.Branch("beamz",&beamz,"beamz/D");
  tree.Branch("phi",&phi,"phi/D");
  tree.Branch("timeRes",&timeRes,"timeRes/D");
  
  test1 = PrtManager::Instance()->GetTest1();
  test2 = PrtManager::Instance()->GetTest2();
  test3 = PrtManager::Instance()->GetTest3();
  beamx = PrtManager::Instance()->GetBeamX();
  beamz = PrtManager::Instance()->GetBeamZ();
  par5 = PrtManager::Instance()->GetPrismStepX();
  par6 = PrtManager::Instance()->GetPrismStepY();
  timeRes = PrtManager::Instance()->GetTimeRes();
  fMethod = PrtManager::Instance()->GetRunType();

  Int_t nEvents = fChain->GetEntries();
  if(end==0) end = nEvents;
  
  std::cout<<"Run started for ["<<start<<","<<end <<"]"<<std::endl;
  Int_t nsHits(0),nsEvents(0),studyId(0), nHits(0), ninfit(1);

  if(start<0) {
    ninfit=abs(start);
    start=0;
  }
    
  for (Int_t ievent=start; ievent<start+nEvents && (events[2]<end || events[4]<end); ievent++){

    if(ievent%1000==0) std::cout<<"Event # "<< ievent << " has "<< nHits <<" hits "<< events[2]<<" "<<events[4]<<std::endl;
    Int_t nhhits(0);
    fChain->GetEntry(ievent);
    nHits = fEvent->GetHitSize();
    if(fEvent->GetType()==1) posz = radiatorL/2.-fEvent->GetPosition().Z();
    else posz = fEvent->GetPosition().Z()-15; // 15 mm - lens thickness  
    
    if(ievent-start==0){
      tree.SetTitle(fEvent->PrintInfo());
      prtangle = fEvent->GetAngle(); // prt_data_info.getAngle();
      phi = fEvent->GetPhi(); //prt_data_info.getPhi(); //
      studyId = fEvent->GetGeometry();      
      mom = fEvent->GetMomentum().Mag();

      std::cout<<"prtangle++  "<<prtangle<< " phi "<<phi<<std::endl;
      
      if(fEvent->GetType()==0){
      	momInBar.RotateY(TMath::Pi()-prtangle*CLHEP::deg+test1);
	momInBar.RotateZ(phi*CLHEP::deg+test2);
      }else{
        momInBar.RotateY(TMath::Pi()-prtangle*CLHEP::deg+test1);
	momInBar.RotateZ(phi*CLHEP::deg+test2);
      }

      if(fVerbose==3){
	cz = momInBar.Unit();
	cz = TVector3(-cz.X(),cz.Y(),cz.Z());
      }
    }
    Double_t momentum = fEvent->GetMomentum().Mag();
    if( fEvent->GetType()==1) momentum /= 1000;
    tofPid = fEvent->GetParticle();
    int pid = prt_get_pid(tofPid);
    if(events[pid]>end) continue;
    
    Double_t angle1(0), angle2(0),sum1(0),sum2(0), sigma(0.010),range(5*sigma),noise(0.4);
    
    fAngleP = acos(sqrt(momentum*momentum+ prt_mass[4]*prt_mass[4])/momentum/1.4738)-0.00; //1.4738 = 370 = 3.35
    fAnglePi= acos(sqrt(momentum*momentum + prt_mass[2]*prt_mass[2])/momentum/1.4738)-0.00; //-0.0014 for 160 25deg

    gF1->SetParameter(0,1);
    gF2->SetParameter(0,1);
    gF1->SetParameter(1,fAngleP);
    gF2->SetParameter(1,fAnglePi);
    gF1->SetParameter(2,sigma);
    gF2->SetParameter(2,sigma);

    
    // if(fMethod==2 && tofPid!=2212) continue;
	
    if(fEvent->GetType()==0){

      Int_t gch, ndirc(0), t2(0), t3h(0), t3v(0),
	str1(0),stl1(0),str2(0),stl2(0);
      Int_t hodo1(0), hodo2(0);
      // if(fabs(fEvent->GetMomentum().Mag()-7)<0.1){
      // 	if( pid==4 && fEvent->GetTest1()<34.4 ) continue;
      // 	if( pid==2 && fEvent->GetTest1()>33.3 ) continue;
      // }
      for(auto h=0; h<nHits; h++) {
      	gch = fEvent->GetHit(h).GetChannel();
	if(gch<prt_maxdircch) ndirc++;

	if(gch==513) t2++;
	if(gch==514) t3h++;
	if(gch==515) t3v++;
	if(gch>=1094 && gch<=1101) hodo1++;
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
      // if(!(t3h && t3v)) continue;
    }
    
    //SearchClusters();
    
    for(Int_t h=0; h<nHits; h++) {
      fHit = fEvent->GetHit(h);
      hitTime = fHit.GetLeadTime();
      if(fEvent->GetType()!=0) hitTime+=fRand.Gaus(0,0.2); // time resol. in case it was not simulated
	
      //======================================== dynamic cuts
      {
	Double_t cut1(7);
	{ //time cuts
	  if(prtangle<=80){
	    if(hitTime<cut1 || hitTime>45) continue;
	    reflected = kTRUE;
	  }else if(prtangle>94){
	    if(hitTime<3 || hitTime>20) continue;
	    reflected = kFALSE;
	  }else{
	    if(hitTime<14)  reflected = kFALSE; //13.5
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

      Int_t pixid=fHit.GetPixelId()-1;
      Int_t mcpid=fHit.GetMcpId();
      Int_t ch = map_mpc[mcpid][pixid];
      Int_t sensorId = 100*mcpid+fHit.GetPixelId();
      if(reflected) lenz = 2*radiatorL - posz;
      else lenz = posz;
      
      if(prt_isBadChannel(ch)) continue;
      Int_t nedge=GetEdge(mcpid, pixid);
      //if(cluster[mcpid][pixid]>8) continue;
      //if(sensorId==1) continue;

      Bool_t isGoodHit(0);      
      Int_t size =fLutNode[sensorId]->Entries();

      for(Int_t i=0; i<size; i++){
	weight = 1; //fLutNode[sensorId]->GetWeight(i);
	dird   = fLutNode[sensorId]->GetEntryCs(i,nedge); // nedge=0
        //dird   = fLutNode[sensorId]->GetEntry(i);
	evtime = fLutNode[sensorId]->GetTime(i);
	Int_t pathid = fLutNode[sensorId]->GetPathId(i);
	Bool_t samepath(false);
	if(pathid==fHit.GetPathInPrizm()) samepath=true;
	//if(fLutNode[sensorId]->GetNRefl(i)!=1 ) continue;
	//if(pathid != 130000 && pathid != 199000) continue;
	//std::cout<<"pathid "<< pathid <<std::endl;
	//if(!samepath) continue;
	
	for(int u=0; u<4; u++){
	  // if((pathid==190000 || pathid==210000) && u == 0) continue; //one from left-right
	  // if((pathid==290000 || pathid==310000) && u == 0) continue; //two from left-right
	  // if((pathid==130000 || pathid==199000) && u == 0) continue; //from up-bottom
	  if(u == 0) dir = dird;
	  if(u == 1) dir.SetXYZ( -dird.X(), dird.Y(), dird.Z());
	  if(u == 2) dir.SetXYZ( dird.X(), -dird.Y(),  dird.Z()); //no need when no divergence in vertical plane
	  if(u == 3) dir.SetXYZ( -dird.X(),-dird.Y(), dird.Z()); //no need when no divergence in vertical plane
	  if(reflected) dir.SetXYZ( dir.X(), dir.Y(), -dir.Z());
	  if(dir.Angle(fnX1) < criticalAngle || dir.Angle(fnY1) < criticalAngle) continue;

	  luttheta = dir.Theta();  
	  if(luttheta > TMath::PiOver2()) luttheta = TMath::Pi()-luttheta;
	  
	  bartime = fabs(lenz/cos(luttheta)/198.);
	  double totaltime = bartime+evtime;
	  //if(fEvent->GetType()==0) totaltime+=0.3;
	  double timediff = totaltime-hitTime;
	  fHist0->Fill(timediff);
	  if(samepath)  fHist0i->Fill(timediff);
	  if(fabs(timediff)>timeRes) continue;
	  
	  fHist3->Fill(fabs(totaltime),hitTime);
	  tangle = momInBar.Angle(dir)-0.002;
	  
	  if(tangle > minChangle && tangle < maxChangle && tangle < 1.85){
	    if(tofPid==211 && fMethod==2) fHistPi->Fill(tangle ,weight);
	    else fHist->Fill(tangle ,weight);
	    
	    if(tofPid==2212) fHistMcp[mcpid]->Fill(tangle ,weight);
	    fHistCh[ch]->Fill(tangle ,weight);
	    
	    if(true && tangle>0.4 && tangle<0.9){
	      sum1 += TMath::Log(gF1->Eval(tangle)+noise);
	      sum2 += TMath::Log(gF2->Eval(tangle)+noise);
	    }
	    
	    // //if(samepath) fHist->Fill(tangle ,weight);
	    if((fRadiator==1 && fabs(tangle-0.815)<0.05) || (fRadiator==2 && fabs(tangle-0.815)<0.2)){
	      isGoodHit=true;
	      fHist2->Fill(totaltime);
	    }
	    if(fVerbose==3){
	      TVector3 rdir = TVector3(-dir.X(),dir.Y(),dir.Z());
	      rdir.RotateUz(cz);	      
	      Double_t lphi = rdir.Phi();
	      Double_t tt =  rdir.Theta();
	      fHist4->Fill(tt*TMath::Sin(lphi),tt*TMath::Cos(lphi));

	      //for cherenckov circle fit
	      gg_gr.SetPoint(gg_i,tt*TMath::Sin(lphi),tt*TMath::Cos(lphi));
	      gg_i++;
	    }
	  }
	}
      }

      fHist1->Fill(hitTime);
      if(isGoodHit){
	nhhits++;
	nsHits++;
	prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
      }
    } 

    fhNph->Fill(nhhits);

    if(tofPid==2212) fhNph_p->Fill(nhhits);
    if(tofPid==211) fhNph_pi->Fill(nhhits);
    
    Double_t sum = sum1-sum2;
    if(sum!=0){
      if(tofPid==2212) hLnDiffP->Fill(sum);
      if(tofPid==211) hLnDiffPi->Fill(sum);
      likelihood=sum;
      events[pid]++;
    }

    // if(fVerbose==1){
    //   prt_canvasAdd("ff",800,400);
    //   gF1->Draw();
    //   gF2->SetLineColor(4);
    //   gF2->Draw("same");
      
    //   prt_waitPrimitive("ff");
    //   prt_canvasDel("ff");
    //   //prt_canvasSave(1,0);
    //   //prt_canvasDel(Form("lh_%d",gg_ind));
    // }
	
    if(fVerbose>0 &&  fMethod==3 && nsEvents%ninfit==0){
      if(nsHits>10){
        // if(tofPid==2212 && sum > 0){
	//   std::cout<<"p  "<<sum1 << "   pi "<<sum2 << "  s "<< sum<<std::endl;
	//   if(fVerbose>0)  if(!FindPeak(cangle,spr, prtangle, tofPid)) continue;
	// }

	FindPeak(cangle,spr, prtangle, tofPid);	
	distPid = FindPdg(momentum,cangle);
	nph = nsHits/(Double_t)ninfit;
	spr = spr*1000;
	trr = spr/sqrt(nph);
	theta = fEvent->GetAngle();
	par3 = fEvent->GetTest1();
	tree.Fill();
      }
      ResetHists();
      nsHits=0;
    }

    //if(++nsEvents>=end) break;
  }

  nnratio = fhNph->GetEntries()/(double)end;
  nnratio_pi = fhNph_pi->GetEntries()/(double)end;
  nnratio_p = fhNph_p->GetEntries()/(double)end;

  std::cout<<"nnratio "<<nnratio<<" "<<end <<"  "<< fhNph->GetEntries()<<std::endl;

  TF1 *ff;
  if(fMethod==2){
    gROOT->SetBatch(1);
    if(fhNph->GetEntries()>20){
      fhNph->Fit("gaus","","MQN",5,120);
      ff = fhNph->GetFunction("gaus");
      nph=ff->GetParameter(1);
      nph_err=ff->GetParError(1);
    }
    //nph = prt_fit(fhNph,40,10,50,1).X();
    gROOT->SetBatch(0);
    FindPeak(cangle,spr, prtangle);
    //nph = nsHits/(Double_t)nsEvents;
    spr = spr*1000;
    trr = spr/sqrt(nph);
    theta = fEvent->GetAngle();
    par3 = fEvent->GetTest1();
    if(fVerbose) std::cout<<Form("SPR=%2.2F N=%2.2f +/- %2.2f",spr,nph,nph_err)<<std::endl;     

    // }else{
    if(fVerbose<2) gROOT->SetBatch(1);
    prt_canvasAdd("r_lhood",800,400);
    prt_normalize(hLnDiffP,hLnDiffPi);
    hLnDiffP->SetLineColor(2);

    Double_t m1,m2,s1,s2; 
    if(hLnDiffP->GetEntries()>10){
      hLnDiffP->Fit("gaus","S");
      ff = hLnDiffP->GetFunction("gaus");
      m1=ff->GetParameter(1);
      s1=ff->GetParameter(2);
    }
    if(hLnDiffPi->GetEntries()>10){
      hLnDiffPi->Fit("gaus","S");
      ff = hLnDiffPi->GetFunction("gaus");
      m2=ff->GetParameter(1);
      s2=ff->GetParameter(2);
    }
    separation = (fabs(m2-m1))/(0.5*(s1+s2));
    std::cout<<"separation "<< separation <<std::endl;

    //gStyle->SetOptFit(0);
    //gStyle->SetOptStat(0);
    
    hLnDiffP->SetName(Form("s_%2.2f",separation));
    hLnDiffP->Draw();
    hLnDiffPi->SetLineColor(4);
    hLnDiffPi->Draw("same");
    if(fVerbose) gROOT->SetBatch(0);
  }
  tree.Fill();
  file.Write();

  if(fVerbose>1) prt_waitPrimitive("r_time","w");
  if(fVerbose>0){
    prt_canvasSave(1,0);
    prt_canvasDel("*");
  }
    
  if(fVerbose) ResetHists(); 
}

Int_t g_num =0;
Bool_t PrtLutReco::FindPeak(Double_t& cangle, Double_t& spr, Double_t a, Int_t tofpdg){
  cangle=0;
  spr=0;
  //  gStyle->SetCanvasPreferGL(kTRUE);
  if(fHist->GetEntries()>20 || fHistPi->GetEntries()>20){
    gROOT->SetBatch(1);
    Int_t nfound = fSpect->Search(fHist,1,"",0.9); //0.6
    if(nfound>0) cangle = fSpect->GetPositionX()[0];
    else cangle =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());
    cangle =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());

    if(cangle>0.85) cangle=0.82;
    fFit->SetParameters(100,cangle,0.010);
    fFit->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");      
    fFit->SetParLimits(0,0.1,1E6);
    fFit->SetParLimits(1,cangle-0.04,cangle+0.04); 
    fFit->SetParLimits(2,0.005,0.016); // width 7-10    
    // fFit->FixParameter(2,0.01); 
    // fFit->FixParameter(3,0); 
    // fFit->FixParameter(4,0);

    Int_t status(0);
    if(fMethod==3) status = fHist->Fit("fgaus","lq","",0.6,1);
    else status =fHist->Fit("fgaus","M","",cangle-0.06,cangle+0.06);    
    Double_t chi = fFit->GetChisquare()/fFit->GetNDF();
    
    // if(fFit->GetParError(1)>0.0035){
    // //   // if(fFit->GetParameter(2)>0.011){
    // //   // if(fabs(chi-1<0.3 ){
    //   spr=0;
    //   cangle=0;
    //   fTest=0;
    //   return false;
    // }else{
    //   fTest=chi;
    // }
    
    cangle = fFit->GetParameter(1);
    spr = fFit->GetParameter(2);
 
    if(fVerbose>1) gROOT->SetBatch(0);

    if(fMethod==2 && fVerbose>0){

      TString nid = "";//Form("_%2.0f",a);
      prt_canvasAdd("r_tangle"+nid,800,400);

      // fFit->SetParLimits(2,0.004,0.008); // width 7-10
      // for(Int_t i=0; i<prt_nmcp; i++){
      // 	prt_canvasAdd(Form("r_tangle_%d",i),800,400);
      // 	fHistMcp[i]->Fit("fgaus","lq","",fAngleP-0.03,fAngleP+0.03);
      // 	std::cout<<"if(mcpid=="<< i<<") tangle += "<<fAngleP-fFit->GetParameter(1)<<";" <<std::endl;	
      // 	fHistMcp[i]->Draw();
      // 	drawTheoryLines();
      // }

      // for(Int_t i=0; i<960; i++){
      // 	prt_canvasAdd(Form("r_tangle_ch_%d",i),800,400);
      // 	fHistCh[i]->Fit("fgaus","lq","",fAngleP-0.03,fAngleP+0.03);
      // 	std::cout<<"if(ch=="<< i<<") tangle += "<<fAngleP-fFit->GetParameter(1)<<";" <<std::endl;	
      // 	fHistCh[i]->Draw();
      // }

      //      TString name = Form("r_tangle_%3.1f",test3);
      fHist->SetTitle(Form("theta %3.1f , TOF PID = %d", a, tofpdg));
      fHist->SetMinimum(0);
      //fHist->Scale(1/fHist->GetMaximum());

      prt_normalize(fHist,fHistPi);
      fHistPi->SetLineColor(4);
      fHist->SetLineColor(1);
      
      fHist->Draw();
      fHistPi->Draw("same");
      // gF1->Draw("same");
      // gF2->Draw("same");
      fHisti->SetLineColor(kRed+2);
      if(fHisti->GetEntries()>5) fHisti->Draw("same");

      drawTheoryLines();
      std::cout<<"here0 "<<std::endl;
      prt_canvasAdd("r_time",800,400);
      prt_normalize(fHist1,fHist2);
      fHist1->SetTitle(Form("theta %3.1f", a));
      fHist1->SetLineColor(2);
      fHist1->Draw();
      fHist2->Draw("same");

      prt_canvasAdd("r_nph"+nid,800,400);
      fhNph->Draw();

      prt_canvasAdd("r_diff"+nid,800,400);
      fHist0->SetTitle(Form("theta %3.1f", a));
      fHist0->Draw();
      fHist0i->SetLineColor(kRed+2);
      if(fHist0i->GetEntries()>5)  fHist0i->Draw("same"); 
       
      // prt_canvasAdd("r_cm"+nid,800,400);
      // fHist3->SetTitle(Form("theta %3.1f", a));
      // fHist3->Draw("colz");
      std::cout<<"here1 "<<std::endl;
      if(false){
	Int_t tmax, max=0;
	for(Int_t m=0; m<prt_nmcp;m++){
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
	for(Int_t m=0; m<prt_nmcp;m++){	  
	  prt_hdigi[m]->Scale(1/(Double_t)max);
	}
      }
 
      prt_drawDigi("m,p,v\n",2018);
      prt_cdigi->SetName("r_hp"+nid);
      prt_canvasAdd(prt_cdigi);
      		 
      if(fVerbose==3){
	TCanvas* c2 = new TCanvas("c2","c2",0,0,800,400);
	c2->Divide(2,1);
	c2->cd(1);
     
	fHist4->SetStats(0);
	fHist4->SetTitle(Form("Calculated from LUT, #theta = %3.1f#circ", a));
	fHist4->Draw("colz");
	Double_t x0(0), y0(0), theta(cangle);
	FitRing(x0,y0,theta);
	TVector3 corr(x0,y0,1-TMath::Sqrt(x0*x0+y0*y0));
	std::cout<<"Tcorr "<< corr.Theta()*1000<< "  Pcorr "<< corr.Phi() <<std::endl;

	TLegend *leg = new TLegend(0.5,0.7,0.85,0.87);
	//      leg->SetFillColor(0);
	//leg->SetFillColorAlpha(0,0.8);
	leg->SetFillStyle(0);
	//leg->SetFillStyle(4000); 
	leg->SetBorderSize(0);
	leg->AddEntry((TObject*)0,Form("Entries %0.0f",fHist4->GetEntries()),"");
	leg->AddEntry((TObject*)0,Form("#Delta#theta_{c} %f [mrad]",corr.Theta()*1000),"");
	leg->AddEntry((TObject*)0,Form("#Delta#varphi_{c} %f [mrad]",corr.Phi()),"");
	leg->Draw();

	TArc *arc = new TArc(x0,y0,theta);
	arc->SetLineColor(kRed);
	arc->SetLineWidth(1);
	arc->SetFillStyle(0);
	arc->Draw();
	gg_i=0;
	gg_gr.Set(0);

	c2->cd(2);
	gStyle->SetOptStat(1110); 
	fHist5->SetTitle(Form("True from MC, #theta = %3.1f#circ", a));
	fHist5->Draw("colz");

	c2->Print(Form("spr/tcorr_%3.1f.png", a));
	c2->Modified();
	c2->Update();
	c2->WaitPrimitive("");
      }
    }
  }

  if(fVerbose<2) gROOT->SetBatch(0);
  
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
  for(Int_t m=0; m<prt_nmcp;m++) prt_hdigi[m]->Reset();
}

TF1 *lFit = new TF1("lgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.6,0.9);
TF1 *lFitPi = new TF1("lgausPi","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.6,0.9);

Double_t PrtLutReco::fillLnDiffPPi(Double_t cangle, Int_t tofPid, Double_t mom){
  if(fHist->GetEntries()>20 ){
    Double_t angle1(0), angle2(0), sigma(0.006),range(0.015);

    // //fHist->Scale(1/fHist->GetMaximum());

    // Double_t d1,d2, sum1(0),sum2(0);
    // Int_t sbin = fHist->FindBin(fAngleP-range);
    // Int_t ebin = fHist->FindBin(fAngleP+range); 
    // // fHist->GetXaxis()->GetNbins()
    // for(Int_t i=sbin; i< ebin; i++){
    //   if(fHist->GetBinContent(i) < 0.01 ) continue;
    //   d1 = gF1->Eval(fHist->GetBinCenter(i))- fHist->GetBinContent(i);
    //   d2 = gF1->Eval(fHist->GetBinCenter(i))- fHist->GetBinContent(i);

    //   std::cout<<"f1 "<< gF1->Eval(fHist->GetBinCenter(i)) << "   f2 "<<gF2->Eval(fHist->GetBinCenter(i)) << "    v "<< fHist->GetBinContent(i) <<std::endl;
    
  //   // if(d1>0) sum1+=TMath::Log(d1);
  //   // if(d2>0) sum2+=TMath::Log(d2);
  //   sum1+=TMath::Log(fabs(d1));
  //   sum2+=TMath::Log(fabs(d2));

  // }
  // Double_t amin(sum1),amin2(sum2); 
  
  // lFit->SetRange(fAngleP-range,fAngleP+range);
  // lFit->FixParameter(0,fFit->GetParameter(0));
  // lFit->FixParameter(1,fAngleP);
  // if(fFit->GetParameter(2)>sigma) sigma=fFit->GetParameter(2);
  // lFit->FixParameter(2,sigma);
  // lFit->FixParameter(3,fFit->GetParameter(3));
  // lFit->FixParameter(4,fFit->GetParameter(4));
    
  lFit->SetRange(fAngleP-range,fAnglePi+range);
  lFit->FixParameter(0,fHist->GetMaximum()-0.5);
  lFit->FixParameter(1,fAngleP);
  lFit->FixParameter(2,0.01);
  lFit->FixParameter(3,0);
  lFit->FixParameter(4,0.5);
 

  fHist->Fit("lgaus","lq","",fAngleP-range,fAnglePi+range);
  Double_t amin,amin2,edm,errdef;
  Int_t nvpar,nparx;
  TVirtualFitter *fitter = TVirtualFitter::Fitter(fHist);
  fitter->GetStats(amin,edm,errdef,nvpar,nparx);
  
  // lFitPi->SetRange(fAnglePi-range,fAnglePi+range);
  // lFitPi->SetLineColor(4);
  // lFitPi->FixParameter(0,fFit->GetParameter(0));
  // lFitPi->FixParameter(1,fAnglePi);
  // lFitPi->FixParameter(2,sigma);
  // lFitPi->FixParameter(3,fFit->GetParameter(3));
  // lFitPi->FixParameter(4,fFit->GetParameter(4));

  lFitPi->SetRange(fAngleP-range,fAnglePi+range);
  lFitPi->SetLineColor(4);
  lFitPi->FixParameter(0,fHist->GetMaximum()-0.5);
  lFitPi->FixParameter(1,fAnglePi);
  lFitPi->FixParameter(2,0.01);
  lFitPi->FixParameter(3,0);
  lFitPi->FixParameter(4,0.5);

  fHist->Fit("lgausPi","lq","",fAngleP-range,fAnglePi+range);
  fitter = TVirtualFitter::Fitter(fHist);
  fitter->GetStats(amin2,edm,errdef,nvpar,nparx);
  
  if(fVerbose) printf("tofPid %04d | %1.4f (%1.4f/%1.4f) likelihood is %1.2f/%1.2f \n",tofPid,cangle,fAngleP,fAnglePi, amin, amin2);
  gg_ind++;

  if(fVerbose==1){
    prt_canvasAdd("ff",800,400);
    //prt_canvasAdd(Form("lh_%d",gg_ind),800,400);
    fHist->SetTitle(Form("%d",tofPid));
    fHist->Draw();
    lFit->SetLineColor(2);
    lFit->Draw("same");
    // gF1->Draw("same");
    // gF2->SetLineColor(4);
    // gF2->Draw("same");
     
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

Double_t PrtLutReco::fillLnDiffPPi2(Double_t cangle, Int_t tofPid, Double_t mom){
 if(fHist->GetEntries()>20 ){
  Double_t angle1(0), angle2(0), sigma(0.006),range(0.03);

  Double_t d1,d2, sum1(0),sum2(0);
  Int_t sbin = fHist->FindBin(fAnglePi-range);
  Int_t ebin = fHist->FindBin(fAngleP+range); 
  for(Int_t i=sbin; i< ebin; i++){
    if(fHist->GetBinContent(i)<1 ) continue;
    d1 = 10*fabs(fHist->GetBinContent(i) *(fAngleP  - fHist->GetBinCenter(i)));
    d2 = 10*fabs(fHist->GetBinContent(i) *(fAnglePi - fHist->GetBinCenter(i)));
    if(d1>0 && d2>0){
      std::cout<<"d1  "<<d1 << "   d2    "<< d2 <<std::endl;
      sum1+=TMath::Log(d1);
      sum2+=TMath::Log(d2);
    }
  }
  
  if(fVerbose) printf("tofPid %04d | %1.4f (%1.4f/%1.4f) likelihood is %1.2f/%1.2f \n",tofPid,cangle,fAngleP,fAnglePi, sum1, sum2);
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

void circleFcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
  Int_t np = gg_gr.GetN();
  f = 0;
  Double_t *x = gg_gr.GetX();
  Double_t *y = gg_gr.GetY();
  for (Int_t i=0;i<np;i++) {
    Double_t u = x[i] + par[0];
    Double_t v = y[i] + par[1];
    Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
    f += dr*dr;  
  }
  std::cout<<"fcn  "<< f<<std::endl;
  
}

void circleFcn2(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
  Int_t np = gg_gr.GetN();
  f = 0;
  Double_t *x = gg_gr.GetX();
  Double_t *y = gg_gr.GetY();
  for (Int_t i=0;i<np;i++){
    Double_t u = x[i] + par[0];
    Double_t v = y[i] + par[1];
    Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
    if(dr>0.07) f += dr*dr; 
    else f += fabs(dr);
  }
}

void PrtLutReco::FitRing(Double_t& x0, Double_t& y0, Double_t& theta){
  TGraph ff_gr;
  Int_t ff_i(0);
  Int_t np = gg_gr.GetN();
  Double_t *x = gg_gr.GetX();
  Double_t *y = gg_gr.GetY();
  for (Int_t i=0;i<np;i++) {
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
  Double_t arglist[1] = {0};
  fitter->ExecuteCommand("MINIMIZE", arglist, 0);

  // fitter->SetFCN(circleFcn2);
  // fitter->ExecuteCommand("MINIMIZE", arglist, 0);

  x0 = fitter->GetParameter(0);
  y0 = fitter->GetParameter(1);
  theta = fitter->GetParameter(2);
}

Int_t PrtLutReco::FindPdg(Double_t mom, Double_t cangle){
  cangle =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());
  
  Int_t pdg[]={211,2212};
  Double_t mass[] = {0.139570,0.9382723};
  Double_t tdiff, diff=100;
  Int_t minid=0;
  for(Int_t i=0; i<2; i++){
    tdiff = fabs(cangle - acos(sqrt(mom*mom + mass[i]*mass[i])/mom/1.46907)); //1.46907 - fused silica
    if(tdiff<diff){
      diff = tdiff;
      minid = i;
    }
  }
  return pdg[minid]; 
}

Int_t PrtLutReco::GetEdge(Int_t mcpid, Int_t pixid){
  // Int_t x(0),y(0), piid(pixid) , nedge(0); //new
  // for(Int_t h=0; h<fEvent->GetHitSize(); h++) {
  // 	Int_t pid=fEvent->GetHit(h).GetPixelId();
  // 	Int_t mid=fEvent->GetHit(h).GetMcpId();
  // 	Double_t tdif=fabs(hitTime-fEvent->GetHit(h).GetLeadTime());
  // 	if(mid!=mcpid || pid==piid || tdif>0.3) continue;
  // 	if(pid==piid-1 && piid%8!=0) y-=1;
  // 	if(pid==piid+1 && piid%8!=7) y+=1;

  // 	if(pid==piid+8 && piid<57) x-=1;
  // 	if(pid==piid-8 && piid>8)  x+=1;
  // }

  Int_t x(0),y(0),piid(pixid+1),nedge(0); //old
  for(Int_t h=0; h<fEvent->GetHitSize(); h++) {
    Int_t pid=fEvent->GetHit(h).GetPixelId();
    Int_t mid=fEvent->GetHit(h).GetMcpId();
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

Int_t getneighbours(Int_t m, Int_t p){
  for(Int_t i=0; i<65; i++) if(p==lneighbours[i]) return -1;
  lneighbours[lsize]=p;
  lsize++;
  for(Int_t t=0; t<65; t++){
    if(mcpdata[m][t]){
      for(Int_t i=0; i<65; i++) if(t==lneighbours[i]) continue;
      if((t==p-1 && p%8!=0) || (t==p+1 && p%8!=7) ||
	 (t==p+8 && p<57) || (t==p-8 && p>8)) getneighbours(m,t);
    }
  }
  return lsize;
}

void getclusters(){
  for(Int_t m=0; m<prt_nmcp; m++){
    for(Int_t p=0; p<65; p++){
      if(mcpdata[m][p])  cluster[m][p] = getneighbours(m,p);
      lsize=0;
      for(Int_t i=0; i<65; i++) lneighbours[i]=0;
    }
  }
}
  
void PrtLutReco::SearchClusters(){

  for(Int_t j=0; j<prt_nmcp; j++){
    for(Int_t i=0; i<65; i++){
      mcpdata[j][i]=0;
      cluster[j][i]=0;
    }
  }
  
  for(Int_t h=0; h<fEvent->GetHitSize(); h++) {
    Int_t mid=fEvent->GetHit(h).GetMcpId();
    Int_t pid=fEvent->GetHit(h).GetPixelId()-1;
    mcpdata[mid][pid]=1;
  }
  getclusters();
}

void PrtLutReco::drawTheoryLines(){
  gPad->Update();
  TLine *line = new TLine(0,0,0,1000);
  line->SetX1(fAngleP);
  line->SetX2(fAngleP);
  line->SetY1(gPad->GetUymin());
  line->SetY2(gPad->GetUymax());
  line->SetLineColor(kRed);
  line->Draw();

  TLine *line1 = new TLine(0,0,0,1000);
  line1->SetX1(fAnglePi);
  line1->SetX2(fAnglePi);
  line1->SetY1(gPad->GetUymin());
  line1->SetY2(gPad->GetUymax());
  line1->SetLineColor(kBlue);
  line1->Draw();
}
