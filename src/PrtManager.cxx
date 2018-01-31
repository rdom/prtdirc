#include "TInterpreter.h"
#include "G4String.hh"
#include "G4SystemOfUnits.hh"

#include "PrtManager.h"
#include "PrtHit.h"
#include "PrtLutNode.h"

PrtManager * PrtManager::fInstance= NULL;

PrtManager::PrtManager(G4String outfile, G4int runtype){
  TString filename = outfile.c_str();
  fOutName = filename; 
  fOutName = fOutName.Remove(fOutName.Last('.'));
  fRunType = runtype;
  
  if(fRunType!=2 && fRunType!=3) fRootFile = new TFile(filename,"RECREATE");

  fEvent = new PrtEvent();
  if(fRunType==0 || fRunType==6){
    fTree = new TTree("data","Prototype hits tree");
    fTree->Branch("PrtEvent", "PrtEvent", &fEvent, 64000, 2);
  }

  if(fRunType==1 || fRunType==5 || fRunType==11){
    fLut = new TClonesArray("PrtLutNode");
    fTree = new TTree("prtlut","Look-up table for the geometrical reconstruction.");
    fTree->Branch("LUT",&fLut,256000,2); 
    Int_t Nnodes = 5000;
    
    TClonesArray &fLuta = *fLut; 
    for (Long64_t n=0; n<Nnodes; n++) {
      new((fLuta)[n]) PrtLutNode(n);
    }    
  }

  // if(fRunType==2){
  //   fTree = new TTree("recodata","Reconstructed info for the prototype");
  //   fTrackInfoArray = new TClonesArray("PrtTrackInfo");
  //   fTree->Branch("PrtTrackInfo",&fTrackInfoArray,256000,2); 
  // }
  
  // fHist = new TH1F("id", "name", 100, 0., 100);

  fPhysList = 0;
  fParticle = 0;
  fMomentum = TVector3(0,0,0);
  fGeometry = 3;
  fAngle = 0;
  fPhi = 0;
  fRadiator=1;
  fRadiatorL=0;
  fRadiatorW=0;
  fRadiatorH=0;
  fShift = 150;
  fTest1 = 0;
  fTest2 = 0;
  fTest3 = 0;
  fLens = 0;
  fMcpLayout = 2014;
  fBeamDimension = 0;

  fPrismStepX=0;
  fPrismStepY=0;
  fRStepX=0;
  fRStepY=0;
  fLStepX=0;
  fLStepY=0;
  fBeamX=0;
  fBeamZ=-1;
  fTime=0;
  fTimeRes=0.2;
  fInfo="";
  fMixPiP=false;

  fnX1 = TVector3(1,0,0);   
  fnY1 = TVector3( 0,1,0);
  fCriticalAngle = asin(1.00028/1.47125);
  
  std::cout<<"PrtManager has been successfully initialized. " <<std::endl;
}

PrtManager* PrtManager::Instance(G4String outfile, G4int runtype){
  if ( !fInstance){
    std::cout<<"Info in (PrtManager::Instance): Making a new instance. "<<std::endl;
    fInstance = new PrtManager(outfile, runtype);
  }
  return fInstance;
}

void PrtManager::AddEvent(PrtEvent event){
  if(fRunType==0 || fRunType==6){
    fEvent = new PrtEvent(event);
    fEvent->SetType(1);
    fEvent->SetPhysList(fPhysList);
    //    fEvent->SetAngle((180*deg-fAngle)/deg);
    fEvent->SetAngle(fAngle);
    fEvent->SetPhi(fPhi);
    
    //fEvent->SetRadiatorL(fRadiatorL);
    //fEvent->SetRadiatorW(fRadiatorW);
    //fEvent->SetRadiatorH(fRadiatorH);
    fEvent->SetParticle(fParticle);
    fEvent->SetMomentum(fMomentum);
    fEvent->SetGeometry(fGeometry);
    fEvent->SetLens(fLens);
    fEvent->SetTest1(fTest1);
    fEvent->SetTest2(fTest2);
    fEvent->SetPrismStepX(fPrismStepX);
    fEvent->SetPrismStepY(fPrismStepY);
    fEvent->SetBeamX(fBeamX);
    fEvent->SetBeamZ(fBeamZ);
    fEvent->SetTimeRes(fTimeRes);
    fEvent->SetTime(fTime);
    fEvent->SetInfo(fInfo);
  }
}


void PrtManager::AddHit(PrtHit hit){
  if(fRunType==0 || fRunType==6){
    if ( fEvent ){
      fEvent->AddHit(hit);
    }else{
      std::cout<<"Event does not exist. Create it first. "<<std::endl;
    }
  }
  if(fRunType==1 || fRunType==5 || fRunType==11){
    if(fMomentum.Angle(fnX1) > fCriticalAngle && fMomentum.Angle(fnY1) > fCriticalAngle){
      Int_t id = 100*hit.GetMcpId() + hit.GetPixelId();
      Double_t time = hit.GetLeadTime();
      if(fRunType==11) time -= fTime;      
      ((PrtLutNode*)(fLut->At(id)))->
	AddEntry(id, fMomentum, hit.GetPathInPrizm(),
		 hit.GetNreflectionsInPrizm(),
		 time,hit.GetLocalPos(),hit.GetDigiPos());
    }
  }
}

void PrtManager::AddTrackInfo(PrtTrackInfo trackinfo){
  // new ((*fTrackInfoArray)[fTrackInfoArray->GetEntriesFast()]) PrtTrackInfo(trackinfo);
}



void PrtManager::Fill(){
  if(fRunType==0 || fRunType==6){
    fTree->Fill();
    fEvent->Clear();
  }
  // if(fRunType==2){
  //   fTree->Fill();
  //   fTrackInfoArray->Clear();
  // }
}

void PrtManager::FillLut(){
  if(fRunType==1 || fRunType==5 || fRunType==11) fTree->Fill();
}
