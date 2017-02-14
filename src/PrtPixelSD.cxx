#include "PrtPixelSD.h"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include <TVector3.h>

#include "PrtEvent.h"
#include "PrtPrizmHit.h"

#include "PrtRunAction.h"
#include "PrtManager.h"

PrtPixelSD::PrtPixelSD( const G4String& name, 
			const G4String& hitsCollectionName,
			G4int nofCells)
  : G4VSensitiveDetector(name){
  collectionName.insert(hitsCollectionName);
}

PrtPixelSD::~PrtPixelSD(){ 
}

void PrtPixelSD::Initialize(G4HCofThisEvent* hce){


  Double_t qe_m[15][64]={{0.63,0.78,0.73,0.69,0.67,0.67,0.69,0.60,0.73,0.82,0.78,0.76,0.75,0.77,0.82,0.78,0.79,0.85,0.78,0.76,0.75,0.77,0.82,0.79,0.81,0.85,0.78,0.76,0.75,0.77,0.83,0.80,0.81,0.84,0.77,0.74,0.74,0.77,0.83,0.79,0.80,0.81,0.76,0.73,0.74,0.76,0.82,0.78,0.80,0.79,0.75,0.74,0.76,0.78,0.83,0.79,0.60,0.76,0.73,0.72,0.73,0.76,0.79,0.70},
		{0.61,0.69,0.70,0.70,0.70,0.69,0.69,0.61,0.77,0.82,0.83,0.83,0.83,0.82,0.81,0.76,0.77,0.83,0.84,0.84,0.84,0.84,0.82,0.77,0.78,0.84,0.85,0.85,0.85,0.84,0.83,0.78,0.78,0.84,0.85,0.85,0.85,0.84,0.83,0.79,0.77,0.83,0.85,0.85,0.85,0.85,0.83,0.80,0.76,0.83,0.84,0.85,0.85,0.84,0.83,0.81,0.71,0.82,0.84,0.84,0.84,0.84,0.83,0.80},
		{0.66,0.66,0.66,0.67,0.66,0.66,0.67,0.62,0.72,0.71,0.72,0.77,0.77,0.77,0.79,0.77,0.72,0.71,0.71,0.74,0.77,0.77,0.79,0.76,0.71,0.71,0.70,0.70,0.75,0.76,0.79,0.75,0.71,0.70,0.70,0.69,0.70,0.73,0.77,0.74,0.71,0.70,0.70,0.69,0.69,0.69,0.72,0.72,0.71,0.71,0.70,0.69,0.69,0.68,0.69,0.68,0.68,0.71,0.70,0.69,0.68,0.68,0.69,0.65},
		{0.79,0.86,0.85,0.84,0.84,0.83,0.82,0.77,0.85,0.90,0.90,0.89,0.89,0.89,0.88,0.84,0.85,0.90,0.90,0.90,0.90,0.89,0.88,0.84,0.85,0.90,0.90,0.90,0.89,0.89,0.88,0.84,0.86,0.91,0.91,0.90,0.89,0.89,0.88,0.83,0.86,0.91,0.90,0.90,0.89,0.89,0.87,0.84,0.86,0.90,0.90,0.89,0.89,0.88,0.87,0.83,0.83,0.86,0.86,0.86,0.86,0.85,0.85,0.82},
		{0.46,0.52,0.52,0.52,0.52,0.52,0.52,0.48,0.75,0.82,0.83,0.83,0.83,0.83,0.82,0.79,0.76,0.83,0.84,0.84,0.84,0.84,0.83,0.80,0.75,0.83,0.84,0.84,0.84,0.84,0.83,0.80,0.76,0.83,0.84,0.84,0.84,0.84,0.83,0.80,0.74,0.83,0.84,0.84,0.84,0.84,0.83,0.79,0.75,0.83,0.84,0.84,0.84,0.84,0.83,0.80,0.74,0.82,0.83,0.83,0.83,0.83,0.82,0.79},
		{0.66,0.70,0.71,0.70,0.71,0.70,0.70,0.56,0.78,0.81,0.81,0.81,0.82,0.82,0.82,0.68,0.78,0.81,0.82,0.81,0.82,0.82,0.82,0.69,0.77,0.81,0.82,0.82,0.82,0.82,0.82,0.69,0.77,0.82,0.82,0.82,0.82,0.82,0.82,0.68,0.77,0.82,0.82,0.82,0.83,0.83,0.82,0.69,0.77,0.82,0.82,0.82,0.83,0.83,0.82,0.71,0.71,0.77,0.75,0.76,0.76,0.76,0.76,0.64},
		{0.50,0.63,0.64,0.64,0.64,0.65,0.65,0.56,0.72,0.90,0.92,0.91,0.89,0.90,0.89,0.74,0.72,0.92,0.93,0.93,0.93,0.93,0.91,0.76,0.71,0.92,0.94,0.94,0.94,0.93,0.92,0.77,0.70,0.92,0.94,0.94,0.94,0.93,0.92,0.77,0.69,0.92,0.94,0.94,0.93,0.93,0.91,0.77,0.66,0.89,0.92,0.92,0.92,0.91,0.89,0.76,0.46,0.74,0.75,0.76,0.77,0.77,0.75,0.63},
		{0.62,0.68,0.68,0.71,0.71,0.70,0.69,0.52,0.70,0.76,0.77,0.80,0.81,0.80,0.79,0.65,0.71,0.76,0.77,0.81,0.81,0.81,0.80,0.67,0.71,0.76,0.76,0.80,0.82,0.81,0.80,0.69,0.70,0.75,0.76,0.78,0.81,0.81,0.80,0.69,0.70,0.75,0.76,0.76,0.80,0.81,0.80,0.70,0.71,0.75,0.75,0.75,0.76,0.79,0.79,0.70,0.57,0.62,0.62,0.62,0.62,0.62,0.63,0.56},
		{0.57,0.69,0.70,0.72,0.73,0.73,0.72,0.66,0.76,0.89,0.90,0.91,0.91,0.91,0.89,0.84,0.77,0.90,0.92,0.93,0.93,0.92,0.91,0.84,0.77,0.91,0.92,0.93,0.93,0.92,0.91,0.85,0.77,0.91,0.92,0.93,0.93,0.92,0.91,0.85,0.76,0.90,0.92,0.93,0.93,0.92,0.91,0.84,0.74,0.89,0.91,0.91,0.92,0.91,0.89,0.83,0.64,0.78,0.79,0.82,0.83,0.83,0.82,0.76},
		{0.74,0.77,0.79,0.82,0.78,0.65,0.78,0.74,0.79,0.81,0.80,0.86,0.82,0.75,0.82,0.70,0.80,0.82,0.81,0.83,0.83,0.69,0.82,0.80,0.81,0.82,0.81,0.81,0.79,0.78,0.87,0.84,0.79,0.82,0.81,0.80,0.80,0.79,0.90,0.73,0.79,0.80,0.81,0.80,0.78,0.82,0.97,0.80,0.82,0.81,0.81,0.80,0.80,0.80,0.93,0.83,0.81,0.83,0.83,0.81,0.79,0.80,0.93,0.81},
		{0.75,0.91,0.90,0.90,0.89,0.87,0.86,0.75,0.87,0.98,0.99,0.99,0.99,0.99,0.98,0.90,0.89,0.99,1.00,1.00,1.00,1.00,0.99,0.92,0.91,0.99,1.00,1.01,1.01,1.01,1.00,0.92,0.92,0.99,1.00,1.01,1.01,1.01,0.99,0.92,0.91,0.98,1.00,1.00,1.00,1.00,0.99,0.91,0.91,0.97,0.98,0.99,0.99,0.99,0.98,0.89,0.65,0.85,0.88,0.91,0.92,0.94,0.94,0.82},
		{0.77,0.85,0.81,0.84,0.83,0.83,0.83,0.72,0.88,0.93,0.92,0.91,0.90,0.90,0.89,0.83,0.87,0.94,0.93,0.91,0.91,0.90,0.90,0.85,0.88,0.94,0.93,0.91,0.91,0.91,0.90,0.86,0.87,0.94,0.93,0.92,0.91,0.91,0.90,0.86,0.87,0.94,0.92,0.91,0.91,0.90,0.90,0.87,0.85,0.93,0.92,0.90,0.90,0.89,0.89,0.86,0.75,0.86,0.84,0.82,0.81,0.80,0.82,0.76},
		{0.67,0.77,0.78,0.78,0.78,0.78,0.78,0.70,0.80,0.92,0.93,0.93,0.93,0.93,0.92,0.83,0.80,0.93,0.94,0.94,0.94,0.94,0.93,0.84,0.79,0.94,0.95,0.95,0.95,0.94,0.93,0.85,0.77,0.94,0.95,0.95,0.95,0.95,0.93,0.85,0.76,0.94,0.95,0.95,0.95,0.94,0.93,0.85,0.75,0.93,0.94,0.94,0.94,0.93,0.92,0.85,0.67,0.84,0.84,0.84,0.84,0.84,0.83,0.75},
		{0.68,0.80,0.79,0.79,0.77,0.76,0.75,0.62,0.76,0.88,0.89,0.89,0.89,0.89,0.88,0.74,0.77,0.89,0.90,0.90,0.90,0.90,0.89,0.75,0.78,0.89,0.90,0.90,0.91,0.90,0.89,0.77,0.78,0.89,0.90,0.90,0.90,0.90,0.89,0.77,0.78,0.88,0.89,0.90,0.90,0.89,0.88,0.78,0.78,0.87,0.88,0.89,0.89,0.88,0.87,0.78,0.64,0.75,0.75,0.76,0.76,0.76,0.75,0.65},
		{0.66,0.80,0.80,0.80,0.80,0.79,0.79,0.66,0.76,0.88,0.87,0.88,0.89,0.89,0.90,0.79,0.76,0.88,0.87,0.88,0.89,0.90,0.90,0.79,0.76,0.87,0.87,0.87,0.89,0.90,0.90,0.79,0.77,0.86,0.86,0.87,0.88,0.89,0.90,0.79,0.75,0.86,0.85,0.86,0.87,0.88,0.90,0.79,0.76,0.85,0.84,0.84,0.86,0.87,0.89,0.78,0.54,0.64,0.64,0.70,0.74,0.76,0.78,0.63}
  };
  
  
  for(Int_t m=0; m<15; m++){
    for(Int_t i=0; i<64; i++){
      fQe_space[m][i]=qe_m[m][i];
    }
  }
  
  
 // TTree *gTree = new TTree("Prt","Prototype hits tree");
 // Event *fHit = 0;
 // gTree->Branch("event", "Event", &event, 64000, 0);

  // // Create hits collection
  // fHitsCollection 
  //   = new B4cCalorHitsCollection(SensitiveDetectorName, collectionName[0]); 

  // // Add this collection in hce
  // G4int hcID 
  //   = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  // hce->AddHitsCollection( hcID, fHitsCollection ); 

  // // Create hits
  // // fNofCells for cells + one more for total sums 
  // for (G4int i=0; i<fNofCells+1; i++ ) {
  //   fHitsCollection->insert(new B4cCalorHit());
  // }
 
  //PrtManager::Instance()->AddEvent(PrtEvent());
}

G4bool PrtPixelSD::ProcessHits(G4Step* step, G4TouchableHistory* hist){  
  // // energy deposit
  // G4double edep = step->GetTotalEnergyDeposit();
  
  // // step length
  // G4double stepLength = 0.;
  // if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
  //   stepLength = step->GetStepLength();
  // }

  // if ( edep==0. && stepLength == 0. ) return false;     

  if(step == 0) return false;
 
  //G4ThreeVector translation = hist->GetTranslation();
  //G4ThreeVector localpos = step->GetPreStepPoint()->GetPhysicalVolume()->GetObjectTranslation();
  G4TouchableHistory* touchable = (G4TouchableHistory*)(step->GetPostStepPoint()->GetTouchable());

  // Get cell id 
  G4int layerNumber = touchable->GetReplicaNumber(0);
  //G4cout<< " PixelId = "<<layerNumber << G4endl;
  G4Track* track = step->GetTrack();
  const G4DynamicParticle* dynParticle = track->GetDynamicParticle();
  G4ParticleDefinition* particle = dynParticle->GetDefinition();  
  G4String ParticleName = particle->GetParticleName();
   
  G4ThreeVector globalpos = step->GetPostStepPoint()->GetPosition();
  G4ThreeVector localpos = touchable->GetHistory()->GetTopTransform().TransformPoint(globalpos);
  G4ThreeVector translation = touchable->GetHistory()->GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0,0,0));
  G4ThreeVector inPrismpos = touchable->GetHistory()->GetTransform( 1 ).TransformPoint(globalpos);
  G4ThreeVector g4mom = track->GetMomentum(); //track->GetVertexMomentumDirection(); //
  G4ThreeVector g4pos = track->GetVertexPosition();
 
  TVector3 globalPos(inPrismpos.x(),inPrismpos.y(),inPrismpos.z());
  
  if(PrtManager::Instance()->GetRunType() == 6){ //focal plane scan
    globalPos = TVector3(globalpos.x(),globalpos.y(),globalpos.z());
  }
  
  TVector3 localPos(localpos.x(),localpos.y(),localpos.z());
  translation=touchable->GetHistory()->GetTransform( 1 ).TransformPoint(translation);
  TVector3 digiPos(translation.x(),translation.y(),translation.z());
  TVector3 momentum(g4mom.x(),g4mom.y(),g4mom.z());
  G4ThreeVector lp = touchable->GetHistory()->GetTransform(1).TransformPoint(g4pos); //pos in wDirc
  TVector3 position(lp.x(),lp.y(),lp.z());
  
  // information form prizm
  G4SDManager* fSDM = G4SDManager::GetSDMpointer();
  G4RunManager* fRM = G4RunManager::GetRunManager();
  G4int collectionID = fSDM->GetCollectionID("PrizmHitsCollection");
  const G4Event* currentEvent = fRM->GetCurrentEvent();
  G4HCofThisEvent* HCofEvent = currentEvent->GetHCofThisEvent();
  PrtPrizmHitsCollection* prizmCol = (PrtPrizmHitsCollection*)(HCofEvent->GetHC(collectionID));
 
  Double_t pathId = 0;
  Int_t refl=0;
  for (G4int i=0;i<prizmCol->entries();i++){
    PrtPrizmHit* phit = (*prizmCol)[i];
    if(phit->GetTrackID()==track->GetTrackID()) {
      refl++;
      pathId += phit->GetNormalId()*1000*refl;
    }
  }

  //std::cout<<"Number of reflections: "<<refl <<std::endl;

  PrtHit hit;
  Int_t mcpid=touchable->GetReplicaNumber(1);
  Int_t pixid = touchable->GetReplicaNumber(0);
  hit.SetMcpId(mcpid);
  hit.SetPixelId(pixid);
  hit.SetGlobalPos(globalPos);
  hit.SetLocalPos(localPos);
  hit.SetDigiPos(digiPos);
  hit.SetPosition(position);
  hit.SetMomentum(momentum);
  if(PrtManager::Instance()->GetRunType()==6){
    G4ThreeVector mominend = step->GetPostStepPoint()->GetMomentum();
    TVector3 mominendv(mominend.x(),mominend.y(),mominend.z());
    hit.SetMomentum(mominendv);
  }
  hit.SetParticleId(track->GetTrackID());
  hit.SetParentParticleId(track->GetParentID());
  hit.SetNreflectionsInPrizm(refl-1);
  hit.SetPathInPrizm(pathId);
  hit.SetCherenkovMC(PrtManager::Instance()->GetCurrentCherenkov());
  // time since track created

  G4double time = step->GetPreStepPoint()->GetLocalTime();
  if(PrtManager::Instance()->GetRunType()==0) time = G4RandGauss::shoot(time,PrtManager::Instance()->GetTimeRes()); //resolution [ns]
  hit.SetLeadTime(time);
  Double_t wavelength = 1.2398/(track->GetMomentum().mag()*1E6)*1000;
  hit.SetTotTime(wavelength); //set photon wavelength
  // time since event created
  // step->GetPreStepPoint()->GetGlobalTime()*1000

  return true;
  Bool_t charge_sharing(true);
  if(PrtManager::Instance()->GetRunType()==0 && PrtManager::Instance()->GetMcpLayout()>=2015 && charge_sharing){
    if(fQe_space[mcpid][pixid]>G4UniformRand()) PrtManager::Instance()->AddHit(hit);
    else charge_sharing=false;
  }else{
    PrtManager::Instance()->AddHit(hit);
  }

  if(PrtManager::Instance()->GetRunType()==0 && PrtManager::Instance()->GetMcpLayout()>=2015 && charge_sharing){
    //charge sharing for 8x8 MCP
    Double_t pixdim(53/16.),chargesig(1),threshold(0.3);
    Double_t x(localPos.x()), y(localPos.y());
    Int_t p(pixid);
    Bool_t ok(false);
    Double_t expd = exp(-(pixdim-fabs(x))/chargesig);
    
    if(x<0 && pixid%8!=1 && expd>G4UniformRand() && expd<threshold){ok=true; p-=1;}
    if(x>0 && pixid%8!=0 && expd>G4UniformRand() && expd<threshold){ok=true; p+=1;}
    expd = exp(-(pixdim-fabs(y))/chargesig);
    if(y<0 && pixid>8    && expd>G4UniformRand() && expd<threshold){ok=true; p-=8;}
    if(y>0 && pixid<57   && expd>G4UniformRand() && expd<threshold){ok=true; p+=8;}
 
    if(ok) {
      hit.SetPixelId(p);
      PrtManager::Instance()->AddHit(hit);
    }
    
  }
      
  return true;
}

void PrtPixelSD::EndOfEvent(G4HCofThisEvent*){ 
  G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  if(eventNumber%1==0 && PrtManager::Instance()->GetRunType()==0) std::cout<<"Event # "<<eventNumber <<std::endl;
  if(eventNumber%1000==0 && PrtManager::Instance()->GetRunType()!=0) std::cout<<"Event # "<<eventNumber <<std::endl;
  PrtManager::Instance()->Fill();
}

