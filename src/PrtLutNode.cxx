#include "PrtLutNode.h"

ClassImp(PrtLutNode)

// -----   Default constructor   -------------------------------------------
PrtLutNode::PrtLutNode(){ 
  fSize = 0;
}

// -----   Standard constructors   -----------------------------------------
PrtLutNode::PrtLutNode(Int_t nodeId):fDetectorId(nodeId),fSize(0){ 
}

void PrtLutNode::AddEntry(Int_t detectorId, TVector3 dir, Double_t pathid, Int_t nrefl, Double_t time, TVector3 hitpos, TVector3 digipos, Double_t weight) {
  fDetectorId = detectorId; 
  fNodeArray.push_back(dir);
  fPathIdArray.push_back(pathid);
  fWeightArray.push_back(weight);
  fNRefl.push_back(nrefl);
  fTimeArray.push_back(time);
  fHitPos.push_back(hitpos);
  fDigiPos = digipos;
  fSize++;
}
