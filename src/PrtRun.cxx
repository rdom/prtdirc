#include "PrtRun.h"

PrtRun::PrtRun()
  : fShortInfo(""), fName(""), fId(0), fRunType(0), fStudy(0), fMc(0), fPhysList(0), fPid(0), fGeometry(0), fLens(0), fRadiator(0), fPmtLayout(0), fTrigger(0), fNpmt(0), fNpix(0), fTheta(0), fPhi(0),
    fMomentum(0), fPrismStepX(0), fPrismStepY(0), fBeamX(0), fBeamZ(0), fBeamSize(0), fTimeSigma(0), fSimOffset(0), fRadiatorL(0), fRadiatorW(0),
    fRadiatorH(0), fTest1(0), fTest2(0), fTest3(0) {}

TString PrtRun::getInfo() {
  TString info = fShortInfo;
  info += Form("Run type %d", fRunType);
  info += Form("Study %d", fStudy);
  info +=  Form("Physics list %d", fPhysList);
  info += Form("Particle  id %d", fPid);
  info += Form("Particle momentum %f", fMomentum);
  info += Form("Geometry id %d", fGeometry);
  info += Form("Lens  id %d", fLens);
  info += Form("Prism step X %f", fPrismStepX);
  info += Form("Prism step Y %f", fPrismStepY);
  info += Form("MCP's time resolution %f", fTimeSigma);
  return info;
}

ClassImp(PrtRun)
