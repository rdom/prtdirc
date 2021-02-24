#include "PrtRun.h"

ClassImp(PrtRun)

// // -----   Default constructor   -------------------------------------------
PrtRun::PrtRun(): frunId(0), fstudyId(0),fmc(0),fphysList(0),fpid(0),fgeometry(0),
  flens(0),fradiator(0),ftrigger(0),ftheta(0),fphi(0),
  fmomentum(TVector3(0,0,0)),fposition(TVector3(0,0,0)),
  fprismStepX(0),fprismStepY(0),fbeamX(0),fbeamZ(0),fbeamSize(0),ftimeSigma(0),
  ftof(0),ftofPi(0),ftofP(0),ftest1(0),ftest2(0){ 
}

void PrtRun::PrintInfo(){
  std::cout<< Form("Run %d \n",frunId)<<std::endl;
  std::cout<< Form("Study %d \n",fstudyId)<<std::endl;
  std::cout<< Form("Physics list %d \n",fphysList)<<std::endl;
  std::cout<< Form("Particle  id %d \n",fpid)<<std::endl;
  std::cout<< Form("Particle momentum %f \n", fmomentum.Mag())<<std::endl;
  std::cout<< Form("Geometry id %d \n", fgeometry)<<std::endl;
  std::cout<< Form("Lens  id %d \n",    flens)<<std::endl;
  std::cout<< Form("Prism step X %f \n",    fprismStepX)<<std::endl;
  std::cout<< Form("Prism step Y %f \n",    fprismStepY)<<std::endl;
  std::cout<< Form("MCP's time resolution %f \n",    ftimeSigma)<<std::endl;
}
