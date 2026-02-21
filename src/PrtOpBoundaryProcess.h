// -----------------------------------------
// PrtOpBoundaryProcess.h
//
// author  : r.dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtOpBoundaryProcess_h
#define PrtOpBoundaryProcess_h

#include "globals.hh"
#include "G4OpBoundaryProcess.hh"

class PrtOpBoundaryProcess : public G4OpBoundaryProcess {
 public:
  PrtOpBoundaryProcess();
  ~PrtOpBoundaryProcess(){};

  // added to satisfy vtable when linked Geant4 lacks G4OpBoundaryProcess::SetInvokeSD
  virtual void SetInvokeSD(G4bool flag) override;


 public:
  G4VParticleChange *PostStepDoIt(const G4Track &aTrack, const G4Step &aStep);

 private:
  int fLensId;
  int fRunType;
};

#endif
