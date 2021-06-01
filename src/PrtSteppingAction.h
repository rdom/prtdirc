// -----------------------------------------
// PrtSteppingAction.h
//
// author  : r.dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtSteppingAction_h
#define PrtSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class PrtSteppingAction : public G4UserSteppingAction {
 public:
  PrtSteppingAction();
  virtual ~PrtSteppingAction();

  // method from the base class
  virtual void UserSteppingAction(const G4Step *);

 private:
  int fRunType;
  int fScintillationCounter;
  int fCerenkovCounter;
  int fEventNumber;
};

#endif
