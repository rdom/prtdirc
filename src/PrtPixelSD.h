// -----------------------------------------
// PrtPixelSD.h
//
// author  : r.dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtPixelSD_h
#define PrtPixelSD_h 1

#include "G4VSensitiveDetector.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

class PrtPixelSD : public G4VSensitiveDetector {
 public:
  PrtPixelSD(const G4String &name, const G4String &hitsCollectionName, int nofCells);
  virtual ~PrtPixelSD();

  // methods from base class
  virtual void Initialize(G4HCofThisEvent *hitCollection);
  virtual bool ProcessHits(G4Step *step, G4TouchableHistory *history);
  virtual void EndOfEvent(G4HCofThisEvent *hitCollection);

 private:
  int fMcpLayout, fRunType;
  int fNofCells;
  double fQe_space[15][64];
  int fMultHit[15][256];
  int fMap_Mpc[15][256];
};

#endif
