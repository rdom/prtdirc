#include "PrtCherenkovProcess.h"

#include "PrtManager.h"

#include "G4Poisson.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

PrtCherenkovProcess::PrtCherenkovProcess(const G4String &processName, G4ProcessType type)
  : G4Cerenkov(processName, type) {

  fQEtype = 0;
  if (PrtManager::Instance()->getRun()->getPmtLayout() == 2030) fQEtype = 1;

  // create a detector efficiency function:
  {
    // quantum efficiency data from Alex Britting, Jan 25, 2011
    // unit is percent
    // first value is at 200 nm, last at 700 nm
    // credible range start around 250nm, >= 280nm to be safe
    double eff_std[502] = {
      0,     231.84, 615.36, 657.4, 258.78, 9839.92, 44.67, 67.87, 51.01, 41.49, 5.36,  49.4,
      2.13,  35.49,  8.66,   5.03,  7.51,   13.27,   18.71, 3.92,  3.66,  8.2,   0.56,  7.68,
      2.87,  10.06,  3.47,   3.39,  6.99,   6.01,    4.92,  6.25,  5.97,  6.92,  8.29,  10.45,
      8.68,  8.6,    9.79,   11.76, 9.53,   10.98,   9.9,   10.97, 11.31, 10.88, 10.78, 12.16,
      12.38, 12.37,  13.04,  12.36, 13.18,  13.7,    13.85, 13.66, 13.98, 14.55, 14.93, 14.82,
      14.97, 14.98,  15.14,  15.35, 15.37,  15.43,   15.49, 15.59, 15.84, 15.84, 15.92, 16.01,
      16.22, 16.41,  16.42,  16.52, 16.86,  17.1,    17.17, 17.22, 17.46, 17.79, 17.99, 18.13,
      18.33, 18.34,  18.53,  18.72, 18.95,  19.02,   19.15, 19.28, 19.45, 19.66, 19.69, 19.77,
      19.73, 19.95,  19.98,  20.17, 20.29,  20.33,   20.37, 20.47, 20.48, 20.57, 20.75, 20.8,
      20.84, 20.86,  20.88,  21.0,  21.06,  21.0,    21.06, 21.06, 21.04, 21.1,  21.14, 21.08,
      21.17, 21.3,   21.38,  21.49, 21.58,  21.69,   21.77, 21.87, 22.02, 22.13, 22.29, 22.35,
      22.45, 22.53,  22.55,  22.64, 22.67,  22.73,   22.74, 22.71, 22.79, 22.76, 22.77, 22.76,
      22.75, 22.78,  22.7,   22.68, 22.72,  22.66,   22.64, 22.7,  22.67, 22.71, 22.67, 22.75,
      22.77, 22.83,  22.84,  22.93, 22.97,  23.0,    23.08, 23.16, 23.27, 23.25, 23.37, 23.44,
      23.49, 23.55,  23.52,  23.58, 23.64,  23.63,   23.58, 23.64, 23.63, 23.62, 23.64, 23.63,
      23.66, 23.59,  23.59,  23.56, 23.58,  23.63,   23.57, 23.66, 23.62, 23.67, 23.64, 23.54,
      23.57, 23.51,  23.53,  23.45, 23.3,   23.41,   23.25, 23.21, 23.08, 23.01, 22.92, 22.9,
      22.76, 22.76,  22.61,  22.53, 22.48,  22.39,   22.29, 22.24, 22.2,  22.12, 22.07, 21.96,
      21.89, 21.87,  21.76,  21.74, 21.58,  21.49,   21.48, 21.37, 21.29, 21.2,  21.17, 21.03,
      20.98, 20.92,  20.85,  20.76, 20.69,  20.58,   20.56, 20.47, 20.37, 20.32, 20.24, 20.13,
      20.08, 19.9,   19.84,  19.77, 19.69,  19.63,   19.51, 19.41, 19.27, 19.06, 19.01, 18.87,
      18.7,  18.49,  18.41,  18.17, 17.98,  17.84,   17.69, 17.5,  17.25, 17.15, 16.98, 16.79,
      16.66, 16.48,  16.32,  16.19, 16.02,  15.88,   15.77, 15.67, 15.5,  15.39, 15.23, 15.09,
      15.04, 14.92,  14.75,  14.7,  14.5,   14.45,   14.34, 14.25, 14.16, 14.13, 14.0,  13.92,
      13.84, 13.76,  13.73,  13.61, 13.54,  13.52,   13.45, 13.41, 13.39, 13.31, 13.22, 13.17,
      13.13, 13.06,  13.2,   13.09, 12.97,  12.92,   12.73, 12.65, 12.4,  12.22, 12.02, 11.79,
      11.59, 11.33,  11.03,  10.68, 10.46,  10.14,   9.88,  9.62,  9.36,  9.14,  8.87,  8.63,
      8.51,  8.24,   8.07,   7.88,  7.77,   7.65,    7.52,  7.35,  7.27,  7.21,  7.1,   6.92,
      6.89,  6.79,   6.74,   6.56,  6.54,   6.5,     6.39,  6.33,  6.25,  6.27,  6.14,  6.06,
      6.04,  6.01,   5.91,   5.89,  5.79,   5.75,    5.75,  5.67,  5.61,  5.51,  5.52,  5.43,
      5.43,  5.34,   5.31,   5.35,  5.23,   5.2,     5.14,  5.11,  5.11,  5.01,  4.98,  4.93,
      4.99,  4.89,   4.82,   4.87,  4.8,    4.7,     4.65,  4.65,  4.61,  4.49,  4.56,  4.44,
      4.42,  4.44,   4.35,   4.35,  4.27,   4.29,    4.19,  4.13,  4.08,  4.02,  4.07,  3.92,
      3.95,  3.88,   3.82,   3.86,  3.74,   3.71,    3.66,  3.72,  3.62,  3.55,  3.56,  3.57,
      3.45,  3.38,   3.36,   3.36,  3.28,   3.25,    3.19,  3.26,  3.13,  3.17,  3.15,  3.04,
      2.98,  2.93,   2.98,   2.9,   2.89,   2.9,     2.81,  2.74,  2.81,  2.68,  2.73,  2.7,
      2.57,  2.58,   2.55,   2.55,  2.37,   2.39,    2.39,  2.44,  2.37,  2.26,  2.27,  2.27,
      2.23,  2.26,   2.14,   2.08,  2.15,   2.06,    2.09,  2.04,  2.0,   1.95,  2.02,  1.87,
      1.9,   1.8,    1.87,   1.85,  1.87,   1.81,    1.86,  1.74,  1.74,  1.63,  1.59,  1.5,
      1.5,   1.44,   1.47,   1.32,  1.24,   1.28,    1.19,  1.21,  1.21,  1.1,   1.1,   1.05,
      1.06,  0.94,   0.92,   0.87,  0.92,   0.81,    0.86,  0.78,  0.77,  0.8,   0.67,  0.7,
      0.81,  0.61,   0.64,   0.71,  0.66,   0.67,    0.68,  0.69,  0.68,  0.73};

    double eff_400[251] = {
      0,    0,    14.0, 14.8, 14.5, 14.9, 14.4, 14.2, 13.9, 14.6, 15.2, 15.7, 16.4, 16.9, 17.5,
      17.7, 18.1, 18.8, 19.3, 19.8, 20.6, 21.4, 22.4, 23.1, 23.6, 24.1, 24.2, 24.6, 24.8, 25.2,
      25.7, 26.5, 27.1, 28.2, 29.0, 29.9, 30.8, 31.1, 31.7, 31.8, 31.6, 31.5, 31.5, 31.3, 31.0,
      30.8, 30.8, 30.4, 30.2, 30.3, 30.2, 30.1, 30.1, 30.1, 29.8, 29.9, 29.8, 29.7, 29.7, 29.7,
      29.8, 29.8, 29.9, 29.9, 29.8, 29.9, 29.8, 29.9, 29.8, 29.7, 29.8, 29.7, 29.8, 29.6, 29.5,
      29.7, 29.7, 29.8, 30.1, 30.4, 31.0, 31.3, 31.5, 31.8, 31.8, 31.9, 32.0, 32.0, 32.0, 32.0,
      32.2, 32.2, 32.1, 31.8, 31.8, 31.8, 31.7, 31.6, 31.6, 31.7, 31.5, 31.5, 31.4, 31.3, 31.3,
      31.2, 30.8, 30.7, 30.5, 30.3, 29.9, 29.5, 29.3, 29.2, 28.6, 28.2, 27.9, 27.8, 27.3, 27.0,
      26.6, 26.1, 25.9, 25.5, 25.0, 24.6, 24.2, 23.8, 23.4, 23.0, 22.7, 22.4, 21.9, 21.4, 21.2,
      20.7, 20.3, 19.8, 19.6, 19.3, 18.9, 18.7, 18.3, 17.9, 17.8, 17.8, 16.7, 16.5, 16.4, 16.0,
      15.6, 15.6, 15.2, 14.9, 14.6, 14.4, 14.1, 13.8, 13.6, 13.3, 13.0, 12.8, 12.6, 12.3, 12.0,
      11.9, 11.7, 11.5, 11.2, 11.1, 10.9, 10.7, 10.4, 10.3, 9.9,  9.8,  9.6,  9.3,  9.1,  9.0,
      8.8,  8.5,  8.3,  8.3,  8.2,  7.9,  7.8,  7.7,  7.5,  7.3,  7.1,  6.9,  6.7,  6.6,  6.3,
      6.2,  6.0,  5.8,  5.7,  5.6,  5.4,  5.2,  5.1,  4.9,  4.8,  4.6,  4.5,  4.4,  4.2,  4.1,
      4.0,  3.8,  3.7,  3.5,  3.3,  3.2,  3.1,  3.0,  2.9,  2.5,  2.4,  2.4,  2.3,  2.3,  2.1,
      1.8,  1.6,  1.5,  1.5,  1.6,  1.8,  1.9,  1.4,  0.8,  0.9,  0.8,  0.7,  0.6,  0.3,  0.3,
      0.5,  0.3,  0.4,  0.3,  0.1,  0.2,  0.1,  0.2,  0.3,  0.0};

    double lambda[2][1000] = {{0}}, fEfficiencyR[2][1000] = {{0}}, fCollectionEff = 0.65;

    for (Int_t i = 1; i < 1000; i++) {
      // Std QE (Photonis 2)
      fCollectionEff = 0.65;
      lambda[0][i] = 200 + i - 1;
      fEfficiencyR[0][i] = (i > 40 && i < 500) ? eff_std[i] * 0.01 * fCollectionEff : 0;

      // HIQ 400 (Photonis 1)
      fCollectionEff = 0.95;
      lambda[1][i] = 180 + (i - 1) * 2;
      fEfficiencyR[1][i] = (i < 251) ? eff_400[i] * 0.01 * fCollectionEff : 0;
    }

    fDetEff[0] = new TGraph(1000, lambda[0], fEfficiencyR[0]);
    fDetEff[1] = new TGraph(1000, lambda[1], fEfficiencyR[1]);
  }
  fDetEff[0]->SetBit(TGraph::kIsSortedX);
  fDetEff[1]->SetBit(TGraph::kIsSortedX);
}

G4VParticleChange *PrtCherenkovProcess::PostStepDoIt(const G4Track &aTrack, const G4Step &aStep)

// This routine is called for each tracking Step of a charged particle
// in a radiator. A Poisson-distributed number of photons is generated
// according to the Cerenkov formula, distributed evenly along the track
// segment and uniformly azimuth w.r.t. the particle direction. The
// parameters are then transformed into the Master Reference System, and
// they are added to the particle change.

{
  //////////////////////////////////////////////////////
  // Should we ensure that the material is dispersive?
  //////////////////////////////////////////////////////

  aParticleChange.Initialize(aTrack);

  const G4DynamicParticle *aParticle = aTrack.GetDynamicParticle();
  const G4Material *aMaterial = aTrack.GetMaterial();

  G4StepPoint *pPreStepPoint = aStep.GetPreStepPoint();
  G4StepPoint *pPostStepPoint = aStep.GetPostStepPoint();

  G4ThreeVector x0 = pPreStepPoint->GetPosition();
  G4ThreeVector p0 = aStep.GetDeltaPosition().unit();
  G4double t0 = pPreStepPoint->GetGlobalTime();

  G4MaterialPropertiesTable *aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();
  if (!aMaterialPropertiesTable) return pParticleChange;

  G4MaterialPropertyVector *Rindex = aMaterialPropertiesTable->GetProperty("RINDEX");
  if (!Rindex) return pParticleChange;

  // particle charge
  const G4double charge = aParticle->GetDefinition()->GetPDGCharge();

  // particle beta
  const G4double beta = (pPreStepPoint->GetBeta() + pPostStepPoint->GetBeta()) / 2.;

  G4double MeanNumberOfPhotons = GetAverageNumberOfPhotons(charge, beta, aMaterial, Rindex);
  
  if (MeanNumberOfPhotons <= 0.0) {

    // return unchanged particle and no secondaries

    aParticleChange.SetNumberOfSecondaries(0);

    return pParticleChange;
  }

  G4double step_length;
  step_length = aStep.GetStepLength();

  MeanNumberOfPhotons = MeanNumberOfPhotons * step_length;

  G4int NumPhotons = (G4int)G4Poisson(MeanNumberOfPhotons);

  if (NumPhotons <= 0) {

    // return unchanged particle and no secondaries

    aParticleChange.SetNumberOfSecondaries(0);

    return pParticleChange;
  }

  ////////////////////////////////////////////////////////////////

  int runid = PrtManager::Instance()->getRun()->getRunType();
  if (runid == 5) {
    NumPhotons = 2000;
  }
  aParticleChange.SetNumberOfSecondaries(NumPhotons);

  if (20) { // fTrackSecondariesFirst
    if (aTrack.GetTrackStatus() == fAlive) aParticleChange.ProposeTrackStatus(fSuspend);
  }

  ////////////////////////////////////////////////////////////////

  // monochromatic photons
  G4double Pmin = 3.18 * 1E-6;
  G4double Pmax = Pmin;

  if (runid == 5) {
    Pmin = Rindex->GetMinLowEdgeEnergy();
    Pmax = Rindex->GetMaxLowEdgeEnergy();
  }

  G4double dp = Pmax - Pmin;

  G4double nMax = Rindex->GetMaxValue();

  G4double BetaInverse = 1. / beta;

  G4double maxCos = BetaInverse / nMax;
  G4double maxSin2 = (1.0 - maxCos) * (1.0 + maxCos);

  const G4double beta1 = pPreStepPoint->GetBeta();
  const G4double beta2 = pPostStepPoint->GetBeta();

  G4double MeanNumberOfPhotons1 = GetAverageNumberOfPhotons(charge, beta1, aMaterial, Rindex);
  G4double MeanNumberOfPhotons2 = GetAverageNumberOfPhotons(charge, beta2, aMaterial, Rindex);

  if (runid == 5) {
    MeanNumberOfPhotons1 = NumPhotons;
    MeanNumberOfPhotons2 = NumPhotons;
  }

  for (G4int i = 0; i < NumPhotons; i++) {

    // Determine photon energy

    G4double rand;
    G4double sampledEnergy, sampledRI;
    G4double cosTheta, sin2Theta;

    // sample an energy
    do {
      rand = G4UniformRand();
      sampledEnergy = Pmin + rand * dp;
      sampledRI = Rindex->Value(sampledEnergy);
      cosTheta = BetaInverse / sampledRI;

      if (cosTheta > 1) {
        std::cout << "Warning - PrtCherenkovProcess:  cosTheta " << cosTheta << std::endl;
        return pParticleChange;
      }

      sin2Theta = (1.0 - cosTheta) * (1.0 + cosTheta);
      rand = G4UniformRand();
 
    } while (rand * maxSin2 > sin2Theta);

    if (runid == 5) {
      double lambda = 197.0 * 2.0 * pi / (sampledEnergy * 1.0E6);
      if (G4UniformRand() > fDetEff[fQEtype]->Eval(lambda)) continue; // pmt's quantum efficiency
    }

    // Generate random position of photon on cone surface
    // defined by Theta

    rand = G4UniformRand();

    G4double phi = twopi * rand; // twopi*PrtManager::Instance()->GetTest()/360.;//twopi*rand;
    G4double sinPhi = std::sin(phi);
    G4double cosPhi = std::cos(phi);

    // calculate x,y, and z components of photon energy
    // (in coord system with primary particle direction
    //  aligned with the z axis)

    G4double sinTheta = std::sqrt(sin2Theta);
    G4double px = sinTheta * cosPhi;
    G4double py = sinTheta * sinPhi;
    G4double pz = cosTheta;

    // Create photon momentum direction vector
    // The momentum direction is still with respect
    // to the coordinate system where the primary
    // particle direction is aligned with the z axis

    G4ParticleMomentum photonMomentum(px, py, pz);

    // Rotate momentum direction back to global reference
    // system

    photonMomentum.rotateUz(p0); // test

    // Determine polarization of new photon

    G4double sx = cosTheta * cosPhi;
    G4double sy = cosTheta * sinPhi;
    G4double sz = -sinTheta;

    G4ThreeVector photonPolarization(sx, sy, sz);

    // Rotate back to original coord system

    photonPolarization.rotateUz(p0);

    // Generate a new photon:

    G4DynamicParticle *aCerenkovPhoton =
      new G4DynamicParticle(G4OpticalPhoton::OpticalPhoton(), photonMomentum);
    aCerenkovPhoton->SetPolarization(photonPolarization.x(), photonPolarization.y(),
                                     photonPolarization.z());

    aCerenkovPhoton->SetKineticEnergy(sampledEnergy);

    // Generate new G4Track object:

    G4double delta, NumberOfPhotons, N;

    do {
      rand = G4UniformRand(); // PrtManager::Instance()->GetTest();//G4UniformRand();
      delta = rand * aStep.GetStepLength();
      NumberOfPhotons = MeanNumberOfPhotons1 - delta *
                                                 (MeanNumberOfPhotons1 - MeanNumberOfPhotons2) /
                                                 aStep.GetStepLength();
      N = G4UniformRand() * std::max(MeanNumberOfPhotons1, MeanNumberOfPhotons2);
    } while (N > NumberOfPhotons);

    G4double deltaTime =
      delta / ((pPreStepPoint->GetVelocity() + pPostStepPoint->GetVelocity()) / 2.);

    G4double aSecondaryTime = t0 + deltaTime;

    G4ThreeVector aSecondaryPosition = x0 + rand * aStep.GetDeltaPosition();

    G4Track *aSecondaryTrack = new G4Track(aCerenkovPhoton, aSecondaryTime, aSecondaryPosition);

    aSecondaryTrack->SetTouchableHandle(aStep.GetPreStepPoint()->GetTouchableHandle());

    aSecondaryTrack->SetParentID(aTrack.GetTrackID());

    aParticleChange.AddSecondary(aSecondaryTrack);
  }

  if (verboseLevel > 0) {
    G4cout << "\n Exiting from G4Cerenkov::DoIt -- NumberOfSecondaries = "
           << aParticleChange.GetNumberOfSecondaries() << G4endl;
  }

  return pParticleChange;
}

G4double PrtCherenkovProcess::GetAverageNumberOfPhotons(const G4double charge, const G4double beta,
                                                        const G4Material *aMaterial,
                                                        G4MaterialPropertyVector *Rindex) const {
  
  const G4double Rfact = 369.81 / (eV * cm);

  if (beta <= 0.0) return 0.0;

  G4double BetaInverse = 1. / beta;

  // Vectors used in computation of Cerenkov Angle Integral:
  // 	- Refraction Indices for the current material
  //	- new G4PhysicsOrderedFreeVector allocated to hold CAI's

  G4int materialIndex = aMaterial->GetIndex();

  // Retrieve the Cerenkov Angle Integrals for this material

  G4PhysicsOrderedFreeVector *CerenkovAngleIntegrals =
    (G4PhysicsOrderedFreeVector *)((*thePhysicsTable)(materialIndex));

  if (!(CerenkovAngleIntegrals->IsFilledVectorExist())) return 0.0;
  
  // Min and Max photon energies
  G4double Pmin = Rindex->GetMinLowEdgeEnergy();
  G4double Pmax = Rindex->GetMaxLowEdgeEnergy();

  // Min and Max Refraction Indices
  G4double nMin = Rindex->GetMinValue();
  G4double nMax = Rindex->GetMaxValue();

  // Max Cerenkov Angle Integral
  G4double CAImax = CerenkovAngleIntegrals->GetMaxValue();

  G4double dp, ge;

  // If n(Pmax) < 1/Beta -- no photons generated

  if (nMax < BetaInverse) {
    dp = 0;
    ge = 0;
  }

  // otherwise if n(Pmin) >= 1/Beta -- photons generated

  else if (nMin > BetaInverse) {
    dp = Pmax - Pmin;
    ge = CAImax;
  }

  // If n(Pmin) < 1/Beta, and n(Pmax) >= 1/Beta, then
  // we need to find a P such that the value of n(P) == 1/Beta.
  // Interpolation is performed by the GetEnergy() and
  // Value() methods of the G4MaterialPropertiesTable and
  // the GetValue() method of G4PhysicsVector.

  else {
    Pmin = Rindex->GetEnergy(BetaInverse);
    dp = Pmax - Pmin;

    // need boolean for current implementation of G4PhysicsVector
    // ==> being phased out
    G4bool isOutRange;
    G4double CAImin = CerenkovAngleIntegrals->GetValue(Pmin, isOutRange);
    ge = CAImax - CAImin;

    if (verboseLevel > 0) {
      G4cout << "CAImin = " << CAImin << G4endl;
      G4cout << "ge = " << ge << G4endl;
    }
  }

  // Calculate number of photons
  G4double NumPhotons =
    Rfact * charge / eplus * charge / eplus * (dp - ge * BetaInverse * BetaInverse);

  return NumPhotons;
}
