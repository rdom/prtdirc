#include <TTree.h>
#include <TFile.h>
#include <TVector3.h>
#include <TInterpreter.h>
#include <TClonesArray.h>
#include "../../prttools/PrtTools.h"

#include "../src/PrtLutNode.h"

void lutmean_cs(TString inFile = "../data/lut.root") {
  TString outFile = inFile.Copy().ReplaceAll(".root", ".cs_avr.root");

  PrtTools t;  
  TClonesArray *fLutNew = new TClonesArray("PrtLutNode");
  TTree *fTreeNew = new TTree("prtlut", "Look-up table for DIRC. Averaged");
  fTreeNew->Branch("LUT", &fLutNew, 256000, 2);  
  
  TClonesArray &fLutaNew = *fLutNew;
  for (Long64_t n = 0; n < t.maxdircch(); n++) {
    new ((fLutaNew)[n]) PrtLutNode(-1);
  }

  TFile *f = TFile::Open(inFile, "READ");
  TTree *tree = (TTree *)f->Get("prtlut");
  TClonesArray *fLut = new TClonesArray("PrtLutNode");
  tree->SetBranchAddress("LUT", &fLut);
  tree->GetEntry(0);

  std::vector<TVector3> vArray[100];
  std::vector<TVector3> lArray[100];
  std::vector<Double_t> tArray[100];
  std::vector<Double_t> pArray;

  Double_t cut(0), weight[9];
  TVector3 dir, pos, sum[9];
  Double_t angle, minangle, pathid, time, sumt;
  PrtLutNode *node;

  for (Int_t inode = 0; inode < fLut->GetEntriesFast(); inode++) {
    if (inode % 100 == 0) cout << "Node # " << inode << endl;
    node = (PrtLutNode *)fLut->At(inode);
    Int_t size = node->Entries();
    if (size < 1) continue;
    for (int i = 0; i < size; i++) {
      dir = node->GetEntry(i);
      pos = node->GetHitPos(i);
      time = node->GetTime(i);
      pathid = node->GetPathId(i);

      bool newid = true;
      for (uint j = 0; j < pArray.size(); j++) {
        if (pathid == pArray[j]) {
          vArray[j].push_back(dir);
          lArray[j].push_back(pos);
          tArray[j].push_back(time);
          newid = false;
        }
      }
      if (newid) {
        vArray[pArray.size()].push_back(dir);
        lArray[pArray.size()].push_back(pos);
        tArray[pArray.size()].push_back(time);
        pArray.push_back(pathid);
      }
    }

    for (uint j = 0; j < pArray.size(); j++) {
      for (Int_t s = 0; s < 9; s++) {
        weight[s] = 0;
        sum[s] = TVector3(0, 0, 0);
      }

      sumt = 0;
      for (uint v = 0; v < vArray[j].size(); v++) {
        sum[0] += vArray[j][v];
        weight[0]++;
        if (lArray[j][v].X() < cut) {
          sum[1] += vArray[j][v];
          weight[1]++;
        }
        if (lArray[j][v].X() < cut && lArray[j][v].Y() > cut) {
          sum[2] += vArray[j][v];
          weight[2]++;
        }
        if (lArray[j][v].Y() > cut) {
          sum[3] += vArray[j][v];
          weight[3]++;
        }
        if (lArray[j][v].X() > cut && lArray[j][v].Y() > cut) {
          sum[4] += vArray[j][v];
          weight[4]++;
        }

        if (lArray[j][v].X() > cut) {
          sum[5] += vArray[j][v];
          weight[5]++;
        }
        if (lArray[j][v].X() > cut && lArray[j][v].Y() < cut) {
          sum[6] += vArray[j][v];
          weight[6]++;
        }
        if (lArray[j][v].Y() < cut) {
          sum[7] += vArray[j][v];
          weight[7]++;
        }
        if (lArray[j][v].X() < cut && lArray[j][v].Y() < cut) {
          sum[8] += vArray[j][v];
          weight[8]++;
        }

        sumt += tArray[j][v];
      }

      if (weight[0] < 5) continue;

      for (int s = 0; s < 9; s++) {
        if (weight[s] == 0) {
          sum[s] = sum[0];
          continue;
        }
        sum[s] *= 1 / weight[s];
        sum[s] = sum[s].Unit();
      }

      sumt *= 1 / weight[0];

      ((PrtLutNode *)(fLutNew->At(inode)))
        ->AddEntry(node->GetDetectorId(), sum[0], pArray[j], node->GetNRefl(0), sumt, node->GetDigiPos(), node->GetDigiPos(), vArray[j].size() / (Double_t)size, sum[1], sum[2],
                   sum[3], sum[4], sum[5], sum[6], sum[7], sum[8], node->GetDigiPos());
    }
    for (int i = 0; i < 100; i++) {
      vArray[i].clear();
      tArray[i].clear();
      lArray[i].clear();
    }
    pArray.clear();
  }

  TFile *fFileNew = TFile::Open(outFile, "RECREATE");
  fTreeNew->Fill();
  fTreeNew->Write();
}
