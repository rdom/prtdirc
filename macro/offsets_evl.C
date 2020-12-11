void offsets_evl(TString in = "~/data/jul18/401/lut_beam_401_20S_cs_avr.root"){
  
  TFile* f = new TFile(in);
  TTree *t = (TTree *) f->Get("prtlut") ;
  TClonesArray* clut = new TClonesArray("PrtLutNode");
  t->SetBranchAddress("LUT",&clut); 
  t->GetEntry(0);

  double dist,cns = 197.0;
  TVector3 dir, dif = TVector3(260,0,980.5); //diffuser positioin
  PrtLutNode *node;
  std::cout<<"clut->GetEntriesFast() "<<clut->GetEntriesFast()<<std::endl;

  TGraph *gr = new TGraph();
  
  for (int c=0; c<clut->GetEntriesFast(); c++){
    node = (PrtLutNode*) clut->At(c);
    int size = node->Entries();    
    if(size<1) continue;

    dir = node->GetHitPos(0);
    dist = (dir-dif).Mag();
    gr->SetPoint(c,c,dist/cns);
  }

  TString out = in.ReplaceAll("_cs_avr.root","_evl.root");  
  TFile ofile(out,"RECREATE");
  gr->SetName("evl_diffuser_1");
  gr->Write();
  ofile.Write();
  ofile.Close();
  std::cout<<"saved in "<<out<<std::endl;
  
    
}
