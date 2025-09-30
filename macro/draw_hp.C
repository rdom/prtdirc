#if defined(__ACLIC__)
#include "../../prttools/PrtTools.h"
#else
R__LOAD_LIBRARY(../build/libPrt.so)
#endif

void draw_hp(TString infile = "../build/hits.root") {

  PrtTools t(infile);
  // double theta = t.run()->getTheta();
  int n = 0;

  TH1F *h1 = new TH1F("h1", "h1", 100, 0, 30);
  TH1F *h2 = new TH1F("h2", "h2", 100, 0, 30);

  while (t.next() && t.i() < 500000) {
    bool bl = 0;

    for (auto hit : t.event()->getHits()) {
      if (hit.getChannel() == 517) bl = true;
    }
    int pid = t.event()->getPid();
    for (auto hit : t.event()->getHits()) {
      int ch = hit.getChannel();
      int pmt = hit.getPmt();
      int pix = hit.getPixel();

      double time = hit.getLeadTime();

      if (pmt == 3 && pix == 32) {
        if (pid == 2) h1->Fill(time);
        else h2->Fill(time);
      }

      if (pid == 2) t.fill_digi(pmt, pix);
    }    
    if (pid == 2) n++;
    // if (n > 5000) break;
  }
  
  auto cdigi = t.draw_digi();
  t.add_canvas(cdigi);
  t.save_canvas("data/drawHP", 1);

  t.add_canvas("pdf");
  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);
  h1->Draw();
  h2->Draw("same");
  
  t.write_string("digi.csv", t.pix_digi("m,p,v\n"));
}
