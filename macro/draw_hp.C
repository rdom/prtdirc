#include "../../prttools/PrtTools.h"

void draw_hp(TString infile = "../build/hits.root") {

  PrtTools t(infile);

  while (t.next() && t.i() < 10000) {
    for (auto hit : t.event()->getHits()) {
      int ch = hit.getChannel();
      int pmt = hit.getPmt();
      int pix = hit.getPixel();
      double time = hit.getLeadTime();

      if (t.pid() == 2) t.fill_digi(pmt, pix);
    }
  }

  auto cdigi = t.draw_digi();
  t.add_canvas(cdigi);
  t.save_canvas("data/drawHP", 0);

  t.write_string("digi.csv", t.pix_digi("m,p,v\n"));  
}
