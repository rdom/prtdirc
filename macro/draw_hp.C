#include "../../prttools/PrtTools.h"

void draw_hp(TString infile = "../build/hits.root") {

  PrtTools t(infile);

  while(t.next()){
    for (auto hit : t.event()->getHits()) {
      int mcp = hit.getPmt();
      int pix = hit.getPixel();
      double time = hit.getLeadTime();
      int ch = hit.getChannel();
            
      if (t.pid() == 2) t.fill_digi(mcp,pix);
    }
  }

  auto cdigi = t.draw_digi();
  t.add_canvas(cdigi);
  t.save_canvas("data/drawHP", 0);
}
