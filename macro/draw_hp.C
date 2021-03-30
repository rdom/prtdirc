#include "../../prttools/PrtTools.h"

void draw_hp(TString infile = "../build/hits.root") {

  PrtTools t;
  if (!t.init_run(infile, 1)) return;

  for (int ievent = 0; ievent < t.entries(); ievent++) {
    t.next(ievent, 1000);
    for (auto hit : t.event()->getHits()) {
      int mcp = hit.getPmt();
      int pix = hit.getPixel();
      double time = hit.getLeadTime();
      int ch = hit.getChannel();
            
      if (mcp > 7) continue;
      if (t.pid() == 2) t.getdigi(mcp)->Fill(pix % 8, pix / 8);
    }
  }

  auto cdigi = t.draw_digi(0, 0);
  t.add_canvas(cdigi);
  t.save_canvas("data/drawHP", 0);
}
