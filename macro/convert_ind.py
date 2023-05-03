import sys
import ROOT
import numpy as np

ROOT.gInterpreter.ProcessLine('#include "../../prttools/PrtTools.h"')
ROOT.gSystem.Load('../build/libPrt.so')

# t = ROOT.PrtTools("/home/drc/data/jul18/403/beam_18215163732S.root") #  ../../prttools/hits.root
# t = ROOT.PrtTools("/home/drc/data/jul18/403/beam_18215164946S.root")


infile = "hits_403.root"
if(len(sys.argv) > 1):
    infile = sys.argv[1] 

t = ROOT.PrtTools(infile)

nhits = 100
timebins = 100
entries = t.entries()
print("running for ", entries, " events")
x_train = np.zeros((entries,nhits,2))
y_train = np.zeros((entries,1))


while t.next() and t.i() < entries :
    i = t.i()
    mom = t.event().getMomentum().Mag();
    theta = t.event().getTof();
    mom_bin = int(mom * 10)
    theta_bin = int(theta * 5)
    
    ind = 0
    for hit in t.event().getHits() :
        pmt = hit.getPmt()
        pix = hit.getPixel()
        time = hit.getLeadTime()        
        ch = int(hit.getChannel())
        
        if time > 30 or time < 0 or pmt > 7 :
            continue

        time_bin = int(4*time)

        # for j in range(20):
        #     if time_bin >= 50:
        #         time_bin -= 50

        if time_bin >= timebins:
            continue

        tx = int(8 * (pmt // 2) + pix % 8);
        ty = int(8 * (pmt % 2) + pix / 8);
        ic = tx * 16 + ty;

        x_train[i,ind, 0] = ic
        x_train[i,ind, 1] = time_bin
        y_train[i] = t.pid()

        ind = ind + 1
        if(ind >= nhits):
            break

        # encode momentum and theta angle
        x_train[i,98, 0] = 512
        x_train[i,98, 1] = mom_bin            
        x_train[i,99, 0] = 513 + int(theta_bin / nhits)
        x_train[i,99, 1] = theta_bin % nhits        
      
outfile = infile.replace(".root","")
np.savez(outfile, x_train= x_train, y_train=y_train)
