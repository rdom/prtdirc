import sys
import ROOT
import numpy as np

ROOT.gInterpreter.ProcessLine('#include "../../prttools/PrtTools.h"')
ROOT.gSystem.Load('../build/libPrt.so')

t = ROOT.PrtTools("/home/drc/data/jul18/403/beam_18215163732S.root") #  ../../prttools/hits.root
# t = ROOT.PrtTools("/home/drc/data/jul18/403/beam_18215164946S.root")

stat = int(sys.argv[1])
ne_train = stat
x_train = np.zeros((ne_train,100,2))
y_train = np.zeros((ne_train,1))

x_test = np.zeros((4000,100,2))
y_test = np.zeros((4000,1))

while t.next() and t.i() < ne_train + 4000 :
    i = t.i()
    ind = 0
    for hit in t.event().getHits() :
        pmt = hit.getPmt()
        pix = hit.getPixel()
        time = hit.getLeadTime()        
        ch = int(hit.getChannel())
        
        if time > 30 or time < 0 or pmt > 7 :
            continue

        time_bin = int(4*time)

        for j in range(20):
            if time_bin >= 50:
                time_bin -= 50

        if time_bin > 50:
            continue

        tx = int(8 * (pmt // 2) + pix % 8);
        ty = int(8 * (pmt % 2) + pix / 8);
        ic = tx * 16 + ty;

        if i < ne_train:
            x_train[i,ind, 0] = ic
            x_train[i,ind, 1] = time_bin
            y_train[i] = t.pid()
        else :
            x_test[i-ne_train,ind,0] = ic
            x_test[i-ne_train,ind,1] = time_bin
            y_test[i-ne_train] = t.pid()

        ind = ind + 1
        if(ind >= 100):
            break


nid = "data_stat_ind_f_" + str(ne_train);
np.savez(nid, x_train= x_train, y_train=y_train, x_test=x_test, y_test= y_test)
