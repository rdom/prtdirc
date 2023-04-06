import sys
import ROOT
import numpy as np

ROOT.gInterpreter.ProcessLine('#include "../../prttools/PrtTools.h"')
# ROOT.gSystem.Load('../../prttools/PrtTools_cxx.so')
ROOT.gSystem.Load('../build/libPrt.so')

t = ROOT.PrtTools("/home/drc/data/jul18/403/beam_18215163732S.root") #  ../../prttools/hits.root
# t = ROOT.PrtTools("/home/drc/data/jul18/403/beam_18215164946S.root")

stat = int(sys.argv[1])
ne_train = stat


# x_train = np.zeros((ne_train,8,64))
# y_train = np.zeros((ne_train,1))
# x_test = np.zeros((4000,8,64))
# y_test = np.zeros((4000,1))


x_train = np.zeros((ne_train,16,32))
y_train = np.zeros((ne_train,1))

x_test = np.zeros((4000,16,32))
y_test = np.zeros((4000,1))

while t.next() and t.i() < ne_train + 4000 :
    i = t.i()
    for hit in t.event().getHits() :
        # tof = e.getTof()
        # tofPi = fEvent->getTofPi()
        # tofP = fEvent->getTofP()
        pmt = hit.getPmt()
        pix = hit.getPixel() - 1
        time = hit.getLeadTime()

        if time > 30 :
            continue

        tx = int(8 * (pmt // 2) + pix % 8);
        ty = int(8 * (pmt % 2) + pix / 8);

        if i < ne_train:
            x_train[i,ty,tx] = time
            y_train[i] = t.pid()
        else :
            x_test[i-ne_train,ty,tx] = time
            y_test[i-ne_train] = t.pid()

        # if i < ne_train:
        #     x_train[i,hit.getPmt(),hit.getPixel()-1] = time
        #     y_train[i] = t.pid()
        # else :
        #     x_test[i-ne_train,hit.getPmt(),hit.getPixel()-1] = time
        #     y_test[i-ne_train] = t.pid()
        

nid = "data_stat_16x32_t_" + str(ne_train);
np.savez(nid, x_train= x_train, y_train=y_train, x_test=x_test, y_test= y_test)
