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


x_train = np.zeros((ne_train,100,3))
y_train = np.zeros((ne_train,1))

x_test = np.zeros((4000,100,3))
y_test = np.zeros((4000,1))

while t.next() and t.i() < ne_train + 4000 :
    i = t.i()
    ind = 0
    for hit in t.event().getHits() :
        # tof = e.getTof()
        # tofPi = fEvent->getTofPi()
        # tofP = fEvent->getTofP()
        pmt = hit.getPmt()
        pix = hit.getPixel() - 1
        time = hit.getLeadTime()
        
        ch = int(hit.getChannel())
        
        if time > 30 :
            continue

        # if t.pid() == 2 :
        #     time = 5
        # else :
        #     time = 10

        time_bin = int(5*time)

        tx = int(8 * (pmt // 2) + pix % 8);
        ty = int(8 * (pmt % 2) + pix / 8);

        if i < ne_train:
            x_train[i,ind, 1] = ch
            x_train[i,ind, 2] = time_bin
            y_train[i] = t.pid()
        else :
            x_test[i-ne_train,ind,1] = ch
            x_test[i-ne_train,ind,2] = time_bin
            y_test[i-ne_train] = t.pid()

        ind = ind + 1
        if(ind >= 100):
            break

b = 0
for e in range(ne_train):
    if b>=32:
        b=0

    for i in range(100):
        x_train[e,i, 0] = b
    b = b + 1

b = 0
for e in range(4000):    
    if b>=32:
        b=0
    for i in range(100):
        x_test[e,i, 0] = b
    b = b + 1
    
nid = "data_stat_ind_" + str(ne_train);
np.savez(nid, x_train= x_train, y_train=y_train, x_test=x_test, y_test= y_test)
