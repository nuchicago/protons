from ROOT import TGraph, gStyle, gROOT, TROOT, TFile, gDirectory, TTree, TCanvas
from array import array
import numpy as np

#filepath  = raw_input("path to primary ID file")
#rootpath = raw_input("path to root data file")

rootpath = '../files/data_anaTree.root'

gROOT.SetBatch(True)
c1 = TCanvas("c1","c1",1000,1000)
c1.Divide(2,1)

#EventDisplayXZ = TH2D("EventDisplayXZ","Beam Selected Event",100,,120,0,48)
#EventDisplayYZ = TH2D("EventDisplayYZ","Beam Selected Event",100,)

eventx = array('d')
eventy = array('d')
eventz = array('d')

#myfile = open(filepath,'r')
rootfile  = TFile(rootpath, "READ")
tree = rootfile.Get("anatree/anatree")
entries = tree.GetEntriesFast()

for entry in tree:
  track_xpos = entry.track_xpos
  track_zpos =track_zpos

for i in range(ntracks_reco):
  for point in len(track_zpos[i])
  EventDisplayXZ[i] = TGraph(len(eventz),eventz, eventx)
  EventDisplayYZ[i] = TGraph(len(eventz),eventx,eventy)