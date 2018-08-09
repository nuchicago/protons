from ROOT import TH1D, TH2D, TCanvas, TGraph, TGraph2D, TFile, TChain
from rootpy.tree import TreeChain
import numpy as np
from root_numpy import root2array , tree2array
import sys
import os.path

"""
command line inputs:
  $ python EventDisplay.py <list of root anatrees> <list of primaries>
"""


if len(sys.argv != 2):
  print "please provide necessary inputs \n $ python EventDisplay.py <list of root anatrees> <list of ProtonXsec primaries> "

else:

  print "#### Open Files"
  entries, primaries = np.loadtxt(sys.argv[1], delimiter = ",", unpack = True, skiprows = 1)

  f = open(sys.argv[0],"r")
  ntuplesList = f.split()


  numNtuples =  ntuples[0]
  

  

  chainList = []
  for line in range(1,numNtuples+1):
    ntupleName = ntuplesList[line]
    print "ntuple " + str(line) + " : " + ntupleName
    ntupleName += "/anatree/anatree"
    chainList.append(ntupleName)
  
  tupleChain =  TreeChain("tupleChain", chainList )

  f.close()

  for jentry in  entries:
    ientry  = tupleChain.LoadEntry(jentry)
    if ientry < 0:
      continue
    nb = tupleChain.GetEntry(jentry)
    track_xpos = tupleChain.track_xpos
    print track_xpos



#~/.virtualenvs/pyrootenv/bin/python rootpy/setup.py install












