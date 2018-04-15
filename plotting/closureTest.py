from ROOT import *
import os
import math
import argparse

from ROOT import TEfficiency
from ROOT import gStyle , TCanvas , TGraphErrors
from array import array


kineticEnergy = []
crossSec      = []
crossSec_el   = []
crossSec_inel = []
zero          = []



fileSlabs  = TFile('../files/MCOutput.root')

hDense   = fileSlabs.Get("hxs")
hSlabs   = fileSlabs.Get("sxs")
#hReco    = fileSlabs.Get("hreco_xs")
#hBad   	 = fileSlabs.Get("hreco_bad")

hDense.SetMarkerColor(kGreen-2)
hDense.SetLineColor(kGreen-2)

hSlabs.SetMarkerColor(kBlue)
hSlabs.SetLineColor(kBlue)

#hReco.SetMarkerColor(kBlack)
#hReco.SetLineColor(kBlack)
#
#hBad.SetMarkerColor(kBlack)
#hBad.SetLineColor(kBlack)

hDense.SetMarkerStyle(22)
hSlabs.SetMarkerStyle(23)
#hReco.SetMarkerStyle(23)
#hBad.SetMarkerStyle(23)

hDense.SetMarkerSize(.72)
hSlabs.SetMarkerSize(.72)
#hReco.SetMarkerSize(.72)
#hBad.SetMarkerSize(.72)

parser = argparse.ArgumentParser()
parser.add_argument("fname"   , nargs='?', default = 'ProtonPlusG4.txt', type = str, help="insert fileName")
args    = parser.parse_args()
fname   = args.fname

def is_number(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


title = ""
with open(fname) as f:
    for fLine in f.readlines():
        w = fLine.split()
        if is_number(w[0]):
            runIn    = int(w[0])
            ke       = float(w[1])
            xs_el       = float(w[2])
            xs_in       = float(w[3])
            xstot       = float(w[4])
            if ke > 1000:
                continue		
            kineticEnergy.append(ke)
            crossSec.append(xstot)
            crossSec_el.append(xs_el)
            crossSec_inel.append(xs_in)
            zero.append(0.)
        else:
            if "for" not in fLine: 
                continue
            title =  fLine[9:]




c1=TCanvas("c1" ,"Data" ,200 ,10 ,500 ,500) #make nice
c1.SetGrid ()
#define some data points . . .
x      = array('f', kineticEnergy )
y      = array('f', crossSec)
y_el   = array('f', crossSec_el)
y_inel = array('f', crossSec_inel)
exl    = array('f', zero)
exr    = array('f', zero)


nPoints=len(x)
# . . . and hand over to TGraphErros object
#gr      = TGraphErrors ( nPoints , x , y     , exl, exr )
#grel    = TGraphErrors ( nPoints , x , y_el  , exl, exr )
grinel  = TGraphErrors ( nPoints , x , y_inel, exl, exr )
grinel.SetTitle("P-Ar Inelastic Cross Section; Kinetic Energy [MeV]; Cross Section [barn]")
grinel . GetXaxis().SetRangeUser(0,1000)
grinel . GetYaxis().SetRangeUser(0,1.3)

grinel . GetYaxis().SetTitleOffset(1.2)

#gr . SetLineWidth(2) ;
#grel . SetLineWidth(2) ;
grinel . SetLineWidth(2) ;

#gr . SetLineColor(kGreen-2) ;
#grel . SetLineColor(kBlue) ;
grinel . SetLineColor(kRed) ;

#gr . SetFillColor(0) ;
#grel . SetFillColor(0) ;
grinel . SetFillColor(0) ;

#gr . Draw ( "APL" ) ;
#grel . Draw ( "PLsame" ) ;
#grinel . Draw ( "PLsame" ) ;
grinel . Draw ( "APL" ) ;

# ~~~ which to draw ~~~
hDense.Draw("same")
hSlabs.Draw("same")
#hReco.Draw("same")
#hBad.Draw("same")

#hIne.Draw("same")

legend = TLegend(.14,.68,.54,.88)


legend.AddEntry(grinel,"G4Prediction: Bertini Cascade")
legend.AddEntry(hDense,"G4 Spt XS: z=0.03cm")
legend.AddEntry(hSlabs,"G4 Sparse Spts XS: z=0.5cm")
#legend.AddEntry(hReco,"Reconstructed XS: z=0.5cm")
#legend.AddEntry(hBad,"Incorrectly Corrected Reconstructed XS: z=0.5cm")

legend.Draw("same")


c1 . Update ()



raw_input()  
