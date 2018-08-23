# protons
Proton-Ar Inelastic XS.

To build and run a test:

mkdir bin
mkdir lib
make clean
make
make
./bin/ProtonXsec.exe jobOptions/jobOptions.wirehit.in


# general info

anatree root file location is specified in files/Data_wire_ntuples.dat
modify this file to use different anaTree.root

# From points.C ----> modules

* MC Analzyer should take the following blocks of code from the root macro:
	* Block A: the geant 4 loop where we fill truth curves
	* Block D: deciding whether or not entries into the num are signal or background
	* Block F: deciding whether or not entries into the denom are signal or background

* Beam Selector will need two modes (or we need two of them)
	* MC version: can grab the quick and dirty version from Block B.
	* Data version: should probably rewrite one. Mass peaks and stuff.

* Event Selector should take from the following:
	* Block C: Numerator! This just needs to know which track to look at (primary track # from Beam Selector).
	* Block E: Denominator! This guy needs to know if there was an interacting candidate, if so where it was (x,y,z).

* Not included in the above mentioned modules is (are?)
	* Block G: puts all this together and then fills a bunch of histograms. Interacting/Incident KE (signal/background). Makes an unfolding matrix. etc.
	* Block H: calculates some cross-sections!

