
#ifndef BEAMSELECTOR_H
#define BEAMSELECTOR_H
#include <vector>
#include <iostream>
#include "../Utilities/UtilityFunctions.h"


class BeamSelector {

 public:

  BeamSelector( );
 
  bool PrimaryTrack(std::vector< std::vector<double> > *track_zpos, int ntracks_reco,
  					 double zPointCutoff, int& reco_primary, double& first_reco_z);

  int  isTPCPrimary(std::vector< std::vector<double> > *track_zpos, int ntracks_reco,
  				 bool mc_mode, double zPointCutoff, int& reco_primary, double& first_reco_z);
  
  std::vector< std::vector<double> > wcTPCMatch(double wc_x, double wc_y, double wc_theta, double wc_phi,
                  std::vector< std::vector<double> > *track_xpos,
                  std::vector< std::vector<double> > *track_ypos,
                  std::vector< std::vector<double> > *track_zpos, int ntracks_reco);

 private:


};

#endif

 
