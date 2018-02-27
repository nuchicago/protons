
#ifndef BEAMSELECTOR_H
#define BEAMSELECTOR_H
#include <vector>
#include <iostream>



class BeamSelector {

 public:

  BeamSelector( );
 
  bool PrimaryTrack(std::vector< std::vector<double> > *track_zpos, int ntracks_reco,
  					 double zPointCutoff, int& reco_primary, double& first_reco_z);

 private:


};

#endif

 
