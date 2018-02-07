
#ifndef BEAMSELECTOR_H
#define BEAMSELECTOR_H
#include <vector>
#include <iostream>



class BeamSelector {

 public:

  BeamSelector( );
 
  bool isProton(std::vector< std::vector<double> > *track_zpos, int ntracks_reco, bool mc_mode);

 private:


};

#endif

 
