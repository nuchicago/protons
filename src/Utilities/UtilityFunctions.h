#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H


#include <vector>
#include <iostream>
#include "math.h"



namespace UtilityFunctions {

  double getTrackTheta( int trackId,
                        std::vector< std::vector<double> > *track_xpos,
                        std::vector< std::vector<double> > *track_ypos,
                        std::vector< std::vector<double> > *track_zpos);

  double getTrackPhi(   int trackId,
                        std::vector< std::vector<double> > *track_xpos,
                        std::vector< std::vector<double> > *track_ypos);

};

#endif