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
  double pointDistance(double x1, double y1, double z1, 
                        double x2, double y2, double z2);
  std::vector<int> zOrderedTrack(std::vector<std::vector<double>> *track_zpos, int trackId,
                                    std::vector<int> *ntrack_hits);
  std::vector<int> zOrderedTrack2(std::vector<std::vector<double>> *track_zpos, int trackId,
                                    std::vector<int> *ntrack_hits);

  std::vector<double> pointProjector(std::vector<double> *point0, std::vector<double> *point1, 
                      double planeCz, double planeC0);

};

#endif
