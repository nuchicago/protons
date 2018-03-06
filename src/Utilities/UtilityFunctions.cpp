#include "UtilityFunctions.h"



namespace UtilityFunctions {

  double getTrackTheta( int trackId,
  						std::vector< std::vector<double> > *track_xpos,
                        std::vector< std::vector<double> > *track_ypos,
                        std::vector< std::vector<double> > *track_zpos){


  	double xDist = (*track_xpos)[trackId][1]-(*track_xpos)[trackId][0];
  	double yDist = (*track_ypos)[trackId][1]-(*track_ypos)[trackId][0];
  	double zDist = (*track_zpos)[trackId][1]-(*track_zpos)[trackId][0];
  	double d = sqrt(pow(xDist,2)+pow(yDist,2));
  	double Theta = atan(d/zDist);

  	return Theta;

  };

  double getTrackPhi( 	int trackId,
  						std::vector< std::vector<double> > *track_xpos,
                        std::vector< std::vector<double> > *track_ypos){


  	double xDist = (*track_xpos)[trackId][0];
  	double yDist = (*track_ypos)[trackId][0];
  	double Phi = atan(yDist/xDist);

  	return Phi;

  };

};
