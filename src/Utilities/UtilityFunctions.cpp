#include "UtilityFunctions.h"



namespace UtilityFunctions {

  double getTrackTheta( int trackId,
  						std::vector< std::vector<double> > *track_xpos,
                        std::vector< std::vector<double> > *track_ypos,
                        std::vector< std::vector<double> > *track_zpos){


  	double xDist = (*track_xpos)[trackId][3]-(*track_xpos)[trackId][0];
  	double yDist = (*track_ypos)[trackId][3]-(*track_ypos)[trackId][0];
  	double zDist = (*track_zpos)[trackId][3]-(*track_zpos)[trackId][0];
  	double d = sqrt(pow(xDist,2)+pow(yDist,2));
  	double Theta = atan(d/zDist);

  	return Theta;

  };

  double getTrackPhi( 	int trackId,
  						std::vector< std::vector<double> > *track_xpos,
                        std::vector< std::vector<double> > *track_ypos){


  	double xDist = (*track_xpos)[trackId][3] - (*track_xpos)[trackId][0];
  	double yDist = (*track_ypos)[trackId][3] - (*track_ypos)[trackId][0];
  	double Phi = atan2(yDist,xDist);


    if(Phi < 0){Phi = Phi + 8 * atan(1);}

  	return Phi;

  };
  double pointDistance(double x1, double y1, double z1, double x2, double y2, double z2){

    double dist = sqrt( pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
    return dist;

      

  }


    

};
