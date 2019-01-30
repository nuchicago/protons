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
  std::vector<int> zOrderedTrack(std::vector<std::vector<double>> *track_zpos, int trackId,
                                    std::vector<int> *ntrack_hits){


    /////////////////////////////////////////////////////////////////////
    //
    //  Orders track by z position
    //  Inputs: track_zpos: vector of vectors with reconsted track z position points
    //          trackId: position of target track in vector of track points
    //          ntrack_hits: vector, number of hit points per reconstructed track
    //
    //  Returns: vector with indices of first, second, last point ordered by z value
    //
    /////////////////////////////////////////////////////////////////////


    double zmin1 = (*track_zpos)[trackId][0]; 
    double zmin2 = (*track_zpos)[trackId][1];
    int index1 = 0;
    int index2 = 1;
    int indexLast = (*ntrack_hits)[trackId] - 1;

    //std::cout << "zmin1 : " << zmin1<< std::endl;
    //std::cout << "zmin2 : " << zmin2<< std::endl;

    if (zmin2 < zmin1){
      double zbuffer  = zmin2; 
      zmin2 = zmin1; 
      zmin1 = zbuffer; 
      int indexbuffer = index2;
      index2 = index1;
      index1 = indexbuffer;  
      //std::cout << "reordering" << std::endl;
      //std::cout << "zmin1 : " << zmin1<< std::endl;
      //std::cout << "zmin2 : " << zmin2<< std::endl;
    }
    
    double zmax = (*track_zpos)[trackId][indexLast];
    
   //std::cout << "zmax : " << zmax<< std::endl;

    for (int ipoint = 2 ; ipoint < (*ntrack_hits)[trackId] ; ipoint++){
      double zval  = (*track_zpos)[trackId][ipoint];
      if (zval < zmin1) {
        zmin2 = zmin1;
        zmin1  = zval;
        index1  =  ipoint;
      }
      else if (zval < zmin2){
        zmin2 = zval;
        index2 =  ipoint;
      }
      if (zval > zmax){
        zmax = zval;
        indexLast = ipoint;
      } 
    }
     

    std::vector<int> trackInd = {index1, index2, indexLast};


    return trackInd;
  }

  std::vector<double> pointProjector(std::vector<double> *point0, std::vector<double> *point1, 
                      double planeCz, double planeC0){


  std::cout << "function parameters:" << planeCz << ", " << planeC0 << std::endl; 
  double t_param = ((*point0)[0] - planeCz * (*point0)[2] - planeC0)/ ( planeCz*((*point1)[2] - (*point0)[2]) - ((*point1)[0] - (*point0)[0]));


  //double t_param = (planeCz * (*point0)[2] - (*point0)[0] + planeC0) / (((*point1)[0] - (*point0)[0]) - planeCz * ((*point1)[2] - (*point0)[2]));
  std::cout << "t parameter : " << t_param << std::endl;
  double x_intercept = (*point0)[0] + t_param * ((*point1)[0] - (*point0)[0]);
  double y_intercept = (*point0)[1] + t_param * ((*point1)[1] - (*point0)[1]);
  double z_intercept = (*point0)[2] + t_param * ((*point1)[2] - (*point0)[2]);

  std::vector<double> output = {x_intercept, y_intercept, z_intercept};

  return output;

  }
};
