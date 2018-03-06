//#############################################################################
//
// BeamSelector.cpp
//
//  A desription of what this class does
//
// January 17, 2018
//
//#############################################################################


#include "BeamSelector.h"




//=============================================================================
// Constructor
//=============================================================================

BeamSelector::BeamSelector(){  

  
}


//=============================================================================
// classsifyEvent()
//=============================================================================
bool BeamSelector::PrimaryTrack(std::vector<std::vector<double>> *track_zpos,int ntracks_reco, 
  double zPointCutoff, int& reco_primary, double& first_reco_z){

      bool print = true;


      

        if(print){
          std::cout<<"Primary Track Lookup"<<std::endl;
          std::cout<<">>>>>>>>>>>>>>>>reco<<<<<<<<<<<<<<<<"<<std::endl;
        }

        for(int rtrack = 0; rtrack < ntracks_reco; rtrack++){
          double z1 = (*track_zpos)[rtrack][0];
          if(print){ 
            std::cout << "first z point: " << z1 << std::endl;
          }
          if(z1 < first_reco_z && z1 < zPointCutoff){
            first_reco_z = z1;
            reco_primary = rtrack;
          }
        }//<---End loop reco tracks
        if(reco_primary == -1){return false;}
        else{
          if(print){std::cout << "earliest primary: " << reco_primary << std::endl;} 
        return true;}
      

}

int BeamSelector::isTPCPrimary(std::vector<std::vector<double>> *track_zpos,int ntracks_reco,  bool mc_mode, 
  double zPointCutoff, int& reco_primary, double& first_reco_z){

      bool print = true;
      if(mc_mode){
        if(print){
          std::cout<<"MC Primary Selection"<<std::endl;
          std::cout<<">>>>>>>>>>>>>>>>reco<<<<<<<<<<<<<<<<"<<std::endl;
        }

        for(int rtrack = 0; rtrack < ntracks_reco; rtrack++) {
          double z1 = (*track_zpos)[rtrack][0];
          if(print){std::cout << "first z point: " << z1 << std::endl;}
          if(z1 < first_reco_z && z1 < zPointCutoff) {
            first_reco_z = z1;
            reco_primary = rtrack;
          }
        }//<---End loop reco tracks
        if(print){std::cout << "earliest primary: " << reco_primary << std::endl;} 
        return reco_primary;
      }

      // ## if data ##
      else {
        std::cout<<"Beam Selection for data not available\n"<<std::endl;
        return -999;
      }

}

std::vector<std::vector<double>> BeamSelector::wcTPCMatch(double wc_x, double wc_y, double wc_theta, double wc_phi,
                  std::vector< std::vector<double> > *track_xpos,
                  std::vector< std::vector<double> > *track_ypos,
                  std::vector< std::vector<double> > *track_zpos, int ntracks_reco){


  std::vector<std::vector<double>> vOut;
  
  for(int mtrack = 0; mtrack < ntracks_reco; mtrack++){

    if((*track_zpos)[mtrack][0] < 4){
      double delX = wc_x - (*track_xpos)[mtrack][0];
      double delY = wc_y - (*track_ypos)[mtrack][0];
      double delTheta = wc_theta - UtilityFunctions::getTrackTheta(mtrack, track_xpos,track_ypos,track_zpos);
      double delPhi = wc_phi - UtilityFunctions::getTrackPhi(mtrack, track_xpos,track_ypos);

      std::vector<double> trackDeltas = {delX, delY, delTheta, delPhi};
      vOut.push_back(trackDeltas);
      }
  }
  return vOut;
}





//#############################################################################
//
// END 
//
//#############################################################################
