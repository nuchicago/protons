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






//#############################################################################
//
// END 
//
//#############################################################################
