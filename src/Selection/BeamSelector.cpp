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
bool BeamSelector::isProton(std::vector<std::vector<double>> *track_zpos,int ntracks_reco,  bool mc_mode, double zPointCutoff){

      bool print = true;


      if(mc_mode){

        if(print){
          std::cout<<"MC Beam Selection"<<std::endl;
          std::cout<<">>>>>>>>>>>>>>>>reco<<<<<<<<<<<<<<<<"<<std::endl;
        }
        int reco_primary = -1;
        double first_reco_z = 99.;
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
      else{
        std::cout<<"Beam Selection for data not available\n"<<std::endl;
      return false;}

}


//#############################################################################
//
// END 
//
//#############################################################################
