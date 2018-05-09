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
  double zPointCutoff, int& reco_primary, double& first_reco_z, int verbose){

      bool print = false;
      if(verbose == 2){print = true;}
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
        // ### adding a data version of this function ###
        // ### right now all it does is make a vector of ###
        // ### track ids for early tracks (first Zpoint < 2cm) ###
        // ### im randomly picking the first one in this vector ###
        // ### later need to do matching on these tracks !!!!!! ###
        std::vector<int> track_id_early_track;
        for(int rtrack = 0; rtrack < ntracks_reco; rtrack++) {
          double z1 = (*track_zpos)[rtrack][0];
          if(z1 < zPointCutoff) {
            track_id_early_track.push_back(rtrack);
          }
        }//<---End loop reco tracks
        // ## returning randomly one of the early tracks ##
        if(track_id_early_track.size()) {
          return track_id_early_track.at(0);
        }
        else {
          return -1;
        }
      }//<-- End if data

}

std::vector<std::vector<double>> BeamSelector::wcTPCMatchPlots(double wc_x, double wc_y, double wc_theta, double wc_phi,
                  std::vector< std::vector<double> > *track_xpos,
                  std::vector< std::vector<double> > *track_ypos,
                  std::vector< std::vector<double> > *track_zpos, int ntracks_reco, double zPointCutoff, int& numEntering){


  std::vector<std::vector<double>> vOut;
  
  for(int mtrack = 0; mtrack < ntracks_reco; mtrack++){

    if((*track_zpos)[mtrack][0] < zPointCutoff){
      numEntering++;
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

std::vector<double> BeamSelector::wcTPCMatch(double wc_x, double wc_y, double wc_theta, double wc_phi,
                  std::vector< std::vector<double> > *track_xpos,
                  std::vector< std::vector<double> > *track_ypos,
                  std::vector< std::vector<double> > *track_zpos, int ntracks_reco, double zPointCutoff, int& MatchedTrack, int& numEntering){

  
  double rValueMin = 999.;
  double delXMin;
  double delYMin;
  double delThetaMin;
  double delPhiMin;
  
  for(int mtrack = 0; mtrack < ntracks_reco; mtrack++){

    if((*track_zpos)[mtrack][0] < zPointCutoff){
      numEntering++;
      double delX = wc_x - (*track_xpos)[mtrack][0];
      double delY = wc_y - (*track_ypos)[mtrack][0];
      double delTheta = wc_theta - UtilityFunctions::getTrackTheta(mtrack, track_xpos,track_ypos,track_zpos);
      double delPhi = wc_phi - UtilityFunctions::getTrackPhi(mtrack, track_xpos,track_ypos);


      double rValue = sqrt(pow(delX,2) + pow(delY,2));
      if(rValue < rValueMin){
        rValueMin = rValue;
        delXMin = delX;
        delYMin = delY;
        MatchedTrack = mtrack;
        delThetaMin = delTheta;
        delPhiMin = delPhi;
      }
    }
  }


  std::vector<double> minVector = {rValueMin,delXMin,delYMin,delThetaMin,delPhiMin};

  return minVector;
}




double BeamSelector::getMCInitialKE(double initial_ke, unsigned int geant_list_size, std::vector<int> *process_primary, 
                                    std::vector<int> *NTrTrajPts, std::vector<std::vector<double>> *MidPosX,
                                    std::vector<std::vector<double>> *MidPosY, std::vector<std::vector<double>> *MidPosZ, 
                                    std::vector<std::vector<double>> *MidPx, std::vector<std::vector<double>> *MidPy,
                                    std::vector<std::vector<double>> *MidPz) {

    double mass = 938.57;
    int first_pt = 0;

    for(unsigned int g4part = 0; g4part < geant_list_size; g4part++) {
      if((*process_primary)[g4part] != 1) {continue;}// skipping non primary g4 ID
      for(int pt = 1; pt < (*NTrTrajPts)[g4part]; pt++){
        double xpos = (*MidPosX)[g4part][pt];
        double ypos = (*MidPosY)[g4part][pt];
        double zpos = (*MidPosZ)[g4part][pt];
    
        double p = sqrt(pow(1000*(*MidPx)[g4part][pt], 2)
                      + pow(1000*(*MidPy)[g4part][pt], 2)
                      + pow(1000*(*MidPz)[g4part][pt], 2));
        double ke = sqrt(pow(mass, 2) + pow(p, 2)) - mass; 
        if( xpos > 0 && xpos < 47.5 && ypos > -20 && ypos < 20 && zpos > 0 && zpos < 90 ) {
          if(first_pt == 0){initial_ke = ke; first_pt++;}
        }//<--End if in tpc
      }//<--End spt loop
    }//<--End g4 loop


    return initial_ke;

}


double BeamSelector::getDataInitialKE(double initial_ke, double wctrk_momentum) {
  
  double mass = 938.57;
  double wc_ke = sqrt(pow(mass, 2) + pow(wctrk_momentum, 2)) - mass;
  initial_ke = wc_ke;


    return initial_ke; 
}


//#############################################################################
//
// END 
//
//#############################################################################
