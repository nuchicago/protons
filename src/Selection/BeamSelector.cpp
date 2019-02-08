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

const double BeamSelector::pi = 3.14159;                // the beloved constant
const double BeamSelector::massProton = 0.938;          // proton mass GeV
const double BeamSelector::massPion = 0.140;            // piplus/minus mass GeV
const double BeamSelector::massElectron = 0.000511;     // electron mass GeV
const double BeamSelector::massKaon = 0.494;            // kplus/kminus mass GeV
const double BeamSelector::c_light = 29.9792458;        // cm/ns - speed of light in vacuum
//const double BeamSelector::tofLength = 665.2;         // cm


//=============================================================================
// Constructor
//=============================================================================

BeamSelector::BeamSelector(){
  
}


//=============================================================================
// classsifyEvent()
//=============================================================================




std::vector<double> BeamSelector::BeamCentering(double wc_x, double wc_y,// double wc_theta, double wc_phi,
                  std::vector< std::vector<double> > *track_xpos,
                  std::vector< std::vector<double> > *track_ypos,
                  std::vector< std::vector<double> > *track_zpos, int ntracks_reco, std::vector<int> *ntrack_hits,
                   double zPointCutoff,
                   int& best_candidate){
          //////////////////////////////////////////////////////////
          //  finds closest track to Wire chamber projection
          //
          //
          //
          //
          //
          //  Returns: vector of doubles made up of
          //      {trackID, 1st Z point index, 2nd Z point index, Last Z point Index,
          //      WC to tpc track start distance in XY, angle between WC track and TPC track}
          /////////////////////////////////////////////////////////


        std::vector<double> centeringMatch;
        double rValueMin = 999.;

        for (int itrack = 0; itrack < ntracks_reco ;  itrack++ ){
          //std::cout << "track : " << itrack  <<"ntrack_hits" << (*ntrack_hits)[itrack]<< std::endl;

          std::vector<int> zIndices = UtilityFunctions::zOrderedTrack2(track_zpos,itrack, ntrack_hits);

          //std::cout << "loading zmin1 : "<< zIndices[0 ]<< std::endl;
          double zmin1 = (*track_zpos)[itrack][zIndices[0]];
          //std::cout << "loading zmin2 : "<< zIndices[1] << std::endl;
          double zmin2 = (*track_zpos)[itrack][zIndices[1]];
          //std::cout << "loading zlast : "<< zIndices[2] << std::endl;
          double zmax = (*track_zpos)[itrack][zIndices[2]];

          
          
          double delXMin;
          double delYMin;
          double alphaMin;
        
          if(zmin1 < zPointCutoff){

             //std::cout << "Calculating delX : " << itrack << std::endl;
            double delX = wc_x - (*track_xpos)[itrack][zIndices[0]];
            double delY = wc_y - (*track_ypos)[itrack][zIndices[0]];
            //double delTheta = wc_theta - UtilityFunctions::getTrackTheta(mtrack, track_xpos,track_ypos,track_zpos);
            //double delPhi;   ########reminder to reimplement angle.


            double alpha = -999.;
            double rValue = sqrt(pow(delX,2) + pow(delY,2));

            

            if(rValue < rValueMin){
                rValueMin = rValue;
                centeringMatch = {1.* itrack, 1.* zIndices[0], 1.* zIndices[1], 1.* zIndices[2], rValue, alpha};
                best_candidate =  itrack;

              }
            }
        //std::cout << "Moving to next tracl \n " << std::endl;
      }
      return centeringMatch;
    }

void BeamSelector::SetMeanXY(double xMean, double yMean){
      xMeanTPCentry = xMean;
      yMeanTPCentry = yMean;
    }

std::vector<double> BeamSelector::BeamMatching(double wc_x, double wc_y, double wc_theta, double wc_phi,
                  std::vector< std::vector<double> > *track_xpos,
                  std::vector< std::vector<double> > *track_ypos,
                  std::vector< std::vector<double> > *track_zpos, int ntracks_reco, std::vector<int> *ntrack_hits,
                   std::vector<double> *track_length, int& best_candidate, std::vector<double> BSoptions){
          //////////////////////////////////////////////////////////
          //  Matches Wire Chamber to TPC and Applies Cuts
          //    to print summary of cuts, use BeamSelector::printSummary
          //
          //
          //
          //
          //  Returns: vector of doubles made up of
          //      {isMatch(1 or 0), trackID, 1st Z point index, 2nd Z point index, Last Z point Index,
          //      ,delX, delY, angle between WC track and TPC track, uniqueMatch(1 or 0)}
          /////////////////////////////////////////////////////////
        //if(BSoptions[0]){std::cout << "Launching Beam Matching" << std::endl;}

        numEventsStart++;
        std::vector<double> matchInfo {0,-1, -1, -1, -1,-1, -1, -1,-1};
        int numMatchesFound = 0;
        int numShortTracks = 0;
        int numPileupTracks = 0;
        int foundMatch = 0;
        double rValueMin = 999.;

        bool zConditionMet = false;
        bool angleConditionMet = false;
        bool circleConditionMet = false;
        bool pileupConditionMet = false;
        bool showerConditionMet = false;

        EnteringTrkStart.clear();
        EnteringTrkSecond.clear();
        EnteringTrkID.clear();
        EnteringTrkEnd.clear();
        EnteringTrkAlpha.clear();
        EnteringTrkRdist.clear();
        EnteringTrkDelX.clear();
        EnteringTrkDelY.clear();
        EnteringTrkNhits.clear();
        EnteringTrkLength.clear();
        
        AllTrkZmin1.clear();
        AllTrkZmin2.clear();
        AllTrkZmax.clear();
        AllTrkAlpha.clear();
        AllTrkRdist.clear();
        AllTrkLength.clear();
        AllTrkNhits.clear();
        AllTrkDelX.clear();
        AllTrkDelY.clear();
        AllTrkTheta.clear();

        AllTrkAlphaEnd.clear();
        AllTrkDirCompare.clear();

        EnteringTrkAlphaEnd.clear();
        EnteringTrkDirCompare.clear();

        CircleTrkAlpha.clear();
        



        for (int itrack = 0; itrack < ntracks_reco ;  itrack++ ){

          //if(BSoptions[0]){std::cout << "ordering Z points" << std::endl;}

          std::vector<int> zIndices = UtilityFunctions::zOrderedTrack2(track_zpos,itrack, ntrack_hits);

          double zmin1 = (*track_zpos)[itrack][zIndices[0]];
          double zmin2 = (*track_zpos)[itrack][zIndices[1]];
          double zmax = (*track_zpos)[itrack][zIndices[2]];

          //if(BSoptions[0]){std::cout << "finished ordering" << std::endl;}

          
          double delXMin;
          double delYMin;
          double alphaMin;

          double trkLength = (*track_length)[itrack];

          if(trkLength < BSoptions[8]){ numShortTracks++;}
          if(zmin1 < BSoptions[6]){numPileupTracks++;}

          if(zmin1 < BSoptions[1]){
            if(!zConditionMet){numZcutoff++;}
            zConditionMet = true;

            EnteringTrkStart.push_back(zIndices[0]);
            EnteringTrkSecond.push_back(zIndices[1]);
            EnteringTrkID.push_back(itrack);
            EnteringTrkEnd.push_back(zIndices[2]);
            EnteringTrkLength.push_back((*track_length)[itrack]);
            EnteringTrkNhits.push_back((*ntrack_hits)[itrack]);


          }


            AllTrkZmin1.push_back(zmin1);
            AllTrkZmin2.push_back(zmin2);
            AllTrkZmax.push_back(zmax);
            AllTrkLength.push_back((*track_length)[itrack]);
            AllTrkNhits.push_back((*ntrack_hits)[itrack]);

            //if(BSoptions[0]){std::cout << "calculating delX" << std::endl;}

             //std::cout << "Calculating delX : " << itrack << std::endl;
            double delX = wc_x - (*track_xpos)[itrack][zIndices[0]];
            double delY = wc_y - (*track_ypos)[itrack][zIndices[0]];

            

              double xDistTheta = (*track_xpos)[itrack][zIndices[1]]-(*track_xpos)[itrack][zIndices[0]];
              double yDistTheta = (*track_ypos)[itrack][zIndices[1]]-(*track_ypos)[itrack][zIndices[0]];
              double zDistTheta = (*track_zpos)[itrack][zIndices[1]]-(*track_zpos)[itrack][zIndices[0]];
              double d = sqrt(pow(xDistTheta,2)+pow(yDistTheta,2));
              double Theta = atan(d/zDistTheta) * (180 / pi);


            double tpc_vec []= { ((*track_xpos)[itrack][zIndices[1]] -(*track_xpos)[itrack][zIndices[0]]),
                              ((*track_ypos)[itrack][zIndices[1]] -(*track_ypos)[itrack][zIndices[0]]),
                              ((*track_zpos)[itrack][zIndices[1]] -(*track_zpos)[itrack][zIndices[0]])};


            double tpc_vec_end []= { ((*track_xpos)[itrack][zIndices[2]] -(*track_xpos)[itrack][zIndices[0]]),
                              ((*track_ypos)[itrack][zIndices[2]] -(*track_ypos)[itrack][zIndices[0]]),
                              ((*track_zpos)[itrack][zIndices[2]] -(*track_zpos)[itrack][zIndices[0]])};
            
            double tpc_size = UtilityFunctions::pointDistance((*track_xpos)[itrack][zIndices[1]],
                                                              (*track_ypos)[itrack][zIndices[1]],
                                                              (*track_zpos)[itrack][zIndices[1]], 
                                                              (*track_xpos)[itrack][zIndices[0]],
                                                              (*track_ypos)[itrack][zIndices[0]],
                                                              (*track_zpos)[itrack][zIndices[0]]);

            double tpc_end_size = UtilityFunctions::pointDistance((*track_xpos)[itrack][zIndices[2]],
                                                              (*track_ypos)[itrack][zIndices[2]],
                                                              (*track_zpos)[itrack][zIndices[2]], 
                                                              (*track_xpos)[itrack][zIndices[0]],
                                                              (*track_ypos)[itrack][zIndices[0]],
                                                              (*track_zpos)[itrack][zIndices[0]]);

            double wc_vec []= { sin(wc_theta)*cos(wc_phi), sin(wc_theta)* sin(wc_phi), cos(wc_theta)};

            double wc_dot_tpc = wc_vec[0]*tpc_vec[0]+wc_vec[1]*tpc_vec[1]+wc_vec[2]*tpc_vec[2];

            double wc_dot_tpc_end = wc_vec[0]*tpc_vec_end[0]+wc_vec[1]*tpc_vec_end[1]+wc_vec[2]*tpc_vec_end[2];

            double tpc_dot_tpc_end = tpc_vec[0]*tpc_vec_end[0]+tpc_vec[1]*tpc_vec_end[1]+tpc_vec[2]*tpc_vec_end[2];

            double alpha =  acos(wc_dot_tpc / (tpc_size)) * (180/pi);  
            double alpha_end = acos(wc_dot_tpc_end / (tpc_end_size)) * (180/pi);
            double dir_vector_angle = acos(tpc_dot_tpc_end / (tpc_end_size * tpc_size)) * (180/pi);

            AllTrkAlphaEnd.push_back(alpha_end);
            AllTrkDirCompare.push_back(dir_vector_angle);
            AllTrkAlpha.push_back(alpha);
            AllTrkDelX.push_back(delX);
            AllTrkDelY.push_back(delY);
            AllTrkTheta.push_back(Theta);

          //if(BSoptions[0]){std::cout << "calculated angles" << std::endl;}



            double rValue = sqrt(pow(delX,2) + pow(delY,2));

            AllTrkRdist.push_back(rValue);

            if(zConditionMet){
            EnteringTrkAlphaEnd.push_back(alpha_end);
            EnteringTrkDirCompare.push_back(dir_vector_angle);
            EnteringTrkAlpha.push_back(alpha);
            EnteringTrkDelX.push_back(delX);
            EnteringTrkDelY.push_back(delY);
            EnteringTrkRdist.push_back(rValue);
            EnteringTrkTheta.push_back(Theta);
            }

            double adjustedR = (sqrt(pow((delX - xMeanTPCentry),2) + pow((delY - yMeanTPCentry),2)));

            if(adjustedR < BSoptions[3]){
              if(!circleConditionMet){numXYdeltaCut++;}
              circleConditionMet = true;
              CircleTrkAlpha.push_back(alpha);

              if(alpha < BSoptions[2]){
                if(!angleConditionMet){numAlphaCut++;}
                angleConditionMet = true;
                

                if (rValue < rValueMin){

                  matchInfo = {1., static_cast <double> (itrack), static_cast <double> (zIndices[0]),
                    static_cast <double> (zIndices[1]),static_cast <double> (zIndices[2]), delX, delY,alpha,-1.};
                  rValueMin = rValue;
                  alphaMin = alpha;
                  delXMin = delX;
                  delYMin = delY;
                  best_candidate = itrack;
                  numMatchesFound++;
                }
              }// end if circle cut
            }//end if alpha cut
          //end if ztpc cut
        
      }//end of track loop


      //if(BSoptions[0]){std::cout << "### Done with track loop\n" << std::endl;}

      matchInfo[8] = static_cast<double>(numMatchesFound);
      if(matchInfo[0]){
        if(BSoptions[9]){
          if(numMatchesFound > 1){matchInfo[0] = 0;}
          else{numUniqueMatch++;}
        }
      }

      PileupTracksBuffer = numPileupTracks;
      ShowerTracksBuffer = numShortTracks;
      if(matchInfo[0]){
        if(BSoptions[4]){
          if(numPileupTracks > BSoptions[5]){matchInfo[0] = 0;}
          else{

            numPileupCut++;
            if(numShortTracks > BSoptions[7]){matchInfo[0] = 0;}
            else{
              numShowerCut++;
                
             } 
          }
        }
      }
      return matchInfo;
    }

void BeamSelector::printSummary(std::vector <double> BSoptions){

  std::cout << "\n------- Beam Selection Results -------\n"<< std::endl;

  std::cout << "Events Passed through Matching: "<< numEventsStart << std::endl; 
  std::cout << "Events with at least one TPC track Z < " << BSoptions[1] << ": "<< numZcutoff << std::endl;
  std::cout << "Events passing circular distance cut r < " << BSoptions[3] << ": "<< numXYdeltaCut << std::endl;
  std::cout << "Events passing alpha angle cut a < "  << BSoptions[2] << ": " << numAlphaCut << std::endl;
  if (BSoptions[9]){std::cout << "Events with unique matched track: " << numUniqueMatch << std::endl;}
  if(BSoptions[4]){
  std::cout << "\n------- Additional Data Selection-------\n"<< std::endl;
  std::cout << "Events passing pileup cut ( > " << BSoptions[5] << " tracks in first "<< BSoptions[6]<<" cm): " << numPileupCut << std::endl;
  std::cout << "Events passing EM shower cut ( "  << BSoptions[7] <<" or more tracks < " << BSoptions[8] <<" cm long): " << numShowerCut << std::endl;
  
  }
}


int BeamSelector::isTPCPrimary(std::vector<std::vector<double>> *track_zpos,int ntracks_reco,bool mc_mode, 
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
                  std::vector< std::vector<double> > *track_zpos, int ntracks_reco,
                  double zPointCutoff, int& numEntering){


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

std::vector<double> BeamSelector::dataTPCPrimary(double wc_x, double wc_y, double wc_theta, double wc_phi,
                  std::vector< std::vector<double> > *track_xpos,
                  std::vector< std::vector<double> > *track_ypos,
                  std::vector< std::vector<double> > *track_zpos, int ntracks_reco, double zPointCutoff,
                   int& MatchedTrack, int& numEntering){

  
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
      double delPhi;
      if (wc_phi < 0){
        delPhi = (wc_phi + 8* atan(1)) - UtilityFunctions::getTrackPhi(mtrack, track_xpos,track_ypos);
      }
      else{
        delPhi = wc_phi  - UtilityFunctions::getTrackPhi(mtrack, track_xpos,track_ypos);

      }


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


double BeamSelector::getDataInitialKE(double initial_ke, double wctrk_momentum, double ParticleMass) {
  
  //double mass = 938.57;
  double wc_ke = sqrt(pow(ParticleMass, 2) + pow(wctrk_momentum, 2)) - ParticleMass;
  double ke_loss = 60;
  initial_ke = wc_ke - ke_loss;


    return initial_ke; 
}

bool BeamSelector::MassCut(double wctrk_momentum, double tofObject, double tofLength, double tofOffset, double& ParticleMass,
  double MassCutMin, double MassCutMax){
  
  //ParticleMass  = wctrk_momentum * sqrt(abs((pow(tofObject * c_light,2))/pow( tofLength ,2) - 1));


// ======================================================================
// Calculate particle mass [MeV/c^2] from beamline momentum and TOF
// (TOF offset calculated externally from best fit to proton band)

  double adjustedToF = tofObject + tofOffset;
  double radical = (pow(adjustedToF,2)*pow(c_light,2))/(pow(tofLength,2)) - 1;

  //if( radical >= 0. ) {ParticleMass = wctrk_momentum * sqrt( radical );}
  //else{ParticleMass = -wctrk_momentum * sqrt( -radical );}

  ParticleMass = wctrk_momentum * sqrt( abs(radical));


  if (ParticleMass < MassCutMax && ParticleMass > MassCutMin){return true;}
  else{return false;}

}


std::vector <double> BeamSelector::backProjections(double wctrk_XFace, double wctrk_YFace, double wctrk_momentum,
                                                    double wctrk_theta, double wctrk_phi){

    //////////////////////////////////////////////////////////
  //  Calculates key beamline variables for data-driven beam Monte Carlo
  //  projecting back X Y at the face of tpc to wc 4 (z = -100)
  //
  //
  //
  //
  //
  //  Returns: vector of doubles made up of
  //      {Xproj_wc4, Yproj_wc4, Z_proj_wc4, theta_xz, theta yz, momentum}
  /////////////////////////////////////////////////////////

  double phi_2pi;

  if (wctrk_phi < 0 ){
    phi_2pi = 2 * pi + wctrk_phi;
  } 
  else{phi_2pi = wctrk_phi;}

  double theta_xz = asin(sin(wctrk_theta) * cos(phi_2pi)) * (180/pi);
  double theta_yz = asin(sin(wctrk_theta) * sin(phi_2pi)) * (180/pi);
  double Yproj_wc4 =  (sin(wctrk_theta) * sin(phi_2pi)) / (cos(wctrk_theta)) * (-100) + wctrk_YFace;
  double Xproj_wc4 =  (sin(wctrk_theta) * cos(phi_2pi)) / (cos(wctrk_theta)) * (-100) + wctrk_XFace;
  double Z_proj_wc4 = -100;

  std::vector<double> output  = { Xproj_wc4, Yproj_wc4, Z_proj_wc4, theta_xz, theta_yz, wctrk_momentum};

  return output;





}



//#############################################################################
//
// END 
//
//#############################################################################
