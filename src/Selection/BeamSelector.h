
#ifndef BEAMSELECTOR_H
#define BEAMSELECTOR_H
#include <vector>
#include <iostream>
#include "../Utilities/UtilityFunctions.h"


class BeamSelector {

 public:

  BeamSelector();
 
  bool PrimaryTrack(std::vector< std::vector<double> > *track_zpos, int ntracks_reco,
  					 double zPointCutoff, int& reco_primary, double& first_reco_z);

  int  isTPCPrimary(std::vector< std::vector<double> > *track_zpos, int ntracks_reco,
  				 bool mc_mode, double zPointCutoff, int& reco_primary, double& first_reco_z, int verbose);
  
  std::vector< std::vector<double> > wcTPCMatchPlots(double wc_x, double wc_y, double wc_theta, double wc_phi,
                  std::vector< std::vector<double> > *track_xpos,
                  std::vector< std::vector<double> > *track_ypos,
                  std::vector< std::vector<double> > *track_zpos, int ntracks_reco, double zPointCutoff,int& numEntering);

  std::vector<double> wcTPCMatch(double wc_x, double wc_y, double wc_theta, double wc_phi,
                  std::vector< std::vector<double> > *track_xpos,
                  std::vector< std::vector<double> > *track_ypos,
                  std::vector< std::vector<double> > *track_zpos, 
                  int ntracks_reco, double zPointCutoff, int& MatchedTrack, int& numEntering);
  

  double getMCInitialKE(double initial_ke, unsigned int geant_list_size, std::vector<int> *process_primary, 
                      std::vector<int> *NTrTrajPts, std::vector<std::vector<double>> *MidPosX,
                      std::vector<std::vector<double>> *MidPosY, std::vector<std::vector<double>> *MidPosZ, 
                      std::vector<std::vector<double>> *MidPx, std::vector<std::vector<double>> *MidPy,
                      std::vector<std::vector<double>> *MidPz); 

  double getDataInitialKE(double initial_ke, double wctrk_momentum);

 private:


};

#endif

 
