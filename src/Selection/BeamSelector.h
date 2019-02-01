#ifndef BEAMSELECTOR_H
#define BEAMSELECTOR_H
#include <vector>
#include <iostream>
#include "../Utilities/UtilityFunctions.h"


class BeamSelector {

 public:

  BeamSelector();

  // counters for command line summary
  double numEventsStart = 0;
  double numMassCut = 0;
  double numtofvalid = 0;
  double numZcutoff = 0;
  double numWCTrack = 0;
  double numXYdeltaCut = 0;
  double numAlphaCut = 0;
  double numPileupCut  = 0;
  double numShowerCut = 0;
  double numUniqueMatch = 0;

  // Beam Centering Parameters
  double xMeanTPCentry = 0;
  double yMeanTPCentry = 0;

  // Beam Matching containers
  std::vector<int> EnteringTrkStart;
  std::vector<int> EnteringTrkSecond;
  std::vector<int> EnteringTrkID;
  std::vector<int> EnteringTrkEnd;
  std::vector<double> EnteringTrkAlpha;
  std::vector<double> EnteringTrkDelX;
  std::vector<double> EnteringTrkDelY;
  std::vector<double> EnteringTrkRdist;
  std::vector<double> EnteringTrkTheta;
  std::vector<double> EnteringTrkNhits;
  std::vector<double> EnteringTrkLength;


  std::vector<double> AllTrkZmax;
  std::vector<double> AllTrkZmin2;
  std::vector<double> AllTrkTheta;
  std::vector<double> AllTrkDelX;
  std::vector<double> AllTrkDelY;
  std::vector<double> AllTrkZmin1;
  std::vector<double> AllTrkAlpha;
  std::vector<double> AllTrkRdist;
  std::vector<double> AllTrkNhits;
  std::vector<double> AllTrkLength;
  


  int ShowerTracksBuffer;
  int PileupTracksBuffer;

  void SetMeanXY(double xMean, double yMean);

  std::vector<double> BeamCentering(double wc_x, double wc_y,// double wc_theta, double wc_phi,
                  std::vector< std::vector<double> > *track_xpos,
                  std::vector< std::vector<double> > *track_ypos,
                  std::vector< std::vector<double> > *track_zpos, int ntracks_reco, std::vector<int> *ntrack_hits,
                   double zPointCutoff,
                   int& best_candidate);


  std::vector<double> BeamMatching(double wc_x, double wc_y, double wc_theta, double wc_phi,
                  std::vector< std::vector<double> > *track_xpos,
                  std::vector< std::vector<double> > *track_ypos,
                  std::vector< std::vector<double> > *track_zpos, int ntracks_reco, std::vector<int> *ntrack_hits,
                  std::vector<double> *track_length ,int& best_candidate, std::vector<double> BSoptions);

  void printSummary(std::vector <double> BSoptions);
 
  int  isTPCPrimary(std::vector< std::vector<double> > *track_zpos, int ntracks_reco,
  				 bool mc_mode, double zPointCutoff, int& reco_primary, double& first_reco_z, int verbose);
  
  std::vector< std::vector<double> > wcTPCMatchPlots(double wc_x, double wc_y, double wc_theta, double wc_phi,
                  std::vector< std::vector<double> > *track_xpos,
                  std::vector< std::vector<double> > *track_ypos,
                  std::vector< std::vector<double> > *track_zpos, int ntracks_reco, double zPointCutoff,int& numEntering);

  std::vector<double> dataTPCPrimary(double wc_x, double wc_y, double wc_theta, double wc_phi,
                  std::vector< std::vector<double> > *track_xpos,
                  std::vector< std::vector<double> > *track_ypos,
                  std::vector< std::vector<double> > *track_zpos, 
                  int ntracks_reco, double zPointCutoff, int& MatchedTrack, int& numEntering);
  

  double getMCInitialKE(double initial_ke, unsigned int geant_list_size, std::vector<int> *process_primary, 
                      std::vector<int> *NTrTrajPts, std::vector<std::vector<double>> *MidPosX,
                      std::vector<std::vector<double>> *MidPosY, std::vector<std::vector<double>> *MidPosZ, 
                      std::vector<std::vector<double>> *MidPx, std::vector<std::vector<double>> *MidPy,
                      std::vector<std::vector<double>> *MidPz); 

  double getDataInitialKE(double initial_ke, double wctrk_momentum, double ParticleMass);


  bool MassCut(double wctrk_momentum, double tofObject, double tofLength, double tofOffset, double& ParticleMass,
  double MassCutMin, double MassCutMax);


  std::vector <double> backProjections(double wctrk_XFace, double wctrk_YFace, double wctrk_momentum,
                                                    double wctrk_theta, double wctrk_phi);
 protected:

  static const double pi;                  // the beloved constant  
  static const double massProton;          // proton mass GeV
  static const double massPion;            // piplus/minus mass GeV
  static const double massElectron;        // electron mass GeV
  static const double massKaon;            // kplus/kminus mass GeV
  static const double c_light;             // cm/ns - speed of light in vacuum
  //static const double tofLength;            // cm distance between wc

 private:


};

#endif

 
