
#ifndef EVENTSELECTOR_H
#define EVENTSELECTOR_H
#include <vector>
#include <iostream>
#include <cmath>

#include "../Utilities/ROOTinclude.h"


class EventSelector {

 public:

  EventSelector( );
 
  int classifyEvent( int type );

  double* findInt(  double* candidate_array, int reco_primary, Int_t &ntracks_reco, std::vector<int>* ntrack_hits,
                    std::vector<std::vector<double>>* track_xpos, std::vector<std::vector<double>>* track_ypos,
                    std::vector<std::vector<double>>* track_zpos, std::vector<double>* track_end_x,
                    std::vector<double>* track_end_y, std::vector<double>* track_end_z,
                    std::vector<int>* col_track_hits, std::vector<std::vector<double>>* col_track_dedx,
                    std::vector<std::vector<double>>* col_track_pitch_hit,
                    std::vector<std::vector<double>>* col_track_x, std::vector<std::vector<double>>* col_track_y,
                    std::vector<std::vector<double>>* col_track_z );

  //void getSlabInfo( std::vector<double> &calo_slab_xpos, std::vector<double> &calo_slab_ypos, 
  //                  std::vector<double> &calo_slab_ypos, std::vector<double> &calo_slab_zpos);
  
  int getSlabInfo(std::vector<double> &calo_slab_xpos, std::vector<double> &calo_slab_ypos,
                                std::vector<double> &calo_slab_zpos, std::vector<double> &calo_slab_KE,
                                int reco_primary, double z2, double initial_ke,
                                std::vector<int>* col_track_hits, std::vector<std::vector<double>>* col_track_dedx,
                                std::vector<std::vector<double>>* col_track_pitch_hit,
                                std::vector<std::vector<double>>* col_track_x, std::vector<std::vector<double>>* col_track_y,
                                std::vector<std::vector<double>>* col_track_z); 






 private:


};

#endif

 
