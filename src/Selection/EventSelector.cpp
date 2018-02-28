//#############################################################################
//
// EventSelector.cpp
//
//  A desription of what this class does
//
// January 17, 2018
//
//#############################################################################


#include "EventSelector.h"


//=============================================================================
// Constructor
//=============================================================================

EventSelector::EventSelector(){  }


//=============================================================================
// classsifyEvent()
//=============================================================================
int EventSelector::classifyEvent( int type ){

  return type;

}

bool EventSelector::findInt(int reco_primary, Int_t &ntracks_reco, std::vector<int>* ntrack_hits, 
                            std::vector<std::vector<double>>* track_xpos, std::vector<std::vector<double>>* track_ypos,
                            std::vector<std::vector<double>>* track_zpos, std::vector<double>* track_end_x,
                            std::vector<double>* track_end_y, std::vector<double>* track_end_z,
                            std::vector<int>* col_track_hits, std::vector<std::vector<double>>* col_track_dedx,
                            std::vector<std::vector<double>>* col_track_pitch_hit,
                            std::vector<std::vector<double>>* col_track_x, std::vector<std::vector<double>>* col_track_y,
                            std::vector<std::vector<double>>* col_track_z) {


  //std::cout<<"\t\tFINDINT: \n";
  //std::cout<<"\t\t\treco_primary: "<<reco_primary<<std::endl;
  //std::cout<<"\t\t\tntracks_reco: "<<ntracks_reco<<std::endl;
  //std::cout<<"\t\t\tntrack_hits[reco_primary]: "<<(*ntrack_hits)[reco_primary]<<std::endl;
  //std::cout<<"\t\t\tfirst trackxpos: "<<(*track_xpos)[reco_primary][0]<<std::endl;
  //std::cout<<"\t\t\tfirst trackypos: "<<(*track_ypos)[reco_primary][0]<<std::endl;
  //std::cout<<"\t\t\tfirst trackzpos: "<<(*track_zpos)[reco_primary][0]<<std::endl;
  //std::cout<<"\t\t\ttrack_end_x: "<<(*track_end_x)[reco_primary]<<std::endl;
  //std::cout<<"\t\t\ttrack_end_y: "<<(*track_end_y)[reco_primary]<<std::endl;
  //std::cout<<"\t\t\ttrack_end_z: "<<(*track_end_z)[reco_primary]<<std::endl;
  //std::cout<<"\t\t\tcol track hits: "<<(*col_track_hits)[reco_primary]<<std::endl;
  //std::cout<<"\t\t\tfirst hit dedx: "<<(*col_track_dedx)[reco_primary][0]<<std::endl;
  //std::cout<<"\t\t\tfirst hit pitch: "<<(*col_track_pitch_hit)[reco_primary][0]<<std::endl;
  //std::cout<<"\t\t\tfirst hit x: "<<(*col_track_x)[reco_primary][0]<<std::endl;
  //return true;




    // ### Inelastic Event Selection ###
    std::vector<std::vector<double>> branches;
    std::vector<std::vector<double>> kinks;
    int missing_bragg = 0;
    double prim_endx = (*track_end_x)[reco_primary];
    double prim_endy = (*track_end_y)[reco_primary];
    double prim_endz = (*track_end_z)[reco_primary];
    // # track loop #
    for(int rtrack = 0; rtrack < ntracks_reco; rtrack++){
      if(rtrack == reco_primary){
        // ## primary spacepoint loop ##
        for(int rspt = 1; rspt < (*ntrack_hits)[rtrack]-1; rspt++){
          double prim_x = (*track_xpos)[rtrack][rspt];
          double prim_y = (*track_ypos)[rtrack][rspt];
          double prim_z = (*track_zpos)[rtrack][rspt];
          double prim_prev_x = (*track_xpos)[rtrack][rspt-1];
          double prim_prev_y = (*track_ypos)[rtrack][rspt-1];
          double prim_prev_z = (*track_zpos)[rtrack][rspt-1];
          double prim_next_x = (*track_xpos)[rtrack][rspt+1];
          double prim_next_y = (*track_ypos)[rtrack][rspt+1];
          double prim_next_z = (*track_zpos)[rtrack][rspt+1];

          // # calculating angle along the track #
          double this_pt [] = {prim_x, prim_y, prim_z};
          double last_pt [] = {prim_prev_x, prim_prev_y, prim_prev_z};
          double next_pt [] = {prim_next_x, prim_next_y, prim_next_z};
          double a_vec [] = {prim_x - prim_prev_x, prim_y - prim_prev_y, prim_z - prim_prev_z};
          double b_vec [] = {prim_next_x - prim_x, prim_next_y - prim_y, prim_next_z - prim_z};
          double a_dot_b = (a_vec[0]*b_vec[0]) + (a_vec[1]*b_vec[1]) + (a_vec[2]*b_vec[2]);
          double mag_a = sqrt( pow(a_vec[0],2) + pow(a_vec[1],2) + pow(a_vec[2],2));
          double mag_b = sqrt( pow(b_vec[0],2) + pow(b_vec[1],2) + pow(b_vec[2],2));
          double denom = mag_a * mag_b;
          double theta = acos(a_dot_b / (mag_a*mag_b));
          if(TMath::IsNaN(theta)){theta = 0.;}
          //RDSptAngle->Fill(theta);
          if(theta>4){
            std::vector<double> kink_tuple = {1.*rspt, theta};
            kinks.push_back(kink_tuple);
          }//<--End if this is a kink
          //std::cout<<"\tp-1: ("<<prim_prev_x<<", "<<prim_prev_y<<", "<<prim_prev_z<<")\n";
          //std::cout<<"\tp:   ("<<prim_x<<", "<<prim_y<<", "<<prim_z<<")\n";
          //std::cout<<"\tp+1: ("<<prim_next_x<<", "<<prim_next_y<<", "<<prim_next_z<<")\n";
          //std::cout<<"\t\ta_vec: ("<<a_vec[0]<<", "<<a_vec[1]<<", "<<a_vec[2]<<")\n";
          //std::cout<<"\t\tb_vec: ("<<b_vec[0]<<", "<<b_vec[1]<<", "<<b_vec[2]<<")\n";
          //std::cout<<"\t\t\ta.b: "<<a_dot_b<<std::endl;
          //std::cout<<"\t\t\tmag(a): "<<mag_a<<std::endl;
          //std::cout<<"\t\t\tmag(b): "<<mag_b<<std::endl;
          //std::cout<<"\t\t\t\tdenom: "<<denom<<std::endl;
          //std::cout<<"\t\t\t\tnum/denom: "<<a_dot_b/denom<<std::endl;
          //std::cout<<"\t\t\t\t\ttheta: "<<theta<<std::endl;
        }//<---End primary reco space point loop

        // ## calo obj loop looking for brag peak? ##
        double end_track_dist = 0;
        double end_track_dedx_sum = 0;
        double end_track_counter = 0;
        for(int calo_pt = (*col_track_hits)[rtrack]; calo_pt > 0; calo_pt--){
          if(end_track_dist > 2.5){continue;}
          end_track_dist += (*col_track_pitch_hit)[rtrack][calo_pt];
          end_track_dedx_sum += (*col_track_dedx)[rtrack][calo_pt];
          end_track_counter++;
        }
        double end_track_dedx_mean = end_track_dedx_sum / end_track_counter;
        if(end_track_dedx_mean < 13){ // need to make this number a variable at some point
          missing_bragg++;
        }//<--End if this track had no bragg peak

      }//<---End if primary  

      if(rtrack != reco_primary){
        //~~~~~~~~~~~~std::cout<<"\n\n\nother tracks!!!!\n";
        // ## for each of these other tracks need to compare the start of this track ##
        // ## to every spacepoint in the primary ##
        // ## and keep record of the closest spacepoint(s) from the primary ##
        // ## (plural if multiple other tracks) ##
        
        // ## start position ##
        double start_other_x = (*track_xpos)[rtrack][0];
        double start_other_y = (*track_ypos)[rtrack][0];
        double start_other_z = (*track_zpos)[rtrack][0];
        
        // ## pushing every branching track to a 2d vector ##
        // ## (track#, closest primary spt, distance to that spt) ##
        double min_dist_prim_spt = 99;
        double closest_prim_spt = -1;

        for(int prim_spt = 0; prim_spt < (*ntrack_hits)[reco_primary]; prim_spt++){
          double prim_xpos = (*track_xpos)[reco_primary][prim_spt];
          double prim_ypos = (*track_ypos)[reco_primary][prim_spt];
          double prim_zpos = (*track_zpos)[reco_primary][prim_spt];
          double dist_prim_spt = sqrt(pow(start_other_x - prim_xpos, 2) + 
                                      pow(start_other_y - prim_ypos, 2) + 
                                      pow(start_other_z - prim_zpos, 2));  
          if(dist_prim_spt < min_dist_prim_spt){
            min_dist_prim_spt = dist_prim_spt;
            closest_prim_spt = prim_spt;
          }
        }//<-End primary spt loop

        std::vector<double> branch_tuple = {1.*rtrack, closest_prim_spt, min_dist_prim_spt};
        branches.push_back(branch_tuple);
      }//<--End if not primary
    }//<---End reco track loop for finding interaction candidates


    // ### Topological Sorting and Identifying Potential Interactions ###
    // ### (still a part of inelastic event selection, just putting the pieces together) ###
    bool interacting_candidate;
    int candidate_spt = 999;
    double candidate_xpos = 999;
    double candidate_ypos = 999;
    double candidate_zpos = 999;
    std::vector<int> potential_interaction_pts;
    // # topology 3, 4
    if(branches.size()){
      int earliest_branch_spt = 9999;
      for(int branch_pt = 0; branch_pt < branches.size(); branch_pt++){
        if(branches[branch_pt][1] < earliest_branch_spt){
          earliest_branch_spt = branches[branch_pt][1];
        }
      }
      potential_interaction_pts.push_back((int)earliest_branch_spt);
    }//<--End if there were branches
    // # topology 1
    if(kinks.size()){
      int earliest_kink_spt = 9999;
      for(int kink_pt = 0; kink_pt < kinks.size(); kink_pt++){
        if(kinks[kink_pt][0] < earliest_kink_spt){
          earliest_kink_spt = kinks[kink_pt][0];
        }
      }
      potential_interaction_pts.push_back((int)earliest_kink_spt);
    }//<--End if there was a kink
    // # getting earliest of topology 1,3,4
    for(int i = 0; i < potential_interaction_pts.size(); i++){
      interacting_candidate = true;
      if(potential_interaction_pts[i] < candidate_spt){
        candidate_spt = potential_interaction_pts[i];
        candidate_xpos = (*track_xpos)[reco_primary][potential_interaction_pts[i]];
        candidate_ypos = (*track_ypos)[reco_primary][potential_interaction_pts[i]];
        candidate_zpos = (*track_zpos)[reco_primary][potential_interaction_pts[i]];
      }
    }
    // # topology 2
    if(!(branches.size() || kinks.size())){
      if(missing_bragg && prim_endz < 88){
        //nTopology2++;
        interacting_candidate = true;
        candidate_spt = (*ntrack_hits)[reco_primary] - 1;
        candidate_xpos = (*col_track_x)[reco_primary][(*col_track_hits)[reco_primary]-1];
        candidate_ypos = (*col_track_y)[reco_primary][(*col_track_hits)[reco_primary]-1];
        candidate_zpos = (*col_track_z)[reco_primary][(*col_track_hits)[reco_primary]-1];
      }
      else{
        interacting_candidate = false;
      }
    }//<-End if no branches or kinks

    // ## gonna need to change this an array at some point
    return interacting_candidate;



}


//#############################################################################
//
// END 
//
//#############################################################################
