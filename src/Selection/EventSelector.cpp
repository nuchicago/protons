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

double* EventSelector::findInt(double* candidate_array, int reco_primary, Int_t &ntracks_reco, std::vector<int>* ntrack_hits, 
                            std::vector<std::vector<double>>* track_xpos, std::vector<std::vector<double>>* track_ypos,
                            std::vector<std::vector<double>>* track_zpos, std::vector<double>* track_end_x,
                            std::vector<double>* track_end_y, std::vector<double>* track_end_z,
                            std::vector<int>* col_track_hits, std::vector<std::vector<double>>* col_track_dedx,
                            std::vector<std::vector<double>>* col_track_pitch_hit,
                            std::vector<std::vector<double>>* col_track_x, std::vector<std::vector<double>>* col_track_y,
                            std::vector<std::vector<double>>* col_track_z,
                            std::vector<double> UIoptions) {


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

    // # return this #
    //double candidate_array[4];

    //unpacking UI options
    double verbose = UIoptions[0];
    double dedxNoBraggMax = UIoptions[1];
    double branchMaxDist =  UIoptions[2];
    double clusterMaxDist = UIoptions[3];

    // emptying containers
    BranchDistVect.clear();
    ClusterDistVect.clear();
    ClusterIDvect.clear();
    // vectors for primary track (in case track is flipped)
    std::vector<double> primary_xpos;
    std::vector<double> primary_ypos;
    std::vector<double> primary_zpos;
    std::vector<double> col_primary_x;
    std::vector<double> col_primary_y;
    std::vector<double> col_primary_z;
    std::vector<double> col_primary_pitch_hit;
    std::vector<double> col_primary_dedx;
    int primary_hits = (*ntrack_hits)[reco_primary];
    int col_primary_hits =  (*col_track_hits)[reco_primary];


    std::vector<int> zIndices = UtilityFunctions::zOrderedTrack2(track_zpos,reco_primary, ntrack_hits);
    if (zIndices[0] != 0){
      for (int i = primary_hits; i > 0; i--){
        primary_xpos.push_back((*track_xpos)[reco_primary][i-1]);
        primary_ypos.push_back((*track_ypos)[reco_primary][i-1]);
        primary_zpos.push_back((*track_zpos)[reco_primary][i-1]);
      }
      for( int i =  col_primary_hits; i > 0 ; i--){
        col_primary_x.push_back((*col_track_x)[reco_primary][i-1]);
        col_primary_y.push_back((*col_track_y)[reco_primary][i-1]);
        col_primary_z.push_back((*col_track_z)[reco_primary][i-1]);
        col_primary_pitch_hit.push_back((*col_track_pitch_hit)[reco_primary][i-1]);
        col_primary_dedx.push_back((*col_track_dedx)[reco_primary][i-1]);
      }
    }
    else{
    primary_xpos =  (*track_xpos)[reco_primary];
    primary_ypos =  (*track_ypos)[reco_primary];
    primary_zpos =  (*track_zpos)[reco_primary];
    col_primary_x = (*col_track_x)[reco_primary];
    col_primary_y = (*col_track_y)[reco_primary];
    col_primary_z = (*col_track_z)[reco_primary];
    col_primary_pitch_hit = (*col_track_pitch_hit)[reco_primary];
    col_primary_dedx = (*col_track_dedx)[reco_primary];
    }


    // ### Inelastic Event Selection ###
    std::vector<std::vector<double>> branches;
    std::vector<std::vector<double>> kinks;
    std::vector<std::vector<double>> kink_angles;
    double no_bragg_dedx_mean;
    int missing_bragg = 0;
    //double prim_endx = (*track_end_x)[reco_primary];
    //double prim_endy = (*track_end_y)[reco_primary];
    double prim_endz = (*track_end_z)[reco_primary];
    // # track loop #
    for(int rtrack = 0; rtrack < ntracks_reco; rtrack++){
      if(rtrack == reco_primary){
        // ## primary spacepoint loop ##
        for(int rspt = 1; rspt < (*ntrack_hits)[rtrack]-1; rspt++){
          double prim_x = primary_xpos[rspt];
          double prim_y = primary_ypos[rspt];
          double prim_z = primary_zpos[rspt];
          double prim_prev_x = primary_xpos[rspt-1];
          double prim_prev_y = primary_ypos[rspt-1];
          double prim_prev_z = primary_zpos[rspt-1];
          double prim_next_x = primary_xpos[rspt+1];
          double prim_next_y = primary_ypos[rspt+1];
          double prim_next_z = primary_zpos[rspt+1];

          // # calculating angle along the track #
          //double this_pt [] = {prim_x, prim_y, prim_z};
          //double last_pt [] = {prim_prev_x, prim_prev_y, prim_prev_z};
          //double next_pt [] = {prim_next_x, prim_next_y, prim_next_z};
          double a_vec [] = {prim_x - prim_prev_x, prim_y - prim_prev_y, prim_z - prim_prev_z};
          double b_vec [] = {prim_next_x - prim_x, prim_next_y - prim_y, prim_next_z - prim_z};
          double a_dot_b = (a_vec[0]*b_vec[0]) + (a_vec[1]*b_vec[1]) + (a_vec[2]*b_vec[2]);
          double mag_a = sqrt( pow(a_vec[0],2) + pow(a_vec[1],2) + pow(a_vec[2],2));
          double mag_b = sqrt( pow(b_vec[0],2) + pow(b_vec[1],2) + pow(b_vec[2],2));
          //double denom = mag_a * mag_b;
          //double theta = acos(a_dot_b / (mag_a*mag_b));
          double theta = acos(a_dot_b / (mag_a*mag_b)) * (180/3.14);
          if(TMath::IsNaN(theta)){theta = 0.;}
          //RDSptAngle->Fill(theta);
          if(theta>6){
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
        for(int calo_pt = col_primary_hits; calo_pt > 0; calo_pt--){
          if(end_track_dist > 2.5){continue;}
          end_track_dist += col_primary_pitch_hit[calo_pt];
          end_track_dedx_sum += col_primary_dedx[calo_pt];
          end_track_counter++;
        }
        double end_track_dedx_mean = end_track_dedx_sum / end_track_counter;
        //dedxNoBraggVector.push_back(end_track_dedx_mean);
        if(end_track_dedx_mean < dedxNoBraggMax){ // need to make this number a variable at some point
          no_bragg_dedx_mean = end_track_dedx_mean;
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

        double end_other_x = (*track_xpos)[rtrack][-1];
        double end_other_y = (*track_ypos)[rtrack][-1];
        double end_other_z = (*track_zpos)[rtrack][-1];
        
        // ## pushing every branching track to a 2d vector ##
        // ## (track#, closest primary spt, distance to that spt) ##


        //double min_dist_prim_spt_start = 99;
        //double min_dist_prim_spt_end = 99;
        //double closest_prim_spt_start = -1;
        //double closest_prim_spt_end = -1;

        
          double prim_xpos = primary_xpos[-1];
          double prim_ypos = primary_ypos[-1];
          double prim_zpos = primary_zpos[-1];
          double dist_prim_spt_start = sqrt(pow(start_other_x - prim_xpos, 2) + 
                                      pow(start_other_y - prim_ypos, 2) + 
                                      pow(start_other_z - prim_zpos, 2));

          //if(dist_prim_spt_start < min_dist_prim_spt_start){
            //min_dist_prim_spt_start = dist_prim_spt_start;
            //closest_prim_spt_start = prim_spt;
          //}
        //<-End primary spt loop - start point

          double dist_prim_spt_end = sqrt(pow(end_other_x - prim_xpos, 2) + 
                                      pow(end_other_y - prim_ypos, 2) + 
                                      pow(end_other_z - prim_zpos, 2));

          //if(dist_prim_spt_end < min_dist_prim_spt_end){
            //min_dist_prim_spt_end = dist_prim_spt_end;
            //closest_prim_spt_end = prim_spt;
          //}
        //<-End primary spt loop - end point



        double closest_prim_spt;
        double min_dist_prim_spt;

        if(dist_prim_spt_start < dist_prim_spt_end){
        min_dist_prim_spt = dist_prim_spt_start;}
        else{
        min_dist_prim_spt = dist_prim_spt_end;
        }



        BranchDistVect.push_back(min_dist_prim_spt);
        if(min_dist_prim_spt < branchMaxDist){
          closest_prim_spt = primary_hits - 1;
          std::vector<double> branch_tuple = {1.*rtrack, closest_prim_spt, min_dist_prim_spt};
          branches.push_back(branch_tuple);}
      }//<--End if not primary
    }//<---End reco track loop for finding interaction candidates


    // ### Topological Sorting and Identifying Potential Interactions ###
    // ### (still a part of inelastic event selection, just putting the pieces together) ###
    int candidate_spt = 999;
    double candidate_xpos = 999;
    double candidate_ypos = 999;
    double candidate_zpos = 999;
    std::vector<int> potential_interaction_pts;
    std::vector<double> potential_interaction_type;
    std::vector<double> potential_interaction_info1;

    // # topology 3, 4
    if(branches.size()){
      int earliest_branch_spt = 9999;
      for(int branch_pt = 0; branch_pt < branches.size(); branch_pt++){
        if(branches[branch_pt][1] < earliest_branch_spt){
          earliest_branch_spt = branches[branch_pt][1];
        }
      }
      potential_interaction_pts.push_back((int)earliest_branch_spt);
      potential_interaction_type.push_back(3.);
      potential_interaction_info1.push_back(1.);
      if(verbose){std::cout<<"earliest branch pt: "<<earliest_branch_spt<<std::endl;}
    }//<--End if there were branches
    // # topology 1
    if(kinks.size()){
      int earliest_kink_spt = 9999;
      double angleValue = -9999;
      for(int kink_pt = 0; kink_pt < kinks.size(); kink_pt++){
        if(kinks[kink_pt][0] < earliest_kink_spt){
          earliest_kink_spt = kinks[kink_pt][0];
          angleValue = kinks[kink_pt][1];
        }
      }
      potential_interaction_pts.push_back((int)earliest_kink_spt);
      potential_interaction_type.push_back(1.);
      potential_interaction_info1.push_back(angleValue);
      if(verbose){std::cout<<"earliest kink pt: "<<earliest_kink_spt<<std::endl;}
    }//<--End if there was a kink
    // # getting earliest of topology 1,3,4
    for(int i = 0; i < potential_interaction_pts.size(); i++){
      candidate_array[0] = 1;
      std::cout<<"candidate_array[0] = 1\n";
      if((*track_zpos)[reco_primary][potential_interaction_pts[i]] < candidate_zpos){
        candidate_spt = potential_interaction_pts[i];
        candidate_xpos = primary_xpos[potential_interaction_pts[i]];
        candidate_ypos = primary_ypos[potential_interaction_pts[i]];
        candidate_zpos = primary_zpos[potential_interaction_pts[i]];
        candidate_array[1] = candidate_xpos;
        candidate_array[2] = candidate_ypos;
        candidate_array[3] = candidate_zpos;
        candidate_array[4] = potential_interaction_type[i];
        //if (potential_interaction_type[i] == 1){}
        candidate_array[5] = potential_interaction_info1[i];
        
        //else{candidate_array[5] = potential_interaction_info1}
      }
    }


    if(candidate_array[0] == 1.){
        double num_branches_t4 = 0.;

        for (int ibranch = 0 ;  ibranch <  ntracks_reco; ibranch++){

          //if(verbose){std::cout << "branches " << ibranch << std::endl;}
          if(ibranch != reco_primary){

            double dist_start  =  sqrt(pow((*track_xpos)[ibranch][0] - candidate_array[1],2)
              + pow((*track_ypos)[ibranch][0] - candidate_array[2],2)
              + pow((*track_zpos)[ibranch][0] - candidate_array[3],2));

            double dist_end  =  sqrt(pow((*track_xpos)[ibranch][(*ntrack_hits)[ibranch] - 1] - candidate_array[1],2)
              + pow((*track_ypos)[ibranch][(*ntrack_hits)[ibranch] - 1] - candidate_array[2],2)
              + pow((*track_zpos)[ibranch][(*ntrack_hits)[ibranch] - 1] - candidate_array[3],2));

            double dist;
            if(dist_start < dist_end){ dist = dist_start;}
            else{dist = dist_end;}
            

            ClusterDistVect.push_back(dist);
            
            if(dist < clusterMaxDist){num_branches_t4 += 1.;
              if(verbose){std::cout << "Branch ID: " << ibranch << std::endl;}
              ClusterIDvect.push_back(ibranch);
            }
          }
        }
        if(verbose){std::cout << "num_branches_type4 = " << num_branches_t4 << std::endl;}
        if (candidate_array[4] == 3.){candidate_array[5] = num_branches_t4;
          if (num_branches_t4 > 1){candidate_array[4] = 4.;}
         }
        //if (num_branches_t4 > 1){
        //  candidate_array[4] = 4.;
        //}
      }

    // # topology 2
    if(!(branches.size() || kinks.size())){
      if(missing_bragg && prim_endz < 88){
        candidate_array[0] = 1;
        candidate_spt = col_primary_hits - 1;
        candidate_xpos = col_primary_x[candidate_spt];
        candidate_ypos = col_primary_y[candidate_spt];
        candidate_zpos = col_primary_z[candidate_spt];
        candidate_array[1] = candidate_xpos;
        candidate_array[2] = candidate_ypos;
        candidate_array[3] = candidate_zpos;
        candidate_array[4] = 2.;
        candidate_array[5] = no_bragg_dedx_mean;
      }
      else{
        candidate_array[0] = 0;
        candidate_array[1] = -1;
        candidate_array[2] = -1;
        candidate_array[3] = -1;
        candidate_array[4] = -1;
        candidate_array[5] = -1;
      }


    }//<-End if no branches or kinks

    /**/



    return candidate_array;

}


int EventSelector::getSlabInfo(std::vector<double> &calo_slab_xpos, std::vector<double> &calo_slab_ypos,
                                std::vector<double> &calo_slab_zpos, std::vector<double> &calo_slab_KE,
                                int reco_primary, double z2, double initial_ke,
                                std::vector<int>* col_track_hits, std::vector<std::vector<double>>* col_track_dedx,
                                std::vector<std::vector<double>>* col_track_pitch_hit,
                                std::vector<std::vector<double>>* col_track_x, std::vector<std::vector<double>>* col_track_y,
                                std::vector<std::vector<double>>* col_track_z) {

    // ### DENOMINATOR ###
    double calo_ke = initial_ke;
    double proj_distance = 0;
    double next_step = 0;
    int ncalo_slab = 1;


    std::vector<double> col_primary_x;
    std::vector<double> col_primary_y;
    std::vector<double> col_primary_z;
    std::vector<double> col_primary_pitch_hit;
    std::vector<double> col_primary_dedx;
    int col_primary_hits =  (*col_track_hits)[reco_primary];


    std::vector<int> zIndices = UtilityFunctions::zOrderedTrack2(col_track_z,reco_primary, col_track_hits);
    if (zIndices[0] != 0){
      
      for( int i =  col_primary_hits; i > 0 ; i--){
        col_primary_x.push_back((*col_track_x)[reco_primary][i-1]);
        col_primary_y.push_back((*col_track_y)[reco_primary][i-1]);
        col_primary_z.push_back((*col_track_z)[reco_primary][i-1]);
        col_primary_pitch_hit.push_back((*col_track_pitch_hit)[reco_primary][i-1]);
        col_primary_dedx.push_back((*col_track_dedx)[reco_primary][i-1]);
      }
    }
    else{
    col_primary_x = (*col_track_x)[reco_primary];
    col_primary_y = (*col_track_y)[reco_primary];
    col_primary_z = (*col_track_z)[reco_primary];
    col_primary_pitch_hit = (*col_track_pitch_hit)[reco_primary];
    col_primary_dedx = (*col_track_dedx)[reco_primary];
    }

    for(int calo_pt = 0; calo_pt < col_primary_hits - 1; calo_pt++){
      proj_distance += col_primary_pitch_hit[calo_pt];
      next_step = proj_distance + col_primary_pitch_hit[calo_pt+1];
      double calo_x = col_primary_x[calo_pt]; 
      double calo_y = col_primary_y[calo_pt]; 
      double calo_z = col_primary_z[calo_pt]; 
      double calo_de = col_primary_dedx[calo_pt]*
                       col_primary_pitch_hit[calo_pt];
      double next_ke = calo_ke - calo_de; 
      double calo_next_x = col_primary_x[calo_pt+1]; 
      double calo_next_y = col_primary_y[calo_pt+1]; 
      double calo_next_z = col_primary_z[calo_pt+1]; 
      if(proj_distance < ncalo_slab*z2 && next_step > ncalo_slab*z2){
        // ## 3d parametrization of 2 calo points ##
        // ## using this to calculate the (x,y,z) of the slab pos ##
        double vx = calo_next_x - calo_x;
        double vy = calo_next_y - calo_y;
        double vz = calo_next_z - calo_z;
        double r = ncalo_slab*z2;
        double A = pow(vx,2) + pow(vy,2) + pow(vz,2);
        double B = 2*(calo_x*vx + calo_y*vy + calo_z*vz);
        double C = pow(proj_distance,2) - pow(r,2);
        double t = (-1*B + pow( pow(B,2) - 4*A*C, .5))/(2*A);
        double calo_slab_x = calo_x + t*vx;
        double calo_slab_y = calo_y + t*vy;
        double calo_slab_z = calo_z + t*vz;
        // ## 1d line between KE values ##
        double calo_ke_slope = (next_ke - calo_ke) / (next_step - proj_distance);
        double calo_slab_ke = calo_ke_slope*(r - proj_distance) + calo_ke;
        //std::cout<<"~~~fillin histos??~~~\n";
        //std::cout<<"number of the calo slab: "<<ncalo_slab<<std::endl;
        //std::cout<<"ncalo_slab*slab size: "<<ncalo_slab*z2<<std::endl;
        //std::cout<<"\t\tprojected distance: "<<proj_distance<<" ke: "<<calo_ke<<std::endl;
        //std::cout<<"\t\tnext step: "<<next_step<<" next ke: "<<next_ke<<std::endl;
        //std::cout<<"\t\t\tcalo pt ("<<calo_x<<", "<<calo_y<<", "<<calo_z<<")\n";
        //std::cout<<"\t\t\tnext calo pt ("<<calo_next_x<<", "<<calo_next_y<<", "<<calo_next_z<<")\n";
        //std::cout<<"\n\tdistance to next slab: "<<r-proj_distance<<std::endl;
        //std::cout<<"\t\tslab pt: ("<<calo_slab_x<<", "<<calo_slab_y<<", "<<calo_slab_z<<")\n";
        //std::cout<<"\t\tslab ke: "<<calo_slab_ke<<std::endl;
        ncalo_slab++;
        calo_slab_xpos.push_back(calo_slab_x);
        calo_slab_ypos.push_back(calo_slab_y);
        calo_slab_zpos.push_back(calo_slab_z);
        calo_slab_KE.push_back(calo_slab_ke);
      }//<--End if this calo obj and the next step surround a slab
      calo_ke -= calo_de;
    }//<---End loop over reco calo objects to get slab information
    return 1;
}


//#############################################################################
//
// END 
//
//#############################################################################
