#define points_cxx
#include "points.h"
#include <TH2.h>
//#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// ### don't go more than 100 char on the width ---------------------------------------------------|
TH1D *hNumPts = new TH1D("hNumPts", "Number of spts inside the TPC", 50, 0, 2000);
TH1D *hDistanceBetweenPoints = new TH1D("hDistanceBetweenPoints", "dist between pts", 25, 0, .1);
TH1D *hfirstx = new TH1D("hfirstx", "dense first x", 100, 0, 50);
TH1D *hfirsty = new TH1D("hfirsty", "dense first x", 80,  0, 40);
TH1D *hfirstz = new TH1D("hfirstz", "dense first z", 50, 0, 10);
TH1D *hdedx = new TH1D("hdedx", "dedx", 500, 0, 50);
TH1D *hintke = new TH1D("hintke", "int ke", 20, 0, 1000);
TH1D *hincke = new TH1D("hincke", "inc ke", 20, 0, 1000);
TH1D *hxs    = new TH1D("hxs",    "xs",     20, 0, 1000);

TH1D *sDistanceBetweenSlabs = new TH1D("sDistanceBetweenSlabs", "distance between slabs", 100, 0,5);
TH1D *sdedx = new TH1D("sdedx", "slab dedx", 100, 0,   50);
TH1D *sxpos = new TH1D("sxpos", "slab xpos", 100, 0,   50);
TH1D *sypos = new TH1D("sypos", "slab ypos", 100, -20, 20);
TH1D *szpos = new TH1D("szpos", "slab zpos", 100, 0,   100);
TH1D *snslb = new TH1D("snslb", "slab numb", 100, 0,   100);
TH1D *sfirstx = new TH1D("sfirstx", "slab first x", 100, 0, 50);
TH1D *sfirsty = new TH1D("sfirsty", "slab first x", 80,  0, 40);
TH1D *sfirstz = new TH1D("sfirstz", "slab first z", 50, 0, 10);

TH1D *sincke = new TH1D("sincke", "slab inc ke", 20, 0, 1000);
TH1D *sintke = new TH1D("sintke", "slab int ke", 20, 0, 1000);
TH1D *sxs    = new TH1D("sxs",    "slab xs",     20, 0, 1000);

TH1D *ratio_int = new TH1D("ratio_int", "interaction ratio", 20, 0, 1000);
TH1D *ratio_inc = new TH1D("ratio_inc", "intcident ratio",   20, 0, 1000);
TH1D *ratio_xs  = new TH1D("ratio_xs",  "xs ratio",          20, 0, 1000);
TH1D *ratio_entries = new TH1D("ratio_entries", "ratio incident entries", 100, 0, 2);
TH2D *ahh = new TH2D("ahh", ":/", 100, 0, 2, 100, 0, 100);

TH1D *RDSptAngle = new TH1D("RDSptAngle", "Angle Between Spts", 1000, 0, 10);




double mass = 938.57;
double z = 0.03;
double z2 = .5;
double rho = 1.3954;                   // ## g/cm3
double molar_mass = 39.95;                    // ## g/mol
double N_A = 6.022 * pow(10, 23);      // ## num/mol

double recip_num_density = molar_mass / (rho * z * N_A); // ## cm2/num
double sparse_recip_num_density = molar_mass / (rho * z2 * N_A); // ## cm2/num
double barn = pow(10, -24);

// ## dumb debug counters ##
int TotSlabEntries = 0;
int TotDenseEntries = 0;
int nG4Interactions = 0;
int nRecoCandidates = 0;
int nRecoSignalEvts = 0;
int nRecoBckgrdEvts = 0;
int nTopology1 = 0;
int nTopology2 = 0;
int nTopology3 = 0;
int nTopology4 = 0;

bool print = true;
bool event_slection_print  = false;
bool geant4_print = false;
bool num_matching_print = false;

void points::Loop()
{
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   // ## event loop ##
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%1000 == 0){std::cout<<"nentries: "<<jentry<<std::endl;}
      if(print){std::cout<<"\n---------------------new event-------------------------"<<std::endl;}
      // ## printout vars ##
      bool dense_int = false;
      bool g4_interaction = false;
      double intx = 0;
      double inty = 0;
      double intz = 0;

      double initial_ke = 0;
      int num_pts_inTPC = 0;
      int first_pt = 0;
      int first_dense_slab = 0;
      double total_dense_dist = 0;
      int num_entries_dense = 0;

      int inslabs = 0;
      std::vector<double> true_slab_xpos;
      std::vector<double> true_slab_ypos;
      std::vector<double> true_slab_zpos;

      for(int g4part = 0; g4part < geant_list_size; g4part++){
        if((*process_primary)[g4part] != 1){continue;}
        if(geant4_print){
          std::cout << "dense spt loop " << (*NTrTrajPts)[g4part] << std::endl;
        }
        // ## dense ##
        for(int pt = 1; pt < (*NTrTrajPts)[g4part]; pt++){
          double xpos = (*MidPosX)[g4part][pt];
          double ypos = (*MidPosY)[g4part][pt];
          double zpos = (*MidPosZ)[g4part][pt];
          double prev_xpos = (*MidPosX)[g4part][pt-1];
          double prev_ypos = (*MidPosY)[g4part][pt-1];
          double prev_zpos = (*MidPosZ)[g4part][pt-1];
    
          double p = sqrt(pow(1000*(*MidPx)[g4part][pt], 2)
                        + pow(1000*(*MidPy)[g4part][pt], 2)
                        + pow(1000*(*MidPz)[g4part][pt], 2));
          double prev_p = sqrt(pow(1000*(*MidPx)[g4part][pt-1], 2)
                             + pow(1000*(*MidPy)[g4part][pt-1], 2)
                             + pow(1000*(*MidPz)[g4part][pt-1], 2));
          double ke = sqrt(pow(mass, 2) + pow(p, 2)) - mass; 
          double prev_ke = sqrt(pow(mass, 2) + pow(prev_p, 2)) - mass; 
          if( xpos > 0 && xpos < 47.5 && ypos > -20 && ypos < 20 && zpos > 0 && zpos < 90 ){
            if(first_pt == 0){initial_ke = ke; first_pt++;}
            num_pts_inTPC++;
            double dist_between_points = sqrt( pow(xpos - prev_xpos, 2) 
                                             + pow(ypos - prev_ypos, 2) 
                                             + pow(zpos - prev_zpos, 2));

            double energy_loss = prev_ke - ke;
            double de_dx = energy_loss / dist_between_points;
            total_dense_dist += dist_between_points;
            if(total_dense_dist < 1.){continue;}
            if(!first_dense_slab){
              first_dense_slab = 1;
              hfirstx->Fill(xpos);
              hfirsty->Fill(ypos);
              hfirstz->Fill(zpos);
            }
            // ### Fill Histos ###
            hDistanceBetweenPoints->Fill(dist_between_points);
            hdedx->Fill(de_dx);
            hincke->Fill(ke);
            num_entries_dense++; TotDenseEntries++;
            for(int nint = 0; nint < (*InteractionPoint).size(); nint++){
              if(pt == (*InteractionPoint)[nint]){
                //std::cout << "\t\tinteraction point" << std::endl;
                if((*InteractionPointType)[nint] == 13){
                  nG4Interactions++;
                  g4_interaction = true;
                  intx = xpos;inty = ypos;intz = zpos;
                  std::cout << "\t\t\tinelastic! " << pt << std::endl;
                  std::cout<<"\t\t\t\tx,y,z: "<<prev_xpos<<", "<<prev_ypos<<", "<<prev_zpos<<std::endl;
                  std::cout<<"\t\t\t\tx,y,z: "<<xpos<<", "<<ypos<<", "<<zpos<<std::endl;
                  dense_int = true;
                  hintke->Fill(prev_ke);
                }
              }
            }
          }//<---End if tpc
        }//<---End spt loop
        // ### using slabs ###
        //int inslabs = 0;
        double total_slab_distance = 0;
        //std::cout << "\tslab loop.." << std::endl;
        for(int slab = 0; slab < (*SlabN).size(); slab++){
          double slab_e = 1000*(*SlabE)[slab];
          double slab_x = (*SlabX)[slab];
          double slab_y = (*SlabY)[slab];
          double slab_z = (*SlabZ)[slab];
          double slab_p = sqrt( pow(1000*(*SlapX)[slab], 2) +
                                pow(1000*(*SlapY)[slab], 2) +
                                pow(1000*(*SlapZ)[slab], 2) );
          double slab_ke = sqrt( pow(slab_p, 2) + pow(mass, 2) ) - mass;
          if(slab_x > 0   && slab_x < 47.5 && 
             slab_y > -20 && slab_y < 20   && 
             slab_z > 0   && slab_z < 90){
            inslabs++;
            if(slab !=0){
              double slab_previous_p = sqrt( pow(1000*(*SlapX)[slab-1], 2) +
                                             pow(1000*(*SlapY)[slab-1], 2) +
                                             pow(1000*(*SlapZ)[slab-1], 2) ); 
              double slab_previous_ke = sqrt(pow(slab_previous_p,2)+pow(mass,2)) - mass;
              double slab_de = abs(slab_ke - slab_previous_ke);
              double slab_dx = sqrt( pow(slab_x - (*SlabX)[slab-1], 2) +
                                     pow(slab_y - (*SlabY)[slab-1], 2) +
                                     pow(slab_z - (*SlabZ)[slab-1], 2));
              double slab_dedx = slab_de / slab_dx;
              total_slab_distance += slab_dx;
              sDistanceBetweenSlabs->Fill(slab_dx);
              sdedx->Fill(slab_dedx);
            }
            else{sfirstx->Fill(slab_x);sfirsty->Fill(slab_y);sfirstz->Fill(slab_z);}
            true_slab_xpos.push_back(slab_x);
            true_slab_ypos.push_back(slab_y);
            true_slab_zpos.push_back(slab_z);
            sxpos->Fill(slab_x);
            sypos->Fill(slab_y);
            szpos->Fill(slab_z);
            // ## getting ke vars ##
            sincke->Fill(slab_ke);
            TotSlabEntries++;
          }//<---End if in tpc
        }//<---End slab loop
        snslb->Fill(inslabs);
        //if(fmod(total_dense_dist,(int)total_dense_dist) < .05){
        if(geant4_print){
          std::cout << "dense total distance: " << total_dense_dist << std::endl;
          std::cout << "slab total distance: " << total_slab_distance << std::endl;
          std::cout << "num slabs: " << inslabs << std::endl;
          std::cout << "\tnum entries dense: " << num_entries_dense << std::endl;
          std::cout << "\tratio: " << (.5*inslabs)/(z*num_entries_dense) << std::endl;
        }
        ratio_entries->Fill((1.*inslabs)/(z*num_entries_dense));
        ahh->Fill((1.*inslabs)/(z*num_entries_dense), total_dense_dist);
        
        // ## interaction debugging :/ ##
        if(dense_int){
          //std::cout<<"\nintdebug"<<std::endl;
          //std::cout << "number of slabs: " << inslabs << std::endl;
          //std::cout << "interaction point:  " << intx              << "\t\t"
          //                                    << inty              << "\t\t"
          //                                    << intz              << "\t\t" << std::endl;
          //std::cout << "last slab position: " << (*SlabX)[inslabs-1] << "\t\t" 
          //                                    << (*SlabY)[inslabs-1] << "\t\t"
          //                                    << (*SlabZ)[inslabs-1] << "\t\t" << std::endl;
          double dist_penult = sqrt( pow(intx - (*SlabX)[inslabs-1], 2) +
                                     pow(inty - (*SlabY)[inslabs-1], 2) +
                                     pow(intz - (*SlabZ)[inslabs-1], 2) );
          if(inslabs && dist_penult < 1){
            double slab_int_p = sqrt( pow(1000*(*SlapX)[inslabs-1], 2) +
                                      pow(1000*(*SlapY)[inslabs-1], 2) +
                                      pow(1000*(*SlapZ)[inslabs-1], 2) );
            double slab_int_ke = sqrt( pow(slab_int_p, 2) + pow(mass, 2)) - mass;
            sintke->Fill(slab_int_ke);   
          }//<---End if inelastic interaction 
        }
      }//<---End geant particle loop
      hNumPts->Fill(num_pts_inTPC);

      // ### Reconstructed Vars Below ###
      // ### Going to redo a lot of the work in the old macro ###
      // ### Hopefully it's clearer! ###
      if(print){
        std::cout<<"\n>>>>>>>>>>>>>>>>reco<<<<<<<<<<<<<<<<"<<std::endl;
      }

      // ### Sloppy version of primary proton selection ###
      // ### can change this when making the modular version ###
      // ### but I need to get a primary particle for now ###
      int reco_primary = -1;
      double first_reco_z = 99.;
      for(int rtrack = 0; rtrack < ntracks_reco; rtrack++){
        double z1 = (*track_zpos)[rtrack][0];
        if(print){
          std::cout << "first z point: " << z1 << std::endl;
        }
        if(z1 < first_reco_z && z1 < 2){
          first_reco_z = z1;
          reco_primary = rtrack;
        }
      }//<---End loop reco tracks
      if(reco_primary == -1){continue;}
      if(print){
        std::cout << "earliest primary: " << reco_primary << std::endl; 
      }



      // ### Inelastic Event Selection ###
      // ### In principle this is a cleaner version of this code block ###
      // ### Also technically more accurate ###
      // ### Changes I am making:
      //            - doing away with the bunch of if statements for number of tracks reco
      //            - catching earlier kinks and fork tracks 
      //            - trying to do all of this and keep the topologies as is! (are?)
      std::vector<std::vector<double>> branches;
      std::vector<std::vector<double>> kinks;
      int missing_bragg = 0;
      double prim_endx = (*track_end_x)[reco_primary];
      double prim_endy = (*track_end_y)[reco_primary];
      double prim_endz = (*track_end_z)[reco_primary];
      // # track loop #
      for(int rtrack = 0; rtrack < ntracks_reco; rtrack++){
        if(rtrack == reco_primary){
          std::cout<<"t1 debug\n";
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
            RDSptAngle->Fill(theta);
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
          std::cout<<"calo obj loop"<<std::endl;
          for(int calo_pt = (*col_track_hits)[rtrack]; calo_pt > 0; calo_pt--){
            if(end_track_dist > 2.5){continue;}
            end_track_dist += (*col_track_pitch_hit)[rtrack][calo_pt];
            end_track_dedx_sum += (*col_track_dedx)[rtrack][calo_pt];
            end_track_counter++;
          }
          double end_track_dedx_mean = end_track_dedx_sum / end_track_counter;
          if(end_track_dedx_mean < 13){ // need to make this number a variable at some point
            missing_bragg++;
          }


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
        }
      }
      // # topology 2
      if(!(branches.size() || kinks.size())){
        if(missing_bragg && prim_endz < 88){
          nTopology2++;
          interacting_candidate = true;
          candidate_spt = (*ntrack_hits)[reco_primary] - 1;
        }
        else{
          interacting_candidate = false;
        }
      }

      // ### MC ONLY ###
      // ### Comparing interaction candidate position to truth info ###
      // ### if the int candidate is within a certain distance of a true inelastic interaction ###
      // ### then, pass a label to the histo filling (singal or background) ###
      if(interacting_candidate){ 
        nRecoCandidates++;
        bool signal;
        double reco_x = (*track_xpos)[reco_primary][candidate_spt];
        double reco_y = (*track_ypos)[reco_primary][candidate_spt];
        double reco_z = (*track_zpos)[reco_primary][candidate_spt];
        std::cout<<"I found a potential interaction here\n";
        std::cout<<"\t\t("<<reco_x<<", "<<reco_y<<", "<<reco_z<<")\n";
        // mark as signal or background ...
        if(g4_interaction){
          std::cout<<"\nthere was an inelastic yeah?"<<std::endl;
          std::cout<<"\t\t("<<intx<<", "<<inty<<", "<<intz<<")\n";
          double reco_x = (*track_xpos)[reco_primary][candidate_spt];
          double reco_y = (*track_ypos)[reco_primary][candidate_spt];
          double reco_z = (*track_zpos)[reco_primary][candidate_spt];
          double dist_reco_g4 = sqrt( pow(reco_x-intx,2)
                                    + pow(reco_y-inty,2)
                                    + pow(reco_z-intz,2) );
          if(dist_reco_g4 < 2){
            std::cout<<"signal event!\n";
            signal = true;
            nRecoSignalEvts++;
            // pass a label to histo filling
          }
          else{
            std::cout<<"background event!\n";
            signal = false;
            nRecoBckgrdEvts++;
            // pass a label to histo filling
          }
        }//<--End if flag an interaction
        else{
          std::cout<<"background event!\n";
          signal = false;
          nRecoBckgrdEvts++;
          // pass a label to histo filling
        }
      }//<--End if flagged an interaction 


      // ### DENOMINATOR ###
      // ### also a cleaner version of this code block ###
      // ### writing it to look more like how it would end up being used in a modular version ###
      double calo_ke = initial_ke;
      std::vector<double> calo_slab_xpos;
      std::vector<double> calo_slab_ypos;
      std::vector<double> calo_slab_zpos;
      std::vector<double> calo_slab_KE;
      double proj_distance = 0;
      double next_step = 0;
      int ncalo_slab = 1;
      for(int calo_pt = 0; calo_pt < (*col_track_hits)[reco_primary] - 1; calo_pt++){
        proj_distance += (*col_track_pitch_hit)[reco_primary][calo_pt];
        next_step = proj_distance + (*col_track_pitch_hit)[reco_primary][calo_pt+1];
        double calo_x = (*col_track_x)[reco_primary][calo_pt]; 
        double calo_y = (*col_track_y)[reco_primary][calo_pt]; 
        double calo_z = (*col_track_z)[reco_primary][calo_pt]; 
        double calo_de = (*col_track_dedx)[reco_primary][calo_pt]*
                         (*col_track_pitch_hit)[reco_primary][calo_pt];
        double next_ke = calo_ke - calo_de; 
        double calo_next_x = (*col_track_x)[reco_primary][calo_pt+1]; 
        double calo_next_y = (*col_track_y)[reco_primary][calo_pt+1]; 
        double calo_next_z = (*col_track_z)[reco_primary][calo_pt+1]; 
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
        }
        calo_ke -= calo_de;
      }
      
      // ### MC ONLY ###
      // ### taking the slabs that will build incident distribution and comparing to mc ###
      // ### then pass labels for signal/background distributions ###
      std::cout<<"\n``````````````````\n";
      //std::cout<<"number of true slabs: "<<inslabs<<std::endl;
      //std::cout<<"number of calo slabs: "<<calo_slab_KE.size()<<std::endl;
      std::vector<int> calo_slab_signal(calo_slab_KE.size(), 0);
      int calo_slab_counter = 0;
      for(int calo_slab = 0; calo_slab < calo_slab_KE.size(); calo_slab++){
        double calo_slab_x = calo_slab_xpos[calo_slab];
        double calo_slab_y = calo_slab_ypos[calo_slab];
        double calo_slab_z = calo_slab_zpos[calo_slab];
        //std::cout<<"\tthis slab pos (x,y,z): ("<<calo_slab_xpos[calo_slab]
        //                                 <<", "<<calo_slab_ypos[calo_slab]
        //                                 <<", "<<calo_slab_zpos[calo_slab]<<")\n";
        double min_dist_slab = 99;
        double min_downstream_slab = 99;
        double min_upstream_slab = 99;
        int closest_downstream_slab = -1;
        int closest_upstream_slab = -1;
        for(int true_slab = 0; true_slab < inslabs; true_slab++){
          double true_slab_x = true_slab_xpos[true_slab];
          double true_slab_y = true_slab_ypos[true_slab];
          double true_slab_z = true_slab_zpos[true_slab];
          double dist_calo_true = sqrt( pow(calo_slab_x - true_slab_x, 2) +
                                        pow(calo_slab_y - true_slab_y, 2) +
                                        pow(calo_slab_z - true_slab_z, 2));
          if(dist_calo_true < min_dist_slab){
            min_dist_slab = dist_calo_true;
          }
          // # true slabs downstream calo #
          if(true_slab_z > calo_slab_z){
            if(dist_calo_true < min_downstream_slab){
              min_downstream_slab = dist_calo_true;
              closest_downstream_slab = true_slab;
            }
          }
          // # true slabs upstream calo #
          if(true_slab_z < calo_slab_z){
            if(dist_calo_true < min_upstream_slab){
              min_upstream_slab = dist_calo_true;
              closest_upstream_slab = true_slab;
            }
          }
        }//<--End loop of true slabs
        //std::cout<<"\t\tdistance to closest true slab: "<<min_dist_slab<<std::endl;
        //std::cout<<"\t\tdistance to closest ds true slab: "<<min_downstream_slab<<std::endl;
        //std::cout<<"\t\tdistance to closest us true slab: "<<min_upstream_slab<<std::endl;
        // ## signal background! ##
        if(closest_downstream_slab != -1 && closest_upstream_slab != -1){
          //std::cout<<"\t\t\tthere's a slab on both side!"<<std::endl;
          //std::cout<<"\t\tdistance to closest true slab: "<<min_dist_slab<<std::endl;
          //std::cout<<"\t\tdistance to closest ds true slab: "<<min_downstream_slab<<std::endl;
          //std::cout<<"\t\tdistance to closest us true slab: "<<min_upstream_slab<<std::endl;
          calo_slab_signal[calo_slab] = 1;
          calo_slab_counter++;
        }//<--End if surrounded!
        else{
          // # check if it's the first slab
            // # check distance to closest downstream true slab
          if(calo_slab == 0){
            //std::cout<<"\t\t\tthis is the first slab :0"<<std::endl;
            calo_slab_signal[calo_slab] = 1;
            calo_slab_counter++;
          }
          if(calo_slab_signal[calo_slab-1] == 1){
            //std::cout<<"\t\t\tthis slab is not surrounded but the previous one was\n";
            if(min_upstream_slab < 1){calo_slab_signal[calo_slab] = 1; calo_slab_counter++;}
          }
        }//<--End if not surrounded
      }
      std::cout<<"number of true slabs: "<<inslabs<<std::endl;
      std::cout<<"number of calo slabs: "<<calo_slab_KE.size()<<std::endl;
      std::cout<<"number of signal calo slabs: "<<calo_slab_counter<<std::endl;

      // ### MORE HISTOGRAM FILLING ###
      // ### Take calo slabs and take interaction point candidates ###
      // ### use to fill histograms appropriately ###


   }//<---End tree event loop


   // ### printouts, diagnostics, counters ###
    std::cout<<"=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\n";
    std::cout<<"Number of G4 Inelastic Events: "<<nG4Interactions<<std::endl;
    std::cout<<"Number of Reco Candidates: "<<nRecoCandidates<<std::endl;
    std::cout<<"\tSignal: "<<nRecoSignalEvts<<std::endl;
    std::cout<<"\tBackground: "<<nRecoBckgrdEvts<<std::endl;
    std::cout<<"\t\tTopology1: "<<nTopology1<<std::endl;
    std::cout<<"\t\tTopology2: "<<nTopology2<<std::endl;
    std::cout<<"\t\tTopology3: "<<nTopology3<<std::endl;
    std::cout<<"\t\tTopology4: "<<nTopology4<<std::endl;

   // ## xs calculation ##
   for(int iBin = 0; iBin < hintke->GetNbinsX(); iBin++){
    if(hincke->GetBinContent(iBin) == 0){continue;}
    double ratio   = 1.*hintke->GetBinContent(iBin) / hincke->GetBinContent(iBin);  
    double temp_xs = ratio * recip_num_density; 
    double xs      = temp_xs / barn;

    double num_err = pow(hintke->GetBinContent(iBin), .5);
    double num = hintke->GetBinContent(iBin);
    if(num == 0){num = 1;}
    double term1 = num_err/num;
    double dem_err = pow(hincke->GetBinContent(iBin), .5);
    double dem = hincke->GetBinContent(iBin);
    if(dem == 0){dem =1;}
    double term2 = dem_err/dem;
    double totalError = temp_xs*pow(pow(term1,2) + pow(term2,2), 0.5)/barn;//*recip_num_density*barn;

    //std::cout <<"xs: " << xs << " +- " << totalError << std::endl;
    hxs->SetBinContent(iBin, xs); 
    hxs->SetBinError(iBin,totalError);
   }
   for(int iBin = 0; iBin < sintke->GetNbinsX(); iBin++){
    if(sincke->GetBinContent(iBin) == 0){continue;}
    double ratio   = 1.*sintke->GetBinContent(iBin) / sincke->GetBinContent(iBin);  
    double temp_xs = ratio * sparse_recip_num_density; 
    double xs      = temp_xs / barn;

    double num_err = pow(sintke->GetBinContent(iBin), .5);
    double num = sintke->GetBinContent(iBin);
    if(num == 0){num = 1;}
    double term1 = num_err/num;
    double dem_err = pow(sincke->GetBinContent(iBin), .5);
    double dem = sincke->GetBinContent(iBin);
    if(dem == 0){dem =1;}
    double term2 = dem_err/dem;
    double totalError = temp_xs*pow(pow(term1,2) + pow(term2,2), 0.5)/barn;//*recip_num_density*barn;

    //std::cout <<"xs: " << xs << " +- " << totalError << std::endl;
    sxs->SetBinContent(iBin, xs); 
    sxs->SetBinError(iBin,totalError);
   }
   // ## ratio plots ##
   for(int iBin = 0; iBin < sintke->GetNbinsX(); iBin++){
    if(hintke->GetBinContent(iBin)){
      double int_ratio = sintke->GetBinContent(iBin) / hintke->GetBinContent(iBin);
      ratio_int->SetBinContent(iBin, int_ratio);
    }
    if(hincke->GetBinContent(iBin)){
      double inc_ratio = (z2*sincke->GetBinContent(iBin)) / (z*hincke->GetBinContent(iBin));
      ratio_inc->SetBinContent(iBin, inc_ratio);
    }
   }
   // ## normalization debuggin ##
   double normalization_ratio = (z2*TotSlabEntries)/(z*TotDenseEntries);
   //std::cout << "++++++++++++++++++++++++++++" << std::endl;
   //std::cout << "Total Dense Entries: " << TotDenseEntries << std::endl;
   //std::cout << "Total Slabs Entries: " << TotSlabEntries  << std::endl;
   //std::cout << "ratio: " << normalization_ratio << std::endl;

   // ## writing histograms ##
   TFile myfile("spacepoints.root", "RECREATE");
   hNumPts->Write();
   hDistanceBetweenPoints->Write();
   hdedx->Write();
   hfirstx->Write();
   hfirsty->Write();
   hfirstz->Write();
   hintke->Write();
   hincke->Write();
   hxs->Write();
   sDistanceBetweenSlabs->Write();
   sdedx->Write();
   sxpos->Write();
   sypos->Write();
   szpos->Write();
   sfirstx->Write();
   sfirsty->Write();
   sfirstz->Write();
   snslb->Write();
   sincke->Write();
   sintke->Write();
   sxs->Write();
   ratio_int->Write();
   ratio_inc->Write();
   ratio_entries->Write();
   ahh->Write();
   RDSptAngle->Write();
}//<---End macro
