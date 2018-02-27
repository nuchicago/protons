#define points_cxx
#include "points.h"
#include <TH2.h>
//#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// ### don't go more than 110 char on the width -------------------------------------------------------------|
// ### ---- this was more just a note for me to try and make things more readable 
// ### ---- didn't really succeed everywhere...


// === Diagnositcs ===
TH1D *hNumPts = new TH1D("hNumPts", "Number of spts inside the TPC", 50, 0, 2000);
TH1D *hDistanceBetweenPoints = new TH1D("hDistanceBetweenPoints", "dist between pts", 25, 0, .1);
TH2D *hdistvske = new TH2D("hdistvske", "dist btwn pts vs ke", 1000, 0, .1, 20, 0, 1000);
TH1D *hfirstx = new TH1D("hfirstx", "dense first x", 100, 0, 50);
TH1D *hfirsty = new TH1D("hfirsty", "dense first x", 80,  0, 40);
TH1D *hfirstz = new TH1D("hfirstz", "dense first z", 50, 0, 10);
TH1D *hdedx = new TH1D("hdedx", "dedx", 500, 0, 50);


// === Deliverables ===
TH1D *hintke = new TH1D("hintke", "int ke", 20, 0, 1000);
TH1D *hincke = new TH1D("hincke", "inc ke", 20, 0, 1000);
TH1D *h2incke = new TH1D("h2incke", "inc ke", 20, 0, 1000);
TH1D *hxs    = new TH1D("hxs",    "xs",     20, 0, 1000);

// === Diagnostics ===
TH1D *mintke = new TH1D("mintke", "int ke", 20, 0, 1000);
TH1D *mincke = new TH1D("mincke", "inc ke", 20, 0, 1000);
TH1D *mxs = new TH1D("mxs", "xs", 20, 0, 1000);


// === Diagnostics ===
TH1D *sDistanceBetweenSlabs = new TH1D("sDistanceBetweenSlabs", "distance between slabs", 100, 0,5);
TH1D *sdedx = new TH1D("sdedx", "slab dedx", 100, 0,   50);
TH1D *sxpos = new TH1D("sxpos", "slab xpos", 100, 0,   50);
TH1D *sypos = new TH1D("sypos", "slab ypos", 100, -20, 20);
TH1D *szpos = new TH1D("szpos", "slab zpos", 100, 0,   100);
TH1D *snslb = new TH1D("snslb", "slab numb", 100, 0,   100);
TH1D *sfirstx = new TH1D("sfirstx", "slab first x", 100, 0, 50);
TH1D *sfirsty = new TH1D("sfirsty", "slab first x", 80,  0, 40);
TH1D *sfirstz = new TH1D("sfirstz", "slab first z", 50, 0, 10);

// === Deliverables ===
TH1D *sincke = new TH1D("sincke", "slab inc ke", 20, 0, 1000);
TH1D *sintke = new TH1D("sintke", "slab int ke", 20, 0, 1000);
TH1D *sxs    = new TH1D("sxs",    "slab xs",     20, 0, 1000);

// === Diagnostics ===
TH2D *slab_vs_dense_intke = new TH2D("slab_vs_dense_intke", "slab vs dense int ke", 20, 0, 1000, 20, 0, 1000);
TH1D *ratio_int = new TH1D("ratio_int", "interaction ratio", 20, 0, 1000);
TH1D *ratio_inc = new TH1D("ratio_inc", "intcident ratio",   20, 0, 1000);
TH1D *ratio_xs  = new TH1D("ratio_xs",  "xs ratio",          20, 0, 1000);
TH1D *ratio_entries = new TH1D("ratio_entries", "ratio incident entries", 100, 0, 2);
TH2D *ahh = new TH2D("ahh", ":/", 100, 0, 2, 100, 0, 100);
TH1D *RDSptAngle = new TH1D("RDSptAngle", "Angle Between Spts", 1000, 0, 10);


// === Deliverables (xs ingredients!) ===
TH1D *hreco_incke = new TH1D("hreco_incke", "incident ke", 20, 0, 1000);
TH1D *hreco_incke_signal = new TH1D("hreco_incke_signal", "incident ke (signal)", 20, 0, 1000);
TH1D *hreco_folded_incke_signal = new TH1D("hreco_folded_incke_signal", "incident ke (signal)", 20, 0, 1000);
TH1D *hreco_unfolded_incke_signal = new TH1D("hreco_unfolded_incke_signal", "incident ke (signal)", 20, 0, 1000);
TH1D *hreco_incke_background = new TH1D("hreco_incke_background", "incident ke (background)", 20, 0, 1000);
TH1D *hreco_intke = new TH1D("hreco_intke", "interacting ke", 20, 0, 1000);
TH1D *hreco_intke_signal = new TH1D("hreco_intke_signal", "interacting ke (signal)", 20, 0, 1000);
TH1D *hreco_folded_intke_signal = new TH1D("hreco_folded_intke_signal", "interacting ke (signal)", 20, 0, 1000);
TH1D *hreco_unfolded_intke_signal = new TH1D("hreco_unfolded_intke_signal", "interacting ke (signal)", 20, 0, 1000);
TH1D *hreco_intke_background = new TH1D("hreco_intke_background", "interacting ke (background)", 20, 0, 1000);
TH1D *hreco_intke_eff = new TH1D("hreco_int_eff", "interacting selection efficiency", 20, 0, 1000);
TH1D *hreco_incke_eff = new TH1D("hreco_inc_eff", "incident selection efficiency", 20, 0, 1000);
TH2D *hreco_unfolding_matrix=new TH2D("hreco_unfolding_matrix","energy unfolding matrix",20,0,1000,20,0,1000);
TH2D *hreco_unfolding_matrix_normalized=new TH2D("hreco_unfolding_matrix_normalized","energy unfolding matrix",20,0,1000,20,0,1000);
TH1D *hreco_xs = new TH1D("hreco_xs", "p-ar inelastic xs", 20, 0, 1000);

// ## constants for xs measurement ##
double mass = 938.57;
// z: dense g4 spts. z2: sparse g4 spts and reco slab size. z3: trying to prove a point about slab size.
double z = 0.03;
double z2 = .5;
double z3 = .1;

double rho = 1.3954;                   // ## g/cm3
double molar_mass = 39.95;                    // ## g/mol
double N_A = 6.022 * pow(10, 23);      // ## num/mol

double recip_num_density = molar_mass / (rho * z * N_A); // ## cm2/num
double dense_recip_num_density = molar_mass / (rho * z3 * N_A); // ## cm2/num
double sparse_recip_num_density = molar_mass / (rho * z2 * N_A); // ## cm2/num
double barn = pow(10, -24);

// ## dumb debug counters ##
// ## these help make sure things are working properly ##
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

// ## flags for whether or not we want to print things
bool print = false;
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
      // ## geant 4 interaction variables ##
      bool dense_int = false;
      bool g4_interaction = false;
      double intx = 0;
      double inty = 0;
      double intz = 0;
      double int_ke = -1;

      bool g4_int = false;
      double g4_intx = -99;
      double g4_inty = -99;
      double g4_intz = -99;

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
      std::vector<double> true_slab_ke;


	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~ Start Code Block A ~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      // ### MC ONLY ###
	  // ### big geant 4 particle loop ###
      for(int g4part = 0; g4part < geant_list_size; g4part++){
        if((*process_primary)[g4part] != 1){continue;}
        if(geant4_print){
          std::cout << "dense spt loop " << (*NTrTrajPts)[g4part] << std::endl;
        }

        for(int nint = 0; nint < (*InteractionPoint).size(); nint++){
          if((*InteractionPointType)[nint] == 13){
            g4_intx = (*MidPosX)[g4part][(*InteractionPoint)[nint]];
            g4_inty = (*MidPosY)[g4part][(*InteractionPoint)[nint]];
            g4_intz = (*MidPosZ)[g4part][(*InteractionPoint)[nint]];
            g4_int = true;
            //std::cout<<"found inelastic outside loop\n";
            //std::cout<<"\t("<<g4_intx<<", "<<g4_inty<<", "<<g4_intz<<")\n";
          }//<--End if the interaction is inelastic
        }//<--End loop over all interactions (could be 0!)



        // ## dense ##
        double g4_proj_dist = 0;
        double prev_g4_proj_dist = 0;
        int ndense_slab = 1;
        double last_dense_slabx = -99;
        double last_dense_slaby = -99;
        double last_dense_slabz = -99;
        double last_dense_slabke = -99;

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
            h2incke->Fill(ke, dist_between_points/z);
            hdistvske->Fill(dist_between_points, ke);
            num_entries_dense++; TotDenseEntries++;
            for(int nint = 0; nint < (*InteractionPoint).size(); nint++){
              if(pt == (*InteractionPoint)[nint]){
                //std::cout << "\t\tinteraction point" << std::endl;
                if((*InteractionPointType)[nint] == 13){
                  nG4Interactions++;
                  g4_interaction = true;
                  intx = xpos;inty = ypos;intz = zpos;
                  int_ke = prev_ke;
                  //std::cout << "\t\t\tinelastic! " << pt << std::endl;
                  //std::cout<<"\t\t\t\tx,y,z: "<<prev_xpos<<", "<<prev_ypos<<", "<<prev_zpos<<std::endl;
                  //std::cout<<"\t\t\t\tx,y,z: "<<xpos<<", "<<ypos<<", "<<zpos<<std::endl;
                  dense_int = true;
                  hintke->Fill(prev_ke);
                }//<--End if the interaction is inelastic
              }//<--End if the point we're looking at is the interaction point
            }//<--End loop over all interactions (could be 0!)

            // ## using .1cm slabs ##
            prev_g4_proj_dist = g4_proj_dist;
            g4_proj_dist += dist_between_points; 
            //std::cout<<"prev spt total distance: "<<prev_g4_proj_dist<<std::endl;
            //std::cout<<"dense spt total distance: "<<g4_proj_dist<<std::endl; 
            //std::cout<<"next slab: "<<ndense_slab*z3<<std::endl;
            if(prev_g4_proj_dist < ndense_slab*z3 && g4_proj_dist > ndense_slab*z3){
              double g4_vx = xpos - prev_xpos;
              double g4_vy = ypos - prev_ypos;
              double g4_vz = zpos - prev_zpos;
              double g4_r = ndense_slab*z3;
              double g4_A = pow(g4_vx,2) + pow(g4_vy,2) + pow(g4_vz,2);  
              double g4_B = 2*(prev_xpos*g4_vx + prev_ypos*g4_vy + prev_zpos*g4_vz);
              double g4_C = pow(prev_g4_proj_dist,2) - pow(g4_r,2);
              double g4_t = (-1*g4_B + pow( pow(g4_B,2) - 4*g4_A*g4_C, .5))/(2*g4_A);
              double dense_slab_x = prev_xpos + g4_t*g4_vx;
              double dense_slab_y = prev_ypos + g4_t*g4_vy;
              double dense_slab_z = prev_zpos + g4_t*g4_vz;
              double dense_ke_slope = (ke - prev_ke) / (g4_proj_dist - prev_g4_proj_dist);
              double dense_slab_ke = dense_ke_slope*(g4_r - prev_g4_proj_dist) + prev_ke;
              last_dense_slabx = dense_slab_x;
              last_dense_slaby = dense_slab_y;
              last_dense_slabz = dense_slab_z;
              last_dense_slabke = dense_slab_ke;
              mincke->Fill(dense_slab_ke);
              //if(g4_int){
              //  double dist_to_int = sqrt(pow(dense_slab_x - g4_intx, 2) + 
              //                            pow(dense_slab_y - g4_inty, 2) + 
              //                            pow(dense_slab_z - g4_intz, 2)); 
              //  //std::cout<<"distance to interaction: "<<dist_to_int<<std::endl;
              //  //std::cout<<"\tp0: ("<<prev_xpos<<", "<<prev_ypos<<", "<<prev_zpos<<")\n";
              //  //std::cout<<"\tke: "<<prev_ke<<std::endl;
              //  //std::cout<<"\tp1: ("<<xpos<<", "<<ypos<<", "<<zpos<<")\n";
              //  //std::cout<<"\tke: "<<ke<<std::endl;
              //  //std::cout<<"\tslab pt: ("<<dense_slab_x<<", "<<dense_slab_y<<", "<<dense_slab_z<<")\n";
              //  //std::cout<<"\tslab ke: "<<dense_slab_ke<<std::endl;
              //  //std::cout<<std::endl;
              //  if(dist_to_int < .115){
              //    std::cout<<"\t\t\t\tfilling interacting histo\n";
              //    mintke->Fill(dense_slab_ke);
              //  }
              //}//<--End if there was an inelastic interaction in this event
              ndense_slab++;
            }//<--End if these two g4 spts surround a .1cm slab
          }//<---End if tpc
        }//<---End spt loop
        if(g4_int){
          double dist_to_int = sqrt(pow(last_dense_slabx - g4_intx, 2) + 
                                    pow(last_dense_slaby - g4_inty, 2) + 
                                    pow(last_dense_slabz - g4_intz, 2)); 
          if(g4_intz > 1 && g4_intz < 90){
            mintke->Fill(last_dense_slabke);
          }
        }//<--End if there was an inelastic interaction in this event
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
          double zero_protection = 0;
          if(slab_x > 0   && slab_x < 47.5 && 
             slab_y > -20 && slab_y < 20   && 
             slab_z > 0   && slab_z < 90){
            inslabs++;
            if(slab !=0){
              double slab_previous_p = sqrt( pow(1000*(*SlapX)[slab-1], 2) +
                                             pow(1000*(*SlapY)[slab-1], 2) +
                                             pow(1000*(*SlapZ)[slab-1], 2) ); 
              double slab_previous_ke = sqrt(pow(slab_previous_p,2)+pow(mass,2)) - mass;
              zero_protection = slab_previous_ke;
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
            if(slab_ke){
              sincke->Fill(slab_ke);
              true_slab_ke.push_back(slab_ke);
            }
            if(!slab_ke){
              sincke->Fill(zero_protection);
              true_slab_ke.push_back(zero_protection);
            }
            TotSlabEntries++;
          }//<---End if in tpc
        }//<---End slab loop
        snslb->Fill(inslabs);
        //if(fmod(total_dense_dist,(int)total_dense_dist) < .05){}
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

          //double dist_penult = sqrt( pow(intx - (*SlabX)[inslabs-1], 2) +
          //                           pow(inty - (*SlabY)[inslabs-1], 2) +
          //                           pow(intz - (*SlabZ)[inslabs-1], 2) );
          //if(inslabs && dist_penult < 1){
          //  double slab_int_p = sqrt( pow(1000*(*SlapX)[inslabs-1], 2) +
          //                            pow(1000*(*SlapY)[inslabs-1], 2) +
          //                            pow(1000*(*SlapZ)[inslabs-1], 2) );
          //  double slab_int_ke = sqrt( pow(slab_int_p, 2) + pow(mass, 2)) - mass;
          //  sintke->Fill(slab_int_ke);   
          //}//<---End if inelastic interaction 
          if(true_slab_ke[inslabs-1]){
            sintke->Fill(true_slab_ke[inslabs-1]);
          }
          else{
            sintke->Fill(true_slab_ke[inslabs-2]);
          }
          slab_vs_dense_intke->Fill(true_slab_ke[inslabs-1], int_ke);
          //if(!true_slab_ke[inslabs-1]){
          //  std::cout<<"dense ke: "<<int_ke<<std::endl;
          //  std::cout<<"slab ke: "<<true_slab_ke[inslabs-1]<<std::endl;
          //  std::cout<<"prev slab ke: "<<true_slab_ke[inslabs-2]<<std::endl;
          //  std::cout<<std::endl;
          //}

        }//<--End if this event had an inelastic interaction


      }//<---End geant particle loop
      hNumPts->Fill(num_pts_inTPC);
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~~ End Code Block A ~~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~ Start Code Block B ~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~~ End Code Block B ~~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~ Start Code Block C ~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
          nTopology2++;
          interacting_candidate = true;
          candidate_spt = (*ntrack_hits)[reco_primary] - 1;
          candidate_xpos = (*col_track_x)[reco_primary][(*col_track_hits)[reco_primary]-1];
          candidate_ypos = (*col_track_y)[reco_primary][(*col_track_hits)[reco_primary]-1];
          candidate_zpos = (*col_track_z)[reco_primary][(*col_track_hits)[reco_primary]-1];
        }
        else{
          interacting_candidate = false;
        }
      }

	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~~ End Code Block C ~~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~ Start Code Block D ~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      // ### MC ONLY ###
      // ### Comparing interaction candidate position to truth info ###
      // ### if the int candidate is within a certain distance of a true inelastic interaction ###
      // ### then, pass a label to the histo filling (singal or background) ###
      bool signal = false;
      if(interacting_candidate){ 
        nRecoCandidates++;
        //bool signal = false;
        //double reco_x = (*track_xpos)[reco_primary][candidate_spt];
        //double reco_y = (*track_ypos)[reco_primary][candidate_spt];
        //double reco_z = (*track_zpos)[reco_primary][candidate_spt];
        double reco_x = candidate_xpos;
        double reco_y = candidate_ypos;
        double reco_z = candidate_zpos;
        //std::cout<<"I found a potential interaction here\n";
        //std::cout<<"\t\t("<<reco_x<<", "<<reco_y<<", "<<reco_z<<")\n";
        // mark as signal or background ...
        if(g4_interaction){
          //std::cout<<"\nthere was an inelastic yeah?"<<std::endl;
          //std::cout<<"\t\t("<<intx<<", "<<inty<<", "<<intz<<")\n";
          double reco_x = (*track_xpos)[reco_primary][candidate_spt];
          double reco_y = (*track_ypos)[reco_primary][candidate_spt];
          double reco_z = (*track_zpos)[reco_primary][candidate_spt];
          double dist_reco_g4 = sqrt( pow(reco_x-intx,2)
                                    + pow(reco_y-inty,2)
                                    + pow(reco_z-intz,2) );
          if(dist_reco_g4 < 2){
            //std::cout<<"signal event!\n";
            signal = true;
            nRecoSignalEvts++;
            // pass a label to histo filling
          }
          else{
            //std::cout<<"background event!\n";
            signal = false;
            nRecoBckgrdEvts++;
            // pass a label to histo filling
          }
        }//<--End if flag an interaction
        else{
          //std::cout<<"background event!\n";
          signal = false;
          nRecoBckgrdEvts++;
          // pass a label to histo filling
        }
      }//<--End if flagged an interaction 

	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~~ End Code Block D ~~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~ Start Code Block E ~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
        }//<--End if this calo obj and the next step surround a slab
        calo_ke -= calo_de;
      }//<---End loop over reco calo objects to get slab information

	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~~ End Code Block E ~~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~ Start Code Block F ~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      
      // ### MC ONLY ###
      // ### taking the slabs that will build incident distribution and comparing to mc ###
      // ### then pass labels for signal/background distributions ###
      //std::cout<<"\n``````````````````\n";
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
            if(min_upstream_slab < .5){calo_slab_signal[calo_slab] = 1; calo_slab_counter++;}
          }
        }//<--End if not surrounded
      }//<---End loop on calo level slabs

	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~~ End Code Block F ~~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~ Start Code Block G ~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      // ### MORE HISTOGRAM FILLING ###
      // ### Take calo slabs and take interaction point candidates ###
      // ### use to fill histograms appropriately ###

      // ## getting which slab will be used for interactions ##
      int calo_int_slab = 999;
      if(interacting_candidate){
        //double int_candidate_x = (*track_xpos)[reco_primary][candidate_spt];
        //double int_candidate_y = (*track_ypos)[reco_primary][candidate_spt];
        //double int_candidate_z = (*track_zpos)[reco_primary][candidate_spt];
        double int_candidate_x = candidate_xpos;
        double int_candidate_y = candidate_ypos;
        double int_candidate_z = candidate_zpos;
        // loop over slabs to find slab closest to interaction candidate
        double min_dist_int = 99;
        for(int calo_slab = 0; calo_slab < calo_slab_KE.size(); calo_slab++){
          double calo_slab_x = calo_slab_xpos[calo_slab];  
          double calo_slab_y = calo_slab_ypos[calo_slab];  
          double calo_slab_z = calo_slab_zpos[calo_slab];  
          double dist_int_slab = sqrt( pow(int_candidate_x - calo_slab_x, 2) 
                                     + pow(int_candidate_y - calo_slab_y, 2) 
                                     + pow(int_candidate_z - calo_slab_z, 2) );
          if(calo_slab_z < int_candidate_z){
            if(dist_int_slab < min_dist_int){
              min_dist_int = dist_int_slab;
              calo_int_slab = calo_slab;
            }
          }//<--End if this slab is upstream of int
        }//<--End calo slab loop
        //std::cout<<"\n]]]]]]]]]]]]]]]\n";
        //std::cout<<"found possible interaction\n";
        //std::cout<<"\tint slab: "<<calo_int_slab<<std::endl;
        //std::cout<<"\tdist btwn int candidate and closest slab: "<<min_dist_int<<std::endl;
      }//<---End if interaction candidate

      // ## incident slabs ## 
      int ninc_entries = 0;
      for(int calo_slab = 1; calo_slab < calo_slab_KE.size(); calo_slab++){
        if(calo_slab > calo_int_slab){continue;}//<--stop after interaction slab 
        //std::cout<<"\tincident entry: "<<std::endl;
        //std::cout<<"\t\tke: "<<calo_slab_KE[calo_slab]<<std::endl;
        //std::cout<<"\t\tsignal? "<<calo_slab_signal[calo_slab]<<std::endl;
        hreco_incke->Fill(calo_slab_KE[calo_slab]);
        if(calo_slab_signal[calo_slab]){
          hreco_incke_signal->Fill(calo_slab_KE[calo_slab]); 
          hreco_unfolding_matrix->Fill(calo_slab_KE[calo_slab], true_slab_ke[calo_slab-1]);
          ninc_entries++;
        }
        if(!calo_slab_signal[calo_slab]){
          hreco_incke_background->Fill(calo_slab_KE[calo_slab]);
        }
        if(calo_slab == calo_int_slab){
          //std::cout<<"\tinteraction energy: "<<calo_slab_KE[calo_slab]<<std::endl;
          // ###### need to decide what to do with interaction pts far away from slab
          // ###### likely a feature of differences between reco spts and calo objs :/
          // ## also need to check on non-physical entries (negative ?)
          // ### should probably also just grab any non terminating protons as interactions...
          hreco_intke->Fill(calo_slab_KE[calo_slab]);
          if(signal){
            hreco_intke_signal->Fill(calo_slab_KE[calo_slab]); 
          }
          if(!signal){
            hreco_intke_background->Fill(calo_slab_KE[calo_slab]);
          }
        }//<-- End if this is the interacting slab
      }//<--End calo slab loop

      //std::cout<<"unfolding debugging\n";
      //std::cout<<"\tnumber of true slabs: "<<inslabs<<std::endl;
      //std::cout<<"\tnumber of inc entries: "<<ninc_entries<<std::endl;
      //for(int calo_slab = 0; calo_slab < calo_slab_KE.size(); calo_slab++){
      //  if(calo_slab > calo_int_slab){continue;}
      //  if(!calo_slab_signal[calo_slab]){continue;}
      //  hreco_unfolding
      //}
      //for(int true_slab = 0; true_slab < inslabs; true_slab++){
      //  //std::cout<<"\t\ttrue pt: ("<<true_slab_xpos[true_slab]<<", " 
      //  //                           <<true_slab_ypos[true_slab]<<", "
      //  //                           <<true_slab_zpos[true_slab]<<")\n";
      //}
      //std::cout<<std::endl;

	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  // ~~~~~~~~~ End Code Block G ~~~~~~~~~~~~
	  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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


	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// ~~~~~~~~ Start Code Block H ~~~~~~~~~~~
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   // ## xs calculation ##
   for(int iBin = 0; iBin < hintke->GetNbinsX(); iBin++){
    if(h2incke->GetBinContent(iBin) == 0){continue;}
    double ratio   = 1.*hintke->GetBinContent(iBin) / h2incke->GetBinContent(iBin);  
    double temp_xs = ratio * recip_num_density; 
    double xs      = temp_xs / barn;

    double num_err = pow(hintke->GetBinContent(iBin), .5);
    double num = hintke->GetBinContent(iBin);
    if(num == 0){num = 1;}
    double term1 = num_err/num;
    double dem_err = pow(h2incke->GetBinContent(iBin), .5);
    double dem = h2incke->GetBinContent(iBin);
    if(dem == 0){dem =1;}
    double term2 = dem_err/dem;
    double totalError = temp_xs*pow(pow(term1,2) + pow(term2,2), 0.5)/barn;//*recip_num_density*barn;

    //std::cout <<"xs: " << xs << " +- " << totalError << std::endl;
    hxs->SetBinContent(iBin, xs); 
    hxs->SetBinError(iBin,totalError);
   }

   for(int iBin = 0; iBin < mintke->GetNbinsX(); iBin++){
    if(mincke->GetBinContent(iBin) == 0){continue;}
    double ratio   = 1.*mintke->GetBinContent(iBin) / mincke->GetBinContent(iBin);  
    double temp_xs = ratio * dense_recip_num_density; 
    double xs      = temp_xs / barn;

    double num_err = pow(mintke->GetBinContent(iBin), .5);
    double num = mintke->GetBinContent(iBin);
    if(num == 0){num = 1;}
    double term1 = num_err/num;
    double dem_err = pow(mincke->GetBinContent(iBin), .5);
    double dem = mincke->GetBinContent(iBin);
    if(dem == 0){dem =1;}
    double term2 = dem_err/dem;
    double totalError = temp_xs*pow(pow(term1,2) + pow(term2,2), 0.5)/barn;//*recip_num_density*barn;

    //std::cout <<"xs: " << xs << " +- " << totalError << std::endl;
    mxs->SetBinContent(iBin, xs); 
    mxs->SetBinError(iBin,totalError);
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



    std::cout<<"xs calculation using observables and corrections...\n";
    // ## folded signal distributions ##
    for(int iBin = 0; iBin < hreco_intke->GetNbinsX(); iBin++){
      double total = hreco_intke->GetBinContent(iBin);
      double background = hreco_intke_background->GetBinContent(iBin);
      hreco_folded_intke_signal->SetBinContent(iBin, total-background); 
    }
    for(int iBin = 0; iBin < hreco_incke->GetNbinsX(); iBin++){
      double total = hreco_incke->GetBinContent(iBin);
      double background = hreco_incke_background->GetBinContent(iBin);
      hreco_folded_incke_signal->SetBinContent(iBin, total-background); 
    }

    // ## unfolding matrix normalization ##
    for(int iBin = 0; iBin <hreco_incke->GetNbinsX(); iBin++){
      int column_total = 0;
      for(int jBin = 0; jBin < hreco_unfolding_matrix->GetNbinsX(); jBin++){
        column_total += hreco_unfolding_matrix->GetBinContent(iBin, jBin);
      }
      if(!column_total){continue;}
      for(int jBin = 0; jBin < hreco_unfolding_matrix_normalized->GetNbinsX(); jBin++){
        double normalized_value = hreco_unfolding_matrix->GetBinContent(iBin, jBin) / column_total;
        hreco_unfolding_matrix_normalized->SetBinContent(iBin, jBin, normalized_value);
      }
    }

    // ## unfolding signal distributions ##
    for(int iBin = 0; iBin < hreco_unfolding_matrix_normalized->GetNbinsX(); iBin++){
      int n_int_entries = hreco_folded_intke_signal->GetBinContent(iBin);
      int n_inc_entries = hreco_folded_incke_signal->GetBinContent(iBin);
      for(int jBin = 0; jBin < hreco_unfolding_matrix_normalized->GetNbinsY(); jBin++){
        double weight = hreco_unfolding_matrix_normalized->GetBinContent(iBin, jBin);
        double int_value = n_int_entries * weight;
        double inc_value = n_inc_entries * weight;
        hreco_unfolded_intke_signal->AddBinContent(jBin, int_value);
        hreco_unfolded_incke_signal->AddBinContent(jBin, inc_value);
      }
    }

    // ## epsilon curves ##
    for(int iBin = 0; iBin < hreco_intke->GetNbinsX(); iBin++){
      if(sintke->GetBinContent(iBin)){
        double int_eff = hreco_unfolded_intke_signal->GetBinContent(iBin) / sintke->GetBinContent(iBin);
        hreco_intke_eff->SetBinContent(iBin, int_eff);
      }
      if(sincke->GetBinContent(iBin)){
        double inc_eff = hreco_unfolded_incke_signal->GetBinContent(iBin) / sincke->GetBinContent(iBin);
        hreco_incke_eff->SetBinContent(iBin, inc_eff);
      }
    }

    // ## putting the pieces together ##
    for(int iBin = 0; iBin < hreco_unfolded_intke_signal->GetNbinsX(); iBin++){
      if(hreco_unfolded_incke_signal->GetBinContent(iBin) == 0){continue;}
      // num:   (N_int - Background_int)*U_ij * 1/eps
      double num = hreco_unfolded_intke_signal->GetBinContent(iBin) / hreco_intke_eff->GetBinContent(iBin);
      if(num == 0){num = 1;}
      // denom: (N_inc - Background_inc)*U_ij * 1/eps
      double dem = hreco_unfolded_incke_signal->GetBinContent(iBin) / hreco_incke_eff->GetBinContent(iBin);
      if(dem == 0){dem =1;}

      // # ratio #
      double ratio = num / dem;
      double temp_xs = ratio * sparse_recip_num_density;
      double xs = temp_xs / barn;
    
      // # error #
      double num_err = pow(num, .5);
      double term1 = num_err/num;
      double dem_err = pow(sincke->GetBinContent(iBin), .5);
      double term2 = dem_err/dem;
      double totalError = temp_xs*pow(pow(term1,2) + pow(term2,2), 0.5)/barn;//*recip_num_density*barn;

      //std::cout <<"xs: " << xs << " +- " << totalError << std::endl;
      hreco_xs->SetBinContent(iBin, xs); 
      hreco_xs->SetBinError(iBin,totalError);
    }

	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// ~~~~~~~~~ End Code Block H ~~~~~~~~~~~~
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



   // ## writing histograms ##
   TFile myfile("spacepoints.root", "RECREATE");
   hNumPts->Write();
   hDistanceBetweenPoints->Write();
   hdistvske->Write();
   hdedx->Write();
   hfirstx->Write();
   hfirsty->Write();
   hfirstz->Write();
   hintke->Write();
   hincke->Write();
   h2incke->Write();
   hxs->Write();
   mintke->Write();
   mincke->Write();
   mxs->Write();
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
   slab_vs_dense_intke->Write();
   ratio_int->Write();
   ratio_inc->Write();
   ratio_entries->Write();
   ahh->Write();
   RDSptAngle->Write();

   hreco_incke->Write();
   hreco_incke_signal->Write();
   hreco_folded_incke_signal->Write();
   hreco_unfolded_incke_signal->Write();
   hreco_incke_background->Write();
   hreco_intke->Write();
   hreco_intke_signal->Write();
   hreco_folded_intke_signal->Write();
   hreco_unfolded_intke_signal->Write();
   hreco_intke_background->Write();
   hreco_intke_eff->Write();
   hreco_incke_eff->Write();
   hreco_unfolding_matrix->Write();
   hreco_unfolding_matrix_normalized->Write();
   hreco_xs->Write();
}//<---End macro
