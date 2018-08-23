#define EventSummaryMC_cxx
#include "EventSummaryMC.h"
#include <TH2.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TPolyMarker3D.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include "Riostream.h"
#include "math.h"
#include "../src/Utilities/UtilityFunctions.h"



void EventSummaryMC::Loop()
{
//   In a ROOT session, you can do:
//      root> .L EventSummaryMC.C
//      root> EventSummaryMC t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   gROOT->SetBatch(true);
   gStyle->SetPalette(53);

    TGraph2D *wirePlotLims1 = new TGraph2D(2);
    wirePlotLims1->SetName("wirePlotLims1");
    wirePlotLims1->SetTitle("");
    wirePlotLims1->SetPoint(0,0,0,0);
    wirePlotLims1->SetPoint(1,240,3100,200);
    wirePlotLims1->SetMarkerStyle(1);

    TGraph2D *wirePlotLims2 = new TGraph2D(2);
    wirePlotLims2->SetName("wirePlotLims2");
    wirePlotLims2->SetTitle("");
    wirePlotLims2->SetPoint(0,240,0,0);
    wirePlotLims2->SetPoint(1,480,3100,200);
    wirePlotLims2->SetMarkerStyle(1);

  TCanvas *c1 = new TCanvas("c1","Canvas",1920,1080);
  


  TPad *pad1 = new TPad("pad1", "Induction Plane",0.0,0.5,0.3,1.0,0);
  TPad *pad2 = new TPad("pad2", "XZ Tracks" ,0.3, 0.5, 0.6, 1.0, 0);
  TPad *pad3 = new TPad("pad3", "Residuals",0.6,0.6,0.85,1.0,0);
  TPad *pad4 = new TPad("pad4", "Relevant Figures",0.85,0.6,1.0,1.0,0);
  TPad *pad5 = new TPad("pad5", "collection Plane",0.0,0.0,0.3,0.5,0);
  TPad *pad6 = new TPad("pad6", "YZ Tracks",0.3,0.0,0.6,0.5,0);
  TPad *pad7 = new TPad("pad7", "Graph 3D",0.6,0.0,1.0,0.6,0);

  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad4->Draw();
  pad5->Draw();
  pad6->Draw();
  pad7->Draw();

   bool isMC = true;

  // stuff for the summary residual plot
  std::vector<double> total_dedx_v;
  std::vector<double> total_res_v;
  TGraph gtotal_res_dedx_item;
  TGraph * gtotal_res_dedx;

  //Opening primaryID file

   std::ifstream idfile("../files/PrimaryID_MC.txt");
   Long64_t fileEntry;
   int reco_primary, isInelastic, Int_type;
   double Int_x, Int_y, Int_z, selCrit, InitialKE, IntKE;

   // looping over event ids in file
   while(idfile >> fileEntry >> reco_primary >> isInelastic >> Int_x >> Int_y >> Int_z >> Int_type >> selCrit >> InitialKE >> IntKE){

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   Long64_t jentry = fileEntry;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;


      double res_buffer = 0;
      int primary_pts = (*col_track_hits)[reco_primary];
      double dedx_graphpts[primary_pts], res_graphpts[primary_pts];


      if (primary_pts > 0){
        double max_dedx = -1;
        //calculating residual range
        for(int ipos = primary_pts - 1; ipos >= 0  ; ipos--){
          double res_range;
          if (ipos == primary_pts -1){res_range = 0;}
          else{
           res_range =  pointDistance(
            (*col_track_x)[reco_primary][ipos],
            (*col_track_y)[reco_primary][ipos],
            (*col_track_z)[reco_primary][ipos],
            (*col_track_x)[reco_primary][ipos + 1],
            (*col_track_y)[reco_primary][ipos + 1],
            (*col_track_z)[reco_primary][ipos + 1]);}

          res_buffer += res_range;
          res_graphpts[ipos] = res_buffer;
          dedx_graphpts[ipos] = (*col_track_dedx)[reco_primary][ipos];
          total_res_v.push_back(res_buffer);
          if (dedx_graphpts[ipos] > max_dedx){ max_dedx = dedx_graphpts[ipos];}
          total_dedx_v.push_back(dedx_graphpts[ipos]);}//end of loop for residual range

        //####### setting up graphs for raw wire hits

        TGraph2D *inductionHits = new TGraph2D(nhits);
        inductionHits->SetName("inductionHits");

        TGraph2D *collectionHits = new TGraph2D(nhits);
        collectionHits->SetName("collectionHits");

        int indGraphEntry = 0;
        int colGraphEntry = 0;
        for (int ihit = 0; ihit < nhits; ihit++){
          if((*hit_wire)[ihit] < 240){
            inductionHits->SetPoint(indGraphEntry,(*hit_wire)[ihit],(*hit_time)[ihit],(*hit_amp)[ihit]);
            indGraphEntry++;
          }
          else{
            collectionHits->SetPoint(colGraphEntry,(*hit_wire)[ihit],(*hit_time)[ihit],(*hit_amp)[ihit]);
            colGraphEntry++;
          }
        }

          //####### creating a TGraph2d for reconstructed track

        int primary_size = (*ntrack_hits)[reco_primary];

        /*
        std::vector<double> ptracks_vx;
        ptracks_vx.reserve(primary_size);
        std::vector<double> ptracks_vy;
        ptracks_vy.reserve(primary_size);
        std::vector<double> ptracks_vz;
        ptracks_vz.reserve(primary_size);
          ptracks_vx.insert(ptracks_vx.end(), (*track_xpos)[reco_primary].begin(), (*track_xpos)[reco_primary].end());
          ptracks_vy.insert(ptracks_vy.end(), (*track_ypos)[reco_primary].begin(), (*track_ypos)[reco_primary].end());
          ptracks_vz.insert(ptracks_vz.end(), (*track_zpos)[reco_primary].begin(), (*track_zpos)[reco_primary].end());

        */

        //std::cout << "ntracks_reco : "<< ntracks_reco << std::endl;
        TGraph2D *Event3dPrimary =  new TGraph2D(primary_size);

        TGraph *EventYZprimary  = new TGraph(primary_size);
        TGraph *EventXZprimary = new TGraph(primary_size);
        Event3dPrimary->SetName("Event3dPrimary");
        EventYZprimary->SetName("EventYZprimary");
        EventXZprimary->SetName("EventXZprimary");

        for (int ipoint = 0 ; ipoint < primary_size;  ipoint++){
          Event3dPrimary->SetPoint(ipoint, (*track_zpos)[reco_primary][ipoint],
            (*track_xpos)[reco_primary][ipoint],(*track_ypos)[reco_primary][ipoint]);
          EventXZprimary->SetPoint(ipoint, (*track_zpos)[reco_primary][ipoint],
            (*track_xpos)[reco_primary][ipoint]);
          EventYZprimary->SetPoint(ipoint, (*track_zpos)[reco_primary][ipoint],
            (*track_ypos)[reco_primary][ipoint]);
        }

        //std::cout << "primary : " << reco_primary << std::endl;
        //std::cout << "track 0 size:" << (*ntrack_hits)[0] << std::endl;
        //std::cout << "track 1 size: " << (*ntrack_hits)[1] << std::endl;

        int secondaries_size = 0;
        for (int itrack = 0; itrack < ntracks_reco; itrack++ ){
          if (itrack != reco_primary){
            int ctrack_size = (*ntrack_hits)[itrack];
            secondaries_size += ctrack_size;
          }
        }

        //std::cout << "secondaries_size: " << secondaries_size << std::endl;
        //####### creating a TGraph2d for all other reconstructed tracks

        TGraph *EventXZother = new TGraph( secondaries_size);
        EventXZother->SetName("EventXZother");
        TGraph *EventYZother = new TGraph( secondaries_size);
        EventYZother->SetName("EventYZother");
        
          
        


        TGraph2D *Event3dOther = new TGraph2D();
        Event3dOther->SetName("Event3dOther");
        if (secondaries_size > 0){
          int graphPtBuffer = 0;
          for (int itrack = 0; itrack < ntracks_reco; itrack++ ){
            if (itrack != reco_primary){
              //std::cout << "which track" << itrack << std::endl;
              
              for(int ipoint = 0; ipoint < (*ntrack_hits)[itrack]; ipoint++){
                //std::cout << ipoint << "point filled" << std::endl;
                Event3dOther->SetPoint(graphPtBuffer,(*track_zpos)[itrack][ipoint],
                  (*track_xpos)[itrack][ipoint],(*track_ypos)[itrack][ipoint]);
                EventYZother->SetPoint(graphPtBuffer,(*track_zpos)[itrack][ipoint],(*track_ypos)[itrack][ipoint]);
                EventXZother->SetPoint(graphPtBuffer,(*track_zpos)[itrack][ipoint],(*track_xpos)[itrack][ipoint]);
                graphPtBuffer++;
              }
              //std::cout << "track num : " << itrack << std::endl;
              //std::cout << "num points: " << (*ntrack_hits)[itrack] << std::endl;
            }
          }
        }


        // ######## legends and wc markers
        TLegend * Legend3d =  new TLegend();
        TLegend * LegendXZ =  new TLegend();
        TLegend * LegendYZ = new TLegend();
        TPolyMarker3D *wcPos3d = new TPolyMarker3D(1);
        TPolyMarker *wcPosYZ = new TPolyMarker(1);
        TPolyMarker *wcPosXZ = new TPolyMarker(1);

        if(!isMC){

          wcPos3d->SetPoint(0,0,wctrk_XFace[0],wctrk_YFace[0]);
          wcPosYZ->SetPoint(0,0,wctrk_YFace[0]);
          wcPosXZ->SetPoint(0,0,wctrk_XFace[0]);
          wcPos3d->SetMarkerStyle(4);
          wcPos3d->SetMarkerColor(kOrange+1);
          wcPosYZ->SetMarkerStyle(4);
          wcPosYZ->SetMarkerColor(kOrange+1);
          wcPosXZ->SetMarkerStyle(4);
          wcPosXZ->SetMarkerColor(kOrange+1);
          //EventYZwc->SetPoint(0,0,wctrk_YFace[0]);
          Legend3d->AddEntry(wcPos3d,"Beam Position","p");
          LegendYZ->AddEntry(wcPosYZ,"Beam Position","p");
          LegendXZ->AddEntry(wcPosXZ,"Beam Position","p");
          
        }

        inductionHits->SetMarkerStyle(7);
        collectionHits->SetMarkerStyle(7);
        if(isInelastic){
          Legend3d->AddEntry(Event3dPrimary,"Primary (Inelastic)", "l");
          LegendYZ->AddEntry(EventYZprimary,"Primary (Inelastic)", "l");
          LegendXZ->AddEntry(EventXZprimary,"Primary (Inelastic)", "l");

        }
        else{
          Legend3d->AddEntry(Event3dPrimary,"Primary (No Inelastic)", "l");
          LegendYZ->AddEntry(EventYZprimary,"Primary (No Inelastic)", "l");
          LegendXZ->AddEntry(EventXZprimary,"Primary (No Inelastic)", "l");
        }

        //using a 2 point Tgraph2D to set limits on wire displays
        pad1->cd();
        inductionHits->SetTitle("");
        wirePlotLims1->GetXaxis()->SetTitle("Wire ");
        wirePlotLims1->GetYaxis()->SetTitle("Time Tick");
        wirePlotLims1->Draw("P");
        inductionHits->Draw("PCOLZ SAME");
        inductionHits->GetHistogram()->SetMinimum(0);
        //inductionHits->GetHistogram()->SetMaximum(180);

        // Setting limits on 2D reconstructed track plots
        EventYZprimary->GetXaxis()->SetLimits(-5,90);                 // along X
        EventYZprimary->GetHistogram()->SetMaximum(25.);   // along          
        EventYZprimary->GetHistogram()->SetMinimum(-25.);  //   Y     
        EventYZother->GetXaxis()->SetLimits(-5.,90);                 // along X
        EventYZother->GetHistogram()->SetMaximum(25.);   // along          
        EventYZother->GetHistogram()->SetMinimum(-25.);  //   Y     

        EventXZprimary->GetXaxis()->SetLimits(-5.,90);                 // along X
        EventXZprimary->GetHistogram()->SetMaximum(55);   // along          
        EventXZprimary->GetHistogram()->SetMinimum(-5);  //   Y     
        EventXZother->GetXaxis()->SetLimits(-5.,90);                 // along X
        EventXZother->GetHistogram()->SetMaximum(55);   // along          
        EventXZother->GetHistogram()->SetMinimum(-5);  //   Y     

        EventXZprimary->SetMarkerStyle(6);
        EventYZprimary->SetMarkerStyle(6);
        

        if (secondaries_size > 0){
          Legend3d->AddEntry(Event3dOther,"Secondaries", "l");
          LegendXZ->AddEntry(EventXZother,"Secondaries", "l");
          LegendYZ->AddEntry(EventYZother,"Secondaries", "l");
        }

        pad1->SetTheta(90.); pad1->SetPhi(0.001);
        pad1->SetFrameFillColor(kAzure+10);

        pad5->cd();
        collectionHits->SetTitle("");
        wirePlotLims2->GetXaxis()->SetTitle("Wire");
        wirePlotLims2->GetYaxis()->SetTitle("Time Tick");
        wirePlotLims2->Draw("P");

        collectionHits->Draw("PCOLZ SAME");
        //gPad->Update();
        //TPaletteAxis *colPalette = 
        //(TPaletteAxis*)collectionHits->GetHistogram->GetListOfFunctions()->FindObject("palette");
        //colPalette->GetAxis()->SetRangeUser(0,200);
        //gPad->Update();
        collectionHits->GetHistogram()->SetMinimum(0);
        //collectionHits->GetHistogram()->SetMaximum(180);
            
        if (secondaries_size > 0){//EventYZother->Draw("psame");}
        }
        if (!isMC){//EventYZwc->Draw("psame");}
        }
        pad5->SetTheta(90.); pad5->SetPhi(0.001);
        pad5->SetFrameFillColor(kAzure+10);

        pad3->cd();
        TF1 *res_fx = new TF1("res_fx","[0]* x**([1])", 0, 90);
        res_fx->SetParameters(17, -0.42);
        res_fx->SetLineColor(4);

        TGraph *gres_dedx = new TGraph(primary_pts,res_graphpts,dedx_graphpts);
        TGraph *gcol_dedx = new TGraph(primary_pts);

        for (int tpoint = 0 ; tpoint < (*ntrack_hits)[reco_primary]; tpoint++){
          gcol_dedx->SetPoint(tpoint,(*col_track_rr)[reco_primary][tpoint] / 1000, (*col_track_dedx)[reco_primary][tpoint]);
        }

        if (max_dedx > 35) {gres_dedx->GetHistogram()->SetMaximum(max_dedx);}
        else{gres_dedx->GetHistogram()->SetMaximum(35);}
        gres_dedx->GetHistogram()->SetMinimum(0);
        gres_dedx->GetXaxis()->SetTitle("Residual Range [cm]");
        gres_dedx->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
        gres_dedx->SetTitle("");
        gres_dedx->SetMarkerStyle(9);
        gcol_dedx->SetMarkerStyle(9);
        gcol_dedx->SetMarkerColor(2);

        //gres_dedx->SetMarkerSize(1);
        gres_dedx->SetMarkerColor(38);

        TLegend *res_legend = new TLegend(0.5,0.7,0.9,0.90);
        res_legend->SetNColumns(2);
        res_legend->AddEntry(res_fx,"#frac{dE}{dx} = A R^{b}","l");
        res_legend->AddEntry((TObject*)0, "A = 17", "");
        res_legend->AddEntry(gres_dedx,"Measured dE/dx","p");
        res_legend->AddEntry((TObject*)0, "b = -0.42", "");

        gres_dedx->Draw("AP");
        //gcol_dedx->Draw("psame");
        res_fx->Draw("lsame");
        res_legend->Draw();

        pad7->cd();
        
        TGraph2D *Event3dTPC = new TGraph2D();
        Event3dTPC->SetName("Event3dTPC");
        Event3dTPC->SetPoint(0,-1,0,-20);
        Event3dTPC->SetPoint(1,90,48,20);
        Event3dTPC->SetTitle("Reconstructed Track");
        Event3dTPC->GetXaxis()->SetTitle("Z [cm]");
        Event3dTPC->GetYaxis()->SetTitle("X [cm]");
        Event3dTPC->GetZaxis()->SetTitle("Y [cm]");
        Event3dTPC->Draw("P");
        Event3dPrimary->Draw("P SAME");
        if (!isMC){wcPos3d->Draw("P0 SAME");}

        TPolyMarker3D *IntPoint = new TPolyMarker3D(1);
        IntPoint->SetPoint(0, Int_z, Int_x, Int_y);
        TPolyMarker *IntPointXZ = new TPolyMarker(1);
        IntPointXZ->SetPoint(0, Int_z, Int_x);
        TPolyMarker *IntPointYZ = new TPolyMarker(1);
        IntPointYZ->SetPoint(0, Int_z, Int_y);
        
        Legend3d->Draw();
        
        if(secondaries_size > 0){
          Event3dOther->SetMarkerStyle(6);
          Event3dOther->SetMarkerColor(4);
          Event3dOther->SetLineColor(4);

          Event3dOther->Draw("P SAME");}


        Event3dPrimary->SetMarkerStyle(6);
        if(isInelastic){
          Event3dPrimary->SetMarkerColor(8);
          Event3dPrimary->SetLineColor(8);
          Legend3d->AddEntry(IntPoint,"interaction point", "p");
          IntPoint->Draw("P SAME");
          IntPoint->SetName("IntPoint");
          IntPoint->SetMarkerStyle(2);



        }
        else{
          Event3dPrimary->SetMarkerColor(2);
          Event3dPrimary->SetLineColor(2);
          EventXZprimary->SetMarkerColor(2);
          EventXZprimary->SetLineColor(2);
          EventYZprimary->SetMarkerColor(2);
          EventYZprimary->SetLineColor(2);
        }
        

        //####### 2D reco track plots

        pad2->cd();


        EventXZprimary->SetTitle("");
        EventXZprimary->GetXaxis()->SetTitle("Z [cm]");
        EventXZprimary->GetYaxis()->SetTitle("X [cm]");
        EventXZprimary->Draw("AP");
        LegendXZ->Draw();

        if (secondaries_size > 0){
          EventXZother->SetMarkerStyle(6);
          EventXZother->SetMarkerColor(4);
          EventXZother->SetLineColor(4);
          EventXZother->Draw("psame");}
        if (!isMC){wcPosXZ->Draw("psame");}

        if(isInelastic){
          EventXZprimary->SetMarkerColor(8);
          EventXZprimary->SetLineColor(8);
          LegendXZ->AddEntry(IntPointXZ,"interaction point", "p");
          IntPointXZ->Draw("P SAME");
          //IntPointXZ->SetName("IntPointXZ");
          IntPointXZ->SetMarkerStyle(2);
        }


        pad6->cd();

        EventYZprimary->SetTitle("");
        EventYZprimary->GetXaxis()->SetTitle("Z [cm]");
        EventYZprimary->GetYaxis()->SetTitle("Y [cm]");
        EventYZprimary->Draw("AP");
        LegendYZ->Draw();

        if (secondaries_size > 0){
          EventYZother->SetMarkerStyle(6);
          EventYZother->SetMarkerColor(4);
          EventYZother->SetLineColor(4);

          EventYZother->Draw("psame");}
        if (!isMC){wcPosYZ->Draw("psame");}
        

        if(isInelastic){

          EventYZprimary->SetMarkerColor(8);
          EventYZprimary->SetLineColor(8);
          LegendYZ->AddEntry(IntPointYZ,"interaction point", "p");
          IntPointYZ->Draw("P SAME");
          //IntPointYZ->SetName("IntPointYZ");
          IntPointYZ->SetMarkerStyle(2);


        }



        pad4->cd();

        

        TPaveText *pt = new TPaveText(.05,.1,.95,.8);

        if(!isMC){
        char beamPtext[50];
        sprintf(beamPtext,"Beam Momentum : %.02f MeV/c", (*wctrk_momentum));}

        char efieldtext[50];
        sprintf(efieldtext,"E field: %.03f kV/cm", (*efield));
        char initialKEtext[50];
        sprintf(initialKEtext,"Initial KE: %.03f MeV", InitialKE);
        char intKEtext[50];
        sprintf(intKEtext,"Interaction KE: %.03f MeV", IntKE);
        char eventtext[50];
        sprintf(eventtext,"MC event: %lld ", jentry);
        char topologytext[50];
        sprintf(topologytext, "Inelastic topology type: %d", Int_type);

        char selCritText[50];

        if(Int_type == 1 ){
          sprintf(selCritText, "Kink Angle = %0.3f deg", selCrit);
        }
        else if(Int_type == 2 ){
          sprintf(selCritText, "Average End dE/dx = %0.3f MeV", selCrit);
        }
        else{sprintf(selCritText, "Number of branches = %0.f", selCrit);

        }



        
        pt->AddText(eventtext);
        pt->AddText(efieldtext);
        if (Int_type > 0) {pt->AddText(topologytext);
        pt->AddText(selCritText);}
        //pt->AddText(beamPtext);

        pt->AddText(initialKEtext);
        pt->AddText(intKEtext);

        pt->Draw();

        // gPad->Modified(); 
        gPad->Update();
        c1->Update();

        char eventdisp_title[100];
        sprintf(eventdisp_title,"../plotting/images/EventSummary/EventSummaryMC%lld.png",jentry);
        
        //char gres_dedx_title[100];
        //sprintf(gres_dedx_title,"%s/dedx_Residuals/gres_dedx%d.png",UI->plotIndividual,event);


        //~gres_dedx;
        c1->Print(eventdisp_title,"png");
        delete gres_dedx;
        delete inductionHits;
        delete collectionHits;
        delete EventXZother;
        delete EventYZother;
        delete EventXZprimary;
        delete EventYZprimary;
        delete Event3dPrimary;
        delete Event3dOther;
        delete Event3dTPC;
        delete IntPoint;
          
        
        } //end if primary has > 0 pts
   } //end while

   idfile.close();
} //end loop()
