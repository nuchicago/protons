#include "PlotModule.h"

PlotModule::PlotModule(){}
PlotModule::PlotModule(TPostScript *psOutputFile){

  


  //char eventdisp_open[100];
  //sprintf(eventdisp_open,"%s(",psOutputFile);
  //c1->Print(eventdisp_open);

}

void PlotModule::CloseSummary(TPostScript *psOutputFile){
  //char eventdisp_close[100];
  //sprintf(eventdisp_close,"%s)",psOutputFile);
  //c1->Print(eventdisp_close);
  psOutputFile->Close();
  std::cout << "Closing PlotModule" << std::endl;
}

void PlotModule::EventSummary( TPostScript *psOutputFile,int run, int subrun, int event, int ntracks_reco,
                    std::vector<double> beamMatchInfo, double *interactionInfo, std::vector<int> ClusterIDvect,
                     double InitialKE,double IntKE, double beamMass,
                    std::vector< std::vector<double> > *track_xpos,
                    std::vector< std::vector<double> > *track_ypos,
                    std::vector< std::vector<double> > *track_zpos,
                    std::vector<int> *ntrack_hits,
                    std::vector< std::vector<double> > *col_track_x,
                    std::vector< std::vector<double> > *col_track_y,
                    std::vector< std::vector<double> > *col_track_z,
                    std::vector<int> *col_track_hits,
                    std::vector<std::vector<double>> *col_track_dedx,
                    std::vector<std::vector<double>> *col_track_rr,
                    double wctrk_XFace, double wctrk_YFace,
                    int nhits,
                    std::vector<double> *hit_time,
                    std::vector<double> *hit_amp,
                    std::vector<double> *hit_wire
                    ){



   std::cout << "Beginning Summary" << std::endl;

   psOutputFile->NewPage();
   gROOT->SetBatch(true);
   gStyle->SetPalette(53);

    TCanvas *c1 = new TCanvas("c1","Canvas",1920,1080);

    TPad * pad1 = new TPad("pad1", "Induction Plane",0.0,0.5,0.3,1.0,0);
    TPad * pad2 = new TPad("pad2", "XZ Tracks" ,0.3, 0.5, 0.6, 1.0, 0);
    TPad * pad3 = new TPad("pad3", "Residuals",0.6,0.6,0.85,1.0,0);
    TPad * pad4 = new TPad("pad4", "Relevant Figures",0.85,0.6,1.0,1.0,0);
    TPad * pad5 = new TPad("pad5", "collection Plane",0.0,0.0,0.3,0.5,0);
    TPad * pad6 = new TPad("pad6", "YZ Tracks",0.3,0.0,0.6,0.5,0);
    TPad * pad7 = new TPad("pad7", "Graph 3D",0.6,0.0,1.0,0.6,0);

    TGraph2D *wirePlotLims1 = new TGraph2D(2);
    wirePlotLims1->SetName("wirePlotLims1");
    wirePlotLims1->SetTitle("");
    wirePlotLims1->SetPoint(0,0,0,0);
    wirePlotLims1->SetPoint(1,240,3100,1000);
    wirePlotLims1->SetMarkerStyle(1);

    TGraph2D *wirePlotLims2 = new TGraph2D(2);
    wirePlotLims2->SetName("wirePlotLims2");
    wirePlotLims2->SetTitle("");
    wirePlotLims2->SetPoint(0,240,0,0);
    wirePlotLims2->SetPoint(1,480,3100,1000);
    wirePlotLims2->SetMarkerStyle(1);




  pad1->Draw();
  pad2->Draw();
  pad3->Draw();
  pad4->Draw();
  pad5->Draw();
  pad6->Draw();
  pad7->Draw();

   bool isMC = false;

  // stuff for the summary residual plot


  

   int reco_primary = static_cast <int>(beamMatchInfo[1]);
   int isInelastic =  static_cast<int>(interactionInfo[0]);
   double Int_x = interactionInfo[1];
   double Int_y = interactionInfo[2];
   double Int_z =  interactionInfo[3];
   double Int_type = interactionInfo[4];
   double selCrit = interactionInfo[5];

      double res_buffer = 0;
      int primary_pts = (*col_track_hits)[reco_primary];
      double dedx_graphpts[primary_pts], res_graphpts[primary_pts];


      if (primary_pts > 0){
        double max_dedx = -1;
        //calculating residual range
        for(int ipos = 0; ipos < primary_pts  ; ipos++){
          double res_range;
          
          res_graphpts[ipos] = (*col_track_rr)[reco_primary][ipos];
          dedx_graphpts[ipos] = (*col_track_dedx)[reco_primary][ipos];
          
          if (dedx_graphpts[ipos] > max_dedx){ max_dedx = dedx_graphpts[ipos];}
          }//end of loop for residual range

        //####### setting up graphs for raw wire hits
        std::cout << "TGraphs for raw hits" << std::endl;
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

        /*std::vector<double> ptracks_vx;
        ptracks_vx.reserve(primary_size);
        std::vector<double> ptracks_vy;
        ptracks_vy.reserve(primary_size);
        std::vector<double> ptracks_vz;
        ptracks_vz.reserve(primary_size);

          ptracks_vx.insert(ptracks_vx.end(), (*track_xpos)[reco_primary].begin(), (*track_xpos)[reco_primary].end());
          ptracks_vy.insert(ptracks_vy.end(), (*track_ypos)[reco_primary].begin(), (*track_ypos)[reco_primary].end());
          ptracks_vz.insert(ptracks_vz.end(), (*track_zpos)[reco_primary].begin(), (*track_zpos)[reco_primary].end());

        //std::cout << "ntracks_reco : "<< ntracks_reco << std::endl;

        */

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
        std::cout << "filled primary 2d and 3d" << std::endl;

        //std::cout << "primary : " << reco_primary << std::endl;
        //std::cout << "track 0 start XYZ:" << (*track_start_x)[0] << "\t" <<(*track_start_y)[0]<< "\t"  << (*track_start_z)[0] << std::endl;
        //std::cout << "track 0 start track_ipos: " << (*track_xpos)[0][0]<< "\t"  << (*track_ypos)[0][0]<< "\t"  << (*track_zpos)[0][0] << std::endl;

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
        TGraph *EventXZBranch = new TGraph( secondaries_size);
        EventXZBranch->SetName("EventXZBranch");
        TGraph *EventYZBranch = new TGraph( secondaries_size);
        EventYZBranch->SetName("EventYZBranch");
        
          
        


        TGraph2D *Event3dOther = new TGraph2D();
        TGraph2D *Event3dBranch = new TGraph2D();

        Event3dOther->SetName("Event3dOther");
        if (secondaries_size > 0){
          int graphPtBuffer = 0;
          int graphPtBufferBranch = 0; 
          for (int itrack = 0; itrack < ntracks_reco; itrack++ ){
            if (itrack != reco_primary){
              bool isBranch = false;
              if(ClusterIDvect.size() > 0){
                for (int j = 0; j < ClusterIDvect.size(); j++){
                  if(ClusterIDvect[j] == itrack ){isBranch = true;}
                }
              } 

              if(isBranch){
                  for(int ipoint = 0; ipoint < (*ntrack_hits)[itrack]; ipoint++){
                    Event3dBranch->SetPoint(graphPtBufferBranch,(*track_zpos)[itrack][ipoint],
                      (*track_xpos)[itrack][ipoint],(*track_ypos)[itrack][ipoint]);
                    EventYZBranch->SetPoint(graphPtBufferBranch,(*track_zpos)[itrack][ipoint],(*track_ypos)[itrack][ipoint]);
                    EventXZBranch->SetPoint(graphPtBufferBranch,(*track_zpos)[itrack][ipoint],(*track_xpos)[itrack][ipoint]);
                    graphPtBufferBranch++;
                  }
              }

              //std::cout << "which track" << itrack << std::endl;
              else{
                for(int ipoint = 0; ipoint < (*ntrack_hits)[itrack]; ipoint++){
                  //std::cout << ipoint << "point filled" << std::endl;
                  Event3dOther->SetPoint(graphPtBuffer,(*track_zpos)[itrack][ipoint],
                    (*track_xpos)[itrack][ipoint],(*track_ypos)[itrack][ipoint]);
                  EventYZother->SetPoint(graphPtBuffer,(*track_zpos)[itrack][ipoint],(*track_ypos)[itrack][ipoint]);
                  EventXZother->SetPoint(graphPtBuffer,(*track_zpos)[itrack][ipoint],(*track_xpos)[itrack][ipoint]);
                  graphPtBuffer++;
                }
              }
              //std::cout << "track num : " << itrack << std::endl;
              ////std::cout << "num points: " << (*ntrack_hits)[itrack] << std::endl;
            }
          }
        }
        std::cout << "Filled other tracks reco" << std::endl;

        // ######## legends and wc markers
        TLegend * Legend3d =  new TLegend();
        TLegend * LegendXZ =  new TLegend();
        TLegend * LegendYZ = new TLegend();
        TPolyMarker3D *wcPos3d = new TPolyMarker3D(1);
        TPolyMarker *wcPosYZ = new TPolyMarker(1);
        TPolyMarker *wcPosXZ = new TPolyMarker(1);

        if(!isMC){

          wcPos3d->SetPoint(0,0,wctrk_XFace,wctrk_YFace);
          wcPosYZ->SetPoint(0,0,wctrk_YFace);
          wcPosXZ->SetPoint(0,0,wctrk_XFace);
          wcPos3d->SetMarkerStyle(4);
          wcPos3d->SetMarkerColor(kOrange+1);
          wcPosYZ->SetMarkerStyle(4);
          wcPosYZ->SetMarkerColor(kOrange+1);
          wcPosXZ->SetMarkerStyle(4);
          wcPosXZ->SetMarkerColor(kOrange+1);
          //EventYZwc->SetPoint(0,0,wctrk_YFace);
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

        std::cout << "legends" << std::endl;

        // ######### using a 2 point Tgraph2D to set limits on wire displays
        pad1->cd();
        inductionHits->SetTitle("");
        inductionHits->GetXaxis()->SetTitle("Wire ");
        inductionHits->GetYaxis()->SetTitle("Time Tick");
        wirePlotLims1->Draw("P");
        inductionHits->Draw("PCOLZ SAME");
        std::cout << "finished drawing pad 1" << std::endl;
        inductionHits->GetHistogram()->SetMinimum(0);
        //inductionHits->GetHistogram()->SetMaximum(180);

        // ######## Setting limits on 2D reconstructed track plots
        std::cout << "2D reco plot formatting" << std::endl;

        EventYZprimary->GetXaxis()->SetLimits(-5,90);                 // along X
        EventYZprimary->GetHistogram()->SetMaximum(25.);   // along          
        EventYZprimary->GetHistogram()->SetMinimum(-25.);  //   Y     
        EventYZother->GetXaxis()->SetLimits(-5.,90);                 // along X
        EventYZother->GetHistogram()->SetMaximum(25.);   // along          
        EventYZother->GetHistogram()->SetMinimum(-25.);  //   Y     
        EventYZBranch->GetXaxis()->SetLimits(-5.,90);                 // along X
        EventYZBranch->GetHistogram()->SetMaximum(55);   // along          
        EventYZBranch->GetHistogram()->SetMinimum(-5);  //   Y    

        EventXZprimary->GetXaxis()->SetLimits(-5.,90);                 // along X
        EventXZprimary->GetHistogram()->SetMaximum(55);   // along          
        EventXZprimary->GetHistogram()->SetMinimum(-5);  //   Y     
        EventXZother->GetXaxis()->SetLimits(-5.,90);                 // along X
        EventXZother->GetHistogram()->SetMaximum(55);   // along          
        EventXZother->GetHistogram()->SetMinimum(-5);  //   Y     
        EventXZBranch->GetXaxis()->SetLimits(-5.,90);                 // along X
        EventXZBranch->GetHistogram()->SetMaximum(55);   // along          
        EventXZBranch->GetHistogram()->SetMinimum(-5);  //   Y    



        EventXZprimary->SetMarkerStyle(6);
        EventYZprimary->SetMarkerStyle(6);

        if (secondaries_size > 0){
          Legend3d->AddEntry(Event3dOther,"Pileup", "l");
          LegendXZ->AddEntry(EventXZother,"Pileup", "l");
          LegendYZ->AddEntry(EventYZother,"Pileup", "l");
        }
        if (ClusterIDvect.size() > 0){
          Legend3d->AddEntry(Event3dBranch,"Branches", "l");
          LegendXZ->AddEntry(EventXZBranch,"Branches", "l");
          LegendYZ->AddEntry(EventYZBranch,"Branches", "l");
        }

        pad1->SetTheta(90.); pad1->SetPhi(0.001);
        pad1->SetFrameFillColor(kAzure+7);

        

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
        pad5->SetFrameFillColor(kAzure+7);

        pad3->cd();
        TF1 *res_fx = new TF1("res_fx","[0]* x**([1])", 0, 90);
        res_fx->SetParameters(17, -0.42);
        res_fx->SetLineColor(4);

        TGraph *gres_dedx = new TGraph(primary_pts,res_graphpts,dedx_graphpts);
        if (max_dedx > 35) {gres_dedx->GetHistogram()->SetMaximum(max_dedx);}
        else{gres_dedx->GetHistogram()->SetMaximum(35);}
        gres_dedx->GetHistogram()->SetMinimum(0);
        gres_dedx->GetXaxis()->SetTitle("Residual Range [cm]");
        gres_dedx->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
        gres_dedx->SetTitle("");
        gres_dedx->SetMarkerStyle(9);
        //gres_dedx->SetMarkerSize(1);
        gres_dedx->SetMarkerColor(38);

        TLegend *res_legend = new TLegend(0.5,0.7,0.9,0.90);
        res_legend->SetNColumns(2);
        res_legend->AddEntry(res_fx,"#frac{dE}{dx} = A R^{b}","l");
        res_legend->AddEntry((TObject*)0, "A = 17", "");
        res_legend->AddEntry(gres_dedx,"Measured dE/dx","p");
        res_legend->AddEntry((TObject*)0, "b = -0.42", "");

        gres_dedx->Draw("AP");
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
        if(ClusterIDvect.size() > 0){
          Event3dBranch->SetMarkerStyle(6);
          Event3dBranch->SetMarkerColor(kViolet);
          Event3dBranch->SetLineColor(kViolet);
          Event3dBranch->Draw("P SAME");}


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
        if (ClusterIDvect.size() > 0){
          EventXZBranch->SetMarkerStyle(6);
          EventXZBranch->SetMarkerColor(kViolet);
          EventXZBranch->SetLineColor(kViolet);
          EventXZBranch->Draw("psame");}

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

        if (ClusterIDvect.size() > 0){
          EventYZBranch->SetMarkerStyle(6);
          EventYZBranch->SetMarkerColor(kViolet);
          EventYZBranch->SetLineColor(kViolet);
          EventYZBranch->Draw("psame");}
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

        
        char beamPtext[50];
        //sprintf(beamPtext,"Beam Momentum : %.02f MeV/c", (*wctrk_momentum));}
        char beamMassText[50];
        sprintf(beamMassText,"Mass : %.02f MeV/c^2", beamMass);


        char efieldtext[50];
        //sprintf(efieldtext,"E field: %.03f kV/cm", (*efield));
        char initialKEtext[50];
        sprintf(initialKEtext,"Initial KE: %.03f MeV", InitialKE);
        char intKEtext[50];
        sprintf(intKEtext,"Interaction KE: %.03f MeV", IntKE);
        char eventtext[50];
        sprintf(eventtext,"Data event: %d ", event);
        char topologytext[50];
        sprintf(topologytext, "Inelastic topology type: %f", Int_type);

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
        //pt->AddText(efieldtext);
        pt->AddText(beamMassText);
        if (Int_type > 0) {
          pt->AddText(topologytext);
          pt-> AddText(selCritText);}
        //pt->AddText(beamPtext);
        pt->AddText(initialKEtext);
        pt->AddText(intKEtext);
        pt->Draw();



        // gPad->Modified(); 
        //gPad->Update();
        std::cout << "Update Canvas" << std::endl;
        c1->Update();

        char eventdisp_title[100];
        sprintf(eventdisp_title,"%s/event%d.png",psOutputFile,event);
        
        //char gres_dedx_title[100];
        //sprintf(gres_dedx_title,"%s/dedx_Residuals/gres_dedx%d.png",UI->plotIndividual,event);

        

        //~gres_dedx;
        //c1->Print(eventdisp_title, "png");
        //delete gres_dedx;
        //delete inductionHits;
        //delete collectionHits;
        //delete EventXZother;
        //delete EventYZother;
        //delete EventXZprimary;
        //delete EventYZprimary;
        //delete Event3dPrimary;
        //delete Event3dOther;
        //delete Event3dTPC;
        //delete IntPoint;
          
        
        } //end if primary has > 0 pts
  

} //end loop()
