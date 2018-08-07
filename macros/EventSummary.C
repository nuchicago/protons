#define EventSummary_cxx
#include "EventSummary.h"
#include <TH2.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TPolyMarker3D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "Riostream.h"
#include "math.h"
#include "../src/Utilities/UtilityFunctions.h"



void EventSummary::Loop()
{
//   In a ROOT session, you can do:
//      root> .L EventSummary.C
//      root> EventSummary t
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

	TCanvas *csplit = new TCanvas("csplit","csplit",1500, 750);
  csplit->Divide(2,2);

	 bool isMC = false;

	// stuff for the summary residual plot
	std::vector<double> total_dedx_v;
  std::vector<double> total_res_v;
  TGraph gtotal_res_dedx_item;
  TGraph * gtotal_res_dedx;

  //Opening primaryID file

	 std::ifstream idfile("../files/PrimaryID.txt");
	 Long64_t fileEntry;
	 int reco_primary, isInelastic;
	 double Int_x, Int_y, Int_z;

	 // looping over event ids in file
	 while(idfile >> fileEntry >> reco_primary >> isInelastic >> Int_x >> Int_y >> Int_z){

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
					total_dedx_v.push_back(dedx_graphpts[ipos]);

				}







				int primary_size = (*ntrack_hits)[reco_primary];

				std::vector<double> ptracks_vx;
				ptracks_vx.reserve(primary_size);
				std::vector<double> ptracks_vy;
				ptracks_vy.reserve(primary_size);
				std::vector<double> ptracks_vz;
				ptracks_vz.reserve(primary_size);

					ptracks_vx.insert(ptracks_vx.end(), (*track_xpos)[reco_primary].begin(), (*track_xpos)[reco_primary].end());
					ptracks_vy.insert(ptracks_vy.end(), (*track_ypos)[reco_primary].begin(), (*track_ypos)[reco_primary].end());
					ptracks_vz.insert(ptracks_vz.end(), (*track_zpos)[reco_primary].begin(), (*track_zpos)[reco_primary].end());



				TGraph *inductionHits = new TGraph(primary_size,
					&ptracks_vz[0], &ptracks_vx[0]);
				inductionHits->SetName("inductionHits");

				TGraph *colectionHits = new TGraph(primary_size,
					&ptracks_vz[0], &ptracks_vy[0]);
				colectionHits->SetName("colectionHits");


				int indGraphEntry = 0;
				int colGraphEntry = 0;
				for (int ihit = 0; ihit < nhits; ihit++){
					if((*hit_wire)[ihit] < 240){
						inductionHits->SetPoint(indGraphEntry,(*hit_wire)[ihit],(*hit_time)[ihit]);
						indGraphEntry++;
					}
					else{
						colectionHits->SetPoint(colGraphEntry,(*hit_wire)[ihit],(*hit_time)[ihit]);
						colGraphEntry++;
					}
				}


				TGraph2D *Event3dPrimary =  new TGraph2D(primary_size,&ptracks_vz[0],&ptracks_vx[0],&ptracks_vy[0]);
				Event3dPrimary->SetName("Event3dPrimary");

				



				int secondaries_size = 0;
				for (int itrack = 0; itrack < ntracks_reco && itrack !=reco_primary; itrack++ ){
					
					int ctrack_size = (*ntrack_hits)[itrack];
					secondaries_size += ctrack_size;
				}

				std::vector<double> otracks_vx;
				otracks_vx.reserve(secondaries_size);
				std::vector<double> otracks_vy;
				otracks_vy.reserve(secondaries_size);
				std::vector<double> otracks_vz;
				otracks_vz.reserve(secondaries_size);


				TGraph2D *Event3dOther = new TGraph2D();
				Event3dOther->SetName("Event3dOther");
				if (secondaries_size > 0){
					int graphPtBuffer = 0;
					for (int itrack = 0; itrack < ntracks_reco && itrack !=reco_primary; itrack++ ){
						otracks_vx.insert(otracks_vx.end(), (*track_xpos)[itrack].begin(), (*track_xpos)[itrack].end());
						otracks_vy.insert(otracks_vy.end(), (*track_ypos)[itrack].begin(), (*track_ypos)[itrack].end());
						otracks_vz.insert(otracks_vz.end(), (*track_zpos)[itrack].begin(), (*track_zpos)[itrack].end());
						
						for(int ipoint = 0; ipoint < (*ntrack_hits)[itrack]; ipoint++){
							//std::cout << ipoint <<"point filled" << std::endl;
							Event3dOther->SetPoint(graphPtBuffer,(*track_zpos)[itrack][ipoint],
								(*track_xpos)[itrack][ipoint],(*track_ypos)[itrack][ipoint]);
							graphPtBuffer++;
						}
					}
				}


				/*
				TGraph *EventXZother = new TGraph( secondaries_size, &otracks_vz[0], &otracks_vx[0]);
				EventXZother->SetName("EventXZother");
				TGraph *EventYZother = new TGraph( secondaries_size, &otracks_vz[0], &otracks_vy[0]);
				EventYZother->SetName("EventYZother");
				TGraph *EventXZwc = new TGraph(1);
				EventXZwc->SetName("EventXZwc");
				TGraph *EventYZwc = new TGraph(1);
				EventYZwc->SetName("EventYZwc"); */

				

				
				//TLegend * Legend2d =  new TLegend(0.15,0.85,0.45,1);
				TLegend * Legend3d =  new TLegend();

				if(!isMC){
					//EventXZwc->SetPoint(0,0,wctrk_XFace[0]);
					//EventYZwc->SetPoint(0,0,wctrk_YFace[0]);
					//Legend2d->AddEntry(//EventXZwc,"Wire Chamber Projection","p");
					
				}

				inductionHits->SetMarkerStyle(7);
				colectionHits->SetMarkerStyle(7);
				if(isInelastic){
					/*inductionHits->SetMarkerColor(3);
					colectionHits->SetMarkerColor(3);
					inductionHits->SetLineColor(3);
					colectionHits->SetLineColor(3);*/
					//Legend2d->AddEntry(colectionHits,"Primary (Inelastic)","l");
					Legend3d->AddEntry(Event3dPrimary,"Primary (Inelastic)", "l");
				}
				else{
					/*inductionHits->SetMarkerColor(2);
					colectionHits->SetMarkerColor(2);
					inductionHits->SetLineColor(2);
					colectionHits->SetLineColor(2);*/
					//Legend2d->AddEntry(colectionHits,"Primary (No Inelastic)","l");
					Legend3d->AddEntry(Event3dPrimary,"Primary (No Inelastic)", "l");
				}
				/*
				EventXZother->SetMarkerStyle(7);
				EventXZother->SetMarkerColor(4);
				EventXZother->SetLineColor(4);

				EventYZother->SetMarkerStyle(7);
				EventYZother->SetMarkerColor(4);
				EventYZother->SetLineColor(4);
				

				//EventXZwc->SetMarkerStyle(9);
				//EventXZwc->SetMarkerColor(kOrange);
				EventYZwc->SetMarkerStyle(9);
				EventYZwc->SetMarkerColor(kOrange); */

				colectionHits->GetXaxis()->SetLimits(240,480);                 // along X
				colectionHits->GetHistogram()->SetMaximum(3100.);   // along          
				colectionHits->GetHistogram()->SetMinimum(0);  //   Y     
				/*
				EventYZother->GetXaxis()->SetLimits(-5.,90);                 // along X
				EventYZother->GetHistogram()->SetMaximum(25.);   // along          
				EventYZother->GetHistogram()->SetMinimum(-25.);  //   Y     */

				inductionHits->GetXaxis()->SetLimits(0,240);                 // along X
				inductionHits->GetHistogram()->SetMaximum(3100);   // along          
				inductionHits->GetHistogram()->SetMinimum(0);  //   Y     

				/*
				EventXZother->GetXaxis()->SetLimits(-5.,90);                 // along X
				EventXZother->GetHistogram()->SetMaximum(55);   // along          
				EventXZother->GetHistogram()->SetMinimum(-5);  //   Y      */

				
				
				csplit->cd(1);
				inductionHits->SetTitle("");
				inductionHits->GetXaxis()->SetTitle("Wire ");
				inductionHits->GetYaxis()->SetTitle("Time Tick");
				if (primary_size > 0){inductionHits->Draw("AP");}
				//Legend2d->Draw();}
				if (secondaries_size > 0){//EventXZother->Draw("psame");
					//Legend2d->AddEntry(EventYZother,"Secondaries","l");
					Legend3d->AddEntry(Event3dOther,"Secondaries", "l");}
				if (!isMC){//EventXZwc->Draw("psame");}
				}

				csplit->cd(3);
				colectionHits->SetTitle("");
				colectionHits->GetXaxis()->SetTitle("Wire");
				colectionHits->GetYaxis()->SetTitle("Time Tick");
				if (primary_size > 0){colectionHits->Draw("AP");
						}
				if (secondaries_size > 0){//EventYZother->Draw("psame");}
				}
				if (!isMC){//EventYZwc->Draw("psame");}
				}


				csplit->cd(2);
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
				gres_dedx->SetMarkerStyle(20);
				gres_dedx->SetMarkerSize(1.5);
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
				csplit->cd(4);
				
				TGraph2D *Event3dTPC = new TGraph2D();
				Event3dTPC->SetName("Event3dTPC");
				Event3dTPC->SetPoint(0,-1,0,-20);
				Event3dTPC->SetPoint(1,90,48,20);
				Event3dTPC->SetTitle("3D Event Display");
				Event3dTPC->GetXaxis()->SetTitle("Z [cm]");
				Event3dTPC->GetYaxis()->SetTitle("X [cm]");
				Event3dTPC->GetZaxis()->SetTitle("Y [cm]");
				Event3dTPC->Draw("P");
				Event3dPrimary->Draw("P SAME");

				TPolyMarker3D *IntPoint = new TPolyMarker3D(1);
				IntPoint->SetPoint(0, Int_z, Int_x, Int_y);
				
				Legend3d->Draw();
				
				if(secondaries_size > 0){
					Event3dOther->SetMarkerStyle(7);
					Event3dOther->SetMarkerColor(4);
					Event3dOther->SetLineColor(4);

					Event3dOther->Draw("P SAME");}


				Event3dPrimary->SetMarkerStyle(7);
				if(isInelastic){
					Event3dPrimary->SetMarkerColor(3);
					Event3dPrimary->SetLineColor(3);
					Legend3d->AddEntry(IntPoint,"interaction point", "p");
					IntPoint->Draw("P SAME");
					IntPoint->SetName("IntPoint");
					IntPoint->SetMarkerStyle(2);

				}
				else{Event3dPrimary->SetMarkerColor(2);
					Event3dPrimary->SetLineColor(2);}
				

				// gPad->Modified(); 
				gPad->Update();

				char eventdisp_title[100];
				sprintf(eventdisp_title,"../plotting/images/EventSummary/EventSummary%d.png",event);
				
				//char gres_dedx_title[100];
				//sprintf(gres_dedx_title,"%s/dedx_Residuals/gres_dedx%d.png",UI->plotIndividual,event);


				//~gres_dedx;
				csplit->Print(eventdisp_title,"png");
				delete gres_dedx;
				delete inductionHits;
				delete colectionHits;
				/*delete EventXZother;
				delete EventYZother;
				delete EventYZwc;
				delete //EventXZwc;*/
				delete Event3dPrimary;
				delete Event3dOther;
				delete Event3dTPC;
				delete IntPoint;
					
				
				} //end if primary has > 0 pts
	 } //end while

	 idfile.close();
} //end loop()
