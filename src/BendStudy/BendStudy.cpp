//#############################################################################
//
// BendStudy.cpp
//
//  Measures the track bending issue in LArIAT data
//
// January 17, 2018
//
//#############################################################################

#include "BendStudy.h"

#include "../Selection/EventSelector.h"
#include "../Selection/BeamSelector.h"



//=============================================================================
// MAIN
//=============================================================================

int main( int argc, char * argv[] ){

  char *jobOptionsFile;

  if (argc == 1 )
  {
    jobOptionsFile    = new char[200];

    cout << endl;
    cout << "Welcome to your favorite LArIAT analysis macro...it's mine too." << endl;
    cout << "First, a question to get to know you and your analysis needs a little better." << endl;
    
    cout << endl;
    cout << "Where is your jobOptions file (full path) : ";
    cin >> jobOptionsFile;
  }
  else if( argc == 2 )
  {
    jobOptionsFile   = argv[1];
  }
  else
  {
    cout << "Program confused by user's input...exiting" << endl;
    exit(0);
  }

  //--------------------------------------------------------------------------
  // Run ROOT in batch mode
  //--------------------------------------------------------------------------

  TROOT theRoot("theRoot","root for XS Analysis");
  gROOT->SetBatch(true);


  //--------------------------------------------------------------------------
  // ROOT Style Settings
  //--------------------------------------------------------------------------
  
  //gStyle->SetOptStat(0000);
  //gStyle->SetOptFit(0000);
  //gStyle->SetOptTitle(0);
  //gStyle->SetFillColor(0);
  //gStyle->SetPadColor(0);
  //gStyle->SetCanvasColor(0);
  //gStyle->SetStatColor(0);
  //gStyle->SetTitleColor(0);
  //gStyle->SetPadBorderMode(0);
  //gStyle->SetFrameBorderMode(0);
  //gStyle->SetCanvasBorderMode(0);
  //gStyle->SetPalette(18,0);
  //gStyle->SetPaperSize(20,26);
  gStyle->SetOptStat("emrou");

  //gRandom->SetSeed(0);



  //--------------------------------------------------------------------------
  // Create a BendStudy object and call main analysis function
  //--------------------------------------------------------------------------

  BendStudy *analyzeBending = new BendStudy( jobOptionsFile );

  analyzeBending->AnalyzeFromNtuples();
  
  return 0;

  
} // end main() 


//=============================================================================
// Constructors
//=============================================================================

BendStudy::BendStudy( ) : LArIATAnalysis( ) { 

  tuple = NULL;



}


BendStudy::BendStudy( char* jobOptionsFile ) : LArIATAnalysis( jobOptionsFile ) { 

  tuple = NULL;

  

  //== Open list of ntuple files 
  if( UI->inputFilesSet ){
    openNtupleFiles( UI->inputFiles, tuple );
  }
  else{
    cout << endl << "#### ERROR : Input files not properly specified!!!!" << endl << endl;
    exit(0);
  }

  //== Open output root file and postscript file
  if( UI->rootOutputFileSet ){
    outputFile = new TFile( UI->rootOutputFile, "RECREATE" );
    //ps = new TPostScript( UI->psOutputFile, 112 );
    //ps->Range(26,18); 
    //psPage = 1; 
  }
  else{
    cout << endl << "#### No output files specified!!!!" << endl << endl;
  }

  //output list of selected primary particle IDs

  if (UI->SelEventListSet){
    
    IDfile.open(UI->SelEventList, ios::trunc);
  }

  //== number of events to process
  if( UI->numEventsToProcessSet )
    numEventsToProcess = UI->numEventsToProcess;
  else
    numEventsToProcess = 100000000;


  
  
  verbose = UI->verbose;
  isMC = UI->isMC;
  applyMassCut =  UI->applyMassCut;
  
}


//=============================================================================
// AnalyzeFromNtuples
//=============================================================================

void BendStudy::AnalyzeFromNtuples(){

  // counters for beam cuts



  double numEventsStart = 0;
  double numMassCut = 0;
  double numtofvalid = 0;


  double numTrackTotal = 0;
  double numMultiHit = 0 ;
  double numZcutoff = 0;
  double numHasCandidate = 0; 
  double numTwoTracks = 0;
  double numzProjShort = 0;
  double numBendZcut = 0;
  double numRightSided = 0;
  double numHasLong = 0 ;
  double numXYZMatching = 0;

  

  double numWCTrack = 0;
  double NumEventsSelList = 0;
  bool intHistFilled = false;


  const double wc_zpos_val = 0.1;
  const double *wc_zpos;
  wc_zpos = &wc_zpos_val;

  double bendZcut  =  UI->zBeamCutoff; //using this leftover input option to select stubby tracks





  // ### some variables that are needed for the xs calc ###
  double z2 = UI->zSlabSize; //<-- slab size. need to move this to jobOptions


  
  // globals move later...
  double mass = 938.57;
  //double z = 0.03;
  //double z2 = 0.5;

  double rho = 1.3954;                   // ## g/cm3
  double molar_mass = 39.95;                    // ## g/mol
  double N_A = 6.022 * pow(10, 23);      // ## num/mol
  
  double sparse_recip_num_density = molar_mass / (rho * z2 * N_A); // ## cm2/num
  double barn = pow(10, -24);



  std::cout << "I calculate the p-Ar cross section! \n";
  
  EventSelector *ES = new EventSelector();

  BeamSelector *BS = new BeamSelector();
  
  bookNtuple( tuple , UI->isMC);
  if (tuple == 0) return;

  Long64_t nentries = tuple->GetEntriesFast();
  TCanvas *c = new TCanvas("c","c",1500, 750);
  //c->Divide(1,2);

  //TCanvas *csplit = new TCanvas("csplit","csplit",1500, 750);
  //csplit->Divide(2,2);

// ## Histograms for entering tracks
  //TH2D *tpcInTracksXY =  new TH2D("tpcInTracksXY","Position of TPC track start XY",
    //200, -100, 100, 200, -100, 100);
  //TH2D *wctrkSelectedXY =  new TH2D("wctrkSelectedXY","Position of Wire Chamber Track",
  //200, -100, 100, 200, -100, 100);
  
  TH1D *ShortStartZ =  new TH1D("ShortStartZ","Earlies track start in Z",200, 0, 90);
  TH2D *ZLengthStart = new TH2D("ZLengthStart","Z projection vs Start Position", 200, 0, 90, 200, 0, 90);
  TH1D *ShortEndZ =  new TH1D("ShortEndZ","Earlies track end in Z",25, 0, 10);




  //TH1D *BadTrackStartZ =  new TH1D("BadTrackStartZ","non-selected track start in Z",10, 0, 4);
  //TH1D *InTrackLength =  new TH1D("InTrackLength","Entering Track Length",100, 0, 100);
  TH1D *ShortLength =  new TH1D("ShortLength","Earliest Entering Track Length",200, 0, 100);
  TH1D *LongLength =  new TH1D("LongLength","Second Entering Track Length",200, 0, 100);
  //TH1D *BadTrackLength =  new TH1D("BadTrackLength","non-selected Entering Track Length",250, 0, 100);

  TH1D *inTracksNumHist = new TH1D("inTracksNumHist","Number of Entering Tracks TPC",100,0,10);



  TH2D * BendingZHist =  new TH2D("BendingZHist"," Bending vs Z_{tpc}",50,0,10,400,-8,8);
  TH2D * BendingXHist =  new TH2D("BendingXHist","Bending vs X_{tpc}",50,0,48,400,-8,8);

  //TH1D * delXYHistPx =  new TH1D("delXYHistPx","tpc to wc delta x",200,-100,100,200,-100,100);
  //TH1D * delXYHistPy =  new TH1D("delXYHistPy","tpc to wc delta x",200,-100,100,200,-100,100);

  //TH2D *BadTrackHist =  new TH2D("BadTrackHist","non-selected tracks",200,-100,100,200,-100,100);
  //TH2D *delBadTrackHist =  new TH2D("delBadTrackHist","non-selected tracks",200,-100,100,200,-100,100);

  //TH1D * zProjTrack_out = new TH1D("zProjTrack_out","",100,0,40);
  //TH1D * zProjTrack_in = new TH1D("zProjTrack_in","",100,0,40);
 // TH1D * InTrackLength_out = new TH1D("InTrackLength_out","",100,0,40);
 // TH1D * InTrackLength_in = new TH1D("InTrackLength_in","",100,0,40);
   //z projections of tracks, in/out of circle cut
  TH1D * zProjTrack_short = new TH1D("zProjTrack_short","Earliest Entering Track Length - Z projection",100,0,90);
  TH1D * zProjTrack_all = new TH1D("zProjTrack_all","Track Length - Z projection",200,0,90);
  //TH1D * InTrackTPCnum_in = new TH1D("InTrackTPCnum_in","",100, 0, 10);
  //TH1D * InTrackTPCnum_out = new TH1D("InTrackTPCnum_out","",100, 0, 10);

  TH1D *BendingProxHist = new TH1D("BendingProxHist","short to long tagged track distance", 50, 0, 10);




  

  
  TH1D *BeamMassHist = new TH1D("BeamMassHist","Beamline particle Mass", 100, 0,3000);
  TH1D *BeamMassCutHist = new TH1D("BeamMassCutHist","Beamline particle Mass - after Cut", 100, 0,3000);
  
  //plots for EventSelection Branches

  //TH1D * BranchDistHist = new TH1D("BranchDistHist", "Inelastic Event Branch Distance",50, 0,25);
  //TH1D * ClusterDistHist = new TH1D("ClusterDistHist", "Additional Branch Distance (Type 4)",50, 0,25);
  // ## plot markers for cuts

  //TLine  *MassMinLine =  new TLine(0,UI->MassCutMin, 100, UI->MassCutMin);
  //TLine  *MassMaxLine =  new TLine(0,UI->MassCutMax, 100, UI->MassCutMax);

  // ## xs histos ##
  //TH1D *hreco_initialKE = new TH1D("hreco_initialKE", "initial ke", 20, 0, 1000);
  //TH1D *hreco_intke = new TH1D("hreco_intke", "int ke", 20, 0, 1000);
  //TH1D *hreco_incke = new TH1D("hreco_incke", "inc ke", 20, 0, 1000);
  //TH1D *hreco_xs = new TH1D("hreco_xs", "P-Ar Inelastic XS", 20, 0, 1000);

  // ## to read in from corrections file
  //TH1D *hreco_intke_background = new TH1D("hreco_intke_background", "int ke background", 20, 0, 1000);
  //TH1D *hreco_incke_background = new TH1D("hreco_incke_background", "inc ke background", 20, 0, 1000);
  //TH2D *hreco_unfolding_matrix_normalized = new TH2D("hreco_unfolding_matrix_normalized", "energy unfolding matrix", 20, 0, 1000, 20, 0, 1000);
  //TH1D *hreco_intke_eff = new TH1D("hreco_intke_eff", "interaction selection efficiency", 20, 0, 1000);
  //TH1D *hreco_incke_eff = new TH1D("hreco_incke_eff", "incident selection efficiency", 20, 0, 1000);

  // ## looping once to find Beamline center ##

  




  // ## event loop ##

  int max_numInEvent = 0;

  if(verbose){std::cout << "Starting main loop" << std::endl;}
  for (Long64_t jentry=0; jentry < numEventsToProcess && jentry < nentries; jentry++){
    
    Long64_t ientry = tuple->LoadTree(jentry); 
    if (ientry < 0){
      continue;}
    Long64_t nb = 0;
    nb = tuple->GetEntry(jentry);
    numEventsStart++;
    if(verbose){printEvent();}
    double ParticleMass;
    int numInEvent = 0;

    std::vector<int> selected_short = {-1,-1,-1,-1};
    std::vector<int> selected_downstream = {-1,-1,-1,-1}; 

    
      if (num_wctracks !=1){if(verbose){std::cout << "no WC track \n" << std::endl;}
        continue;}
      else{

          numWCTrack++;
          if(tofObject[0] < 0){
            if(verbose){std::cout << "no valid ToF \n" << std::endl;}
            continue;}
          numtofvalid++;
          ParticleMass = -9999999. ;

          bool isProton = BS->MassCut(wctrk_momentum[0], tofObject[0], ParticleMass, UI->MassCutMin, UI->MassCutMax);

          BeamMassHist->Fill(ParticleMass);
          if(!isProton){
            if(applyMassCut){continue;}
          }

          if (isProton){
            BeamMassCutHist->Fill(ParticleMass);
            numMassCut++;
          }

          //############## my crappy match finding function (I'll pack it somewhere later)


        //std::vector<std::vector<double>> zIndexes; //vector w/ vectors of indices for point z min, second min, z max.
        std::vector<std::vector<double>> matchCandidates; 
        std::vector<double> centeringMatch;

        double numZtracks = 0;
        double numCylinder = 0; 
        for (int itrack = 0; itrack < ntracks_reco ;  itrack++ ){
          numTrackTotal++;

          if ((*ntrack_hits)[itrack] < 2) {continue;}
          numMultiHit++;


          std::vector<int> zIndices = UtilityFunctions::zOrderedTrack(track_zpos,itrack,ntrack_hits);

          int index1 = zIndices[0];
          int index2 = zIndices[1];
          int indexLast = zIndices[2];
          double zmin1 = (*track_zpos)[itrack][zIndices[0]];
          double zmin2 = (*track_zpos)[itrack][zIndices[1]];
          double zmax = (*track_zpos)[itrack][zIndices[2]];


          double zProjLength = zmax - zmin1;
          zProjTrack_all->Fill(zProjLength);
          ShortLength->Fill((*track_length)[itrack]);

          ZLengthStart->Fill(zmin1, zProjLength);

          if(zProjLength > UI->zProjCut){continue;}
          numzProjShort++;

          ShortStartZ->Fill(zmin1);

          

          if (zmin1 > bendZcut){continue;}
          numBendZcut++;

          ShortEndZ->Fill(zmax);
          if(UI->BendDirFilter){
            if ((*track_xpos)[itrack][index1] < (*track_xpos)[itrack][indexLast]){continue;}
          }
          numRightSided++;

          

          double XYZMin = 99.;
          int DownstreamID = -1;
          std::vector<int>  downstreamIndices;

          bool found_downstream  = false;

          for (int otrack = 0 ;  otrack < ntracks_reco; otrack++){
            if (otrack == itrack){continue;}

            std::vector<int> otrackIndices = UtilityFunctions::zOrderedTrack(track_zpos,otrack,ntrack_hits);
            if ((*track_zpos)[otrack][otrackIndices[0]] < UI->zTPCCutoff){

              double delX = (*track_xpos)[itrack][indexLast] - (*track_xpos)[otrack][otrackIndices[0]];
              double delY = (*track_ypos)[itrack][indexLast] - (*track_ypos)[otrack][otrackIndices[0]];
              double delZ = (*track_zpos)[itrack][indexLast] - (*track_zpos)[otrack][otrackIndices[0]];

              double delMatch =  sqrt( pow(delX,2) + pow(delY,2) + pow(delZ,2));

              if (delMatch < XYZMin){
                XYZMin = delMatch;
                DownstreamID = otrack;
                downstreamIndices = otrackIndices;
                found_downstream = true;
              }
            }
          }

          if (!found_downstream){continue;}
          numHasLong++;

          BendingProxHist->Fill(XYZMin);
          if (XYZMin > UI->branchMaxDist){continue;}
          numXYZMatching++;
          


          selected_short[numInEvent] = itrack;
          selected_downstream[numInEvent] = DownstreamID;
          numInEvent++;

          LongLength->Fill((*track_length)[DownstreamID]);

          double zmaxFit = (*track_zpos)[DownstreamID][downstreamIndices[0]] + 4;

          TGraph *BendGraph =  new TGraph(); // X as a function of Z

          for (int ipoint = 0; ipoint < (*ntrack_hits)[DownstreamID] ; ipoint++){
            BendGraph->SetPoint(ipoint,(*track_zpos)[DownstreamID][ipoint],(*track_xpos)[DownstreamID][ipoint]);
          }
          
          TF1 *FitFz = new TF1("FitFz","[0]+[1]*x",0, zmaxFit);
          FitFz->SetParameters(22,1);
          c->cd(1);

          //TGraph *ShortGraph =  new TGraph();
          /*
          for (int ipoint = 0; ipoint < (*ntrack_hits)[Short_ID] ; ipoint++){
            ShortGraph->SetPoint(ipoint,(*track_zpos)[Short_ID][ipoint],(*track_xpos)[Short_ID][ipoint]);
          }
          ShortGraph->SetMarkerStyle(17);
          ShortGraph->SetMarkerColor(3);
          BendGraph->SetMarkerStyle(17);
          BendGraph->SetMarkerColor(4);*/
          

          //BendGraph->GetXaxis()->SetLimits(-1,UI->zTPCCutoff);
          //BendGraph->GetHistogram()->SetMinimum((wctrk_XFace[0] - xMeanTPCentry - UI->rCircleCut));
          //BendGraph->GetHistogram()->SetMaximum((wctrk_XFace[0] - xMeanTPCentry + UI->rCircleCut));
          //BendGraph->Draw("AP");
          //BendGraph->GetXaxis()->SetTitle("Z[cm]");
          //BendGraph->GetYaxis()->SetTitle("X[cm]");


          BendGraph->Fit("FitFz","R");

          /*ShortGraph->Draw("P SAME");

          TGraph *OthersGraph = new TGraph();
          
          if (matchCandidates.size() > 2){
            int graphPtBuffer = 0;
            for (int itrack = 0; itrack < ntracks_reco; itrack++ ){
              if (itrack != Short_ID && itrack != DownstreamID){
                //std::cout << "which track" << itrack << std::endl;
                
                for(int ipoint = 0; ipoint < (*ntrack_hits)[itrack]; ipoint++){
                  //std::cout << ipoint << "point filled" << std::endl;
                  
                  OthersGraph->SetPoint(graphPtBuffer,(*track_zpos)[itrack][ipoint],
                    (*track_ypos)[itrack][ipoint]);
                  graphPtBuffer++;
                }
              }
            }
          }
          OthersGraph->SetMarkerStyle(17);
          OthersGraph->SetMarkerColor(12);
          OthersGraph->Draw("P SAME");*/

          //c->Update();

          
          //char graph_title[100];
          //sprintf(graph_title,"plotting/images/Bending/Graph%d.png",event);
          //c->Update();
          //c->Print(graph_title,"png");




          for(int ipoint = 0; ipoint < (*ntrack_hits)[itrack]; ipoint++){
            double xval  = (*track_xpos)[itrack][ipoint];
            double zval = (*track_zpos)[itrack][ipoint];
            double projX = FitFz->Eval(zval);
            double deltaX = xval - projX;
            BendingZHist->Fill(zval,deltaX);
            BendingXHist->Fill(projX,deltaX);
          }

          delete BendGraph;
          //delete ShortGraph

          
        }//end of track loop
        if (numInEvent > max_numInEvent){
        max_numInEvent = numInEvent;}
        if(UI->SelEventListSet && numInEvent > 0){IDfile << jentry << "\t" << ParticleMass << "\t" << selected_short[0] << "\t"<< selected_downstream[0] << 
        "\t" << selected_short[1]<< "\t"<<  selected_downstream[1]  << "\t" << selected_short[2]<< "\t"<< selected_downstream[2]  <<
         "\t" << selected_short[3]<< "\t"<< selected_downstream[3] <<"\n";}
      }//events with wc track

      


    }//end of event loop


    

 

    std::cout << "\n------- Cut Results -------\n"<< std::endl;

    std::cout << "Total events processed: "<< numEventsStart << std::endl;
    std::cout << "Events with only one wc track: "<< numWCTrack << std::endl;
    std::cout << "Events with valid ToF value: "<< numtofvalid << std::endl;
    std::cout << "Events passing mass cut: "<< numMassCut << std::endl;

    std::cout << "\n------- Per Track Stats -------\n"<< std::endl;
    std::cout << "Tracks processed: " << numTrackTotal << std::endl;
    std::cout << "tracks with more than 1 hit" << numMultiHit << std::endl;
    std::cout << "Tracks with Short track Z proj  < " << UI->zProjCut << ": " << numzProjShort << std::endl;
    std::cout << "Tracks with Short track start < "<< UI->zBeamCutoff << ": " << numBendZcut << std::endl;
    std::cout << "Tracks with track tilt to the Right in XZ: "<< numRightSided  << std::endl;
    std::cout << "Tracks with Long track found :"<< numHasLong << std::endl;
    std::cout << "Tracks with Long track matched within r < " << UI->branchMaxDist  << ": "<< numXYZMatching << std::endl; 

    std::cout << "Max number of tracks selected in an event : " << max_numInEvent << std::endl; 





  






  if (!(isMC)){
    
    if(UI->rootOutputFileSet){
      if(verbose){std::cout << "Writing to outputFile" << std::endl;}
      outputFile->cd();
        
      ShortLength->Write();
      LongLength->Write();
      ShortStartZ->Write();
      ShortEndZ->Write();
      BeamMassHist->Write();
      zProjTrack_short->Write();
      zProjTrack_all->Write();

      BendingXHist->Write();
      BendingZHist->Write();
      BendingProxHist->Write();
      ZLengthStart->Write();


    }
  }
  if(UI->SelEventListSet){IDfile.close();}

}



