void XsecPlots(){

  gROOT->SetBatch(true);

  TFile* XsecFile = new TFile("../files/XsecOutput.root","READ");
  gStyle->SetOptStat("emr");

  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  TH2D *delXYHist = (TH2D*)XsecFile->Get("delXYHist");
  TEllipse *XYCut =  new TEllipse(delXYHist->GetMean(1),delXYHist->GetMean(2),5,0);
  XYCut->SetFillColorAlpha(0,0.0);
  XYCut->SetLineColor(2);
  XYCut->SetLineWidth(5);

  c1->cd(1);
  
  delXYHist->GetXaxis()->SetTitle("X_{wc} - X_{tpc} [cm]");
  delXYHist->GetYaxis()->SetTitle("Y_{wc} - Y_{tpc} [cm]");
  delXYHist->GetXaxis()->SetRangeUser(-48,48);
  delXYHist->GetYaxis()->SetRangeUser(-48,48);
  delXYHist->SetTitle("Wire Chamber - TPC position difference");
  delXYHist->Draw("COLZ");
  XYCut->Draw("same");

  c1->Print("images/delXYHist.png","png");

  TH1D *delThetaHist = (TH1D*)XsecFile->Get("delThetaHist");
  
  delThetaHist->GetXaxis()->SetTitle("#theta_{wc} = #theta_{tpc} [rad]");
  delThetaHist->SetTitle("Wire Chamber - TPC #theta difference");
  delThetaHist->Draw("COLZ");

  c1->Print("images/delThetaHist.png","png");

  TH1D *inTracksNumHist = (TH1D*)XsecFile->Get("inTracksNumHist");
  
  inTracksNumHist->GetXaxis()->SetTitle("Number of Tracks");
  
  inTracksNumHist->GetXaxis()->SetRangeUser(0,8);
  inTracksNumHist->SetTitle("Number of Tracks starting at Z < 4cm");
  inTracksNumHist->Draw("COLZ");

  c1->Print("images/inTracksNumHist.png","png");

  TH1D *wctrkNumHist = (TH1D*)XsecFile->Get("wctrkNumHist");
  
  wctrkNumHist->GetXaxis()->SetTitle("Number of Wire Chamber Tracks");
  wctrkNumHist->GetXaxis()->SetRangeUser(0,8);
  wctrkNumHist->SetTitle("Number of Wire Chamber Tracks per Event");
  wctrkNumHist->Draw("COLZ");

  c1->Print("images/wctrkNumHist.png","png");

  TH1D *delPhiHist = (TH1D*)XsecFile->Get("delPhiHist");
  
  delPhiHist->GetXaxis()->SetTitle("#phi_{wc} - #phi_{tpc} [rad]");
  delPhiHist->GetYaxis()->SetTitle("Number of Tracks");
  delPhiHist->SetTitle("Wire Chamber - TPC #phi difference");
  delPhiHist->Draw("COLZ");

  c1->Print("images/delPhiHist.png","png");

  TH1D *BeamMomentum = (TH1D*)XsecFile->Get("BeamMomentum");
  
  BeamMomentum->GetXaxis()->SetTitle("Momentum [MeV/c]");
  BeamMomentum->GetYaxis()->SetTitle("Number of Tracks");
  BeamMomentum->SetTitle("Beam Momentum");
  BeamMomentum->Draw("COLZ");

  c1->Print("images/BeamMomentum.png","png");

  TH1D *BeamToF = (TH1D*)XsecFile->Get("BeamToF");
  
  BeamToF->GetXaxis()->SetTitle("ToF [ns]");
  //BeamToF->GetXaxis()->SetRangeUser(0,200);
  BeamToF->GetYaxis()->SetTitle("Number of Tracks");
  BeamToF->SetTitle("Time of Flight");
  BeamToF->Draw("COLZ");

  c1->Print("images/BeamToF.png","png");

  TH1D *hreco_initialKE = (TH1D*)XsecFile->Get("hreco_initialKE");
  
  hreco_initialKE->GetXaxis()->SetTitle("KE [MeV]");
  hreco_initialKE->GetXaxis()->SetRangeUser(0,1000);

  hreco_initialKE->SetTitle("Initial KE");
  hreco_initialKE->Draw("COLZ");

  c1->Print("images/hreco_initialKE.png","png");

  TH1D *hreco_incke = (TH1D*)XsecFile->Get("hreco_incke");
  
  hreco_incke->GetXaxis()->SetTitle("KE [MeV]");
  hreco_incke->GetXaxis()->SetRangeUser(0,1000);

  hreco_incke->SetTitle("Incident KE");
  hreco_incke->Draw("COLZ");

  c1->Print("images/hreco_incke.png","png");

  TH1D *hreco_intke = (TH1D*)XsecFile->Get("hreco_intke");
  
  hreco_intke->GetXaxis()->SetTitle("KE [MeV]");
  hreco_intke->GetXaxis()->SetRangeUser(0,1000);
 
  hreco_intke->SetTitle("Interacting KE");
  hreco_intke->Draw("COLZ");

  c1->Print("images/hreco_intke.png","png");



  TH2D *tpcInTracksXY = (TH2D*)XsecFile->Get("tpcInTracksXY");
  
  tpcInTracksXY->GetXaxis()->SetTitle("X [cm]");
  tpcInTracksXY->GetYaxis()->SetTitle("Y[cm]");
  tpcInTracksXY->Draw("COLZ");
  tpcInTracksXY->GetXaxis()->SetRangeUser(-0,47.5);
  tpcInTracksXY->GetYaxis()->SetRangeUser(-20,20);

  c1->Print("images/tpcInTracksXY.png","png");


  TH1D *tpcInTracksZ = (TH1D*)XsecFile->Get("tpcInTracksZ");
  
  tpcInTracksZ->GetXaxis()->SetTitle("Z [cm]");
  tpcInTracksZ->GetYaxis()->SetTitle("");
  tpcInTracksZ->Draw("COLZ");
  c1->Print("images/tpcInTracksZ.png","png");


  TH1D *InTrackLength = (TH1D*)XsecFile->Get("InTrackLength");
  
  InTrackLength->GetXaxis()->SetTitle("length [cm]");
  InTrackLength->GetYaxis()->SetTitle("");
  InTrackLength->Draw("COLZ");
  c1->Print("images/InTrackLength.png","png");


  TH2D *wctrkPositionXY = (TH2D*)XsecFile->Get("wctrkPositionXY");
  
  wctrkPositionXY->GetXaxis()->SetTitle("X [cm]");
  wctrkPositionXY->GetYaxis()->SetTitle("Y[cm]");
  wctrkPositionXY->Draw("COLZ");
  wctrkPositionXY->GetXaxis()->SetRangeUser(0,47.5);
  wctrkPositionXY->GetYaxis()->SetRangeUser(-20,20);

  c1->Print("images/wctrkPositionXY.png","png");

  TH2D *wctrkSelectedXY = (TH2D*)XsecFile->Get("wctrkSelectedXY");
  
  wctrkSelectedXY->GetXaxis()->SetTitle("X [cm]");
  wctrkSelectedXY->GetYaxis()->SetTitle("Y [cm]");
  wctrkSelectedXY->Draw("COLZ");
  wctrkSelectedXY->GetXaxis()->SetRangeUser(0,47.5);
  wctrkSelectedXY->GetYaxis()->SetRangeUser(-20,20);

  c1->Print("images/wctrkSelectedXY.png","png");


  TH2D *BadTrackHist = (TH2D*)XsecFile->Get("BadTrackHist");
  
  BadTrackHist->GetXaxis()->SetTitle("X_{tpc} [cm]");
  BadTrackHist->GetYaxis()->SetTitle("Y_{tpc} [cm]");
  BadTrackHist->Draw("COLZ");
  BadTrackHist->GetXaxis()->SetRangeUser(0,48);
  BadTrackHist->GetYaxis()->SetRangeUser(-20,20);

  c1->Print("images/BadTrackHist.png","png");


  TH2D *delBadTrackHist = (TH2D*)XsecFile->Get("delBadTrackHist");
  
  delBadTrackHist->GetXaxis()->SetTitle("X_{wc} - X_{tpc} [cm]");
  delBadTrackHist->GetYaxis()->SetTitle("Y_{wc} - Y_{tpc} [cm]");
  delBadTrackHist->Draw("COLZ");
  delBadTrackHist->GetXaxis()->SetRangeUser(-48,48);
  delBadTrackHist->GetYaxis()->SetRangeUser(-48,48);

  c1->Print("images/delBadTrackHist.png","png");

  TH2D *PileupHist = (TH2D*)XsecFile->Get("PileupHist");
  
  PileupHist->GetXaxis()->SetTitle("X_{tpc} [cm]");
  PileupHist->GetYaxis()->SetTitle("Y_{tpc} [cm]");
  PileupHist->Draw("COLZ");
  PileupHist->GetXaxis()->SetRangeUser(0,47.5);
  PileupHist->GetYaxis()->SetRangeUser(-20,20);

  c1->Print("images/PileupHist.png","png");

  TH2D *delPileupHist = (TH2D*)XsecFile->Get("delPileupHist");
  
  delPileupHist->GetXaxis()->SetTitle("X_{wc} - X_{tpc} [cm]");
  delPileupHist->GetYaxis()->SetTitle("Y_{wc} - Y_{tpc} [cm]");
  delPileupHist->Draw("COLZ");
  delPileupHist->GetXaxis()->SetRangeUser(-48,48);
  delPileupHist->GetYaxis()->SetRangeUser(-48,48);

  c1->Print("images/delPileupHist.png","png");


  TH2D *delXYHist_pfX = (TH2D*)delXYHist->ProfileX("delXYHist_pfX");
  delXYHist_pfX->SetTitle("X projection of wc-tpc difference");
  delXYHist_pfX->GetXaxis()->SetTitle("X_{wc} - X_{tpc} [cm]");
  delXYHist_pfX->Draw("e");

  c1->Print("images/delXYHist_pfX.png","png");
  TH2D *delXYHist_pfY = (TH2D*)delXYHist->ProfileY("delXYHist_pfY");
  delXYHist_pfY->SetTitle("Y Projection of wc-tpc Difference");
  delXYHist_pfY->Draw("e");
  c1->Print("images/delXYHist_pfY.png","png");

  TH1D *BadTrackLength = (TH1D*)XsecFile->Get("BadTrackLength");
  
  BadTrackLength->GetXaxis()->SetTitle("length [cm]");
  BadTrackLength->GetYaxis()->SetTitle("");
  BadTrackLength->Draw("COLZ");
  c1->Print("images/BadTrackLength.png","png");

  TH1D *BadTrackStartZ = (TH1D*)XsecFile->Get("BadTrackStartZ");
  
  BadTrackStartZ->GetXaxis()->SetTitle("length [cm]");
  BadTrackStartZ->SetTitle("non-selected track start");
  BadTrackStartZ->Draw("COLZ");
  c1->Print("images/BadTrackStartZ.png","png");

  TH1D *PrimaryStartZ = (TH1D*)XsecFile->Get("PrimaryStartZ");
  
  PrimaryStartZ->GetXaxis()->SetTitle("length [cm]");
  PrimaryStartZ->GetYaxis()->SetTitle("");
  PrimaryStartZ->Draw("COLZ");
  c1->Print("images/PrimaryStartZ.png","png");

  TH1D *PrimaryLength = (TH1D*)XsecFile->Get("PrimaryLength");
  
  PrimaryLength->GetXaxis()->SetTitle("length [cm]");
  PrimaryLength->GetYaxis()->SetTitle("");
  PrimaryLength->Draw("COLZ");
  c1->Print("images/PrimaryLength.png","png");

  TH1D *BeamMassHist = (TH1D*)XsecFile->Get("BeamMassHist");
  
  BeamMassHist->GetXaxis()->SetTitle("Mass [MeV/c^{2}]");
  BeamMassHist->GetYaxis()->SetTitle("");
  BeamMassHist->Draw("COLZ");
  BeamMassHist->GetXaxis()->SetRangeUser(-0,3000);

  TLine *MassMaxLine = new TLine(0,1200,100,1200);
  MassMaxLine->SetLineColor(1);
  MassMaxLine->Draw("lsame");

  TLine *MassMinLine = new TLine(0,700, 100, 700);
  MassMinLine->SetLineColor(1);
  MassMinLine->Draw("lsame");

  c1->Print("images/BeamMassHist.png","png");

  TH1D *BeamMassCutHist = (TH1D*)XsecFile->Get("BeamMassCutHist");
  
  BeamMassCutHist->GetXaxis()->SetTitle("Mass [MeV/c^{2}]");
  BeamMassCutHist->GetYaxis()->SetTitle("");
  BeamMassCutHist->Draw("COLZ");
  BeamMassCutHist->GetXaxis()->SetRangeUser(-0,3000);
  c1->Print("images/BeamMassCutHist.png","png");

  TH2D *tofMomentHist = (TH2D*)XsecFile->Get("tofMomentHist");
  
  tofMomentHist->GetXaxis()->SetTitle("Momentum [MeV/c]");
  tofMomentHist->GetYaxis()->SetTitle("Time of Flight");
  tofMomentHist->Draw("COLZ");
  c1->Print("images/tofMomentHist.png","png");

  TH1D *numTracksSelHist = (TH1D*)XsecFile->Get("numTracksSelHist");

  numTracksSelHist->GetXaxis()->SetTitle("Number of Tracks");
  numTracksSelHist->Draw("");
  c1->Print("images/numTracksSelHist.png","png");

  TH1D *tpcPhiHist = (TH1D*)XsecFile->Get("tpcPhiHist");

  tpcPhiHist->GetXaxis()->SetTitle("#phi_{tpc}");
  tpcPhiHist->Draw("");
  c1->Print("images/tpcPhiHist.png","png");

  TH1D *wcPhiHist = (TH1D*)XsecFile->Get("wcPhiHist");

  wcPhiHist->GetXaxis()->SetTitle("#phi_{wc}");
  wcPhiHist->Draw("");
  c1->Print("images/wcPhiHist.png","png");

  TH1D *tpcThetaHist = (TH1D*)XsecFile->Get("tpcThetaHist");

  tpcThetaHist->GetXaxis()->SetTitle("#theta_{tpc}");
  tpcThetaHist->Draw("");
  c1->Print("images/tpcThetaHist.png","png");

  TH1D *wcThetaHist = (TH1D*)XsecFile->Get("wcThetaHist");

  wcThetaHist->GetXaxis()->SetTitle("#theta_{wc}");
  wcThetaHist->Draw("");
  c1->Print("images/wcThetaHist.png","png");

  TH1D *primary_dedx = (TH1D*)XsecFile->Get("primary_dedx");

  primary_dedx->GetXaxis()->SetTitle("dE/dx");
  primary_dedx->Draw("");
  c1->Print("images/primary_dedx.png","png");
}