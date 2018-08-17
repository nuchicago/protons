void XsecPlotsMC(){

  gROOT->SetBatch(true);

  TFile* XsecFile = new TFile("../files/XSecOutputMC.root","READ");
  gStyle->SetOptStat("emrou");

  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);

  c1->cd(1);

  TH1D *hreco_initialKE = (TH1D*)XsecFile->Get("hreco_initialKE");
  
  hreco_initialKE->GetXaxis()->SetTitle("KE [MeV]");
  hreco_initialKE->GetXaxis()->SetRangeUser(0,1000);

  hreco_initialKE->SetTitle("Initial KE");
  hreco_initialKE->Draw("COLZ");

  c1->Print("images/MC/hreco_initialKE.png","png");

  TH1D *hreco_incke = (TH1D*)XsecFile->Get("hreco_incke");
  
  hreco_incke->GetXaxis()->SetTitle("KE [MeV]");
  hreco_incke->GetXaxis()->SetRangeUser(0,1000);

  hreco_incke->SetTitle("Incident KE");
  hreco_incke->Draw("COLZ");

  c1->Print("images/MC/hreco_incke.png","png");

  TH1D *hreco_intke = (TH1D*)XsecFile->Get("hreco_intke");
  
  hreco_intke->GetXaxis()->SetTitle("KE [MeV]");
  hreco_intke->GetXaxis()->SetRangeUser(0,1000);
 
  hreco_intke->SetTitle("Interacting KE");
  hreco_intke->Draw("COLZ");

  c1->Print("images/MC/hreco_intke.png","png");

  /*
  TH1D *tpcInTracksZ = (TH1D*)XsecFile->Get("tpcInTracksZ");
  
  tpcInTracksZ->GetXaxis()->SetTitle("Z [cm]");
  tpcInTracksZ->GetYaxis()->SetTitle("");
  tpcInTracksZ->Draw("COLZ");
  c1->Print("images/MC/tpcInTracksZ.png","png");


  TH1D *InTrackLength = (TH1D*)XsecFile->Get("InTrackLength");
  
  InTrackLength->GetXaxis()->SetTitle("length [cm]");
  InTrackLength->GetYaxis()->SetTitle("");
  InTrackLength->Draw("COLZ");
  c1->Print("images/MC/InTrackLength.png","png");




  TH1D *PrimaryStartZ = (TH1D*)XsecFile->Get("PrimaryStartZ");
  
  PrimaryStartZ->GetXaxis()->SetTitle("length [cm]");
  PrimaryStartZ->GetYaxis()->SetTitle("");
  PrimaryStartZ->Draw("COLZ");
  c1->Print("images/MC/PrimaryStartZ.png","png");

  TH1D *PrimaryLength = (TH1D*)XsecFile->Get("PrimaryLength");
  
  PrimaryLength->GetXaxis()->SetTitle("length [cm]");
  PrimaryLength->GetYaxis()->SetTitle("");
  PrimaryLength->Draw("COLZ");
  c1->Print("images/MC/PrimaryLength.png","png");


  TH1D *primary_dedx = (TH1D*)XsecFile->Get("primary_dedx");

  primary_dedx->GetXaxis()->SetTitle("dE/dx");
  primary_dedx->Draw("");
  c1->Print("images/MC/primary_dedx.png","png");*/

  TH1D *BranchDistHist = (TH1D*)XsecFile->Get("BranchDistHist");
  
  BranchDistHist->GetXaxis()->SetTitle("r [cm]");
  BranchDistHist->GetYaxis()->SetTitle("");
  BranchDistHist->Draw("");
  c1->Print("images/MC/BranchDistHist.png","png");

  TH1D *ClusterDistHist = (TH1D*)XsecFile->Get("ClusterDistHist");
  
  ClusterDistHist->GetXaxis()->SetTitle("r [cm]");
  ClusterDistHist->GetYaxis()->SetTitle("");
  ClusterDistHist->Draw("");
  c1->Print("images/MC/ClusterDistHist.png","png");
}