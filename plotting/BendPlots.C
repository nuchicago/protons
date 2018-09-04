void BendPlots(){

  gROOT->SetBatch(true);

  TFile* XsecFile = new TFile("../files/BendOutput.root","READ");
  gStyle->SetOptStat("emrou");

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
  c1->Print("images/Bending/delXYHist.png","png");

  
  TH1D *zProjTrack_short = (TH1D*)XsecFile->Get("zProjTrack_short");

  zProjTrack_short->SetTitle("Track Z Projection - Selected Tracks");
  zProjTrack_short->GetXaxis()->SetTitle("Z projection [cm]");
  zProjTrack_short->GetXaxis()->SetRangeUser(0,15);

  zProjTrack_short->Draw();

  c1->Print("images/Bending/zProjShortTracks.png","png");

  TH2D *BendingZHist = (TH2D*)XsecFile->Get("BendingZHist");

  BendingZHist->GetXaxis()->SetTitle("Z_{tpc} [cm]");
  //BendingZHist->GetXaxis()->SetRangeUser(0,10);
  //BendingZHist->GetYaxis()->SetRangeUser(0,10);
  BendingZHist->GetYaxis()->SetTitle("#Delta X [cm]");
  BendingZHist->Draw("COLZ");

  c1->Print("images/Bending/BendingZHist.png","png");

  TH2D *BendingXHist = (TH2D*)XsecFile->Get("BendingXHist");

  BendingXHist->GetXaxis()->SetTitle("X_{tpc} [cm]");
  BendingXHist->GetYaxis()->SetTitle("#Delta X [cm]");
  //BendingXHist->GetYaxis()->SetRangeUser(0,10);
  BendingXHist->Draw("COLZ");

  c1->Print("images/Bending/BendingXHist.png","png");
 
  TH1D *BendingProxHist = (TH1D*)XsecFile->Get("BendingProxHist");

  //BendingProxHist->SetTitle("Track Z Projection - Selected Tracks");
  BendingProxHist->GetXaxis()->SetTitle("Track Distance [cm]");
  BendingProxHist->Draw();

  c1->Print("images/Bending/BendingProxHist.png","png");
 
  TH1D *ShortLength = (TH1D*)XsecFile->Get("ShortLength");
  
  ShortLength->GetXaxis()->SetTitle("length [cm]");
  ShortLength->GetXaxis()->SetRangeUser(0,15);
  ShortLength->GetYaxis()->SetTitle("");
  ShortLength->Draw("COLZ");
  c1->Print("images/Bending/ShortLength.png","png");

  TH1D *ShortStartZ = (TH1D*)XsecFile->Get("ShortStartZ");
  
  ShortStartZ->GetXaxis()->SetTitle("Z [cm]");
  //ShortStartZ->GetXaxis()->SetRangeUser(0,4);
  ShortStartZ->GetYaxis()->SetTitle("");
  ShortStartZ->Draw("COLZ");
  c1->Print("images/Bending/ShortStartZ.png","png");
    
  TH1D *ShortEndZ = (TH1D*)XsecFile->Get("ShortEndZ");
  
  ShortEndZ->GetXaxis()->SetTitle("Z [cm]");
  //ShortEndZ->GetXaxis()->SetRangeUser(0,4);
  ShortEndZ->GetYaxis()->SetTitle("");
  ShortEndZ->Draw("COLZ");
  c1->Print("images/Bending/ShortEndZ.png","png");

  /*THStack * InTrackTPCnum_st = new THStack("InTrackTPCnum_st","");

  TH1D *InTrackTPCnum_in = (TH1D*)XsecFile->Get("InTrackTPCnum_in");
  TH1D *InTrackTPCnum_out = (TH1D*)XsecFile->Get("InTrackTPCnum_out");
  InTrackTPCnum_in->SetFillColor(30);
  InTrackTPCnum_out->SetFillColor(46);

  TLegend *InTrackNumLeg = new TLegend(0.7,0.85,0.9,0.9);
  InTrackNumLeg->AddEntry(InTrackTPCnum_in, "Tracks inside circle cut","f");\
  InTrackNumLeg->AddEntry(InTrackTPCnum_out, "Tracks outside circle cut","f");

  InTrackTPCnum_in->GetXaxis()->SetTitle("Track Length[cm]");
  InTrackTPCnum_in->SetTitle("Entering Track Length ");

  InTrackTPCnum_st->Add(InTrackTPCnum_in);
  InTrackTPCnum_st->Add(InTrackTPCnum_out);
  InTrackTPCnum_st->Draw();
  InTrackNumLeg->Draw("same");
  c1->Print("images/Data/InTrackTPCnum_stack.png","png");*/

 



}