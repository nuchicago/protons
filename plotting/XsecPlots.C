void XsecPlots(){

  TFile* XsecFile = new TFile("../files/XsecOutput.root","READ");
  gStyle->SetOptStat(1111);

  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  TH2D *delXYHist = (TH2D*)XsecFile->Get("delXYHist");
  TEllipse *XYCut =  new TEllipse(delXYHist->GetMean(1),delXYHist->GetMean(2),5,0);
  XYCut->SetFillColorAlpha(0,0.0);
  XYCut->SetLineColor(2);
  XYCut->SetLineWidth(5);

  c1->cd(1);
  
  delXYHist->GetXaxis()->SetTitle("X_{wc} - X_{tpc} [cm]");
  delXYHist->GetYaxis()->SetTitle("Y_{wc} - Y_{tpc} [cm]");
  delXYHist->GetXaxis()->SetRangeUser(-45,45);
  delXYHist->GetYaxis()->SetRangeUser(-45,45);
  delXYHist->SetTitle("Wire Chamber - TPC position difference");
  delXYHist->Draw("COLZ");
  XYCut->Draw("same");

  c1->Print("images/delXYHist.png","png");

  TH1D *delThetaHist = (TH1D*)XsecFile->Get("delThetaHist");
  
  delThetaHist->GetXaxis()->SetTitle("#theta_{wc} = #theta_{tpc} [rad]");
  delThetaHist->GetYaxis()->SetTitle("Number of Events");
  delThetaHist->SetTitle("Wire Chamber - TPC #theta difference");
  delThetaHist->Draw("COLZ");

  c1->Print("images/delThetaHist.png","png");

  TH1D *inTracksNumHist = (TH1D*)XsecFile->Get("inTracksNumHist");
  
  inTracksNumHist->GetXaxis()->SetTitle("Number of Tracks");
  inTracksNumHist->GetYaxis()->SetTitle("Number of Events");
  inTracksNumHist->GetXaxis()->SetRangeUser(0,8);
  inTracksNumHist->SetTitle("Number of Tracks starting at Z < 4cm");
  inTracksNumHist->Draw("COLZ");

  c1->Print("images/inTracksNumHist.png","png");

    TH1D *wctrkNumHist = (TH1D*)XsecFile->Get("wctrkNumHist");
  
  wctrkNumHist->GetXaxis()->SetTitle("Number of Wire Chamber Tracks");
  wctrkNumHist->GetYaxis()->SetTitle("Number of Events");
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

  TH2D *tpcInTracksXY = (TH2D*)XsecFile->Get("tpcInTracksXY");
  
  tpcInTracksXY->GetXaxis()->SetTitle("X [cm]");
  tpcInTracksXY->GetYaxis()->SetTitle("Y[cm]");
  tpcInTracksXY->Draw("COLZ");
  tpcInTracksXY->GetXaxis()->SetRangeUser(0,47.5);
  tpcInTracksXY->GetYaxis()->SetRangeUser(-20,20);

  c1->Print("images/tpcInTracksXY.png","png");

  TH2D *wctrkPositionXY = (TH2D*)XsecFile->Get("wctrkPositionXY");
  
  wctrkPositionXY->GetXaxis()->SetTitle("X [cm]");
  wctrkPositionXY->GetYaxis()->SetTitle("Y[cm]");
  wctrkPositionXY->Draw("COLZ");
  wctrkPositionXY->GetXaxis()->SetRangeUser(0,47.5);
  wctrkPositionXY->GetYaxis()->SetRangeUser(-20,20);

  c1->Print("images/wctrkPositionXY.png","png");

  TH2D *wctrkSelectedXY = (TH2D*)XsecFile->Get("wctrkSelectedXY");
  
  wctrkSelectedXY->GetXaxis()->SetTitle("X [cm]");
  wctrkSelectedXY->GetYaxis()->SetTitle("Y[cm]");
  wctrkSelectedXY->Draw("COLZ");
  wctrkSelectedXY->GetXaxis()->SetRangeUser(0,47.5);
  wctrkSelectedXY->GetYaxis()->SetRangeUser(-20,20);

  c1->Print("images/wctrkSelectedXY.png","png");


  TH2D *delXYHist_pfX = (TH2D*)delXYHist->ProfileX("delXYHist_pfX");
  delXYHist_pfX->SetTitle("X projection of wc-tpc difference");
  delXYHist->GetXaxis()->SetTitle("X_{wc} - X_{tpc} [cm]");
  delXYHist_pfX->Draw("e");

  c1->Print("images/delXYHist_pfX.png","png");
  TH2D *delXYHist_pfY = (TH2D*)delXYHist->ProfileY("delXYHist_pfY");
  delXYHist_pfY->SetTitle("Y Projection of wc-tpc Difference");
  delXYHist_pfY->Draw("e");
  c1->Print("images/delXYHist_pfY.png","png");


}