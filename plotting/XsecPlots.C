void XsecPlots(){





  TFile* XsecFile = new TFile("../files/XsecOutput.root","READ");


  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  TH2D *delXYHist = (TH2D*)XsecFile->Get("delXYHist");
  //TEllipse *XYCut =  new TEllipse(delXYHist->GetMean(1),delXYHist->GetMean(2),5,0);
  TEllipse *XYCut =  new TEllipse(0,0,5,0);
  XYCut->SetFillColorAlpha(0,0.0);
  XYCut->SetLineColor(2);
  XYCut->SetLineWidth(5);

  c1->cd(1);
  
  delXYHist->GetXaxis()->SetTitle("X [cm]");
  delXYHist->GetYaxis()->SetTitle("Y [cm]");
  delXYHist->SetTitle("Wire Chamber - TPC position difference");
  delXYHist->Draw("COLZ");
  XYCut->Draw("same");


  c1->Print("delXYHist.png","png");

  TH1D *delThetaHist = (TH1D*)XsecFile->Get("delThetaHist");
  
  delThetaHist->GetXaxis()->SetTitle("#theta [rad]");
  delThetaHist->GetYaxis()->SetTitle("Number of Tracks");
  delThetaHist->SetTitle("Wire Chamber - TPC #theta difference");
  delThetaHist->Draw("COLZ");

  c1->Print("delThetaHist.png","png");


  TH1D *delPhiHist = (TH1D*)XsecFile->Get("delPhiHist");
  
  delPhiHist->GetXaxis()->SetTitle("#phi [rad]");
  delPhiHist->GetYaxis()->SetTitle("Number of Tracks");
  delPhiHist->SetTitle("Wire Chamber - TPC #phi difference");
  delPhiHist->Draw("COLZ");

  c1->Print("delPhiHist.png","png");





}