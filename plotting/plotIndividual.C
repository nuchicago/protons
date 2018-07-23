void plotIndividual(){

  gROOT->SetBatch(true);

  TFile* XsecFile = new TFile("../files/XsecOutput.root","READ");
  gStyle->SetOptStat("emr");

  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);

  c1->cd(1);


  TGraph *gtotal_res_dedx = (TGraph*)XsecFile->Get("gtotal_res_dedx");
          gtotal_res_dedx->GetXaxis()->SetTitle("Residual Range [cm]");
          gtotal_res_dedx->GetYaxis()->SetTitle("dE/dx [MeV/cm]");
          gtotal_res_dedx->SetTitle("");
          gtotal_res_dedx->SetMarkerStyle(6);
          //gtotal_res_dedx->SetMarkerSize(1);
          gtotal_res_dedx->SetMarkerColor(1);
          

  
  gtotal_res_dedx->GetHistogram()->SetMaximum(60.); 
  gtotal_res_dedx->GetHistogram()->SetMinimum(0.);  
  gtotal_res_dedx->Draw("AP");
  c1->Update();
  
  c1->Print("images/gtotal_res_dedx.png","png");
  

}