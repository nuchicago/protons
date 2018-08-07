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
  
  TF1 *res_fx = new TF1("res_fx","[0]* x**([1])", 0, 90);
  res_fx->SetParameters(17, -0.42);
  res_fx->SetLineColor(2);

  TLegend *res_legend = new TLegend(0.5,0.7,0.9,0.90);
  res_legend->SetNColumns(2);
  res_legend->AddEntry(res_fx,"#left( #frac{dE}{dx} #right)_{hyp} = A R^{b}","l");
  res_legend->AddEntry((TObject*)0, "A = 17", "");
  res_legend->AddEntry(gtotal_res_dedx,"#left( #frac{dE}{dx} #right)_{reco}","p");
  res_legend->AddEntry((TObject*)0, "b = -0.42", "");

  res_fx->Draw("lsame");
  res_legend->Draw();


  c1->Update();
  
  c1->Print("images/gtotal_res_dedx.png","png");
  

}