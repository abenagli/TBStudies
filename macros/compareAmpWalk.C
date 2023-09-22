void compareAmpWalk()
{
  TFile* inFile = TFile::Open("plots/2.2/drawCTRMatrix_2.2.root","READ");

  TCanvas* c = new TCanvas();
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.2,-0.4,0.9,0.4) );
  hPad -> SetTitle(";amplitude [V];t-t_{0} [ns]");
  hPad -> Draw();
  
  
  TProfile* p_all_50um = new TProfile("p_all_50um","",200,0.,1.);
  int jj = 0;
  for(int ii = 1; ii <= 16; ++ii)
  {
    if( ii == 1 || ii == 2 || ii == 3 || ii == 4 || ii == 15 || ii == 16 ) continue;
    
    TProfile* p = (TProfile*)( inFile->Get(Form("p_time_vs_amp_CH%02d",ii)) );
    
    float refVal = 0.;
    for(int bin = 1; bin <= p->GetNbinsX(); ++bin)
    {
      float binCenter  = p -> GetBinCenter(bin);
      float binContent = p -> GetBinContent(bin);
      if( binCenter >= 0.5 )
      {
        refVal = binContent;
        break;
      }
    }

    TGraph* g = new TGraph();
    for(int bin = 1; bin <= p->GetNbinsX(); ++bin)
    {
      float binCenter  = p -> GetBinCenter(bin);
      float binContent = p -> GetBinContent(bin);
      int point = g->GetN();
      if( fabs(binContent-refVal) < 2. )
      {
        g -> SetPoint(point,binCenter,binContent-refVal);
        p_all_50um -> Fill(binCenter,binContent-refVal);
      }
    }
    
    g -> SetMarkerSize(0.7);
    g -> SetMarkerColor(51+jj*5);
    g -> SetLineColor(51+jj*5);
    g -> Draw("PL,same");

    ++jj;
  }

  p_all_50um -> SetMarkerSize(0.5);
  p_all_50um -> SetMarkerColor(kBlack);
  p_all_50um -> Draw("same");
  
  TF1* fitFunc_corrAmp_50um = new TF1(Form("fitFunc_corrAmp_50um"),"[0]*log([1]*x)+[2]",0.,1000.);
  fitFunc_corrAmp_50um -> SetParameters(-0.2,0.00000001,-0.5);
  p_all_50um -> Fit(fitFunc_corrAmp_50um,"QNRS","",0.5,1.);
  fitFunc_corrAmp_50um -> Draw("same");
  fitFunc_corrAmp_50um -> SetLineStyle(7);
  fitFunc_corrAmp_50um -> SetLineWidth(2);
  fitFunc_corrAmp_50um -> SetLineColor(kBlack);


  
  // TProfile* p_all_25um = new TProfile("p_all_25um","",200,0.,1.);
  // jj = 0;
  // for(int ii = 1; ii <= 16; ++ii)
  // {
  //   if( ii != 1 && ii != 2 && ii != 15 && ii != 16 ) continue;
    
  //   TProfile* p = (TProfile*)( inFile->Get(Form("p_time_vs_amp_CH%02d",ii)) );
    
  //   float refVal = 0.;
  //   for(int bin = 1; bin <= p->GetNbinsX(); ++bin)
  //   {
  //     float binCenter  = p -> GetBinCenter(bin);
  //     float binContent = p -> GetBinContent(bin);
  //     if( binCenter >= 0.2 )
  //     {
  //       refVal = binContent;
  //       break;
  //     }
  //   }

  //   TGraph* g = new TGraph();
  //   for(int bin = 1; bin <= p->GetNbinsX(); ++bin)
  //   {
  //     float binCenter  = p -> GetBinCenter(bin);
  //     float binContent = p -> GetBinContent(bin);
  //     int point = g->GetN();
  //     if( fabs(binContent-refVal) < 2. )
  //     {
  //       g -> SetPoint(point,binCenter,binContent-refVal);
  //       p_all_25um -> Fill(binCenter,binContent-refVal);
  //     }
  //   }
    
  //   g -> SetMarkerSize(0.7);
  //   g -> SetMarkerColor(41+jj*3);
  //   g -> SetLineColor(41+jj*3);
  //   g -> Draw("PL,same");

  //   ++jj;
  // }

  // p_all_25um -> SetMarkerSize(0.5);
  // p_all_25um -> SetMarkerColor(kBlack);
  // p_all_25um -> Draw("same");
}
