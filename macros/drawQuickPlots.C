void drawAvgAmpRMS(const std::string& inFileName)
{
  TFile* inFile = TFile::Open(inFileName.c_str(),"READ");

  TProfile2D* p2 = (TProfile2D*)( inFile->Get("p2_amp") );
  p2 -> Draw("COLZ");

  TH1F* h1_avgAmp_50um = new TH1F("h1_avgAmp_50um","",50,0.,1.);
  TH1F* h1_avgAmp_25um = new TH1F("h1_avgAmp_25um","",50,0.,1.);
  
  for(int binx = 1; binx <= p2->GetNbinsX(); ++binx)
    for(int biny = 1; biny <= p2->GetNbinsX(); ++biny)
    {
      float content = p2->GetBinContent(binx,biny);

      if( biny == 4 || ( (biny == 1) && (binx == 3 || binx == 4) ) )
        h1_avgAmp_25um -> Fill( content );
      else
        h1_avgAmp_50um -> Fill( content );        
    }

  new TCanvas();
  
  h1_avgAmp_50um -> SetLineColor(kRed);
  h1_avgAmp_50um -> SetLineWidth(2);
  
  h1_avgAmp_25um -> SetLineColor(kBlue);
  h1_avgAmp_25um -> SetLineWidth(2);
  
  h1_avgAmp_25um -> Draw();
  h1_avgAmp_50um -> Draw("same");

  TLatex* latexLabel_50um = new TLatex(0.40,0.80,Form("RMS / mean = %.1f %%",h1_avgAmp_50um->GetRMS()/h1_avgAmp_50um->GetMean()*100));
  latexLabel_50um -> SetNDC();
  latexLabel_50um -> SetTextFont(42);
  latexLabel_50um -> SetTextSize(0.05);
  latexLabel_50um -> SetTextColor(kRed);
  latexLabel_50um -> Draw("same");
  
  TLatex* latexLabel_25um = new TLatex(0.40,0.70,Form("RMS / mean = %.1f %%",h1_avgAmp_25um->GetRMS()/h1_avgAmp_25um->GetMean()*100));
  latexLabel_25um -> SetNDC();
  latexLabel_25um -> SetTextFont(42);
  latexLabel_25um -> SetTextSize(0.05);
  latexLabel_25um -> SetTextColor(kBlue);
  latexLabel_25um -> Draw("same");
}




void drawCTRVsAvgAmp(const std::string& inFileName)
{
  TFile* inFile = TFile::Open(inFileName.c_str(),"READ");

  TGraphErrors* graph_25um = new TGraphErrors();
  TGraphErrors* graph_50um = new TGraphErrors();
  
  for(int ii = 1; ii <= 16; ++ii)
  {
    if( ii == 3 || ii == 4 ) continue;
    
    TH1F* histo = (TH1F*)( inFile->Get(Form("h_CTR_ampCorr_CH%02d",ii)));
    histo -> Fit("gaus","QLS+");
    TF1* fitFunc = (TF1*)( histo->GetFunction("gaus") );

    TH1F* histo2 = (TH1F*)( inFile->Get(Form("h_amp_cut_CH%02d",ii)));

    if( ii == 1 || ii == 2 || ii == 3 || ii == 4 || ii == 15 || ii == 16 )
    {
      graph_25um -> SetPoint(graph_25um->GetN(),histo2->GetMean(),fitFunc->GetParameter(2)*1000.);
      graph_25um -> SetPointError(graph_25um->GetN()-1,histo2->GetMeanError(),fitFunc->GetParError(2)*1000.);
    }
    else
    {
      graph_50um -> SetPoint(graph_50um->GetN(),histo2->GetMean(),fitFunc->GetParameter(2)*1000.);
      graph_50um -> SetPointError(graph_50um->GetN()-1,histo2->GetMeanError(),fitFunc->GetParError(2)*1000.);
    }
  }

  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,40.,0.6,80.) );
  hPad -> SetTitle(";avg. amplitude [V];#sigma_{t} [ps]");
  hPad -> Draw();
  
  graph_50um -> SetMarkerColor(kRed);
  graph_50um -> SetLineColor(kRed);
  graph_50um -> Draw("P,same");
  
  graph_25um -> SetMarkerColor(kBlue);
  graph_25um -> SetLineColor(kBlue);
  graph_25um -> Draw("P,same");
}
