#include "interface/FitUtils.h"



/*** find effective sigma ***/
void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction)
{
  float integralMax = fraction * histo->Integral();
  
  // find first and last bin with non-null content
  int N = histo -> GetNbinsX();
  int M1 = 1;
  int M2 = 1;
  for(int bin1 = 1; bin1 <= N; ++bin1)
  {
    if( histo->GetBinContent(bin1) > 0. && M1 == 1 ) M1 = bin1;
    if( histo->GetBinContent(bin1) > 0. )            M2 = bin1;
  }

  std::map<int,float> binCenters;
  std::map<int,float> binContents;
  std::map<int,float> binIntegrals;
  for(int bin1 = M1; bin1 <= M2; ++bin1)
  {
    if( histo->GetBinContent(bin1) == 0 ) continue;

    binCenters[bin1] = histo->GetBinCenter(bin1);
    binContents[bin1] = histo->GetBinContent(bin1);

    for(int bin2 = M1; bin2 <= bin1; ++bin2)
      binIntegrals[bin1] += histo->GetBinContent(bin2);
  }
  
  float min = 0.;
  float max = 0.;
  float delta = 999999.;
  for(std::map<int,float>::const_iterator mapIt1 = binIntegrals.begin(); mapIt1 != binIntegrals.end(); ++mapIt1)
  {
    for(std::map<int,float>::const_iterator mapIt2 = ++binIntegrals.begin(); mapIt2 != binIntegrals.end(); ++mapIt2)
    {
      if( (mapIt2->second-mapIt1->second) < integralMax ) continue;

      float tmpMin = binCenters[mapIt1->first];
      float tmpMax = binCenters[mapIt2->first];

      if( (tmpMax-tmpMin) < delta )
      {
        delta = (tmpMax - tmpMin);
        min = tmpMin;
        max = tmpMax;
      }

      break;
    }
  }
  
  TH1F* smallHisto = (TH1F*)( histo->Clone("smallHisto") );
  for(int bin = 1; bin <= smallHisto->GetNbinsX(); ++bin)
  {
    if( smallHisto->GetBinCenter(bin) < min )
      smallHisto -> SetBinContent(bin,0);

    if( smallHisto->GetBinCenter(bin) > max )
      smallHisto -> SetBinContent(bin,0);
  }
  smallHisto -> SetFillColor(kYellow);

  float mean = smallHisto -> GetMean();
  float meanErr = smallHisto -> GetMeanError();
  float RMS = smallHisto -> GetRMS();
  float RMSErr = smallHisto -> GetRMSError();

  ret[0] = mean;
  ret[1] = meanErr;
  ret[2] = RMS;
  ret[3] = RMSErr;
  ret[4] = min;
  ret[5] = max;

  std::cout << "mean: " << mean << "   min: " << min << "   max: " << max << std::endl;
}



CTRResult drawCTRPlot(TH1F* histo, const std::string& label, const int& rebin, const int& isMCP0, const int& isMCP1, const float& MCPIntrinsic,
                      const std::string& title, TLatex* latexLabel,
                      TH1F* histo_center, const std::string& center_label, TH1F* histo_border, const std::string& border_label)
{
  CTRResult result;
  
  float* vals = new float[6];
  FindSmallestInterval(vals,histo,0.68);
  histo = (TH1F*)( histo->Clone() );
  histo -> Scale(1./histo->Integral());
  histo -> Rebin(rebin);
  histo -> SetMarkerStyle(20);
  histo -> SetMarkerColor(kBlack);
  histo -> SetLineColor(kBlack);
  histo -> SetMaximum(1.5*histo->GetMaximum());
  histo -> SetTitle(title.c_str());
  histo -> Draw("PE");
  
  float norm = histo -> GetMaximum();
  float mean = vals[0];
  float min = vals[4];
  float max = vals[5];
  float delta = max-min;
  float sigma = 0.5*delta;
  float effSigma = sigma;
  min = min - 2.*delta;
  max = max + 2.*delta;
  
  histo -> GetXaxis() -> SetRangeUser(min,max);  
  
  TLatex* latex;
  if( (isMCP0 && !isMCP1) || (!isMCP0 && isMCP1) )
  {
    latex = new TLatex(0.20,0.85,Form("#splitline{%s:}{#sigma_{single}^{eff} = %.1f ps}",
                                      label.c_str(),sqrt(effSigma*effSigma - MCPIntrinsic*MCPIntrinsic)*1000.));
    result.effSigma = sqrt(effSigma*effSigma - MCPIntrinsic*MCPIntrinsic)*1000.;
  }
  else
  {
    latex = new TLatex(0.20,0.85,Form("#splitline{%s:}{#sigma_{single}^{eff} = %.1f ps}",label.c_str(),fabs(effSigma*1000)/sqrt(2)));
    result.effSigma = fabs(effSigma*1000)/sqrt(2);
  }
  latex -> SetNDC();
  latex -> SetTextFont(42);
  latex -> SetTextSize(0.05);
  latex -> SetTextColor(kBlack);
  latex -> Draw("same");
  
  // gaus fit
  std::string gausName(Form("fitFunc_gaus_corr"));
  TF1* fitFunc_gaus_corr = new TF1(gausName.c_str(),"gaus(0)",mean-1.5*sigma,mean+1.5*sigma);
  fitFunc_gaus_corr -> SetNpx(10000);
  fitFunc_gaus_corr -> SetParameters(norm,mean,sigma);
  fitFunc_gaus_corr -> SetLineColor(kBlack);
  histo -> Fit(gausName.c_str(),"QNRS");
  fitFunc_gaus_corr -> Draw("same");
  
  norm = fitFunc_gaus_corr -> GetParameter(0);
  mean = fitFunc_gaus_corr -> GetParameter(1);
  sigma = fitFunc_gaus_corr -> GetParameter(2);
  float sigmaErr = fitFunc_gaus_corr -> GetParError(2);
  
  TLatex* latex_gaus;
  if( (isMCP0 && !isMCP1) || (!isMCP0 && isMCP1) )
  {
    latex_gaus = new TLatex(0.20,0.75,Form("#sigma_{single}^{gaus} = (%.1f #pm %.1f) ps",
                                           sqrt(sigma*sigma - MCPIntrinsic*MCPIntrinsic)*1000.,
                                           fabs(sigmaErr*1000)));
    result.gausSigma = sqrt(sigma*sigma - MCPIntrinsic*MCPIntrinsic)*1000.;
    result.gausSigmaErr = sigmaErr*1000.;
  }
  else
  {
    latex_gaus = new TLatex(0.20,0.75,Form("#sigma_{single}^{gaus} = (%.1f #pm %.1f) ps",fabs(sigma*1000)/sqrt(2),fabs(sigmaErr*1000)/sqrt(2)));
    result.gausSigma = fabs(sigma*1000)/sqrt(2);
    result.gausSigmaErr = sigmaErr*1000.;
  }
  latex_gaus -> SetNDC();
  latex_gaus -> SetTextFont(42);
  latex_gaus -> SetTextSize(0.05);
  latex_gaus -> SetTextColor(kBlack);
  latex_gaus -> Draw("same");  
  
  
  if( histo_center)
  {
    vals = new float[6];
    FindSmallestInterval(vals,histo_center,0.68);
    histo_center = (TH1F*)( histo_center->Clone() );
    histo_center -> Scale(1./histo_center->Integral());
    histo_center -> Rebin(rebin);
    histo_center -> SetMarkerStyle(22);
    histo_center -> SetMarkerSize(0.7);
    histo_center -> SetMarkerColor(kRed);
    histo_center -> SetLineColor(kRed);
    histo_center -> Draw("PE,same");
    
    norm = histo -> GetMaximum();
    mean = vals[0];
    min = vals[4];
    max = vals[5];
    delta = max-min;
    sigma = 0.5*delta;
    effSigma = sigma;
    min = min - 2.*delta;
    max = max + 2.*delta;
    
    TLatex* latex_center;
    if( (isMCP0 && !isMCP1) || (!isMCP0 && isMCP1) )
      latex_center = new TLatex(0.75,0.40,Form("#splitline{%s:}{#sigma_{single}^{eff} = %.1f ps}",
                                               center_label.c_str(),sqrt(effSigma*effSigma - MCPIntrinsic*MCPIntrinsic)*1000.));
    else
      latex_center = new TLatex(0.75,0.40,Form("#splitline{%s:}{#sigma_{single}^{eff} = %.1f ps}",center_label.c_str(),fabs(effSigma*1000)/sqrt(2)));
    latex_center -> SetNDC();
    latex_center -> SetTextFont(42);
    latex_center -> SetTextSize(0.03);
    latex_center -> SetTextColor(kRed);
    latex_center -> Draw("same");  
  }
  
  
  if( histo_border)
  {
    vals = new float[6];
    FindSmallestInterval(vals,histo_border,0.68);
    histo_border -> Scale(1./histo_border->Integral());
    histo_border = (TH1F*)( histo_border->Clone() );
    histo_border -> Rebin(rebin);
    histo_border -> SetMarkerStyle(23);
    histo_border -> SetMarkerSize(0.7);
    histo_border -> SetMarkerColor(kBlue);
    histo_border -> SetLineColor(kBlue);
    histo_border -> Draw("PE,same");
    
    norm = histo -> GetMaximum();
    mean = vals[0];
    min = vals[4];
    max = vals[5];
    delta = max-min;
    sigma = 0.5*delta;
    effSigma = sigma;
    min = min - 2.*delta;
    max = max + 2.*delta;
    
    TLatex* latex_border;
    if( (isMCP0 && !isMCP1) || (!isMCP0 && isMCP1) )
      latex_border = new TLatex(0.75,0.30,Form("#splitline{%s:}{#sigma_{single}^{eff} = %.1f ps}",
                                               border_label.c_str(),sqrt(effSigma*effSigma - MCPIntrinsic*MCPIntrinsic)*1000.));
    else
      latex_border = new TLatex(0.75,0.30,Form("#splitline{%s:}{#sigma_{single}^{eff} = %.1f ps}",border_label.c_str(),fabs(effSigma*1000)/sqrt(2)));
    latex_border -> SetNDC();
    latex_border -> SetTextFont(42);
    latex_border -> SetTextSize(0.03);
    latex_border -> SetTextColor(kBlue);
    latex_border -> Draw("same");  
  }
  
  
  latexLabel -> Draw("same");
  
  
  gPad -> Update();

  return result;
}



void PrintCanvas(TCanvas* c, CfgManager& opts, const std::string& ch, const std::string& plotDir, const std::string& label)
{
  std::string shortLabel = opts.GetOpt<std::string>(Form("%s.shortLabel",ch.c_str()));
  
  c -> Print(Form("%s/c_%s_%s.png",plotDir.c_str(),label.c_str(),shortLabel.c_str()));
  c -> Print(Form("%s/c_%s_%s.pdf",plotDir.c_str(),label.c_str(),shortLabel.c_str()));
}



void DrawHistogram(CfgManager& opts,
                   const std::string& ch, TH1F* histo,
                   const std::string& title,
                   const float& xMin, const float& xMax, const int& rebin, const bool& logy,
                   TLatex* latexLabel, std::vector<TLine*>* lines)
{
  std::string shortLabel = opts.GetOpt<std::string>(Form("%s.shortLabel",ch.c_str()));
  
  gPad -> SetGridx();
  if( logy ) gPad -> SetLogy();
  
  histo -> SetTitle(title.c_str());
  histo -> Rebin(rebin);
  histo -> SetLineColor(kRed);
  histo -> SetFillColor(kRed);
  histo -> SetFillStyle(3003);
  histo -> GetXaxis() -> SetRangeUser(xMin,xMax);
  histo -> Draw("hist");

  if( latexLabel != NULL ) latexLabel -> Draw("same");

  if( lines != NULL )
  {
    for(unsigned int ii = 0; ii < lines->size(); ++ii)
    {
      lines->at(ii) -> SetLineStyle(2);
      lines->at(ii) -> SetLineWidth(3);
      lines->at(ii) -> Draw("same");
    }
  }

  gPad -> Update();
}



void DrawHistogram2D(CfgManager& opts,
                     const std::string& ch, TH2F* histo2,
                     const std::string& title,
                     const float& xMin, const float& xMax, const float& yMin, const float& yMax, const float& zMin, const float& zMax,
                     const bool& logx, const bool& logy, const std::string& drawOpt,
                     TLatex* latexLabel)
{
  std::string shortLabel = opts.GetOpt<std::string>(Form("%s.shortLabel",ch.c_str()));
  
  gPad -> SetGridx();
  if( logx ) gPad -> SetLogy();
  if( logy ) gPad -> SetLogy();
  
  histo2 -> SetTitle(title.c_str());
  histo2 -> SetLineColor(kRed);
  histo2 -> SetFillColor(kRed);
  histo2 -> SetFillStyle(3003);
  histo2 -> GetXaxis() -> SetRangeUser(xMin,xMax);
  histo2 -> GetYaxis() -> SetRangeUser(yMin,yMax);
  histo2 -> SetMinimum(zMin);
  histo2 -> SetMaximum(zMax);
  histo2 -> Draw(drawOpt.c_str());
  
  if( latexLabel != NULL ) latexLabel -> Draw("same");
  
  TH1F* h_spread = new TH1F("h_spread","",100,zMin,zMax);
  for(int ii = 1; ii <= histo2->GetNbinsX(); ++ii)
    for(int jj = 1; jj <= histo2->GetNbinsY(); ++jj)
    {
      float val = histo2 -> GetBinContent(ii,jj);
      h_spread -> Fill(val);
    }
  TLatex* latex_mean = new TLatex(0.19,0.83,Form("mean = %.2e",h_spread->GetMean()));  
  latex_mean -> SetNDC();
  latex_mean -> SetTextFont(42);
  latex_mean -> SetTextSize(0.04);
  latex_mean -> Draw("same");
  TLatex* latex_rms = new TLatex(0.19,0.78,Form("rms = %.3e",h_spread->GetRMS()));
  latex_rms -> Draw("same");
  latex_rms -> SetNDC();
  latex_rms -> SetTextFont(42);
  latex_rms -> SetTextSize(0.04);
  delete h_spread;
  
  if( latexLabel != NULL ) latexLabel -> Draw("same");
  
  gPad -> Update();
}



void DrawProfile(CfgManager& opts,
                 const std::string& ch, TProfile* prof,
                 const std::string& title,
                 const float& xMin, const float& xMax, const float& yMin, const float& yMax,
                 const int& color, const std::string& drawOpt,
                 TLatex* latexLabel, TF1* func)
{
  std::string shortLabel = opts.GetOpt<std::string>(Form("%s.shortLabel",ch.c_str()));
  
  gPad -> SetGridy();
  
  prof -> SetTitle(title.c_str());
  prof -> SetMarkerColor(color);
  prof -> SetLineColor(color);
  prof -> SetMarkerSize(0.3);
  prof -> GetXaxis() -> SetRangeUser(xMin,xMax);
  prof -> SetMinimum(yMin);
  prof -> SetMaximum(yMax);
  prof -> Draw(drawOpt.c_str());

  if( latexLabel != NULL ) latexLabel -> Draw("same");
  
  if( func != NULL )
  {
    func -> SetLineStyle(2);
    func -> SetLineWidth(2);
    func -> SetLineColor(kRed);
    func -> Draw("same");
  }
  
  gPad -> Update();
}



void DrawProfile2D(CfgManager& opts,
                   const std::string& ch, TProfile2D* prof2,
                   const std::string& title,
                   const float& xMin, const float& xMax, const float& yMin, const float& yMax, const float& zMin, const float& zMax,
                   TLatex* latexLabel, std::vector<TLine*>* lines)
{
  std::string shortLabel = opts.GetOpt<std::string>(Form("%s.shortLabel",ch.c_str()));
  
  prof2 -> SetTitle(title.c_str());
  prof2 -> GetXaxis() -> SetRangeUser(xMin,xMax);
  prof2 -> GetYaxis() -> SetRangeUser(yMin,yMax);
  prof2 -> SetMinimum(zMin);
  prof2 -> SetMaximum(zMax);
  prof2 -> Draw("COLZ");
  
  TH1F* h_spread = new TH1F("h_spread","",100,zMin,zMax);
  for(int ii = 1; ii <= prof2->GetNbinsX(); ++ii)
    for(int jj = 1; jj <= prof2->GetNbinsY(); ++jj)
    {
      float val = prof2 -> GetBinContent(ii,jj);
      h_spread -> Fill(val);
    }
  TLatex* latex_mean = new TLatex(0.19,0.83,Form("mean = %.1e",h_spread->GetMean()));  
  latex_mean -> SetNDC();
  latex_mean -> SetTextFont(42);
  latex_mean -> SetTextSize(0.04);
  latex_mean -> Draw("same");
  TLatex* latex_rms = new TLatex(0.19,0.78,Form("rms = %.2e",h_spread->GetRMS()));
  latex_rms -> Draw("same");
  latex_rms -> SetNDC();
  latex_rms -> SetTextFont(42);
  latex_rms -> SetTextSize(0.04);
  delete h_spread;
  
  if( latexLabel != NULL ) latexLabel -> Draw("same");
  
  if( lines != NULL )
  {
    for(unsigned int ii = 0; ii < lines->size(); ++ii)
      lines->at(ii) -> Draw("same");
  }
  
  gPad -> Update();
}
