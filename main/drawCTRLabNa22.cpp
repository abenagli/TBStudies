#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <vector>
#include <map>

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"



float GetFWHM(const TF1* func, const float& xMin, const float& xMax)
{
  float yMax = func->GetMaximum();
  float yMax2 = 0.5 * yMax;

  float left  = -1.;
  float right = -1.;
  for(int ii = 0; ii < 10000; ++ii)
  {
    double xx = xMin + (xMax-xMin)/10000.*ii;
    if( func->Eval(xx) > yMax2 && left == -1.) left = xx;
    if( func->Eval(xx) < yMax2 && left != -1. && right == -1.) { right = xx; break; }
  }

  return (right-left);
}



int main(int argc, char** argv)
{
  setTDRStyle();

  
  if( argc < 2 )
  {
    std::cout << ">>> drawCTRLabNa22::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  
  
  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  

  
  //--- get parameters
  std::string conf = opts.GetOpt<std::string>("Input.conf");
  
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  plotDir += "/" + conf;
  system(Form("mkdir -p %s",plotDir.c_str()));
  
  std::vector<std::string> channels = opts.GetOpt<std::vector<std::string> >("Channels.channels");
  
  
  
  //-------------------------
  // SINGLE P.E. CALIBRATION
  //-------------------------
  
  //--- open input files
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");
  
  std::vector<std::string> runs_peCalib = opts.GetOpt<std::vector<std::string> >("Input.runs_peCalib");
  TChain* h4_peCalib = new TChain("h4","h4");
  for(auto run : runs_peCalib)
  {
    std::string fileName = Form("%s/lab5015_%d.root",inputDir.c_str(),atoi(run.c_str()));
    std::cout << ">>> Adding flle " << fileName << std::endl;
    h4_peCalib -> Add(fileName.c_str());
  }
  
  
  
  //--- define tree branch addresses
  std::map<std::string,int> channelIds;
  std::map<std::string,int> timeMethodIds;
  
  for(auto ch : channels)
  {
    std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
    if( timeCh != "NULL" ) h4_peCalib -> SetBranchAddress(timeCh.c_str(), &channelIds[timeCh.c_str()]);
  }
  
  float b_rms;
  float* time    = new float[1000];
  float* time_slope = new float[1000];
  float* amp_max = new float[1000];
  float* chi2_max = new float[1000];
  h4_peCalib -> SetBranchAddress("time",    time);
  h4_peCalib -> SetBranchAddress("amp_max", amp_max);
  h4_peCalib -> SetBranchAddress("chi2_max",chi2_max);
  
  int nEntries = h4_peCalib->GetEntries();
  std::cout << ">>> Events read: " << nEntries << std::endl;
  
  
  
  //------------------
  // define histograms
  
  TFile* outFile = TFile::Open(Form("%s/drawCTRLabNa22.root",plotDir.c_str()),"RECREATE");
  outFile -> cd();
  
  std::map<std::string,TH1F*> h_peCalib;
  std::map<std::string,TH1F*> h_peCalib_chi2Cut;  
  for(auto ch : channels)
  {
    h_peCalib[ch]         = new TH1F(Form("h_peCalib_%s",        ch.c_str()),"",8000,0.,1000.);
    h_peCalib_chi2Cut[ch] = new TH1F(Form("h_peCalib_chi2Cut_%s",ch.c_str()),"",8000,0.,1000.);
  }
  


  //-----------------------
  // 1st loop over events
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    h4_peCalib -> GetEntry(entry);
    
    
    for(auto ch: channels)
    {
      // fill amplitude and time plots
      std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
      
      float amp = amp_max[channelIds[timeCh]] / 4096. * 1000.;
      float chi2 = chi2_max[channelIds[timeCh]];
      h_peCalib[ch] -> Fill( amp );
      if( chi2 < 1.5 ) h_peCalib_chi2Cut[ch] -> Fill( amp );
    }
  }
  
  
  //--- draw plots
  TCanvas* c;
  TH1F* histo;
  TH1F* histoCorr;
  
  std::map<std::string,double> peCalib_higGain;
  std::map<std::string,double> peCalibErr_higGain;
  
  std::map<std::string,TLatex*> latexLabels;
  for(auto ch: channels)
  {
    std::vector<std::string> labelVec = opts.GetOpt<std::vector<std::string> >(Form("%s.label",ch.c_str()));
    std::string label = "";
    for(auto it: labelVec) label += it + " ";
    TLatex* latexLabel = new TLatex(0.16,0.96,Form("%s",label.c_str()));
    latexLabel -> SetNDC();
    latexLabel -> SetTextFont(42);
    latexLabel -> SetTextSize(0.03);
    latexLabels[ch] = latexLabel;
  }
  
  c = new TCanvas("c_peCalib","c_peCalib",2000,2000);
  c -> Divide(2,2);
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    c -> cd(index);
    gPad -> SetLogy();
    histo = h_peCalib_chi2Cut[ch];
    histo -> Scale( 1./histo->Integral() );
    histo -> SetTitle(Form(";SiPM%d amplitude (high gain) [mV];event fraction",index));
    histo -> SetLineWidth(2);
    histo -> SetLineColor(kRed+1);
    histo -> SetFillColor(kRed+1);
    histo -> SetFillStyle(3003);
    histo -> GetXaxis() -> SetRangeUser(0.,40.);
    histo -> Draw("hist");
    
    latexLabels[ch] -> Draw("same");

    TF1* g1pe = new TF1("g1pe","gaus",1.,200.); g1pe -> SetParameters(0.01,2.5,0.7);
    g1pe -> SetNpx(10000);
    g1pe -> SetLineColor(kBlack);
    g1pe -> SetLineStyle(7);
    histo -> Fit(g1pe,"QNRS+","",2.,5.);
    g1pe -> Draw("same");

    TF1* gbkg = new TF1("gbkg","gaus",1.,200.); gbkg -> SetParameters(0.01,10.,10.);
    gbkg -> SetNpx(10000);
    gbkg -> SetLineColor(kBlack);
    gbkg -> SetLineStyle(7);
    histo -> Fit(gbkg,"QNRS+","",2.,50.);
    // gbkg -> Draw("same");
    
    TF1* total_temp = new TF1("total","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)",1.,200.);
    total_temp -> SetNpx(10000);
    total_temp -> SetParameter(0,0.01);
    total_temp -> SetParameter(1,gbkg->GetParameter(1));
    total_temp -> SetParameter(2,gbkg->GetParameter(2));
    for(int ipe = 1; ipe <= 6; ++ipe)
    {
      total_temp -> SetParameter(0+3*(ipe),0.01);
      total_temp -> FixParameter(1+3*(ipe),g1pe->GetParameter(1)*ipe);
      total_temp -> FixParameter(2+3*(ipe),g1pe->GetParameter(2));
    }
    histo -> Fit(total_temp,"QNRS+");
    
    TF1* total = new TF1("total","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)",1.,200.);
    total -> SetNpx(10000);
    total -> SetLineColor(kBlue);
    total -> SetLineWidth(2);
    for(int ipar = 0; ipar < total_temp->GetNpar(); ++ipar)
    {
      total -> SetParameter(ipar,total_temp->GetParameter(ipar));
    }
    histo -> Fit(total,"QNRS+");
    total -> Draw("same");

    for(int ipe = 1; ipe <= 5; ++ipe)
    {
      TLatex* latex = new TLatex(0.45,0.85-0.05*(ipe-1),Form("%d p.e. = (%6.3f #pm %.3f) mV",ipe,total->GetParameter(4+3*(ipe-1)),total->GetParError(4+3*(ipe-1))));
      latex -> SetNDC();
      latex -> SetTextFont(82);
      latex -> SetTextSize(0.03);
      latex -> Draw("same");
    }

    
    c -> cd(index+2);
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,6.,6*g1pe->GetParameter(1)) );
    hPad -> SetTitle(Form(";n_{p.e.};SiPM%d amplitude (high gain) [mV]",index));
    hPad -> Draw();

    TGraphErrors* g = new TGraphErrors();
    for(int ipe = 1; ipe <= 5; ++ipe)
    {
      g -> SetPoint(ipe-1,ipe,total->GetParameter(4+3*(ipe-1)));
      g -> SetPointError(ipe-1,0.,total->GetParError(4+3*(ipe-1)));
    }
    g -> SetMarkerStyle(20);
    g -> SetMarkerSize(1.5);
    g -> Draw("P,same");
    
    TF1* calib = new TF1("calib","[0]*x",0.,6.);
    calib -> SetNpx(10000);
    calib -> SetLineColor(kRed+1);
    g -> Fit(calib,"QNRS+","",1.5,6.);
    calib -> Draw("same");
    peCalib_higGain[ch] = calib->GetParameter(0);
    peCalibErr_higGain[ch] = calib->GetParError(0);
    
    TLatex* latex = new TLatex(0.20,0.85,Form("(%.3f #pm %.3f) mV/p.e.",calib->GetParameter(0),calib->GetParError(0)));
    latex -> SetNDC();
    latex -> SetTextFont(82);
    latex -> SetTextSize(0.04);
    latex -> Draw("same");
    
    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();
  }
  c -> Print(Form("%s/c_peCalib.png",plotDir.c_str()));
  c -> Print(Form("%s/c_peCalib.pdf",plotDir.c_str()));
  
  
  
  //-----------------------
  // AMPLITUDE CALIBRATION
  //-----------------------
  
  std::vector<std::string> runs_ampCalib = opts.GetOpt<std::vector<std::string> >("Input.runs_ampCalib");
  TChain* h4_ampCalib = new TChain("h4","h4");
  for(auto run : runs_ampCalib)
  {
    std::string fileName = Form("%s/lab5015_%d.root",inputDir.c_str(),atoi(run.c_str()));
    std::cout << ">>> Adding flle " << fileName << std::endl;
    h4_ampCalib -> Add(fileName.c_str());
  }
  
  
  
  //--- define tree branch addresses
  channelIds.clear();
  timeMethodIds.clear();
  
  for(auto ch : channels)
  {
    std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
    if( ampCh != "NULL" ) h4_ampCalib -> SetBranchAddress(ampCh.c_str(), &channelIds[ampCh.c_str()]);
    
    std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
    if( timeCh != "NULL" ) h4_ampCalib -> SetBranchAddress(timeCh.c_str(), &channelIds[timeCh.c_str()]);
  }
  
  h4_ampCalib -> SetBranchAddress("time",    time);
  h4_ampCalib -> SetBranchAddress("amp_max", amp_max);
  h4_ampCalib -> SetBranchAddress("chi2_max",chi2_max);
  
  nEntries = h4_ampCalib->GetEntries();
  std::cout << ">>> Events read: " << nEntries << std::endl;
  
  
  //------------------
  // define histograms
  std::map<std::string,TH1F*> h_ampCalib;
  std::map<std::string,TH1F*> h_ampCalib_chi2Cut;
  std::map<std::string,TH1F*> h_ampCalib_chi2_lowGain;
  std::map<std::string,TH1F*> h_ampCalib_chi2_higGain;
  for(auto ch : channels)
  {
    h_ampCalib[ch]              = new TH1F(Form("h_ampCalib_%s",             ch.c_str()),"",1000,0.,100.);
    h_ampCalib_chi2Cut[ch]      = new TH1F(Form("h_ampCalib_chi2Cut_%s",     ch.c_str()),"",1000,0.,100.);
    h_ampCalib_chi2_lowGain[ch] = new TH1F(Form("h_ampCalib_chi2_lowGain_%s",ch.c_str()),"",1000,0.,100.);
    h_ampCalib_chi2_higGain[ch] = new TH1F(Form("h_ampCalib_chi2_higGain_%s",ch.c_str()),"",1000,0.,100.);
  }
  


  //-----------------------
  // 1st loop over events
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    h4_ampCalib -> GetEntry(entry);
    
    
    for(auto ch: channels)
    {
      // fill amplitude and time plots
      std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
      std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
      
      float amp     = amp_max[channelIds[ampCh]]  / 4096. * 1000.;
      float ampTime = amp_max[channelIds[timeCh]] / 4096. * 1000.;
      float chi2     = chi2_max[channelIds[ampCh]];
      float chi2Time = chi2_max[channelIds[timeCh]];

      h_ampCalib_chi2_lowGain[ch] -> Fill(chi2);
      h_ampCalib_chi2_higGain[ch] -> Fill(chi2Time);
      
      h_ampCalib[ch] -> Fill( ampTime/amp );
      if( chi2 < 4. && chi2Time < 4. ) h_ampCalib_chi2Cut[ch] -> Fill( ampTime/amp );
    }
  }
  

  //--- legend
  TLegend* legend = new TLegend(0.60,0.70,0.90,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(42);
  
  //--- draw plots
  std::map<std::string,double> peCalib_lowGain;
  std::map<std::string,double> peCalibErr_lowGain;
  
  c = new TCanvas("c_ampCalib","c_ampCalib",2000,2000);
  c -> Divide(2,2);

  int chIt = 0;
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    c -> cd(index);
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(15.,0.1,100.,h_ampCalib[ch]->Integral()) );
    hPad -> SetTitle(Form(";SiPM%d (high gain / low gain) ratio;event fraction",index));
    hPad -> Draw();
    gPad -> SetLogy();
    histo = h_ampCalib[ch];
    histo -> SetLineWidth(1);
    histo -> SetLineColor(17);
    histo -> SetFillColor(17);
    histo -> SetFillStyle(3003);
    histo -> Draw("hist,same");
    histo = h_ampCalib_chi2Cut[ch];
    histo -> SetLineWidth(2);
    histo -> SetLineColor(12);
    histo -> SetFillColor(12);
    histo -> SetFillStyle(3001);
    histo -> Draw("hist,same");
    
    latexLabels[ch] -> Draw("same");
    
    TF1* total = new TF1("total","gaus(0)",1.,200.);
    total -> SetNpx(10000);
    total -> SetLineColor(kBlue);
    total -> SetLineWidth(2);
    histo -> Fit(total,"QNRS+","",20.,40.);
    total -> Draw("same");
    peCalib_lowGain[ch] = peCalib_higGain[ch]/total->GetParameter(1);
    peCalibErr_lowGain[ch] = peCalib_lowGain[ch]*sqrt( pow(peCalibErr_higGain[ch]/peCalib_higGain[ch],2) + pow(total->GetParError(1)/total->GetParameter(1),2) );
  
    TLatex* latex = new TLatex(0.20,0.85,Form("R = (%.2f #pm %.2f)",total->GetParameter(1),total->GetParError(1)));
    latex -> SetNDC();
    latex -> SetTextFont(82);
    latex -> SetTextSize(0.04);
    latex -> Draw("same");

    TLatex* latex2 = new TLatex(0.20,0.75,Form("(%.4f #pm %.4f) mV/p.e.",peCalib_lowGain[ch],peCalibErr_lowGain[ch]));
    latex2 -> SetNDC();
    latex2 -> SetTextFont(82);
    latex2 -> SetTextSize(0.04);
    latex2 -> Draw("same");

    
    c -> cd(index+2);
    histo = h_ampCalib_chi2_higGain[ch];
    histo -> Scale( 1./histo->Integral() );
    histo -> SetTitle(Form(";SiPM%d amplitude fit #chi_{2} / N_{d.o.f.};event fraction",index));
    histo -> SetLineWidth(2);
    histo -> SetLineColor(kRed+1);
    histo -> SetFillColor(kRed+1);
    histo -> SetFillStyle(3003);
    histo -> GetXaxis() -> SetRangeUser(0.,15.);
    histo -> Draw("hist");
    if( chIt == 0 ) legend -> AddEntry(histo,"high gain channel","F");
    
    histo = h_ampCalib_chi2_lowGain[ch];
    histo -> Scale( 1./histo->Integral() );
    histo -> SetLineWidth(2);
    histo -> SetLineColor(kGreen-5);
    histo -> SetFillColor(kGreen-5);
    histo -> SetFillStyle(3003);
    histo -> Draw("hist,same");
    if( chIt == 0 ) legend -> AddEntry(histo,"low gain channel","F");
    
    latexLabels[ch] -> Draw("same");
    legend -> Draw("same");
    
    TLine* line_chi2Cut = new TLine(4.,h_ampCalib_chi2_higGain[ch]->GetMinimum(),4.,2.*h_ampCalib_chi2_higGain[ch]->GetMaximum());
    line_chi2Cut -> SetLineWidth(2);
    line_chi2Cut -> SetLineStyle(7);
    line_chi2Cut -> SetLineColor(kBlack);
    line_chi2Cut -> Draw("same");
    
    gPad -> SetLogy();
    gPad -> Update();

    ++chIt;
  }
  c -> Print(Form("%s/c_ampCalib.png",plotDir.c_str()));
  c -> Print(Form("%s/c_ampCalib.pdf",plotDir.c_str()));
  
  
  
  //----------------
  // SINGLES SPECTRA
  //----------------
  
  std::vector<std::string> runs_singles = opts.GetOpt<std::vector<std::string> >("Input.runs_singles");
  TChain* h4_singles = new TChain("h4","h4");
  for(auto run : runs_singles)
  {
    std::string fileName = Form("%s/lab5015_%d.root",inputDir.c_str(),atoi(run.c_str()));
    std::cout << ">>> Adding flle " << fileName << std::endl;
    h4_singles -> Add(fileName.c_str());
  }
  
  
  
  //--- define tree branch addresses
  channelIds.clear();
  timeMethodIds.clear();
  
  for(auto ch : channels)
  {
    std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
    if( ampCh != "NULL" ) h4_singles -> SetBranchAddress(ampCh.c_str(), &channelIds[ampCh.c_str()]);
    
    std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
    if( timeCh != "NULL" ) h4_singles -> SetBranchAddress(timeCh.c_str(), &channelIds[timeCh.c_str()]);
  }
  
  h4_singles -> SetBranchAddress("time",    time);
  h4_singles -> SetBranchAddress("amp_max", amp_max);
  h4_singles -> SetBranchAddress("chi2_max",chi2_max);
  
  nEntries = h4_singles->GetEntries();
  std::cout << ">>> Events read: " << nEntries << std::endl;
  
  
  //------------------
  // define histograms
  std::map<std::string,TH1F*> h_singles;
  std::map<std::string,TH1F*> h_singles_chi2Cut;
  std::map<std::string,TH1F*> h_singles_chi2_lowGain;
  std::map<std::string,TH1F*> h_singles_chi2_higGain;
  for(auto ch : channels)
  {
    h_singles[ch]      = new TH1F(Form("h_singles_%s",      ch.c_str()),"",500,0.,4000.);
    h_singles_chi2Cut[ch] = new TH1F(Form("h_singles_chi2Cut_%s", ch.c_str()),"",500,0.,4000.);
    h_singles_chi2_lowGain[ch] = new TH1F(Form("h_singles_chi2_lowGain_%s",ch.c_str()),"",1000,0.,100.);
  }
  


  //-----------------------
  // 1st loop over events
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    h4_singles -> GetEntry(entry);
    
    
    for(auto ch: channels)
    {
      // fill amplitude and time plots
      std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
      
      float amp  = amp_max[channelIds[ampCh]]  / 4096. * 1000.;
      float chi2 = chi2_max[channelIds[ampCh]];
      
      h_singles_chi2_lowGain[ch] -> Fill( chi2 );
      
      h_singles[ch] -> Fill( amp/peCalib_lowGain[ch] );
      if( chi2 < 5. ) h_singles_chi2Cut[ch] -> Fill( amp/peCalib_lowGain[ch] );
    }
  }
  
  
  //--- draw plots
  c = new TCanvas("c_singles","c_singles",2000,2000);
  c -> Divide(2,2);
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    c -> cd(index);
    
    histo = h_singles[ch];
    histo -> SetTitle(Form(";SiPM%d amplitude (low gain) [n_{p.e.}];events",index));
    histo -> SetLineWidth(2);
    histo -> SetLineColor(kGreen-5);
    histo -> SetFillColor(kGreen-5);
    histo -> SetFillStyle(3003);
    histo -> GetXaxis() -> SetRangeUser(250.,3500.);
    histo -> Draw("hist");

    histo = h_singles_chi2Cut[ch];
    histo -> SetLineWidth(2);
    histo -> SetLineColor(kGreen-8);
    histo -> SetFillColor(kGreen-8);
    histo -> SetFillStyle(3001);
    histo -> Draw("hist,same");
    
    latexLabels[ch] -> Draw("same");

    histo = h_singles[ch];
    float xMin  = opts.GetOpt<float>(Form("%s.ampMin", ch.c_str()));
    float xMax  = opts.GetOpt<float>(Form("%s.ampMax", ch.c_str()));
    TF1* total = new TF1("total","gaus(0)+gaus(3)+gaus(6)",xMin,xMax);
    total -> SetParameters(1000.,0.5*(xMin+xMax),40.,100.,0.5*(xMin+xMax),80.,100.,0.5*(xMin+xMax),120.);
    total -> SetNpx(10000);
    total -> SetLineColor(kBlue);
    total -> SetLineWidth(2);
    histo -> Fit(total,"QNRS+");
    total -> Draw("same");
    
    TLatex* latex = new TLatex(0.50,0.85,Form("photopeak: %.0f p.e.",total->GetMaximumX()));
    latex -> SetNDC();
    latex -> SetTextFont(82);
    latex -> SetTextSize(0.04);
    latex -> Draw("same");

    float FWHM = GetFWHM(total,xMin,xMax);
    TLatex* latex2 = new TLatex(0.50,0.75,Form("FWHM: %.1f%%",FWHM/total->GetMaximumX()*100.));
    latex2 -> SetNDC();
    latex2 -> SetTextFont(82);
    latex2 -> SetTextSize(0.04);
    latex2 -> Draw("same");

    TLine* line_FWHM = new TLine(total->GetMaximumX()-0.5*FWHM,0.5*total->GetMaximum(),total->GetMaximumX()+0.5*FWHM,0.5*total->GetMaximum());
    line_FWHM -> SetLineColor(kBlue);
    line_FWHM -> SetLineWidth(2);
    line_FWHM -> Draw("same");
    
    gPad -> SetLogy();
    
    c -> cd(index+2);
    histo = h_singles_chi2_lowGain[ch];
    histo -> Scale( 1./histo->Integral() );
    histo -> SetTitle(Form(";SiPM%d #chi_{2} / N_{d.o.f.} (low gain amp. fit);event fraction",index));
    histo -> SetLineWidth(2);
    histo -> SetLineColor(kGreen-5);
    histo -> SetFillColor(kGreen-5);
    histo -> SetFillStyle(3003);
    histo -> GetXaxis() -> SetRangeUser(0.,30.);
    histo -> Draw("hist");

    TLine* line_chi2Cut = new TLine(5.,histo->GetMinimum(),5.,2.*histo->GetMaximum());
    line_chi2Cut -> SetLineWidth(2);
    line_chi2Cut -> SetLineStyle(7);
    line_chi2Cut -> SetLineColor(kBlack);
    line_chi2Cut -> Draw("same");
    
    latexLabels[ch] -> Draw("same");

    gPad -> SetLogy();
    gPad -> Update();
  }
  c -> Print(Form("%s/c_singles.png",plotDir.c_str()));
  c -> Print(Form("%s/c_singles.pdf",plotDir.c_str()));
  
  
  
  //---------------------
  // COINCIDENCES SPECTRA
  //---------------------
  
  std::vector<std::string> runs_coincidences = opts.GetOpt<std::vector<std::string> >("Input.runs_coincidences");
  TChain* h4_coincidences = new TChain("h4","h4");
  for(auto run : runs_coincidences)
  {
    std::string fileName = Form("%s/lab5015_%d.root",inputDir.c_str(),atoi(run.c_str()));
    std::cout << ">>> Adding flle " << fileName << std::endl;
    h4_coincidences -> Add(fileName.c_str());
  }
  
  
  
  //--- define tree branch addresses
  channelIds.clear();
  timeMethodIds.clear();

  std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("Input.timeMethods"));
  
  for(auto ch : channels)
  {
    std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
    if( ampCh != "NULL" ) h4_coincidences -> SetBranchAddress(ampCh.c_str(), &channelIds[ampCh.c_str()]);
    
    std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
    if( timeCh != "NULL" ) h4_coincidences -> SetBranchAddress(timeCh.c_str(), &channelIds[timeCh.c_str()]);
    
    for(auto timeMethod: timeMethods)
      if( timeMethod != "NULL" ) h4_coincidences -> SetBranchAddress(timeMethod.c_str(), &timeMethodIds[timeMethod.c_str()]);
  }
  
  h4_coincidences -> SetBranchAddress("b_rms",    &b_rms);
  h4_coincidences -> SetBranchAddress("time",      time);
  h4_coincidences -> SetBranchAddress("time_slope",time_slope);
  h4_coincidences -> SetBranchAddress("amp_max",   amp_max);
  h4_coincidences -> SetBranchAddress("chi2_max",  chi2_max);
  
  nEntries = h4_coincidences->GetEntries();
  std::cout << ">>> Events read: " << nEntries << std::endl;
  
  
  //------------------
  // define histograms
  std::map<std::string,TH1F*> h_coincidences;
  std::map<std::string,TH1F*> h_coincidences_chi2Cut;
  std::map<std::string,TH1F*> h_coincidences_chi2_lowGain;
  std::map<std::string,TH1F*> h_coincidences_chi2_higGain;

  TH1F* h_ampRatio = new TH1F("h_ampRatio","",1000,0.,2.);
  
  for(auto ch : channels)
  {
    h_coincidences[ch]      = new TH1F(Form("h_coincidences_%s",      ch.c_str()),"",500,0.,4000.);
    h_coincidences_chi2Cut[ch] = new TH1F(Form("h_coincidences_chi2Cut_%s", ch.c_str()),"",500,0.,4000.);
    h_coincidences_chi2_lowGain[ch] = new TH1F(Form("h_coincidences_chi2_lowGain_%s",ch.c_str()),"",1000,0.,100.);
  }
  
  std::map<std::string,std::map<std::string,TH1F*> > h_time;  
  for(auto ch: channels)
    for(auto timeMethod : timeMethods)
    {
      h_time[ch][timeMethod] = new TH1F(Form("h_time_%s_%s",ch.c_str(),timeMethod.c_str()),"",1000,0.,200.);
    }

  std::map<std::string,TH1F*> h_CTR_raw;
  std::map<std::string,TProfile*> p_CTR_vs_ampRatio;
  for(auto timeMethod : timeMethods)
  {
    h_CTR_raw[timeMethod] = new TH1F(Form("h_CTR_raw_%s",timeMethod.c_str()),"",10000,-50.,50.);
    p_CTR_vs_ampRatio[timeMethod] = new TProfile(Form("p_CTR_vs_ampRatio_%s",timeMethod.c_str()),"",1000,0.,2.);
  }
  
  
  
  //-----------------------
  // 1st loop over events
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    h4_coincidences -> GetEntry(entry);
    
    
    for(auto ch: channels)
    {
      // fill amplitude and time plots
      std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
      std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
      float ampMin  = opts.GetOpt<float>(Form("%s.ampMin", ch.c_str()));
      float ampMax  = opts.GetOpt<float>(Form("%s.ampMax", ch.c_str()));
      
      float amp  = amp_max[channelIds[ampCh]]  / 4096. * 1000. / peCalib_lowGain[ch];
      float chi2 = chi2_max[channelIds[ampCh]];

      for(auto timeMethod : timeMethods)
      {
        float tim = time[channelIds[timeCh]+timeMethodIds[timeMethod]];
        
        if( amp >= ampMin  && amp < ampMax )
        {
          h_time[ch][timeMethod] -> Fill( tim );
        }
      }
      
      h_coincidences_chi2_lowGain[ch] -> Fill( chi2 );
      
      h_coincidences[ch] -> Fill( amp );
      if( chi2 < 5. ) h_coincidences_chi2Cut[ch] -> Fill( amp );
    }
    
    
    for(auto timeMethod : timeMethods)
    {
      std::string ch1 = channels.at(0);
      std::string ch2 = channels.at(1);
      
      std::string ampCh1 = opts.GetOpt<std::string>(Form("%s.ampCh",ch1.c_str()));
      std::string ampCh2 = opts.GetOpt<std::string>(Form("%s.ampCh",ch2.c_str()));
      float ampMin1 = opts.GetOpt<float>(Form("%s.ampMin", ch1.c_str()));
      float ampMax1 = opts.GetOpt<float>(Form("%s.ampMax", ch1.c_str()));
      float ampMin2 = opts.GetOpt<float>(Form("%s.ampMin", ch2.c_str()));
      float ampMax2 = opts.GetOpt<float>(Form("%s.ampMax", ch2.c_str()));
      
      float amp1  = amp_max[channelIds[ampCh1]]  / 4096. * 1000. / peCalib_lowGain[ch1];
      float amp2  = amp_max[channelIds[ampCh2]]  / 4096. * 1000. / peCalib_lowGain[ch2];
      if( (amp1 < ampMin1) || (amp1 > ampMax1) ) continue;
      if( (amp2 < ampMin2) || (amp2 > ampMax2) ) continue;
      
      std::string timeCh1 = opts.GetOpt<std::string>(Form("%s.timeCh",ch1.c_str()));
      std::string timeCh2 = opts.GetOpt<std::string>(Form("%s.timeCh",ch2.c_str()));
      float timeMin1 = opts.GetOpt<float>(Form("%s.timeMin", ch1.c_str()));
      float timeMax1 = opts.GetOpt<float>(Form("%s.timeMax", ch1.c_str()));
      float timeMin2 = opts.GetOpt<float>(Form("%s.timeMin", ch2.c_str()));
      float timeMax2 = opts.GetOpt<float>(Form("%s.timeMax", ch2.c_str()));

      float tim1  = time[channelIds[timeCh1]+timeMethodIds[timeMethod]];
      float tim2  = time[channelIds[timeCh2]+timeMethodIds[timeMethod]];
      if( (tim1 < timeMin1) || (tim1 > timeMax1) ) continue;
      if( (tim2 < timeMin2) || (tim2 > timeMax2) ) continue;

      h_ampRatio -> Fill( amp2/amp1 );
      h_CTR_raw[timeMethod] -> Fill( tim2-tim1 );
    }
  }
  
  
  //--- find CTR ranges
  std::map<std::string,float> CTRRanges_mean;
  std::map<std::string,float> CTRRanges_sigma;
  for(auto timeMethod: timeMethods)
  {
    histo = h_CTR_raw[timeMethod];
    float* vals = new float[6];
    FindSmallestInterval(vals,histo,0.68);
    
    float mean = vals[0];
    float min = vals[4];
    float max = vals[5];
    float delta = max-min;
    float sigma = 0.5*delta;
    CTRRanges_mean[timeMethod] = mean;
    CTRRanges_sigma[timeMethod] = sigma;
    std::cout << ">>> timeMethod: " << timeMethod << "   CTR mean: " << CTRRanges_mean[timeMethod] << "   CTR sigma: " << CTRRanges_sigma[timeMethod] << std::endl;
    
    outFile -> cd();
    histo -> Write();
  }

  
  //--- draw plots
  c = new TCanvas("c_coincidences","c_coincidences",2000,2000);
  c -> Divide(2,2);
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    float xMin  = opts.GetOpt<float>(Form("%s.ampMin", ch.c_str()));
    float xMax  = opts.GetOpt<float>(Form("%s.ampMax", ch.c_str()));

    if( index < 0 ) continue;
    
    c -> cd(index);
    
    histo = h_coincidences[ch];
    histo -> SetTitle(Form(";SiPM%d amplitude (low gain) [n_{p.e.}];events",index));
    histo -> SetLineWidth(2);
    histo -> SetLineColor(kGreen-5);
    histo -> SetFillColor(kGreen-5);
    histo -> SetFillStyle(3003);
    histo -> GetXaxis() -> SetRangeUser(250.,3500.);
    histo -> Draw("hist");

    histo = h_coincidences_chi2Cut[ch];
    histo -> SetLineWidth(2);
    histo -> SetLineColor(kGreen-8);
    histo -> SetFillColor(kGreen-8);
    histo -> SetFillStyle(3001);
    histo -> Draw("hist,same");

    TLine* line_ampMin = new TLine(xMin,histo->GetMinimum(),xMin,2.*histo->GetMaximum());
    line_ampMin -> SetLineWidth(2);
    line_ampMin -> SetLineStyle(7);
    line_ampMin -> SetLineColor(kBlack);
    line_ampMin -> Draw("same");
    
    TLine* line_ampMax = new TLine(xMax,histo->GetMinimum(),xMax,2.*histo->GetMaximum());
    line_ampMax -> SetLineWidth(2);
    line_ampMax -> SetLineStyle(7);
    line_ampMax -> SetLineColor(kBlack);
    line_ampMax -> Draw("same");
    
    latexLabels[ch] -> Draw("same");
    
    histo = h_coincidences[ch];
    TF1* total = new TF1("total","gaus(0)+gaus(3)+gaus(6)",xMin,xMax);
    total -> SetParameters(1000.,0.5*(xMin+xMax),40.,100.,0.5*(xMin+xMax),80.,100.,0.5*(xMin+xMax),120.);
    total -> SetNpx(10000);
    total -> SetLineColor(kBlue);
    total -> SetLineWidth(2);
    histo -> Fit(total,"QNRS+");
    total -> Draw("same");
    
    TLatex* latex = new TLatex(0.50,0.85,Form("photopeak: %.0f p.e.",total->GetMaximumX()));
    latex -> SetNDC();
    latex -> SetTextFont(82);
    latex -> SetTextSize(0.04);
    latex -> Draw("same");

    float FWHM = GetFWHM(total,xMin,xMax);
    TLatex* latex2 = new TLatex(0.50,0.75,Form("FWHM: %.1f%%",FWHM/total->GetMaximumX()*100.));
    latex2 -> SetNDC();
    latex2 -> SetTextFont(82);
    latex2 -> SetTextSize(0.04);
    latex2 -> Draw("same");

    TLine* line_FWHM = new TLine(total->GetMaximumX()-0.5*FWHM,0.5*total->GetMaximum(),total->GetMaximumX()+0.5*FWHM,0.5*total->GetMaximum());
    line_FWHM -> SetLineColor(kBlue);
    line_FWHM -> SetLineWidth(2);
    line_FWHM -> Draw("same");
    
    gPad -> SetLogy();
    
    c -> cd(index+2);
    int timeMethodIt = 0;
    for(auto timeMethod : timeMethods)
    {
      histo = h_time[ch][timeMethod];
      histo -> SetTitle(Form(";SiPM%d time (high gain) [ns];event fraction",index));
      histo -> SetLineWidth(2);
      histo -> SetLineColor(kRed-timeMethodIt);
      histo -> SetFillColor(kRed-timeMethodIt);
      histo -> SetFillStyle(3003);
      histo -> GetXaxis() -> SetRangeUser(0.,200.);
      if( timeMethodIt == 0 ) histo -> Draw("hist");
      else                    histo -> Draw("hist,same");
      ++timeMethodIt;
    }
    latexLabels[ch] -> Draw("same");
    
    gPad -> SetLogy();
    gPad -> Update();
  }
  c -> Print(Form("%s/c_coincidences.png",plotDir.c_str()));
  c -> Print(Form("%s/c_coincidences.pdf",plotDir.c_str()));
  
  
  
  //-----------------------
  // 2nd loop over events
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 2nd loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    h4_coincidences -> GetEntry(entry);
    
    
    for(auto timeMethod : timeMethods)
    {
      std::string ch1 = channels.at(0);
      std::string ch2 = channels.at(1);
      
      std::string ampCh1 = opts.GetOpt<std::string>(Form("%s.ampCh",ch1.c_str()));
      std::string ampCh2 = opts.GetOpt<std::string>(Form("%s.ampCh",ch2.c_str()));
      float ampMin1 = opts.GetOpt<float>(Form("%s.ampMin", ch1.c_str()));
      float ampMax1 = opts.GetOpt<float>(Form("%s.ampMax", ch1.c_str()));
      float ampMin2 = opts.GetOpt<float>(Form("%s.ampMin", ch2.c_str()));
      float ampMax2 = opts.GetOpt<float>(Form("%s.ampMax", ch2.c_str()));
      
      float amp1  = amp_max[channelIds[ampCh1]]  / 4096. * 1000. / peCalib_lowGain[ch1];
      float amp2  = amp_max[channelIds[ampCh2]]  / 4096. * 1000. / peCalib_lowGain[ch2];
      if( (amp1 < ampMin1) || (amp1 > ampMax1) ) continue;
      if( (amp2 < ampMin2) || (amp2 > ampMax2) ) continue;
      
      std::string timeCh1 = opts.GetOpt<std::string>(Form("%s.timeCh",ch1.c_str()));
      std::string timeCh2 = opts.GetOpt<std::string>(Form("%s.timeCh",ch2.c_str()));
      float timeMin1 = opts.GetOpt<float>(Form("%s.timeMin", ch1.c_str()));
      float timeMax1 = opts.GetOpt<float>(Form("%s.timeMax", ch1.c_str()));
      float timeMin2 = opts.GetOpt<float>(Form("%s.timeMin", ch2.c_str()));
      float timeMax2 = opts.GetOpt<float>(Form("%s.timeMax", ch2.c_str()));

      float tim1  = time[channelIds[timeCh1]+timeMethodIds[timeMethod]];
      float tim2  = time[channelIds[timeCh2]+timeMethodIds[timeMethod]];
      if( (tim1 < timeMin1) || (tim1 > timeMax1) ) continue;
      if( (tim2 < timeMin2) || (tim2 > timeMax2) ) continue;

      if( (tim2-tim1) < CTRRanges_mean[timeMethod]-3.*CTRRanges_sigma[timeMethod] ) continue;
      if( (tim2-tim1) > CTRRanges_mean[timeMethod]+3.*CTRRanges_sigma[timeMethod] ) continue;
      p_CTR_vs_ampRatio[timeMethod] -> Fill( amp2/amp1,tim2-tim1 );
    }
  }

  
  c = new TCanvas("c_CTR_vs_ampRatio","c_CTR_vs_ampRatio",2000,2000);
  c -> Divide(3,3);

  TProfile* prof;
  std::map<std::string,TF1*> fitFunc_corrAmpRatio;

  int timeMethodIt = 0;
  for(auto timeMethod : timeMethods)
  {
    ++timeMethodIt;
    c -> cd(timeMethodIt);
    
    prof = p_CTR_vs_ampRatio[timeMethod];
    prof -> Rebin(5);
    prof -> SetTitle(Form(";amplitude_{SiPM2} / amplitude_{SiPM1};#Deltat_{1,2} [ns]"));
    prof -> SetMarkerColor(kBlack);
    prof -> SetMarkerSize(1.);
    prof -> GetXaxis() -> SetRangeUser(0.5,2.);
    prof -> GetYaxis() -> SetRangeUser(CTRRanges_mean[timeMethod]-2.*CTRRanges_sigma[timeMethod],CTRRanges_mean[timeMethod]+2.*CTRRanges_sigma[timeMethod]);
    prof -> Draw();

    fitFunc_corrAmpRatio[timeMethod] = new TF1(Form("fitFunc_corrAmpRatio_%s",timeMethod.c_str()),"pol4",0.,2.);
    prof -> Fit(fitFunc_corrAmpRatio[timeMethod],"QNRS+");
    fitFunc_corrAmpRatio[timeMethod] -> Draw("same");

    float thr1 = atof((timeMethod.substr(3,4)).c_str()) / 4096. * 1000. / peCalib_higGain[channels.at(0)];
    float thr2 = atof((timeMethod.substr(3,4)).c_str()) / 4096. * 1000. / peCalib_higGain[channels.at(1)];
    TLatex* latex = new TLatex(0.20,0.85,Form("SiPM1: thr. for time rec. = %.1f",thr1));
    latex -> SetNDC();
    latex -> SetTextFont(82);
    latex -> SetTextSize(0.04);
    latex -> Draw("same");
    TLatex* latex2 = new TLatex(0.20,0.80,Form("SiPM2: thr. for time rec. = %.1f",thr2));
    latex2 -> SetNDC();
    latex2 -> SetTextFont(82);
    latex2 -> SetTextSize(0.04);
    latex2 -> Draw("same");
    
    latexLabels[channels.at(0)] -> Draw("same");
    
    gPad -> Update();
  }
  c -> Print(Form("%s/c_timeWalkCorrection.png",plotDir.c_str()));
  c -> Print(Form("%s/c_timeWalkCorrection.pdf",plotDir.c_str()));



  //-----------------------
  // 3rd loop over events

  std::map<std::string,TH1F*> h_CTR;
  std::map<std::string,TH1F*> h_CTR_ampRatioCorr;

  std::map<std::string,std::map<std::string,TH1F*> > h_noise;
  
  for(auto timeMethod : timeMethods)
  {
    h_CTR[timeMethod] = new TH1F(Form("h_CTR_%s",timeMethod.c_str()),"",1000,CTRRanges_mean[timeMethod]-5.*CTRRanges_sigma[timeMethod],CTRRanges_mean[timeMethod]+5.*CTRRanges_sigma[timeMethod]);
    h_CTR_ampRatioCorr[timeMethod] = new TH1F(Form("h_CTR_ampRatioCorr_%s",timeMethod.c_str()),"",1000,CTRRanges_mean[timeMethod]-5.*CTRRanges_sigma[timeMethod],CTRRanges_mean[timeMethod]+5.*CTRRanges_sigma[timeMethod]);
  }
  for(auto ch : channels)
    for(auto timeMethod : timeMethods)
    {    
      h_noise[ch][timeMethod] = new TH1F(Form("h_noise_%s_%s",ch.c_str(),timeMethod.c_str()),"",10000,0.,10.);
    }
  
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 3rd loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    h4_coincidences -> GetEntry(entry);
    
    
    for(auto timeMethod : timeMethods)
    {
      std::string ch1 = channels.at(0);
      std::string ch2 = channels.at(1);
      
      std::string ampCh1 = opts.GetOpt<std::string>(Form("%s.ampCh",ch1.c_str()));
      std::string ampCh2 = opts.GetOpt<std::string>(Form("%s.ampCh",ch2.c_str()));
      float ampMin1 = opts.GetOpt<float>(Form("%s.ampMin", ch1.c_str()));
      float ampMax1 = opts.GetOpt<float>(Form("%s.ampMax", ch1.c_str()));
      float ampMin2 = opts.GetOpt<float>(Form("%s.ampMin", ch2.c_str()));
      float ampMax2 = opts.GetOpt<float>(Form("%s.ampMax", ch2.c_str()));
      
      float amp1  = amp_max[channelIds[ampCh1]]  / 4096. * 1000. / peCalib_lowGain[ch1];
      float amp2  = amp_max[channelIds[ampCh2]]  / 4096. * 1000. / peCalib_lowGain[ch2];
      if( (amp1 < ampMin1) || (amp1 > ampMax1) ) continue;
      if( (amp2 < ampMin2) || (amp2 > ampMax2) ) continue;
      
      std::string timeCh1 = opts.GetOpt<std::string>(Form("%s.timeCh",ch1.c_str()));
      std::string timeCh2 = opts.GetOpt<std::string>(Form("%s.timeCh",ch2.c_str()));
      float timeMin1 = opts.GetOpt<float>(Form("%s.timeMin", ch1.c_str()));
      float timeMax1 = opts.GetOpt<float>(Form("%s.timeMax", ch1.c_str()));
      float timeMin2 = opts.GetOpt<float>(Form("%s.timeMin", ch2.c_str()));
      float timeMax2 = opts.GetOpt<float>(Form("%s.timeMax", ch2.c_str()));

      float tim1  = time[channelIds[timeCh1]+timeMethodIds[timeMethod]];
      float tim2  = time[channelIds[timeCh2]+timeMethodIds[timeMethod]];
      if( (tim1 < timeMin1) || (tim1 > timeMax1) ) continue;
      if( (tim2 < timeMin2) || (tim2 > timeMax2) ) continue;

      h_noise[ch1][timeMethod] -> Fill( b_rms/time_slope[channelIds[timeCh1]+timeMethodIds[timeMethod]] );
      h_noise[ch2][timeMethod] -> Fill( b_rms/time_slope[channelIds[timeCh2]+timeMethodIds[timeMethod]] );
      
      h_CTR[timeMethod] -> Fill( tim2-tim1 );
      h_CTR_ampRatioCorr[timeMethod] -> Fill( (tim2-tim1) - fitFunc_corrAmpRatio[timeMethod]->Eval(amp2/amp1) + fitFunc_corrAmpRatio[timeMethod]->Eval(h_ampRatio->GetMean()) );
    }
  }


  c = new TCanvas("c_CTR","c_CTR",2000,2000);
  c -> Divide(3,3);

  TGraphErrors* g_CTR_effSigma = new TGraphErrors();
  TGraphErrors* g_CTR_gaus = new TGraphErrors();
  TGraphErrors* g_noise = new TGraphErrors();
  
  timeMethodIt = 0;
  for(auto timeMethod : timeMethods)
  {
    ++timeMethodIt;
    c -> cd(timeMethodIt);
    
    histo = h_CTR[timeMethod];
    histoCorr = h_CTR_ampRatioCorr[timeMethod];
    
    if( histo->Integral() < 100 ) continue;
    CTRResult tRes = drawCTRPlot(histoCorr,"amp. walk corr.",5,false,false,-1,";#Deltat_{1,2} [ns]",latexLabels[channels.at(0)],histo,"raw",NULL,"");
    
    float thr1 = atof((timeMethod.substr(3,4)).c_str()) / 4096. * 1000. / peCalib_higGain[channels.at(0)];
    float thr2 = atof((timeMethod.substr(3,4)).c_str()) / 4096. * 1000. / peCalib_higGain[channels.at(1)];
    TLatex* latex = new TLatex(0.20,0.55,Form("SiPM1 thr.: %.1f",thr1));
    latex -> SetNDC();
    latex -> SetTextFont(82);
    latex -> SetTextSize(0.04);
    latex -> Draw("same");
    TLatex* latex2 = new TLatex(0.20,0.50,Form("SiPM2 thr.: %.1f",thr2));
    latex2 -> SetNDC();
    latex2 -> SetTextFont(82);
    latex2 -> SetTextSize(0.04);
    latex2 -> Draw("same");
    
    gPad -> Update();

    g_CTR_effSigma -> SetPoint(g_CTR_effSigma->GetN(),0.5*(thr1+thr2),tRes.effSigma);
    g_CTR_effSigma -> SetPointError(g_CTR_effSigma->GetN()-1,0.,50.);
    
    g_CTR_gaus -> SetPoint(g_CTR_gaus->GetN(),0.5*(thr1+thr2),tRes.gausSigma);
    g_CTR_gaus -> SetPointError(g_CTR_gaus->GetN()-1,0.,tRes.gausSigmaErr);
    
    g_noise -> SetPoint(g_noise->GetN(),0.5*(thr1+thr2),1000.*0.5*(h_noise[channels.at(0)][timeMethod]->GetMean()+h_noise[channels.at(1)][timeMethod]->GetMean()));
    g_noise -> SetPointError(g_noise->GetN()-1,0.,1000.*0.5*(h_noise[channels.at(0)][timeMethod]->GetMeanError()+h_noise[channels.at(1)][timeMethod]->GetMeanError()));
  }
  c -> Print(Form("%s/c_CTR.png",plotDir.c_str()));
  c -> Print(Form("%s/c_CTR.pdf",plotDir.c_str()));
  
  
  c = new TCanvas("c_CTR_vs_thr","c_CTR_vs_thr",2000,2000);
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.1,0.,500.,2500.) );
  hPad -> SetTitle(Form(";threshold [n_{p.e.}];#sigma_{t} [ps]"));
  hPad -> Draw();
  
  g_CTR_effSigma -> SetLineColor(kRed);
  g_CTR_effSigma -> SetMarkerColor(kRed);
  g_CTR_effSigma -> SetMarkerStyle(20);
  g_CTR_effSigma -> SetMarkerSize(2);
  g_CTR_effSigma -> Draw("P,same");
  
  g_CTR_gaus -> SetLineColor(kBlue);
  g_CTR_gaus -> SetMarkerColor(kBlue);
  g_CTR_gaus -> SetMarkerStyle(20);
  g_CTR_gaus -> SetMarkerSize(2);
  g_CTR_gaus -> Draw("P,same");
  
  g_noise -> SetLineColor(kBlack);
  g_noise -> SetLineStyle(7);
  g_noise -> SetLineWidth(2);
  g_noise -> Draw("L,same");

  TLegend* legend2 = new TLegend(0.60,0.60,0.90,0.90);
  legend2 -> SetFillColor(0);
  legend2 -> SetFillStyle(1000);  
  legend2 -> SetTextFont(42);
  legend2 -> AddEntry(g_CTR_effSigma,"eff. sigma","P");
  legend2 -> AddEntry(g_CTR_gaus,    "gaus. fit", "P");
  legend2 -> AddEntry(g_noise,       "noise contrib. (#sigma_{V} / dV/dt)","L");
  legend2 -> Draw("same");
  
  gPad -> SetLogx();
  gPad -> Update();
  
  c -> Print(Form("%s/c_CTR_vs_thr.png",plotDir.c_str()));
  c -> Print(Form("%s/c_CTR_vs_thr.pdf",plotDir.c_str()));
  
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
