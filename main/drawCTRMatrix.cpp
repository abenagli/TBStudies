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
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"



int main(int argc, char** argv)
{
  setTDRStyle();

  
  if( argc < 2 )
  {
    std::cout << ">>> drawCTRMatrix::usage:   " << argv[0] << " configFile.cfg" << std::endl;
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
  
  
  
  //--- open input files
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  int runMin = opts.GetOpt<int>("Input.runMin");
  int runMax = opts.GetOpt<int>("Input.runMax");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");
  
  TChain* h4 = new TChain("h4","h4");
  for(int run = runMin; run <= runMax; ++run)
  {
    std::string fileName = Form("%s/%d.root",inputDir.c_str(),run);
    std::cout << ">>> Adding flle " << fileName << std::endl;
    h4 -> Add(fileName.c_str());
  }
  

  
  //--- define tree branch addresses
  std::map<std::string,int> channelIds;
  std::map<std::string,int> timeMethodIds;
  
  for(auto ch : channels)
  {
    std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
    if( ampCh != "NULL" ) h4 -> SetBranchAddress(ampCh.c_str(), &channelIds[ampCh.c_str()]);
    
    std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
    if( timeCh != "NULL" ) h4 -> SetBranchAddress(timeCh.c_str(), &channelIds[timeCh.c_str()]);
    
    std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));    
    for(auto timeMethod: timeMethods)
      if( timeMethod != "NULL" ) h4 -> SetBranchAddress(timeMethod.c_str(), &timeMethodIds[timeMethod.c_str()]);
  }
  
  float* time    = new float[1000];
  float* amp_max = new float[1000];
  h4 -> SetBranchAddress("time",   time);
  h4 -> SetBranchAddress("amp_max",amp_max);
  
  int nEntries = h4->GetEntries();
  std::cout << ">>> Events read: " << nEntries << std::endl;
  
  
  
  //------------------
  // define histograms
  
  TFile* outFile = TFile::Open(Form("%s/drawCTRMatrix_%s.root",plotDir.c_str(),conf.c_str()),"RECREATE");
  outFile -> cd();

  TProfile2D* p2_amp = new TProfile2D("p2_amp","",4,-0.5,3.5,4,-0.5,3.5);
  
  std::map<std::string,TH1F*> h_amp;
  std::map<std::string,TH1F*> h_amp_noXTalk;
  std::map<std::string,TH1F*> h_amp_cut;
  std::map<std::string,TH1F*> h_amp_others;
  std::map<std::string,TH1F*> h_rt;
  std::map<std::string,TH1F*> h_rt_cut;
  std::map<std::string,TH1F*> h_time;
  std::map<std::string,TH1F*> h_duration;
  std::map<std::string,TH2F*> h2_duration_vs_amp;
  std::map<std::string,TProfile2D*> p2_lightSharing;
  std::map<std::string,TH1F*> h_CTR_raw;
  std::map<std::string,TProfile*> p_time_vs_amp;
  std::map<std::string,TProfile*> p_time_vs_rt;
  TH2F* h2_CTR_ampCorr_effSigma = new TH2F("h2_CTR_ampCorr_effSigma","",4,-0.5,3.5,4,-0.5,3.5);
  TH2F* h2_CTR_ampCorr_gausFit  = new TH2F("h2_CTR_ampCorr_gausFit", "",4,-0.5,3.5,4,-0.5,3.5);
  TH2F* h2_CTR_RTCorr_effSigma = new TH2F("h2_CTR_RTCorr_effSigma","",4,-0.5,3.5,4,-0.5,3.5);
  TH2F* h2_CTR_RTCorr_gausFit  = new TH2F("h2_CTR_RTCorr_gausFit", "",4,-0.5,3.5,4,-0.5,3.5);
  
  for(auto ch : channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
    std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
    
    h_amp[ch]         = new TH1F(Form("h_amp_%s",        ch.c_str()),"",4000,0.,1.);
    h_amp_noXTalk[ch] = new TH1F(Form("h_amp_noXTalk_%s",ch.c_str()),"",4000,0.,1.);
    h_amp_cut[ch]     = new TH1F(Form("h_amp_cut_%s",    ch.c_str()),"",4000,0.,1.);

    h_rt[ch]         = new TH1F(Form("h_rt_%s",    ch.c_str()),"",2500,0.,100.);
    h_rt_cut[ch]     = new TH1F(Form("h_rt_cut_%s",ch.c_str()),"",2500,0.,100.);
    
    h_time[ch] = new TH1F(Form("h_time_%s",ch.c_str()),"",100,0.,200.);
    h_duration[ch] = new TH1F(Form("h_duration_%s",ch.c_str()),"",100,0.,200.);
    
    h2_duration_vs_amp[ch] = new TH2F(Form("h_duration_vs_amp_%s",ch.c_str()),"",500,0.,1.,200,0.,200.);
    
    if( index < 0 ) continue;

    for(auto ch2 : channels)
    {
      std::string otherLabel = ch + "_" + ch2;
      h_amp_others[otherLabel] = new TH1F(Form("h_amp_%s",otherLabel.c_str()),"",4000,0.,1.);      
    }

    p2_lightSharing[ch] = new TProfile2D(Form("p2_lightSharing_%s",ch.c_str()),"",4,-0.5,3.5,4,-0.5,3.5);
    
    p_time_vs_amp[ch] = new TProfile(Form("p_time_vs_amp_%s",ch.c_str()),"",200,0.,1.);
    p_time_vs_rt[ch] = new TProfile(Form("p_time_vs_rt_%s",ch.c_str()),"",500,0.,100.);
    h_CTR_raw[ch] = new TH1F(Form("h_CTR_raw_%s",ch.c_str()),"",10000,0.,50.);
  }
  
  
  
  //-----------------------
  // 1st loop over events
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    h4 -> GetEntry(entry);
    
    
    for(auto ch: channels)
    {
      // fill amplitude and time plots
      int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
      std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
      std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
      std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));

      if( ampCh == "NULL" ) continue;
      float amp = amp_max[channelIds[ampCh]] / 4096.;
      h_amp[ch] -> Fill( amp );
      
      float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
      float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));

//       for(auto ch2: channels)
//       {
//         if( index < 0 ) continue;
//         std::string otherLabel = ch + "_" + ch2;
        
//         std::string ampChTemp  = opts.GetOpt<std::string>(Form("%s.ampCh", ch2.c_str()));
//         float ampMinTemp = opts.GetOpt<float>(Form("%s.ampMin", ch2.c_str()));
//         float ampTemp = amp_max[channelIds[ampChTemp]] / 4096.;
        
//         if( ampTemp > ampMinTemp )
//           h_amp_others[otherLabel] -> Fill( amp );
//       }

      bool isXTalkOk = true;
      std::vector<std::string> vetoCh = opts.GetOpt<std::vector<std::string> >(Form("%s.vetoCh",ch.c_str()));
      for(auto ch2: vetoCh)
      {
        if( ch2 == "NULL" ) continue;
        std::string ampChTemp  = opts.GetOpt<std::string>(Form("%s.ampCh", ch2.c_str()));
        float ampVeto = opts.GetOpt<float>(Form("%s.ampVeto",ch2.c_str()));
        float ampTemp = amp_max[channelIds[ampChTemp]] / 4096.;
        if( ampTemp > ampVeto)
        {
          isXTalkOk = false;
          // break;
        }
      }
      if( !isXTalkOk ) continue;
      
      h_amp_noXTalk[ch] -> Fill( amp );
      
      if( amp < ampMin || amp > ampMax ) continue;
      if( isnan(amp) ) continue;
      
      h_amp_cut[ch] -> Fill( amp );
      
      if( index > 0 )
      {
        p2_amp -> Fill( (index-1)%4,3-floor((index-1)/4.),amp );
        
        for(auto ch2 : channels)
        {
          int indexTemp  = opts.GetOpt<int>(Form("%s.index", ch2.c_str()));
          std::string ampChTemp  = opts.GetOpt<std::string>(Form("%s.ampCh", ch2.c_str()));
          float ampTemp = amp_max[channelIds[ampChTemp]] / 4096.;
          p2_lightSharing[ch] -> Fill( (indexTemp-1)%4,3-floor((indexTemp-1)/4.),ampTemp );
        }
      }
      
      if( timeCh == "NULL" ) continue;
      float tim  = time[channelIds[timeCh]+timeMethodIds[timeMethods.at(0)]];
      float tim2 = timeMethods.at(1) != "NULL" ? time[channelIds[timeCh]+timeMethodIds[timeMethods.at(1)]] : 0.;
      h_time[ch] -> Fill( tim );
      h_duration[ch] -> Fill( tim2-tim );      
      h2_duration_vs_amp[ch] -> Fill( amp,(tim2-tim) );
      
      if( index < 0 ) continue;

      float tim50 = time[channelIds[ampCh]+timeMethodIds[timeMethods.at(2)]];
      float tim1000 = time[channelIds[ampCh]+timeMethodIds[timeMethods.at(3)]];
      h_rt[ch] -> Fill( tim1000-tim50 );
      
      // get reference channel
      std::string refCh  = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
      std::string ampChRef  = opts.GetOpt<std::string>(Form("%s.ampCh", refCh.c_str()));
      std::string timeChRef = opts.GetOpt<std::string>(Form("%s.timeCh",refCh.c_str()));
      std::vector<std::string> timeMethodsRef = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",refCh.c_str()));
      float ampMinRef = opts.GetOpt<float>(Form("%s.ampMin",refCh.c_str())); 
      float ampMaxRef = opts.GetOpt<float>(Form("%s.ampMax",refCh.c_str()));
      
      float ampRef = amp_max[channelIds[ampChRef]] / 4096.;
      float timRef = time[channelIds[timeChRef]+timeMethodIds[timeMethodsRef.at(0)]];
      
      if( ampRef < ampMinRef || ampRef > ampMaxRef ) continue;
      if( isnan(ampRef) ) continue;
      
      if( isinf(timRef) ) continue;
      if( isnan(timRef) ) continue;
      if( timRef < 0. || timRef > 40. ) continue;
      
      if( isinf(tim) ) continue;
      if( isnan(tim) ) continue;
      if( tim < 0. || tim > 200. ) continue;

      h_rt_cut[ch] -> Fill( tim1000-tim50 );
      
      h_CTR_raw[ch] -> Fill( tim-timRef );
    }
  }
  std::cout << "\n>>> end 1st loop" << std::endl;
  
  
  
  //--- draw plots
  TCanvas* c;
  TH1F* histo;
  TH1F* histoCorr;
  TH2F* histo2;
  
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
  
  
  c = new TCanvas("c_amp_SiPM","c_amp_SiPM",2000,2000);
  c -> Divide(4,4);
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string ampCh = opts.GetOpt<std::string>(Form("%s.ampCh",ch.c_str()));
    float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
    float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));
    if( index < 0 ) continue;
    
    c -> cd(index);
    gPad -> SetLogy();
    histo = h_amp[ch];
    histo -> SetTitle(Form(";SiPM%d amplitude [mV];entries",index));
    histo -> SetLineColor(kRed);
    histo -> SetFillColor(kRed);
    histo -> SetFillStyle(3003);
    histo -> GetXaxis() -> SetRangeUser(0.,1.1*ampMax);
    histo -> Draw();

    histo = h_amp_noXTalk[ch];
    histo -> SetLineColor(kOrange);
    histo -> SetFillColor(kOrange);
    histo -> SetFillStyle(1);
    histo -> Draw("same");
    
    TLine* line_min = new TLine(ampMin,histo->GetMinimum(),ampMin,histo->GetMaximum());
    line_min -> SetLineStyle(2);
    line_min -> Draw("same");
    TLine* line_max = new TLine(ampMax,histo->GetMinimum(),ampMax,histo->GetMaximum());
    line_max -> SetLineStyle(2);
    line_max -> Draw("same");

    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();
  }
  c -> Print(Form("%s/c_matrix_amp_SiPM.png",plotDir.c_str()));
  c -> Print(Form("%s/c_matrix_amp_SiPM.pdf",plotDir.c_str()));
  
  c = new TCanvas("c_time_SiPM","c_time_SiPM",2000,2000);
  c -> Divide(4,4);
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
    if( index < 0 ) continue;
    
    c -> cd(index);
    gPad -> SetLogy();
    histo = h_time[ch];
    histo -> SetTitle(Form(";SiPM%d time [ns];entries",index));
    histo -> SetLineColor(kRed);
    histo -> SetFillColor(kRed);
    histo -> SetFillStyle(3003);
    histo -> GetXaxis() -> SetRangeUser(0.,200.);
    histo -> Draw();

    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();
  }
  c -> Print(Form("%s/c_matrix_time_SiPM.png",plotDir.c_str()));
  c -> Print(Form("%s/c_matrix_time_SiPM.pdf",plotDir.c_str()));
  
  c = new TCanvas("c_duration_SiPM","c_duration_SiPM",2000,2000);
  c -> Divide(4,4);
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
    if( index < 0 ) continue;
    
    c -> cd(index);
    gPad -> SetLogy();
    histo = h_duration[ch];
    histo -> SetTitle(Form(";SiPM%d duration [ns];entries",index));
    histo -> SetLineColor(kRed);
    histo -> SetFillColor(kRed);
    histo -> SetFillStyle(3003);
    histo -> GetXaxis() -> SetRangeUser(0.,200.);
    histo -> Draw();

    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();
  }
  c -> Print(Form("%s/c_matrix_duration_SiPM.png",plotDir.c_str()));
  c -> Print(Form("%s/c_matrix_duration_SiPM.pdf",plotDir.c_str()));
  
  c = new TCanvas("c_duration_vs_amp_SiPM","c_duration_vs_amp_SiPM",2000,2000);
  c -> Divide(4,4);
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
    if( index < 0 ) continue;
    
    c -> cd(index);
    histo2 = h2_duration_vs_amp[ch];
    histo2 -> SetTitle(Form(";SiPM%d amp [mV];SiPM%d duration [ns];entries",index,index));
    histo2 -> GetXaxis() -> SetRangeUser(0.,200.);
    histo2 -> Draw("COLZ");

    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();
  }
  c -> Print(Form("%s/c_matrix_duration_vs_amp_SiPM.png",plotDir.c_str()));
  c -> Print(Form("%s/c_matrix_duration_vs_amp_SiPM.pdf",plotDir.c_str()));
  
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string ampCh = opts.GetOpt<std::string>(Form("%s.ampCh",ch.c_str()));
    std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
    float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
    float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));
    
    if( index > 0 ) continue;
    
    c = new TCanvas(Form("c_amp_%s",ch.c_str()),Form("c_amp_%s",ch.c_str()),2000,2000);
    gPad -> SetLogy();
    histo = h_amp[ch];
    histo -> SetTitle(Form(";%s amplitude [mV];entries",ch.c_str()));
    histo -> SetLineColor(kBlue);
    histo -> SetFillColor(kBlue);
    histo -> SetFillStyle(3003);
    histo -> Draw();
    
    TLine* line_min = new TLine(ampMin,histo->GetMinimum(),ampMin,histo->GetMaximum());
    line_min -> SetLineStyle(2);
    line_min -> Draw("same");
    TLine* line_max = new TLine(ampMax,histo->GetMinimum(),ampMax,histo->GetMaximum());
    line_max -> SetLineStyle(2);
    line_max -> Draw("same");
    
    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();
    c -> Print(Form("%s/c_amp_%s.png",plotDir.c_str(),ch.c_str()));
    c -> Print(Form("%s/c_amp_%s.pdf",plotDir.c_str(),ch.c_str()));
    
    c = new TCanvas(Form("c_time_%s",ch.c_str()),Form("c_time_%s",ch.c_str()),2000,2000);
    gPad -> SetLogy();
    histo = h_time[ch];
    histo -> SetTitle(Form(";%s time [ns];entries",ch.c_str()));
    histo -> SetLineColor(kBlue);
    histo -> SetFillColor(kBlue);
    histo -> SetFillStyle(3003);
    histo -> Draw();

    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();
    c -> Print(Form("%s/c_time_%s.png",plotDir.c_str(),ch.c_str()));
    c -> Print(Form("%s/c_time_%s.pdf",plotDir.c_str(),ch.c_str()));
  }
  
  c = new TCanvas("c_avgAmp","c_avgAmp",2000,2000);
  c -> SetRightMargin(0.25);
  c -> SetBottomMargin(0.25);
  p2_amp -> GetZaxis() -> SetTitleOffset(1.50);
  p2_amp -> SetTitle(";;;average amplitude [mV]");
  p2_amp -> Draw("COLZ");
  c -> Print(Form("%s/c_matrix_avgAmp_SiPM.png",plotDir.c_str()));
  c -> Print(Form("%s/c_matrix_avgAmp_SiPM.pdf",plotDir.c_str()));
  
  
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string ampCh = opts.GetOpt<std::string>(Form("%s.ampCh",ch.c_str()));
    if( index < 0 ) continue;
    
    c = new TCanvas(Form("c_lightSharing_%s",ch.c_str()),Form("c_lightSharing_%s",ch.c_str()),2000,2000);
    c -> SetRightMargin(0.25);
    c -> SetBottomMargin(0.25);
    p2_lightSharing[ch] -> Scale(1./h_amp_cut[ch]->GetMean());
    p2_lightSharing[ch] -> GetZaxis() -> SetTitleOffset(1.50);
    p2_lightSharing[ch] -> SetTitle(";;;norm. average amplitude");
    p2_lightSharing[ch] -> Draw("COLZ,text");
    
    latexLabels[ch] -> Draw("same");
    
    c -> Print(Form("%s/c_matrix_lightSharing_%s.png",plotDir.c_str(),ch.c_str()));
    c -> Print(Form("%s/c_matrix_lightSharing_%s.pdf",plotDir.c_str(),ch.c_str()));
  }
  
  
  c = new TCanvas("c_rt_SiPM","c_rt_SiPM",2000,2000);
  c -> Divide(4,4);
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    c -> cd(index);
    gPad -> SetLogy();
    histo = h_rt[ch];
    histo -> SetTitle(Form(";SiPM%d rise time [ns];entries",index));
    histo -> SetLineColor(kRed);
    histo -> SetFillColor(kRed);
    histo -> SetFillStyle(3003);
    histo -> Draw();
    
    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();
  }
  c -> Print(Form("%s/c_matrix_rt_SiPM.png",plotDir.c_str()));
  c -> Print(Form("%s/c_matrix_rt_SiPM.pdf",plotDir.c_str()));
  
  
  //--- find CTR ranges
  std::map<std::string,float> CTRRanges_mean;
  std::map<std::string,float> CTRRanges_sigma;
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
    if( index < 0 ) continue;
    
    histo = h_CTR_raw[ch];
    float* vals = new float[6];
    FindSmallestInterval(vals,histo,0.68);
    
    float mean = vals[0];
    float min = vals[4];
    float max = vals[5];
    float delta = max-min;
    float sigma = 0.5*delta;
    CTRRanges_mean[ch] = mean;
    CTRRanges_sigma[ch] = sigma;
    std::cout << ">>> ch: " << index << "   CTR mean: " << CTRRanges_mean[ch] << "   CTR sigma: " << CTRRanges_sigma[ch] << std::endl;
    
    outFile -> cd();
    histo -> Write();
  }
  
  
  
  //--- 2nd loop over events
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 2nd loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    h4 -> GetEntry(entry);
    
    for(auto ch: channels)
    {
      // fill amplitude and time plots
      int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
      if( index < 0 ) continue;
      
      std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
      std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
      std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
      float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
      float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));

      if( ampCh == "NULL" ) continue;
      float amp = amp_max[channelIds[ampCh]] / 4096.;
      
      bool isXTalkOk = true;
      std::vector<std::string> vetoCh = opts.GetOpt<std::vector<std::string> >(Form("%s.vetoCh",ch.c_str()));
      for(auto ch2: vetoCh)
      {
        if( ch2 == "NULL" ) continue;
        std::string ampChTemp  = opts.GetOpt<std::string>(Form("%s.ampCh", ch2.c_str()));
        float ampVeto = opts.GetOpt<float>(Form("%s.ampVeto",ch2.c_str()));
        float ampTemp = amp_max[channelIds[ampChTemp]] / 4096.;
        if( ampTemp > ampVeto)
        {
          isXTalkOk = false;
          // break;
        }
      }
      if( !isXTalkOk ) continue;
      
      if( timeCh == "NULL" ) continue;
      float tim = time[channelIds[timeCh]+timeMethodIds[timeMethods.at(0)]];
      float tim50 = time[channelIds[ampCh]+timeMethodIds[timeMethods.at(2)]];
      float tim1000 = time[channelIds[ampCh]+timeMethodIds[timeMethods.at(3)]];
      
      if( amp < ampMin || amp > ampMax ) continue;
      if( isnan(amp) ) continue;
      if( isinf(tim) ) continue;
      if( isnan(tim) ) continue;
      if( tim < 0. || tim > 200. ) continue;  
      
      std::string refCh  = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
      std::string ampChRef  = opts.GetOpt<std::string>(Form("%s.ampCh", refCh.c_str()));
      std::string timeChRef = opts.GetOpt<std::string>(Form("%s.timeCh",refCh.c_str()));
      std::vector<std::string> timeMethodsRef = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",refCh.c_str()));
      float ampMinRef = opts.GetOpt<float>(Form("%s.ampMin",refCh.c_str())); 
      float ampMaxRef = opts.GetOpt<float>(Form("%s.ampMax",refCh.c_str()));
      
      float ampRef = amp_max[channelIds[ampChRef]] / 4096.;
      float timRef = time[channelIds[timeChRef]+timeMethodIds[timeMethodsRef.at(0)]];
      
      if( ampRef < ampMinRef || ampRef > ampMaxRef ) continue;
      if( isnan(ampRef) ) continue;
      if( isinf(timRef) ) continue;
      if( isnan(timRef) ) continue;
      if( timRef < 0. || timRef > 40. ) continue;
      
      if( (tim-timRef) < (CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch]) ) continue;
      if( (tim-timRef) > (CTRRanges_mean[ch]+3.*CTRRanges_sigma[ch]) ) continue;
      
      p_time_vs_amp[ch] -> Fill( amp,tim-timRef );
      p_time_vs_rt[ch] -> Fill( tim1000-tim50,tim-timRef );
    }
  }
  std::cout << "\n>>> end 2nd loop" << std::endl;
  
  
  
  //--- draw plots
  std::map<std::string,TF1*> fitFunc_corrAmp;
  
  c = new TCanvas("c_time_vs_amp_SiPM","c_time_vs_amp_SiPM",2000,2000);
  c -> Divide(4,4);
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;

    std::string fitFunc = opts.GetOpt<std::string>(Form("%s.fitFunc", ch.c_str()));
    
    c -> cd(index);
    TProfile* prof = p_time_vs_amp[ch];
    if( prof->GetEntries() < 30 ) continue;
    
    prof -> SetTitle(Form(";SiPM%d amplitude [mV];t_{SiPM} - t_{MCP} [ns]",index));
    prof -> SetMarkerColor(kBlack);
    prof -> SetMarkerSize(0.3);
    prof -> GetYaxis() -> SetRangeUser(CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+3.*CTRRanges_sigma[ch]);
    prof -> Draw();
    
    fitFunc_corrAmp[ch] = new TF1(Form("fitFunc_corrAmp_%s",ch.c_str()),fitFunc.c_str(),0.,1000.);
    prof -> Fit(fitFunc_corrAmp[ch],"QNRS+");
    
    // fitFunc_corrAmp[ch] = new TF1(Form("fitFunc_corrAmp_%s",ch.c_str()),"pol6",0.,1000.);
    // if( ch == "CH01" || ch == "CH02" || ch == "CH03" || ch == "CH04" || ch == "CH15" || ch == "CH16" )
    //   fitFunc_corrAmp[ch] -> SetParameters(1.08567,-11.4711,48.2153,-124.762,207.922,-203.278,84.6884);
    // else
    //   fitFunc_corrAmp[ch] -> SetParameters(-9.41201,125.669,-643.62,1689.11,-2436.31,1841.14,-570.996);
    
    fitFunc_corrAmp[ch] -> Draw("same");
    
    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();
    
    // outFile -> cd();
    // prof -> Write();
    // fitFunc_corrAmp[ii] -> Write();
  }
  c -> Print(Form("%s/c_matrix_time_vs_amp.png",plotDir.c_str()));
  c -> Print(Form("%s/c_matrix_time_vs_amp.pdf",plotDir.c_str()));
  
  
  
  std::map<std::string,TF1*> fitFunc_corrRT;
  
  c = new TCanvas("c_time_vs_rt_SiPM","c_time_vs_rt_SiPM",2000,2000);
  c -> Divide(4,4);
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;

    std::string fitFuncRT = opts.GetOpt<std::string>(Form("%s.fitFuncRT", ch.c_str()));
    
    c -> cd(index);
    TProfile* prof = p_time_vs_rt[ch];
    if( prof->GetEntries() < 30 ) continue;
    
    prof -> SetTitle(Form(";SiPM%d rise time [ns];t_{SiPM} - t_{MCP} [ns]",index));
    prof -> SetMarkerColor(kBlack);
    prof -> SetMarkerSize(0.3);
    prof -> GetYaxis() -> SetRangeUser(CTRRanges_mean[ch]-7.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+5.*CTRRanges_sigma[ch]);
    prof -> Draw();
    
    fitFunc_corrRT[ch] = new TF1(Form("fitFunc_corrRT_%s",ch.c_str()),fitFuncRT.c_str(),0.,1000.);
    prof -> Fit(fitFunc_corrRT[ch],"QNRS+");
    fitFunc_corrRT[ch] -> Draw("same");

    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();
    
    // outFile -> cd();
    // prof -> Write();
    // fitFunc_corrAmp[ii] -> Write();
  }
  c -> Print(Form("%s/c_matrix_time_vs_rt.png",plotDir.c_str()));
  c -> Print(Form("%s/c_matrix_time_vs_rt.pdf",plotDir.c_str()));  
  

  //--- 3rd loop over events
  std::map<std::string,TH1F*> h_CTR;
  std::map<std::string,TH1F*> h_CTR_ampCorr;
  std::map<std::string,TH1F*> h_CTR_RTCorr;
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    h_CTR[ch]         = new TH1F(Form("h_CTR_%s",ch.c_str()),        "",1000,CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+5.*CTRRanges_sigma[ch]);
    h_CTR_ampCorr[ch] = new TH1F(Form("h_CTR_ampCorr_%s",ch.c_str()),"",1000,CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+5.*CTRRanges_sigma[ch]);
    h_CTR_RTCorr[ch]  = new TH1F(Form("h_CTR_RTCorr_%s",ch.c_str()), "",1000,CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+5.*CTRRanges_sigma[ch]);
  }
  
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 3rd loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    h4 -> GetEntry(entry);
    
    for(auto ch: channels)
    {
      // fill amplitude and time plots
      int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
      if( index < 0 ) continue;
      
      std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
      std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
      std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
      float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
      float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));
      
      if( ampCh == "NULL" ) continue;
      float amp = amp_max[channelIds[ampCh]] / 4096.;
      
      bool isXTalkOk = true;
      std::vector<std::string> vetoCh = opts.GetOpt<std::vector<std::string> >(Form("%s.vetoCh",ch.c_str()));
      for(auto ch2: vetoCh)
      {
        if( ch2 == "NULL" ) continue;
        std::string ampChTemp  = opts.GetOpt<std::string>(Form("%s.ampCh", ch2.c_str()));
        float ampVeto = opts.GetOpt<float>(Form("%s.ampVeto",ch2.c_str()));
        float ampTemp = amp_max[channelIds[ampChTemp]] / 4096.;
        if( ampTemp > ampVeto)
        {
          isXTalkOk = false;
          // break;
        }
      }
      if( !isXTalkOk ) continue;
      
      if( timeCh == "NULL" ) continue;
      float tim = time[channelIds[timeCh]+timeMethodIds[timeMethods.at(0)]];
      float tim50 = time[channelIds[ampCh]+timeMethodIds[timeMethods.at(2)]];
      float tim1000 = time[channelIds[ampCh]+timeMethodIds[timeMethods.at(3)]];
      
      if( amp < ampMin || amp > ampMax ) continue;
      if( isnan(amp) ) continue;
      if( isinf(tim) ) continue;
      if( isnan(tim) ) continue;
      if( tim < 0. || tim > 200. ) continue;  
      
      std::string refCh  = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
      std::string ampChRef  = opts.GetOpt<std::string>(Form("%s.ampCh", refCh.c_str()));
      std::string timeChRef = opts.GetOpt<std::string>(Form("%s.timeCh",refCh.c_str()));
      std::vector<std::string> timeMethodsRef = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",refCh.c_str()));
      float ampMinRef = opts.GetOpt<float>(Form("%s.ampMin",refCh.c_str())); 
      float ampMaxRef = opts.GetOpt<float>(Form("%s.ampMax",refCh.c_str()));
      
      float ampRef = amp_max[channelIds[ampChRef]] / 4096.;
      float timRef = time[channelIds[timeChRef]+timeMethodIds[timeMethodsRef.at(0)]];
      
      if( ampRef < ampMinRef || ampRef > ampMaxRef ) continue;
      if( isnan(ampRef) ) continue;
      if( isinf(timRef) ) continue;
      if( isnan(timRef) ) continue;
      if( timRef < 0. || timRef > 40. ) continue;

      if( (tim-timRef) < (CTRRanges_mean[ch]-3.*CTRRanges_sigma[ch]) ) continue;
      if( (tim-timRef) > (CTRRanges_mean[ch]+3.*CTRRanges_sigma[ch]) ) continue;
      
      float CTR = tim - timRef;
      float CTR_ampCorr = CTR - fitFunc_corrAmp[ch]->Eval(amp) + fitFunc_corrAmp[ch]->Eval( h_amp_cut[ch]->GetMean() );
      float CTR_RTCorr = CTR - fitFunc_corrRT[ch]->Eval(tim1000-tim50) + fitFunc_corrRT[ch]->Eval( h_rt_cut[ch]->GetMean() );
      h_CTR[ch] -> Fill( CTR );
      h_CTR_ampCorr[ch] -> Fill( CTR_ampCorr );
      h_CTR_RTCorr[ch] -> Fill( CTR_RTCorr );
    }
  }
  std::cout << "\n>>> end 3rd loop" << std::endl;
  
  
  
  //--- draw plots
  int rebin = opts.GetOpt<int>("Input.rebin");
  TLatex* latexLabel12 = new TLatex(0.13,0.97,Form("%s",""));
  
  
  c = new TCanvas("c_CTR","c_CTR",2000,2000);
  c -> Divide(4,4);
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    c -> cd(index);
    histo = h_CTR[ch];
    histoCorr = h_CTR_ampCorr[ch];
    
    if( histo->Integral() < 100 ) continue;
    
    CTRResult tRes = drawCTRPlot(histoCorr,"amp. walk corr.",rebin,false,true,0.020,";t_{SiPM} - t_{MCP} [ns]",latexLabel12,histo,"raw",NULL,"");
    
    latexLabels[ch] -> Draw("same");
    
    h2_CTR_ampCorr_effSigma -> Fill( (index-1)%4,3-floor((index-1)/4.),tRes.effSigma );
    h2_CTR_ampCorr_gausFit  -> Fill( (index-1)%4,3-floor((index-1)/4.),tRes.gausSigma );
    
    gPad -> Update();
  }
  c -> Print(Form("%s/c_matrix_tRes.png",plotDir.c_str()));
  c -> Print(Form("%s/c_matrix_tRes.pdf",plotDir.c_str()));
  
  
  c = new TCanvas("c_CTR_RTCorr","c_CTR_RTCorr",2000,2000);
  c -> Divide(4,4);
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    c -> cd(index);
    histo = h_CTR[ch];
    histoCorr = h_CTR_RTCorr[ch];

    if( histo->Integral() < 100 ) continue;
    
    CTRResult tRes = drawCTRPlot(histoCorr,"rise time corr.",rebin,false,true,0.020,";t_{SiPM} - t_{MCP} [ns]",latexLabel12,histo,"raw",NULL,"");
    
    latexLabels[ch] -> Draw("same");
    
    h2_CTR_RTCorr_effSigma -> Fill( (index-1)%4,3-floor((index-1)/4.),tRes.effSigma );
    h2_CTR_RTCorr_gausFit  -> Fill( (index-1)%4,3-floor((index-1)/4.),tRes.gausSigma );
    
    gPad -> Update();
  }
  c -> Print(Form("%s/c_matrix_tRes_RTCorr.png",plotDir.c_str()));
  c -> Print(Form("%s/c_matrix_tRes_RTCorr.pdf",plotDir.c_str()));
  

  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
