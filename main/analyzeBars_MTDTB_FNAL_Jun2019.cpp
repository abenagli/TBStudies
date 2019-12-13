#include "interface/FitUtils.h"
#include "interface/TreeUtils.h"
#include "interface/SetTDRStyle.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRandom3.h"






float FindXMaximum(TH1F* histo, const float& xMin, const float& xMax)
{
  float max = -999999999.;
  int binMax = -1;
  for(int bin = 1; bin <= histo->GetNbinsX(); ++bin)
  {
    if( histo->GetBinCenter(bin) < xMin ) continue;
    if( histo->GetBinCenter(bin) > xMax ) continue;
    if( histo->GetBinContent(bin) > max ) { max = histo->GetBinContent(bin); binMax = bin; };
  }
  return histo->GetBinCenter(binMax);
}

struct EventSingle
{
  std::string stepLabel;
  std::string ch1;
  std::string ch2;
  std::string label1;
  std::string label2;
  std::string label12;
  float tot1;
  float tot2;
  long long time1;
  long long time2;
};

struct EventCoinc
{
  std::string stepLabel;
  std::string ch1;
  std::string ch2;
  std::string ch3;
  std::string ch4;
  std::string label1;
  std::string label2;
  std::string label3;
  std::string label4;
  std::string label1234;
  float tot1;
  float tot2;
  float tot3;
  float tot4;
  long long time1;
  long long time2;
  long long time3;
  long long time4;
};






int main(int argc, char** argv)
{
  setTDRStyle();
  
  
  if( argc < 2 )
  {
    std::cout << ">>> analyzeBars::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  
  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  

  
  //--- get parameters
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("cp /afs/cern.ch/user/a/abenagli/public/index.php %s",plotDir.c_str()));
  
  
  //--- open files and make the tree chain
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  std::string runs = opts.GetOpt<std::string>("Input.runs");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");
  TChain* tree = new TChain("data","data");
  
  std::stringstream ss(runs);
  std::string token;
  while( std::getline(ss,token,',') )
    {
      std::stringstream ss2(token);
      std::string token2;
      int runMin = -1;
      int runMax = -1;
      while( std::getline(ss2,token2,'-') )
	{
	  if( runMin != -1 && runMax == -1 ) runMax = atoi(token2.c_str());
	  if( runMin == -1 ) runMin = atoi(token2.c_str());
	}
      for(int run = runMin; run <= runMax; ++run)
	{
	  //std::string fileName = Form("%s/*%04d*.root",inputDir.c_str(),run);
	  std::string fileName = Form("%s/*%d*.root",inputDir.c_str(),run);
	  std::cout << ">>> Adding flle " << fileName << std::endl;
	  tree -> Add(fileName.c_str());
	}
    }
  
  
  //--- define channels
  std::vector<std::string> channels = opts.GetOpt<std::vector<std::string> >("Channels.channels");
  
  int pairsMode = opts.GetOpt<int>("Channels.pairsMode");
  std::vector<std::string> pairs = opts.GetOpt<std::vector<std::string> >("Channels.pairs");
  std::vector<std::pair<std::string,std::string> > pairsVec;
  for(unsigned int ii = 0; ii < pairs.size()/2; ++ii)
  {
    pairsVec.push_back(std::make_pair(pairs.at(0+ii*2),pairs.at(1+ii*2)));
  }
  
  std::vector<std::string> bars = opts.GetOpt<std::vector<std::string> >("Channels.bars");
  std::vector<std::pair<std::pair<std::string,std::string>,std::pair<std::string,std::string> > > barsVec;
  for(unsigned int ii = 0; ii < bars.size()/4; ++ii)
  {
    barsVec.push_back( std::make_pair(std::make_pair(bars.at(0+ii*4),bars.at(1+ii*4)),std::make_pair(bars.at(2+ii*4),bars.at(3+ii*4))) );
  }
  
  
  //--- get cuts per bar / Vov
  std::map<std::string,std::map<float,float> > cut_totAcc;
  for(auto ch :  channels)
  {
    std::vector<float> Vovs    = opts.GetOpt<std::vector<float> >(Form("%s.Vovs",ch.c_str()));
    std::vector<float> totMins = opts.GetOpt<std::vector<float> >(Form("%s.totMins",ch.c_str()));
    int iter = 0;
    for(auto Vov : Vovs)
    {
      cut_totAcc[ch][Vov] = totMins.at(iter);
      ++iter;
    }
  }
  std::map<std::string,float> cut_totMin;
  std::map<std::string,float> cut_totMax;
  
  
  //--- define branches
  float step1, step2;
  std::vector<float>* tot = 0;
  std::vector<long long>* time = 0;
  std::map<std::string,int> channelIdx;
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("step1",1); tree -> SetBranchAddress("step1",&step1);
  tree -> SetBranchStatus("step2",1); tree -> SetBranchAddress("step2",&step2);
  tree -> SetBranchStatus("tot",  1); tree -> SetBranchAddress("tot",  &tot);
  tree -> SetBranchStatus("time", 1); tree -> SetBranchAddress("time", &time);
  for(auto ch : channels)
  {
    tree -> SetBranchStatus(Form("%s",ch.c_str()), 1);  tree -> SetBranchAddress(Form("%s",ch.c_str()), &channelIdx[ch]);
  }
  
  
  //--- define histograms
  std::map<std::string,int> VovLabels;
  std::map<std::string,int> thLabels;
  std::vector<std::string> stepLabels;
  std::map<std::string,float> map_Vovs;
  std::map<std::string,float> map_ths;
  
  std::map<std::string,TH1F*> h1_tot;
  std::map<std::string,TH1F*> h1_tot_cut;
  std::map<std::string,TH2F*> h2_tot_corr;
  std::map<std::string,TH1F*> h1_totRatio;
  
  std::map<std::string,TH1F*> h1_deltaT_raw;
  std::map<std::string,TH1F*> h1_deltaT;
  std::map<std::string,TProfile*> p1_deltaT_vs_totRatio;

  std::map<std::string,TH1F*> h1_deltaT_totCorr;  
  
  
  
  
  //------------------------
  //--- 1st loop over events
  std::map<std::string,std::vector<EventSingle> > eventsSingle;
  std::map<std::string,std::vector<EventCoinc> > eventsCoinc;
  
  int nEntries = tree->GetEntries();
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    tree -> GetEntry(entry);
    if( entry%10000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
    
    step2 = float(int(step2/100)-1);
    std::string VovLabel = Form("Vov%.1f",step1);
    std::string thLabel = Form("th%02.0f",step2);
    std::string stepLabel = Form("Vov%.1f_th%02.0f",step1,step2);
    
        
    //--- create histograms, if needed
    for(auto ch : channels)
    {
      std::string label = Form("%s_%s",ch.c_str(),stepLabel.c_str());
      
      if( h1_tot[label] == NULL )
      {
        h1_tot[label] = new TH1F(Form("h1_tot_%s",label.c_str()),"",2000,0.,1000.);
        h1_tot_cut[label] = new TH1F(Form("h1_tot_cut_%s",label.c_str()),"",2000,0.,1000.);
        
        VovLabels[VovLabel] += 1;
        thLabels[thLabel] += 1;
        stepLabels.push_back(stepLabel);
        map_Vovs[stepLabel] = step1;
        map_ths[stepLabel] = step2;
      }
    }
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
      
      if( h2_tot_corr[label12] == NULL )
      {
        h2_tot_corr[label12] = new TH2F(Form("h2_tot_corr_%s",label12.c_str()),"",200,0.,1000.,200,0.,1000.);
        h1_totRatio[label12] = new TH1F(Form("h1_totRatio_%s",label12.c_str()),"",1000,0.,5.);
        
        h1_deltaT_raw[label12] = new TH1F(Form("h1_deltaT_raw_%s",label12.c_str()),"",1000,35000.,65000.);
        h1_deltaT[label12] = new TH1F(Form("h1_deltaT_%s",label12.c_str()),"",1000,35000.,65000.);
        p1_deltaT_vs_totRatio[label12] = new TProfile(Form("p1_deltaT_vs_totRatio_%s",label12.c_str()),"",1000,0.,5.);
        
        h1_deltaT_totCorr[label12] = new TH1F(Form("h1_deltaT_totCorr_%s",label12.c_str()),"",1000,35000.,65000.);
      }
    }
    for(auto bar : barsVec)
    {
      std::string ch1 = bar.first.first;
      std::string ch2 = bar.first.second;
      std::string ch3 = bar.second.first;
      std::string ch4 = bar.second.second;
      
      std::string label1234 = Form("%s+%s-%s-%s+%s",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str(),stepLabel.c_str());
      
      if( h1_deltaT[label1234] == NULL )
      { 
        h1_totRatio[label1234] = new TH1F(Form("h1_totRatio_%s",label1234.c_str()),"",1000,0.,5.);
        h2_tot_corr[label1234] = new TH2F(Form("h2_tot_corr_%s",label1234.c_str()),"",200,0.,1000.,200,0.,1000.);
        h1_deltaT_raw[label1234] = new TH1F(Form("h1_deltaT_raw_%s",label1234.c_str()),"",1000,-5000.,5000.);
        h1_deltaT[label1234] = new TH1F(Form("h1_deltaT_%s",label1234.c_str()),"",250,-5000.,5000.);
        
        p1_deltaT_vs_totRatio[label1234] = new TProfile(Form("p1_deltaT_vs_totRatio_%s",label1234.c_str()),"",1000,0.,5.);
        
        h1_deltaT_totCorr[label1234] = new TH1F(Form("h1_deltaT_totCorr_%s",label1234.c_str()),"",250,-5000.,5000.);
      }
    }
    
    
    //--- fill histograms
    for(auto ch1 : channels)
    {
      std::string label1 = Form("%s_%s",ch1.c_str(),stepLabel.c_str());
      int idx1 = channelIdx[ch1];
      
      float tot1 = idx1 >= 0 ? tot->at(idx1)/1000. : -1.;
      
      if( idx1 >= 0 ) h1_tot[label1] -> Fill( tot1 );
    }
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      std::string label1 = Form("%s_%s",ch1.c_str(),stepLabel.c_str());
      std::string label2 = Form("%s_%s",ch2.c_str(),stepLabel.c_str());
      std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
      
      int idx1 = channelIdx[ch1];
      int idx2 = channelIdx[ch2];
      
      float tot1 = idx1 >= 0 ? tot->at(idx1)/1000. : -1.;
      float tot2 = idx2 >= 0 ? tot->at(idx2)/1000. : -1.;
      // long long time1 = idx1 >= 0 ? time->at(idx1) : -1.;
      // long long time2 = idx2 >= 0 ? time->at(idx2) : -1.;
      long long time1 = idx1 >= 0 ? time->at(idx1)+tot1*1000. : -1.;
      long long time2 = idx2 >= 0 ? time->at(idx2)+tot2*1000. : -1.;
      
      // std::cout << "delta: " << time2-time1 << "   tot1: " << tot1*1000 << "   tot2: " << tot2*1000 << std::endl;
      
      if( idx1 >= 0 && idx2 >= 0 )
      {
        h2_tot_corr[label12] -> Fill( tot1,tot2 );
        
        EventSingle anEvent;
        anEvent.stepLabel = stepLabel;
        anEvent.ch1 = ch1;
        anEvent.ch2 = ch2;
        anEvent.label1 = label1;
        anEvent.label2 = label2;
        anEvent.label12 = label12;
        anEvent.tot1 = tot1;
        anEvent.tot2 = tot2;
        anEvent.time1 = time1;
        anEvent.time2 = time2;
        eventsSingle[label12].push_back(anEvent);
      }
    }
    for(auto bar : barsVec)
    {
      std::string ch1 = bar.first.first;
      std::string ch2 = bar.first.second;
      std::string label1 = Form("%s_%s",ch1.c_str(),stepLabel.c_str());
      std::string label2 = Form("%s_%s",ch2.c_str(),stepLabel.c_str());
      int idx1 = channelIdx[ch1];
      int idx2 = channelIdx[ch2];
      
      std::string ch3 = bar.second.first;
      std::string ch4 = bar.second.second;
      std::string label3 = Form("%s_%s",ch3.c_str(),stepLabel.c_str());
      std::string label4 = Form("%s_%s",ch4.c_str(),stepLabel.c_str());
      int idx3 = channelIdx[ch3];
      int idx4 = channelIdx[ch4];
      
      std::string label1234 = Form("%s+%s-%s-%s+%s",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str(),stepLabel.c_str());
      
      float tot1 = idx1 >= 0 ? tot->at(idx1)/1000. : -1.;
      float tot2 = idx2 >= 0 ? tot->at(idx2)/1000. : -1.;
      float tot3 = idx3 >= 0 ? tot->at(idx3)/1000. : -1.;
      float tot4 = idx4 >= 0 ? tot->at(idx4)/1000. : -1.;
      // long long time1 = idx1 >= 0 ? time->at(idx1) : -1.;
      // long long time2 = idx2 >= 0 ? time->at(idx2) : -1.;
      // long long time3 = idx3 >= 0 ? time->at(idx3) : -1.;
      // long long time4 = idx4 >= 0 ? time->at(idx4) : -1.;
      long long time1 = idx1 >= 0 ? time->at(idx1)+tot1*1000. : -1.;
      long long time2 = idx2 >= 0 ? time->at(idx2)+tot2*1000. : -1.;
      long long time3 = idx3 >= 0 ? time->at(idx3)+tot3*1000. : -1.;
      long long time4 = idx4 >= 0 ? time->at(idx4)+tot4*1000. : -1.;
      
      if( idx1 >= 0 && idx2 >= 0 && idx3 >= 0 && idx4 >= 0 )
      {
        h2_tot_corr[label1234] -> Fill( 0.5*(tot1+tot2),0.5*(tot3+tot4) );
        
        EventCoinc anEvent;
        anEvent.stepLabel = stepLabel;
        anEvent.ch1 = ch1;
        anEvent.ch2 = ch2;
        anEvent.ch3 = ch3;
        anEvent.ch4 = ch4;
        anEvent.label1 = label1;
        anEvent.label2 = label2;
        anEvent.label3 = label3;
        anEvent.label4 = label4;
        anEvent.label1234 = label1234;
        anEvent.tot1 = tot1;
        anEvent.tot2 = tot2;
        anEvent.tot3 = tot3;
        anEvent.tot4 = tot4;
        anEvent.time1 = time1;
        anEvent.time2 = time2;
        anEvent.time3 = time3;
        anEvent.time4 = time4;
        eventsCoinc[label1234].push_back(anEvent);
      }
    }
  }
  std::cout << std::endl;
  
  std::vector<std::string>::iterator iter;
  iter = std::unique(stepLabels.begin(),stepLabels.end());
  stepLabels.resize( std::distance(stepLabels.begin(),iter) );  
  
  std::sort(stepLabels.begin(),stepLabels.end());
  
  
  
  
  //------------------
  //--- draw 1st plots
  TCanvas* c;
  float* vals = new float[6];
  TLatex* latex;
  
  std::map<std::string,TGraphErrors*> g_tot_vs_th;
  std::map<std::string,TGraphErrors*> g_tot_vs_Vov;
  
  for(auto stepLabel : stepLabels)
  {
    float Vov = map_Vovs[stepLabel];
    float th = map_ths[stepLabel];
    std::string VovLabel(Form("Vov%.1f",Vov));
    std::string thLabel(Form("th%02.0f",th));
    
    
    //--------------------------------------------------------
    
    
    for(auto ch : channels)
    {
      std::string label = ch + "_" + stepLabel;
      
      c = new TCanvas(Form("c_tot_%s",label.c_str()),Form("c_tot_%s",label.c_str()));
      // gPad -> SetLogy();
      
      h1_tot[label] -> SetTitle(";ToT [ns];entries");
      h1_tot[label] -> SetLineColor(kRed);
      h1_tot[label] -> Draw();
      float max1 = FindXMaximum(h1_tot[label],cut_totAcc[ch][Vov],1000.);
      h1_tot[label] -> GetXaxis() -> SetRangeUser(0.1*max1,2.5*max1);
      TF1* fitFunc1 = new TF1("fitFunc1","pol2",max1-0.075*max1,max1+0.075*max1);
      h1_tot[label] -> Fit(fitFunc1,"QNRS+");
      fitFunc1 -> SetLineColor(kRed);
      fitFunc1 -> SetLineWidth(3);
      fitFunc1 -> Draw("same");
      cut_totMin[ch+"_"+stepLabel] = 0.9*fitFunc1->GetMaximumX();
      cut_totMax[ch+"_"+stepLabel] = 1.1*fitFunc1->GetMaximumX();
      TLine* line_totMin1 = new TLine(cut_totMin[ch+"_"+stepLabel],h1_tot[label]->GetMinimum(),cut_totMin[ch+"_"+stepLabel],h1_tot[label]->GetMaximum());
      line_totMin1 -> SetLineColor(kRed);
      line_totMin1 -> Draw("same");
      TLine* line_totMax1 = new TLine(cut_totMax[ch+"_"+stepLabel],h1_tot[label]->GetMinimum(),cut_totMax[ch+"_"+stepLabel],h1_tot[label]->GetMaximum());
      line_totMax1 -> SetLineColor(kRed);
      line_totMax1 -> Draw("same");
      latex = new TLatex(0.65,0.85,Form("%s",ch.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");      
      c -> Print(Form("%s/c_tot_%s.png",plotDir.c_str(),label.c_str()));
      c -> Print(Form("%s/c_tot_%s.pdf",plotDir.c_str(),label.c_str()));
      
      
      if( g_tot_vs_th[ch+"_"+VovLabel] == NULL )
        g_tot_vs_th[ch+"_"+VovLabel] = new TGraphErrors();
      
      if( g_tot_vs_Vov[ch+"_"+thLabel] == NULL )
        g_tot_vs_Vov[ch+"_"+thLabel] = new TGraphErrors();
      
      g_tot_vs_th[ch+"_"+VovLabel] -> SetPoint(g_tot_vs_th[ch+"_"+VovLabel]->GetN(),63-th,fitFunc1->GetMaximumX());
      g_tot_vs_th[ch+"_"+VovLabel] -> SetPointError(g_tot_vs_th[ch+"_"+VovLabel]->GetN()-1,0.,0.);
      
      g_tot_vs_Vov[ch+"_"+thLabel] -> SetPoint(g_tot_vs_Vov[ch+"_"+thLabel]->GetN(),Vov,fitFunc1->GetMaximumX());
      g_tot_vs_Vov[ch+"_"+thLabel] -> SetPointError(g_tot_vs_Vov[ch+"_"+thLabel]->GetN()-1,0.,0.);
    }
    
    
    //--------------------------------------------------------
    
    
    for(auto pair : pairsVec)
    {  
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
      
      c = new TCanvas(Form("c_tot_corr_%s",label12.c_str()),Form("c_tot_corr_%s",label12.c_str()));
      gPad -> SetLogz();
      
      h2_tot_corr[label12] -> SetTitle(Form(";%s ToT [ns];%s ToT [ns]",ch1.c_str(),ch2.c_str()));
      h2_tot_corr[label12] -> Draw("colz");
      
      c -> Print(Form("%s/c_tot_corr_%s.png",plotDir.c_str(),label12.c_str()));
      c -> Print(Form("%s/c_tot_corr_%s.pdf",plotDir.c_str(),label12.c_str()));
    }
    
    
    //--------------------------------------------------------
    
    
    for(auto bar : barsVec)
    {
      std::string ch1 = bar.first.first;
      std::string ch2 = bar.first.second;
      std::string label1 = Form("%s_%s",ch1.c_str(),stepLabel.c_str());
      std::string label2 = Form("%s_%s",ch2.c_str(),stepLabel.c_str());
      
      std::string ch3 = bar.second.first;
      std::string ch4 = bar.second.second;
      std::string label3 = Form("%s_%s",ch3.c_str(),stepLabel.c_str());
      std::string label4 = Form("%s_%s",ch4.c_str(),stepLabel.c_str());
      
      std::string label1234 = Form("%s+%s-%s-%s+%s",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str(),stepLabel.c_str());
      
      c = new TCanvas(Form("c_tot_corr_%s",label1234.c_str()),Form("c_tot_corr_%s",label1234.c_str()));
      gPad -> SetLogz();
      
      h2_tot_corr[label1234] -> SetTitle(Form(";%s+%s ToT [ns];%s+%s ToT [ns]",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str()));
      h2_tot_corr[label1234] -> Draw("colz");
      
      c -> Print(Form("%s/c_tot_corr_%s.png",plotDir.c_str(),label1234.c_str()));
      c -> Print(Form("%s/c_tot_corr_%s.pdf",plotDir.c_str(),label1234.c_str()));
    }
  }
  
  
  //--------------------------------------------------------
  
  
  for(auto pair : pairsVec)
  {
    std::string ch1 = pair.first;
    std::string ch2 = pair.second;
    
    c = new TCanvas(Form("c_tot_vs_th_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_tot_vs_th_%s-%s",ch1.c_str(),ch2.c_str()));
    // gPad -> SetLogy();
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,500.) );
    hPad -> SetTitle(";63 - threshold [DAC];ToT [ns]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    int iter = 0;
    for(auto mapIt : VovLabels)
    {
      std::string label1 = ch1+"_"+mapIt.first;
      std::string label2 = ch2+"_"+mapIt.first;
      TGraph* g_tot1 = g_tot_vs_th[label1];
      TGraph* g_tot2 = g_tot_vs_th[label2];
      
      g_tot1 -> SetLineColor(1+iter);
      g_tot1 -> SetMarkerColor(1+iter);
      g_tot1 -> SetMarkerStyle(20);
      g_tot1 -> Draw("PL,same");
      
      g_tot2 -> SetLineColor(1+iter);
      g_tot2 -> SetMarkerColor(1+iter);
      g_tot2 -> SetMarkerStyle(25);
      g_tot2 -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_tot_vs_th_%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
    c -> Print(Form("%s/c_tot_vs_th_%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
    
    
    c = new TCanvas(Form("c_tot_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_tot_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()));
    // gPad -> SetLogy();
    
    hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,500.) );
    hPad -> SetTitle(";V_{ov} [V];ToT [ns]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    iter = 0;
    for(auto mapIt : thLabels)
    {
      std::string label1 = ch1+"_"+mapIt.first;
      std::string label2 = ch2+"_"+mapIt.first;
      TGraph* g_tot1 = g_tot_vs_Vov[label1];
      TGraph* g_tot2 = g_tot_vs_Vov[label2];
      
      g_tot1 -> SetLineColor(1+iter);
      g_tot1 -> SetMarkerColor(1+iter);
      g_tot1 -> SetMarkerStyle(20);
      g_tot1 -> Draw("PL,same");
      
      g_tot2 -> SetLineColor(1+iter);
      g_tot2 -> SetMarkerColor(1+iter);
      g_tot2 -> SetMarkerStyle(25);
      g_tot2 -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_tot_vs_Vov_%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
    c -> Print(Form("%s/c_tot_vs_Vov_%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
  }
  
  
  
  
  //------------------------
  //--- 2nd loop over events
  std::map<std::string,std::vector<EventSingle> > eventsSingle2;
  std::map<std::string,std::vector<EventCoinc> > eventsCoinc2;
  
  for(auto mapIt : eventsSingle)
  {
    std::string label = mapIt.first;
    
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 2nd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      EventSingle anEvent = mapIt.second.at(entry);
      
      if( anEvent.tot1 > cut_totMin[anEvent.ch1+"_"+anEvent.stepLabel] && anEvent.tot1 < cut_totMax[anEvent.ch1+"_"+anEvent.stepLabel] ) h1_tot_cut[anEvent.label1] -> Fill( anEvent.tot1 );
      if( anEvent.tot2 > cut_totMin[anEvent.ch2+"_"+anEvent.stepLabel] && anEvent.tot2 < cut_totMax[anEvent.ch2+"_"+anEvent.stepLabel] ) h1_tot_cut[anEvent.label2] -> Fill( anEvent.tot2 );
      
      if( anEvent.tot1 > cut_totMin[anEvent.ch1+"_"+anEvent.stepLabel] &&
          anEvent.tot2 > cut_totMin[anEvent.ch2+"_"+anEvent.stepLabel] &&
          anEvent.tot1 < cut_totMax[anEvent.ch1+"_"+anEvent.stepLabel] &&
          anEvent.tot2 < cut_totMax[anEvent.ch2+"_"+anEvent.stepLabel] )
      {
        h1_totRatio[anEvent.label12] -> Fill( anEvent.tot2 / anEvent.tot1 );
        
        h1_deltaT_raw[anEvent.label12] -> Fill( anEvent.time2-anEvent.time1 );
        
        EventSingle anEvent2;
        anEvent2.stepLabel = anEvent.stepLabel;
        anEvent2.ch1 = anEvent.ch1;
        anEvent2.ch2 = anEvent.ch2;
        anEvent2.label1 = anEvent.label1;
        anEvent2.label2 = anEvent.label2;
        anEvent2.label12 = anEvent.label12;
        anEvent2.tot1 = anEvent.tot1;
        anEvent2.tot2 = anEvent.tot2;
        anEvent2.time1 = anEvent.time1;
        anEvent2.time2 = anEvent.time2;
        eventsSingle2[anEvent.label12].push_back(anEvent2);
      }
    }
    std::cout << std::endl;
  }
  
  for(auto mapIt : eventsCoinc)
  {
    std::string label = mapIt.first;
    
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 2nd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      EventCoinc anEvent = mapIt.second.at(entry);
      
      if( anEvent.tot1 > cut_totMin[anEvent.ch1+"_"+anEvent.stepLabel] && 
          anEvent.tot2 > cut_totMin[anEvent.ch2+"_"+anEvent.stepLabel] && 
          anEvent.tot3 > cut_totMin[anEvent.ch3+"_"+anEvent.stepLabel] && 
          anEvent.tot4 > cut_totMin[anEvent.ch4+"_"+anEvent.stepLabel] &&
          anEvent.tot1 < cut_totMax[anEvent.ch1+"_"+anEvent.stepLabel] && 
          anEvent.tot2 < cut_totMax[anEvent.ch2+"_"+anEvent.stepLabel] && 
          anEvent.tot3 < cut_totMax[anEvent.ch3+"_"+anEvent.stepLabel] && 
          anEvent.tot4 < cut_totMax[anEvent.ch4+"_"+anEvent.stepLabel] )
      {
        h1_totRatio[anEvent.label1234] -> Fill( (anEvent.tot3+anEvent.tot4) / (anEvent.tot1+anEvent.tot2) );
        
        h1_deltaT_raw[anEvent.label1234] -> Fill( 0.5*(anEvent.time3+anEvent.time4) - 0.5*(anEvent.time1+anEvent.time2) );
        
        EventCoinc anEvent2;
        anEvent2.stepLabel = anEvent.stepLabel;
        anEvent2.ch1 = anEvent.ch1;
        anEvent2.ch2 = anEvent.ch2;
        anEvent2.ch3 = anEvent.ch3;
        anEvent2.ch4 = anEvent.ch4;
        anEvent2.label1 = anEvent.label1;
        anEvent2.label2 = anEvent.label2;
        anEvent2.label3 = anEvent.label3;
        anEvent2.label4 = anEvent.label4;
        anEvent2.label1234 = anEvent.label1234;
        anEvent2.tot1 = anEvent.tot1;
        anEvent2.tot2 = anEvent.tot2;
        anEvent2.tot3 = anEvent.tot3;
        anEvent2.tot4 = anEvent.tot4;
        anEvent2.time1 = anEvent.time1;
        anEvent2.time2 = anEvent.time2;
        anEvent2.time3 = anEvent.time3;
        anEvent2.time4 = anEvent.time4;
        eventsCoinc2[anEvent.label1234].push_back(anEvent2);
      }
    }
    std::cout << std::endl;
  }
  
  
  
  
  //------------------
  //--- draw 2nd plots
  std::map<std::string,float> CTRMeans;
  std::map<std::string,float> CTRSigmas;
  
  for(auto stepLabel : stepLabels)
  {
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      std::string label1 = ch1 + "_" + stepLabel;
      std::string label2 = ch2 + "_" + stepLabel;
      std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
      
      
      //--------------------------------------------------------
      
      
      FindSmallestInterval(vals,h1_deltaT_raw[label12],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans[label12] = mean;
      CTRSigmas[label12] = effSigma;
      
      
      //--------------------------------------------------------
      
      
      c = new TCanvas(Form("c_totRatio_%s",label12.c_str()),Form("c_totRatio_%s",label12.c_str()));
      // gPad -> SetLogy();
      
      h1_totRatio[label12] -> SetTitle(Form(";%s ToT / %s ToT;entries",ch1.c_str(),ch2.c_str()));
      h1_totRatio[label12] -> Draw("colz");
      
      TF1* fitFunc = new TF1(Form("fitFunc_totRatio_%s",label12.c_str()),"gaus",
                             h1_totRatio[label12]->GetMean()-1.5*h1_totRatio[label12]->GetRMS(),
                             h1_totRatio[label12]->GetMean()+1.5*h1_totRatio[label12]->GetRMS());
      fitFunc -> SetParameters(1,h1_totRatio[label12]->GetMean(),h1_totRatio[label12]->GetRMS());
      h1_totRatio[label12] -> Fit(fitFunc,"QNRSL+");
      fitFunc -> SetLineColor(kRed);
      fitFunc -> SetLineWidth(1);
      fitFunc -> Draw("same");
      
      latex = new TLatex(0.55,0.85,Form("#sigma = %.1f %%",100.*fitFunc->GetParameter(2)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");
      
      c -> Print(Form("%s/c_totRatio_%s.png",plotDir.c_str(),label12.c_str()));
      c -> Print(Form("%s/c_totRatio_%s.pdf",plotDir.c_str(),label12.c_str()));
    }
    
    
    for(auto bar : barsVec)
    {
      std::string ch1 = bar.first.first;
      std::string ch2 = bar.first.second;
      std::string label1 = Form("%s_%s",ch1.c_str(),stepLabel.c_str());
      std::string label2 = Form("%s_%s",ch2.c_str(),stepLabel.c_str());

      std::string ch3 = bar.second.first;
      std::string ch4 = bar.second.second;
      std::string label3 = Form("%s_%s",ch3.c_str(),stepLabel.c_str());
      std::string label4 = Form("%s_%s",ch4.c_str(),stepLabel.c_str());

      std::string label1234 = Form("%s+%s-%s-%s+%s",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str(),stepLabel.c_str());
      
      
      //--------------------------------------------------------
      
      
      FindSmallestInterval(vals,h1_deltaT_raw[label1234],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans[label1234] = mean;
      CTRSigmas[label1234] = effSigma;
      
      
      //--------------------------------------------------------
      
      
      c = new TCanvas(Form("c_totRatio_%s",label1234.c_str()),Form("c_totRatio_%s",label1234.c_str()));
      gPad -> SetLogy();
      
      h1_totRatio[label1234] -> SetTitle(Form(";%s+%s ToT / %s+%s ToT;entries",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str()));
      h1_totRatio[label1234] -> Draw("colz");
      
      c -> Print(Form("%s/c_totRatio_%s.png",plotDir.c_str(),label1234.c_str()));
      c -> Print(Form("%s/c_totRatio_%s.pdf",plotDir.c_str(),label1234.c_str()));
    }
  }
  
  
  
  
  //------------------------
  //--- 3rd loop over events
  for(auto mapIt : eventsSingle)
  {
    std::string label = mapIt.first;
    
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 3rd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      EventSingle anEvent = mapIt.second.at(entry);
      
      if( anEvent.tot1 > cut_totMin[anEvent.ch1+"_"+anEvent.stepLabel] &&
          anEvent.tot2 > cut_totMin[anEvent.ch2+"_"+anEvent.stepLabel] &&
          anEvent.tot1 < cut_totMax[anEvent.ch1+"_"+anEvent.stepLabel] &&
          anEvent.tot2 < cut_totMax[anEvent.ch2+"_"+anEvent.stepLabel] )
      {
        float timeLow = CTRMeans[anEvent.label12] - 2.* CTRSigmas[anEvent.label12];
        float timeHig = CTRMeans[anEvent.label12] + 2.* CTRSigmas[anEvent.label12];
        float timeDiff = anEvent.time2 - anEvent.time1;
        h1_deltaT[anEvent.label12] -> Fill( timeDiff );
        if( ( timeDiff > timeLow ) &&
            ( timeDiff < timeHig ) )
          p1_deltaT_vs_totRatio[anEvent.label12] -> Fill( anEvent.tot2/anEvent.tot1,anEvent.time2-anEvent.time1 );    
      }
    }
    std::cout << std::endl;
  }
  
  for(auto mapIt : eventsCoinc)
  {
    std::string label = mapIt.first;
    
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 3rd loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
      EventCoinc anEvent = mapIt.second.at(entry);
      
      if( anEvent.tot1 > cut_totMin[anEvent.ch1+"_"+anEvent.stepLabel] && 
          anEvent.tot2 > cut_totMin[anEvent.ch2+"_"+anEvent.stepLabel] && 
          anEvent.tot3 > cut_totMin[anEvent.ch3+"_"+anEvent.stepLabel] && 
          anEvent.tot4 > cut_totMin[anEvent.ch4+"_"+anEvent.stepLabel] &&
          anEvent.tot1 < cut_totMax[anEvent.ch1+"_"+anEvent.stepLabel] && 
          anEvent.tot2 < cut_totMax[anEvent.ch2+"_"+anEvent.stepLabel] && 
          anEvent.tot3 < cut_totMax[anEvent.ch3+"_"+anEvent.stepLabel] && 
          anEvent.tot4 < cut_totMax[anEvent.ch4+"_"+anEvent.stepLabel] )
      {
        float timeLow = CTRMeans[anEvent.label1234] - 2.* CTRSigmas[anEvent.label1234];
        float timeHig = CTRMeans[anEvent.label1234] + 2.* CTRSigmas[anEvent.label1234];
        
        float timeComb = 0.5*(anEvent.time3+anEvent.time4) - 0.5*(anEvent.time1+anEvent.time2);
        h1_deltaT[anEvent.label1234] -> Fill( timeComb );
        if( ( timeComb > timeLow ) &&
            ( timeComb < timeHig ) )
          p1_deltaT_vs_totRatio[anEvent.label1234] -> Fill( (anEvent.tot3+anEvent.tot4)/(anEvent.tot1+anEvent.tot2),timeComb );
      }
    }
    std::cout << std::endl;
  }
  
  
  
  
  //------------------
  //--- draw 3rd plots
  std::map<std::string,TF1*> fitFunc_totCorr;
  
  for(auto stepLabel : stepLabels)
  {
    float Vov = map_Vovs[stepLabel];
    float th = map_ths[stepLabel];
    std::string VovLabel(Form("Vov%.1f",Vov));
    std::string thLabel(Form("th%02.0f",th));
    
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      std::string label1 = ch1 + "_" + stepLabel;
      std::string label2 = ch2 + "_" + stepLabel;
      std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
      
      
      //--------------------------------------------------------
      
      
      c = new TCanvas(Form("c_deltaT_vs_totRatio_%s",label12.c_str()),Form("c_deltaT_vs_totRatio_%s",label12.c_str()));
      
      p1_deltaT_vs_totRatio[label12] -> SetTitle(Form(";%s ToT / %s ToT;#Deltat [ps]",ch2.c_str(),ch1.c_str()));
      p1_deltaT_vs_totRatio[label12] -> GetXaxis() -> SetRangeUser(h1_totRatio[label12]->GetMean()-3.*h1_totRatio[label12]->GetRMS(),
                                                                   h1_totRatio[label12]->GetMean()+3.*h1_totRatio[label12]->GetRMS());
      p1_deltaT_vs_totRatio[label12] -> Draw("");
      
      float fitXMin = h1_totRatio[label12]->GetMean() - 2.*h1_totRatio[label12]->GetRMS();
      float fitXMax = h1_totRatio[label12]->GetMean() + 2.*h1_totRatio[label12]->GetRMS();
      fitFunc_totCorr[label12] = new TF1(Form("fitFunc_totCorr_%s",label12.c_str()),"pol4",fitXMin,fitXMax);
      p1_deltaT_vs_totRatio[label12] -> Fit(fitFunc_totCorr[label12],"QNRS+");
      fitFunc_totCorr[label12] -> SetLineColor(kRed);
      fitFunc_totCorr[label12] -> SetLineWidth(2);
      fitFunc_totCorr[label12] -> Draw("same");
      
      c -> Print(Form("%s/c_deltaT_vs_totRatio_%s.png",plotDir.c_str(),label12.c_str()));
      c -> Print(Form("%s/c_deltaT_vs_totRatio_%s.pdf",plotDir.c_str(),label12.c_str()));
    }
    
    
    for(auto bar : barsVec)
    {
      std::string ch1 = bar.first.first;
      std::string ch2 = bar.first.second;
      std::string label1 = Form("%s_%s",ch1.c_str(),stepLabel.c_str());
      std::string label2 = Form("%s_%s",ch2.c_str(),stepLabel.c_str());

      std::string ch3 = bar.second.first;
      std::string ch4 = bar.second.second;
      std::string label3 = Form("%s_%s",ch3.c_str(),stepLabel.c_str());
      std::string label4 = Form("%s_%s",ch4.c_str(),stepLabel.c_str());

      std::string label1234 = Form("%s+%s-%s-%s+%s",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str(),stepLabel.c_str());
      
      
      //--------------------------------------------------------
      
      
      c = new TCanvas(Form("c_deltaT_vs_totRatio_%s",label1234.c_str()),Form("c_deltaT_vs_totRatio_%s",label1234.c_str()));
      
      p1_deltaT_vs_totRatio[label1234] -> SetTitle(Form(";%s+%s ToT / %s+%s ToT;#Deltat [ps]",ch3.c_str(),ch4.c_str(),ch1.c_str(),ch2.c_str()));
      p1_deltaT_vs_totRatio[label1234] -> GetXaxis() -> SetRangeUser(h1_totRatio[label1234]->GetMean()-3.*h1_totRatio[label1234]->GetRMS(),
                                                                     h1_totRatio[label1234]->GetMean()+3.*h1_totRatio[label1234]->GetRMS());
      p1_deltaT_vs_totRatio[label1234] -> Draw("");
      
      float fitXMin = h1_totRatio[label1234]->GetMean() - 2.*h1_totRatio[label1234]->GetRMS();
      float fitXMax = h1_totRatio[label1234]->GetMean() + 2.*h1_totRatio[label1234]->GetRMS();
      fitFunc_totCorr[label1234] = new TF1(Form("fitFunc_totCorr_%s",label1234.c_str()),"pol4",fitXMin,fitXMax);
      p1_deltaT_vs_totRatio[label1234] -> Fit(fitFunc_totCorr[label1234],"QNRS+");
      fitFunc_totCorr[label1234] -> SetLineColor(kRed);
      fitFunc_totCorr[label1234] -> SetLineWidth(2);
      fitFunc_totCorr[label1234] -> Draw("same");
      
      c -> Print(Form("%s/c_deltaT_vs_totRatio_%s.png",plotDir.c_str(),label1234.c_str()));
      c -> Print(Form("%s/c_deltaT_vs_totRatio_%s.pdf",plotDir.c_str(),label1234.c_str()));
    }
  }
  
  
  
  
  //------------------------
  //--- 4th loop over events
  for(auto mapIt : eventsSingle2)
  {
    std::string label = mapIt.first;
    
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 4th loop (" << label << "): reading entry " << entry << " / " << nEntries << "\r" << std::flush;
      EventSingle anEvent = mapIt.second.at(entry);
      
      float deltaT = anEvent.time2-anEvent.time1;
      float totCorr = fitFunc_totCorr[label]->Eval(anEvent.tot2/anEvent.tot1) - fitFunc_totCorr[label]->Eval(h1_tot_cut[anEvent.label2]->GetMean()/h1_tot_cut[anEvent.label1]->GetMean());
      h1_deltaT_totCorr[label] -> Fill( deltaT - totCorr );
    }
    std::cout << std::endl;
  }
  
  for(auto mapIt : eventsCoinc2)
  {
    std::string label = mapIt.first;
    
    nEntries = mapIt.second.size();
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%1000 == 0 ) std::cout << ">>> 4th loop (" << label << "): reading entry " << entry << " / " << nEntries << "\r" << std::flush;
      EventCoinc anEvent = mapIt.second.at(entry);
      
      float deltaT = 0.5*(anEvent.time3+anEvent.time4) - 0.5*(anEvent.time1+anEvent.time2);
      float totCorr = fitFunc_totCorr[label]->Eval((anEvent.tot3+anEvent.tot4)/(anEvent.tot1+anEvent.tot2)) - fitFunc_totCorr[label]->Eval((h1_tot_cut[anEvent.label3]->GetMean()+h1_tot_cut[anEvent.label4]->GetMean())/(h1_tot_cut[anEvent.label1]->GetMean()+h1_tot_cut[anEvent.label2]->GetMean()));
      h1_deltaT_totCorr[label] -> Fill( deltaT - totCorr );
    }
    std::cout << std::endl;
  }
  
  
  
  //------------------
  //--- draw 4th plots
  std::map<std::string,TGraphErrors*> g_tRes_effSigma_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_gaus_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_effSigma_vs_Vov;
  std::map<std::string,TGraphErrors*> g_tRes_gaus_vs_Vov;
  
  std::map<std::string,TGraphErrors*> g_tRes_totCorr_effSigma_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_totCorr_gaus_vs_th;
  std::map<std::string,TGraphErrors*> g_tRes_totCorr_effSigma_vs_Vov;
  std::map<std::string,TGraphErrors*> g_tRes_totCorr_gaus_vs_Vov;
  
  for(auto stepLabel : stepLabels)
  {
    float Vov = map_Vovs[stepLabel];
    float th = map_ths[stepLabel];
    std::string VovLabel(Form("Vov%.1f",Vov));
    std::string thLabel(Form("th%02.0f",th));
    
    for(auto pair : pairsVec)
    {
      std::string ch1 = pair.first;
      std::string ch2 = pair.second;
      std::string label1 = ch1 + "_" + stepLabel;
      std::string label2 = ch2 + "_" + stepLabel;
      std::string label12 = Form("%s-%s_%s",ch1.c_str(),ch2.c_str(),stepLabel.c_str());
      
      
      //--------------------------------------------------------
      
      c = new TCanvas(Form("c_deltaT_totCorr_%s",label12.c_str()),Form("c_deltaT_totCorr_%s",label12.c_str()));
      
      h1_deltaT_totCorr[label12] -> GetXaxis() -> SetRangeUser(h1_deltaT_totCorr[label12]->GetMean()-5.*h1_deltaT_totCorr[label12]->GetRMS(),
                                                               h1_deltaT_totCorr[label12]->GetMean()+5.*h1_deltaT_totCorr[label12]->GetRMS());          
      h1_deltaT_totCorr[label12] -> SetTitle(Form(";tot-corrected #Deltat [ps];entries"));
      h1_deltaT_totCorr[label12] -> SetLineWidth(2);
      h1_deltaT_totCorr[label12] -> SetLineColor(kBlue);
      h1_deltaT_totCorr[label12] -> SetMarkerColor(kBlue);
      h1_deltaT_totCorr[label12] -> Draw("");
      
      float fitXMin = CTRMeans[label12] - 1.5*CTRSigmas[label12];
      float fitXMax = CTRMeans[label12] + 1.5*CTRSigmas[label12];
      TF1* fitFunc = new TF1(Form("fitFunc_totCorr_%s",label12.c_str()),"gaus",fitXMin,fitXMax);
      fitFunc -> SetParameters(1,h1_deltaT_totCorr[label12]->GetMean(),h1_deltaT_totCorr[label12]->GetRMS());
      h1_deltaT_totCorr[label12] -> Fit(fitFunc,"QNRSL+");
      fitFunc -> SetLineColor(kBlue);
      fitFunc -> SetLineWidth(2);
      fitFunc -> Draw("same");
      
      FindSmallestInterval(vals,h1_deltaT_totCorr[label12],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      
      latex = new TLatex(0.55,0.85,Form("time walk corr. #splitline{#sigma_{CTR}^{eff} = %.1f ps}{#sigma_{CTR}^{gaus} = %.1f ps}",effSigma,fitFunc->GetParameter(2)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlue);
      latex -> Draw("same");
      
      if( g_tRes_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel] == NULL )
      {
        g_tRes_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel] = new TGraphErrors();
        g_tRes_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel] = new TGraphErrors();
        g_tRes_totCorr_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel] = new TGraphErrors();
        g_tRes_totCorr_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel] = new TGraphErrors();        
      }
      if( g_tRes_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel] == NULL )
      {
        g_tRes_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel] = new TGraphErrors();
        g_tRes_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel] = new TGraphErrors();
        g_tRes_totCorr_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel] = new TGraphErrors();
        g_tRes_totCorr_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel] = new TGraphErrors();
      } 
      
      if( pairsMode == 1 )
      {
        g_tRes_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPoint(g_tRes_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN(),63-th,effSigma/sqrt(2));
        g_tRes_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPointError(g_tRes_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN()-1,0.,5.);
        g_tRes_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPoint(g_tRes_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN(),63-th,fitFunc->GetParameter(2)/sqrt(2));
        g_tRes_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPointError(g_tRes_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN()-1,0.,fitFunc->GetParError(2)/sqrt(2));
        
        g_tRes_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPoint(g_tRes_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN(),Vov,effSigma/sqrt(2));
        g_tRes_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPointError(g_tRes_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN()-1,0.,5.);
        g_tRes_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPoint(g_tRes_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN(),Vov,fitFunc->GetParameter(2)/sqrt(2));
        g_tRes_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPointError(g_tRes_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN()-1,0.,fitFunc->GetParError(2)/sqrt(2));
      }
      if( pairsMode == 2 )
      {
        g_tRes_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPoint(g_tRes_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN(),63-th,effSigma/2);
        g_tRes_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPointError(g_tRes_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN()-1,0.,5.);
        g_tRes_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPoint(g_tRes_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN(),63-th,fitFunc->GetParameter(2)/2);
        g_tRes_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPointError(g_tRes_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN()-1,0.,fitFunc->GetParError(2)/2.);
        
        g_tRes_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPoint(g_tRes_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN(),Vov,effSigma/2);
        g_tRes_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPointError(g_tRes_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN()-1,0.,5.);
        g_tRes_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPoint(g_tRes_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN(),Vov,fitFunc->GetParameter(2)/2);
        g_tRes_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPointError(g_tRes_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN()-1,0.,fitFunc->GetParError(2)/2.);
      }
      
      
      
      h1_deltaT[label12] -> SetLineWidth(2);
      h1_deltaT[label12] -> SetLineColor(kRed);
      h1_deltaT[label12] -> SetMarkerColor(kRed);
      h1_deltaT[label12] -> Draw("same");
      
      fitXMin = CTRMeans[label12] - 1.5*CTRSigmas[label12];
      fitXMax = CTRMeans[label12] + 1.5*CTRSigmas[label12];
      fitFunc = new TF1(Form("fitFunc_%s",label12.c_str()),"gaus",fitXMin,fitXMax);
      fitFunc -> SetParameters(1,h1_deltaT[label12]->GetMean(),h1_deltaT[label12]->GetRMS());
      h1_deltaT[label12] -> Fit(fitFunc,"QNRSL+");
      fitFunc -> SetLineColor(kRed);
      fitFunc -> SetLineWidth(2);
      fitFunc -> Draw("same");
      
      FindSmallestInterval(vals,h1_deltaT[label12],0.68);
      mean = vals[0];
      min = vals[4];
      max = vals[5];
      delta = max-min;
      sigma = 0.5*delta;
      effSigma = sigma;
      
      latex = new TLatex(0.55,0.65,Form("#splitline{raw #sigma_{CTR}^{eff} = %.1f ps}{raw #sigma_{CTR}^{gaus} = %.1f ps}",effSigma,fitFunc->GetParameter(2)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");
      
      if( pairsMode == 1 )
      {
        g_tRes_totCorr_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPoint(g_tRes_totCorr_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN(),63-th,effSigma/sqrt(2));
        g_tRes_totCorr_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPointError(g_tRes_totCorr_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN()-1,0.,5.);
        g_tRes_totCorr_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPoint(g_tRes_totCorr_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN(),63-th,fitFunc->GetParameter(2)/sqrt(2));
        g_tRes_totCorr_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPointError(g_tRes_totCorr_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN()-1,0.,fitFunc->GetParError(2)/sqrt(2));
        
        g_tRes_totCorr_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPoint(g_tRes_totCorr_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN(),Vov,effSigma/sqrt(2));
        g_tRes_totCorr_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPointError(g_tRes_totCorr_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN()-1,0.,5.);
        g_tRes_totCorr_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPoint(g_tRes_totCorr_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN(),Vov,fitFunc->GetParameter(2)/sqrt(2));
        g_tRes_totCorr_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPointError(g_tRes_totCorr_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN()-1,0.,fitFunc->GetParError(2)/sqrt(2));
      }
      if( pairsMode == 2 )
      {
        g_tRes_totCorr_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPoint(g_tRes_totCorr_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN(),63-th,effSigma/2.);
        g_tRes_totCorr_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPointError(g_tRes_totCorr_effSigma_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN()-1,0.,5.);
        g_tRes_totCorr_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPoint(g_tRes_totCorr_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN(),63-th,fitFunc->GetParameter(2)/2.);
        g_tRes_totCorr_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel] -> SetPointError(g_tRes_totCorr_gaus_vs_th[ch1+"-"+ch2+"_"+VovLabel]->GetN()-1,0.,fitFunc->GetParError(2)/2.);
        
        g_tRes_totCorr_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPoint(g_tRes_totCorr_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN(),Vov,effSigma/2.);
        g_tRes_totCorr_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPointError(g_tRes_totCorr_effSigma_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN()-1,0.,5.);
        g_tRes_totCorr_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPoint(g_tRes_totCorr_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN(),Vov,fitFunc->GetParameter(2)/2.);
        g_tRes_totCorr_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel] -> SetPointError(g_tRes_totCorr_gaus_vs_Vov[ch1+"-"+ch2+"_"+thLabel]->GetN()-1,0.,fitFunc->GetParError(2)/2.);
      }
      
      c -> Print(Form("%s/c_deltaT_totCorr_%s.png",plotDir.c_str(),label12.c_str()));
      c -> Print(Form("%s/c_deltaT_totCorr_%s.pdf",plotDir.c_str(),label12.c_str()));
    }
    
    
    for(auto bar : barsVec)
    {
      std::string ch1 = bar.first.first;
      std::string ch2 = bar.first.second;
      std::string label1 = Form("%s_%s",ch1.c_str(),stepLabel.c_str());
      std::string label2 = Form("%s_%s",ch2.c_str(),stepLabel.c_str());
      
      std::string ch3 = bar.second.first;
      std::string ch4 = bar.second.second;
      std::string label3 = Form("%s_%s",ch3.c_str(),stepLabel.c_str());
      std::string label4 = Form("%s_%s",ch4.c_str(),stepLabel.c_str());
      
      std::string label1234 = Form("%s+%s-%s-%s+%s",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str(),stepLabel.c_str());
      
      
      //--------------------------------------------------------
      
      
      c = new TCanvas(Form("c_deltaT_totCorr_%s",label1234.c_str()),Form("c_deltaT_totCorr_%s",label1234.c_str()));
      // gPad -> SetLogy();
      
      h1_deltaT_totCorr[label1234] -> SetTitle(Form(";tot-corrected #Deltat [ps];entries"));
      h1_deltaT_totCorr[label1234] -> SetMarkerColor(kBlue);
      h1_deltaT_totCorr[label1234] -> SetLineColor(kBlue);
      h1_deltaT_totCorr[label1234] -> Draw("");
      
      float fitXMin = CTRMeans[label1234] - 1.*CTRSigmas[label1234];
      float fitXMax = CTRMeans[label1234] + 1.*CTRSigmas[label1234];
      TF1* fitFunc = new TF1(Form("fitFunc_totCorr_%s",label1234.c_str()),"gaus",fitXMin,fitXMax);
      fitFunc -> SetParameters(1,h1_deltaT_totCorr[label1234]->GetMean(),h1_deltaT_totCorr[label1234]->GetRMS());
      h1_deltaT_totCorr[label1234] -> Fit(fitFunc,"QNRSL+");
      fitFunc -> SetLineColor(kBlue-1);
      fitFunc -> SetLineWidth(2);
      fitFunc -> Draw("same");
      
      FindSmallestInterval(vals,h1_deltaT_totCorr[label1234],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      
      latex = new TLatex(0.55,0.85,Form("#splitline{time walk corr. #sigma_{CTR}^{eff} = %.1f ps}{time walk corr. #sigma_{CTR}^{gaus} = %.1f ps}",effSigma,fitFunc->GetParameter(2)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlue);
      latex -> Draw("same");
      
      if( g_tRes_effSigma_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel] == NULL )
      {
        g_tRes_effSigma_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel] = new TGraphErrors();
        g_tRes_gaus_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel] = new TGraphErrors();
        g_tRes_totCorr_effSigma_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel] = new TGraphErrors();
        g_tRes_totCorr_gaus_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel] = new TGraphErrors();
      }
      if( g_tRes_effSigma_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel] == NULL )
      {          
        g_tRes_effSigma_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel] = new TGraphErrors();
        g_tRes_gaus_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel] = new TGraphErrors();
        g_tRes_totCorr_effSigma_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel] = new TGraphErrors();
        g_tRes_totCorr_gaus_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel] = new TGraphErrors();
      }
      
      g_tRes_effSigma_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel] -> SetPoint(g_tRes_effSigma_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel]->GetN(),63-th,effSigma/sqrt(2.));
      g_tRes_effSigma_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel] -> SetPointError(g_tRes_effSigma_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel]->GetN()-1,0.,5.);
      g_tRes_gaus_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel] -> SetPoint(g_tRes_gaus_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel]->GetN(),63-th,fitFunc->GetParameter(2)/sqrt(2.));
      g_tRes_gaus_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel] -> SetPointError(g_tRes_gaus_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel]->GetN()-1,0.,fitFunc->GetParError(2));
      
      g_tRes_effSigma_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel] -> SetPoint(g_tRes_effSigma_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel]->GetN(),Vov,effSigma/sqrt(2.));
      g_tRes_effSigma_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel] -> SetPointError(g_tRes_effSigma_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel]->GetN()-1,0.,5.);
      g_tRes_gaus_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel] -> SetPoint(g_tRes_gaus_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel]->GetN(),Vov,fitFunc->GetParameter(2)/sqrt(2.));
      g_tRes_gaus_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel] -> SetPointError(g_tRes_gaus_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel]->GetN()-1,0.,fitFunc->GetParError(2));
      
      
      h1_deltaT[label1234] -> SetMarkerColor(kRed);
      h1_deltaT[label1234] -> SetLineColor(kRed);
      h1_deltaT[label1234] -> Draw("same");
      
      fitXMin = CTRMeans[label1234] - 1.*CTRSigmas[label1234];
      fitXMax = CTRMeans[label1234] + 1.*CTRSigmas[label1234];
      fitFunc = new TF1(Form("fitFunc_%s",label1234.c_str()),"gaus",fitXMin,fitXMax);
      fitFunc -> SetParameters(1,h1_deltaT[label1234]->GetMean(),h1_deltaT[label1234]->GetRMS());
      h1_deltaT[label1234] -> Fit(fitFunc,"QNRSL+");
      fitFunc -> SetLineColor(kRed-1);
      fitFunc -> SetLineWidth(2);
      fitFunc -> Draw("same");
      
      FindSmallestInterval(vals,h1_deltaT[label1234],0.68);
      mean = vals[0];
      min = vals[4];
      max = vals[5];
      delta = max-min;
      sigma = 0.5*delta;
      effSigma = sigma;
      
      latex = new TLatex(0.55,0.65,Form("#splitline{raw #sigma_{CTR}^{eff} = %.1f ps}{raw #sigma_{CTR}^{gaus} = %.1f ps}",effSigma,fitFunc->GetParameter(2)));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");
      
      g_tRes_totCorr_effSigma_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel] -> SetPoint(g_tRes_totCorr_effSigma_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel]->GetN(),63-th,effSigma/sqrt(2.));
      g_tRes_totCorr_effSigma_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel] -> SetPointError(g_tRes_totCorr_effSigma_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel]->GetN()-1,0.,5.);
      g_tRes_totCorr_gaus_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel] -> SetPoint(g_tRes_totCorr_gaus_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel]->GetN(),63-th,fitFunc->GetParameter(2)/sqrt(2.));
      g_tRes_totCorr_gaus_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel] -> SetPointError(g_tRes_totCorr_gaus_vs_th[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+VovLabel]->GetN()-1,0.,fitFunc->GetParError(2));
      
      g_tRes_totCorr_effSigma_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel] -> SetPoint(g_tRes_totCorr_effSigma_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel]->GetN(),Vov,effSigma/sqrt(2.));
      g_tRes_totCorr_effSigma_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel] -> SetPointError(g_tRes_totCorr_effSigma_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel]->GetN()-1,0.,5.);
      g_tRes_totCorr_gaus_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel] -> SetPoint(g_tRes_totCorr_gaus_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel]->GetN(),Vov,fitFunc->GetParameter(2)/sqrt(2.));
      g_tRes_totCorr_gaus_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel] -> SetPointError(g_tRes_totCorr_gaus_vs_Vov[ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+thLabel]->GetN()-1,0.,fitFunc->GetParError(2));
      
      c -> Print(Form("%s/c_deltaT_totCorr_%s.png",plotDir.c_str(),label1234.c_str()));
      c -> Print(Form("%s/c_deltaT_totCorr_%s.pdf",plotDir.c_str(),label1234.c_str()));
    }
  }
  
  
  //--------------------------------------------------------
  
  
  for(auto pair : pairsVec)
  {
    std::string ch1 = pair.first;
    std::string ch2 = pair.second;
    
    c = new TCanvas(Form("c_tRes_vs_th_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_tRes_vs_th_%s-%s",ch1.c_str(),ch2.c_str()));
    // gPad -> SetLogy();
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,500.) );
    if( pairsMode == 1 )
      hPad -> SetTitle(";63 - threshold [DAC];#sigma_{t_{diff}} / #sqrt{2} [ps]");
    if( pairsMode == 2 )
      hPad -> SetTitle(";63 - threshold [DAC];#sigma_{t_{diff}} / 2 [ps]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    int iter = 0;
    for(auto mapIt : VovLabels)
    {
      std::string label = ch1+"-"+ch2+"_"+mapIt.first;
      TGraph* g_effSigma = g_tRes_effSigma_vs_th[label];
      TGraph* g_gaus = g_tRes_gaus_vs_th[label];
      TGraph* g_totCorr_effSigma = g_tRes_totCorr_effSigma_vs_th[label];
      TGraph* g_totCorr_gaus = g_tRes_totCorr_gaus_vs_th[label];
      
      g_effSigma -> SetLineColor(1+iter);
      g_effSigma -> SetMarkerColor(1+iter);
      g_effSigma -> SetMarkerStyle(20);
      // g_effSigma -> Draw("PL,same");
      
      g_gaus -> SetLineColor(1+iter);
      g_gaus -> SetMarkerColor(1+iter);
      g_gaus -> SetMarkerStyle(20);
      g_gaus -> Draw("PL,same");
      
      g_totCorr_effSigma -> SetLineColor(1+iter);
      g_totCorr_effSigma -> SetMarkerColor(1+iter);
      g_totCorr_effSigma -> SetMarkerStyle(21);
      // g_totCorr_effSigma -> Draw("PL,same");
      
      g_totCorr_gaus -> SetLineColor(1+iter);
      g_totCorr_gaus -> SetMarkerColor(1+iter);
      g_totCorr_gaus -> SetMarkerStyle(25);
      g_totCorr_gaus -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_tRes_vs_th_%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
    c -> Print(Form("%s/c_tRes_vs_th_%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
  }
  
  for(auto bar : barsVec)
  {
    std::string ch1 = bar.first.first;
    std::string ch2 = bar.first.second;
    std::string ch3 = bar.second.first;
    std::string ch4 = bar.second.second;
    
    c = new TCanvas(Form("c_tRes_vs_th_%s+%s-%s+%s",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str()),Form("c_tRes_vs_th_%s+%s-%s+%s",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str()));
    // gPad -> SetLogy();
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(-1.,0.,64.,500.) );
    hPad -> SetTitle(";63 - threshold [DAC];#sigma_{CTR} / #sqrt{2} [ps]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    int iter = 0;
    for(auto mapIt : VovLabels)
    {
      std::string label = ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+mapIt.first;
      TGraph* g_effSigma = g_tRes_effSigma_vs_th[label];
      TGraph* g_gaus = g_tRes_gaus_vs_th[label];
      TGraph* g_totCorr_effSigma = g_tRes_totCorr_effSigma_vs_th[label];
      TGraph* g_totCorr_gaus = g_tRes_totCorr_gaus_vs_th[label];
      
      g_effSigma -> SetLineColor(1+iter);
      g_effSigma -> SetMarkerColor(1+iter);
      g_effSigma -> SetMarkerStyle(20);
      // g_effSigma -> Draw("PL,same");
      
      g_gaus -> SetLineColor(1+iter);
      g_gaus -> SetMarkerColor(1+iter);
      g_gaus -> SetMarkerStyle(20);
      g_gaus -> Draw("PL,same");
      
      g_totCorr_effSigma -> SetLineColor(1+iter);
      g_totCorr_effSigma -> SetMarkerColor(1+iter);
      g_totCorr_effSigma -> SetMarkerStyle(21);
      // g_totCorr_effSigma -> Draw("PL,same");
      
      g_totCorr_gaus -> SetLineColor(1+iter);
      g_totCorr_gaus -> SetMarkerColor(1+iter);
      g_totCorr_gaus -> SetMarkerStyle(25);
      g_totCorr_gaus -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_tRes_vs_th_%s+%s_%s+%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str()));
    c -> Print(Form("%s/c_tRes_vs_th_%s+%s_%s+%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str()));
  }
  
  
  //--------------------------------------------------------
  
  
  for(auto pair : pairsVec)
  {
    std::string ch1 = pair.first;
    std::string ch2 = pair.second;
    
    c = new TCanvas(Form("c_tRes_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()),Form("c_tRes_vs_Vov_%s-%s",ch1.c_str(),ch2.c_str()));
    // gPad -> SetLogy();
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,500.) );
    hPad -> SetTitle(";V_{ov} [V];#sigma_{t_{diff}} / 2 [ps]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    int iter = 0;
    for(auto mapIt : thLabels)
    {
      std::string label = ch1+"-"+ch2+"_"+mapIt.first;
      TGraph* g_effSigma = g_tRes_effSigma_vs_Vov[label];
      TGraph* g_gaus = g_tRes_gaus_vs_Vov[label];
      TGraph* g_totCorr_effSigma = g_tRes_totCorr_effSigma_vs_Vov[label];
      TGraph* g_totCorr_gaus = g_tRes_totCorr_gaus_vs_Vov[label];
      
      g_effSigma -> SetLineColor(1+iter);
      g_effSigma -> SetMarkerColor(1+iter);
      g_effSigma -> SetMarkerStyle(20);
      // g_effSigma -> Draw("PL,same");
      
      g_gaus -> SetLineColor(1+iter);
      g_gaus -> SetMarkerColor(1+iter);
      g_gaus -> SetMarkerStyle(20);
      g_gaus -> Draw("PL,same");
      
      g_totCorr_effSigma -> SetLineColor(1+iter);
      g_totCorr_effSigma -> SetMarkerColor(1+iter);
      g_totCorr_effSigma -> SetMarkerStyle(21);
      // g_totCorr_effSigma -> Draw("PL,same");
      
      g_totCorr_gaus -> SetLineColor(1+iter);
      g_totCorr_gaus -> SetMarkerColor(1+iter);
      g_totCorr_gaus -> SetMarkerStyle(25);
      g_totCorr_gaus -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_tRes_vs_Vov_%s-%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
    c -> Print(Form("%s/c_tRes_vs_Vov_%s-%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str()));
  }
  
  for(auto bar : barsVec)
  {
    std::string ch1 = bar.first.first;
    std::string ch2 = bar.first.second;
    std::string ch3 = bar.second.first;
    std::string ch4 = bar.second.second;
    
    c = new TCanvas(Form("c_tRes_vs_Vov_%s+%s-%s+%s",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str()),Form("c_tRes_vs_Vov_%s+%s-%s+%s",ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str()));
    // gPad -> SetLogy();
    
    TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,10.,500.) );
    hPad -> SetTitle(";V_{ov} [V];#sigma_{CTR} / #sqrt{2} [ps]");
    hPad -> Draw();
    gPad -> SetGridy();
    
    int iter = 0;
    for(auto mapIt : thLabels)
    {
      std::string label = ch1+"+"+ch2+"-"+ch3+"+"+ch4+"_"+mapIt.first;
      TGraph* g_effSigma = g_tRes_effSigma_vs_Vov[label];
      TGraph* g_gaus = g_tRes_gaus_vs_Vov[label];
      TGraph* g_totCorr_effSigma = g_tRes_totCorr_effSigma_vs_Vov[label];
      TGraph* g_totCorr_gaus = g_tRes_totCorr_gaus_vs_Vov[label];
      
      g_effSigma -> SetLineColor(1+iter);
      g_effSigma -> SetMarkerColor(1+iter);
      g_effSigma -> SetMarkerStyle(20);
      // g_effSigma -> Draw("PL,same");
      
      g_gaus -> SetLineColor(1+iter);
      g_gaus -> SetMarkerColor(1+iter);
      g_gaus -> SetMarkerStyle(20);
      g_gaus -> Draw("PL,same");
      
      g_totCorr_effSigma -> SetLineColor(1+iter);
      g_totCorr_effSigma -> SetMarkerColor(1+iter);
      g_totCorr_effSigma -> SetMarkerStyle(21);
      // g_totCorr_effSigma -> Draw("PL,same");
      
      g_totCorr_gaus -> SetLineColor(1+iter);
      g_totCorr_gaus -> SetMarkerColor(1+iter);
      g_totCorr_gaus -> SetMarkerStyle(25);
      g_totCorr_gaus -> Draw("PL,same");
      
      latex = new TLatex(0.55,0.85-0.04*iter,Form("%s",mapIt.first.c_str()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack+iter);
      latex -> Draw("same");
      
      ++iter;
    }
    
    c -> Print(Form("%s/c_tRes_vs_Vov_%s+%s_%s+%s.png",plotDir.c_str(),ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str()));
    c -> Print(Form("%s/c_tRes_vs_Vov_%s+%s_%s+%s.pdf",plotDir.c_str(),ch1.c_str(),ch2.c_str(),ch3.c_str(),ch4.c_str()));
  }
  
  
}
