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
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"



int main(int argc, char** argv)
{
  setTDRStyle();
  
  if( argc < 2 )
  {
    std::cout << ">>> drawAlignment::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  
  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);

  
  
  //--- get parameters
  std::string label = opts.GetOpt<std::string>("Input.label");
  
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  plotDir += "/" + label;
  system(Form("mkdir -p %s",plotDir.c_str()));
  
  std::vector<std::string> channels = opts.GetOpt<std::vector<std::string> >("Channels.channels");

  int doTracking = opts.GetOpt<int>("Options.doTracking");
  int hodoPlane  = opts.GetOpt<int>("Options.hodoPlane");
  int nFibresMin = opts.GetOpt<int>("Options.nFibresMin");
  int nFibresMax = opts.GetOpt<int>("Options.nFibresMax");
  

  
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
  TreeVars tv;
  InitTreeVars(h4,tv,opts);
  
  int nEntries = h4->GetEntries();
  std::cout << ">>> Events read: " << nEntries << std::endl;  
  
  
  
  //------------------
  // define histograms
  
  TFile* outFile = TFile::Open(Form("%s/drawAlignment_%s.root",plotDir.c_str(),label.c_str()),"RECREATE");
  outFile -> cd();
  
  TH1F* h1_hodo1X  = new TH1F("h_hodo1X", "",128,-32.,32.);
  TH1F* h1_hodo1Y  = new TH1F("h_hodo1Y", "",128,-32.,32.);
  TH2F* h2_hodo1XY = new TH2F("h_hodo1XY","",128,-32.,32.,128,-32.,32.);
  
  TH1F* h1_hodo2X  = new TH1F("h_hodo2X", "",128,-32.,32.);
  TH1F* h1_hodo2Y  = new TH1F("h_hodo2Y", "",128,-32.,32.);
  TH2F* h2_hodo2XY = new TH2F("h_hodo2XY","",128,-32.,32.,128,-32.,32.);
  
  TH1F* h1_hodo12X = new TH1F("h_hodo12X","",256,-32.,32.);
  TH1F* h1_hodo12Y = new TH1F("h_hodo12Y","",256,-32.,32.);
  
  std::map<std::string,TH1F*> h1_amp;
  std::map<std::string,TProfile2D*> p2_eff_vs_XY;
  
  std::map<std::string,TH1F*> h1_eff_hodo1;
  std::map<std::string,TH1F*> h1_eff_hodo2;
  std::map<std::string,int> n_pass_hodo1;
  std::map<std::string,int> n_pass_hodo2;
  std::map<std::string,int> n_tot_hodo1;
  std::map<std::string,int> n_tot_hodo2;
  for(auto ch : channels)
  {
    std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
    
    h1_amp[ch] = new TH1F(Form("h_amp_%s",ch.c_str()),"",4000,0.,1.);
    p2_eff_vs_XY[ch] = new TProfile2D(Form("p2_eff_vs_XY_%s",ch.c_str()),"",128,-32.,32.,128,-32.,32.);
    
    h1_eff_hodo1[ch] = new TH1F(Form("h1_eff_hodo1_%s",ch.c_str()),"",1000,0.,1.);
    h1_eff_hodo2[ch] = new TH1F(Form("h1_eff_hodo2_%s",ch.c_str()),"",1000,0.,1.);
  }
  
  

  //--------------------------------
  // define outfile with good spills

  std::ofstream goodSpillList_new("goodSpills_new.txt",std::ios::out);
  
  std::ifstream goodSpillList("goodSpills.txt",std::ios::in);
  std::map<std::pair<int,int>,bool> goodSpills;
  std::string line;
  while(1)
  {
    getline(goodSpillList,line,'\n');
    if( !goodSpillList.good() ) break;
    if( line.at(0) == '#' ) continue;
    // std::cout << "Reading line " << line << std::endl;

    std::stringstream ss(line);
    int run, spill;
    ss >> run >> spill;
    goodSpills[std::make_pair(run,spill)] = true;
  }


  
  //-----------------------
  // 1st loop over events
  int oldSpill = -1;
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    h4 -> GetEntry(entry);
    // if( !goodSpills[std::make_pair(tv.run,tv.spill)] ) continue;
    
    
    if( debugMode )
    {
      std::cout << ">>>>>> entry: " << entry << "   event: " << tv.event << "   spill: " << tv.spill << "   run: " << tv.run << std::endl;
    }
    
    if( tv.hodoX[0] < -64. ) continue;
    if( tv.hodoY[0] < -64. ) continue;
    if( tv.hodoX[1] < -64. ) continue;
    if( tv.hodoY[1] < -64. ) continue;
    
    if( (tv.nFibresOnX[0] >= 2 && tv.nFibresOnX[0] <= 3) && (tv.nFibresOnY[0] >= 2 && tv.nFibresOnY[0] <= 3) )
    {
      h1_hodo1X -> Fill( tv.hodoX[0] );
      h1_hodo1Y -> Fill( tv.hodoY[0] );
      h2_hodo1XY -> Fill( tv.hodoX[0],tv.hodoY[0] );
    }

    if( (tv.nFibresOnX[1] >= 2 && tv.nFibresOnX[1] <= 3) && (tv.nFibresOnY[1] >= 2 && tv.nFibresOnY[1] <= 3) )
    {
      h1_hodo2X -> Fill( tv.hodoX[1] );
      h1_hodo2Y -> Fill( tv.hodoY[1] );
      h2_hodo2XY -> Fill( tv.hodoX[1],tv.hodoY[1] );
    }
    
    if( (tv.nFibresOnX[0] >= 2 && tv.nFibresOnX[0] <= 3) && (tv.nFibresOnY[0] >= 2 && tv.nFibresOnY[0] <= 3) &&
        (tv.nFibresOnX[1] >= 2 && tv.nFibresOnX[1] <= 3) && (tv.nFibresOnY[1] >= 2 && tv.nFibresOnY[1] <= 3) )
    {
      h1_hodo12X -> Fill( tv.hodoX[1]-tv.hodoX[0] );
      h1_hodo12Y -> Fill( tv.hodoY[1]-tv.hodoY[0] );
    }
    
    
    if( tv.spill != oldSpill )
    {
      bool goodSpill = false;
      for(auto ch: channels)
      {
        if( n_tot_hodo1[ch] == 0 ) continue;
        
        h1_eff_hodo1[ch] -> Fill( 1.*n_pass_hodo1[ch]/n_tot_hodo1[ch] );
        h1_eff_hodo2[ch] -> Fill( 1.*n_pass_hodo2[ch]/n_tot_hodo2[ch] );
        
        if( 1.*n_pass_hodo1[ch]/n_tot_hodo1[ch] > 0.85 )
          goodSpill = true;
      }
      
      for(auto ch: channels)
      {
        n_pass_hodo1[ch] = 0;
        n_tot_hodo1[ch] = 0;
        n_pass_hodo2[ch] = 0;
        n_tot_hodo2[ch] = 0;
      }

      if( goodSpill )
        goodSpillList_new << tv.run << " " << oldSpill << "\n";
      oldSpill = tv.spill;
    }
    
    
    for(auto ch: channels)
    {
      // reconstruct position
      ReconstructHodoPosition(tv,opts,ch,bool(doTracking),hodoPlane);
      
      
      // fill amplitude and time plots
      std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
      
      if( ampCh == "NULL" ) continue;
      float amp = tv.amp_max[tv.channelIds[ampCh]] / 4096.;
      h1_amp[ch] -> Fill( amp );
      
      if( amp < 0. || amp > 1. ) continue;
      
      if( (tv.nFibresOnX[0] < nFibresMin || tv.nFibresOnX[0] > nFibresMax) ) continue;
      if( (tv.nFibresOnX[1] < nFibresMin || tv.nFibresOnX[1] > nFibresMax) ) continue;
      if( (tv.nFibresOnY[0] < nFibresMin || tv.nFibresOnY[0] > nFibresMax) ) continue;
      if( (tv.nFibresOnY[1] < nFibresMin || tv.nFibresOnY[1] > nFibresMax) ) continue;
      
      float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str()));
      float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));
      
      float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str())); 
      float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
      float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
      float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
      
      if( amp > ampMin && amp < ampMax )
      {
        p2_eff_vs_XY[ch] -> Fill( tv.beamX,tv.beamY,1. );
        
        if( tv.beamX > xtalXMin && tv.beamX < xtalXMax && tv.beamY > xtalYMin && tv.beamY < xtalYMax )
        {
          n_tot_hodo1[ch] += 1;
          n_pass_hodo1[ch] += 1;
        }
        if( (tv.hodoX[1]-3.50) > xtalXMin && (tv.hodoX[1]-3.50) < xtalXMax && (tv.hodoY[1]-0.12) > xtalYMin && (tv.hodoY[1]-0.12) < xtalYMax )
        {
          n_tot_hodo2[ch] += 1;
          n_pass_hodo2[ch] += 1;
        }
      }
      else
      {
        p2_eff_vs_XY[ch] -> Fill( tv.beamX,tv.beamY,0. );
        
        if( tv.beamX > xtalXMin && tv.beamX < xtalXMax && tv.beamY > xtalYMin && tv.beamY < xtalYMax )
        {
          n_tot_hodo1[ch] += 1;
        }
        if( (tv.hodoX[1]-3.50) > xtalXMin && (tv.hodoX[1]-3.50) < xtalXMax && (tv.hodoY[1]-0.12) > xtalYMin && (tv.hodoY[1]-0.12) < xtalYMax )
        {
          n_tot_hodo2[ch] += 1;
        }        
      }
    }
  }
  std::cout << "\n>>> end 1st loop" << std::endl;
  
  
  
  //--- draw plots
  TCanvas* c;
  TH1F* histo;
  TH1F* histoCorr;
  TH2F* histo2;
  TProfile2D* prof2;


  
  c = new TCanvas(Form("hodo1XY"),Form("hodo1XY"));
  
  histo2 = h2_hodo1XY;
  histo2 -> SetTitle(Form(";hodo1 X;hodo1 Y"));
  histo2 -> Draw("COLZ");
  
  gPad -> Update();
  
  c -> Print(Form("%s/c_hodo1XY.png",plotDir.c_str()));
  c -> Print(Form("%s/c_hodo1XY.pdf",plotDir.c_str()));
  
  
  
  c = new TCanvas(Form("hodo2XY"),Form("hodo2XY"));
  
  histo2 = h2_hodo2XY;
  histo2 -> SetTitle(Form(";hodo2 X;hodo2 Y"));
  histo2 -> Draw("COLZ");
  
  gPad -> Update();
  
  c -> Print(Form("%s/c_hodo2XY.png",plotDir.c_str()));
  c -> Print(Form("%s/c_hodo2XY.pdf",plotDir.c_str()));  
  
  
  
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
  
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string label = opts.GetOpt<std::string>(Form("%s.label",ch.c_str()));
    
    c = new TCanvas(Form("c_amp_%s",label.c_str()),Form("c_amp_%s",label.c_str()));
    
    std::string ampCh = opts.GetOpt<std::string>(Form("%s.ampCh",ch.c_str()));
    float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
    float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));
    float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str())); 
    float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
    float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
    float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
    
    TLine* line_xtalXLow = new TLine(xtalXMin,xtalYMin,xtalXMax,xtalYMin);
    line_xtalXLow -> SetLineStyle(2);
    line_xtalXLow -> SetLineWidth(3);
    line_xtalXLow -> SetLineColor(kMagenta);
    TLine* line_xtalXHig = new TLine(xtalXMin,xtalYMax,xtalXMax,xtalYMax);
    line_xtalXHig -> SetLineStyle(2);
    line_xtalXHig -> SetLineWidth(3);
    line_xtalXHig -> SetLineColor(kMagenta);
    
    TLine* line_xtalYLow = new TLine(xtalXMin,xtalYMin,xtalXMin,xtalYMax);
    line_xtalYLow -> SetLineStyle(2);
    line_xtalYLow -> SetLineWidth(3);
    line_xtalYLow -> SetLineColor(kMagenta);
    TLine* line_xtalYHig = new TLine(xtalXMax,xtalYMin,xtalXMax,xtalYMax);
    line_xtalYHig -> SetLineStyle(2);
    line_xtalYHig -> SetLineWidth(3);
    line_xtalYHig -> SetLineColor(kRed);
    
    if( index < 0 ) continue;
    
    gPad -> SetLogy();
    histo = h1_amp[ch];
    histo -> Scale(1./histo->Integral());
    histo -> SetTitle(Form(";amplitude [mV];entries"));
    histo -> SetLineColor(kRed);
    histo -> SetFillColor(kRed);
    histo -> SetFillStyle(3003);
    histo -> GetXaxis() -> SetRangeUser(0.,1.1*ampMax);
    histo -> Draw("hist");
    
    TLine* line_min = new TLine(ampMin,histo->GetMinimum(),ampMin,histo->GetMaximum());
    line_min -> SetLineStyle(2);
    line_min -> Draw("same");
    TLine* line_max = new TLine(ampMax,histo->GetMinimum(),ampMax,histo->GetMaximum());
    line_max -> SetLineStyle(2);
    line_max -> Draw("same");
    
    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();
    
    c -> Print(Form("%s/c_amp_%s.png",plotDir.c_str(),label.c_str()));
    c -> Print(Form("%s/c_amp_%s.pdf",plotDir.c_str(),label.c_str()));
    
    
    c = new TCanvas(Form("c_eff_vs_XY_%s",label.c_str()),Form("c_eff_vs_XY_%s",label.c_str()));
    gPad -> SetLogz();
    
    prof2 = p2_eff_vs_XY[ch];
    prof2 -> SetTitle(Form(";beam X [mm];beam Y [mm]"));
    prof2 -> SetMinimum(0.001);
    prof2 -> SetMaximum(1.1);
    prof2 -> Draw("COLZ");
    
    latexLabels[ch] -> Draw("same");
    
    line_xtalXLow -> Draw("same");
    line_xtalXHig -> Draw("same");
    line_xtalYLow -> Draw("same");
    line_xtalYHig -> Draw("same");
    
    gPad -> Update();
    
    c -> Print(Form("%s/c_eff_vs_XY_%s.png",plotDir.c_str(),label.c_str()));
    c -> Print(Form("%s/c_eff_vs_XY_%s.pdf",plotDir.c_str(),label.c_str()));
  }
  
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;

  goodSpillList_new.close();
}
