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



struct Event
{
  std::map<std::string,float> amps;
  std::map<std::string,float> times;
  float x;
  float y;
};



int main(int argc, char** argv)
{
  setTDRStyle();

  
  if( argc < 2 )
  {
    std::cout << ">>> drawCTRSingles::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  
  TRandom3 r;
  
  
  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  

  
  //--- get parameters
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  system(Form("mkdir -p %s",plotDir.c_str()));
  
  std::vector<std::string> channels = opts.GetOpt<std::vector<std::string> >("Channels.channels");
  
  int doTracking = opts.GetOpt<int>("Options.doTracking");
  int hodoPlane  = opts.GetOpt<int>("Options.hodoPlane");
  int nFibresMin = opts.GetOpt<int>("Options.nFibresMin");
  int nFibresMax = opts.GetOpt<int>("Options.nFibresMax");
  
  float energySmearing = opts.GetOpt<float>("Options.energySmearing");
  
  
  //--- open input files
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  std::string runs = opts.GetOpt<std::string>("Input.runs");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");
  std::string treeName = opts.GetOpt<std::string>("Input.treeName");
  
  TChain* h4 = new TChain(treeName.c_str(),treeName.c_str());
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
      // std::string fileName = Form("%s/%d/*.root",inputDir.c_str(),run);
      std::string fileName = Form("%s/*%d*.root",inputDir.c_str(),run);
      std::cout << ">>> Adding flle " << fileName << std::endl;
      h4 -> Add(fileName.c_str());
    }
  }
  
  std::string goodRunsFile = opts.GetOpt<std::string>("Input.goodRunsFile");
  std::ifstream goodRunList(goodRunsFile.c_str(),std::ios::in);
  std::map<int,bool> goodRuns;
  std::string line;
  while(1)
  {
    getline(goodRunList,line,'\n');
    if( !goodRunList.good() ) break;
    if( line.at(0) == '#' ) continue;
    // std::cout << "Reading line " << line << std::endl;

    std::stringstream ss(line);
    int run, spill;
    ss >> run;
    goodRuns[run] = true;
  }
  

  
  //--- define tree branch addresses
  TreeVars tv;
  InitTreeVars(h4,tv,opts);
  //InitTreeVarsFNAL(h4,tv,opts);
  
  int nEntries = h4->GetEntries();
  std::cout << ">>> Events read: " << nEntries << std::endl;
  
  
  
  //------------------
  // define histograms

  std::string conf = opts.GetOpt<std::string>("Input.conf");
  TFile* outFile = TFile::Open(Form("%s/drawCTRSingles_%s.root",plotDir.c_str(),conf.c_str()),"RECREATE");
  outFile -> cd();

  TTree* outTree = new TTree("results","results");

  std::map<std::string,float> t_CTR_effSigma;
  std::map<std::string,float> t_CTR_gausSigma;
  std::map<std::string,float> t_CTR_gausSigmaErr;
  std::map<std::string,float> t_CTR_ampCorr_effSigma;
  std::map<std::string,float> t_CTR_ampCorr_gausSigma;
  std::map<std::string,float> t_CTR_ampCorr_gausSigmaErr;
  std::map<std::string,float> t_CTR_ampCorr_posCorr_effSigma;
  std::map<std::string,float> t_CTR_ampCorr_posCorr_gausSigma;
  std::map<std::string,float> t_CTR_ampCorr_posCorr_gausSigmaErr;
  
  for(auto ch : channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    outTree -> Branch(Form("CTR_effSigma_%s",ch.c_str()),&t_CTR_effSigma[ch]);
    outTree -> Branch(Form("CTR_gausSigma_%s",ch.c_str()),&t_CTR_gausSigma[ch]);
    outTree -> Branch(Form("CTR_gausSigmaErr_%s",ch.c_str()),&t_CTR_gausSigmaErr[ch]);
    outTree -> Branch(Form("CTR_ampCorr_effSigma_%s",ch.c_str()),&t_CTR_ampCorr_effSigma[ch]);
    outTree -> Branch(Form("CTR_ampCorr_gausSigma_%s",ch.c_str()),&t_CTR_ampCorr_gausSigma[ch]);
    outTree -> Branch(Form("CTR_ampCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_ampCorr_gausSigmaErr[ch]);
    outTree -> Branch(Form("CTR_ampCorr_posCorr_effSigma_%s",ch.c_str()),&t_CTR_ampCorr_posCorr_effSigma[ch]);
    outTree -> Branch(Form("CTR_ampCorr_posCorr_gausSigma_%s",ch.c_str()),&t_CTR_ampCorr_posCorr_gausSigma[ch]);
    outTree -> Branch(Form("CTR_ampCorr_posCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_ampCorr_posCorr_gausSigmaErr[ch]);
  }
  
  std::map<std::string,TProfile*> p_eff_vs_X;
  std::map<std::string,TProfile*> p_eff_vs_Y;
  std::map<std::string,TProfile2D*> p2_eff_vs_XY;
  
  std::map<std::string,TH1F*> h_amp;
  std::map<std::string,TH1F*> h_amp_cut;
  std::map<std::string,TProfile*> p_amp_vs_X;
  std::map<std::string,TProfile*> p_amp_vs_Y;
  std::map<std::string,TProfile2D*> p2_amp_vs_XY;

  std::map<std::string,TH1F*> h_b_rms;  
  std::map<std::string,TH1F*> h_time;
  
  std::map<std::string,TProfile*> p_time_vs_amp;
  std::map<std::string,TProfile*> p_time_ampCorr_vs_amp;
  
  std::map<std::string,TProfile*> p_time_ampCorr_vs_X;
  std::map<std::string,TProfile*> p_time_ampCorr_vs_Y;
  std::map<std::string,TProfile2D*> p2_time_ampCorr_vs_XY;
  std::map<std::string,TProfile*> p_time_ampCorr_posCorr_vs_X;
  std::map<std::string,TProfile*> p_time_ampCorr_posCorr_vs_Y;
  std::map<std::string,TProfile2D*> p2_time_ampCorr_posCorr_vs_XY;
  
  std::map<std::string,TH1F*> h_occupancyX;
  std::map<std::string,TH1F*> h_occupancyY;  
  std::map<std::string,TH2F*> h2_occupancyXY;
  
  std::map<std::string,TH1F*> h_chessOccupancyX;
  std::map<std::string,TH1F*> h_chessOccupancyY;
  std::map<std::string,TH2F*> h2_chessOccupancyXY;
  std::map<std::string,TProfile*> map_p_time_vs_amp;
  
  std::map<std::string,TH1F*> h_CTR_raw;
  
  for(auto ch : channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str()));
    float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
    float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
    float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
    
    std::vector<std::string> ampChs = opts.GetOpt<std::vector<std::string> >(Form("%s.ampCh", ch.c_str()));
    for( auto ampCh : ampChs )
    {
      p_eff_vs_X[ch+"_"+ampCh] = new TProfile(Form("p_eff_vs_X_%s",(ch+"_"+ampCh).c_str()),"",100,-10.,40.);
      p_eff_vs_Y[ch+"_"+ampCh] = new TProfile(Form("p_eff_vs_Y_%s",(ch+"_"+ampCh).c_str()),"",100,0.,50.);
      p2_eff_vs_XY[ch+"_"+ampCh] = new TProfile2D(Form("p2_eff_vs_XY_%s",(ch+"_"+ampCh).c_str()),"",100,-10.,40.,100,0.,50.);
      
      h_amp[ch+"_"+ampCh]       = new TH1F(Form("h_amp_%s",      (ch+"_"+ampCh).c_str()),"",1000,0.,1.);
      h_amp_cut[ch+"_"+ampCh]   = new TH1F(Form("h_amp_cut_%s",  (ch+"_"+ampCh).c_str()),"",1000,0.,1.);
      p_amp_vs_X[ch+"_"+ampCh] = new TProfile(Form("p_amp_vs_X_%s",(ch+"_"+ampCh).c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
      p_amp_vs_Y[ch+"_"+ampCh] = new TProfile(Form("p_amp_vs_Y_%s",(ch+"_"+ampCh).c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
      p2_amp_vs_XY[ch+"_"+ampCh] = new TProfile2D(Form("p2_amp_vs_XY_%s",(ch+"_"+ampCh).c_str()),"",int(1*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax,int(1*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    }
    
    std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
    std::vector<std::string> timeChs = opts.GetOpt<std::vector<std::string> >(Form("%s.timeCh",ch.c_str()));
    int timeIt = 0;
    for( auto timeCh : timeChs )
    {
      float ampMin = opts.GetOpt<std::vector<float> >(Form("%s.ampMin",ch.c_str())).at(timeIt);
      float ampMax = opts.GetOpt<std::vector<float> >(Form("%s.ampMax",ch.c_str())).at(timeIt);
      
      std::cout << "creating histogram " << "h_b_rms[" << ch << "+_+" << timeCh << "]" << std::endl;
      h_b_rms[ch+"_"+timeCh] = new TH1F(Form("h_b_rms_%s",(ch+"_"+timeCh).c_str()),"",100000,0.,0.5);
      
      for( auto timeMethod : timeMethods )
      {
        std::string label = ch+"_"+timeCh+"_"+timeMethod;
        
        h_time[label] = new TH1F(Form("h_time_%s",(label).c_str()),"",200,0.,200.);
      
        p_time_vs_amp[label] = new TProfile(Form("p_time_vs_amp_%s",(label).c_str()),"",400,ampMin,ampMax);
        p_time_ampCorr_vs_amp[label] = new TProfile(Form("p_time_ampCorr_vs_amp_%s",(label).c_str()),"",400,ampMin,ampMax);
        
        p_time_ampCorr_vs_X[label] = new TProfile(Form("p_time_ampCorr_vs_X_%s",label.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
        p_time_ampCorr_vs_Y[label] = new TProfile(Form("p_time_ampCorr_vs_Y_%s",label.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
        p2_time_ampCorr_vs_XY[label] = new TProfile2D(Form("p2_time_ampCorr_vs_XY_%s",label.c_str()),"",int(1*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax,int(1*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
      }
      ++timeIt;
    }
    for( auto timeMethod : timeMethods )
    {
      std::string label = ch+"_"+timeMethod+"__avg-ref";
      p_time_ampCorr_vs_X[label] = new TProfile(Form("p_time_ampCorr_vs_X_%s",label.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
      p_time_ampCorr_vs_Y[label] = new TProfile(Form("p_time_ampCorr_vs_Y_%s",label.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
      p2_time_ampCorr_vs_XY[label] = new TProfile2D(Form("p2_time_ampCorr_vs_XY_%s",label.c_str()),"",int(1*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax,int(1*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
      p_time_ampCorr_posCorr_vs_X[label] = new TProfile(Form("p_time_ampCorr_posCorr_vs_X_%s",label.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
      p_time_ampCorr_posCorr_vs_Y[label] = new TProfile(Form("p_time_ampCorr_posCorr_vs_Y_%s",label.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
      p2_time_ampCorr_posCorr_vs_XY[label] = new TProfile2D(Form("p2_time_ampCorr_posCorr_vs_XY_%s",label.c_str()),"",int(1*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax,int(1*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
      
      label = ch+"_"+timeMethod+"__diff-ref";
      p_time_ampCorr_vs_X[label] = new TProfile(Form("p_time_ampCorr_vs_X_%s",label.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
      p_time_ampCorr_vs_Y[label] = new TProfile(Form("p_time_ampCorr_vs_Y_%s",label.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
      p2_time_ampCorr_vs_XY[label] = new TProfile2D(Form("p2_time_ampCorr_vs_XY_%s",label.c_str()),"",int(1*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax,int(1*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
      p_time_ampCorr_posCorr_vs_X[label] = new TProfile(Form("p_time_ampCorr_posCorr_vs_X_%s",label.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
      p_time_ampCorr_posCorr_vs_Y[label] = new TProfile(Form("p_time_ampCorr_posCorr_vs_Y_%s",label.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
      p2_time_ampCorr_posCorr_vs_XY[label] = new TProfile2D(Form("p2_time_ampCorr_posCorr_vs_XY_%s",label.c_str()),"",int(1*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax,int(1*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    }
  
    for( auto timeMethod : timeMethods )
    {
      if( index >= 0 )
      {
        h_CTR_raw[ch+"_"+timeMethod+"__1-ref"] = new TH1F(Form("h_CTR_raw_%s_%s",ch.c_str(),(+"_"+timeMethod+"__1-ref").c_str()),"",20000,-200.,200.);
        h_CTR_raw[ch+"_"+timeMethod+"__2-ref"] = new TH1F(Form("h_CTR_raw_%s_%s",ch.c_str(),(+"_"+timeMethod+"__2-ref").c_str()),"",20000,-200.,200.);
        h_CTR_raw[ch+"_"+timeMethod+"__avg-ref"] = new TH1F(Form("h_CTR_raw_%s_%s",ch.c_str(),(+"_"+timeMethod+"__avg-ref").c_str()),"",20000,-200.,200.);
        h_CTR_raw[ch+"_"+timeMethod+"__diff"] = new TH1F(Form("h_CTR_raw_%s_%s",ch.c_str(),(+"_"+timeMethod+"__diff").c_str()),"",20000,-200.,200.);
      }
    }
    
    h_occupancyX[ch]   = new TH1F(Form("h_occupancyX_%s",ch.c_str()),  "",1.*(xtalXMax-xtalXMin),xtalXMin,xtalXMax);
    h_occupancyY[ch]   = new TH1F(Form("h_occupancyY_%s",ch.c_str()),  "",1.*(xtalYMax-xtalYMin),xtalYMin,xtalYMax);    
    h2_occupancyXY[ch] = new TH2F(Form("h2_occupancyXY_%s",ch.c_str()),"",1.*(xtalXMax-xtalXMin),xtalXMin,xtalXMax,1.*(xtalYMax-xtalYMin),xtalYMin,xtalYMax);
          
    int nBinsX = opts.GetOpt<int>(Form("%s.nBinsX",ch.c_str()));
    int nBinsY = opts.GetOpt<int>(Form("%s.nBinsY",ch.c_str()));

    if( nBinsX > 0 && nBinsY > 0 )
    {
      h_chessOccupancyX[ch]   = new TH1F(Form("h_chessOccupancyX_%s",ch.c_str()),   "",nBinsX,xtalXMin,xtalXMax);
      h_chessOccupancyY[ch]   = new TH1F(Form("h_chessOccupancyY_%s",ch.c_str()),  "",nBinsY,xtalYMin,xtalYMax);
      h2_chessOccupancyXY[ch] = new TH2F(Form("h2_chessOccupancyXY_%s",ch.c_str()),"",nBinsX,xtalXMin,xtalXMax,nBinsY,xtalYMin,xtalYMax);  
      
      for(int binX = 0; binX < nBinsX; ++binX)
        for(int binY = 0; binY < nBinsY; ++binY)
        {
          std::string label(Form("%s_%d-%d",ch.c_str(),binX,binY));
        map_p_time_vs_amp[label] = new TProfile(Form("p_time_vs_amp_%s",label.c_str()),"",50,0.,1.);
        }
    }
  }
  
  TCanvas* c;
  TH1F* histo;
  TH1F* histoCorr;
  TH1F* histoCorr2;
  TProfile* prof;
  
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
  
  
  
  std::vector<Event> events;
  
  //-----------------------
  // 1st loop over events
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    h4 -> GetEntry(entry);
    if( !goodRuns[tv.run] ) continue;
    
    
    std::map<std::string,float> amps;
    std::map<std::string,float> times;
    
    for(auto ch: channels)
    {
      // reconstruct position
      ReconstructHodoPosition(tv,opts,ch,bool(doTracking),hodoPlane);
      
      
      if( doTracking )
      {
        if( (tv.nFibresOnX[0] < nFibresMin || tv.nFibresOnX[0] > nFibresMax) ) continue;
        if( (tv.nFibresOnX[1] < nFibresMin || tv.nFibresOnX[1] > nFibresMax) ) continue;
        if( (tv.nFibresOnY[0] < nFibresMin || tv.nFibresOnY[0] > nFibresMax) ) continue;
        if( (tv.nFibresOnY[1] < nFibresMin || tv.nFibresOnY[1] > nFibresMax) ) continue;
      }
      else
      {
        if( tv.nTracks != 1 ) continue;
      }
      
      
      // fill amplitude and time plots
      std::vector<std::string> ampChs = opts.GetOpt<std::vector<std::string> >(Form("%s.ampCh", ch.c_str()));
      std::vector<std::string> timeChs = opts.GetOpt<std::vector<std::string> >(Form("%s.timeCh", ch.c_str()));
      std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
      float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str())); 
      float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
      float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
      float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
      
      int ampIt = -1;
      for( auto ampCh : ampChs )
      {
        ++ampIt;
        
        float amp = tv.amp_max[tv.channelIds[ampCh]] / 4096.;
        h_amp[ch+"_"+ampCh] -> Fill( amp );
        
        float ampMin = opts.GetOpt<std::vector<float> >(Form("%s.ampMin",ch.c_str())).at(ampIt);
        float ampMax = opts.GetOpt<std::vector<float> >(Form("%s.ampMax",ch.c_str())).at(ampIt);
        
        if( amp > ampMin && amp < ampMax)
        {
          if( tv.beamY > xtalYMin && tv.beamY < xtalYMax ) p_eff_vs_X[ch+"_"+ampCh] -> Fill( tv.beamX,1. );
          if( tv.beamX > xtalXMin && tv.beamX < xtalXMax ) p_eff_vs_Y[ch+"_"+ampCh] -> Fill( tv.beamY,1. );
          p2_eff_vs_XY[ch+"_"+ampCh] -> Fill( tv.beamX,tv.beamY,1. );
          
          if( tv.beamY > xtalYMin && tv.beamY < xtalYMax ) p_amp_vs_X[ch+"_"+ampCh] -> Fill( tv.beamX,amp );
          if( tv.beamX > xtalXMin && tv.beamX < xtalXMax ) p_amp_vs_Y[ch+"_"+ampCh] -> Fill( tv.beamY,amp );
          p2_amp_vs_XY[ch+"_"+ampCh] -> Fill( tv.beamX,tv.beamY,amp );
        }
        else
        {
          if( tv.beamY > xtalYMin && tv.beamY < xtalYMax ) p_eff_vs_X[ch+"_"+ampCh] -> Fill( tv.beamX,0. );
          if( tv.beamX > xtalXMin && tv.beamX < xtalXMax ) p_eff_vs_Y[ch+"_"+ampCh] -> Fill( tv.beamY,0. );
          p2_eff_vs_XY[ch+"_"+ampCh] -> Fill( tv.beamX,tv.beamY,0. );
        }
        
        if( tv.beamX < xtalXMin || tv.beamX > xtalXMax ) continue;
        if( tv.beamY < xtalYMin || tv.beamY > xtalYMax ) continue;
        
        if( amp < ampMin || amp > ampMax ) continue;
        if( isnan(amp) ) continue;
        
        h_amp_cut[ch+"_"+ampCh] -> Fill( amp );
        amps[ch+"_"+ampCh] = amp;
        
        std::string timeCh = timeChs.at(ampIt);
        float timeMin = opts.GetOpt<float>(Form("%s.timeMin",ch.c_str()));
        float timeMax = opts.GetOpt<float>(Form("%s.timeMax",ch.c_str()));
        
        float b_rms = tv.b_rms[tv.channelIds[timeCh]] / 4096.;
        h_b_rms[ch+"_"+timeCh] -> Fill( b_rms );
        if( b_rms > 0.002 ) continue;
        
        for( auto timeMethod : timeMethods )
        {          
          float tim = tv.time[tv.channelIds[timeCh]+tv.timeMethodIds[timeMethod]];
          h_time[ch+"_"+timeCh+"_"+timeMethod] -> Fill( tim );
          
          if( isinf(tim) ) continue;
          if( isnan(tim) ) continue;
          if( tim < timeMin || tim > timeMax ) continue;
          
          times[ch+"_"+timeCh+"_"+timeMethod] = tim;
        }
      }
      
    } // loop over channels
    
    Event aEvent;
    aEvent.amps = amps;
    aEvent.times = times;
    aEvent.x = tv.beamX;
    aEvent.y = tv.beamY;
    events.push_back( aEvent );
    
    
    for(auto ch: channels)
    {
      int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
      std::vector<std::string> ampChs = opts.GetOpt<std::vector<std::string> >(Form("%s.ampCh", ch.c_str()));
      std::vector<std::string> timeChs = opts.GetOpt<std::vector<std::string> >(Form("%s.timeCh", ch.c_str()));
      std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
      
      if( index < 0 ) continue;
      
      std::string time1Ch = timeChs.at(0);
      std::string time2Ch = timeChs.at(1);
      
      std::string refCh = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
      std::string timeChRef = opts.GetOpt<std::string>(Form("%s.timeCh",refCh.c_str()));
      std::string timeMethodRef = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",refCh.c_str())).at(0);
      
      for( auto timeMethod : timeMethods )
      {
        std::string label1 = ch+"_"+time1Ch+"_"+timeMethod;
        std::string label2 = ch+"_"+time2Ch+"_"+timeMethod;
        std::string labelRef = refCh+"_"+timeChRef+"_"+timeMethodRef;
        
        if( (times.count(label1) > 0) && (times.count(labelRef) > 0) )
          h_CTR_raw[ch+"_"+timeMethod+"__1-ref"] -> Fill( times[label1]-times[labelRef] );
        
        if( (times.count(label2) > 0) && (times.count(labelRef) > 0) )
          h_CTR_raw[ch+"_"+timeMethod+"__2-ref"] -> Fill( times[label2]-times[labelRef] );
        
        if( (times.count(label1) > 0) && (times.count(label2) > 0) && (times.count(labelRef) > 0) )
          h_CTR_raw[ch+"_"+timeMethod+"__avg-ref"] -> Fill( 0.5*(times[label1]+times[label2])-times[labelRef] );
        
        if( (times.count(label1) > 0) && (times.count(label2) > 0) )
          h_CTR_raw[ch+"_"+timeMethod+"__diff"] -> Fill( 0.5*(times[label2]-times[label1]) );
      }
    }
  }
  std::cout << "\n>>> end 1st loop" << std::endl;
  
  
  
  //--- draw 1st plots
  for(auto ch: channels)
  {
    std::vector<std::string> ampChs = opts.GetOpt<std::vector<std::string> >(Form("%s.ampCh", ch.c_str()));
    std::vector<std::string> timeChs = opts.GetOpt<std::vector<std::string> >(Form("%s.timeCh", ch.c_str()));
    std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
    int ampIt = 0;
    for( auto ampCh : ampChs )
    {    
      float ampMin = opts.GetOpt<std::vector<float> >(Form("%s.ampMin",ch.c_str())).at(ampIt);
      float ampMax = opts.GetOpt<std::vector<float> >(Form("%s.ampMax",ch.c_str())).at(ampIt);
      float timeMin = opts.GetOpt<float>(Form("%s.timeMin",ch.c_str())); 
      float timeMax = opts.GetOpt<float>(Form("%s.timeMax",ch.c_str()));
      float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str())); 
      float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
      float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
      float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
      float ampLow = h_amp_cut[ch+"_"+ampCh]->GetMean()-1.*h_amp_cut[ch+"_"+ampCh]->GetRMS();
      float ampHig = h_amp_cut[ch+"_"+ampCh]->GetMean()+1.*h_amp_cut[ch+"_"+ampCh]->GetRMS();
      
      TLine* line_min_amp = new TLine(ampMin,h_amp[ch+"_"+ampCh]->GetMinimum(),ampMin,h_amp[ch+"_"+ampCh]->GetMaximum());
      TLine* line_max_amp = new TLine(ampMax,h_amp[ch+"_"+ampCh]->GetMinimum(),ampMax,h_amp[ch+"_"+ampCh]->GetMaximum());
      std::vector<TLine*> lines_amp;
      lines_amp.push_back(line_min_amp);
      lines_amp.push_back(line_max_amp);
      
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
      
      std::vector<TLine*> lines_xtal;
      lines_xtal.push_back(line_xtalXLow);
      lines_xtal.push_back(line_xtalXHig);
      lines_xtal.push_back(line_xtalYLow);
      lines_xtal.push_back(line_xtalYHig);
      
            
      c = new TCanvas();
      std::vector<TH1F*> h_amps;
      h_amps.push_back(h_amp[ch+"_"+ampCh]);
      h_amps.push_back(h_amp_cut[ch+"_"+ampCh]);
      DrawHistogram(opts,ch,h_amps, ";amplitude [V];entries", 0.,1.1*ampMax,4,true,latexLabels[ch],&lines_amp,true,ampMin,0.6); 
      PrintCanvas(c,opts,ch,plotDir,Form("1_amp_%d",ampIt));
      delete c;
      
      std::string timeCh = timeChs.at(ampIt);
      
      c = new TCanvas();
      DrawHistogram(opts,ch,h_b_rms[ch+"_"+timeCh], ";baseline RMS [V];entries", 0.00001,0.1,1,true,latexLabels[ch]);
      gPad -> SetLogx();
      PrintCanvas(c,opts,ch,plotDir,Form("1_b_rms_%d_%s",ampIt,timeCh.c_str()));
      delete c;
      
      for( auto timeMethod : timeMethods )
      {
        TLine* line_min_time = new TLine(timeMin,h_time[ch+"_"+timeCh+"_"+timeMethod]->GetMinimum(),timeMin,h_time[ch+"_"+timeCh+"_"+timeMethod]->GetMaximum());
        TLine* line_max_time = new TLine(timeMax,h_time[ch+"_"+timeCh+"_"+timeMethod]->GetMinimum(),timeMax,h_time[ch+"_"+timeCh+"_"+timeMethod]->GetMaximum());
        std::vector<TLine*> lines_time;
        lines_time.push_back(line_min_time);
        lines_time.push_back(line_max_time);
        
        c = new TCanvas();
        DrawHistogram(opts,ch,h_time[ch+"_"+timeCh+"_"+timeMethod], ";time [ns];entries", 0.,200., 1,true,latexLabels[ch],&lines_time);
        PrintCanvas(c,opts,ch,plotDir,Form("1_time_%d_%s",ampIt,timeMethod.c_str()));
        delete c;
      }
      
      c = new TCanvas("c","c",2000,600);
      c -> Divide(3,1);
      c->cd(1); DrawProfile  (opts,ch,p_eff_vs_X[ch+"_"+ampCh],  ";hodoscope x [mm];efficiency",                 -10.,40.,0.001,1.1,kBlack,"",latexLabels[ch]);
      c->cd(2); DrawProfile  (opts,ch,p_eff_vs_Y[ch+"_"+ampCh],  ";hodoscope y [mm];efficiency",                   0.,50.,0.001,1.1,kBlack,"",latexLabels[ch]);
      c->cd(3); DrawProfile2D(opts,ch,p2_eff_vs_XY[ch+"_"+ampCh],";hodoscope x [mm];hodoscope y [mm];efficiency",-10.,40.,0.,50.,0.001,1.1,latexLabels[ch],&lines_xtal);
      PrintCanvas(c,opts,ch,plotDir,Form("1_eff_vs_XY_%d",ampIt)); delete c;
      
      c = new TCanvas("c","c",2000,600);
      c -> Divide(3,1);
      c->cd(1); DrawProfile  (opts,ch,p_amp_vs_X[ch+"_"+ampCh],  ";hodoscope x [mm];amplitude [V]",                 xtalXMin,xtalXMax,ampLow,ampHig,kBlack,"",latexLabels[ch]);
      c->cd(2); DrawProfile  (opts,ch,p_amp_vs_Y[ch+"_"+ampCh],  ";hodoscope y [mm];amplitude [V]",                 xtalYMin,xtalYMax,ampLow,ampHig,kBlack,"",latexLabels[ch]);
      c->cd(3); DrawProfile2D(opts,ch,p2_amp_vs_XY[ch+"_"+ampCh],";hodoscope x [mm];hodoscope y [mm];amplitude [V]",xtalXMin,xtalXMax,xtalYMin,xtalYMax,ampLow,ampHig,latexLabels[ch]);    
      PrintCanvas(c,opts,ch,plotDir,Form("1_amp_vs_XY_%d",ampIt)); delete c;
      
      ++ampIt;
    }
  }
  
  
  
  //--- find CTR ranges
  std::map<std::string,float> CTRRanges_mean;
  std::map<std::string,float> CTRRanges_sigma;
  for(auto mapIt : h_CTR_raw)
  {
    histo = mapIt.second;
    float* vals = new float[6];
    FindSmallestInterval(vals,histo,0.68);
    
    float mean = vals[0];
    float min = vals[4];
    float max = vals[5];
    float delta = max-min;
    float sigma = 0.5*delta;
    CTRRanges_mean[mapIt.first] = mean;
    CTRRanges_sigma[mapIt.first] = sigma;
    std::cout << ">>> " << mapIt.first << "   CTR mean: " << CTRRanges_mean[mapIt.first] << "   CTR sigma: " << CTRRanges_sigma[mapIt.first] << std::endl;
    
    outFile -> cd();
    histo -> Write();
  }
  
  
  
  
  
  
  //------------------------
  //--- 2nd loop over events
  int entry = 0;
  nEntries = events.size();
  for(auto event: events)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 2nd loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;    
    ++entry;
    
    std::map<std::string,float> amps = event.amps;
    std::map<std::string,float> times = event.times;
    
    
    for(auto ch: channels)
    {
      int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
      if( index < 0 ) continue;
      
      std::vector<std::string> ampChs = opts.GetOpt<std::vector<std::string> >(Form("%s.ampCh", ch.c_str()));
      std::vector<std::string> timeChs = opts.GetOpt<std::vector<std::string> >(Form("%s.timeCh",ch.c_str()));      
      std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
      
      std::string refCh = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
      std::string timeChRef = opts.GetOpt<std::string>(Form("%s.timeCh",refCh.c_str()));
      std::string timeMethodRef = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",refCh.c_str())).at(0);
      
      int timeIt = 0;
      for(auto timeCh : timeChs)
      {
        std::string ampCh = ampChs.at(timeIt);
        for(auto timeMethod : timeMethods)
        {
          std::string label = ch+"_"+timeCh+"_"+timeMethod;
          std::string labelRef = refCh+"_"+timeChRef+"_"+timeMethodRef;
          
          float timeLow = CTRRanges_mean[ch+"_"+timeMethod+Form("__%d-ref",timeIt+1)]-2.*CTRRanges_sigma[ch+"_"+timeMethod+Form("__%d-ref",timeIt+1)];
          float timeHig = CTRRanges_mean[ch+"_"+timeMethod+Form("__%d-ref",timeIt+1)]+2.*CTRRanges_sigma[ch+"_"+timeMethod+Form("__%d-ref",timeIt+1)];
          
          if( (times[label]-times[labelRef]) < timeLow ||
              (times[label]-times[labelRef]) > timeHig )
            continue;
          
          if( (times.count(label) > 0) && (times.count(labelRef) > 0) )
            p_time_vs_amp[label] -> Fill( amps[ch+"_"+ampCh],times[label]-times[labelRef] );
        }
        
        ++timeIt;
      }
    }
  }
  std::cout << "\n>>> end 2nd loop" << std::endl;
  
  
  
  //--- draw 2nd plots
  std::map<std::string,TF1*> fitFunc_corrAmp;
  
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    std::vector<std::string> ampChs = opts.GetOpt<std::vector<std::string> >(Form("%s.ampCh",ch.c_str()));      
    std::vector<std::string> timeChs = opts.GetOpt<std::vector<std::string> >(Form("%s.timeCh",ch.c_str()));      
    std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
    
    int timeIt = 0;
    for(auto timeCh : timeChs)
    {
      std::string ampCh = ampChs.at(timeIt);
      float ampMin = opts.GetOpt<std::vector<float> >(Form("%s.ampMin",ch.c_str())).at(timeIt);
      float ampMax = opts.GetOpt<std::vector<float> >(Form("%s.ampMax",ch.c_str())).at(timeIt);
      
      for(auto timeMethod : timeMethods)
      {
        std::string label = ch+"_"+timeCh+"_"+timeMethod;
        prof = p_time_vs_amp[label];
        
        float timeLow = CTRRanges_mean[ch+"_"+timeMethod+Form("__%d-ref",timeIt+1)]-4.*CTRRanges_sigma[ch+"_"+timeMethod+Form("__%d-ref",timeIt+1)];
        float timeHig = CTRRanges_mean[ch+"_"+timeMethod+Form("__%d-ref",timeIt+1)]+3.*CTRRanges_sigma[ch+"_"+timeMethod+Form("__%d-ref",timeIt+1)];
        
        std::string fitFunc = opts.GetOpt<std::string>(Form("%s.fitFunc",ch.c_str()));
        float fitXMin = opts.GetOpt<float>(Form("%s.fitFunc",ch.c_str()),1);
        float fitXMax = opts.GetOpt<float>(Form("%s.fitFunc",ch.c_str()),2);
        int nPar = opts.OptExist(Form("%s.fitFunc",ch.c_str()),3) ? opts.GetOpt<int>(Form("%s.fitFunc",ch.c_str()),3) : 0;
        
        fitFunc_corrAmp[label] = new TF1(Form("fitFunc_corrAmp_%s",label.c_str()),fitFunc.c_str(),ampMin,ampMax);
        for(int iPar = 0; iPar < nPar; ++iPar)
          fitFunc_corrAmp[label] -> SetParameter(iPar,opts.GetOpt<float>(Form("%s.fitFunc",ch.c_str()),4+iPar));
        
        prof -> Fit(fitFunc_corrAmp[label],"QNRS+","",fitXMin,fitXMax);
        
        c = new TCanvas();
        DrawProfile(opts,ch,p_time_vs_amp[label],";amplitude [V];#Delta_{t} [ns]",ampMin,ampMax,timeLow,timeHig,kBlack,"",latexLabels[ch],fitFunc_corrAmp[label],h_amp_cut[ch+"_"+ampCh]->GetMean());
        PrintCanvas(c,opts,ch,plotDir,Form("2_time_vs_amp_%d_%s",timeIt,timeMethod.c_str()));
        delete c;
      }
      
      ++timeIt;
    }
  }
  
  
  
  
  
  
  //------------------------
  //--- 3rd loop over events
  std::map<std::string,TH1F*> h_CTR;
  std::map<std::string,TH1F*> h_CTR_ampCorr;
  for(auto mapIt : h_CTR_raw)
  {
    h_CTR[mapIt.first] = new TH1F(Form("h_CTR_%s",mapIt.first.c_str()),"",1000,CTRRanges_mean[mapIt.first]-5.*CTRRanges_sigma[mapIt.first],CTRRanges_mean[mapIt.first]+5.*CTRRanges_sigma[mapIt.first]);
    h_CTR_ampCorr[mapIt.first] = new TH1F(Form("h_CTR_ampCorr_%s",mapIt.first.c_str()),"",1000,CTRRanges_mean[mapIt.first]-5.*CTRRanges_sigma[mapIt.first],CTRRanges_mean[mapIt.first]+5.*CTRRanges_sigma[mapIt.first]);
  }
  
  std::map<std::string,TF1*> fitFunc_corrPos;
  
  entry = 0;
  for(auto event: events)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 3rd loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;    
    ++entry;
    
    std::map<std::string,float> amps = event.amps;
    std::map<std::string,float> times = event.times;
    
    
    for(auto ch: channels)
    {
      int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
      if( index < 0 ) continue;
      
      std::vector<std::string> ampChs = opts.GetOpt<std::vector<std::string> >(Form("%s.ampCh", ch.c_str()));
      std::vector<std::string> timeChs = opts.GetOpt<std::vector<std::string> >(Form("%s.timeCh",ch.c_str()));      
      std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
      
      std::string ampCh1 = ampChs.at(0);
      std::string ampCh2 = ampChs.at(1);
      std::string timeCh1 = timeChs.at(0);
      std::string timeCh2 = timeChs.at(1);
      
      std::string refCh = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
      std::string timeChRef = opts.GetOpt<std::string>(Form("%s.timeCh",refCh.c_str()));
      std::string timeMethodRef = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",refCh.c_str())).at(0);
      
      for(auto timeMethod : timeMethods)
      {
        std::string label1 = ch+"_"+timeCh1+"_"+timeMethod;
        std::string label2 = ch+"_"+timeCh2+"_"+timeMethod;
        std::string labelRef = refCh+"_"+timeChRef+"_"+timeMethodRef;
        
        float timeLow1 = CTRRanges_mean[ch+"_"+timeMethod+"__1-ref"]-4.*CTRRanges_sigma[ch+"_"+timeMethod+"__1-ref"];
        float timeHig1 = CTRRanges_mean[ch+"_"+timeMethod+"__1-ref"]+4.*CTRRanges_sigma[ch+"_"+timeMethod+"__1-ref"];
        float timeLow2 = CTRRanges_mean[ch+"_"+timeMethod+"__2-ref"]-4.*CTRRanges_sigma[ch+"_"+timeMethod+"__2-ref"];
        float timeHig2 = CTRRanges_mean[ch+"_"+timeMethod+"__2-ref"]+4.*CTRRanges_sigma[ch+"_"+timeMethod+"__2-ref"];
        float timeLowAvg = CTRRanges_mean[ch+"_"+timeMethod+"__avg-ref"]-4.*CTRRanges_sigma[ch+"_"+timeMethod+"__avg-ref"];
        float timeHigAvg = CTRRanges_mean[ch+"_"+timeMethod+"__avg-ref"]+4.*CTRRanges_sigma[ch+"_"+timeMethod+"__avg-ref"];
        
        if( (times.count(label1) > 0) && (times.count(labelRef) > 0) )
        {
          float CTR = times[label1]-times[labelRef];

          if( CTR >= timeLow1 && CTR <= timeHig1 )
          {
            h_CTR[ch+"_"+timeMethod+"__1-ref"] -> Fill( CTR );
            
            float CTR_ampCorr = CTR - fitFunc_corrAmp[label1]->Eval(amps[ch+"_"+ampCh1]) + fitFunc_corrAmp[label1]->Eval( h_amp_cut[ch+"_"+ampCh1]->GetMean() );
            h_CTR_ampCorr[ch+"_"+timeMethod+"__1-ref"] -> Fill( CTR_ampCorr );
            
            p_time_ampCorr_vs_amp[label1] -> Fill( amps[ch+"_"+ampCh1],CTR_ampCorr );
            p_time_ampCorr_vs_X[label1] -> Fill( event.x,CTR_ampCorr );
            p_time_ampCorr_vs_Y[label1] -> Fill( event.y,CTR_ampCorr );
            p2_time_ampCorr_vs_XY[label1] -> Fill( event.x,event.y,CTR_ampCorr );
          }
        }
        if( (times.count(label2) > 0) && (times.count(labelRef) > 0) )
        {
          float CTR = times[label2]-times[labelRef];
          
          if( CTR >= timeLow2 && CTR <= timeHig2 )
          {
            h_CTR[ch+"_"+timeMethod+"__2-ref"] -> Fill( CTR );
            
            float CTR_ampCorr = CTR - fitFunc_corrAmp[label2]->Eval(amps[ch+"_"+ampCh2]) + fitFunc_corrAmp[label2]->Eval( h_amp_cut[ch+"_"+ampCh2]->GetMean() );
            h_CTR_ampCorr[ch+"_"+timeMethod+"__2-ref"] -> Fill( CTR_ampCorr );
            
            p_time_ampCorr_vs_amp[label2] -> Fill( amps[ch+"_"+ampCh2],CTR_ampCorr );
            p_time_ampCorr_vs_X[label2] -> Fill( event.x,CTR_ampCorr );
            p_time_ampCorr_vs_Y[label2] -> Fill( event.y,CTR_ampCorr );
            p2_time_ampCorr_vs_XY[label2] -> Fill( event.x,event.y,CTR_ampCorr );
          }
        }
        if( (times.count(label1) > 0) && (times.count(label2) > 0) && (times.count(labelRef) > 0) )
        {
          float CTR = 0.5*(times[label1]+times[label2])-times[labelRef];
            
          if( CTR >= timeLowAvg && CTR <= timeHigAvg )
          {
            h_CTR[ch+"_"+timeMethod+"__avg-ref"] -> Fill( CTR );
            
            float CTR_ampCorr = CTR - 0.5*fitFunc_corrAmp[label1]->Eval(amps[ch+"_"+ampCh1]) + 0.5*fitFunc_corrAmp[label1]->Eval( h_amp_cut[ch+"_"+ampCh1]->GetMean() )
                                    - 0.5*fitFunc_corrAmp[label2]->Eval(amps[ch+"_"+ampCh2]) + 0.5*fitFunc_corrAmp[label2]->Eval( h_amp_cut[ch+"_"+ampCh2]->GetMean() );
            h_CTR_ampCorr[ch+"_"+timeMethod+"__avg-ref"] -> Fill( CTR_ampCorr );
            
            p_time_ampCorr_vs_X[ch+"_"+timeMethod+"__avg-ref"] -> Fill( event.x,CTR_ampCorr );
            p_time_ampCorr_vs_Y[ch+"_"+timeMethod+"__avg-ref"] -> Fill( event.y,CTR_ampCorr );
            p2_time_ampCorr_vs_XY[ch+"_"+timeMethod+"__avg-ref"] -> Fill( event.x,event.y,CTR_ampCorr );
          }
        }
        
      }
    }
  }
  std::cout << "\n>>> end 3rd loop" << std::endl;
  
  
  
  //--- draw 3rd plots
  std::map<std::string,CTRResult> ctrResults_ampCorr;
  int rebin = opts.GetOpt<int>("Input.rebin");
  TLatex* latexLabel12 = new TLatex(0.13,0.97,Form("%s",""));
  
  std::map<std::string,TGraphErrors*> g_tRes_vs_thr;
  
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::vector<std::string> ampChs = opts.GetOpt<std::vector<std::string> >(Form("%s.ampCh",ch.c_str()));      
    std::vector<std::string> timeChs = opts.GetOpt<std::vector<std::string> >(Form("%s.timeCh",ch.c_str()));    
    std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
    
    float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str()));
    float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
    
    if( index < 0 ) continue;
    
    int timeIt = 0;
    for(auto timeCh : timeChs)
    {
      std::string ampCh = ampChs.at(timeIt);
      float ampMin = opts.GetOpt<std::vector<float> >(Form("%s.ampMin",ch.c_str())).at(timeIt);
      float ampMax = opts.GetOpt<std::vector<float> >(Form("%s.ampMax",ch.c_str())).at(timeIt);
      
      for(auto timeMethod : timeMethods)
      {
        std::string label = ch+"_"+timeCh+"_"+timeMethod;
        
        float timeLow = CTRRanges_mean[ch+"_"+timeMethod+Form("__%d-ref",timeIt+1)]-4.*CTRRanges_sigma[ch+"_"+timeMethod+Form("__%d-ref",timeIt+1)];
        float timeHig = CTRRanges_mean[ch+"_"+timeMethod+Form("__%d-ref",timeIt+1)]+3.*CTRRanges_sigma[ch+"_"+timeMethod+Form("__%d-ref",timeIt+1)];
        
        c = new TCanvas();
        DrawProfile(opts,ch,p_time_ampCorr_vs_amp[label],";amplitude [V];#Delta_{t} [ns]",ampMin,ampMax,timeLow,timeHig,kBlack,"",latexLabels[ch]);
        PrintCanvas(c,opts,ch,plotDir,Form("3_time_ampCorr_vs_amp_%d_%s",timeIt,timeMethod.c_str()));
        delete c;
      }
      
      ++timeIt;
    }
    
    
    for(auto timeMethod : timeMethods)
    {
      std::string label1 = ch+"_"+timeChs.at(0)+"_"+timeMethod;
      std::string label2 = ch+"_"+timeChs.at(1)+"_"+timeMethod;
      
      float timeLow1 = CTRRanges_mean[ch+"_"+timeMethod+Form("__%d-ref",1)]-4.*CTRRanges_sigma[ch+"_"+timeMethod+Form("__%d-ref",1)];
      float timeHig1 = CTRRanges_mean[ch+"_"+timeMethod+Form("__%d-ref",1)]+3.*CTRRanges_sigma[ch+"_"+timeMethod+Form("__%d-ref",1)];
      float timeLow2 = CTRRanges_mean[ch+"_"+timeMethod+Form("__%d-ref",2)]-4.*CTRRanges_sigma[ch+"_"+timeMethod+Form("__%d-ref",2)];
      float timeHig2 = CTRRanges_mean[ch+"_"+timeMethod+Form("__%d-ref",2)]+3.*CTRRanges_sigma[ch+"_"+timeMethod+Form("__%d-ref",2)];
      
      
      fitFunc_corrPos[ch+"_"+timeMethod+"__avg-ref"] = new TF1(Form("fitFunc_corrPos_%s",(ch+"_"+timeMethod+"__avg-ref").c_str()),"pol2",xtalXMin,xtalXMax);
      p_time_ampCorr_vs_X[ch+"_"+timeMethod+"__avg-ref"] -> Fit(fitFunc_corrPos[ch+"_"+timeMethod+"__avg-ref"],"QNRS+","",xtalXMin,xtalXMax);
      fitFunc_corrPos[ch+"_"+timeMethod+"__avg-ref"] -> SetLineColor(kBlack);
      
      std::vector<TProfile*> profs;
      profs.push_back(p_time_ampCorr_vs_X[label1]);
      profs.push_back(p_time_ampCorr_vs_X[label2]);
      profs.push_back(p_time_ampCorr_vs_X[ch+"_"+timeMethod+"__avg-ref"]);
      std::vector<int> colors;
      colors.push_back(kRed);
      colors.push_back(kBlue);
      colors.push_back(kBlack);
      
      c = new TCanvas();
      DrawProfile(opts,ch,profs,";x [mm];#Delta_{t} [ns]",xtalXMin,xtalXMax,std::min(timeLow1,timeLow2),std::max(timeHig1,timeHig2),colors,"",latexLabels[ch]);
      fitFunc_corrPos[ch+"_"+timeMethod+"__avg-ref"] -> Draw("same");
      PrintCanvas(c,opts,ch,plotDir,Form("3_time_ampCorr_vs_X_%s",timeMethod.c_str()));
      delete c;
    }
    
    
    float intrinsic = opts.GetOpt<float>(Form("%s.intrinsic",ch.c_str()));
    std::string shortLabel = opts.GetOpt<std::string>(Form("%s.shortLabel",ch.c_str()));    
    
    g_tRes_vs_thr[ch+"_ampCorr"] = new TGraphErrors();
    
    for(auto timeMethod : timeMethods)
    {
      c = new TCanvas(Form("c_CTR_ampCorr"),Form("c_CTR_ampCorr"));
      
      TH1F* histoAvg = h_CTR_ampCorr[ch+"_"+timeMethod+"__avg-ref"];
      TH1F* histo1 = h_CTR_ampCorr[ch+"_"+timeMethod+"__1-ref"];
      TH1F* histo2 = h_CTR_ampCorr[ch+"_"+timeMethod+"__2-ref"];
      
      if( intrinsic < 0 )
      {
        ctrResults_ampCorr[ch+"_"+timeMethod+"__1-ref"]   = drawCTRPlot(histo1,"amp. walk corr.",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
        ctrResults_ampCorr[ch+"_"+timeMethod+"__2-ref"]   = drawCTRPlot(histo2,"amp. walk corr.",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
        ctrResults_ampCorr[ch+"_"+timeMethod+"__avg-ref"] = drawCTRPlot(histoAvg,"amp. walk corr.",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,histo1,"1",histo2,"2");
      }      
      else
      {
        ctrResults_ampCorr[ch+"_"+timeMethod+"__1-ref"] = drawCTRPlot(histo1,"amp. walk corr.",rebin,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
        ctrResults_ampCorr[ch+"_"+timeMethod+"__2-ref"] = drawCTRPlot(histo2,"amp. walk corr.",rebin,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
        ctrResults_ampCorr[ch+"_"+timeMethod+"__avg-ref"] = drawCTRPlot(histoAvg,"amp. walk corr.",rebin,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,histo1,"1",histo2,"2");
      }
      
      // t_CTR_ampCorr_effSigma[ch] = ctrResults.effSigma;
      // t_CTR_ampCorr_gausSigma[ch] = ctrResults.gausSigma;
      // t_CTR_ampCorr_gausSigmaErr[ch] = ctrResults.gausSigmaErr;
      
      latexLabels[ch] -> Draw("same");
      
      gPad -> Update();
      
      c -> Print(Form("%s/c_3_tRes_ampCorr_%s_%s.png",plotDir.c_str(),timeMethod.c_str(),shortLabel.c_str()));
      c -> Print(Form("%s/c_3_tRes_ampCorr_%s_%s.pdf",plotDir.c_str(),timeMethod.c_str(),shortLabel.c_str()));
      
      
      std::string tm(timeMethod);
      double x = atof((tm.erase(0,3)).c_str())/4096.;      
      g_tRes_vs_thr[ch+"_ampCorr"] -> SetPoint(g_tRes_vs_thr[ch+"_ampCorr"]->GetN(),x,ctrResults_ampCorr[ch+"_"+timeMethod+"__avg-ref"].gausSigma);
      g_tRes_vs_thr[ch+"_ampCorr"] -> SetPointError(g_tRes_vs_thr[ch+"_ampCorr"]->GetN()-1,0.,ctrResults_ampCorr[ch+"_"+timeMethod+"__avg-ref"].gausSigmaErr);
      
      delete c;
    }
  }
  
  
  
  
  
  
  //------------------------
  //--- 4th loop over events
  std::map<std::string,TH1F*> h_CTR_ampCorr_posCorr;
  std::map<std::string,TH1F*> h_CTR_ampCorr_posCorr_weighted;
  for(auto mapIt : h_CTR_raw)
  {
    h_CTR_ampCorr_posCorr[mapIt.first] = new TH1F(Form("h_CTR_ampCorr_posCorr_%s",mapIt.first.c_str()),"",1000,CTRRanges_mean[mapIt.first]-5.*CTRRanges_sigma[mapIt.first],CTRRanges_mean[mapIt.first]+5.*CTRRanges_sigma[mapIt.first]);
    h_CTR_ampCorr_posCorr_weighted[mapIt.first] = new TH1F(Form("h_CTR_ampCorr_posCorr_weighted_%s",mapIt.first.c_str()),"",1000,CTRRanges_mean[mapIt.first]-5.*CTRRanges_sigma[mapIt.first],CTRRanges_mean[mapIt.first]+5.*CTRRanges_sigma[mapIt.first]);
  }
  
  entry = 0;
  for(auto event: events)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 4th loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;    
    ++entry;
    
    std::map<std::string,float> amps = event.amps;
    std::map<std::string,float> times = event.times;
    
    
    for(auto ch: channels)
    {
      int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
      if( index < 0 ) continue;
      
      std::vector<std::string> ampChs = opts.GetOpt<std::vector<std::string> >(Form("%s.ampCh", ch.c_str()));
      std::vector<std::string> timeChs = opts.GetOpt<std::vector<std::string> >(Form("%s.timeCh",ch.c_str()));      
      std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
      
      float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str()));
      float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
      
      std::string ampCh1 = ampChs.at(0);
      std::string ampCh2 = ampChs.at(1);
      std::string timeCh1 = timeChs.at(0);
      std::string timeCh2 = timeChs.at(1);
      
      std::string refCh = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
      std::string timeChRef = opts.GetOpt<std::string>(Form("%s.timeCh",refCh.c_str()));
      std::string timeMethodRef = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",refCh.c_str())).at(0);
      
      for(auto timeMethod : timeMethods)
      {
        std::string label1 = ch+"_"+timeCh1+"_"+timeMethod;
        std::string label2 = ch+"_"+timeCh2+"_"+timeMethod;
        std::string labelRef = refCh+"_"+timeChRef+"_"+timeMethodRef;
        
        float timeLow1 = CTRRanges_mean[ch+"_"+timeMethod+"__1-ref"]-4.*CTRRanges_sigma[ch+"_"+timeMethod+"__1-ref"];
        float timeHig1 = CTRRanges_mean[ch+"_"+timeMethod+"__1-ref"]+4.*CTRRanges_sigma[ch+"_"+timeMethod+"__1-ref"];
        float timeLow2 = CTRRanges_mean[ch+"_"+timeMethod+"__2-ref"]-4.*CTRRanges_sigma[ch+"_"+timeMethod+"__2-ref"];
        float timeHig2 = CTRRanges_mean[ch+"_"+timeMethod+"__2-ref"]+4.*CTRRanges_sigma[ch+"_"+timeMethod+"__2-ref"];
        float timeLowAvg = CTRRanges_mean[ch+"_"+timeMethod+"__avg-ref"]-4.*CTRRanges_sigma[ch+"_"+timeMethod+"__avg-ref"];
        float timeHigAvg = CTRRanges_mean[ch+"_"+timeMethod+"__avg-ref"]+4.*CTRRanges_sigma[ch+"_"+timeMethod+"__avg-ref"];
        
        if( (times.count(label1) > 0) && (times.count(label2) > 0) && (times.count(labelRef) > 0) )
        {
          float CTR = 0.5*(times[label1]+times[label2])-times[labelRef];
          
          if( CTR >= timeLowAvg && CTR <= timeHigAvg )
          {
            h_CTR[ch+"_"+timeMethod+"__avg-ref"] -> Fill( CTR );
            
            float CTR_ampCorr = CTR - 0.5*fitFunc_corrAmp[label1]->Eval(amps[ch+"_"+ampCh1]) + 0.5*fitFunc_corrAmp[label1]->Eval( h_amp_cut[ch+"_"+ampCh1]->GetMean() )
                                    - 0.5*fitFunc_corrAmp[label2]->Eval(amps[ch+"_"+ampCh2]) + 0.5*fitFunc_corrAmp[label2]->Eval( h_amp_cut[ch+"_"+ampCh2]->GetMean() );
            
            float CTR_ampCorr_posCorr = CTR_ampCorr - fitFunc_corrPos[ch+"_"+timeMethod+"__avg-ref"]->Eval(event.x) + fitFunc_corrPos[ch+"_"+timeMethod+"__avg-ref"]->Eval(0.5*(xtalXMin+xtalXMax));
            
            h_CTR_ampCorr_posCorr[ch+"_"+timeMethod+"__avg-ref"] -> Fill( CTR_ampCorr_posCorr );
            
            float sigma1 = ctrResults_ampCorr[ch+"_"+timeMethod+"__1-ref"].gausSigma;
            float sigma2 = ctrResults_ampCorr[ch+"_"+timeMethod+"__2-ref"].gausSigma;
            float time1 = times[label1] - fitFunc_corrAmp[label1]->Eval(amps[ch+"_"+ampCh1]) + fitFunc_corrAmp[label1]->Eval( h_amp_cut[ch+"_"+ampCh1]->GetMean() );
            float time2 = times[label2] - fitFunc_corrAmp[label2]->Eval(amps[ch+"_"+ampCh2]) + fitFunc_corrAmp[label2]->Eval( h_amp_cut[ch+"_"+ampCh2]->GetMean() );
            float CTR_ampCorr_posCorr_weighted = ( time1/pow(sigma1,2) + time2/pow(sigma2,2) ) / ( 1/pow(sigma1,2) + 1/pow(sigma2,2) ) - times[labelRef] - fitFunc_corrPos[ch+"_"+timeMethod+"__avg-ref"]->Eval(event.x) + fitFunc_corrPos[ch+"_"+timeMethod+"__avg-ref"]->Eval(0.5*(xtalXMin+xtalXMax));
            
            h_CTR_ampCorr_posCorr_weighted[ch+"_"+timeMethod+"__avg-ref"] -> Fill( CTR_ampCorr_posCorr_weighted  );
          }
        }
        
      }
    }
  }
  std::cout << "\n>>> end 4th loop" << std::endl;
  
  
  
  //--- draw 4th plots
  std::map<std::string,CTRResult> ctrResults_ampCorr_posCorr;
  std::map<std::string,CTRResult> ctrResults_ampCorr_posCorr_weighted;

  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::vector<std::string> ampChs = opts.GetOpt<std::vector<std::string> >(Form("%s.ampCh",ch.c_str()));      
    std::vector<std::string> timeChs = opts.GetOpt<std::vector<std::string> >(Form("%s.timeCh",ch.c_str()));    
    std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
    
    if( index < 0 ) continue;
    
    float intrinsic = opts.GetOpt<float>(Form("%s.intrinsic",ch.c_str()));
    std::string shortLabel = opts.GetOpt<std::string>(Form("%s.shortLabel",ch.c_str()));    
    
    g_tRes_vs_thr[ch+"_ampCorr_posCorr"] = new TGraphErrors();
    g_tRes_vs_thr[ch+"_ampCorr_posCorr_weighted"] = new TGraphErrors();
        
    for(auto timeMethod : timeMethods)
    {
      c = new TCanvas(Form("c_CTR_ampCorr_posCorr"),Form("c_CTR_ampCorr_posCorr"));
      
      TH1F* histo_ampCorr_posCorr = h_CTR_ampCorr_posCorr[ch+"_"+timeMethod+"__avg-ref"];
      TH1F* histo_ampCorr = h_CTR_ampCorr[ch+"_"+timeMethod+"__avg-ref"];
      
      if( intrinsic < 0 )
        ctrResults_ampCorr_posCorr[ch+"_"+timeMethod+"__avg-ref"] = drawCTRPlot(histo_ampCorr_posCorr,"amp. walk + pos. corr.",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,histo_ampCorr,"amp. walk corr.",NULL,"");
      else
        ctrResults_ampCorr_posCorr[ch+"_"+timeMethod+"__avg-ref"] = drawCTRPlot(histo_ampCorr_posCorr,"amp. walk + pos. corr.",rebin,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,histo_ampCorr,"amp. walk corr.",NULL,"");
      // t_CTR_ampCorr_effSigma[ch] = ctrResults.effSigma;
      // t_CTR_ampCorr_gausSigma[ch] = ctrResults.gausSigma;
      // t_CTR_ampCorr_gausSigmaErr[ch] = ctrResults.gausSigmaErr;
      
      latexLabels[ch] -> Draw("same");
      
      gPad -> Update();
      
      c -> Print(Form("%s/c_4_tRes_ampCorr_posCorr_%s_%s.png",plotDir.c_str(),timeMethod.c_str(),shortLabel.c_str()));
      c -> Print(Form("%s/c_4_tRes_ampCorr_posCorr_%s_%s.pdf",plotDir.c_str(),timeMethod.c_str(),shortLabel.c_str()));
      
      delete c;
      
      
      c = new TCanvas(Form("c_CTR_ampCorr_posCorr_weighted"),Form("c_CTR_ampCorr_posCorr_weighted"));
      
      TH1F* histo_ampCorr_posCorr_weighted = h_CTR_ampCorr_posCorr_weighted[ch+"_"+timeMethod+"__avg-ref"];
      
      if( intrinsic < 0 )
        ctrResults_ampCorr_posCorr_weighted[ch+"_"+timeMethod+"__avg-ref"] = drawCTRPlot(histo_ampCorr_posCorr_weighted,"amp. walk + pos. corr. weighted",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
      else
        ctrResults_ampCorr_posCorr_weighted[ch+"_"+timeMethod+"__avg-ref"] = drawCTRPlot(histo_ampCorr_posCorr_weighted,"amp. walk + pos. corr. weighted",rebin,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
      // t_CTR_ampCorr_effSigma[ch] = ctrResults.effSigma;
      // t_CTR_ampCorr_gausSigma[ch] = ctrResults.gausSigma;
      // t_CTR_ampCorr_gausSigmaErr[ch] = ctrResults.gausSigmaErr;
      
      latexLabels[ch] -> Draw("same");
      
      gPad -> Update();
      
      c -> Print(Form("%s/c_4_tRes_ampCorr_posCorr_weighted_%s_%s.png",plotDir.c_str(),timeMethod.c_str(),shortLabel.c_str()));
      c -> Print(Form("%s/c_4_tRes_ampCorr_posCorr_weighted_%s_%s.pdf",plotDir.c_str(),timeMethod.c_str(),shortLabel.c_str()));
      
      delete c;
      
      
      std::string tm(timeMethod);
      double x = atof((tm.erase(0,3)).c_str())/4096.;
      g_tRes_vs_thr[ch+"_ampCorr_posCorr"] -> SetPoint(g_tRes_vs_thr[ch+"_ampCorr_posCorr"]->GetN(),x,ctrResults_ampCorr_posCorr[ch+"_"+timeMethod+"__avg-ref"].gausSigma);
      g_tRes_vs_thr[ch+"_ampCorr_posCorr"] -> SetPointError(g_tRes_vs_thr[ch+"_ampCorr_posCorr"]->GetN()-1,0.,ctrResults_ampCorr_posCorr[ch+"_"+timeMethod+"__avg-ref"].gausSigmaErr);
      
      g_tRes_vs_thr[ch+"_ampCorr_posCorr_weighted"] -> SetPoint(g_tRes_vs_thr[ch+"_ampCorr_posCorr_weighted"]->GetN(),x,ctrResults_ampCorr_posCorr_weighted[ch+"_"+timeMethod+"__avg-ref"].gausSigma);
      g_tRes_vs_thr[ch+"_ampCorr_posCorr_weighted"] -> SetPointError(g_tRes_vs_thr[ch+"_ampCorr_posCorr_weighted"]->GetN()-1,0.,ctrResults_ampCorr_posCorr_weighted[ch+"_"+timeMethod+"__avg-ref"].gausSigmaErr);
    }
    
    
    c = new TCanvas();
    std::vector<TGraph*> graphs;
    graphs.push_back(g_tRes_vs_thr[ch+"_ampCorr"]);
    graphs.push_back(g_tRes_vs_thr[ch+"_ampCorr_posCorr"]);
    graphs.push_back(g_tRes_vs_thr[ch+"_ampCorr_posCorr_weighted"]);
    std::vector<int> colors;
    colors.push_back(kBlack);
    colors.push_back(kOrange);
    colors.push_back(kViolet);
    DrawGraph(opts,ch,graphs,";threshold [mV];#sigma_{t} [ps]",0.001,1.,20.,100.,colors,"PL",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,Form("4_timeRes_vs_thr"));
    delete c;
  }
  
  
  
  /*
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
    float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));
    float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str())); 
    float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
    float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
    float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
    float timeLow = CTRRanges_mean[ch]-2.*CTRRanges_sigma[ch];
    float timeHig = CTRRanges_mean[ch]+1.*CTRRanges_sigma[ch];
    
    TLine* line_xtalXLow = new TLine(0.5*(xtalXMin+xtalXMax)-1.5,xtalYMin,0.5*(xtalXMin+xtalXMax)-1.5,xtalYMax);
    line_xtalXLow -> SetLineStyle(2);
    line_xtalXLow -> SetLineWidth(3);
    line_xtalXLow -> SetLineColor(kBlue);
    TLine* line_xtalXHig = new TLine(0.5*(xtalXMin+xtalXMax)+1.5,xtalYMin,0.5*(xtalXMin+xtalXMax)+1.5,xtalYMax);
    line_xtalXHig -> SetLineStyle(2);
    line_xtalXHig -> SetLineWidth(3);
    line_xtalXHig -> SetLineColor(kBlue);
    
    TLine* line_xtalYLow = new TLine(xtalXMin,0.5*(xtalYMin+xtalYMax)-1.5,xtalXMax,0.5*(xtalYMin+xtalYMax)-1.5);
    line_xtalYLow -> SetLineStyle(2);
    line_xtalYLow -> SetLineWidth(3);
    line_xtalYLow -> SetLineColor(kRed);
    TLine* line_xtalYHig = new TLine(xtalXMin,0.5*(xtalYMin+xtalYMax)+1.5,xtalXMax,0.5*(xtalYMin+xtalYMax)+1.5);
    line_xtalYHig -> SetLineStyle(2);
    line_xtalYHig -> SetLineWidth(3);
    line_xtalYHig -> SetLineColor(kRed);
    
    std::vector<TLine*> lines_xtal;
    lines_xtal.push_back(line_xtalXLow);
    lines_xtal.push_back(line_xtalXHig);
    lines_xtal.push_back(line_xtalYLow);
    lines_xtal.push_back(line_xtalYHig);

    
    c = new TCanvas();
    DrawProfile(opts,ch,p_time_ampCorr_vs_amp[ch],";amplitude [V];amp.-corrected #Delta_{t} [ns]",ampMin,ampMax,timeLow,timeHig,kBlack,"",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,"3_time_ampCorr_vs_amp");
    delete c;
    
    c = new TCanvas("c","c",2000,600);
    c -> Divide(3,1);
    c->cd(1); DrawProfile  (opts,ch,p_time_ampCorr_vs_X[ch],    ";hodoscope x [mm];amp.-corrected #Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(1); DrawProfile  (opts,ch,p_time_cut_ampCorr_vs_X[ch],";hodoscope x [mm];amp.-corrected #Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kRed,  "same",latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_ampCorr_vs_Y[ch],    ";hodoscope y [mm];amp.-corrected #Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_cut_ampCorr_vs_Y[ch],";hodoscope y [mm];amp.-corrected #Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlue, "same",latexLabels[ch]);
    c->cd(3); DrawProfile2D(opts,ch,p2_time_ampCorr_vs_XY[ch],  ";hodoscope x [mm];hodoscope y [mm];amp.-corrected #Delta_{t} [ns]",xtalXMin,xtalXMax,xtalYMin,xtalYMax,timeLow,timeHig,latexLabels[ch],&lines_xtal);
    PrintCanvas(c,opts,ch,plotDir,"3_time_ampCorr_vs_XY");
    delete c;
  }
  
  */
  /*
  //------------------------
  //--- 4th loop over events
  int rebin = opts.GetOpt<int>("Input.rebin");
  TLatex* latexLabel12 = new TLatex(0.13,0.97,Form("%s",""));

  
  std::map<std::string,TH1F*> h_CTR_ampCorr_posCorr;
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    h_CTR_ampCorr_posCorr[ch] = new TH1F(Form("h_CTR_ampCorr_posCorr_%s",ch.c_str()),"",1000,CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+5.*CTRRanges_sigma[ch]);
    h_CTR_ampCorr_posCorr[ch] -> Sumw2();
  }
  
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 4th loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    h4 -> GetEntry(entry);
    // if( !goodSpills[std::make_pair(tv.run,tv.spill)] ) continue;
    
    
    for(auto ch: channels)
    {
      int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
      float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str()));
      float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
      float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
      float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
      if( index < 0 ) continue;
      
      // reconstruct position
      ReconstructHodoPosition(tv,opts,ch,bool(doTracking),hodoPlane);      
      
      // default event selection
      AnalysisVars av;
      if( !AcceptEvent(av,tv,opts,ch,nFibresMin,nFibresMax,(CTRRanges_mean[ch]-3.*CTRRanges_sigma[ch]),(CTRRanges_mean[ch]+3.*CTRRanges_sigma[ch])) ) continue;
      if( energySmearing > 0. ) av.amp = av.amp * r.Gaus(1.,energySmearing);
      
      float CTR = av.time - av.timeRef;
      float CTR_ampCorr = CTR - fitFunc_corrAmp[ch]->Eval(av.amp) + fitFunc_corrAmp[ch]->Eval( h_amp_cut[ch]->GetMean() );
      // float CTR_ampCorr_posCorr = CTR_ampCorr - p2_time_ampCorr_vs_XY[ch]->GetBinContent(p2_time_ampCorr_vs_XY[ch]->FindBin(tv.beamX,tv.beamY))
      //                                         + p2_time_ampCorr_vs_XY[ch]->GetBinContent(p2_time_ampCorr_vs_XY[ch]->FindBin(0.5*(xtalXMin+xtalXMax)+2,0.5*(xtalYMin+xtalYMax)+2));
      float CTR_ampCorr_posCorr = CTR_ampCorr - p_time_ampCorr_vs_X[ch]->GetBinContent(p_time_ampCorr_vs_X[ch]->FindBin(tv.beamX))
                                              + p_time_ampCorr_vs_X[ch]->GetBinContent(p_time_ampCorr_vs_X[ch]->FindBin(0.5*(xtalXMin+xtalXMax)+2))
                                              - p_time_ampCorr_vs_Y[ch]->GetBinContent(p_time_ampCorr_vs_Y[ch]->FindBin(tv.beamY))
                                              + p_time_ampCorr_vs_Y[ch]->GetBinContent(p_time_ampCorr_vs_Y[ch]->FindBin(0.5*(xtalYMin+xtalYMax)+2));
      
      p_time_ampCorr_posCorr_vs_X[ch] -> Fill( tv.beamX,CTR_ampCorr_posCorr );
      p_time_ampCorr_posCorr_vs_Y[ch] -> Fill( tv.beamY,CTR_ampCorr_posCorr );
      p2_time_ampCorr_posCorr_vs_XY[ch] -> Fill( tv.beamX,tv.beamY,CTR_ampCorr_posCorr );
      
      int binX = h_occupancyX[ch] -> Fill(tv.beamX);
      int binY = h_occupancyY[ch] -> Fill(tv.beamY);      
      if( binX >= 0 && binY >= 0 )
      {
        float weight = h2_occupancyXY[ch] -> GetBinContent(binX,binY) / nSelectedEntries[ch];
        
        h_CTR_ampCorr_posCorr[ch] -> Fill( CTR_ampCorr_posCorr,weight );
      }
    }
  }
  std::cout << "\n>>> end 4th loop" << std::endl;
  
  
  
  //--- draw 4th plots
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;

    float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
    float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));
    float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str())); 
    float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
    float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
    float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
    float timeLow = CTRRanges_mean[ch]-2.5*CTRRanges_sigma[ch];
    float timeHig = CTRRanges_mean[ch]+1.5*CTRRanges_sigma[ch];
    
    TLine* line_xtalXLow = new TLine(0.5*(xtalXMin+xtalXMax)-1.5,xtalYMin,0.5*(xtalXMin+xtalXMax)-1.5,xtalYMax);
    line_xtalXLow -> SetLineStyle(2);
    line_xtalXLow -> SetLineWidth(3);
    line_xtalXLow -> SetLineColor(kBlue);
    TLine* line_xtalXHig = new TLine(0.5*(xtalXMin+xtalXMax)+1.5,xtalYMin,0.5*(xtalXMin+xtalXMax)+1.5,xtalYMax);
    line_xtalXHig -> SetLineStyle(2);
    line_xtalXHig -> SetLineWidth(3);
    line_xtalXHig -> SetLineColor(kBlue);
    
    TLine* line_xtalYLow = new TLine(xtalXMin,0.5*(xtalYMin+xtalYMax)-1.5,xtalXMax,0.5*(xtalYMin+xtalYMax)-1.5);
    line_xtalYLow -> SetLineStyle(2);
    line_xtalYLow -> SetLineWidth(3);
    line_xtalYLow -> SetLineColor(kRed);
    TLine* line_xtalYHig = new TLine(xtalXMin,0.5*(xtalYMin+xtalYMax)+1.5,xtalXMax,0.5*(xtalYMin+xtalYMax)+1.5);
    line_xtalYHig -> SetLineStyle(2);
    line_xtalYHig -> SetLineWidth(3);
    line_xtalYHig -> SetLineColor(kRed);
    
    std::vector<TLine*> lines_xtal;
    lines_xtal.push_back(line_xtalXLow);
    lines_xtal.push_back(line_xtalXHig);
    lines_xtal.push_back(line_xtalYLow);
    lines_xtal.push_back(line_xtalYHig);
    
    
    c = new TCanvas("c","c",2000,600);
    c -> Divide(3,1);
    c->cd(1); DrawProfile  (opts,ch,p_time_ampCorr_posCorr_vs_X[ch],    ";hodoscope x [mm];amp.+pos.-corrected #Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(1); DrawProfile  (opts,ch,p_time_cut_ampCorr_posCorr_vs_X[ch],";hodoscope x [mm];amp.+pos.-corrected #Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kRed,  "same",latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_ampCorr_posCorr_vs_Y[ch],    ";hodoscope y [mm];amp.+pos.-corrected #Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_cut_ampCorr_posCorr_vs_Y[ch],";hodoscope y [mm];amp.+pos.-corrected #Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlue, "same",latexLabels[ch]);
    c->cd(3); DrawProfile2D(opts,ch,p2_time_ampCorr_posCorr_vs_XY[ch],  ";hodoscope x [mm];hodoscope y [mm];amp.+pos.-corrected #Delta_{t} [ns]",xtalXMin,xtalXMax,xtalYMin,xtalYMax,timeLow,timeHig,latexLabels[ch],&lines_xtal);
    PrintCanvas(c,opts,ch,plotDir,"4_time_ampCorr_posCorr_vs_XY");
    delete c;
  }
  
  
  
  //--- draw CTR plots
  CTRResult ctrResults;
  
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string shortLabel = opts.GetOpt<std::string>(Form("%s.shortLabel",ch.c_str()));
    float intrinsic = opts.GetOpt<float>(Form("%s.intrinsic",ch.c_str()));
    std::string refCh  = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
    std::string shortLabelRef = opts.GetOpt<std::string>(Form("%s.shortLabel",refCh.c_str()));
    if( index < 0 ) continue;
    
    c = new TCanvas(Form("c_CTR_%s",shortLabel.c_str()),Form("c_CTR_%s",shortLabel.c_str()));
    
    histo = h_CTR[ch];
    
    if( intrinsic < 0 )
      ctrResults = drawCTRPlot(histo,"raw",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
    else
      ctrResults = drawCTRPlot(histo,"raw",rebin,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
    t_CTR_effSigma[ch] = ctrResults.effSigma;
    t_CTR_gausSigma[ch] = ctrResults.gausSigma;
    t_CTR_gausSigmaErr[ch] = ctrResults.gausSigmaErr;
    
    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();

    c -> Print(Form("%s/c_4_tRes_%s.png",plotDir.c_str(),shortLabel.c_str()));
    c -> Print(Form("%s/c_4_tRes_%s.pdf",plotDir.c_str(),shortLabel.c_str()));
  }
  
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string shortLabel = opts.GetOpt<std::string>(Form("%s.shortLabel",ch.c_str()));
    float intrinsic = opts.GetOpt<float>(Form("%s.intrinsic",ch.c_str()));
    std::string refCh  = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
    std::string shortLabelRef = opts.GetOpt<std::string>(Form("%s.shortLabel",refCh.c_str()));
    if( index < 0 ) continue;
    
    c = new TCanvas(Form("c_CTR_ampCorr_%s",shortLabel.c_str()),Form("c_CTR_%s",shortLabel.c_str()));

    histo = h_CTR[ch];
    histoCorr = h_CTR_ampCorr[ch];
    
    if( intrinsic < 0 )
      ctrResults = drawCTRPlot(histoCorr,"amp. walk corr.",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,histo,"raw",NULL,"");
    else
      ctrResults = drawCTRPlot(histoCorr,"amp. walk corr.",rebin,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,histo,"raw",NULL,"");
    t_CTR_ampCorr_effSigma[ch] = ctrResults.effSigma;
    t_CTR_ampCorr_gausSigma[ch] = ctrResults.gausSigma;
    t_CTR_ampCorr_gausSigmaErr[ch] = ctrResults.gausSigmaErr;
    
    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();

    c -> Print(Form("%s/c_4_tRes_ampCorr_%s.png",plotDir.c_str(),shortLabel.c_str()));
    c -> Print(Form("%s/c_4_tRes_ampCorr_%s.pdf",plotDir.c_str(),shortLabel.c_str()));
  }
  
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string shortLabel = opts.GetOpt<std::string>(Form("%s.shortLabel",ch.c_str()));
    float intrinsic = opts.GetOpt<float>(Form("%s.intrinsic",ch.c_str()));
    std::string refCh  = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
    std::string shortLabelRef = opts.GetOpt<std::string>(Form("%s.shortLabel",refCh.c_str()));
    if( index < 0 ) continue;
    
    c = new TCanvas(Form("c_CTR_ampCorr_posCorr_%s",shortLabel.c_str()),Form("c_CTR_%s",shortLabel.c_str()));
    
    histo = h_CTR[ch];
    histoCorr = h_CTR_ampCorr[ch];
    histoCorr2 = h_CTR_ampCorr_posCorr[ch];
    
    if( intrinsic < 0 )
      ctrResults = drawCTRPlot(histoCorr2,"amp. walk + pos. corr.",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,histoCorr,"amp. corr",histo,"raw");
    else
      ctrResults = drawCTRPlot(histoCorr2,"amp. walk + pos. corr.",rebin,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,histoCorr,"amp. corr",histo,"raw");
    t_CTR_ampCorr_posCorr_effSigma[ch] = ctrResults.effSigma;
    t_CTR_ampCorr_posCorr_gausSigma[ch] = ctrResults.gausSigma;
    t_CTR_ampCorr_posCorr_gausSigmaErr[ch] = ctrResults.gausSigmaErr;
    
    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();

    c -> Print(Form("%s/c_4_tRes_ampCorr_posCorr_%s.png",plotDir.c_str(),shortLabel.c_str()));
    c -> Print(Form("%s/c_4_tRes_ampCorr_posCorr_%s.pdf",plotDir.c_str(),shortLabel.c_str()));
  }
  
  for(auto ch: channels)
  {
    int nBinsX = opts.GetOpt<int>(Form("%s.nBinsX",ch.c_str()));
    int nBinsY = opts.GetOpt<int>(Form("%s.nBinsY",ch.c_str()));
    float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str())); 
    float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
    float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
    float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
    float intrinsic = opts.GetOpt<float>(Form("%s.intrinsic",ch.c_str()));
    
    if( nBinsX <= 0 || nBinsY <= 0 ) continue;
    
    TH2F* h2_CTR_ampCorr_chess = new TH2F(Form("h2_CTR_ampCorr_chess_%s",ch.c_str()),"",nBinsX,xtalXMin,xtalXMax,nBinsY,xtalYMin,xtalYMax);
    
    for(int binX = 0; binX < nBinsX; ++binX)
      for(int binY = 0; binY < nBinsY; ++binY)
      {
        std::string label(Form("%s_%d-%d",ch.c_str(),binX,binY));
        histo =  map_h_CTR_ampCorr[label];
        if( histo->GetEntries() < 100 ) continue;

        CTRResult result;        
        if( intrinsic < 0 )
          result = drawCTRPlot(histo,"amp. walk corr.",1,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
        else
          result = drawCTRPlot(histo,"amp. walk corr.",1,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
        
        h2_CTR_ampCorr_chess -> SetBinContent( binX+1,binY+1,result.effSigma );
      }
    
    c = new TCanvas();
    DrawHistogram2D(opts,ch,h2_CTR_ampCorr_chess,";hodoscope x[mm];hodoscope y [mm];time resolution [ps]",xtalXMin,xtalXMax,xtalYMin,xtalYMax,30.,70.,false,false,"colz,text",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,"4_tRes_chess_ampCorr");
    delete c;    
  }
  
  outTree -> Fill();
  */
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
