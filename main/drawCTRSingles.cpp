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
    std::cout << ">>> drawCTRSingles::usage:   " << argv[0] << " configFile.cfg" << std::endl;
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
    // std::string fileName = Form("%s/%d/*.root",inputDir.c_str(),run);
    std::string fileName = Form("%s/%d.root",inputDir.c_str(),run);
    std::cout << ">>> Adding flle " << fileName << std::endl;
    h4 -> Add(fileName.c_str());
  }

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
  

  
  //--- define tree branch addresses
  TreeVars tv;
  InitTreeVars(h4,tv,opts);
  
  int nEntries = h4->GetEntries();
  std::cout << ">>> Events read: " << nEntries << std::endl;
  
  
  
  //------------------
  // define histograms

  std::string conf = opts.GetOpt<std::string>("Input.conf");
  TFile* outFile = TFile::Open(Form("%s/drawCTRSingles_%s.root",plotDir.c_str(),conf.c_str()),"RECREATE");
  outFile -> cd();

  TTree* outTree = new TTree("results","results");

  std::map<std::string,float> t_Vbias;
  std::map<std::string,float> t_NINOthr;
  std::map<std::string,float> t_CTR_effSigma;
  std::map<std::string,float> t_CTR_gausSigma;
  std::map<std::string,float> t_CTR_gausSigmaErr;
  std::map<std::string,float> t_CTR_ampCorr_effSigma;
  std::map<std::string,float> t_CTR_ampCorr_gausSigma;
  std::map<std::string,float> t_CTR_ampCorr_gausSigmaErr;
  std::map<std::string,float> t_CTR_ampCorr_posCorr_effSigma;
  std::map<std::string,float> t_CTR_ampCorr_posCorr_gausSigma;
  std::map<std::string,float> t_CTR_ampCorr_posCorr_gausSigmaErr;
  std::map<std::string,float> t_CTR_RTCorr_effSigma;
  std::map<std::string,float> t_CTR_RTCorr_gausSigma;
  std::map<std::string,float> t_CTR_RTCorr_gausSigmaErr;
  std::map<std::string,float> t_CTR_RTCorr_posCorr_effSigma;
  std::map<std::string,float> t_CTR_RTCorr_posCorr_gausSigma;
  std::map<std::string,float> t_CTR_RTCorr_posCorr_gausSigmaErr;
  std::map<std::string,float> t_CTR_distAmpCorr_effSigma;
  std::map<std::string,float> t_CTR_distAmpCorr_gausSigma;
  std::map<std::string,float> t_CTR_distAmpCorr_gausSigmaErr;
  std::map<std::string,float> t_CTR_distRTCorr_effSigma;
  std::map<std::string,float> t_CTR_distRTCorr_gausSigma;
  std::map<std::string,float> t_CTR_distRTCorr_gausSigmaErr;
  
  for(auto ch : channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    outTree -> Branch(Form("Vbias_%s",ch.c_str()),&t_Vbias[ch]);
    outTree -> Branch(Form("NINOthr_%s",ch.c_str()),&t_NINOthr[ch]);
    outTree -> Branch(Form("CTR_effSigma_%s",ch.c_str()),&t_CTR_effSigma[ch]);
    outTree -> Branch(Form("CTR_gausSigma_%s",ch.c_str()),&t_CTR_gausSigma[ch]);
    outTree -> Branch(Form("CTR_gausSigmaErr_%s",ch.c_str()),&t_CTR_gausSigmaErr[ch]);
    outTree -> Branch(Form("CTR_ampCorr_effSigma_%s",ch.c_str()),&t_CTR_ampCorr_effSigma[ch]);
    outTree -> Branch(Form("CTR_ampCorr_gausSigma_%s",ch.c_str()),&t_CTR_ampCorr_gausSigma[ch]);
    outTree -> Branch(Form("CTR_ampCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_ampCorr_gausSigmaErr[ch]);
    outTree -> Branch(Form("CTR_ampCorr_posCorr_effSigma_%s",ch.c_str()),&t_CTR_ampCorr_posCorr_effSigma[ch]);
    outTree -> Branch(Form("CTR_ampCorr_posCorr_gausSigma_%s",ch.c_str()),&t_CTR_ampCorr_posCorr_gausSigma[ch]);
    outTree -> Branch(Form("CTR_ampCorr_posCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_ampCorr_posCorr_gausSigmaErr[ch]);
    outTree -> Branch(Form("CTR_RTCorr_effSigma_%s",ch.c_str()),&t_CTR_RTCorr_effSigma[ch]);
    outTree -> Branch(Form("CTR_RTCorr_gausSigma_%s",ch.c_str()),&t_CTR_RTCorr_gausSigma[ch]);
    outTree -> Branch(Form("CTR_RTCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_RTCorr_gausSigmaErr[ch]);
    outTree -> Branch(Form("CTR_RTCorr_posCorr_effSigma_%s",ch.c_str()),&t_CTR_RTCorr_posCorr_effSigma[ch]);
    outTree -> Branch(Form("CTR_RTCorr_posCorr_gausSigma_%s",ch.c_str()),&t_CTR_RTCorr_posCorr_gausSigma[ch]);
    outTree -> Branch(Form("CTR_RTCorr_posCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_RTCorr_posCorr_gausSigmaErr[ch]);
    outTree -> Branch(Form("CTR_distAmpCorr_effSigma_%s",ch.c_str()),&t_CTR_distAmpCorr_effSigma[ch]);
    outTree -> Branch(Form("CTR_distAmpCorr_gausSigma_%s",ch.c_str()),&t_CTR_distAmpCorr_gausSigma[ch]);
    outTree -> Branch(Form("CTR_distAmpCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_distAmpCorr_gausSigmaErr[ch]);
    outTree -> Branch(Form("CTR_distRTCorr_effSigma_%s",ch.c_str()),&t_CTR_distRTCorr_effSigma[ch]);
    outTree -> Branch(Form("CTR_distRTCorr_gausSigma_%s",ch.c_str()),&t_CTR_distRTCorr_gausSigma[ch]);
    outTree -> Branch(Form("CTR_distRTCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_distRTCorr_gausSigmaErr[ch]);
  }
  
  std::map<std::string,TProfile*> p_eff_vs_X;
  std::map<std::string,TProfile*> p_eff_vs_Y;
  std::map<std::string,TProfile2D*> p2_eff_vs_XY;

  std::map<std::string,TH1F*> h_amp;
  std::map<std::string,TH1F*> h_amp_cut;
  std::map<std::string,TProfile*> p_amp_vs_X;
  std::map<std::string,TProfile*> p_amp_vs_Y;
  std::map<std::string,TProfile2D*> p2_amp_vs_XY;
  
  std::map<std::string,TH1F*> h_RT;
  std::map<std::string,TH1F*> h_RT_cut;

  std::map<std::string,TH1F*> h_time;
  std::map<std::string,TH1F*> h_tot;
  std::map<std::string,TH2F*> h2_tot_vs_amp;
  
  std::map<std::string,TProfile*> p_time_vs_amp;
  std::map<std::string,TProfile*> p_time_vs_RT;
  std::map<std::string,TProfile*> p_time_vs_dist;
  std::map<std::string,TProfile2D*> p2_time_vs_distAmp;
  std::map<std::string,TProfile2D*> p2_time_vs_distRT;
  
  std::map<std::string,TProfile*> p_time_ampCorr_vs_amp;
  std::map<std::string,TProfile*> p_time_ampCorr_posCorr_vs_amp;
  
  std::map<std::string,TProfile*> p_time_RTCorr_vs_RT;
  std::map<std::string,TProfile*> p_time_RTCorr_posCorr_vs_RT;
  
  std::map<std::string,TProfile*> p_time_vs_X;
  std::map<std::string,TProfile*> p_time_cut_vs_X;
  std::map<std::string,TProfile*> p_time_vs_Y;
  std::map<std::string,TProfile*> p_time_cut_vs_Y;
  std::map<std::string,TProfile2D*> p2_time_vs_XY;
  std::map<std::string,TProfile*> p_time_ampCorr_vs_X;
  std::map<std::string,TProfile*> p_time_cut_ampCorr_vs_X;
  std::map<std::string,TProfile*> p_time_ampCorr_vs_Y;
  std::map<std::string,TProfile*> p_time_cut_ampCorr_vs_Y;
  std::map<std::string,TProfile2D*> p2_time_ampCorr_vs_XY;
  std::map<std::string,TProfile*> p_time_ampCorr_posCorr_vs_X;
  std::map<std::string,TProfile*> p_time_cut_ampCorr_posCorr_vs_X;
  std::map<std::string,TProfile*> p_time_ampCorr_posCorr_vs_Y;
  std::map<std::string,TProfile*> p_time_cut_ampCorr_posCorr_vs_Y;
  std::map<std::string,TProfile2D*> p2_time_ampCorr_posCorr_vs_XY;
  std::map<std::string,TProfile*> p_time_RTCorr_vs_X;
  std::map<std::string,TProfile*> p_time_cut_RTCorr_vs_X;
  std::map<std::string,TProfile*> p_time_RTCorr_vs_Y;
  std::map<std::string,TProfile*> p_time_cut_RTCorr_vs_Y;
  std::map<std::string,TProfile2D*> p2_time_RTCorr_vs_XY;
  std::map<std::string,TProfile*> p_time_RTCorr_posCorr_vs_X;
  std::map<std::string,TProfile*> p_time_cut_RTCorr_posCorr_vs_X;
  std::map<std::string,TProfile*> p_time_RTCorr_posCorr_vs_Y;
  std::map<std::string,TProfile*> p_time_cut_RTCorr_posCorr_vs_Y;
  std::map<std::string,TProfile2D*> p2_time_RTCorr_posCorr_vs_XY;
  std::map<std::string,TProfile*> p_time_distAmpCorr_vs_X;
  std::map<std::string,TProfile*> p_time_cut_distAmpCorr_vs_X;
  std::map<std::string,TProfile*> p_time_distAmpCorr_vs_Y;
  std::map<std::string,TProfile*> p_time_cut_distAmpCorr_vs_Y;
  std::map<std::string,TProfile2D*> p2_time_distAmpCorr_vs_XY;
  
  std::map<std::string,TProfile*> p_time_distAmpCorr_vs_amp;
  std::map<std::string,TProfile*> p_time_distAmpCorr_vs_dist;
  
  std::map<std::string,TProfile*> p_time_distRTCorr_vs_amp;
  std::map<std::string,TProfile*> p_time_distRTCorr_vs_RT;
  std::map<std::string,TProfile*> p_time_distRTCorr_vs_dist;
  
  std::map<std::string,TH1F*> h_chessOccupancyX;
  std::map<std::string,TH1F*> h_chessOccupancyY;
  std::map<std::string,TH2F*> h2_chessOccupancyXY;
  std::map<std::string,TProfile*> map_p_time_vs_amp;
  
  std::map<std::string,TH1F*> h_CTR_raw;
  
  for(auto ch : channels)
  {
    float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
    float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));
    float rtMin = opts.GetOpt<float>(Form("%s.rtMin",ch.c_str())); 
    float rtMax = opts.GetOpt<float>(Form("%s.rtMax",ch.c_str()));
    float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str()));
    float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
    float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
    float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
    
    p_eff_vs_X[ch] = new TProfile(Form("p_eff_vs_X_%s",ch.c_str()),"",80,-20.,20.);
    p_eff_vs_Y[ch] = new TProfile(Form("p_eff_vs_Y_%s",ch.c_str()),"",80,-20.,20.);
    p2_eff_vs_XY[ch] = new TProfile2D(Form("p2_eff_vs_XY_%s",ch.c_str()),"",80,-20.,20.,80,-20.,20.);
    
    h_amp[ch]     = new TH1F(Form("h_amp_%s",    ch.c_str()),"",4000,0.,1.);
    h_amp_cut[ch] = new TH1F(Form("h_amp_cut_%s",ch.c_str()),"",4000,0.,1.);
    p_amp_vs_X[ch] = new TProfile(Form("p_amp_vs_X_%s",ch.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
    p_amp_vs_Y[ch] = new TProfile(Form("p_amp_vs_Y_%s",ch.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p2_amp_vs_XY[ch] = new TProfile2D(Form("p2_amp_vs_XY_%s",ch.c_str()),"",int(1*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax,int(1*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    
    h_RT[ch]     = new TH1F(Form("h_RT_%s",    ch.c_str()),"",2500,0.,100.);
    h_RT_cut[ch] = new TH1F(Form("h_RT_cut_%s",ch.c_str()),"",2500,0.,100.);
    
    h_time[ch] = new TH1F(Form("h_time_%s",ch.c_str()),"",200,0.,200.);
    h_tot[ch] = new TH1F(Form("h_tot_%s",ch.c_str()),"",100,0.,200.);
    h2_tot_vs_amp[ch] = new TH2F(Form("h_tot_vs_amp_%s",ch.c_str()),"",500,0.,1.,200,0.,200.);
    
    p_time_vs_amp[ch] = new TProfile(Form("p_time_vs_amp_%s",ch.c_str()),"",100,ampMin,ampMax);
    p_time_ampCorr_vs_amp[ch] = new TProfile(Form("p_time_ampCorr_vs_amp_%s",ch.c_str()),"",100,ampMin,ampMax);
    p_time_ampCorr_posCorr_vs_amp[ch] = new TProfile(Form("p_time_ampCorr_posCorr_vs_amp_%s",ch.c_str()),"",100,ampMin,ampMax);
    
    p_time_vs_RT[ch] = new TProfile(Form("p_time_vs_RT_%s",ch.c_str()),"",250,0.,100.);
    p_time_RTCorr_vs_RT[ch] = new TProfile(Form("p_time_RTCorr_vs_RT_%s",ch.c_str()),"",250,0.,100.);
    p_time_RTCorr_posCorr_vs_RT[ch] = new TProfile(Form("p_time_RTCorr_posCorr_vs_RT_%s",ch.c_str()),"",250,0.,100.);
    
    p_time_vs_X[ch] = new TProfile(Form("p_time_vs_X_%s",ch.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
    p_time_cut_vs_X[ch] = new TProfile(Form("p_time_cut_vs_X_%s",ch.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
    p_time_vs_Y[ch] = new TProfile(Form("p_time_vs_Y_%s",ch.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p_time_cut_vs_Y[ch] = new TProfile(Form("p_time_cut_vs_Y_%s",ch.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p2_time_vs_XY[ch] = new TProfile2D(Form("p2_time_vs_XY_%s",ch.c_str()),"",int(1*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax,int(1*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p_time_ampCorr_vs_X[ch] = new TProfile(Form("p_time_ampCorr_vs_X_%s",ch.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
    p_time_cut_ampCorr_vs_X[ch] = new TProfile(Form("p_time_cut_ampCorr_vs_X_%s",ch.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
    p_time_ampCorr_vs_Y[ch] = new TProfile(Form("p_time_ampCorr_vs_Y_%s",ch.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p_time_cut_ampCorr_vs_Y[ch] = new TProfile(Form("p_time_cut_ampCorr_vs_Y_%s",ch.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p2_time_ampCorr_vs_XY[ch] = new TProfile2D(Form("p2_time_ampCorr_vs_XY_%s",ch.c_str()),"",int(1*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax,int(1*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p_time_ampCorr_posCorr_vs_X[ch] = new TProfile(Form("p_time_ampCorr_posCorr_vs_X_%s",ch.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
    p_time_cut_ampCorr_posCorr_vs_X[ch] = new TProfile(Form("p_time_cut_ampCorr_posCorr_vs_X_%s",ch.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
    p_time_ampCorr_posCorr_vs_Y[ch] = new TProfile(Form("p_time_ampCorr_posCorr_vs_Y_%s",ch.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p_time_cut_ampCorr_posCorr_vs_Y[ch] = new TProfile(Form("p_time_cut_ampCorr_posCorr_vs_Y_%s",ch.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p2_time_ampCorr_posCorr_vs_XY[ch] = new TProfile2D(Form("p2_time_ampCorr_posCorr_vs_XY_%s",ch.c_str()),"",int(1*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax,int(1*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p_time_RTCorr_vs_X[ch] = new TProfile(Form("p_time_RTCorr_vs_X_%s",ch.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
    p_time_cut_RTCorr_vs_X[ch] = new TProfile(Form("p_time_cut_RTCorr_vs_X_%s",ch.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
    p_time_RTCorr_vs_Y[ch] = new TProfile(Form("p_time_RTCorr_vs_Y_%s",ch.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p_time_cut_RTCorr_vs_Y[ch] = new TProfile(Form("p_time_cut_RTCorr_vs_Y_%s",ch.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p2_time_RTCorr_vs_XY[ch] = new TProfile2D(Form("p2_time_RTCorr_vs_XY_%s",ch.c_str()),"",int(1*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax,int(1*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p_time_RTCorr_posCorr_vs_X[ch] = new TProfile(Form("p_time_RTCorr_posCorr_vs_X_%s",ch.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
    p_time_cut_RTCorr_posCorr_vs_X[ch] = new TProfile(Form("p_time_cut_RTCorr_posCorr_vs_X_%s",ch.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
    p_time_RTCorr_posCorr_vs_Y[ch] = new TProfile(Form("p_time_RTCorr_posCorr_vs_Y_%s",ch.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p_time_cut_RTCorr_posCorr_vs_Y[ch] = new TProfile(Form("p_time_cut_RTCorr_posCorr_vs_Y_%s",ch.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p2_time_RTCorr_posCorr_vs_XY[ch] = new TProfile2D(Form("p2_time_RTCorr_posCorr_vs_XY_%s",ch.c_str()),"",int(1*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax,int(1*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p_time_distAmpCorr_vs_X[ch] = new TProfile(Form("p_time_distAmpCorr_vs_X_%s",ch.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
    p_time_cut_distAmpCorr_vs_X[ch] = new TProfile(Form("p_time_cut_distAmpCorr_vs_X_%s",ch.c_str()),"",int(2*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax);
    p_time_distAmpCorr_vs_Y[ch] = new TProfile(Form("p_time_distAmpCorr_vs_Y_%s",ch.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p_time_cut_distAmpCorr_vs_Y[ch] = new TProfile(Form("p_time_cut_distAmpCorr_vs_Y_%s",ch.c_str()),"",int(2*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    p2_time_distAmpCorr_vs_XY[ch] = new TProfile2D(Form("p2_time_distAmpCorr_vs_XY_%s",ch.c_str()),"",int(1*(xtalXMax-xtalXMin)),xtalXMin,xtalXMax,int(1*(xtalYMax-xtalYMin)),xtalYMin,xtalYMax);
    
    p2_time_vs_distAmp[ch] = new TProfile2D(Form("p2_time_vs_distAmp_%s",ch.c_str()),"",20,0.,10.,20,ampMin,ampMax);
    p_time_distAmpCorr_vs_amp[ch] = new TProfile(Form("p_time_distAmpCorr_vs_amp_%s",ch.c_str()),"",100,ampMin,ampMax);
    p_time_distAmpCorr_vs_dist[ch] = new TProfile(Form("p_time_distAmpCorr_vs_dist_%s",ch.c_str()),"",100,0.,10.);
    
    p2_time_vs_distRT[ch] = new TProfile2D(Form("p2_time_vs_distRT_%s",ch.c_str()),"",20,0.,10.,20,rtMin,rtMax);
    p_time_distRTCorr_vs_RT[ch] = new TProfile(Form("p_time_distRTCorr_vs_RT_%s",ch.c_str()),"",100,rtMin,rtMax);
    p_time_distRTCorr_vs_amp[ch] = new TProfile(Form("p_time_distRTCorr_vs_amp_%s",ch.c_str()),"",100,ampMin,ampMax);
    p_time_distRTCorr_vs_dist[ch] = new TProfile(Form("p_time_distRTCorr_vs_dist_%s",ch.c_str()),"",100,0.,10.);
    
    h_CTR_raw[ch] = new TH1F(Form("h_CTR_raw_%s",ch.c_str()),"",10000,-50.,50.);
    
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
  
  
  
  //-----------------------
  // 1st loop over events
  if( maxEntries > 0 ) nEntries = maxEntries;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    h4 -> GetEntry(entry);
    // if( !goodSpills[std::make_pair(tv.run,tv.spill)] ) continue;

    
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
        if( (tv.nFibresOnX[hodoPlane] < nFibresMin || tv.nFibresOnX[hodoPlane] > nFibresMax) ) continue;
        if( (tv.nFibresOnY[hodoPlane] < nFibresMin || tv.nFibresOnY[hodoPlane] > nFibresMax) ) continue;
        if( fabs((tv.hodoX[1]-tv.hodoX[0]-3.50)) > 2. ) continue;
        if( fabs((tv.hodoY[1]-tv.hodoY[0]-0.12)) > 2. ) continue;
      }
      
      // fill amplitude and time plots
      int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
      std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
      std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
      std::string trgCh  = opts.GetOpt<std::string>(Form("%s.trgCh",ch.c_str()));
      std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
      
      if( ampCh == "NULL" ) continue;
      float amp = tv.amp_max[tv.channelIds[ampCh]] / 4096.;
      h_amp[ch] -> Fill( amp );
      
      float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
      float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));
      float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str())); 
      float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
      float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
      float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
      
      if( amp > ampMin )
      {
        if( tv.beamY > xtalYMin && tv.beamY < xtalYMax) p_eff_vs_X[ch] -> Fill( tv.beamX,1. );
        if( tv.beamX > xtalXMin && tv.beamX < xtalXMax) p_eff_vs_Y[ch] -> Fill( tv.beamY,1. );
        p2_eff_vs_XY[ch] -> Fill( tv.beamX,tv.beamY,1. );
        
        if( tv.beamY > xtalYMin && tv.beamY < xtalYMax) p_amp_vs_X[ch] -> Fill( tv.beamX,amp );
        if( tv.beamX > xtalXMin && tv.beamX < xtalXMax) p_amp_vs_Y[ch] -> Fill( tv.beamY,amp );
        p2_amp_vs_XY[ch] -> Fill( tv.beamX,tv.beamY,amp );
      }
      else
      {
        if( tv.beamY > xtalYMin && tv.beamY < xtalYMax) p_eff_vs_X[ch] -> Fill( tv.beamX,0. );
        if( tv.beamX > xtalXMin && tv.beamX < xtalXMax) p_eff_vs_Y[ch] -> Fill( tv.beamY,0. );
        p2_eff_vs_XY[ch] -> Fill( tv.beamX,tv.beamY,0. );
      }
      
      if( amp < ampMin || amp > ampMax ) continue;
      if( isnan(amp) ) continue;
      
      h_amp_cut[ch] -> Fill( amp );
      
      if( timeCh == "NULL" ) continue;
      float timTrg = trgCh != "NULL" ? tv.time[tv.channelIds[trgCh]+tv.timeMethodIds["LED"]] : 0;
      float tim  = tv.time[tv.channelIds[timeCh]+tv.timeMethodIds[timeMethods.at(0)]];
      float tim2 = timeMethods.at(1) != "NULL" ? tv.time[tv.channelIds[timeCh]+tv.timeMethodIds[timeMethods.at(1)]] : 0.;
      h_time[ch] -> Fill( tim );

      float timeMin = opts.GetOpt<float>(Form("%s.timeMin",ch.c_str())); 
      float timeMax = opts.GetOpt<float>(Form("%s.timeMax",ch.c_str()));
      if( tim > timeMin && tim < timeMax )
      {
        h_tot[ch] -> Fill( tim2-tim );      
        h2_tot_vs_amp[ch] -> Fill( amp,(tim2-tim) );
      }
      
      if( index < 0 ) continue;
      
      
      float tim50 = timeMethods.at(2) != "NULL" ? tv.time[tv.channelIds[ampCh]+tv.timeMethodIds[timeMethods.at(2)]] : 0.;
      float tim1000 = timeMethods.at(3) != "NULL" ? tv.time[tv.channelIds[ampCh]+tv.timeMethodIds[timeMethods.at(3)]] : 0.;
      h_RT[ch] -> Fill( tim1000-tim50 );
      
      // get reference channel
      std::string refCh  = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
      std::string ampChRef  = opts.GetOpt<std::string>(Form("%s.ampCh", refCh.c_str()));
      std::string timeChRef = opts.GetOpt<std::string>(Form("%s.timeCh",refCh.c_str()));
      std::string trgChRef  = opts.GetOpt<std::string>(Form("%s.trgCh",refCh.c_str()));
      std::vector<std::string> timeMethodsRef = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",refCh.c_str()));
      float ampMinRef = opts.GetOpt<float>(Form("%s.ampMin",refCh.c_str())); 
      float ampMaxRef = opts.GetOpt<float>(Form("%s.ampMax",refCh.c_str()));
      float timeMinRef = opts.GetOpt<float>(Form("%s.timeMin",refCh.c_str())); 
      float timeMaxRef = opts.GetOpt<float>(Form("%s.timeMax",refCh.c_str()));
      
      float ampRef = tv.amp_max[tv.channelIds[ampChRef]] / 4096.;
      float timTrgRef = trgChRef != "NULL" ? tv.time[tv.channelIds[trgChRef]+tv.timeMethodIds["LED"]] : 0;
      float timRef = tv.time[tv.channelIds[timeChRef]+tv.timeMethodIds[timeMethodsRef.at(0)]];
      
      if( ampRef < ampMinRef || ampRef > ampMaxRef ) continue;
      if( isnan(ampRef) ) continue;
      
      if( isinf(timRef) ) continue;
      if( isnan(timRef) ) continue;
      if( timRef < timeMinRef || timRef > timeMaxRef ) continue;
      
      if( isinf(tim) ) continue;
      if( isnan(tim) ) continue;
      if( tim < timeMin || tim > timeMax ) continue;
      
      h_RT_cut[ch] -> Fill( tim1000-tim50 );
      h_CTR_raw[ch] -> Fill( (tim-timTrg)-(timRef-timTrgRef) );
    }
  }
  std::cout << "\n>>> end 1st loop" << std::endl;
  
  
  
  //--- draw 1st plots
  for(auto ch: channels)
  {
    float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
    float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));
    float timeMin = opts.GetOpt<float>(Form("%s.timeMin",ch.c_str())); 
    float timeMax = opts.GetOpt<float>(Form("%s.timeMax",ch.c_str()));
    float totMin = opts.GetOpt<float>(Form("%s.totMin",ch.c_str())); 
    float totMax = opts.GetOpt<float>(Form("%s.totMax",ch.c_str()));
    float rtMin = opts.GetOpt<float>(Form("%s.rtMin",ch.c_str())); 
    float rtMax = opts.GetOpt<float>(Form("%s.rtMax",ch.c_str()));
    float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str())); 
    float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
    float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
    float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
    float ampLow = h_amp_cut[ch]->GetMean()-2.*h_amp_cut[ch]->GetRMS();
    float ampHig = h_amp_cut[ch]->GetMean()+2.*h_amp_cut[ch]->GetRMS();
    
    TLine* line_min_amp = new TLine(ampMin,h_amp[ch]->GetMinimum(),ampMin,h_amp[ch]->GetMaximum());
    TLine* line_max_amp = new TLine(ampMax,h_amp[ch]->GetMinimum(),ampMax,h_amp[ch]->GetMaximum());
    std::vector<TLine*> lines_amp;
    lines_amp.push_back(line_min_amp);
    lines_amp.push_back(line_max_amp);
    
    TLine* line_min_time = new TLine(timeMin,h_time[ch]->GetMinimum(),timeMin,h_time[ch]->GetMaximum());
    TLine* line_max_time = new TLine(timeMax,h_time[ch]->GetMinimum(),timeMax,h_time[ch]->GetMaximum());
    std::vector<TLine*> lines_time;
    lines_time.push_back(line_min_time);
    lines_time.push_back(line_max_time);
    
    TLine* line_min_tot = new TLine(totMin,h_tot[ch]->GetMinimum(),totMin,h_tot[ch]->GetMaximum());
    TLine* line_max_tot = new TLine(totMax,h_tot[ch]->GetMinimum(),totMax,h_tot[ch]->GetMaximum());
    std::vector<TLine*> lines_tot;
    lines_tot.push_back(line_min_tot);
    lines_tot.push_back(line_max_tot);
    
    TLine* line_min_rt = new TLine(rtMin,h_RT[ch]->GetMinimum(),rtMin,h_RT[ch]->GetMaximum());
    TLine* line_max_rt = new TLine(rtMax,h_RT[ch]->GetMinimum(),rtMax,h_RT[ch]->GetMaximum());
    std::vector<TLine*> lines_rt;
    lines_rt.push_back(line_min_rt);
    lines_rt.push_back(line_max_rt);
    
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
    
    
    c = new TCanvas(); DrawHistogram(opts,ch,h_amp[ch],     ";amplitude [V];entries", 0.,1.1*ampMax,4,true,latexLabels[ch],&lines_amp);  PrintCanvas(c,opts,ch,plotDir,"1_amp");      delete c;
    c = new TCanvas(); DrawHistogram(opts,ch,h_time[ch],    ";time [ns];entries",     0.,200.,      1,true,latexLabels[ch],&lines_time); PrintCanvas(c,opts,ch,plotDir,"1_time");     delete c;
    c = new TCanvas(); DrawHistogram(opts,ch,h_tot[ch],     ";ToT [ns];entries",      0.,200.,      1,true,latexLabels[ch],&lines_tot);  PrintCanvas(c,opts,ch,plotDir,"1_tot");      delete c;
    c = new TCanvas(); DrawHistogram(opts,ch,h_RT[ch],      ";rise time [ns];entries",0.,100.,      1,true,latexLabels[ch],&lines_rt);   PrintCanvas(c,opts,ch,plotDir,"1_RT");       delete c;
    
    c = new TCanvas();
    DrawHistogram2D(opts,ch,h2_tot_vs_amp[ch],";amplitude [V];ToT [ns];entries",ampMin,ampMax,80.,160.,h2_tot_vs_amp[ch]->GetMinimum(),h2_tot_vs_amp[ch]->GetMaximum(),false,false,"colz",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,"1_tot_vs_amp");
    delete c;

    c = new TCanvas("c","c",2000,600);
    c -> Divide(3,1);
    c->cd(1); DrawProfile  (opts,ch,p_eff_vs_X[ch],  ";hodoscope x [mm];efficiency",                 -30.,30.,0.001,1.1,kBlack,"",latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_eff_vs_Y[ch],  ";hodoscope y [mm];efficiency",                 -30.,30.,0.001,1.1,kBlack,"",latexLabels[ch]);
    c->cd(3); DrawProfile2D(opts,ch,p2_eff_vs_XY[ch],";hodoscope x [mm];hodoscope y [mm];efficiency",-30.,30.,-30.,30.,0.001,1.1,latexLabels[ch],&lines_xtal);
    PrintCanvas(c,opts,ch,plotDir,"1_eff_vs_XY"); delete c;

    c = new TCanvas("c","c",2000,600);
    c -> Divide(3,1);
    c->cd(1); DrawProfile  (opts,ch,p_amp_vs_X[ch],  ";hodoscope x [mm];amplitude [V]",                 xtalXMin,xtalXMax,ampLow,ampHig,kBlack,"",latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_amp_vs_Y[ch],  ";hodoscope y [mm];amplitude [V]",                 xtalYMin,xtalYMax,ampLow,ampHig,kBlack,"",latexLabels[ch]);
    c->cd(3); DrawProfile2D(opts,ch,p2_amp_vs_XY[ch],";hodoscope x [mm];hodoscope y [mm];amplitude [V]",xtalXMin,xtalXMax,xtalYMin,xtalYMax,ampLow,ampHig,latexLabels[ch]);    
    PrintCanvas(c,opts,ch,plotDir,"1_amp_vs_XY"); delete c;
    
    // if( index > 0 )
    // {
    //   TF1* fitFunc_landau = new TF1(Form("fitFunc_landau_%s",shortLabel.c_str()),"[0]*TMath::Landau(x,[1],[2])",ampMin,ampMax);
    //   fitFunc_landau -> SetParameters(0.01,histo->GetMean(),0.1);
    //   histo -> Fit(fitFunc_landau,"QNRS+");
    //   fitFunc_landau -> SetLineColor(kRed+1);
    //   fitFunc_landau -> SetLineWidth(2);
    //   fitFunc_landau -> Draw("same");

    //   TLatex* latexLandau = new TLatex(0.25,0.80,Form("Landau fit MPV: %.2f V",fitFunc_landau->GetParameter(1)));
    //   latexLandau -> SetNDC();
    //   latexLandau -> SetTextFont(42);
    //   latexLandau -> SetTextSize(0.04);
    //   latexLandau -> SetTextColor(kRed+1);
    //   latexLandau -> Draw("same");
    // }
  }


  
  //--- find CTR ranges
  std::map<std::string,float> CTRRanges_mean;
  std::map<std::string,float> CTRRanges_sigma;
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
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
      if( !AcceptEvent(av,tv,opts,ch,nFibresMin,nFibresMax,(CTRRanges_mean[ch]-4.*CTRRanges_sigma[ch]),(CTRRanges_mean[ch]+4.*CTRRanges_sigma[ch])) ) continue;

      std::string Vbias = opts.GetOpt<std::string>(Form("%s.Vbias", ch.c_str()));
      std::string NINOthr = opts.GetOpt<std::string>(Form("%s.NINOthr", ch.c_str()));
      t_Vbias[ch] = tv.VbiasVals[Vbias];
      t_NINOthr[ch] = tv.NINOthrVals[NINOthr];
      
      p_time_vs_amp[ch] -> Fill( av.amp,av.time-av.timeRef );
      p_time_vs_RT[ch] -> Fill( av.rt,av.time-av.timeRef );
      p2_time_vs_distAmp[ch] -> Fill( av.dist,av.amp,av.time-av.timeRef );
      p2_time_vs_distRT[ch] -> Fill( av.dist,av.rt,av.time-av.timeRef );
      
      p_time_vs_X[ch] -> Fill( tv.beamX,av.time-av.timeRef );
      if( (tv.beamY >= 0.5*(xtalYMin+xtalYMax)-1.5) && (tv.beamY <= 0.5*(xtalYMin+xtalYMax)+1.5) )
        p_time_cut_vs_X[ch] -> Fill( tv.beamX,av.time-av.timeRef );
      p_time_vs_Y[ch] -> Fill( tv.beamY,av.time-av.timeRef );
      if( (tv.beamX >= 0.5*(xtalXMin+xtalXMax)-1.5) && (tv.beamX <= 0.5*(xtalXMin+xtalXMax)+1.5) )
        p_time_cut_vs_Y[ch] -> Fill( tv.beamY,av.time-av.timeRef );
      p2_time_vs_XY[ch] -> Fill( tv.beamX,tv.beamY,av.time-av.timeRef );

      int binX = h_chessOccupancyX[ch] -> Fill(tv.beamX);
      int binY = h_chessOccupancyY[ch] -> Fill(tv.beamY);
      h2_chessOccupancyXY[ch] -> Fill(tv.beamX,tv.beamY);
      
      std::string label(Form("%s_%d-%d",ch.c_str(),binX-1,binY-1));
      map_p_time_vs_amp[label] -> Fill( av.amp,av.time-av.timeRef );
    }
  }
  std::cout << "\n>>> end 2nd loop" << std::endl;
  
  
  
  //--- draw 2nd plots
  std::map<std::string,TF1*> fitFunc_corrAmp;
  std::map<std::string,TF1*> fitFunc_corrRT;
  std::map<std::string,TF1*> fitFunc_corrDist;
  
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
    float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));
    float rtMin = opts.GetOpt<float>(Form("%s.rtMin",ch.c_str())); 
    float rtMax = opts.GetOpt<float>(Form("%s.rtMax",ch.c_str()));
    float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str())); 
    float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
    float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
    float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
    float timeLow = CTRRanges_mean[ch]-3.*CTRRanges_sigma[ch];
    float timeHig = CTRRanges_mean[ch]+3.*CTRRanges_sigma[ch];
    
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
    
    
    prof = p_time_vs_amp[ch];
    
    std::string fitFunc = opts.GetOpt<std::string>(Form("%s.fitFunc",ch.c_str()));
    float fitXMin = opts.GetOpt<float>(Form("%s.fitFunc",ch.c_str()),1);
    float fitXMax = opts.GetOpt<float>(Form("%s.fitFunc",ch.c_str()),2);
    int nPar = opts.OptExist(Form("%s.fitFunc",ch.c_str()),3) ? opts.GetOpt<int>(Form("%s.fitFunc",ch.c_str()),3) : 0;
    
    fitFunc_corrAmp[ch] = new TF1(Form("fitFunc_corrAmp_%s",ch.c_str()),fitFunc.c_str(),ampMin,ampMax);
    for(int iPar = 0; iPar < nPar; ++iPar)
      fitFunc_corrAmp[ch] -> SetParameter(iPar,opts.GetOpt<float>(Form("%s.fitFunc",ch.c_str()),4+iPar));
    
    prof -> Fit(fitFunc_corrAmp[ch],"QNRS+","",fitXMin,fitXMax);
    
    c = new TCanvas();
    DrawProfile(opts,ch,p_time_vs_amp[ch],";amplitude [V];#Delta_{t} [ns]",ampMin,ampMax,timeLow,timeHig,kBlack,"",latexLabels[ch],fitFunc_corrAmp[ch]);
    PrintCanvas(c,opts,ch,plotDir,"2_time_vs_amp");
    delete c;
    
    
    prof = p_time_vs_RT[ch];
    
    std::string fitFuncRT = opts.GetOpt<std::string>(Form("%s.fitFuncRT",ch.c_str()));
    float fitXMinRT = opts.GetOpt<float>(Form("%s.fitFuncRT",ch.c_str()),1);
    float fitXMaxRT = opts.GetOpt<float>(Form("%s.fitFuncRT",ch.c_str()),2);
    int nParRT = opts.OptExist(Form("%s.fitFuncRT",ch.c_str()),3) ? opts.GetOpt<int>(Form("%s.fitFuncRT",ch.c_str()),3) : 0;
    
    fitFunc_corrRT[ch] = new TF1(Form("fitFunc_corrRT_%s",ch.c_str()),fitFuncRT.c_str(),rtMin,rtMax);
    for(int iPar = 0; iPar < nParRT; ++iPar)
      fitFunc_corrRT[ch] -> SetParameter(iPar,opts.GetOpt<float>(Form("%s.fitFuncRT",ch.c_str()),4+iPar));
    
    prof -> Fit(fitFunc_corrRT[ch],"QNRS+","",fitXMinRT,fitXMaxRT);
    
    c = new TCanvas();
    DrawProfile(opts,ch,p_time_vs_RT[ch],";rise time [ns];#Delta_{t} [ns]",rtMin,rtMax,timeLow,timeHig,kBlack,"",latexLabels[ch],fitFunc_corrRT[ch]);
    PrintCanvas(c,opts,ch,plotDir,"2_time_vs_RT");
    delete c;
    
    
    int nBinsX = opts.GetOpt<int>(Form("%s.nBinsX",ch.c_str()));
    int nBinsY = opts.GetOpt<int>(Form("%s.nBinsY",ch.c_str()));
    for(int binX = 0; binX < nBinsX; ++binX)
      for(int binY = 0; binY < nBinsY; ++binY)
      {
        std::string label(Form("%s_%d-%d",ch.c_str(),binX,binY));
        prof = map_p_time_vs_amp[label];
        
        fitFunc_corrAmp[label] = new TF1(Form("fitFunc_corrAmp_%s",label.c_str()),fitFunc.c_str(),ampMin,ampMax);
        for(int iPar = 0; iPar < nPar; ++iPar)
          fitFunc_corrAmp[label] -> SetParameter(iPar,opts.GetOpt<float>(Form("%s.fitFunc",label.c_str()),4+iPar));
        
        prof -> Fit(fitFunc_corrAmp[label],"QNRS+","",fitXMin,fitXMax);        
      }
    
    
    c = new TCanvas();
    DrawProfile2D(opts,ch,p2_time_vs_distAmp[ch],";distance from center [mm];amplitude [V];#Delta_{t} [ns]",0.,10.,ampMin,ampMax,timeLow,timeHig,latexLabels[ch]); 
    PrintCanvas(c,opts,ch,plotDir,"2_time_vs_distAmp");
    delete c;

    
    c = new TCanvas("c","c",2000,600);
    c -> Divide(3,1);
    c->cd(1); DrawProfile  (opts,ch,p_time_vs_X[ch],    ";hodoscope x [mm];#Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(1); DrawProfile  (opts,ch,p_time_cut_vs_X[ch],";hodoscope x [mm];#Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kRed,  "same",latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_vs_Y[ch],    ";hodoscope y [mm];#Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_cut_vs_Y[ch],";hodoscope y [mm];#Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlue, "same",latexLabels[ch]);
    c->cd(3); DrawProfile2D(opts,ch,p2_time_vs_XY[ch],  ";hodoscope x [mm];hodoscope y [mm];#Delta_{t} [ns]",xtalXMin,xtalXMax,xtalYMin,xtalYMax,timeLow,timeHig,latexLabels[ch],&lines_xtal);
    PrintCanvas(c,opts,ch,plotDir,"2_time_vs_XY");
    delete c;
    
    
    c = new TCanvas();
    DrawProfile2D(opts,ch,p2_time_vs_distRT[ch],";distance from center [mm];rise time [ns];#Delta_{t} [ns]",0.,10.,rtMin,rtMax,timeLow,timeHig,latexLabels[ch]); 
    PrintCanvas(c,opts,ch,plotDir,"2_time_vs_distRT");
    delete c;
    
    
    c = new TCanvas();
    DrawHistogram2D(opts,ch,h2_chessOccupancyXY[ch],";hodoscope x [mm];hodoscope y [mm];entries",xtalXMin,xtalXMax,xtalYMin,xtalYMax,h2_chessOccupancyXY[ch]->GetMinimum(),h2_chessOccupancyXY[ch]->GetMaximum(),false,false,"colz",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,"2_chessOccupancyXY");
    delete c;    
  }
  
  
  
  //------------------------
  //--- 3rd loop over events
  std::map<std::string,TH1F*> h_CTR;
  std::map<std::string,TH1F*> h_CTR_ampCorr;
  std::map<std::string,TH1F*> h_CTR_RTCorr;
  std::map<std::string,TH1F*> h_CTR_distAmpCorr;
  std::map<std::string,TH1F*> h_CTR_distRTCorr;
  
  std::map<std::string,TH1F*> map_h_CTR_ampCorr;
  
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    h_CTR[ch]         = new TH1F(Form("h_CTR_%s",ch.c_str()),        "",1000,CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+5.*CTRRanges_sigma[ch]);
    h_CTR_ampCorr[ch] = new TH1F(Form("h_CTR_ampCorr_%s",ch.c_str()),"",1000,CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+5.*CTRRanges_sigma[ch]);
    h_CTR_RTCorr[ch]  = new TH1F(Form("h_CTR_RTCorr_%s",ch.c_str()), "",1000,CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+5.*CTRRanges_sigma[ch]);
    h_CTR_distAmpCorr[ch] = new TH1F(Form("h_CTR_distAmpCorr_%s",ch.c_str()),"",1000,CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+5.*CTRRanges_sigma[ch]);
    h_CTR_distRTCorr[ch] = new TH1F(Form("h_CTR_distRTCorr_%s",ch.c_str()),"",1000,CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+5.*CTRRanges_sigma[ch]);
    
    int nBinsX = opts.GetOpt<int>(Form("%s.nBinsX",ch.c_str()));
    int nBinsY = opts.GetOpt<int>(Form("%s.nBinsY",ch.c_str()));
    for(int binX = 0; binX < nBinsX; ++binX)
      for(int binY = 0; binY < nBinsY; ++binY)
      {
        std::string label(Form("%s_%d-%d",ch.c_str(),binX,binY));    
        map_h_CTR_ampCorr[label] = new TH1F(Form("map_CTR_ampCorr_%s",label.c_str()),"",1000,CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+5.*CTRRanges_sigma[ch]);
      }
  }
  
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 3rd loop: reading entry " << entry << " / " << nEntries << "\r" << std::flush;
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
      
      float CTR = av.time - av.timeRef;
      float CTR_ampCorr = CTR - fitFunc_corrAmp[ch]->Eval(av.amp) + fitFunc_corrAmp[ch]->Eval( h_amp_cut[ch]->GetMean() );
      float CTR_RTCorr = CTR - fitFunc_corrRT[ch]->Eval(av.rt) + fitFunc_corrRT[ch]->Eval( h_RT_cut[ch]->GetMean() );
      float CTR_distAmpCorr = CTR - p2_time_vs_distAmp[ch]->GetBinContent(p2_time_vs_distAmp[ch]->FindBin(av.dist,av.amp)) + p2_time_vs_distAmp[ch]->GetBinContent(p2_time_vs_distAmp[ch]->FindBin(3.,h_amp_cut[ch]->GetMean()));
      float CTR_distRTCorr = CTR - p2_time_vs_distRT[ch]->GetBinContent(p2_time_vs_distRT[ch]->FindBin(av.dist,av.rt)) + p2_time_vs_distRT[ch]->GetBinContent(p2_time_vs_distRT[ch]->FindBin(3.,h_RT_cut[ch]->GetMean()));
      
      p_time_ampCorr_vs_amp[ch] -> Fill( av.amp,CTR_ampCorr );
      p_time_ampCorr_vs_X[ch] -> Fill( tv.beamX,CTR_ampCorr );
      if( (tv.beamY >= 0.5*(xtalYMin+xtalYMax)-1.5) && (tv.beamY <= 0.5*(xtalYMin+xtalYMax)+1.5) )
        p_time_cut_ampCorr_vs_X[ch] -> Fill( tv.beamX,CTR_ampCorr );
      p_time_ampCorr_vs_Y[ch] -> Fill( tv.beamY,CTR_ampCorr );
      if( (tv.beamX >= 0.5*(xtalXMin+xtalXMax)-1.5) && (tv.beamX <= 0.5*(xtalXMin+xtalXMax)+1.5) )
        p_time_cut_ampCorr_vs_Y[ch] -> Fill( tv.beamY,CTR_ampCorr );
      p2_time_ampCorr_vs_XY[ch] -> Fill( tv.beamX,tv.beamY,CTR_ampCorr );

      p_time_RTCorr_vs_RT[ch] -> Fill( av.rt,CTR_RTCorr );
      p_time_RTCorr_vs_X[ch] -> Fill( tv.beamX,CTR_RTCorr );
      if( (tv.beamY >= 0.5*(xtalYMin+xtalYMax)-1.5) && (tv.beamY <= 0.5*(xtalYMin+xtalYMax)+1.5) )
        p_time_cut_RTCorr_vs_X[ch] -> Fill( tv.beamX,CTR_RTCorr );
      p_time_RTCorr_vs_Y[ch] -> Fill( tv.beamY,CTR_RTCorr );
      if( (tv.beamX >= 0.5*(xtalXMin+xtalXMax)-1.5) && (tv.beamX <= 0.5*(xtalXMin+xtalXMax)+1.5) )
        p_time_cut_RTCorr_vs_Y[ch] -> Fill( tv.beamY,CTR_RTCorr );
      p2_time_RTCorr_vs_XY[ch] -> Fill( tv.beamX,tv.beamY,CTR_RTCorr );
      
      p_time_distAmpCorr_vs_amp[ch] -> Fill( av.amp,CTR_distAmpCorr );
      p_time_distAmpCorr_vs_dist[ch] -> Fill( av.dist,CTR_distAmpCorr );
      p_time_distAmpCorr_vs_X[ch] -> Fill( tv.beamX,CTR_distAmpCorr );
      if( (tv.beamY >= 0.5*(xtalYMin+xtalYMax)-1.5) && (tv.beamY <= 0.5*(xtalYMin+xtalYMax)+1.5) )
        p_time_cut_distAmpCorr_vs_X[ch] -> Fill( tv.beamX,CTR_distAmpCorr );
      p_time_distAmpCorr_vs_Y[ch] -> Fill( tv.beamY,CTR_distAmpCorr );
      if( (tv.beamX >= 0.5*(xtalXMin+xtalXMax)-1.5) && (tv.beamX <= 0.5*(xtalXMin+xtalXMax)+1.5) )
        p_time_cut_distAmpCorr_vs_Y[ch] -> Fill( tv.beamY,CTR_distAmpCorr );
      p2_time_distAmpCorr_vs_XY[ch] -> Fill( tv.beamX,tv.beamY,CTR_distAmpCorr );

      p_time_distRTCorr_vs_RT[ch] -> Fill( av.rt,CTR_distRTCorr );
      p_time_distRTCorr_vs_amp[ch] -> Fill( av.amp,CTR_distRTCorr );
      p_time_distRTCorr_vs_dist[ch] -> Fill( av.dist,CTR_distRTCorr );
      
      h_CTR[ch] -> Fill( CTR );
      h_CTR_ampCorr[ch] -> Fill( CTR_ampCorr );
      h_CTR_RTCorr[ch] -> Fill( CTR_RTCorr );
      h_CTR_distAmpCorr[ch] -> Fill( CTR_distAmpCorr );
      h_CTR_distRTCorr[ch] -> Fill( CTR_distRTCorr );
      
      
      int binX = h_chessOccupancyX[ch] -> Fill(tv.beamX);
      int binY = h_chessOccupancyY[ch] -> Fill(tv.beamY);
      std::string label(Form("%s_%d-%d",ch.c_str(),binX-1,binY-1));
      // map_h_CTR_ampCorr[label] -> Fill( CTR-fitFunc_corrAmp[label]->Eval(av.amp)+fitFunc_corrAmp[label]->Eval(h_amp_cut[ch]->GetMean()) );
      map_h_CTR_ampCorr[label] -> Fill( CTR-fitFunc_corrAmp[ch]->Eval(av.amp)+fitFunc_corrAmp[ch]->Eval(h_amp_cut[ch]->GetMean()) );
    }
  }
  std::cout << "\n>>> end 3rd loop" << std::endl;
  
  
  
  //--- draw 3rd plots
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;

    float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
    float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));
    float rtMin = opts.GetOpt<float>(Form("%s.rtMin",ch.c_str())); 
    float rtMax = opts.GetOpt<float>(Form("%s.rtMax",ch.c_str()));
    float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str())); 
    float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
    float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
    float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
    float timeLow = CTRRanges_mean[ch]-3.*CTRRanges_sigma[ch];
    float timeHig = CTRRanges_mean[ch]+3.*CTRRanges_sigma[ch];
    
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

    
    c = new TCanvas();
    DrawProfile(opts,ch,p_time_RTCorr_vs_RT[ch],";rise time [ns];r.t.-corrected #Delta_{t} [ns]",rtMin,rtMax,timeLow,timeHig,kBlack,"",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,"3_time_RTCorr_vs_RT");
    delete c;
    
    c = new TCanvas("c","c",2000,600);
    c -> Divide(3,1);
    c->cd(1); DrawProfile  (opts,ch,p_time_RTCorr_vs_X[ch],    ";hodoscope x [mm];r.t.-corrected #Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(1); DrawProfile  (opts,ch,p_time_cut_RTCorr_vs_X[ch],";hodoscope x [mm];r.t.-corrected #Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kRed,  "same",latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_RTCorr_vs_Y[ch],    ";hodoscope y [mm];r.t.-corrected #Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_cut_RTCorr_vs_Y[ch],";hodoscope y [mm];r.t.-corrected #Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlue, "same",latexLabels[ch]);
    c->cd(3); DrawProfile2D(opts,ch,p2_time_RTCorr_vs_XY[ch],  ";hodoscope x [mm];hodoscope y [mm];r.t.-corrected #Delta_{t} [ns]",xtalXMin,xtalXMax,xtalYMin,xtalYMax,timeLow,timeHig,latexLabels[ch],&lines_xtal);
    PrintCanvas(c,opts,ch,plotDir,"3_time_RTCorr_vs_XY");
    delete c;
    
    
    c = new TCanvas();
    DrawProfile(opts,ch,p_time_distAmpCorr_vs_amp[ch],";amplitude [V];dist.&amp.-corrected #Delta_{t} [ns]",ampMin,ampMax,timeLow,timeHig,kBlack,"",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,"3_time_distAmpCorr_vs_amp");
    delete c;
    
    c = new TCanvas();    
    DrawProfile(opts,ch,p_time_distAmpCorr_vs_dist[ch],";distance from center [mm];dist&amp.-corrected #Delta_{t} [ns]",0.,10.,timeLow,timeHig,kBlack,"",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,"3_time_distAmpCorr_vs_dist");
    delete c;
    
    c = new TCanvas("c","c",2000,600);
    c -> Divide(3,1);
    c->cd(1); DrawProfile  (opts,ch,p_time_distAmpCorr_vs_X[ch],    ";hodoscope x [mm];dist.&amp.-corrected #Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(1); DrawProfile  (opts,ch,p_time_cut_distAmpCorr_vs_X[ch],";hodoscope x [mm];dist.&amp.-corrected #Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kRed,  "same",latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_distAmpCorr_vs_Y[ch],    ";hodoscope y [mm];dist.&amp.-corrected #Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_cut_distAmpCorr_vs_Y[ch],";hodoscope y [mm];dist.&amp.-corrected #Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlue, "same",latexLabels[ch]);
    c->cd(3); DrawProfile2D(opts,ch,p2_time_distAmpCorr_vs_XY[ch],  ";hodoscope x [mm];hodoscope y [mm];dist.&amp.-corrected #Delta_{t} [ns]",xtalXMin,xtalXMax,xtalYMin,xtalYMax,timeLow,timeHig,latexLabels[ch],&lines_xtal);
    PrintCanvas(c,opts,ch,plotDir,"3_time_distAmpCorr_vs_XY");
    delete c;
    
    
    c = new TCanvas();
    DrawProfile(opts,ch,p_time_distRTCorr_vs_RT[ch],";rise time [ns];dist.&r.t.-corrected #Delta_{t} [ns]",0.,100.,timeLow,timeHig,kBlack,"",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,"3_time_distRTCorr_vs_RT");
    delete c;
    
    c = new TCanvas();
    DrawProfile(opts,ch,p_time_distRTCorr_vs_amp[ch],";amplitude [V];dist.&r.t.-corrected #Delta_{t} [ns]",ampMin,ampMax,timeLow,timeHig,kBlack,"",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,"3_time_distRTCorr_vs_amp");
    delete c;
    
    c = new TCanvas();    
    DrawProfile(opts,ch,p_time_distRTCorr_vs_dist[ch],";distance from center [mm];dist&r.t.-corrected #Delta_{t} [ns]",0.,10.,timeLow,timeHig,kBlack,"",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,"3_time_distRTCorr_vs_dist");
    delete c;    
  }
  
  
  
  //------------------------
  //--- 4th loop over events
  int rebin = opts.GetOpt<int>("Input.rebin");
  TLatex* latexLabel12 = new TLatex(0.13,0.97,Form("%s",""));

  
  std::map<std::string,TH1F*> h_CTR_ampCorr_posCorr;
  std::map<std::string,TH1F*> h_CTR_RTCorr_posCorr;
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    if( index < 0 ) continue;
    
    h_CTR_ampCorr_posCorr[ch] = new TH1F(Form("h_CTR_ampCorr_posCorr_%s",ch.c_str()),"",1000,CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+5.*CTRRanges_sigma[ch]);
    h_CTR_RTCorr_posCorr[ch]  = new TH1F(Form("h_CTR_RTCorr_posCorr_%s",ch.c_str()), "",1000,CTRRanges_mean[ch]-5.*CTRRanges_sigma[ch],CTRRanges_mean[ch]+5.*CTRRanges_sigma[ch]);
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
      
      float CTR = av.time - av.timeRef;
      float CTR_ampCorr = CTR - fitFunc_corrAmp[ch]->Eval(av.amp) + fitFunc_corrAmp[ch]->Eval( h_amp_cut[ch]->GetMean() );
      float CTR_RTCorr = CTR - fitFunc_corrRT[ch]->Eval(av.rt) + fitFunc_corrRT[ch]->Eval( h_RT_cut[ch]->GetMean() );
      float CTR_ampCorr_posCorr = CTR_ampCorr - p2_time_ampCorr_vs_XY[ch]->GetBinContent(p2_time_ampCorr_vs_XY[ch]->FindBin(tv.beamX,tv.beamY))
                                              + p2_time_ampCorr_vs_XY[ch]->GetBinContent(p2_time_ampCorr_vs_XY[ch]->FindBin(0.5*(xtalXMin+xtalXMax)+2,0.5*(xtalYMin+xtalYMax)+2));
      float CTR_RTCorr_posCorr = CTR_RTCorr - p2_time_RTCorr_vs_XY[ch]->GetBinContent(p2_time_RTCorr_vs_XY[ch]->FindBin(tv.beamX,tv.beamY))
                                             + p2_time_RTCorr_vs_XY[ch]->GetBinContent(p2_time_RTCorr_vs_XY[ch]->FindBin(0.5*(xtalXMin+xtalXMax)+2,0.5*(xtalYMin+xtalYMax)+2));
      
      p_time_ampCorr_posCorr_vs_amp[ch] -> Fill( av.amp,CTR_ampCorr_posCorr );
      p_time_ampCorr_posCorr_vs_X[ch] -> Fill( tv.beamX,CTR_ampCorr_posCorr );
      p_time_ampCorr_posCorr_vs_Y[ch] -> Fill( tv.beamY,CTR_ampCorr_posCorr );
      p2_time_ampCorr_posCorr_vs_XY[ch] -> Fill( tv.beamX,tv.beamY,CTR_ampCorr_posCorr );
      
      p_time_RTCorr_posCorr_vs_RT[ch] -> Fill( av.rt,CTR_RTCorr_posCorr );
      p_time_RTCorr_posCorr_vs_X[ch] -> Fill( tv.beamX,CTR_RTCorr_posCorr );
      p_time_RTCorr_posCorr_vs_Y[ch] -> Fill( tv.beamY,CTR_RTCorr_posCorr );
      p2_time_RTCorr_posCorr_vs_XY[ch] -> Fill( tv.beamX,tv.beamY,CTR_RTCorr_posCorr );
      
      h_CTR_ampCorr_posCorr[ch] -> Fill( CTR_ampCorr_posCorr );
      h_CTR_RTCorr_posCorr[ch] -> Fill( CTR_RTCorr_posCorr );
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
    float timeLow = CTRRanges_mean[ch]-3.*CTRRanges_sigma[ch];
    float timeHig = CTRRanges_mean[ch]+3.*CTRRanges_sigma[ch];
    
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
    DrawProfile(opts,ch,p_time_ampCorr_posCorr_vs_amp[ch],";amplitude [V];amp.+pos.-corrected #Delta_{t} [ns]",ampMin,ampMax,timeLow,timeHig,kBlack,"",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,"4_time_ampCorr_posCorr_vs_amp");
    delete c;

    c = new TCanvas("c","c",2000,600);
    c -> Divide(3,1);
    c->cd(1); DrawProfile  (opts,ch,p_time_ampCorr_posCorr_vs_X[ch],    ";hodoscope x [mm];amp.+pos.-corrected #Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(1); DrawProfile  (opts,ch,p_time_cut_ampCorr_posCorr_vs_X[ch],";hodoscope x [mm];amp.+pos.-corrected #Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kRed,  "same",latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_ampCorr_posCorr_vs_Y[ch],    ";hodoscope y [mm];amp.+pos.-corrected #Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_cut_ampCorr_posCorr_vs_Y[ch],";hodoscope y [mm];amp.+pos.-corrected #Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlue, "same",latexLabels[ch]);
    c->cd(3); DrawProfile2D(opts,ch,p2_time_ampCorr_posCorr_vs_XY[ch],  ";hodoscope x [mm];hodoscope y [mm];amp.+pos.-corrected #Delta_{t} [ns]",xtalXMin,xtalXMax,xtalYMin,xtalYMax,timeLow,timeHig,latexLabels[ch],&lines_xtal);
    PrintCanvas(c,opts,ch,plotDir,"4_time_ampCorr_posCorr_vs_XY");
    delete c;
    
    
    c = new TCanvas();
    DrawProfile(opts,ch,p_time_RTCorr_posCorr_vs_RT[ch],";rise time [ns];r.t.+pos.-corrected #Delta_{t} [ns]",0.,40.,timeLow,timeHig,kBlack,"",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,"4_time_RTCorr_posCorr_vs_RT");
    delete c;
    
    c = new TCanvas("c","c",2000,600);
    c -> Divide(3,1);
    c->cd(1); DrawProfile  (opts,ch,p_time_RTCorr_posCorr_vs_X[ch],    ";hodoscope x [mm];r.t.+pos.-corrected #Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(1); DrawProfile  (opts,ch,p_time_cut_RTCorr_posCorr_vs_X[ch],";hodoscope x [mm];r.t.+pos.-corrected #Delta_{t} [ns]",                 xtalXMin,xtalXMax,timeLow,timeHig,kRed,  "same",latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_RTCorr_posCorr_vs_Y[ch],    ";hodoscope y [mm];r.t.+pos.-corrected #Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlack,"",    latexLabels[ch]);
    c->cd(2); DrawProfile  (opts,ch,p_time_cut_RTCorr_posCorr_vs_Y[ch],";hodoscope y [mm];r.t.+pos.-corrected #Delta_{t} [ns]",                 xtalYMin,xtalYMax,timeLow,timeHig,kBlue, "same",latexLabels[ch]);
    c->cd(3); DrawProfile2D(opts,ch,p2_time_RTCorr_posCorr_vs_XY[ch],  ";hodoscope x [mm];hodoscope y [mm];r.t.+pos.-corrected #Delta_{t} [ns]",xtalXMin,xtalXMax,xtalYMin,xtalYMax,timeLow,timeHig,latexLabels[ch],&lines_xtal);
    PrintCanvas(c,opts,ch,plotDir,"4_time_RTCorr_posCorr_vs_XY");
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
    
    histo = h_CTR_ampCorr[ch];
    
    if( intrinsic < 0 )
      ctrResults = drawCTRPlot(histo,"amp. walk corr.",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
    else
      ctrResults = drawCTRPlot(histo,"amp. walk corr.",rebin,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
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
    
    c = new TCanvas(Form("c_CTR_RTCorr_%s",shortLabel.c_str()),Form("c_CTR_%s",shortLabel.c_str()));
    
    histo = h_CTR_RTCorr[ch];
    
    if( intrinsic < 0 )
      ctrResults = drawCTRPlot(histo,"r.t. corr.",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
    else
      ctrResults = drawCTRPlot(histo,"r.t. corr.",rebin,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
    t_CTR_RTCorr_effSigma[ch] = ctrResults.effSigma;
    t_CTR_RTCorr_gausSigma[ch] = ctrResults.gausSigma;
    t_CTR_RTCorr_gausSigmaErr[ch] = ctrResults.gausSigmaErr;
    
    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();

    c -> Print(Form("%s/c_4_tRes_RTCorr_%s.png",plotDir.c_str(),shortLabel.c_str()));
    c -> Print(Form("%s/c_4_tRes_RTCorr_%s.pdf",plotDir.c_str(),shortLabel.c_str()));
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
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string shortLabel = opts.GetOpt<std::string>(Form("%s.shortLabel",ch.c_str()));
    float intrinsic = opts.GetOpt<float>(Form("%s.intrinsic",ch.c_str()));
    std::string refCh  = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
    std::string shortLabelRef = opts.GetOpt<std::string>(Form("%s.shortLabel",refCh.c_str()));
    if( index < 0 ) continue;
    
    c = new TCanvas(Form("c_CTR_RTCorr_posCorr_%s",shortLabel.c_str()),Form("c_CTR_%s",shortLabel.c_str()));
    
    histo = h_CTR[ch];
    histoCorr = h_CTR_RTCorr[ch];
    histoCorr2 = h_CTR_RTCorr_posCorr[ch];
    
    if( intrinsic < 0 )
      ctrResults = drawCTRPlot(histoCorr2,"r.t. corr. + pos. corr.",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,histoCorr,"r.t. corr",histo,"raw");
    else
      ctrResults = drawCTRPlot(histoCorr2,"r.t. corr. + pos. corr.",rebin,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,histoCorr,"r.t. corr",histo,"raw");
    t_CTR_RTCorr_posCorr_effSigma[ch] = ctrResults.effSigma;
    t_CTR_RTCorr_posCorr_gausSigma[ch] = ctrResults.gausSigma;
    t_CTR_RTCorr_posCorr_gausSigmaErr[ch] = ctrResults.gausSigmaErr;
    
    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();

    c -> Print(Form("%s/c_4_tRes_RTCorr_posCorr_%s.png",plotDir.c_str(),shortLabel.c_str()));
    c -> Print(Form("%s/c_4_tRes_RTCorr_posCorr_%s.pdf",plotDir.c_str(),shortLabel.c_str()));
  }
  
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string shortLabel = opts.GetOpt<std::string>(Form("%s.shortLabel",ch.c_str()));
    float intrinsic = opts.GetOpt<float>(Form("%s.intrinsic",ch.c_str()));
    std::string refCh  = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
    std::string shortLabelRef = opts.GetOpt<std::string>(Form("%s.shortLabel",refCh.c_str()));
    if( index < 0 ) continue;
    
    c = new TCanvas(Form("c_CTR_distAmpCorr_%s",shortLabel.c_str()),Form("c_CTR_distAmpCorr_%s",shortLabel.c_str()));
    
    histo = h_CTR[ch];
    histoCorr = h_CTR_ampCorr[ch];
    histoCorr2 = h_CTR_distAmpCorr[ch];
    
    if( intrinsic < 0 )
      ctrResults = drawCTRPlot(histoCorr2,"amp. walk & dist. corr.",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,histoCorr,"amp. corr",histo,"raw");
    else
      ctrResults = drawCTRPlot(histoCorr2,"amp. walk & dist. corr.",rebin,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,histoCorr,"amp. corr",histo,"raw");
    t_CTR_distAmpCorr_effSigma[ch] = ctrResults.effSigma;
    t_CTR_distAmpCorr_gausSigma[ch] = ctrResults.gausSigma;
    t_CTR_distAmpCorr_gausSigmaErr[ch] = ctrResults.gausSigmaErr;
    
    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();
    
    c -> Print(Form("%s/c_4_tRes_distAmpCorr_%s.png",plotDir.c_str(),shortLabel.c_str()));
    c -> Print(Form("%s/c_4_tRes_distAmpCorr_%s.pdf",plotDir.c_str(),shortLabel.c_str()));
  }
  
  for(auto ch: channels)
  {
    int index = opts.GetOpt<int>(Form("%s.index", ch.c_str()));
    std::string shortLabel = opts.GetOpt<std::string>(Form("%s.shortLabel",ch.c_str()));
    float intrinsic = opts.GetOpt<float>(Form("%s.intrinsic",ch.c_str()));
    std::string refCh  = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
    std::string shortLabelRef = opts.GetOpt<std::string>(Form("%s.shortLabel",refCh.c_str()));
    if( index < 0 ) continue;
    
    c = new TCanvas(Form("c_CTR_distRTCorr_%s",shortLabel.c_str()),Form("c_CTR_distRTCorr_%s",shortLabel.c_str()));
    
    histo = h_CTR[ch];
    histoCorr = h_CTR_ampCorr[ch];
    histoCorr2 = h_CTR_distRTCorr[ch];
    
    if( intrinsic < 0 )
      ctrResults = drawCTRPlot(histoCorr2,"amp. walk & r.t. corr.",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,histoCorr,"amp. corr",histo,"raw");
    else
      ctrResults = drawCTRPlot(histoCorr2,"amp. walk & r.t. corr.",rebin,false,true,intrinsic,Form(";#Deltat [ns]"),latexLabel12,histoCorr,"amp. corr",histo,"raw");
    t_CTR_distRTCorr_effSigma[ch] = ctrResults.effSigma;
    t_CTR_distRTCorr_gausSigma[ch] = ctrResults.gausSigma;
    t_CTR_distRTCorr_gausSigmaErr[ch] = ctrResults.gausSigmaErr;
    
    latexLabels[ch] -> Draw("same");
    
    gPad -> Update();
    
    c -> Print(Form("%s/c_4_tRes_distRTCorr_%s.png",plotDir.c_str(),shortLabel.c_str()));
    c -> Print(Form("%s/c_4_tRes_distRTCorr_%s.pdf",plotDir.c_str(),shortLabel.c_str()));
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
    
    TH2F* h2_CTR_ampCorr_chess = new TH2F(Form("h2_CTR_ampCorr_chess_%s",ch.c_str()),"",nBinsX,xtalXMin,xtalXMax,nBinsY,xtalYMin,xtalYMax);
    
    for(int binX = 0; binX < nBinsX; ++binX)
      for(int binY = 0; binY < nBinsY; ++binY)
      {
        std::string label(Form("%s_%d-%d",ch.c_str(),binX,binY));
        histo =  map_h_CTR_ampCorr[label];
        if( histo->GetEntries() < 100 ) continue;

        CTRResult result;        
        if( intrinsic < 0 )
          result = drawCTRPlot(histo,"amp. walk corr.",rebin,false,false,-1.,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
        else
          result = drawCTRPlot(histo,"amp. walk corr.",rebin,false,false,intrinsic,Form(";#Deltat [ns]"),latexLabel12,NULL,"",NULL,"");
        
        h2_CTR_ampCorr_chess -> SetBinContent( binX+1,binY+1,result.effSigma );
      }
    
    c = new TCanvas();
    DrawHistogram2D(opts,ch,h2_CTR_ampCorr_chess,";hodoscope x[mm];hodoscope y [mm];time resolution [ps]",xtalXMin,xtalXMax,xtalYMin,xtalYMax,20.,70.,false,false,"colz,text",latexLabels[ch]);
    PrintCanvas(c,opts,ch,plotDir,"4_tRes_chess_ampCorr");
    delete c;    
  }
  
  outTree -> Fill();

  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
