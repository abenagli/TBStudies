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
#include "TGraph.h"



int main(int argc, char** argv)
{
  setTDRStyle();

  
  if( argc < 2 )
  {
    std::cout << ">>> testTLU::usage:   " << argv[0] << " configFile.cfg" << std::endl;
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
  
  
  
  //--- open input files
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  int runMin = opts.GetOpt<int>("Input.runMin");
  int runMax = opts.GetOpt<int>("Input.runMax");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");

  long long int refTime_h4 = opts.GetOpt<long long int>("Input.refTimeH4");
  long long int refTime_tlu = opts.GetOpt<long long int>("Input.refTimeTLU");
  
  TChain* h4 = new TChain("h4","h4");
  for(int run = runMin; run <= runMax; ++run)
  {
    std::string fileName = Form("%s/%d.root",inputDir.c_str(),run);
    std::cout << ">>> Adding flle " << fileName << std::endl;
    h4 -> Add(fileName.c_str());
  }
  
  std::string inputFileTLU = opts.GetOpt<std::string>("Input.inputFileTLU");
  TFile* inFileTLU = TFile::Open(inputFileTLU.c_str(),"READ");

  TTree* tlu = (TTree*)( inFileTLU->Get("tracks") );
  

  
  //--- define tree branch addresses
  unsigned int run;
  unsigned int spill;
  long long int h4_timestamp;
  float* time    = new float[1000];
  float* amp_max = new float[1000];
  int AMPBAR01;
  int AMPMAT09;
  int AMPMAT10;
  int AMPMAT11;
  int AMPMAT12;
  int AMPMAT13;
  int AMPMAT14;
  int AMPMAT15;
  int AMPMAT16;
  h4 -> SetBranchAddress("timestamp",&h4_timestamp);
  h4 -> SetBranchAddress("run",       &run);
  h4 -> SetBranchAddress("spill",     &spill);
  h4 -> SetBranchAddress("time",       time);
  h4 -> SetBranchAddress("amp_max",    amp_max);
  h4 -> SetBranchAddress("AMPBAR1",   &AMPBAR01);
  h4 -> SetBranchAddress("AMPMAT9",   &AMPMAT09);
  h4 -> SetBranchAddress("AMPMAT10",  &AMPMAT10);
  h4 -> SetBranchAddress("AMPMAT11",  &AMPMAT11);
  h4 -> SetBranchAddress("AMPMAT12",  &AMPMAT12);
  h4 -> SetBranchAddress("AMPMAT13",  &AMPMAT13);
  h4 -> SetBranchAddress("AMPMAT14",  &AMPMAT14);
  h4 -> SetBranchAddress("AMPMAT15",  &AMPMAT15);
  h4 -> SetBranchAddress("AMPMAT16",  &AMPMAT16);
      
  long long int tlu_timestamp;
  int euEvt;
  std::vector<int>* iden = new std::vector<int>;
  std::vector<double>* xPos = new std::vector<double>;
  std::vector<double>* yPos = new std::vector<double>;
  std::vector<double>* dxdz = new std::vector<double>;
  std::vector<double>* dydz = new std::vector<double>;
  tlu -> SetBranchAddress("timeStamp",&tlu_timestamp);
  tlu -> SetBranchAddress("euEvt",&euEvt);
  tlu -> SetBranchAddress("iden", &iden);
  tlu -> SetBranchAddress("xPos", &xPos);
  tlu -> SetBranchAddress("yPos", &yPos);
  tlu -> SetBranchAddress("dxdz", &dxdz);
  tlu -> SetBranchAddress("dydz", &dydz);
  
  int nEntries_h4  = h4->GetEntries();
  int nEntries_tlu = tlu->GetEntries();
  std::cout << ">>> Events read (h4): "  << nEntries_h4 << std::endl;
  std::cout << ">>> Events read (tlu): " << nEntries_tlu << std::endl;  
  
  
  
  //------------------
  // define histograms
  
  TFile* outFile = TFile::Open(Form("%s/testTLU.root",plotDir.c_str()),"RECREATE");
  outFile -> cd();
  
  TH1F* histo = new TH1F("histo","",10000000,-100000000,100000000);

  TProfile* prof = new TProfile("prof","",9299,1232800,1242099);
  TGraph* g = new TGraph();

  TProfile2D* p2_amp_AMPBAR01 = new TProfile2D("p2_amp_AMPBAR01","",100,-20.,20.,100,-20.,20.);
    
  TProfile2D* p2_amp_AMPMAT09 = new TProfile2D("p2_amp_AMPMAT09","",100,-20.,20.,100,-20.,20.);
  TProfile2D* p2_amp_AMPMAT10 = new TProfile2D("p2_amp_AMPMAT10","",100,-20.,20.,100,-20.,20.);
  TProfile2D* p2_amp_AMPMAT11 = new TProfile2D("p2_amp_AMPMAT11","",100,-20.,20.,100,-20.,20.);
  TProfile2D* p2_amp_AMPMAT12 = new TProfile2D("p2_amp_AMPMAT12","",100,-20.,20.,100,-20.,20.);
  
  TProfile2D* p2_amp_AMPMAT13 = new TProfile2D("p2_amp_AMPMAT13","",100,-20.,20.,100,-20.,20.);
  TProfile2D* p2_amp_AMPMAT14 = new TProfile2D("p2_amp_AMPMAT14","",100,-20.,20.,100,-20.,20.);
  TProfile2D* p2_amp_AMPMAT15 = new TProfile2D("p2_amp_AMPMAT15","",100,-20.,20.,100,-20.,20.);
  TProfile2D* p2_amp_AMPMAT16 = new TProfile2D("p2_amp_AMPMAT16","",100,-20.,20.,100,-20.,20.);

  
  //-----------------------
  // 1st loop over events
  if( maxEntries > 0 ) nEntries_h4 = maxEntries;
  int entry_tlu = 0;
  for(int entry = 4; entry < nEntries_h4; ++entry)
  {
    if( entry%1000 == 0 ) std::cout << ">>> 1st loop: reading entry " << entry << " / " << nEntries_h4 << "\r" << std::flush;
    h4 -> GetEntry(entry);
    tlu -> GetEntry(entry_tlu);
    
    //long long int diff = ((tlu_timestamp-refTime_tlu)/384.02523) - (h4_timestamp-refTime_h4);
    long long int diff = ((tlu_timestamp-refTime_tlu)/384.02527) - (h4_timestamp-refTime_h4);
    
    if( fabs(diff) > 20000. )
    {
      std::cout << std::endl;
      std::cout << "run: " << run << "   spill: " << spill << std::endl;
      std::cout << "h4 entry: "  << entry     << "   h4 timestamp: " << h4_timestamp-refTime_h4 << std::endl;
      std::cout << "tlu entry: " << entry_tlu << "   tlu event: " << euEvt << "   tlu timestamp: " << std::fixed << std::setprecision(0) << (tlu_timestamp-refTime_tlu)/384.025 << std::endl;
      std::cout << "diff: " << diff << std::endl;
    }
        
    histo -> Fill( diff );
    
    if( fabs(diff) > 200000. )
    {
      entry -= 1;
      ++entry_tlu;
      if( entry_tlu >= nEntries_tlu )
        break;
      
      continue;
    }
    
    g -> SetPoint(g->GetN(),h4_timestamp-refTime_h4,(tlu_timestamp-refTime_tlu)/384.025);

    float x = -99.;
    float y = -99.;
    for(int ii = 0; ii < int(iden->size()/1); ++ii)
    {
      x = xPos->at(1+1*ii-1) + dxdz->at(1+1*ii-1)*830.;
      y = yPos->at(1+1*ii-1) + dydz->at(1+1*ii-1)*830.;
      
      float amp;
      amp = amp_max[AMPBAR01]; p2_amp_AMPBAR01 -> Fill(x,y,amp);
      amp = amp_max[AMPMAT09]; p2_amp_AMPMAT09 -> Fill(x,y,amp);
      amp = amp_max[AMPMAT10]; p2_amp_AMPMAT10 -> Fill(x,y,amp);
      amp = amp_max[AMPMAT11]; p2_amp_AMPMAT11 -> Fill(x,y,amp);
      amp = amp_max[AMPMAT12]; p2_amp_AMPMAT12 -> Fill(x,y,amp);
      amp = amp_max[AMPMAT13]; p2_amp_AMPMAT13 -> Fill(x,y,amp);
      amp = amp_max[AMPMAT14]; p2_amp_AMPMAT14 -> Fill(x,y,amp);
      amp = amp_max[AMPMAT15]; p2_amp_AMPMAT15 -> Fill(x,y,amp);
      amp = amp_max[AMPMAT16]; p2_amp_AMPMAT16 -> Fill(x,y,amp);
    }

    prof -> Fill( run*100+spill,diff );
    
    ++entry_tlu;
  }
  std::cout << "\n>>> end 1st loop" << std::endl;
  
  
  g -> Write("g");
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
