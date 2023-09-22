#include "interface/AnalysisUtils.h"
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
#include "TEfficiency.h"



int main(int argc, char** argv)
{
  setTDRStyle();

  if( argc < 2 )
  {
    std::cout << ">>> geometryPlots::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  

  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  

  
  //--- get parameters
  std::string inFileName = opts.GetOpt<std::string>("Input.fileName");
  std::string label = opts.GetOpt<std::string>("Input.label");
  
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  system(Form("mkdir -p %s",plotDir.c_str()));
  
  
  
  //--- get tree
  TFile* inFile = TFile::Open(inFileName.c_str(),"READ");
  TTree* tree = (TTree*)( inFile->Get("FTLDumpHits/hits_tree") );

  std::vector<int>*   gen_idx = new std::vector<int>;
  std::vector<float>* gen_pt  = new std::vector<float>;
  std::vector<float>* gen_eta = new std::vector<float>;
  std::vector<float>* gen_phi = new std::vector<float>;
  
  std::vector<std::vector<float> >* simHits_time = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* simHits_energy = new std::vector<std::vector<float> >;
  std::vector<float>* simHits_energySum = new std::vector<float>;
  std::vector<float>* simHits_n = new std::vector<float>;

  std::vector<std::vector<float> >* tracks_pt  = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* tracks_eta = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* tracks_phi = new std::vector<std::vector<float> >;
  
  std::vector<std::vector<float> >* recHits_time = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* recHits_energy = new std::vector<std::vector<float> >;
  std::vector<float>* recHits_energySum = new std::vector<float>;
  std::vector<float>* recHits_n = new std::vector<float>;
  
  tree -> SetBranchAddress("idx", &gen_idx);
  tree -> SetBranchAddress("pt",  &gen_pt);
  tree -> SetBranchAddress("eta", &gen_eta);
  tree -> SetBranchAddress("phi", &gen_phi);

  tree -> SetBranchAddress("simHits_time",      &simHits_time);
  tree -> SetBranchAddress("simHits_energy",    &simHits_energy);
  tree -> SetBranchAddress("simHits_energySum", &simHits_energySum);
  tree -> SetBranchAddress("simHits_n",         &simHits_n);

  tree -> SetBranchAddress("tracks_pt",  &tracks_pt);
  tree -> SetBranchAddress("tracks_eta", &tracks_eta);
  tree -> SetBranchAddress("tracks_phi", &tracks_phi);
    
  tree -> SetBranchAddress("recHits_time",      &recHits_time);
  tree -> SetBranchAddress("recHits_energy",    &recHits_energy);
  tree -> SetBranchAddress("recHits_energySum", &recHits_energySum);
  tree -> SetBranchAddress("recHits_n",         &recHits_n);
  
  
  
  //--- open output file
  TFile* outFile = TFile::Open(Form("%s/geometryPlots_%s.root",plotDir.c_str(),label.c_str()),"RECREATE");
  outFile -> cd();  

  TH1F* h1_track_DR = new TH1F("h1_track_DR","",1000,0.,1.);
  
  TH1F* h1_simHits_n = new TH1F("h1_simHits_n","",20,-0.5,19.5);
  TH1F* h1_simHits_energy = new TH1F("h1_simHits_energy","",1000,0.,100.);
  TH1F* h1_simHits_energySum = new TH1F("h1_simHits_energySum","",1000,0.,100.);
  TProfile* p1_simHits_n_vs_eta = new TProfile("p1_simHits_n_vs_eta","",50,0.,1.5);
  TProfile* p1_simHits_n_vs_phi = new TProfile("p1_simHits_n_vs_phi","",720,-3.15,3.15);
  TProfile* p1_simHits_energy_vs_eta = new TProfile("p1_simHits_energy_vs_eta","",50,0.,1.5);
  TProfile* p1_simHits_energy_vs_phi = new TProfile("p1_simHits_energy_vs_phi","",720,-3.15,3.15);
  TProfile* p1_simHits_energySum_vs_eta = new TProfile("p1_simHits_energySum_vs_eta","",50,0.,1.5);
  TProfile* p1_simHits_energySum_vs_phi = new TProfile("p1_simHits_energySum_vs_phi","",720,-3.15,3.15);
  TProfile* p1_simHits_time_vs_eta = new TProfile("p1_simHits_time_vs_eta","",50,0.,1.5);
  TProfile* p1_simHits_time_vs_phi = new TProfile("p1_simHits_time_vs_phi","",720,-3.15,3.15);
  
  TH1F* h1_recHits_n = new TH1F("h1_recHits_n","",20,-0.5,19.5);
  TH1F* h1_recHits_energy = new TH1F("h1_recHits_energy","",1000,0.,100.);
  TH1F* h1_recHits_energySum = new TH1F("h1_recHits_energySum","",1000,0.,100.);
  TProfile* p1_recHits_n_vs_eta = new TProfile("p1_recHits_n_vs_eta","",50,0.,1.5);
  TProfile* p1_recHits_n_vs_phi = new TProfile("p1_recHits_n_vs_phi","",720,-3.15,3.15);
  TProfile* p1_recHits_energy_vs_eta = new TProfile("p1_recHits_energy_vs_eta","",50,0.,1.5);
  TProfile* p1_recHits_energy_vs_phi = new TProfile("p1_recHits_energy_vs_phi","",720,-3.15,3.15);
  TProfile* p1_recHits_energySum_vs_eta = new TProfile("p1_recHits_energySum_vs_eta","",50,0.,1.5);
  TProfile* p1_recHits_energySum_vs_phi = new TProfile("p1_recHits_energySum_vs_phi","",720,-3.15,3.15);
  TProfile* p1_recHits_time_vs_eta = new TProfile("p1_recHits_time_vs_eta","",50,0.,1.5);
  TProfile* p1_recHits_time_vs_phi = new TProfile("p1_recHits_time_vs_phi","",720,-3.15,3.15);
  
  TEfficiency* p1_eff_vs_eta = new TEfficiency("p1_eff_vs_eta","",100,0.,1.5);
  TEfficiency* p1_eff_vs_phi = new TEfficiency("p1_eff_vs_phi","",720,-3.15,3.15);
  
  
  
  //--- loop over events
  int nEntries = tree->GetEntries();
  for(int entry = 0; entry < nEntries; ++entry)
  {
    tree -> GetEntry(entry);
    
    int idx = 0;

    
    if( tracks_pt->at(idx).size() == 0 ) continue;
    
    float DR = DeltaR(gen_eta->at(idx),gen_phi->at(idx),tracks_eta->at(idx).at(0),tracks_phi->at(idx).at(0));
    h1_track_DR -> Fill( DR );
    if( DR > 0.03 ) continue;

    
    for(unsigned int jj = 0; jj < simHits_energy->at(idx).size(); ++jj)
    {
      // std::cout << "energy: " << simHits_energy->at(idx).at(jj) << std::endl;
      h1_simHits_energy -> Fill( 1000.*simHits_energy->at(idx).at(jj) );
      p1_simHits_energy_vs_eta -> Fill( gen_eta->at(idx),1000.*simHits_energy->at(idx).at(jj) );
      p1_simHits_energy_vs_phi -> Fill( gen_phi->at(idx),1000.*simHits_energy->at(idx).at(jj) );
      if( 1000.*simHits_energy->at(idx).at(jj) > 1. )
      {
        p1_simHits_time_vs_eta -> Fill( gen_eta->at(idx),1000.*simHits_time->at(idx).at(jj) );
        p1_simHits_time_vs_phi -> Fill( gen_phi->at(idx),1000.*simHits_time->at(idx).at(jj) );
      }
    }
    h1_simHits_n -> Fill( simHits_n->at(idx) );
    p1_simHits_n_vs_eta -> Fill( gen_eta->at(idx),simHits_n->at(idx) );
    p1_simHits_n_vs_phi -> Fill( gen_phi->at(idx),simHits_n->at(idx) );
    h1_simHits_energySum -> Fill( 1000.*simHits_energySum->at(idx) );
    p1_simHits_energySum_vs_eta -> Fill( gen_eta->at(idx),1000.*simHits_energySum->at(idx) );
    p1_simHits_energySum_vs_phi -> Fill( gen_phi->at(idx),1000.*simHits_energySum->at(idx) );
    
    
    for(unsigned int jj = 0; jj < recHits_energy->at(idx).size(); ++jj)
    {
      // std::cout << "energy: " << recHits_energy->at(idx).at(jj) << std::endl;
      h1_recHits_energy -> Fill( recHits_energy->at(idx).at(jj) );
      p1_recHits_energy_vs_eta -> Fill( gen_eta->at(idx),recHits_energy->at(idx).at(jj) );
      p1_recHits_energy_vs_phi -> Fill( gen_phi->at(idx),recHits_energy->at(idx).at(jj) );
      if( recHits_energy->at(idx).at(jj) > 1. )
      {
        p1_recHits_time_vs_eta -> Fill( gen_eta->at(idx),recHits_time->at(idx).at(jj) );
        p1_recHits_time_vs_phi -> Fill( gen_phi->at(idx),recHits_time->at(idx).at(jj) );
      }
    }
    h1_recHits_n -> Fill( recHits_n->at(idx) );
    p1_recHits_n_vs_eta -> Fill( gen_eta->at(idx),recHits_n->at(idx) );
    p1_recHits_n_vs_phi -> Fill( gen_phi->at(idx),recHits_n->at(idx) );
    h1_recHits_energySum -> Fill( recHits_energySum->at(idx) );
    p1_recHits_energySum_vs_eta -> Fill( gen_eta->at(idx),recHits_energySum->at(idx) );
    p1_recHits_energySum_vs_phi -> Fill( gen_phi->at(idx),recHits_energySum->at(idx) );
    
    
    p1_eff_vs_eta -> Fill( recHits_energySum->at(idx)>1.,gen_eta->at(idx) );
    if( gen_eta->at(idx)< 1.4 )
      p1_eff_vs_phi -> Fill( recHits_energySum->at(idx)>1.,gen_phi->at(idx) );
  }
  
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;

}

