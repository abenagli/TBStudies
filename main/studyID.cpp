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
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "Math/Vector3D.h"



int main(int argc, char** argv)
{
  setTDRStyle();

  if( argc < 2 )
  {
    std::cout << ">>> studyID::usage:   " << argv[0] << " configFile.cfg" << std::endl;
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
  
  std::vector<float> ptRanges = opts.GetOpt<std::vector<float> >("Options.ptRanges");
  float ptRangesMin;
  float ptRangesMax;
  
  std::vector<float> etaRanges = opts.GetOpt<std::vector<float> >("Options.etaRanges");
  float etaRangesMin;
  float etaRangesMax;
  
  std::vector<float> EthrVals = opts.GetOpt<std::vector<float> >("Options.EthrVals");
  
  int nEtaBins = opts.GetOpt<int>("Options.nEtaBins");
  float etaMin = opts.GetOpt<float>("Options.etaMin");
  float etaMax = opts.GetOpt<float>("Options.etaMax");
  
  int nPhiBins = opts.GetOpt<int>("Options.nPhiBins");
  float phiMin = opts.GetOpt<float>("Options.phiMin");
  float phiMax = opts.GetOpt<float>("Options.phiMax");
  
  int nPtBins = opts.GetOpt<int>("Options.nPtBins");
  float ptMin = opts.GetOpt<float>("Options.ptMin");
  float ptMax = opts.GetOpt<float>("Options.ptMax");
  
  std::string plotDir = opts.GetOpt<std::string>("Input.plotDir");
  system(Form("mkdir -p %s",plotDir.c_str()));
  
  
  
  //--- get tree
  TFile* inFile = TFile::Open(inFileName.c_str(),"READ");
  TTree* tree = (TTree*)( inFile->Get("FTLDumpHits/DumpHits") );
  tree -> SetBranchStatus("*",0);
  
  std::vector<float>* tracks_pt  = new std::vector<float>;
  std::vector<float>* tracks_eta = new std::vector<float>;
  std::vector<float>* tracks_phi = new std::vector<float>;
  std::vector<float>* tracks_eta_atBTL = new std::vector<float>;
  std::vector<float>* tracks_phi_atBTL = new std::vector<float>;
  std::vector<int>* tracks_isHighPurity = new std::vector<int>;
  std::vector<float>* tracks_mcMatch_genPt = new std::vector<float>;
  std::vector<float>* tracks_mcMatch_DR = new std::vector<float>;
  tree -> SetBranchStatus("track_pt",           1); tree -> SetBranchAddress("track_pt",           &tracks_pt);
  tree -> SetBranchStatus("track_eta",          1); tree -> SetBranchAddress("track_eta",          &tracks_eta);
  tree -> SetBranchStatus("track_phi",          1); tree -> SetBranchAddress("track_phi",          &tracks_phi);
  tree -> SetBranchStatus("track_eta_atBTL",    1); tree -> SetBranchAddress("track_eta_atBTL",    &tracks_eta_atBTL);
  tree -> SetBranchStatus("track_phi_atBTL",    1); tree -> SetBranchAddress("track_phi_atBTL",    &tracks_phi_atBTL);
  tree -> SetBranchStatus("track_isHighPurity", 1); tree -> SetBranchAddress("track_isHighPurity", &tracks_isHighPurity);
  tree -> SetBranchStatus("track_mcMatch_genPt",1); tree -> SetBranchAddress("track_mcMatch_genPt",&tracks_mcMatch_genPt);
  tree -> SetBranchStatus("track_mcMatch_DR",   1); tree -> SetBranchAddress("track_mcMatch_DR",   &tracks_mcMatch_DR);
  
  std::vector<int>* recHits_det = new std::vector<int>;
  std::vector<float>* recHits_time = new std::vector<float>;
  std::vector<float>* recHits_energy = new std::vector<float>;
  std::vector<int>* recHits_module = new std::vector<int>;
  std::vector<int>* recHits_modType = new std::vector<int>;
  tree -> SetBranchStatus("recHits_det",1);       tree -> SetBranchAddress("recHits_det",       &recHits_det);
  tree -> SetBranchStatus("recHits_time",1);      tree -> SetBranchAddress("recHits_time",      &recHits_time);
  tree -> SetBranchStatus("recHits_energy",1);    tree -> SetBranchAddress("recHits_energy",    &recHits_energy);
  tree -> SetBranchStatus("recHits_module",1);    tree -> SetBranchAddress("recHits_module",    &recHits_module);
  tree -> SetBranchStatus("recHits_modType",1);   tree -> SetBranchAddress("recHits_modType",   &recHits_modType);
  
  std::vector<std::vector<int> >*   matchedRecHits_det = new std::vector<std::vector<int> >;
  std::vector<std::vector<float> >* matchedRecHits_time = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* matchedRecHits_energy = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* matchedRecHits_energyCorr = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* matchedRecHits_track_DR = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* matchedRecHits_track_dist = new std::vector<std::vector<float> >;
  std::vector<std::vector<int> >*   matchedRecHits_modType = new std::vector<std::vector<int> >;
  tree -> SetBranchStatus("matchedRecHits_det",1);        tree -> SetBranchAddress("matchedRecHits_det",        &matchedRecHits_det);
  tree -> SetBranchStatus("matchedRecHits_time",1);       tree -> SetBranchAddress("matchedRecHits_time",       &matchedRecHits_time);
  tree -> SetBranchStatus("matchedRecHits_energy",1);     tree -> SetBranchAddress("matchedRecHits_energy",     &matchedRecHits_energy);
  tree -> SetBranchStatus("matchedRecHits_energyCorr",1); tree -> SetBranchAddress("matchedRecHits_energyCorr", &matchedRecHits_energyCorr);
  tree -> SetBranchStatus("matchedRecHits_track_DR",1);   tree -> SetBranchAddress("matchedRecHits_track_DR",   &matchedRecHits_track_DR);
  tree -> SetBranchStatus("matchedRecHits_track_dist",1); tree -> SetBranchAddress("matchedRecHits_track_dist", &matchedRecHits_track_dist);
  tree -> SetBranchStatus("matchedRecHits_modType",1);    tree -> SetBranchAddress("matchedRecHits_modType",    &matchedRecHits_modType);
  
  
  
  //--- open output file
  TFile* outFile = TFile::Open(Form("%s/IDPlots_%s.root",plotDir.c_str(),label.c_str()),"RECREATE");
  outFile -> cd();  
  
  TH1F* h1_tracks_n = new TH1F("h1_tracks_n","",250,-0.5,249.5);
  TH1F* h1_tracks_pt = new TH1F("h1_tracks_pt","",nPtBins,ptMin,ptMax);
  TH1F* h1_tracks_eta = new TH1F("h1_tracks_eta","",nEtaBins,etaMin,etaMax);
  TH1F* h1_tracks_phi = new TH1F("h1_tracks_phi","",nPhiBins,phiMin,phiMax);
  TH1F* h1_tracks_eta_atBTL = new TH1F("h1_tracks_eta_atBTL","",nEtaBins,etaMin,etaMax);
  TH1F* h1_tracks_phi_atBTL = new TH1F("h1_tracks_phi_atBTL","",nPhiBins,phiMin,phiMax);
  
  TH1F* h1_recHits_energy = new TH1F("h1_recHit_energy","",1000,0.,100.);
  TH1F* h1_recHits_time = new TH1F("h1_recHit_time","",500,-10.,40.);
  
  TH1F* h1_matchedRecHit_energySumCorr = new TH1F("h1_matchedRecHit_energySumCorr","",1000,0.,100.);
  TH1F* h1_matchedRecHit_energySum = new TH1F("h1_matchedRecHit_energySum","",1000,0.,100.);
  TH1F* h1_matchedRecHit_energy = new TH1F("h1_matchedRecHit_energy","",1000,0.,100.);
  TH1F* h1_matchedRecHit_time = new TH1F("h1_matchedRecHit_time","",500,-10.,40.);
  TH1F* h1_matchedRecHit_track_DR = new TH1F("h1_matchedRecHit_track_DR","",10000,0.,10.);
  TH1F* h1_matchedRecHit_track_dist = new TH1F("h1_matchedRecHit_track_dist","",10000,0.,1000.);

  TProfile* p1_matchedRecHit_time_vs_eta = new TProfile("p1_matchedRecHit_time_vs_eta","",nEtaBins,etaMin,etaMax);
  
  TH1F* h1_nEntries_ptRanges = new TH1F("h1_nEntries_ptRanges","",ptRanges.size()-1,ptRanges.data());
  TH1F* h1_nEntries_etaRanges = new TH1F("h1_nEntries_etaRanges","",etaRanges.size()-1,etaRanges.data());
  /*
  for(auto Ethr : EthrVals)
  {
    p1_matchedRecHit_eff_vs_pt[Ethr] = new TEfficiency(Form("p1_matchedRecHit_eff_vs_pt__Ethr%.1fMeV",Ethr),"",nPtBins,ptMin,ptMax);
    
    for(unsigned int ii = 0; ii < ptRanges.size()-1; ++ii)
    {
      ptRangesMin = ptRanges.at(ii);
      ptRangesMax = ptRanges.at(ii+1);
      
      h1_matchedRecHit_totEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TH1F(Form("h1_matchedRecHit_totEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",1000.,0.,500.);
      
      p1_matchedRecHit_eff_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TEfficiency(Form("p1_matchedRecHit_eff_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
      p1_matchedRecHit_n_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_n_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
      p1_matchedRecHit_totEnergy_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_totEnergy_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
      p1_matchedRecHit_avgEnergy_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_avgEnergy_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
      p1_matchedRecHit_maxEnergy_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_maxEnergy_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
      p1_matchedRecHit_maxOverTotEnergy_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_maxOverTotEnergy_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
      p1_matchedRecHit_timeRes_vs_eta_totEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_timeRes_vs_eta__totEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
      p1_matchedRecHit_timeRes_vs_eta_maxEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_timeRes_vs_eta__maxEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
      p1_matchedRecHit_timeRes_vs_eta_sumEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_timeRes_vs_eta__sumEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
      p1_matchedRecHit_timeRes_stat_vs_eta_totEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_timeRes_stat_vs_eta__totEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
      p1_matchedRecHit_timeRes_DCR_vs_eta_totEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_timeRes_DCR_vs_eta__totEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
      
      p1_matchedRecHit_eff_vs_phi[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TEfficiency(Form("p1_matchedRecHit_eff_vs_phi__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nPhiBins,phiMin,phiMax);
      p1_matchedRecHit_eff_vs_phiFold[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TEfficiency(Form("p1_matchedRecHit_eff_vs_phiFold__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",4*nPhiBins,phiMin,phiMax);
      p1_matchedRecHit_n_vs_phi[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_n_vs_phi__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nPhiBins,phiMin,phiMax);
      p1_matchedRecHit_totEnergy_vs_phi[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_totEnergy_vs_phi__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nPhiBins,phiMin,phiMax);
      p1_matchedRecHit_avgEnergy_vs_phi[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_avgEnergy_vs_phi__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nPhiBins,phiMin,phiMax);
    }
    
    ptRangesMin = 0.8;
    ptRangesMax = ptRanges.at(ptRanges.size()-1);

    h1_matchedRecHit_totEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TH1F(Form("h1_matchedRecHit_totEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",1000.,0.,500.);
      
    p1_matchedRecHit_eff_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TEfficiency(Form("p1_matchedRecHit_eff_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
    p1_matchedRecHit_n_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_n_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
    p1_matchedRecHit_totEnergy_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_totEnergy_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
    p1_matchedRecHit_avgEnergy_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_avgEnergy_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
    p1_matchedRecHit_maxEnergy_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_maxEnergy_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
    p1_matchedRecHit_maxOverTotEnergy_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_maxOverTotEnergy_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
    p1_matchedRecHit_timeRes_vs_eta_totEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_timeRes_vs_eta__totEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
    p1_matchedRecHit_timeRes_vs_eta_maxEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_timeRes_vs_eta__maxEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
    p1_matchedRecHit_timeRes_vs_eta_sumEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_timeRes_vs_eta__sumEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
    p1_matchedRecHit_timeRes_stat_vs_eta_totEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_timeRes_stat_vs_eta__totEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
    p1_matchedRecHit_timeRes_DCR_vs_eta_totEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_timeRes_DCR_vs_eta__totEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nEtaBins,etaMin,etaMax);
    
    p1_matchedRecHit_eff_vs_phi[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TEfficiency(Form("p1_matchedRecHit_eff_vs_phi__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nPhiBins,phiMin,phiMax);
    p1_matchedRecHit_eff_vs_phiFold[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TEfficiency(Form("p1_matchedRecHit_eff_vs_phiFold__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",4*nPhiBins,phiMin,phiMax);
    p1_matchedRecHit_n_vs_phi[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_n_vs_phi__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nPhiBins,phiMin,phiMax);
    p1_matchedRecHit_totEnergy_vs_phi[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_totEnergy_vs_phi__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nPhiBins,phiMin,phiMax);
    p1_matchedRecHit_avgEnergy_vs_phi[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedRecHit_avgEnergy_vs_phi__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",nPhiBins,phiMin,phiMax);
  }
  */
  
  
  //--- loop over events
  int nEntries = tree->GetEntries();
  // nEntries = 10000;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%100 == 0 ) std::cout << ">>> reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    
    tree -> GetEntry(entry);
    
    
    
    //------------------------------
    //--- fill matched recHits plots
    if( debugMode ) std::cout << ">>> fill matched recHits plots" << std::endl;
    int nGoodTracks = 0;
    for(unsigned int trackIt = 0; trackIt < tracks_pt->size(); ++trackIt)
    {
      float pt = tracks_pt->at(trackIt);
      float eta = tracks_eta->at(trackIt);
      float eta_atBTL = tracks_eta_atBTL->at(trackIt);
      float phi_atBTL = tracks_phi_atBTL->at(trackIt);
      int isHighPurity = 1;//tracks_isHighPurity->at(trackIt);
      float genPt = tracks_mcMatch_genPt->at(trackIt);
      float DR = tracks_mcMatch_DR->at(trackIt);
      
      if( eta_atBTL < -100. ) continue;
      if( isHighPurity != 1 ) continue;
      if( DR > 0.01 ) continue;
      if( fabs(pt/genPt-1.) > 0.05 ) continue;
      
      ++nGoodTracks;
      h1_tracks_pt -> Fill( tracks_pt->at(trackIt) );
      h1_tracks_eta -> Fill( fabs(tracks_eta->at(trackIt)) );
      h1_tracks_phi -> Fill( tracks_phi->at(trackIt) );
      h1_tracks_eta_atBTL -> Fill( fabs(tracks_eta_atBTL->at(trackIt)) );
      h1_tracks_phi_atBTL -> Fill( tracks_phi_atBTL->at(trackIt) );
      
      int bin = h1_nEntries_ptRanges -> Fill( pt );
      if( bin < 1 || bin > int(ptRanges.size()) ) continue;

      float energySum = 0.;
      float energySumCorr = 0.;
      float energyMax = -999.;
      for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(trackIt)).size(); ++recHitIt)
      {
        if( (matchedRecHits_det->at(trackIt)).at(recHitIt) != 1 ) continue;
        float recHitE = (matchedRecHits_energy->at(trackIt)).at(recHitIt);
        
        h1_matchedRecHit_energy -> Fill( recHitE );
        h1_matchedRecHit_time -> Fill( (matchedRecHits_time->at(trackIt)).at(recHitIt) );
        
        h1_matchedRecHit_track_DR -> Fill( (matchedRecHits_track_DR->at(trackIt)).at(recHitIt) );
        h1_matchedRecHit_track_dist -> Fill( (matchedRecHits_track_dist->at(trackIt)).at(recHitIt) );
        
        energySum += recHitE;
        if( (matchedRecHits_modType->at(trackIt)).at(recHitIt) == 1 )
        {
          energySumCorr += (matchedRecHits_energyCorr->at(trackIt)).at(recHitIt);
          if( recHitE > energyMax )
          {
            energyMax = recHitE;
          }
        }
        if( (matchedRecHits_modType->at(trackIt)).at(recHitIt) == 2 )
        {
          energySumCorr += (matchedRecHits_energyCorr->at(trackIt)).at(recHitIt) * 3.75/3.00;
          if( recHitE > energyMax )
          {
            energyMax = recHitE;
          }
        }
        if( (matchedRecHits_modType->at(trackIt)).at(recHitIt) == 3 )
        {
          energySumCorr += (matchedRecHits_energyCorr->at(trackIt)).at(recHitIt) * 3.75/2.40;
          if( recHitE > energyMax )
          {
            energyMax = recHitE;
          }
        }
        
        p1_matchedRecHit_time_vs_eta -> Fill( fabs(eta),(matchedRecHits_time->at(trackIt)).at(recHitIt) );
      }
      h1_matchedRecHit_energySum -> Fill( energySum );
      h1_matchedRecHit_energySumCorr -> Fill( energySumCorr );
        
    } // end loop over tracks
    
    h1_tracks_n -> Fill( nGoodTracks );
    
  } // end loop over events
  
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;

}
