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
    std::cout << ">>> studyHits::usage:   " << argv[0] << " configFile.cfg" << std::endl;
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
  
  int nCrystalsPerMatrix = opts.GetOpt<int>("Input.nCrystalsPerMatrix");
  int nMatricesPerMod = opts.GetOpt<int>("Input.nMatricesPerMod");
  int nModsPerRU = opts.GetOpt<int>("Input.nModsPerRU");
  int nModsPerType = opts.GetOpt<int>("Input.nModsPerType");
  int nRUs = opts.GetOpt<int>("Input.nRUs");
  int nTrays = opts.GetOpt<int>("Input.nTrays");
  
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
  
  float sigma_clock = opts.GetOpt<float>("Options.sigma_clock");
  float sigma_digi  = opts.GetOpt<float>("Options.sigma_digi");
  float sigma_ele   = opts.GetOpt<float>("Options.sigma_ele");
  float LY  = opts.GetOpt<float>("Options.LY");
  float LCE = opts.GetOpt<float>("Options.LCE");
  float PDE = opts.GetOpt<float>("Options.PDE");
  float DCR = opts.GetOpt<float>("Options.DCR");
  float sigma_red = 1.;
  std::size_t found = label.find("bar");
  if( found != std::string::npos) sigma_red = sqrt(2.);

  TF1* f_fluence_vs_eta = new TF1("f_fluence_vs_eta","pol2",0.,1.5);
  f_fluence_vs_eta -> SetParameters(1.7,0.0474763,0.109945);
  
  std::string plotDir = opts.GetOpt<std::string>("Input.plotDir");
  plotDir += std::string(Form("%s_PDE%.2f_DCR%.0fGHz",label.c_str(),PDE,DCR));
  system(Form("mkdir -p %s",plotDir.c_str()));
  
  
  
  //--- get tree
  TFile* inFile = TFile::Open(inFileName.c_str(),"READ");
  //TTree* tree = (TTree*)( inFile->Get("FTLDumpHits/hits_tree") );
  TTree* tree = (TTree*)( inFile->Get("FTLDumpHits/DumpHits") );
  // TTree* tree = (TTree*)( inFile->Get("DumpHits") );
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

  std::vector<int>* simHits_det = new std::vector<int>;
  std::vector<float>* simHits_time = new std::vector<float>;
  std::vector<float>* simHits_energy = new std::vector<float>;
  std::vector<int>* simHits_ieta = new std::vector<int>;
  std::vector<int>* simHits_iphi = new std::vector<int>;  
  std::vector<int>* simHits_module = new std::vector<int>;
  std::vector<int>* simHits_modType = new std::vector<int>;
  std::vector<float>* simHits_entry_local_x = new std::vector<float>;
  std::vector<float>* simHits_entry_local_y = new std::vector<float>;
  std::vector<float>* simHits_entry_local_z = new std::vector<float>;
  std::vector<float>* simHits_exit_local_x = new std::vector<float>;
  std::vector<float>* simHits_exit_local_y = new std::vector<float>;
  std::vector<float>* simHits_exit_local_z = new std::vector<float>;
  tree -> SetBranchStatus("simHits_det",1);       tree -> SetBranchAddress("simHits_det",       &simHits_det);
  tree -> SetBranchStatus("simHits_time",1);      tree -> SetBranchAddress("simHits_time",      &simHits_time);
  tree -> SetBranchStatus("simHits_energy",1);    tree -> SetBranchAddress("simHits_energy",    &simHits_energy);
  tree -> SetBranchStatus("simHits_ieta",1);      tree -> SetBranchAddress("simHits_ieta",      &simHits_ieta);
  tree -> SetBranchStatus("simHits_iphi",1);      tree -> SetBranchAddress("simHits_iphi",      &simHits_iphi);
  tree -> SetBranchStatus("simHits_module",1);    tree -> SetBranchAddress("simHits_module",    &simHits_module);
  tree -> SetBranchStatus("simHits_modType",1);   tree -> SetBranchAddress("simHits_modType",   &simHits_modType);
  // tree -> SetBranchStatus("simHits_entry_local_x",1); tree -> SetBranchAddress("simHits_entry_local_x",&simHits_entry_local_x);
  // tree -> SetBranchStatus("simHits_entry_local_y",1); tree -> SetBranchAddress("simHits_entry_local_y",&simHits_entry_local_y);
  // tree -> SetBranchStatus("simHits_entry_local_z",1); tree -> SetBranchAddress("simHits_entry_local_z",&simHits_entry_local_z);
  // tree -> SetBranchStatus("simHits_exit_local_x",1); tree -> SetBranchAddress("simHits_exit_local_x",&simHits_exit_local_x);
  // tree -> SetBranchStatus("simHits_exit_local_y",1); tree -> SetBranchAddress("simHits_exit_local_y",&simHits_exit_local_y);
  // tree -> SetBranchStatus("simHits_exit_local_z",1); tree -> SetBranchAddress("simHits_exit_local_z",&simHits_exit_local_z);
  
  std::vector<std::vector<int> >*   matchedSimHits_det = new std::vector<std::vector<int> >;
  std::vector<std::vector<float> >* matchedSimHits_time = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* matchedSimHits_energy = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* matchedSimHits_entry_local_x = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* matchedSimHits_entry_local_y = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* matchedSimHits_entry_local_z = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* matchedSimHits_ieta = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* matchedSimHits_iphi = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* matchedSimHits_track_RDphi = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* matchedSimHits_track_Dz = new std::vector<std::vector<float> >;
  tree -> SetBranchStatus("matchedSimHits_det",          1); tree -> SetBranchAddress("matchedSimHits_det",          &matchedSimHits_det);
  tree -> SetBranchStatus("matchedSimHits_time",         1); tree -> SetBranchAddress("matchedSimHits_time",         &matchedSimHits_time);
  tree -> SetBranchStatus("matchedSimHits_energy",       1); tree -> SetBranchAddress("matchedSimHits_energy",       &matchedSimHits_energy);
  tree -> SetBranchStatus("matchedSimHits_entry_local_x",1); tree -> SetBranchAddress("matchedSimHits_entry_local_x",&matchedSimHits_entry_local_x);
  tree -> SetBranchStatus("matchedSimHits_entry_local_y",1); tree -> SetBranchAddress("matchedSimHits_entry_local_y",&matchedSimHits_entry_local_y);
  tree -> SetBranchStatus("matchedSimHits_entry_local_z",1); tree -> SetBranchAddress("matchedSimHits_entry_local_z",&matchedSimHits_entry_local_z);
  tree -> SetBranchStatus("matchedSimHits_ieta",         1); tree -> SetBranchAddress("matchedSimHits_ieta",         &matchedSimHits_ieta);
  tree -> SetBranchStatus("matchedSimHits_iphi",         1); tree -> SetBranchAddress("matchedSimHits_iphi",         &matchedSimHits_iphi);
  tree -> SetBranchStatus("matchedSimHits_track_RDphi",  1); tree -> SetBranchAddress("matchedSimHits_track_RDphi",  &matchedSimHits_track_RDphi);
  tree -> SetBranchStatus("matchedSimHits_track_Dz",     1); tree -> SetBranchAddress("matchedSimHits_track_Dz",     &matchedSimHits_track_Dz);
  
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
  TFile* outFile = TFile::Open(Form("%s/hitsPlots_%s.root",plotDir.c_str(),label.c_str()),"RECREATE");
  outFile -> cd();  
  
  TH1F* h1_tracks_n = new TH1F("h1_tracks_n","",250,-0.5,249.5);
  TH1F* h1_tracks_pt = new TH1F("h1_tracks_pt","",500,0.,20.);
  TH1F* h1_tracks_eta = new TH1F("h1_tracks_eta","",nEtaBins,etaMin,etaMax);
  TH1F* h1_tracks_phi = new TH1F("h1_tracks_phi","",nPhiBins,phiMin,phiMax);
  TH1F* h1_tracks_eta_atBTL = new TH1F("h1_tracks_eta_atBTL","",nEtaBins,etaMin,etaMax);
  TH1F* h1_tracks_phi_atBTL = new TH1F("h1_tracks_phi_atBTL","",nPhiBins,phiMin,phiMax);
  
  TH1F* h1_simHits_energy = new TH1F("h1_simHit_energy","",1000,0.,100.);
  TH1F* h1_simHits_time = new TH1F("h1_simHit_time","",500,-10.,40.);

  std::map<std::pair<int,int>,TH1F*> h1_simHit_time_perChannel;
  std::map<std::pair<int,int>,TH1F*> h1_simHit_timeRMS_perChannel;
  
  TH1F* h1_matchedSimHit_track_RDphi = new TH1F("h1_matchedSimHit_track_RDphi","",10000,-1.,1.);
  TH1F* h1_matchedSimHit_track_Dz = new TH1F("h1_matchedSimHit_track_Dz","",10000,-5.,5.);
  std::map<std::pair<float,float>,std::map<std::pair<float,float>,TH1F*> > h1_matchedSimHit_track_RDphi_pt_eta;
  std::map<std::pair<float,float>,std::map<std::pair<float,float>,TH1F*> > h1_matchedSimHit_track_Dz_pt_eta;
  
  std::map<float,std::map<std::pair<float,float>,TProfile*> > p1_matchedSimHit_totEnergy_vs_local_x_pt;
  std::map<float,std::map<std::pair<float,float>,TProfile*> > p1_matchedSimHit_totEnergy_vs_local_y_pt;
  std::map<float,std::map<std::pair<float,float>,TProfile*> > p1_matchedSimHit_maxEnergy_vs_local_x_pt;
  std::map<float,std::map<std::pair<float,float>,TProfile*> > p1_matchedSimHit_maxEnergy_vs_local_y_pt;
  std::map<float,std::map<std::pair<float,float>,TProfile*> > p1_matchedSimHit_maxOverTotEnergy_vs_local_x_pt;
  std::map<float,std::map<std::pair<float,float>,TProfile*> > p1_matchedSimHit_maxOverTotEnergy_vs_local_y_pt;
  
  std::map<float,std::map<std::pair<float,float>,TProfile*> > p1_matchedSimHit_timeRes_vs_local_x_totEnergy;
  std::map<float,std::map<std::pair<float,float>,TProfile*> > p1_matchedSimHit_timeRes_vs_local_x_maxEnergy;
  std::map<float,std::map<std::pair<float,float>,TProfile*> > p1_matchedSimHit_timeRes_vs_local_x_sumEnergy;
  std::map<float,std::map<std::pair<float,float>,TProfile*> > p1_matchedSimHit_timeRes_vs_local_y_totEnergy;
  std::map<float,std::map<std::pair<float,float>,TProfile*> > p1_matchedSimHit_timeRes_vs_local_y_maxEnergy;
  std::map<float,std::map<std::pair<float,float>,TProfile*> > p1_matchedSimHit_timeRes_vs_local_y_sumEnergy;
  
  for(unsigned int ii = 0; ii < ptRanges.size()-1; ++ii)
  {
    ptRangesMin = ptRanges.at(ii);
    ptRangesMax = ptRanges.at(ii+1);

    for(auto Ethr : EthrVals)
    {
      if( label == "tile" )
      {
        p1_matchedSimHit_totEnergy_vs_local_x_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_totEnergy_vs_local_x__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-0.69,0.69);
        p1_matchedSimHit_totEnergy_vs_local_y_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_totEnergy_vs_local_y__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-0.69,0.69);
        p1_matchedSimHit_maxEnergy_vs_local_x_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_maxEnergy_vs_local_x__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-0.69,0.69);
        p1_matchedSimHit_maxEnergy_vs_local_y_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_maxEnergy_vs_local_y__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-0.69,0.69);
        p1_matchedSimHit_maxOverTotEnergy_vs_local_x_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_maxOverTotEnergy_vs_local_x__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-0.69,0.69);
        p1_matchedSimHit_maxOverTotEnergy_vs_local_y_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_maxOverTotEnergy_vs_local_y__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-0.69,0.69);
      }
      if( label == "barphi" )
      {
        p1_matchedSimHit_totEnergy_vs_local_x_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_totEnergy_vs_local_x__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-3.,3.);
        p1_matchedSimHit_totEnergy_vs_local_y_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_totEnergy_vs_local_y__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-0.18,0.18);
        p1_matchedSimHit_maxEnergy_vs_local_x_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_maxEnergy_vs_local_x__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-3.,3.);
        p1_matchedSimHit_maxEnergy_vs_local_y_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_maxEnergy_vs_local_y__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-0.18,0.18);
        p1_matchedSimHit_maxOverTotEnergy_vs_local_x_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_maxOverTotEnergy_vs_local_x__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-3.,3.);
        p1_matchedSimHit_maxOverTotEnergy_vs_local_y_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_maxOverTotEnergy_vs_local_y__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-0.18,0.18);
      }
      if( label == "barzflat" )
      {
        p1_matchedSimHit_totEnergy_vs_local_x_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_totEnergy_vs_local_x__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-0.18,0.18);
        p1_matchedSimHit_totEnergy_vs_local_y_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_totEnergy_vs_local_y__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-3.,3.);
        p1_matchedSimHit_maxEnergy_vs_local_x_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_maxEnergy_vs_local_x__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-0.18,0.18);
        p1_matchedSimHit_maxEnergy_vs_local_y_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_maxEnergy_vs_local_y__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-3.,3.);
        p1_matchedSimHit_maxOverTotEnergy_vs_local_x_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_maxOverTotEnergy_vs_local_x__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-0.18,0.18);
        p1_matchedSimHit_maxOverTotEnergy_vs_local_y_pt[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] = new TProfile(Form("p1_matchedSimHit_maxOverTotEnergy_vs_local_y__pt%04.1f-%04.1f__Ethr%.1fMeV",ptRangesMin,ptRangesMax,Ethr),"",200,-3.,3.);
      }
    }
    
    for(unsigned int jj = 0; jj < etaRanges.size()-1; ++jj)
    {
      etaRangesMin = etaRanges.at(jj);
      etaRangesMax = etaRanges.at(jj+1);    
      
      h1_matchedSimHit_track_RDphi_pt_eta[std::make_pair(ptRangesMin,ptRangesMax)][std::make_pair(etaRangesMin,etaRangesMax)] = new TH1F(Form("h1_matchedSimHit_track_RDphi__pt%04.1f-%04.1f__eta%02.1f-%02.1f",ptRangesMin,ptRangesMax,etaRangesMin,etaRangesMax),"",1000,-1.,1.);
      h1_matchedSimHit_track_Dz_pt_eta[std::make_pair(ptRangesMin,ptRangesMax)][std::make_pair(etaRangesMin,etaRangesMax)] = new TH1F(Form("h1_matchedSimHit_track_Dz__pt%04.1f-%04.1f__eta%02.1f-%02.1f",ptRangesMin,ptRangesMax,etaRangesMin,etaRangesMax),"",1000,-5.,5.);
    }
  }
  
  std::map<float,std::map<int,int> > simHits_n_vs_RU_energyCut;
  std::map<float,std::map<int,int> > simHits_n_vs_ieta_energyCut;
  std::map<float,std::map<int,int> > recHits_n_vs_RU_energyCut;
  
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

  std::map<float,std::map<std::pair<float,float>,TH1F*> > h1_matchedRecHit_totEnergy;

  double xAxis_Ethr[101];
  xAxis_Ethr[0] = pow(10.,-2.005);
  for(int jj = 1; jj <= 100; ++jj)
    xAxis_Ethr[jj] = pow(10.,-2.005+3./100.*jj);
  TEfficiency* p1_matchedRecHit_eff_vs_Ethr = new TEfficiency("p1_matchedRecHit_eff_vs_Ethr","",100,xAxis_Ethr);
  TProfile* p1_matchedRecHit_timeRes_vs_Ethr_totEnergy = new TProfile("p1_matchedRecHit_timeRes_vs_Ethr_totEnergy","",100,xAxis_Ethr);
  TProfile* p1_matchedRecHit_timeRes_vs_Ethr_maxEnergy = new TProfile("p1_matchedRecHit_timeRes_vs_Ethr_maxEnergy","",100,xAxis_Ethr);
  TProfile* p1_matchedRecHit_timeRes_vs_Ethr_sumEnergy = new TProfile("p1_matchedRecHit_timeRes_vs_Ethr_sumEnergy","",100,xAxis_Ethr);
  
  std::map<float,std::map<std::pair<float,float>,TEfficiency*> > p1_matchedRecHit_eff_vs_eta;
  std::map<float,std::map<std::pair<float,float>,TProfile*> >    p1_matchedRecHit_n_vs_eta;
  std::map<float,std::map<std::pair<float,float>,TProfile*> >    p1_matchedRecHit_totEnergy_vs_eta;
  std::map<float,std::map<std::pair<float,float>,TProfile*> >    p1_matchedRecHit_avgEnergy_vs_eta;
  std::map<float,std::map<std::pair<float,float>,TProfile*> >    p1_matchedRecHit_maxEnergy_vs_eta;
  std::map<float,std::map<std::pair<float,float>,TProfile*> >    p1_matchedRecHit_maxOverTotEnergy_vs_eta;
  std::map<float,std::map<std::pair<float,float>,TProfile*> >    p1_matchedRecHit_timeRes_vs_eta_totEnergy;
  std::map<float,std::map<std::pair<float,float>,TProfile*> >    p1_matchedRecHit_timeRes_vs_eta_maxEnergy;
  std::map<float,std::map<std::pair<float,float>,TProfile*> >    p1_matchedRecHit_timeRes_vs_eta_sumEnergy;
  std::map<float,std::map<std::pair<float,float>,TProfile*> >    p1_matchedRecHit_timeRes_stat_vs_eta_totEnergy;
  std::map<float,std::map<std::pair<float,float>,TProfile*> >    p1_matchedRecHit_timeRes_DCR_vs_eta_totEnergy;
  
  std::map<float,std::map<std::pair<float,float>,TEfficiency*> > p1_matchedRecHit_eff_vs_phi;
  std::map<float,std::map<std::pair<float,float>,TEfficiency*> > p1_matchedRecHit_eff_vs_phiFold;
  std::map<float,std::map<std::pair<float,float>,TProfile*> >    p1_matchedRecHit_n_vs_phi;
  std::map<float,std::map<std::pair<float,float>,TProfile*> >    p1_matchedRecHit_totEnergy_vs_phi;
  std::map<float,std::map<std::pair<float,float>,TProfile*> >    p1_matchedRecHit_avgEnergy_vs_phi;

  std::map<float,TEfficiency*> p1_matchedRecHit_eff_vs_pt;
  
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
  
  
  
  //--- loop over events
  int nEntries = tree->GetEntries();
  // nEntries = 10000;
  for(int entry = 0; entry < nEntries; ++entry)
  {
    if( entry%100 == 0 ) std::cout << ">>> reading entry " << entry << " / " << nEntries << "\r" << std::flush;
    
    tree -> GetEntry(entry);
    

    
    
    
    
    //--------------------------
    //--- fill all simHits plots
    if( debugMode ) std::cout << ">>> fill all simHits plots" << std::endl;
    std::map<std::pair<int,int>,int> nTracks_perSimHit;
    std::map<std::pair<int,int>,float> energy_perSimHit;
    std::map<std::pair<int,int>,int> RU_perSimHit;
    bool write = false;
    for(unsigned int simHitIt = 0; simHitIt < simHits_energy->size(); ++simHitIt)
    {
      if( simHits_det->at(simHitIt) != 1 ) continue;
      int ieta = simHits_ieta->at(simHitIt);
      int iphi = simHits_iphi->at(simHitIt);
      int RU = ((simHits_module->at(simHitIt)-1) + (simHits_modType->at(simHitIt)-1)*nModsPerType) / nModsPerRU;
      float energy = simHits_energy->at(simHitIt);
      float time = simHits_time->at(simHitIt);
      
      // std::pair<int,int> map_key(ieta,iphi);
      // if( h1_simHit_time_perChannel[map_key] == NULL )
      // {
      //   h1_simHit_time_perChannel[map_key] = new TH1F(Form("h1_simHit_time_perChannel_%d-%d",ieta,iphi),"",5000,-5.,45.);
      //   h1_simHit_time_perChannel[map_key] -> SetDirectory(0);
      // }
      // h1_simHit_time_perChannel[map_key] -> Fill( time,energy );
      
      energy_perSimHit[std::make_pair(ieta,iphi)] += energy;
      RU_perSimHit[std::make_pair(ieta,iphi)] = RU;
    }
    
    for(int jj = 0; jj < 200; ++jj)
    {
      float cut = pow(10.,-2.+3./200.*jj);
      
      for(auto mapIt : energy_perSimHit)
      {
        int RU = RU_perSimHit[mapIt.first];
        float energy = energy_perSimHit[mapIt.first];
        if( energy > cut )
        {
          ++simHits_n_vs_RU_energyCut[cut][RU];
          ++simHits_n_vs_ieta_energyCut[cut][mapIt.first.first];
        }
      }
    }
    
    for(unsigned int simHitIt = 0; simHitIt < simHits_energy->size(); ++simHitIt)
    {
      h1_simHits_energy -> Fill( simHits_energy->at(simHitIt) );
      h1_simHits_time -> Fill( simHits_time->at(simHitIt) );
    }

    
    //------------------------------
    //--- fill matched simHits plots
    if( debugMode ) std::cout << ">>> fill matched simHits plots" << std::endl;
    for(unsigned int trackIt = 0; trackIt < tracks_pt->size(); ++trackIt)
    {
      float pt = tracks_pt->at(trackIt);
      float eta = tracks_eta->at(trackIt);
      float eta_atBTL = tracks_eta_atBTL->at(trackIt);
      int isHighPurity = 1;//tracks_isHighPurity->at(trackIt);
      float genPt = tracks_mcMatch_genPt->at(trackIt);
      float DR = tracks_mcMatch_DR->at(trackIt);
      
      if( eta_atBTL < -100. ) continue;
      if( isHighPurity != 1 ) continue;
      if( DR > 0.01 ) continue;
      if( fabs(pt/genPt-1.) > 0.05 ) continue;
      
      int ptBin = h1_nEntries_ptRanges -> Fill( pt );
      if( ptBin < 1 || ptBin > int(ptRanges.size()) ) continue;
      
      int etaBin = h1_nEntries_etaRanges -> Fill( fabs(eta) );
      if( etaBin < 1 || etaBin > int(etaRanges.size()) ) continue;
      
      for(unsigned int simHitIt = 0; simHitIt < (matchedSimHits_energy->at(trackIt)).size(); ++simHitIt)
      {
        h1_matchedSimHit_track_RDphi -> Fill( (matchedSimHits_track_RDphi->at(trackIt)).at(simHitIt) );
        h1_matchedSimHit_track_Dz -> Fill( (matchedSimHits_track_Dz->at(trackIt)).at(simHitIt) );
        
        h1_matchedSimHit_track_RDphi_pt_eta[std::make_pair(ptRanges.at(ptBin-1),ptRanges.at(ptBin))][std::make_pair(etaRanges.at(etaBin-1),etaRanges.at(etaBin))] -> Fill( (matchedSimHits_track_RDphi->at(trackIt)).at(simHitIt) );
        h1_matchedSimHit_track_Dz_pt_eta[std::make_pair(ptRanges.at(ptBin-1),ptRanges.at(ptBin))][std::make_pair(etaRanges.at(etaBin-1),etaRanges.at(etaBin))] -> Fill( (matchedSimHits_track_Dz->at(trackIt)).at(simHitIt) );
      }
      
      std::map<std::pair<int,int>,float> energy_perRecHit;
      int entry_ieta = 999;
      int entry_iphi = 999;
      float entry_x = 999.;
      float entry_y = 999.;
      float entry_z_min = 999.;
      for(unsigned int simHitIt = 0; simHitIt < (matchedSimHits_energy->at(trackIt)).size(); ++simHitIt)
      {
        if( (matchedSimHits_det->at(trackIt)).at(simHitIt) != 1 ) continue;
        int ieta = (matchedSimHits_ieta->at(trackIt)).at(simHitIt);
        int iphi = (matchedSimHits_iphi->at(trackIt)).at(simHitIt);
        float energy = (matchedSimHits_energy->at(trackIt)).at(simHitIt);
        
        energy_perRecHit[std::make_pair(ieta,iphi)] += energy;
        
        if( (matchedSimHits_entry_local_z->at(trackIt)).at(simHitIt) < entry_z_min )
        {
          entry_z_min = (matchedSimHits_entry_local_z->at(trackIt)).at(simHitIt);
          entry_x = (matchedSimHits_entry_local_x->at(trackIt)).at(simHitIt);
          entry_y = (matchedSimHits_entry_local_y->at(trackIt)).at(simHitIt);
          entry_ieta = ieta;
          entry_iphi = iphi;
        }
      }
      
      for(auto Ethr : EthrVals)
      {
        float energy_tot = 0.;
        float energy_max = -999;
        for(auto mapIt : energy_perRecHit)
        {
          if( mapIt.second < Ethr ) continue;
          if( mapIt.second > energy_max ) energy_max = mapIt.second;
          energy_tot += mapIt.second;
        }
        
        // if( energy_tot <= 0. ) continue;
        
        // p1_matchedSimHit_totEnergy_vs_local_x_pt[Ethr][std::make_pair(ptRanges.at(ptBin-1),ptRanges.at(ptBin))] -> Fill( entry_x,energy_tot );
        // p1_matchedSimHit_totEnergy_vs_local_y_pt[Ethr][std::make_pair(ptRanges.at(ptBin-1),ptRanges.at(ptBin))] -> Fill( entry_y,energy_tot );
        // p1_matchedSimHit_maxEnergy_vs_local_x_pt[Ethr][std::make_pair(ptRanges.at(ptBin-1),ptRanges.at(ptBin))] -> Fill( entry_x,energy_max );
        // p1_matchedSimHit_maxEnergy_vs_local_y_pt[Ethr][std::make_pair(ptRanges.at(ptBin-1),ptRanges.at(ptBin))] -> Fill( entry_y,energy_max );
        // p1_matchedSimHit_maxOverTotEnergy_vs_local_x_pt[Ethr][std::make_pair(ptRanges.at(ptBin-1),ptRanges.at(ptBin))] -> Fill( entry_x,energy_max/energy_tot );
        // p1_matchedSimHit_maxOverTotEnergy_vs_local_y_pt[Ethr][std::make_pair(ptRanges.at(ptBin-1),ptRanges.at(ptBin))] -> Fill( entry_y,energy_max/energy_tot );
      }
    }
    
    
    //--------------------------
    //--- fill all recHits plots
    if( debugMode ) std::cout << ">>> fill all recHits plots" << std::endl;
    for(int jj = 0; jj < 200; ++jj)
    {
      float cut = pow(10.,-2.+3./200.*jj);
      
      for(unsigned int recHitIt = 0; recHitIt < recHits_energy->size(); ++recHitIt)
      {
        if( recHits_det->at(recHitIt) != 1 ) continue;
        
        int RU = ((recHits_module->at(recHitIt)-1) + (recHits_modType->at(recHitIt)-1)*nModsPerType) / nModsPerRU;
        if( recHits_energy->at(recHitIt) > cut ) ++recHits_n_vs_RU_energyCut[cut][RU];
      }
    }
    
    for(unsigned int recHitIt = 0; recHitIt < recHits_energy->size(); ++recHitIt)
    {
      if( recHits_det->at(recHitIt) != 1 ) continue;
      
      h1_recHits_energy -> Fill( recHits_energy->at(recHitIt) );
      h1_recHits_time -> Fill( recHits_time->at(recHitIt) );
    }
    
    
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

      float phi_atBTL_fold = -999.;
      int iTray = int( (fabs(phi_atBTL) - (0.0415-0.006)) / (2.*((2.*3.14159-4.*0.0415)/36.)) );
      phi_atBTL_fold = fabs(phi_atBTL) - 2.*((2.*3.14159-4.*0.0415)/36.) * iTray - 0.006;
      
      ++nGoodTracks;
      h1_tracks_pt -> Fill( tracks_pt->at(trackIt) );
      h1_tracks_eta -> Fill( fabs(tracks_eta->at(trackIt)) );
      h1_tracks_phi -> Fill( tracks_phi->at(trackIt) );
      h1_tracks_eta_atBTL -> Fill( fabs(tracks_eta_atBTL->at(trackIt)) );
      h1_tracks_phi_atBTL -> Fill( tracks_phi_atBTL->at(trackIt) );
      
      int bin = h1_nEntries_ptRanges -> Fill( pt );
      if( bin < 1 || bin > int(ptRanges.size()) ) continue;

      float DCRCorr = 1.;
      float fluenceCorr = 1.;
      
      float energyMax = -999999;
      float energySum = 0.;
      float energySumCorr = 0.;
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
            DCRCorr = 3.15*3.75/9.;
            fluenceCorr = f_fluence_vs_eta->Eval(fabs(tracks_eta_atBTL->at(trackIt))) / f_fluence_vs_eta->Eval(1.45);
          }
        }
        if( (matchedRecHits_modType->at(trackIt)).at(recHitIt) == 2 )
        {
          energySumCorr += (matchedRecHits_energyCorr->at(trackIt)).at(recHitIt) * 3.75/3.00;
          if( recHitE > energyMax )
          {
            energyMax = recHitE;
            DCRCorr = 3.15*3.00/9.;
            fluenceCorr = f_fluence_vs_eta->Eval(fabs(tracks_eta_atBTL->at(trackIt))) / f_fluence_vs_eta->Eval(1.45);
          }
        }
        if( (matchedRecHits_modType->at(trackIt)).at(recHitIt) == 3 )
        {
          energySumCorr += (matchedRecHits_energyCorr->at(trackIt)).at(recHitIt) * 3.75/2.40;
          if( recHitE > energyMax )
          {
            energyMax = recHitE;
            DCRCorr = 3.15*2.40/9.;
            fluenceCorr = f_fluence_vs_eta->Eval(fabs(tracks_eta_atBTL->at(trackIt))) / f_fluence_vs_eta->Eval(1.45);
          }
        }
        
        p1_matchedRecHit_time_vs_eta -> Fill( fabs(eta),(matchedRecHits_time->at(trackIt)).at(recHitIt) );
      }
      h1_matchedRecHit_energySum -> Fill( energySum );
      h1_matchedRecHit_energySumCorr -> Fill( energySumCorr );
      

      // plot vs Ethr
      for(int jj = 0; jj < 100; ++jj)
      {
        float cut = pow(10.,-2.+3./100.*jj);

        int matchedRecHit_n = 0;
        float matchedRecHit_totEnergy = 0.;
        float matchedRecHit_maxEnergy = -999.;
        std::vector<float> Npe_sumEnergy;
        
        for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(trackIt)).size(); ++recHitIt)
        {
          if( (matchedRecHits_energy->at(trackIt)).at(recHitIt) > cut )
          {
            ++matchedRecHit_n;
            matchedRecHit_totEnergy += (matchedRecHits_energy->at(trackIt)).at(recHitIt);
            if( (matchedRecHits_energy->at(trackIt)).at(recHitIt) > matchedRecHit_maxEnergy )
              matchedRecHit_maxEnergy = (matchedRecHits_energy->at(trackIt)).at(recHitIt);
            Npe_sumEnergy.push_back( LY*(matchedRecHits_energy->at(trackIt)).at(recHitIt)*LCE*PDE );
          }
        }
        
        if( matchedRecHit_totEnergy > 2.75 ) p1_matchedRecHit_eff_vs_Ethr -> Fill( 1.,cut );
        else                                 p1_matchedRecHit_eff_vs_Ethr -> Fill( 0.,cut );
        
        if( matchedRecHit_totEnergy < 2.75 ) continue;
        
        float Npe_totEnergy = LY*matchedRecHit_totEnergy*LCE*PDE;
        float sigmaT_stat_totEnergy = 35.4*sqrt(12480/Npe_totEnergy);
        float sigmaT_DCR_totEnergy = 32.*sqrt(DCR*DCRCorr*fluenceCorr/20.)*(9000./Npe_totEnergy);
        float sigmaT_totEnergy = sqrt( pow(sigma_clock,2.) + pow(sigma_digi/sigma_red,2.) + pow(sigma_ele/sigma_red,2.) + pow(sigmaT_stat_totEnergy/sigma_red,2.) + pow(sigmaT_DCR_totEnergy/sigma_red,2.) );
        
        float Npe_maxEnergy = LY*matchedRecHit_maxEnergy*LCE*PDE;
        float sigmaT_stat_maxEnergy = 35.4*sqrt(12480/Npe_maxEnergy);
        float sigmaT_DCR_maxEnergy = 32.*sqrt(DCR*DCRCorr*fluenceCorr/20.)*(9000./Npe_maxEnergy);
        float sigmaT_maxEnergy = sqrt( pow(sigma_clock,2.) + pow(sigma_digi/sigma_red,2.) + pow(sigma_ele/sigma_red,2.) + pow(sigmaT_stat_maxEnergy/sigma_red,2.) + pow(sigmaT_DCR_maxEnergy/sigma_red,2.) );
        
        float sigmaT_sumEnergy = 0.;
        for(auto vecIt : Npe_sumEnergy)
        {
          float Npe_temp = vecIt;
          float sigmaT_stat_temp = 35.4*sqrt(12480/Npe_temp);
          float sigmaT_DCR_temp = 32.*sqrt(DCR*DCRCorr*fluenceCorr/20.)*(9000./Npe_temp);
          float sigmaT_temp = sqrt( pow(sigma_clock,2.) + pow(sigma_digi/sigma_red,2.) + pow(sigma_ele/sigma_red,2.) + pow(sigmaT_stat_temp/sigma_red,2.) + pow(sigmaT_DCR_temp/sigma_red,2.) );
          sigmaT_sumEnergy += 1. / pow(sigmaT_temp,2.);
        }
        sigmaT_sumEnergy = sqrt( 1./sigmaT_sumEnergy );
        
        p1_matchedRecHit_timeRes_vs_Ethr_totEnergy -> Fill( cut,sigmaT_totEnergy );
        p1_matchedRecHit_timeRes_vs_Ethr_maxEnergy -> Fill( cut,sigmaT_maxEnergy );
        p1_matchedRecHit_timeRes_vs_Ethr_sumEnergy -> Fill( cut,sigmaT_sumEnergy );
      }
      
      
      // loop over Ethr
      for(auto Ethr : EthrVals)
      {
        int matchedRecHit_n = 0;
        float matchedRecHit_totEnergy = 0.;
        float matchedRecHit_maxEnergy = -999.;
        std::vector<float> Npe_sumEnergy;
        
        for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(trackIt)).size(); ++recHitIt)
        {
          if( (matchedRecHits_energy->at(trackIt)).at(recHitIt) > Ethr )
          {
            ++matchedRecHit_n;
            matchedRecHit_totEnergy += (matchedRecHits_energy->at(trackIt)).at(recHitIt);
            if( (matchedRecHits_energy->at(trackIt)).at(recHitIt) > matchedRecHit_maxEnergy )
              matchedRecHit_maxEnergy = (matchedRecHits_energy->at(trackIt)).at(recHitIt);
            Npe_sumEnergy.push_back( LY*(matchedRecHits_energy->at(trackIt)).at(recHitIt)*LCE*PDE );
          }
        }
        
        float Npe_totEnergy = LY*matchedRecHit_totEnergy*LCE*PDE;
        float sigmaT_stat_totEnergy = 35.4*sqrt(12480/Npe_totEnergy);
        float sigmaT_DCR_totEnergy = 32.*sqrt(DCR*DCRCorr*fluenceCorr/20.)*(9000./Npe_totEnergy);
        float sigmaT_totEnergy = sqrt( pow(sigma_clock,2.) + pow(sigma_digi/sigma_red,2.) + pow(sigma_ele/sigma_red,2.) + pow(sigmaT_stat_totEnergy/sigma_red,2.) + pow(sigmaT_DCR_totEnergy/sigma_red,2.) );
        
        float Npe_maxEnergy = LY*matchedRecHit_maxEnergy*LCE*PDE;
        float sigmaT_stat_maxEnergy = 35.4*sqrt(12480/Npe_maxEnergy);
        float sigmaT_DCR_maxEnergy = 32.*sqrt(DCR*DCRCorr*fluenceCorr/20.)*(9000./Npe_maxEnergy);
        float sigmaT_maxEnergy = sqrt( pow(sigma_clock,2.) + pow(sigma_digi/sigma_red,2.) + pow(sigma_ele/sigma_red,2.) + pow(sigmaT_stat_maxEnergy/sigma_red,2.) + pow(sigmaT_DCR_maxEnergy/sigma_red,2.) );
        
        float sigmaT_sumEnergy = 0.;
        for(auto vecIt : Npe_sumEnergy)
        {
          float Npe_temp = vecIt;
          float sigmaT_stat_temp = 35.4*sqrt(12480/Npe_temp);
          float sigmaT_DCR_temp = 32.*sqrt(DCR*DCRCorr*fluenceCorr/20.)*(9000./Npe_temp);
          float sigmaT_temp = sqrt( pow(sigma_clock,2.) + pow(sigma_digi/sigma_red,2.) + pow(sigma_ele/sigma_red,2.) + pow(sigmaT_stat_temp/sigma_red,2.) + pow(sigmaT_DCR_temp/sigma_red,2.) );
          sigmaT_sumEnergy += 1. / pow(sigmaT_temp,2.);
        }
        sigmaT_sumEnergy = sqrt( 1./sigmaT_sumEnergy );
        
        
        // all pt
        p1_matchedRecHit_n_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( eta_atBTL,matchedRecHit_n );
        if( matchedRecHit_totEnergy > 2.75 ) p1_matchedRecHit_eff_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( 1.,fabs(eta_atBTL) );
        else                                 p1_matchedRecHit_eff_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( 0.,fabs(eta_atBTL) );
        
        p1_matchedRecHit_n_vs_phi[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( phi_atBTL,matchedRecHit_n );
        if( matchedRecHit_totEnergy > 2.75 ) p1_matchedRecHit_eff_vs_phi[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( 1.,phi_atBTL );
        else                                 p1_matchedRecHit_eff_vs_phi[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( 0.,phi_atBTL );
        if( matchedRecHit_totEnergy > 2.75 ) p1_matchedRecHit_eff_vs_phiFold[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( 1.,phi_atBTL_fold );
        else                                 p1_matchedRecHit_eff_vs_phiFold[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( 0.,phi_atBTL_fold );
        
        if( matchedRecHit_totEnergy > 2.75 ) p1_matchedRecHit_eff_vs_pt[Ethr] -> Fill( 1.,pt );
        else                                 p1_matchedRecHit_eff_vs_pt[Ethr] -> Fill( 0.,pt );
                
        h1_matchedRecHit_totEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( matchedRecHit_totEnergy );
        
        if( matchedRecHit_totEnergy > 2.75 )
        {
          p1_matchedRecHit_totEnergy_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( fabs(eta_atBTL),matchedRecHit_totEnergy );
          p1_matchedRecHit_avgEnergy_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( fabs(eta_atBTL),matchedRecHit_totEnergy/matchedRecHit_n );
          p1_matchedRecHit_maxEnergy_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( fabs(eta_atBTL),matchedRecHit_maxEnergy );
          p1_matchedRecHit_maxOverTotEnergy_vs_eta[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( fabs(eta_atBTL),matchedRecHit_maxEnergy/matchedRecHit_totEnergy );
          
          p1_matchedRecHit_timeRes_vs_eta_totEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( fabs(eta_atBTL),sigmaT_totEnergy );
          p1_matchedRecHit_timeRes_vs_eta_maxEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( fabs(eta_atBTL),sigmaT_maxEnergy );
          p1_matchedRecHit_timeRes_vs_eta_sumEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( fabs(eta_atBTL),sigmaT_sumEnergy );

          p1_matchedRecHit_timeRes_stat_vs_eta_totEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( fabs(eta_atBTL),sigmaT_stat_totEnergy/sigma_red );
          p1_matchedRecHit_timeRes_DCR_vs_eta_totEnergy[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( fabs(eta_atBTL),sigmaT_DCR_totEnergy/sigma_red );
          
          p1_matchedRecHit_totEnergy_vs_phi[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( phi_atBTL,matchedRecHit_totEnergy );
          p1_matchedRecHit_avgEnergy_vs_phi[Ethr][std::make_pair(ptRangesMin,ptRangesMax)] -> Fill( phi_atBTL,matchedRecHit_totEnergy/matchedRecHit_n );
        }

        
        // per pt-bin
        p1_matchedRecHit_n_vs_eta[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( eta_atBTL,matchedRecHit_n );
        if( matchedRecHit_totEnergy > 2.75 ) p1_matchedRecHit_eff_vs_eta[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( 1.,fabs(eta_atBTL) );
        else                                 p1_matchedRecHit_eff_vs_eta[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( 0.,fabs(eta_atBTL) );
        
        p1_matchedRecHit_n_vs_phi[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( phi_atBTL,matchedRecHit_n );
        if( matchedRecHit_n > 0 ) p1_matchedRecHit_eff_vs_phi[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( 1.,phi_atBTL );
        else                      p1_matchedRecHit_eff_vs_phi[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( 0.,phi_atBTL );
        if( matchedRecHit_n > 0 ) p1_matchedRecHit_eff_vs_phiFold[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( 1.,phi_atBTL_fold );
        else                      p1_matchedRecHit_eff_vs_phiFold[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( 0.,phi_atBTL_fold );
        
        h1_matchedRecHit_totEnergy[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( matchedRecHit_totEnergy );
        
        if( matchedRecHit_totEnergy > 2.75 )
        {
          p1_matchedRecHit_totEnergy_vs_eta[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( fabs(eta_atBTL),matchedRecHit_totEnergy );
          p1_matchedRecHit_avgEnergy_vs_eta[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( fabs(eta_atBTL),matchedRecHit_totEnergy/matchedRecHit_n );
          p1_matchedRecHit_maxEnergy_vs_eta[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( fabs(eta_atBTL),matchedRecHit_maxEnergy );
          p1_matchedRecHit_maxOverTotEnergy_vs_eta[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( fabs(eta_atBTL),matchedRecHit_maxEnergy/matchedRecHit_totEnergy );

          p1_matchedRecHit_timeRes_vs_eta_totEnergy[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( fabs(eta_atBTL),sigmaT_totEnergy );
          p1_matchedRecHit_timeRes_vs_eta_maxEnergy[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( fabs(eta_atBTL),sigmaT_maxEnergy );
          p1_matchedRecHit_timeRes_vs_eta_sumEnergy[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( fabs(eta_atBTL),sigmaT_sumEnergy );
          
          p1_matchedRecHit_totEnergy_vs_phi[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( phi_atBTL,matchedRecHit_totEnergy );
          p1_matchedRecHit_avgEnergy_vs_phi[Ethr][std::make_pair(ptRanges.at(bin-1),ptRanges.at(bin))] -> Fill( phi_atBTL,matchedRecHit_totEnergy/matchedRecHit_n );
        }
      } // end loop over Ethr
      
    } // end loop over tracks
    
    h1_tracks_n -> Fill( nGoodTracks );

    
    // if( entry%200 == 0 && entry > 0 )
    // {
    //   for(auto map_key : h1_simHit_time_perChannel )
    //   {
    //     int ieta = map_key.first.first;
    //     int iphi = map_key.first.second;

    //     if( h1_simHit_timeRMS_perChannel[map_key.first] == NULL )
    //     {
    //       h1_simHit_timeRMS_perChannel[map_key.first] = new TH1F(Form("h1_simHit_timeRMS_perChannel_%d-%d",ieta,iphi),"",1000,0.,20.);
    //       // h1_simHit_timeRMS_perChannel[map_key.first] -> SetDirectory(0);
    //     }
        
    //     h1_simHit_timeRMS_perChannel[map_key.first] -> Fill( map_key.second->GetRMS() );
        
    //     map_key.second -> Reset("ICESM");
    //   }
    // }

    
  } // end loop over events
  
  
  
  //-------------------------
  // tracker resolution plots
  for(unsigned int ii = 0; ii < ptRanges.size()-1; ++ii)
  {
    ptRangesMin = ptRanges.at(ii);
    ptRangesMax = ptRanges.at(ii+1);
    
    for(unsigned int jj = 0; jj < etaRanges.size()-1; ++jj)
    {
      etaRangesMin = etaRanges.at(jj);
      etaRangesMax = etaRanges.at(jj+1);    

      TH1F* histo;
      TF1* func;
      
      histo = h1_matchedSimHit_track_RDphi_pt_eta[std::make_pair(ptRangesMin,ptRangesMax)][std::make_pair(etaRangesMin,etaRangesMax)];
      histo -> Fit("gaus","QLS+","",0.-1.*histo->GetRMS(),0.+1.*histo->GetRMS());
      func = (TF1*)( histo->GetFunction("gaus") );
      func -> SetLineWidth(3);
      
      histo = h1_matchedSimHit_track_Dz_pt_eta[std::make_pair(ptRangesMin,ptRangesMax)][std::make_pair(etaRangesMin,etaRangesMax)];
      histo -> Fit("gaus","QLS+","",0.-1.*histo->GetRMS(),0.+1.*histo->GetRMS());
      func = (TF1*)( histo->GetFunction("gaus") );
      func -> SetLineWidth(3);
    }
  }
  
  
  
  outFile -> cd();
  
  for(int iRU = 0; iRU < nRUs; ++iRU)
  {
    TGraph* g_simHits_n_vs_RU_energyCut = new TGraph();
    TGraph* g_simHits_PU200_n_vs_RU_energyCut = new TGraph();
    TGraph* g_simHits_occ_vs_RU_energyCut = new TGraph();
    TGraph* g_simHits_PU200_occ_vs_RU_energyCut = new TGraph();

    TGraph* g_recHits_n_vs_RU_energyCut = new TGraph();
    TGraph* g_recHits_PU200_n_vs_RU_energyCut = new TGraph();
    TGraph* g_recHits_occ_vs_RU_energyCut = new TGraph();
    TGraph* g_recHits_PU200_occ_vs_RU_energyCut = new TGraph();
    
    for(int jj = 0; jj < 200; ++jj)
    {
      float cut = pow(10.,-2.+3./200.*jj);
      
      g_simHits_n_vs_RU_energyCut -> SetPoint(jj,cut,1.*simHits_n_vs_RU_energyCut[cut][iRU]/nEntries);
      g_simHits_PU200_n_vs_RU_energyCut -> SetPoint(jj,cut,200.*simHits_n_vs_RU_energyCut[cut][iRU]/nEntries);
      g_simHits_occ_vs_RU_energyCut -> SetPoint(jj,cut,1.*simHits_n_vs_RU_energyCut[cut][iRU]/nEntries/(nCrystalsPerMatrix*nMatricesPerMod*nModsPerType*3*nTrays/nRUs));
      g_simHits_PU200_occ_vs_RU_energyCut -> SetPoint(jj,cut,200.*simHits_n_vs_RU_energyCut[cut][iRU]/nEntries/(nCrystalsPerMatrix*nMatricesPerMod*nModsPerType*3*nTrays/nRUs));
      
      g_recHits_n_vs_RU_energyCut -> SetPoint(jj,cut,1.*recHits_n_vs_RU_energyCut[cut][iRU]/nEntries);
      g_recHits_PU200_n_vs_RU_energyCut -> SetPoint(jj,cut,200.*recHits_n_vs_RU_energyCut[cut][iRU]/nEntries);
      g_recHits_occ_vs_RU_energyCut -> SetPoint(jj,cut,1.*recHits_n_vs_RU_energyCut[cut][iRU]/nEntries/(nCrystalsPerMatrix*nMatricesPerMod*nModsPerType*3*nTrays/nRUs));
      g_recHits_PU200_occ_vs_RU_energyCut -> SetPoint(jj,cut,200.*recHits_n_vs_RU_energyCut[cut][iRU]/nEntries/(nCrystalsPerMatrix*nMatricesPerMod*nModsPerType*3*nTrays/nRUs));
    }
    
    g_simHits_n_vs_RU_energyCut -> Write(Form("g_simHits_n_vs_RU_energyCut_RU%d",iRU));
    g_simHits_PU200_n_vs_RU_energyCut -> Write(Form("g_simHits_PU200_n_vs_RU_energyCut_RU%d",iRU));
    g_simHits_occ_vs_RU_energyCut -> Write(Form("g_simHits_occ_vs_RU_energyCut_RU%d",iRU));
    g_simHits_PU200_occ_vs_RU_energyCut -> Write(Form("g_simHits_PU200_occ_vs_RU_energyCut_RU%d",iRU));
    
    g_recHits_n_vs_RU_energyCut -> Write(Form("g_recHits_n_vs_RU_energyCut_RU%d",iRU));
    g_recHits_PU200_n_vs_RU_energyCut -> Write(Form("g_recHits_PU200_n_vs_RU_energyCut_RU%d",iRU));
    g_recHits_occ_vs_RU_energyCut -> Write(Form("g_recHits_occ_vs_RU_energyCut_RU%d",iRU));
    g_recHits_PU200_occ_vs_RU_energyCut -> Write(Form("g_recHits_PU200_occ_vs_RU_energyCut_RU%d",iRU));
  }
  
  // int n_ieta = nZCrystalsPerMod*nModsPerType*3;
  // int n_iphi = nPhiCrystalsPerMod*nTrays/2;

  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;

}
