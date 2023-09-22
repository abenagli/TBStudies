#include "interface/SetTDRStyle.h"

#include <iostream>
#include <vector>
#include <map>

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"



int main(int argc, char** argv)
{
  setTDRStyle();
  
  if( argc < 2 )
  {
    std::cout << ">>> drawMapping.exe::usage:   ./bin/drawMapping.exe   crysLayout [1=default,2=thicker,3=thickest]" << std::endl;
    return -1;
  }
  int crysLayout = atoi(argv[1]);
  

  TChain* tree = new TChain("FTLDumpHits/DumpHits","FTLDumpHits/DumpHits");
  
  if( crysLayout == 1 ) tree -> Add("/Users/abenagli/Work/TIMING/TBStudies/data/TDR/CMSSW_13_0_0_pre1/ntuple_MinBias_14TeV_default_MinBias_14TeV.root");
  if( crysLayout == 2 ) tree -> Add("/Users/abenagli/Work/TIMING/TBStudies/data/TDR/CMSSW_13_0_0_pre1/ntuple_MinBias_14TeV_thicker_MinBias_14TeV.root");
  if( crysLayout == 3 ) tree -> Add("/Users/abenagli/Work/TIMING/TBStudies/data/TDR/CMSSW_13_0_0_pre1/ntuple_MinBias_14TeV_thickest_MinBias_14TeV.root");
  
  
  tree -> SetBranchStatus("*",0);
  
  std::vector<float>* tracks_pt  = new std::vector<float>;
  std::vector<float>* tracks_eta = new std::vector<float>;
  std::vector<float>* tracks_phi = new std::vector<float>;
  std::vector<float>* tracks_mcMatch_genPt = new std::vector<float>;
  std::vector<float>* tracks_mcMatch_DR = new std::vector<float>;
  tree -> SetBranchStatus("track_pt",           1); tree -> SetBranchAddress("track_pt",           &tracks_pt);
  tree -> SetBranchStatus("track_eta_atBTL",    1); tree -> SetBranchAddress("track_eta_atBTL",    &tracks_eta);
  tree -> SetBranchStatus("track_phi",          1); tree -> SetBranchAddress("track_phi",          &tracks_phi);
  tree -> SetBranchStatus("track_mcMatch_genPt",1); tree -> SetBranchAddress("track_mcMatch_genPt",&tracks_mcMatch_genPt);
  tree -> SetBranchStatus("track_mcMatch_DR",   1); tree -> SetBranchAddress("track_mcMatch_DR",   &tracks_mcMatch_DR);
  
  std::vector<std::vector<int> >*   matchedRecHits_det = new std::vector<std::vector<int> >;
  std::vector<std::vector<float> >* matchedRecHits_energy = new std::vector<std::vector<float> >;
  std::vector<std::vector<int> >* matchedRecHits_runit = new std::vector<std::vector<int> >;
  std::vector<std::vector<int> >* matchedRecHits_modType = new std::vector<std::vector<int> >;
  std::vector<std::vector<int> >* matchedRecHits_module = new std::vector<std::vector<int> >;
  tree -> SetBranchStatus("matchedRecHits_det",1);     tree -> SetBranchAddress("matchedRecHits_det",     &matchedRecHits_det);
  tree -> SetBranchStatus("matchedRecHits_energy",1);  tree -> SetBranchAddress("matchedRecHits_energy",  &matchedRecHits_energy);
  tree -> SetBranchStatus("matchedRecHits_runit",1);   tree -> SetBranchAddress("matchedRecHits_runit",   &matchedRecHits_runit);
  tree -> SetBranchStatus("matchedRecHits_modType",1); tree -> SetBranchAddress("matchedRecHits_modType", &matchedRecHits_modType);
  tree -> SetBranchStatus("matchedRecHits_module",1);  tree -> SetBranchAddress("matchedRecHits_module",  &matchedRecHits_module);
  
  
  //--- open output file
  TFile* outFile;
  if( crysLayout == 1 ) outFile = TFile::Open(Form("./plots/channelsPlots_default.root"),  "RECREATE");
  if( crysLayout == 2 ) outFile = TFile::Open(Form("./plots/channelsPlots_thicker.root"),  "RECREATE");
  if( crysLayout == 3 ) outFile = TFile::Open(Form("./plots/channelsPlots_thickest.root"), "RECREATE");
  outFile -> cd();
  
  TProfile* p_RUId_vs_eta = new TProfile("p_RUId_vs_eta","",10000,0.,2.);
  TProfile* p_runit_vs_eta = new TProfile("p_runit_vs_eta","",10000,0.,2.);
  TProfile* p_modType_vs_eta = new TProfile("p_modType_vs_eta","",10000,0.,2.);
  TProfile* p_moduleId_vs_eta = new TProfile("p_moduleId_vs_eta","",10000,0.,2.);
  
  
  
  //--- loop over events
  int nEntries = tree -> GetEntries();
  for(int entry = 0; entry < 50000; ++entry)
  {
    if( entry%10000 == 0 ) std::cout << ">>> reading entry " << entry << " / " << nEntries << "\r" << std::endl;
    tree -> GetEntry(entry);
    
    for(unsigned int trackIt = 0; trackIt < tracks_pt->size(); ++trackIt)
    {
      float pt = tracks_pt->at(trackIt);
      float feta = fabs(tracks_eta->at(trackIt));
      float genPt = tracks_mcMatch_genPt->at(trackIt);
      float DR = tracks_mcMatch_DR->at(trackIt);
      
      if( DR > 0.01 ) continue;
      if( fabs(pt/genPt-1.) > 0.05 ) continue;
      
      
      int seedIt = -1;
      float seedE = -1.;
      float totEnergy = 0;
      for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(trackIt)).size(); ++recHitIt)
      {
        if( (matchedRecHits_det->at(trackIt)).at(recHitIt) != 1 ) continue;
        float recHitE = (matchedRecHits_energy->at(trackIt)).at(recHitIt);
        if( recHitE < 1. ) continue;
        
        if( recHitE > seedE )
        {
          seedE = recHitE;
          seedIt = recHitIt;
        }
        
        totEnergy += recHitE;
      }
      
      if( totEnergy < 1. ) continue;
      
      int modType = (matchedRecHits_modType->at(trackIt)).at(seedIt);
      int runit  = (matchedRecHits_runit->at(trackIt)).at(seedIt);
      int RUId = crysLayout != 2 ? runit + (modType-1)*2 : (runit-1) + (modType-1)*3 + 1;
      int moduleId = (matchedRecHits_module->at(trackIt)).at(seedIt)-1;

      bool keepCluster = true;
      for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(trackIt)).size(); ++recHitIt)
      {
        if( (matchedRecHits_det->at(trackIt)).at(recHitIt) != 1 ) continue;
        float recHitE = (matchedRecHits_energy->at(trackIt)).at(recHitIt);
        if( recHitE < 1. ) continue;

        int modType_temp = (matchedRecHits_modType->at(trackIt)).at(recHitIt);
        int runit_temp  = (matchedRecHits_runit->at(trackIt)).at(recHitIt);
        int RUId_temp = crysLayout != 2 ? runit_temp + (modType_temp-1)*2 : (runit_temp-1) + (modType_temp-1)*3 + 1;
        
        if( RUId_temp != RUId ) keepCluster = false;
      }

      if( !keepCluster ) continue;
      
      p_RUId_vs_eta -> Fill(feta,RUId);
      p_moduleId_vs_eta -> Fill(feta,int(moduleId/3));
      p_runit_vs_eta -> Fill(feta,runit);
      p_modType_vs_eta -> Fill(feta,modType);
    }
  }
  
  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
  return 0;
}
