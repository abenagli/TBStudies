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
#include "TMVA/Reader.h"



int main(int argc, char** argv)
{
  setTDRStyle();

  if( argc < 2 )
  {
    std::cout << ">>> studyPhotonID::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  

  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  
  
  
  //--- get parameters
  std::vector<std::string> inFileNames = opts.GetOpt<std::vector<std::string> >("Input.fileNames");
    
  std::vector<float> ptRanges = opts.GetOpt<std::vector<float> >("Options.ptRanges");
  float ptRangesMin;
  float ptRangesMax;
  
  std::vector<float> etaRanges = opts.GetOpt<std::vector<float> >("Options.etaRanges");
  float etaRangesMin;
  float etaRangesMax;
  
  
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
  std::string plotLabel = opts.GetOpt<std::string>("Input.plotLabel");  

  
  
  //--- open output file
  TFile* outFile = TFile::Open(Form("%s/IDPlots_%s.root",plotDir.c_str(),plotLabel.c_str()),"RECREATE");
  outFile -> cd();  
  
  float eventId = 0.;
  float isSig;
  float isEB;
  float mva;
  float nRecHits;
  float energySum;
  float energySumCorr;
  float energySeed;
  float energyRatio;
  float sieie;
  float sipip;
  TTree* outTree = new TTree("outTree","outTree");
  outTree -> Branch("eventId",    &eventId,        "eventId/F");
  outTree -> Branch("isSig",      &isSig,            "isSig/F");
  outTree -> Branch("isEB",       &isEB,              "isEB/F");
  outTree -> Branch("mva",        &mva,                "mva/F");
  outTree -> Branch("nRecHits",   &nRecHits,      "nRecHits/I");
  outTree -> Branch("energySum",  &energySum,    "energySum/F");
  outTree -> Branch("energySeed", &energySeed,  "energySeed/F");
  outTree -> Branch("energyRatio",&energyRatio,"energyRatio/F");
  outTree -> Branch("sieie",      &sieie,            "sieie/F");
  outTree -> Branch("sipip",      &sipip,            "sipip/F");
  
  
  
  //-----
  // TMVA
  std::map<std::string,float*> varMap;
  varMap["eventId"] = &eventId;
  varMap["isSig"] = &isSig;
  varMap["isEB"] = &isEB;
  varMap["mva"] = &mva;
  varMap["nRecHits"] = &nRecHits;
  varMap["energySum"] = &energySum;
  varMap["energySumCorr"] = &energySumCorr;
  varMap["energySeed"] = &energySeed;
  varMap["energyRatio"] = &energyRatio;
  varMap["sieie"] = &sieie;
  varMap["sipip"] = &sipip;
  
  std::vector<std::string> diphoMVA_labels = opts.GetOpt<std::vector<std::string> >("Input.diphoMVA_labels");
  std::map<std::string,TMVA::Reader*> diphoMVAReaders;
  std::map<std::string,std::string> diphoMVA_methods;
  for(unsigned int ii = 0; ii < diphoMVA_labels.size(); ++ii)
  {
    std::string diphoMVA_label = diphoMVA_labels.at(ii);

    diphoMVA_methods[diphoMVA_label] = opts.GetOpt<std::string>(Form("Input.%s.method",diphoMVA_label.c_str()));
    std::string weightsFile = opts.GetOpt<std::string>(Form("Input.%s.weightsFile",diphoMVA_label.c_str()));
    std::vector<std::string> inputVariables = opts.GetOpt<std::vector<std::string> >(Form("Input.%s.inputVariables",diphoMVA_label.c_str()));

    diphoMVAReaders[diphoMVA_label] = new TMVA::Reader( "!Color:!Silent" );

    for(unsigned int jj = 0; jj < inputVariables.size(); ++jj)
    {
      std::string inputVariable = inputVariables.at(jj);
      diphoMVAReaders[diphoMVA_label] -> AddVariable(inputVariable.c_str(),varMap[inputVariable.c_str()]);
    }

    diphoMVAReaders[diphoMVA_label] -> BookMVA( diphoMVA_methods[diphoMVA_label],weightsFile.c_str() );
  }


  
  //--- loop over samples
  std::map<std::string,std::map<float,int> > ROC_nEvents_EB;
  std::map<std::string,std::map<float,int> > ROC2_nEvents_EB;

  std::map<std::string,std::map<float,int> > ROC_nEvents_EE;
  std::map<std::string,std::map<float,int> > ROC2_nEvents_EE;
  
  for(int inFileIt = 0; inFileIt < int(inFileNames.size()/2); ++inFileIt)
  {
    std::string inFileName = inFileNames.at(2*inFileIt+0);
    std::string label = inFileNames.at(2*inFileIt+1);
    std::cout << "\n\n>>> processing " << label << std::endl;

    
    //--- get tree
    TFile* inFile = TFile::Open(inFileName.c_str(),"READ");
    TTree* tree = (TTree*)( inFile->Get("FTLDumpPhotons/DumpPhotons") );
    tree -> SetBranchStatus("*",0);
    
    std::vector<float>* photons_pt  = new std::vector<float>;
    std::vector<float>* photons_eta = new std::vector<float>;
    std::vector<float>* photons_phi = new std::vector<float>;
    std::vector<float>* photons_mva = new std::vector<float>;
    std::vector<float>* photons_mcMatch_genPt = new std::vector<float>;
    std::vector<float>* photons_mcMatch_DR = new std::vector<float>;
    tree -> SetBranchStatus("photons_pt",           1); tree -> SetBranchAddress("photons_pt",           &photons_pt);
    tree -> SetBranchStatus("photons_eta",          1); tree -> SetBranchAddress("photons_eta",          &photons_eta);
    tree -> SetBranchStatus("photons_phi",          1); tree -> SetBranchAddress("photons_phi",          &photons_phi);
    tree -> SetBranchStatus("photons_mva",          1); tree -> SetBranchAddress("photons_mva",          &photons_mva);
    tree -> SetBranchStatus("photons_mcMatch_genPt",1); tree -> SetBranchAddress("photons_mcMatch_genPt",&photons_mcMatch_genPt);
    tree -> SetBranchStatus("photons_mcMatch_DR",   1); tree -> SetBranchAddress("photons_mcMatch_DR",   &photons_mcMatch_DR);
    
    std::vector<std::vector<int> >*   matchedRecHits_det = new std::vector<std::vector<int> >;
    std::vector<std::vector<float> >* matchedRecHits_time = new std::vector<std::vector<float> >;
    std::vector<std::vector<float> >* matchedRecHits_energy = new std::vector<std::vector<float> >;
    std::vector<std::vector<float> >* matchedRecHits_energyCorr = new std::vector<std::vector<float> >;
    std::vector<std::vector<float> >* matchedRecHits_photon_DR = new std::vector<std::vector<float> >;
    std::vector<std::vector<int> >*   matchedRecHits_modType = new std::vector<std::vector<int> >;
    std::vector<std::vector<int> >*   matchedRecHits_ieta = new std::vector<std::vector<int> >;
    std::vector<std::vector<int> >*   matchedRecHits_iphi = new std::vector<std::vector<int> >;
    tree -> SetBranchStatus("matchedRecHits_det",1);           tree -> SetBranchAddress("matchedRecHits_det",           &matchedRecHits_det);
    tree -> SetBranchStatus("matchedRecHits_time",1);          tree -> SetBranchAddress("matchedRecHits_time",          &matchedRecHits_time);
    tree -> SetBranchStatus("matchedRecHits_energy",1);        tree -> SetBranchAddress("matchedRecHits_energy",        &matchedRecHits_energy);
    tree -> SetBranchStatus("matchedRecHits_energyCorr",1);    tree -> SetBranchAddress("matchedRecHits_energyCorr",    &matchedRecHits_energyCorr);
    tree -> SetBranchStatus("matchedRecHits_photon_DR",1);   tree -> SetBranchAddress("matchedRecHits_photon_DR",   &matchedRecHits_photon_DR);
    tree -> SetBranchStatus("matchedRecHits_modType",1);       tree -> SetBranchAddress("matchedRecHits_modType",       &matchedRecHits_modType);
    tree -> SetBranchStatus("matchedRecHits_ieta",1);          tree -> SetBranchAddress("matchedRecHits_ieta",          &matchedRecHits_ieta);
    tree -> SetBranchStatus("matchedRecHits_iphi",1);          tree -> SetBranchAddress("matchedRecHits_iphi",          &matchedRecHits_iphi);
    

    outFile -> cd();
    
    TH1F* h1_photons_n = new TH1F(Form("h1_%s_photons_n",label.c_str()),"",250,-0.5,249.5);
    TH1F* h1_photons_pt = new TH1F(Form("h1_%s_photons_pt",label.c_str()),"",nPtBins,ptMin,ptMax);
    TH1F* h1_photons_eta = new TH1F(Form("h1_%s_photons_eta",label.c_str()),"",nEtaBins,etaMin,etaMax);
    TH1F* h1_photons_phi = new TH1F(Form("h1_%s_photons_phi",label.c_str()),"",nPhiBins,phiMin,phiMax);
    
    TH1F* h1_recHits_energy = new TH1F(Form("h1_%s_recHit_energy",label.c_str()),"",1000,0.,100.);
    TH1F* h1_recHits_time = new TH1F(Form("h1_%s_recHit_time",label.c_str()),"",500,-10.,40.);
    
    TH1F* h1_matchedRecHit_n_BTL = new TH1F(Form("h1_%s_matchedRecHit_n_BTL",label.c_str()),"",100,-0.5,99.5);
    TH1F* h1_matchedRecHit_energySumCorr_BTL = new TH1F(Form("h1_%s_matchedRecHit_energySumCorr_BTL",label.c_str()),"",1000,0.,100.);
    TH1F* h1_matchedRecHit_energySum_BTL = new TH1F(Form("h1_%s_matchedRecHit_energySum_BTL",label.c_str()),"",1000,0.,100.);
    TH1F* h1_matchedRecHit_energySeed_BTL = new TH1F(Form("h1_%s_matchedRecHit_energySeed_BTL",label.c_str()),"",1000,0.,100.);
    TH1F* h1_matchedRecHit_energyRatio_BTL = new TH1F(Form("h1_%s_matchedRecHit_energyRatio_BTL",label.c_str()),"",100.,0.,1.);
    TH1F* h1_matchedRecHit_energy_BTL = new TH1F(Form("h1_%s_matchedRecHit_energy_BTL",label.c_str()),"",1000,0.,100.);
    TH1F* h1_matchedRecHit_time_BTL = new TH1F(Form("h1_%s_matchedRecHit_time_BTL",label.c_str()),"",500,-10.,40.);
    TH1F* h1_matchedRecHit_photon_DR_BTL = new TH1F(Form("h1_%s_matchedRecHit_photon_DR_BTL",label.c_str()),"",10000,0.,10.);
    TH1F* h1_matchedRecHit_sieie_BTL = new TH1F(Form("h1_%s_matchedRecHit_sieie_BTL",label.c_str()),"",1000,0.,100.);
    TH1F* h1_matchedRecHit_sipip_BTL = new TH1F(Form("h1_%s_matchedRecHit_sipip_BTL",label.c_str()),"",1000,0.,100.);
    TProfile2D* p2_matchedRecHit_energy_vs_ieta_iphi_BTL = new TProfile2D(Form("p2_%s_matchedRecHit_energy_vs_ieta_iphi_BTL",label.c_str()),"",5,-2.5,2.5,5,-2.5,2.5);
    
    TH1F* h1_matchedRecHit_n_ETL = new TH1F(Form("h1_%s_matchedRecHit_n_ETL",label.c_str()),"",100,-0.5,99.5);
    TH1F* h1_matchedRecHit_energySumCorr_ETL = new TH1F(Form("h1_%s_matchedRecHit_energySumCorr_ETL",label.c_str()),"",1000,0.,10.);
    TH1F* h1_matchedRecHit_energySum_ETL = new TH1F(Form("h1_%s_matchedRecHit_energySum_ETL",label.c_str()),"",1000,0.,10.);
    TH1F* h1_matchedRecHit_energySeed_ETL = new TH1F(Form("h1_%s_matchedRecHit_energySeed_ETL",label.c_str()),"",1000,0.,10.);
    TH1F* h1_matchedRecHit_energyRatio_ETL = new TH1F(Form("h1_%s_matchedRecHit_energyRatio_ETL",label.c_str()),"",100.,0.,1.);
    TH1F* h1_matchedRecHit_energy_ETL = new TH1F(Form("h1_%s_matchedRecHit_energy_ETL",label.c_str()),"",1000,0.,10.);
    TH1F* h1_matchedRecHit_time_ETL = new TH1F(Form("h1_%s_matchedRecHit_time_ETL",label.c_str()),"",500,-10.,40.);
    TH1F* h1_matchedRecHit_photon_DR_ETL = new TH1F(Form("h1_%s_matchedRecHit_photon_DR_ETL",label.c_str()),"",10000,0.,10.);
    TH1F* h1_matchedRecHit_sieie_ETL = new TH1F(Form("h1_%s_matchedRecHit_sieie_ETL",label.c_str()),"",1000,0.,10.);
    TH1F* h1_matchedRecHit_sipip_ETL = new TH1F(Form("h1_%s_matchedRecHit_sipip_ETL",label.c_str()),"",1000,0.,10.);
    
    TProfile* p1_matchedRecHit_time_vs_eta = new TProfile(Form("p1_%s_matchedRecHit_time_vs_eta",label.c_str()),"",nEtaBins,etaMin,etaMax);
    
    TH1F* h1_nEntries_ptRanges = new TH1F(Form("h1_%s_nEntries_ptRanges",label.c_str()),"",ptRanges.size()-1,ptRanges.data());
    TH1F* h1_nEntries_etaRanges = new TH1F(Form("h1_%s_nEntries_etaRanges",label.c_str()),"",etaRanges.size()-1,etaRanges.data());
    
    
    
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
      for(unsigned int eleIt = 0; eleIt < photons_pt->size(); ++eleIt)
      {
        float pt = photons_pt->at(eleIt);
        float eta = photons_eta->at(eleIt);
        
        mva = photons_mva->at(eleIt);
        if( mva < -1 ) continue;
        
        isSig = label=="ele" ? 1 : 0;
        
        if( pt < 20. ) continue;
        
        ++nGoodTracks;
        h1_photons_pt -> Fill( photons_pt->at(eleIt) );
        h1_photons_eta -> Fill( fabs(photons_eta->at(eleIt)) );
        h1_photons_phi -> Fill( photons_phi->at(eleIt) );
        
        int ptBin = h1_nEntries_ptRanges -> Fill( pt );
        if( ptBin < 1 || ptBin > int(ptRanges.size()) ) continue;

        int etaBin = h1_nEntries_etaRanges -> Fill( fabs(eta) );
        if( etaBin < 1 || etaBin > int(etaRanges.size()) ) continue;

        
        // BTL
        if( fabs(eta) < 1.5 )
        {
          isEB = 1;

          nRecHits = 0;
          energySum = 0.;
          energySumCorr = 0.;
          energySeed = -1.;
          int ietaSeed = 0;
          int iphiSeed = 0;
          float ietaAvg = 0;
          float iphiAvg = 0;
          for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(eleIt)).size(); ++recHitIt)
          {
            if( (matchedRecHits_det->at(eleIt)).at(recHitIt) != 1 ) continue;
            
            float recHitE = (matchedRecHits_energy->at(eleIt)).at(recHitIt);
            h1_matchedRecHit_energy_BTL -> Fill( recHitE );
            
            if( recHitE < 1. ) continue;

            float DR = (matchedRecHits_photon_DR->at(eleIt)).at(recHitIt);
            h1_matchedRecHit_photon_DR_BTL -> Fill( DR );
            if( DR > 0.03 ) continue;

            ++nRecHits;
            
            if( recHitE > energySeed )
            {
              energySeed = recHitE;
              ietaSeed = (matchedRecHits_ieta->at(eleIt)).at(recHitIt);
              iphiSeed = (matchedRecHits_iphi->at(eleIt)).at(recHitIt);
            }
            
            h1_matchedRecHit_time_BTL -> Fill( (matchedRecHits_time->at(eleIt)).at(recHitIt) );

            p1_matchedRecHit_time_vs_eta -> Fill( fabs(eta),(matchedRecHits_time->at(eleIt)).at(recHitIt) );
            
            energySum += recHitE;
            if( (matchedRecHits_modType->at(eleIt)).at(recHitIt) == 1 )
              energySumCorr += (matchedRecHits_energyCorr->at(eleIt)).at(recHitIt);
            if( (matchedRecHits_modType->at(eleIt)).at(recHitIt) == 2 )
              energySumCorr += (matchedRecHits_energyCorr->at(eleIt)).at(recHitIt) * 3.75/3.00;
            if( (matchedRecHits_modType->at(eleIt)).at(recHitIt) == 3 )
              energySumCorr += (matchedRecHits_energyCorr->at(eleIt)).at(recHitIt) * 3.75/2.40;

            ietaAvg += (matchedRecHits_ieta->at(eleIt)).at(recHitIt)*recHitE;
            iphiAvg += (matchedRecHits_iphi->at(eleIt)).at(recHitIt)*recHitE;
          }
          energyRatio = energySum > 0. ? energySeed/energySum : -1.;

          if( energySum <= 0. ) continue;
          h1_matchedRecHit_n_BTL -> Fill ( nRecHits );          
          h1_matchedRecHit_energySum_BTL -> Fill( energySum );
          h1_matchedRecHit_energySumCorr_BTL -> Fill( energySumCorr );
          h1_matchedRecHit_energySeed_BTL -> Fill( energySeed );
          h1_matchedRecHit_energyRatio_BTL -> Fill( energySeed/energySum );
          
          sieie = 0.;
          sipip = 0.;
          ietaAvg /= energySum;
          iphiAvg /= energySum;

          for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(eleIt)).size(); ++recHitIt)
          {
            if( (matchedRecHits_det->at(eleIt)).at(recHitIt) != 1 ) continue;

            int ieta = (matchedRecHits_ieta->at(eleIt)).at(recHitIt);
            int iphi = (matchedRecHits_iphi->at(eleIt)).at(recHitIt);
            float recHitE = (matchedRecHits_energy->at(eleIt)).at(recHitIt);
            if( recHitE < 1. ) continue;

            float DR = (matchedRecHits_photon_DR->at(eleIt)).at(recHitIt);
            if( DR > 0.03 ) continue;
            
            sieie += recHitE * pow(ieta-ietaAvg,2);
            sipip += recHitE * pow(iphi-iphiAvg,2);
            
            p2_matchedRecHit_energy_vs_ieta_iphi_BTL -> Fill(iphi-iphiSeed,ieta-ietaSeed,recHitE/energySum);
          }
          sieie = energySum > 0. ? sqrt(sieie/energySum) : -1.;
          sipip = energySum > 0. ? sqrt(sipip/energySum) : -1.;
          
          h1_matchedRecHit_sieie_BTL -> Fill(sieie);
          h1_matchedRecHit_sipip_BTL -> Fill(sipip);
          
          
          // evaluate MVA
          for(unsigned int ii = 0; ii < diphoMVA_labels.size(); ++ii)
          {
            std::string diphoMVA_label = diphoMVA_labels.at(ii);
            float mva_new = diphoMVAReaders[diphoMVA_label] -> EvaluateMVA(diphoMVA_methods[diphoMVA_label].c_str());
            
            for(float val = -1.; val <= 1.; val += 0.01)
            {
              if( fabs(eta) < 1.5 )
              {
                if( photons_mva->at(eleIt) >= val )
                  ROC_nEvents_EB[label][val] += 1;
                if( mva_new >= val )
                  ROC2_nEvents_EB[label][val] += 1;
              }
            }
          }
          
          outTree -> Fill();
        } // BTL
        
        
        
        // ETL
        if( fabs(eta) > 1.5 )
        {
          isEB = 0;
          
          nRecHits = 0;
          energySum = 0.;
          energySumCorr = 0.;
          energySeed = -1.;
          int ietaSeed = 0;
          int iphiSeed = 0;
          float ietaAvg = 0;
          float iphiAvg = 0;
          for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(eleIt)).size(); ++recHitIt)
          {
            if( (matchedRecHits_det->at(eleIt)).at(recHitIt) != 2 ) continue;
            
            float recHitE = (matchedRecHits_energy->at(eleIt)).at(recHitIt);
            h1_matchedRecHit_energy_ETL -> Fill( recHitE );
            
            if( recHitE < 0. ) continue;

            float DR = (matchedRecHits_photon_DR->at(eleIt)).at(recHitIt);
            h1_matchedRecHit_photon_DR_ETL -> Fill( DR );
            if( DR > 0.02 ) continue;

            ++nRecHits;
            
            if( recHitE > energySeed )
            {
              energySeed = recHitE;
              ietaSeed = (matchedRecHits_ieta->at(eleIt)).at(recHitIt);
              iphiSeed = (matchedRecHits_iphi->at(eleIt)).at(recHitIt);
            }
            
            h1_matchedRecHit_time_ETL -> Fill( (matchedRecHits_time->at(eleIt)).at(recHitIt) );

            p1_matchedRecHit_time_vs_eta -> Fill( fabs(eta),(matchedRecHits_time->at(eleIt)).at(recHitIt) );
            
            energySum += recHitE;
            energySumCorr += recHitE;

            ietaAvg += (matchedRecHits_ieta->at(eleIt)).at(recHitIt)*recHitE;
            iphiAvg += (matchedRecHits_iphi->at(eleIt)).at(recHitIt)*recHitE;
          }
          energyRatio = energySum > 0. ? energySeed/energySum : -1.;
          
          h1_matchedRecHit_n_ETL -> Fill ( nRecHits );          
          h1_matchedRecHit_energySum_ETL -> Fill( energySum );
          h1_matchedRecHit_energySumCorr_ETL -> Fill( energySumCorr );
          h1_matchedRecHit_energySeed_ETL -> Fill( energySeed );
          h1_matchedRecHit_energyRatio_ETL -> Fill( energySeed/energySum );
          
          sieie = 0.;
          sipip = 0.;
          ietaAvg /= energySum;
          iphiAvg /= energySum;

          for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(eleIt)).size(); ++recHitIt)
          {
            if( (matchedRecHits_det->at(eleIt)).at(recHitIt) != 2 ) continue;

            int ieta = (matchedRecHits_ieta->at(eleIt)).at(recHitIt);
            int iphi = (matchedRecHits_iphi->at(eleIt)).at(recHitIt);
            float recHitE = (matchedRecHits_energy->at(eleIt)).at(recHitIt);
            if( recHitE < 0. ) continue;

            float DR = (matchedRecHits_photon_DR->at(eleIt)).at(recHitIt);
            if( DR > 0.02 ) continue;
            
            sieie += recHitE * pow(ieta-ietaAvg,2);
            sipip += recHitE * pow(iphi-iphiAvg,2);
          }
          sieie = energySum > 0. ? sqrt(sieie/energySum) : -1.;
          sipip = energySum > 0. ? sqrt(sipip/energySum) : -1.;
          
          h1_matchedRecHit_sieie_ETL -> Fill(sieie);
          h1_matchedRecHit_sipip_ETL -> Fill(sipip);

          
          // evaluate MVA
          for(unsigned int ii = 0; ii < diphoMVA_labels.size(); ++ii)
          {
            std::string diphoMVA_label = diphoMVA_labels.at(ii);
            float mva_new = diphoMVAReaders[diphoMVA_label] -> EvaluateMVA(diphoMVA_methods[diphoMVA_label].c_str());
            
            for(float val = -1.; val <= 1.; val += 0.01)
            {
              if( fabs(eta) > 1.5 )
              {
                if( photons_mva->at(eleIt) >= val )
                  ROC_nEvents_EE[label][val] += 1;
                if( mva_new >= val )
                  ROC2_nEvents_EE[label][val] += 1;
              }
            }
          }
          
          outTree -> Fill();          
        } // ETL

        ++eventId;
        
      } // end loop over tracks
      
      h1_photons_n -> Fill( nGoodTracks );
      
    } // end loop over events
    
  } // end loop over samples

  

  //--- make ROC plots
  TGraph* g_ROC_EB = new TGraph();
  
  std::map<float,int> ROC_nEvents_EB_ele = ROC_nEvents_EB["ele"];
  std::map<float,int> ROC_nEvents_EB_pi  = ROC_nEvents_EB["pi"];
  for(std::map<float,int>::const_iterator mapIt = ROC_nEvents_EB_ele.begin(); mapIt != ROC_nEvents_EB_ele.end(); ++mapIt)
  {
    float val = mapIt->first;
    g_ROC_EB -> SetPoint(g_ROC_EB->GetN(),1.*ROC_nEvents_EB_ele[val]/ROC_nEvents_EB_ele[-1.],1.*ROC_nEvents_EB_pi[val]/ROC_nEvents_EB_pi[-1.]);
  }

  g_ROC_EB -> Write("g_ROC_EB");

  
  TGraph* g_ROC2_EB = new TGraph();
    
  std::map<float,int> ROC2_nEvents_EB_ele = ROC2_nEvents_EB["ele"];
  std::map<float,int> ROC2_nEvents_EB_pi  = ROC2_nEvents_EB["pi"];
  for(std::map<float,int>::const_iterator mapIt = ROC2_nEvents_EB_ele.begin(); mapIt != ROC2_nEvents_EB_ele.end(); ++mapIt)
  {
    float val = mapIt->first;
    g_ROC2_EB -> SetPoint(g_ROC2_EB->GetN(),1.*ROC2_nEvents_EB_ele[val]/ROC2_nEvents_EB_ele[-1.],1.*ROC2_nEvents_EB_pi[val]/ROC2_nEvents_EB_pi[-1.]);
  }

  g_ROC2_EB -> Write("g_ROC2_EB");


  
  TGraph* g_ROC_EE = new TGraph();
  
  std::map<float,int> ROC_nEvents_EE_ele = ROC_nEvents_EE["ele"];
  std::map<float,int> ROC_nEvents_EE_pi  = ROC_nEvents_EE["pi"];
  for(std::map<float,int>::const_iterator mapIt = ROC_nEvents_EE_ele.begin(); mapIt != ROC_nEvents_EE_ele.end(); ++mapIt)
  {
    float val = mapIt->first;
    g_ROC_EE -> SetPoint(g_ROC_EE->GetN(),1.*ROC_nEvents_EE_ele[val]/ROC_nEvents_EE_ele[-1.],1.*ROC_nEvents_EE_pi[val]/ROC_nEvents_EE_pi[-1.]);
  }
  
  g_ROC_EE -> Write("g_ROC_EE");
  
  TGraph* g_ROC2_EE = new TGraph();
  
  std::map<float,int> ROC2_nEvents_EE_ele = ROC2_nEvents_EE["ele"];
  std::map<float,int> ROC2_nEvents_EE_pi  = ROC2_nEvents_EE["pi"];
  for(std::map<float,int>::const_iterator mapIt = ROC2_nEvents_EE_ele.begin(); mapIt != ROC2_nEvents_EE_ele.end(); ++mapIt)
  {
    float val = mapIt->first;
    g_ROC2_EE -> SetPoint(g_ROC2_EE->GetN(),1.*ROC2_nEvents_EE_ele[val]/ROC2_nEvents_EE_ele[-1.],1.*ROC2_nEvents_EE_pi[val]/ROC2_nEvents_EE_pi[-1.]);
  }

  g_ROC2_EE -> Write("g_ROC2_EE");
  

  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
}
