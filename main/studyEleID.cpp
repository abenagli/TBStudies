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
    std::cout << ">>> studyEleID::usage:   " << argv[0] << " configFile.cfg" << std::endl;
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
  float fabsEta;
  
  TTree* outTree = new TTree("outTree","outTree");
  outTree -> Branch("eventId",    &eventId,        "eventId/F");
  outTree -> Branch("isSig",      &isSig,            "isSig/F");
  outTree -> Branch("isEB",       &isEB,              "isEB/F");
  outTree -> Branch("mva",        &mva,                "mva/F");
  outTree -> Branch("nRecHits",   &nRecHits,      "nRecHits/F");
  outTree -> Branch("energySum",  &energySum,    "energySum/F");
  outTree -> Branch("energySeed", &energySeed,  "energySeed/F");
  outTree -> Branch("energyRatio",&energyRatio,"energyRatio/F");
  outTree -> Branch("sieie",      &sieie,            "sieie/F");
  outTree -> Branch("sipip",      &sipip,            "sipip/F");
  outTree -> Branch("fabsEta",    &fabsEta,        "fabsEta/F");  
  
  float ele_fbrem;
  float ele_gsfHits;
  float ele_gsfChi2;
  float ele_ctfHits;
  float ele_ctfChi2;
  float ele_sieie;
  float ele_sipip;
  float ele_circ;
  float ele_r9;
  float ele_etaw;
  float ele_phiw;
  float ele_hoe;
  float ele_eop;
  float ele_eseedopout;
  float ele_detain;
  float ele_dphiin;
  float ele_trackreliso;
  outTree -> Branch("ele_fbrem",      &ele_fbrem,            "ele_fbrem/F");
  outTree -> Branch("ele_gsfHits",    &ele_gsfHits,        "ele_gsfHits/F");
  outTree -> Branch("ele_gsfChi2",    &ele_gsfChi2,        "ele_gsfChi2/F");
  outTree -> Branch("ele_ctfHits",    &ele_ctfHits,        "ele_ctfHits/F");
  outTree -> Branch("ele_ctfChi2",    &ele_ctfChi2,        "ele_ctfChi2/F");
  outTree -> Branch("ele_sieie",      &ele_sieie,            "ele_sieie/F");
  outTree -> Branch("ele_sipip",      &ele_sipip,            "ele_sipip/F");
  outTree -> Branch("ele_circ",       &ele_circ,              "ele_circ/F");
  outTree -> Branch("ele_r9",         &ele_r9,                  "ele_r9/F");
  outTree -> Branch("ele_etaw",       &ele_etaw,              "ele_etaw/F");
  outTree -> Branch("ele_phiw",       &ele_phiw,              "ele_phiw/F");
  outTree -> Branch("ele_hoe",        &ele_hoe,                "ele_hoe/F");
  outTree -> Branch("ele_eop",        &ele_eop,                "ele_eop/F");
  outTree -> Branch("ele_eseedopout", &ele_eseedopout,      "ele_eopout/F");
  outTree -> Branch("ele_detain",     &ele_detain,          "ele_detain/F");
  outTree -> Branch("ele_dphiin",     &ele_dphiin,          "ele_dphiin/F");
  outTree -> Branch("ele_trackreliso",&ele_trackreliso,"ele_trackreliso/F");
  
  
  
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
  varMap["fabsEta"] = &fabsEta;
  varMap["ele_fbrem"] = &ele_fbrem;
  varMap["ele_gsfHits"] = &ele_gsfHits;
  varMap["ele_gsfChi2"] = &ele_gsfChi2;
  varMap["ele_ctfHits"] = &ele_ctfHits;
  varMap["ele_ctfChi2"] = &ele_ctfChi2;
  varMap["ele_sieie"] = &ele_sieie;
  varMap["ele_sipip"] = &ele_sipip;
  varMap["ele_circ"] = &ele_circ;
  varMap["ele_r9"] = &ele_r9;
  varMap["ele_etaw"] = &ele_etaw;
  varMap["ele_phiw"] = &ele_phiw;
  varMap["ele_hoe"] = &ele_hoe;
  varMap["ele_eop"] = &ele_eop;
  varMap["ele_eseedopout"] = &ele_eseedopout;
  varMap["ele_detain"] = &ele_detain;
  varMap["ele_dphiin"] = &ele_dphiin;
  varMap["ele_trackreliso"] = &ele_trackreliso;
  
  std::vector<std::string> MVA_labels = opts.GetOpt<std::vector<std::string> >("Input.MVA_labels");
  std::map<std::string,TMVA::Reader*> MVAReaders;
  std::map<std::string,std::string> MVA_methods;
  for(unsigned int ii = 0; ii < MVA_labels.size(); ++ii)
  {
    std::string MVA_label = MVA_labels.at(ii);

    MVA_methods[MVA_label] = opts.GetOpt<std::string>(Form("Input.%s.method",MVA_label.c_str()));
    std::string weightsFile = opts.GetOpt<std::string>(Form("Input.%s.weightsFile",MVA_label.c_str()));
    std::vector<std::string> inputVariables = opts.GetOpt<std::vector<std::string> >(Form("Input.%s.inputVariables",MVA_label.c_str()));

    MVAReaders[MVA_label] = new TMVA::Reader( "!Color:!Silent" );

    for(unsigned int jj = 0; jj < inputVariables.size(); ++jj)
    {
      std::string inputVariable = inputVariables.at(jj);
      MVAReaders[MVA_label] -> AddVariable(inputVariable.c_str(),varMap[inputVariable.c_str()]);
    }

    MVAReaders[MVA_label] -> BookMVA( MVA_methods[MVA_label],weightsFile.c_str() );
  }


  
  //--- loop over samples
  std::vector<std::string> MVALabels_EB;
  MVALabels_EB.push_back("eleID_MTDOnly_EB");
  MVALabels_EB.push_back("eleID_eleOnlyNoIso_EB");
  MVALabels_EB.push_back("eleID_eleOnlyIso_EB");
  MVALabels_EB.push_back("eleID_eleMTDNoIso_EB");
  MVALabels_EB.push_back("eleID_eleMTDIso_EB");
  
  std::vector<std::string> MVALabels_EE;
  MVALabels_EE.push_back("eleID_MTDOnly_EE");
  
  std::map<std::string,std::map<std::string,std::map<float,int> > > ROC_nEvents_EB;
  std::map<std::string,std::map<std::string,std::map<float,int> > > ROC_nEvents_EE;
  
  
  for(int inFileIt = 0; inFileIt < int(inFileNames.size()/2); ++inFileIt)
  {
    std::string inFileName = inFileNames.at(2*inFileIt+0);
    std::string label = inFileNames.at(2*inFileIt+1);
    std::cout << "\n\n>>> processing " << label << std::endl;

    
    //--- get tree
    TChain* tree = new TChain("FTLDumpElectrons/DumpElectrons","FTLDumpElectrons/DumpElectrons");
    tree -> Add((inFileName+"*.root").c_str());
    
    tree -> SetBranchStatus("*",0);
    
    std::vector<float>* electrons_pt  = new std::vector<float>;
    std::vector<float>* electrons_eta = new std::vector<float>;
    std::vector<float>* electrons_phi = new std::vector<float>;
    std::vector<float>* electrons_mva = new std::vector<float>;
    std::vector<float>* electrons_t = new std::vector<float>;
    std::vector<float>* electrons_pathLength = new std::vector<float>;
    std::vector<float>* electrons_t_atBTL = new std::vector<float>;
    std::vector<float>* electrons_eta_atBTL = new std::vector<float>;
    std::vector<float>* electrons_phi_atBTL = new std::vector<float>;
    std::vector<float>* electrons_eta_atETL = new std::vector<float>;
    std::vector<float>* electrons_phi_atETL = new std::vector<float>;
    std::vector<float>* electrons_mcMatch_genPt = new std::vector<float>;
    std::vector<float>* electrons_mcMatch_DR = new std::vector<float>;
    std::vector<float>* electrons_trackIso = new std::vector<float>;
    tree -> SetBranchStatus("electrons_pt",           1); tree -> SetBranchAddress("electrons_pt",           &electrons_pt);
    tree -> SetBranchStatus("electrons_eta",          1); tree -> SetBranchAddress("electrons_eta",          &electrons_eta);
    tree -> SetBranchStatus("electrons_phi",          1); tree -> SetBranchAddress("electrons_phi",          &electrons_phi);
    tree -> SetBranchStatus("electrons_mva",          1); tree -> SetBranchAddress("electrons_mva",          &electrons_mva);
    tree -> SetBranchStatus("electrons_t",            1); tree -> SetBranchAddress("electrons_t",            &electrons_t);
    tree -> SetBranchStatus("electrons_pathLength",   1); tree -> SetBranchAddress("electrons_pathLength",   &electrons_pathLength);
    tree -> SetBranchStatus("electrons_t_atBTL",      1); tree -> SetBranchAddress("electrons_t_atBTL",      &electrons_t_atBTL);
    tree -> SetBranchStatus("electrons_eta_atBTL",    1); tree -> SetBranchAddress("electrons_eta_atBTL",    &electrons_eta_atBTL);
    tree -> SetBranchStatus("electrons_phi_atBTL",    1); tree -> SetBranchAddress("electrons_phi_atBTL",    &electrons_phi_atBTL);
    tree -> SetBranchStatus("electrons_eta_atETL",    1); tree -> SetBranchAddress("electrons_eta_atETL",    &electrons_eta_atETL);
    tree -> SetBranchStatus("electrons_phi_atETL",    1); tree -> SetBranchAddress("electrons_phi_atETL",    &electrons_phi_atETL);
    tree -> SetBranchStatus("electrons_mcMatch_genPt",1); tree -> SetBranchAddress("electrons_mcMatch_genPt",&electrons_mcMatch_genPt);
    tree -> SetBranchStatus("electrons_mcMatch_DR",   1); tree -> SetBranchAddress("electrons_mcMatch_DR",   &electrons_mcMatch_DR);
    tree -> SetBranchStatus("electrons_trackIso",     1); tree -> SetBranchAddress("electrons_trackIso",     &electrons_trackIso);
    
    std::vector<float>* electrons_fbrem = new std::vector<float>;
    std::vector<float>* electrons_gsfHits = new std::vector<float>;
    std::vector<float>* electrons_gsfChi2 = new std::vector<float>;
    std::vector<float>* electrons_ctfHits = new std::vector<float>;
    std::vector<float>* electrons_ctfChi2 = new std::vector<float>;
    std::vector<float>* electrons_dEtaIn = new std::vector<float>;
    std::vector<float>* electrons_dPhiIn = new std::vector<float>;
    std::vector<float>* electrons_eop = new std::vector<float>;
    std::vector<float>* electrons_eSeedOPout = new std::vector<float>;
    std::vector<float>* electrons_sigmaEtaEta = new std::vector<float>;
    std::vector<float>* electrons_sigmaIetaIeta = new std::vector<float>;
    std::vector<float>* electrons_sigmaIphiIphi = new std::vector<float>;
    std::vector<float>* electrons_e1x5 = new std::vector<float>;
    std::vector<float>* electrons_e2x5Max = new std::vector<float>;
    std::vector<float>* electrons_e5x5 = new std::vector<float>;
    std::vector<float>* electrons_r9 = new std::vector<float>;
    std::vector<float>* electrons_etaWidth = new std::vector<float>;
    std::vector<float>* electrons_phiWidth = new std::vector<float>;
    std::vector<float>* electrons_hcalOverEcal = new std::vector<float>;
    tree -> SetBranchStatus("electrons_fbrem",1); tree -> SetBranchAddress("electrons_fbrem",&electrons_fbrem);
    tree -> SetBranchStatus("electrons_gsfHits",1); tree -> SetBranchAddress("electrons_gsfHits",&electrons_gsfHits);
    tree -> SetBranchStatus("electrons_gsfChi2",1); tree -> SetBranchAddress("electrons_gsfChi2",&electrons_gsfChi2);
    tree -> SetBranchStatus("electrons_ctfHits",1); tree -> SetBranchAddress("electrons_ctfHits",&electrons_ctfHits);
    tree -> SetBranchStatus("electrons_ctfChi2",1); tree -> SetBranchAddress("electrons_ctfChi2",&electrons_ctfChi2);
    tree -> SetBranchStatus("electrons_dEtaIn",1); tree -> SetBranchAddress("electrons_dEtaIn",&electrons_dEtaIn);
    tree -> SetBranchStatus("electrons_dPhiIn",1); tree -> SetBranchAddress("electrons_dPhiIn",&electrons_dPhiIn);
    tree -> SetBranchStatus("electrons_eop",1); tree -> SetBranchAddress("electrons_eop",&electrons_eop);
    tree -> SetBranchStatus("electrons_eSeedOPout",1); tree -> SetBranchAddress("electrons_eSeedOPout",&electrons_eSeedOPout);
    tree -> SetBranchStatus("electrons_sigmaEtaEta",1); tree -> SetBranchAddress("electrons_sigmaEtaEta",&electrons_sigmaEtaEta);
    tree -> SetBranchStatus("electrons_sigmaIetaIeta",1); tree -> SetBranchAddress("electrons_sigmaIetaIeta",&electrons_sigmaIetaIeta);
    tree -> SetBranchStatus("electrons_sigmaIphiIphi",1); tree -> SetBranchAddress("electrons_sigmaIphiIphi",&electrons_sigmaIphiIphi);
    tree -> SetBranchStatus("electrons_e1x5",1); tree -> SetBranchAddress("electrons_e1x5",&electrons_e1x5);
    tree -> SetBranchStatus("electrons_e2x5Max",1); tree -> SetBranchAddress("electrons_e2x5Max",&electrons_e2x5Max);
    tree -> SetBranchStatus("electrons_e5x5",1); tree -> SetBranchAddress("electrons_e5x5",&electrons_e5x5);
    tree -> SetBranchStatus("electrons_r9",1); tree -> SetBranchAddress("electrons_r9",&electrons_r9);
    tree -> SetBranchStatus("electrons_etaWidth",1); tree -> SetBranchAddress("electrons_etaWidth",&electrons_etaWidth);
    tree -> SetBranchStatus("electrons_phiWidth",1); tree -> SetBranchAddress("electrons_phiWidth",&electrons_phiWidth);
    tree -> SetBranchStatus("electrons_hcalOverEcal",1); tree -> SetBranchAddress("electrons_hcalOverEcal",&electrons_hcalOverEcal);
    
    std::vector<std::vector<int> >*   matchedRecHits_det = new std::vector<std::vector<int> >;
    std::vector<std::vector<float> >* matchedRecHits_time = new std::vector<std::vector<float> >;
    std::vector<std::vector<float> >* matchedRecHits_energy = new std::vector<std::vector<float> >;
    std::vector<std::vector<float> >* matchedRecHits_energyCorr = new std::vector<std::vector<float> >;
    std::vector<std::vector<float> >* matchedRecHits_electron_DR = new std::vector<std::vector<float> >;
    std::vector<std::vector<int> >*   matchedRecHits_modType = new std::vector<std::vector<int> >;
    std::vector<std::vector<int> >*   matchedRecHits_ieta = new std::vector<std::vector<int> >;
    std::vector<std::vector<int> >*   matchedRecHits_iphi = new std::vector<std::vector<int> >;
    tree -> SetBranchStatus("matchedRecHits_det",1);           tree -> SetBranchAddress("matchedRecHits_det",           &matchedRecHits_det);
    tree -> SetBranchStatus("matchedRecHits_time",1);          tree -> SetBranchAddress("matchedRecHits_time",          &matchedRecHits_time);
    tree -> SetBranchStatus("matchedRecHits_energy",1);        tree -> SetBranchAddress("matchedRecHits_energy",        &matchedRecHits_energy);
    tree -> SetBranchStatus("matchedRecHits_energyCorr",1);    tree -> SetBranchAddress("matchedRecHits_energyCorr",    &matchedRecHits_energyCorr);
    tree -> SetBranchStatus("matchedRecHits_electron_DR",1);   tree -> SetBranchAddress("matchedRecHits_electron_DR",   &matchedRecHits_electron_DR);
    tree -> SetBranchStatus("matchedRecHits_modType",1);       tree -> SetBranchAddress("matchedRecHits_modType",       &matchedRecHits_modType);
    tree -> SetBranchStatus("matchedRecHits_ieta",1);          tree -> SetBranchAddress("matchedRecHits_ieta",          &matchedRecHits_ieta);
    tree -> SetBranchStatus("matchedRecHits_iphi",1);          tree -> SetBranchAddress("matchedRecHits_iphi",          &matchedRecHits_iphi);
    

    outFile -> cd();
    
    TH1F* h1_electrons_n = new TH1F(Form("h1_%s_electrons_n",label.c_str()),"",250,-0.5,249.5);
    TH1F* h1_electrons_pt = new TH1F(Form("h1_%s_electrons_pt",label.c_str()),"",nPtBins,ptMin,ptMax);
    TH1F* h1_electrons_eta = new TH1F(Form("h1_%s_electrons_eta",label.c_str()),"",nEtaBins,etaMin,etaMax);
    TH1F* h1_electrons_phi = new TH1F(Form("h1_%s_electrons_phi",label.c_str()),"",nPhiBins,phiMin,phiMax);
    TH1F* h1_electrons_eta_atBTL = new TH1F(Form("h1_%s_electrons_eta_atBTL",label.c_str()),"",nEtaBins,etaMin,etaMax);
    TH1F* h1_electrons_phi_atBTL = new TH1F(Form("h1_%s_electrons_phi_atBTL",label.c_str()),"",nPhiBins,phiMin,phiMax);
    TH1F* h1_electrons_eta_atETL = new TH1F(Form("h1_%s_electrons_eta_atETL",label.c_str()),"",nEtaBins,etaMin,etaMax);
    TH1F* h1_electrons_phi_atETL = new TH1F(Form("h1_%s_electrons_phi_atETL",label.c_str()),"",nPhiBins,phiMin,phiMax);
    TH1F* h1_electrons_t_BTL = new TH1F(Form("h1_%s_electrons_t_BTL",label.c_str()),"",1000.,-1.,1.);
    TH1F* h1_electrons_pathLength_BTL = new TH1F(Form("h1_%s_electrons_pathLength_BTL",label.c_str()),"",1000.,0.,500.);
    TH1F* h1_electrons_t_ETL = new TH1F(Form("h1_%s_electrons_t_ETL",label.c_str()),"",1000.,-1.,1.);
    TH1F* h1_electrons_pathLength_ETL = new TH1F(Form("h1_%s_electrons_pathLength_ETL",label.c_str()),"",1000.,0.,500.);
    TH1F* h1_electrons_trackRelIso_BTL = new TH1F(Form("h1_%s_electrons_trackRelIso_BTL",label.c_str()),"",10000.,0.,10.);
    TH1F* h1_electrons_trackRelIso_ETL = new TH1F(Form("h1_%s_electrons_trackRelIso_ETL",label.c_str()),"",10000.,0.,10.);
    
    TH1F* h1_recHits_energy = new TH1F(Form("h1_%s_recHit_energy",label.c_str()),"",1000,0.,100.);
    TH1F* h1_recHits_time = new TH1F(Form("h1_%s_recHit_time",label.c_str()),"",500,-10.,40.);
    
    TH1F* h1_electrons_mva_BTL = new TH1F(Form("h1_%s_electrons_mva_BTL",label.c_str()),"",1000,-1.,1.);
    TH1F* h1_matchedRecHit_n_BTL = new TH1F(Form("h1_%s_matchedRecHit_n_BTL",label.c_str()),"",100,-0.5,99.5);
    TH1F* h1_matchedRecHit_energySumCorr_BTL = new TH1F(Form("h1_%s_matchedRecHit_energySumCorr_BTL",label.c_str()),"",1000,0.,100.);
    TH1F* h1_matchedRecHit_energySum_BTL = new TH1F(Form("h1_%s_matchedRecHit_energySum_BTL",label.c_str()),"",1000,0.,100.);
    TH1F* h1_matchedRecHit_energySeed_BTL = new TH1F(Form("h1_%s_matchedRecHit_energySeed_BTL",label.c_str()),"",1000,0.,100.);
    TH1F* h1_matchedRecHit_energyRatio_BTL = new TH1F(Form("h1_%s_matchedRecHit_energyRatio_BTL",label.c_str()),"",200.,0.,2.);
    TH1F* h1_matchedRecHit_energy_BTL = new TH1F(Form("h1_%s_matchedRecHit_energy_BTL",label.c_str()),"",1000,0.,100.);
    TH1F* h1_matchedRecHit_time_BTL = new TH1F(Form("h1_%s_matchedRecHit_time_BTL",label.c_str()),"",2000,-10.,10.);
    TH1F* h1_matchedRecHit_electron_DR_BTL = new TH1F(Form("h1_%s_matchedRecHit_electron_DR_BTL",label.c_str()),"",10000,0.,10.);
    TH1F* h1_matchedRecHit_sieie_BTL = new TH1F(Form("h1_%s_matchedRecHit_sieie_BTL",label.c_str()),"",1000,0.,100.);
    TH1F* h1_matchedRecHit_sipip_BTL = new TH1F(Form("h1_%s_matchedRecHit_sipip_BTL",label.c_str()),"",1000,0.,100.);
    TProfile2D* p2_matchedRecHit_energy_vs_ieta_iphi_BTL = new TProfile2D(Form("p2_%s_matchedRecHit_energy_vs_ieta_iphi_BTL",label.c_str()),"",5,-2.5,2.5,5,-2.5,2.5);

    TH1F* h1_electrons_mva_ETL = new TH1F(Form("h1_%s_electrons_mva_ETL",label.c_str()),"",1000,-1.,1.);
    TH1F* h1_matchedRecHit_n_ETL = new TH1F(Form("h1_%s_matchedRecHit_n_ETL",label.c_str()),"",100,-0.5,99.5);
    TH1F* h1_matchedRecHit_energySumCorr_ETL = new TH1F(Form("h1_%s_matchedRecHit_energySumCorr_ETL",label.c_str()),"",1000,0.,10.);
    TH1F* h1_matchedRecHit_energySum_ETL = new TH1F(Form("h1_%s_matchedRecHit_energySum_ETL",label.c_str()),"",1000,0.,10.);
    TH1F* h1_matchedRecHit_energySeed_ETL = new TH1F(Form("h1_%s_matchedRecHit_energySeed_ETL",label.c_str()),"",1000,0.,10.);
    TH1F* h1_matchedRecHit_energyRatio_ETL = new TH1F(Form("h1_%s_matchedRecHit_energyRatio_ETL",label.c_str()),"",200.,0.,2.);
    TH1F* h1_matchedRecHit_energy_ETL = new TH1F(Form("h1_%s_matchedRecHit_energy_ETL",label.c_str()),"",1000,0.,10.);
    TH1F* h1_matchedRecHit_time_ETL = new TH1F(Form("h1_%s_matchedRecHit_time_ETL",label.c_str()),"",2000,-10.,10.);
    TH1F* h1_matchedRecHit_electron_DR_ETL = new TH1F(Form("h1_%s_matchedRecHit_electron_DR_ETL",label.c_str()),"",10000,0.,10.);
    TH1F* h1_matchedRecHit_sieie_ETL = new TH1F(Form("h1_%s_matchedRecHit_sieie_ETL",label.c_str()),"",1000,0.,1.);
    TH1F* h1_matchedRecHit_sipip_ETL = new TH1F(Form("h1_%s_matchedRecHit_sipip_ETL",label.c_str()),"",1000,0.,1.);
    TProfile2D* p2_matchedRecHit_energy_vs_ieta_iphi_ETL = new TProfile2D(Form("p2_%s_matchedRecHit_energy_vs_ieta_iphi_ETL",label.c_str()),"",5,-2.5,2.5,5,-2.5,2.5);
    
    TProfile* p1_matchedRecHit_time_vs_eta = new TProfile(Form("p1_%s_matchedRecHit_time_vs_eta",label.c_str()),"",nEtaBins,etaMin,etaMax);
    
    TH1F* h1_nEntries_ptRanges = new TH1F(Form("h1_%s_nEntries_ptRanges",label.c_str()),"",ptRanges.size()-1,ptRanges.data());
    TH1F* h1_nEntries_etaRanges = new TH1F(Form("h1_%s_nEntries_etaRanges",label.c_str()),"",etaRanges.size()-1,etaRanges.data());
    
    
    
    //--- loop over events
    int nEntries = tree->GetEntries();
    // nEntries = 10000;
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%100 == 0 ) std::cout << ">>> reading entry " << entry << " / " << nEntries << "\r" << std::flush;
      if( entry%2 != 0 ) continue;
      
      tree -> GetEntry(entry);

      
      //------------------------------
      //--- fill matched recHits plots
      if( debugMode ) std::cout << ">>> fill matched recHits plots" << std::endl;
      int nGoodTracks = 0;
      for(unsigned int eleIt = 0; eleIt < electrons_pt->size(); ++eleIt)
      {
        float pt = electrons_pt->at(eleIt);
        float eta = electrons_eta->at(eleIt);
        float eta_atBTL = electrons_eta_atBTL->at(eleIt);
        float phi_atBTL = electrons_phi_atBTL->at(eleIt);
        float eta_atETL = electrons_eta_atETL->at(eleIt);
        float phi_atETL = electrons_phi_atETL->at(eleIt);
        float mcMatch_DR = electrons_mcMatch_DR->at(eleIt);
        float t = electrons_t->at(eleIt);
        float t_atBTL = electrons_t_atBTL->at(eleIt);
        float pathLength = electrons_pathLength->at(eleIt);
        
        ele_fbrem = electrons_fbrem->at(eleIt);
        ele_gsfHits = electrons_gsfHits->at(eleIt);
        ele_gsfChi2 = electrons_gsfChi2->at(eleIt);
        ele_ctfHits = electrons_ctfHits->at(eleIt);
        ele_ctfChi2 = electrons_ctfChi2->at(eleIt);
        ele_sieie = electrons_sigmaIetaIeta->at(eleIt);
        ele_sipip = electrons_sigmaIphiIphi->at(eleIt);
        ele_circ = (electrons_e5x5->at(eleIt) != 0.) ? 1.-electrons_e1x5->at(eleIt)/electrons_e5x5->at(eleIt) : -1.;
        ele_r9 = electrons_r9->at(eleIt);
        ele_etaw = electrons_etaWidth->at(eleIt);
        ele_phiw = electrons_phiWidth->at(eleIt);
        ele_hoe = electrons_hcalOverEcal->at(eleIt);
        ele_eop = electrons_eop->at(eleIt);
        ele_eseedopout = electrons_eSeedOPout->at(eleIt);
        ele_detain = electrons_dEtaIn->at(eleIt);
        ele_dphiin = electrons_dPhiIn->at(eleIt);
        ele_trackreliso = electrons_trackIso->at(eleIt)/pt;
        if( t == -999. ) t_atBTL += 999.;
        
        fabsEta = fabs(eta);
        mva = electrons_mva->at(eleIt);
        if( mva < -1 ) continue;
        // if( eta_atBTL > -100. && mva <  0.6 ) continue;
        // if( eta_atETL > -100. && mva < -0.2 ) continue;
        
        isSig = (label=="prompt" || label=="promptNoPU") ? 1 : 0;
        
        if( pt < 20. ) continue;
        if( isSig && mcMatch_DR > 0.05) continue;
        if( !isSig && mcMatch_DR < 100. ) continue;
        
        if( eta_atBTL > -100. && ele_r9 > -100. )
        {
          h1_electrons_t_BTL -> Fill( electrons_t->at(eleIt));
          h1_electrons_pathLength_BTL -> Fill( electrons_pathLength->at(eleIt) );
          h1_electrons_trackRelIso_BTL -> Fill( electrons_trackIso->at(eleIt)/pt );
        }
        else
        {
          h1_electrons_t_ETL -> Fill( electrons_t->at(eleIt));
          h1_electrons_pathLength_ETL -> Fill( electrons_pathLength->at(eleIt));          
          h1_electrons_trackRelIso_ETL -> Fill( electrons_trackIso->at(eleIt)/pt );
        }
        
        // if( pathLength < 0. ) continue;
        // if( t < -100. ) continue;
        // if( electrons_trackIso->at(eleIt)/pt < 0.10 ) continue;
        
        ++nGoodTracks;
        h1_electrons_pt -> Fill( electrons_pt->at(eleIt) );
        h1_electrons_eta -> Fill( fabs(electrons_eta->at(eleIt)) );
        h1_electrons_phi -> Fill( electrons_phi->at(eleIt) );
        
        int ptBin = h1_nEntries_ptRanges -> Fill( pt );
        if( ptBin < 1 || ptBin > int(ptRanges.size()) ) continue;

        int etaBin = h1_nEntries_etaRanges -> Fill( fabs(eta) );
        if( etaBin < 1 || etaBin > int(etaRanges.size()) ) continue;
        
        // BTL
        if( eta_atBTL > -100. && ele_r9 > -100. )
        {
          isEB = 1;
          
          h1_electrons_eta_atBTL -> Fill( fabs(electrons_eta_atBTL->at(eleIt)) );
          h1_electrons_phi_atBTL -> Fill( electrons_phi_atBTL->at(eleIt) );
          h1_electrons_mva_BTL -> Fill( electrons_mva->at(eleIt) );
          
          energySeed = 0.;
          float timeSeed = 0.;
          int ietaSeed = 0;
          int iphiSeed = 0;
          for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(eleIt)).size(); ++recHitIt)
          {
            if( (matchedRecHits_det->at(eleIt)).at(recHitIt) != 1 ) continue;
            
            float DR = (matchedRecHits_electron_DR->at(eleIt)).at(recHitIt);
            if( DR > 0.03 ) continue;
            
            float recHitE = (matchedRecHits_energy->at(eleIt)).at(recHitIt);
            float recHitTime = (matchedRecHits_time->at(eleIt)).at(recHitIt);
            
            if( recHitE > energySeed )
            {
              energySeed = recHitE;
              timeSeed = recHitTime;
              ietaSeed = (matchedRecHits_ieta->at(eleIt)).at(recHitIt);
              iphiSeed = (matchedRecHits_iphi->at(eleIt)).at(recHitIt);
            }
          }
          
          nRecHits = 0.;
          energySum = 0.;
          energySumCorr = 0.;
          float ietaAvg = 0;
          float iphiAvg = 0;
          for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(eleIt)).size(); ++recHitIt)
          {
            if( (matchedRecHits_det->at(eleIt)).at(recHitIt) != 1 ) continue;
            
            float recHitE = (matchedRecHits_energy->at(eleIt)).at(recHitIt);
            h1_matchedRecHit_energy_BTL -> Fill( recHitE );
            
            if( recHitE < 0. ) continue;

            float recHitTime = (matchedRecHits_time->at(eleIt)).at(recHitIt);
            h1_matchedRecHit_time_BTL -> Fill( recHitTime-timeSeed );
            p1_matchedRecHit_time_vs_eta -> Fill( fabs(eta),recHitTime );
            
            if( fabs(recHitTime-timeSeed) > 0.1 ) continue;
            
            float DR = (matchedRecHits_electron_DR->at(eleIt)).at(recHitIt);
            h1_matchedRecHit_electron_DR_BTL -> Fill( DR );
            if( DR > 0.03 ) continue;

            nRecHits += 1.;
            
            energySum += recHitE;
            
            ietaAvg += (matchedRecHits_ieta->at(eleIt)).at(recHitIt)*recHitE;
            iphiAvg += (matchedRecHits_iphi->at(eleIt)).at(recHitIt)*recHitE;
          }
          energyRatio = energySum > 0. ? energySeed/energySum : -1.;
          
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
            float recHitTime = (matchedRecHits_time->at(eleIt)).at(recHitIt);
            if( recHitE < 0. ) continue;
            if( fabs(recHitTime-timeSeed) > 0.1 ) continue;
            
            float DR = (matchedRecHits_electron_DR->at(eleIt)).at(recHitIt);
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
          for(float val = -1.; val <= 1.; val += 0.005)
          {
            if( eta_atBTL > -100 && eta_atETL < -100 )
            {
              if( electrons_mva->at(eleIt) >= val )
                ROC_nEvents_EB["old"][label][val] += 1;
            }
          }
          
          for( auto MVA_label : MVALabels_EB )
          {
            float mva_new = MVAReaders[MVA_label] -> EvaluateMVA(MVA_methods[MVA_label].c_str());
            
            for(float val = -1.; val <= 1.; val += 0.005)
            {
              if( eta_atBTL > -100 && eta_atETL < -100 )
              {
                if( mva_new >= val )
                  ROC_nEvents_EB[MVA_label][label][val] += 1;
              }
            }
          }
          
          outTree -> Fill();
        } // BTL
        
        
        
        // ETL
        if( eta_atETL > -100. && ele_r9 < -100. )
        {
          isEB = 0;
          
          h1_electrons_eta_atETL -> Fill( fabs(electrons_eta_atETL->at(eleIt)) );
          h1_electrons_phi_atETL -> Fill( electrons_phi_atETL->at(eleIt) );
          h1_electrons_mva_ETL -> Fill( electrons_mva->at(eleIt) );
          
          energySeed = 0.;
          float timeSeed = 0.;
          int ietaSeed = 0;
          int iphiSeed = 0;
          for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(eleIt)).size(); ++recHitIt)
          {
            if( (matchedRecHits_det->at(eleIt)).at(recHitIt) != 2 ) continue;
            
            float DR = (matchedRecHits_electron_DR->at(eleIt)).at(recHitIt);
            if( DR > 0.02 ) continue;
            
            float recHitE = (matchedRecHits_energy->at(eleIt)).at(recHitIt);
            float recHitTime = (matchedRecHits_time->at(eleIt)).at(recHitIt);
            
            if( recHitE > energySeed )
            {
              energySeed = recHitE;
              timeSeed = recHitTime;
              ietaSeed = (matchedRecHits_ieta->at(eleIt)).at(recHitIt);
              iphiSeed = (matchedRecHits_iphi->at(eleIt)).at(recHitIt);
            }
          }
          
          nRecHits = 0.;
          energySum = 0.;
          energySumCorr = 0.;
          float ietaAvg = 0;
          float iphiAvg = 0;
          for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(eleIt)).size(); ++recHitIt)
          {
            if( (matchedRecHits_det->at(eleIt)).at(recHitIt) != 2 ) continue;
            
            float recHitE = (matchedRecHits_energy->at(eleIt)).at(recHitIt);
            h1_matchedRecHit_energy_ETL -> Fill( recHitE );
            
            if( recHitE < 0. ) continue;

            float recHitTime = (matchedRecHits_time->at(eleIt)).at(recHitIt);
            h1_matchedRecHit_time_ETL -> Fill( recHitTime-timeSeed );
            p1_matchedRecHit_time_vs_eta -> Fill( fabs(eta),recHitTime );
            
            if( fabs(recHitTime-timeSeed) > 0.1 ) continue;
            
            float DR = (matchedRecHits_electron_DR->at(eleIt)).at(recHitIt);
            h1_matchedRecHit_electron_DR_ETL -> Fill( DR );
            if( DR > 0.02 ) continue;

            nRecHits += 1.;
            
            energySum += recHitE;
            
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
            float recHitTime = (matchedRecHits_time->at(eleIt)).at(recHitIt);
            if( recHitE < 0. ) continue;
            if( fabs(recHitTime-timeSeed) > 0.1 ) continue;
            
            float DR = (matchedRecHits_electron_DR->at(eleIt)).at(recHitIt);
            if( DR > 0.02 ) continue;
            
            sieie += recHitE * pow(ieta-ietaAvg,2);
            sipip += recHitE * pow(iphi-iphiAvg,2);
            
            p2_matchedRecHit_energy_vs_ieta_iphi_ETL -> Fill(iphi-iphiSeed,ieta-ietaSeed,recHitE/energySum);
          }
          sieie = energySum > 0. ? sqrt(sieie/energySum) : -1.;
          sipip = energySum > 0. ? sqrt(sipip/energySum) : -1.;
          
          h1_matchedRecHit_sieie_ETL -> Fill(sieie);
          h1_matchedRecHit_sipip_ETL -> Fill(sipip);
          
          
          // evaluate MVA
          for(float val = -1.; val <= 1.; val += 0.005)
          {
            if( eta_atETL > -100 && eta_atBTL < -100 )
            {
              if( electrons_mva->at(eleIt) >= val )
                ROC_nEvents_EE["old"][label][val] += 1;
            }
          }
          
          for( auto MVA_label : MVALabels_EE )
          {
            float mva_new = MVAReaders[MVA_label] -> EvaluateMVA(MVA_methods[MVA_label].c_str());
            
            for(float val = -1.; val <= 1.; val += 0.005)
            {
              if( eta_atETL > -100 && eta_atBTL < -100 )
              {
                if( mva_new >= val )
                  ROC_nEvents_EE[MVA_label][label][val] += 1;
              }
            }
          }
          
          outTree -> Fill();
        } // ETL
        
        
        
        ++eventId;
        
      } // end loop over tracks
      
      h1_electrons_n -> Fill( nGoodTracks );
      
    } // end loop over events
    
  } // end loop over samples
  
  
  
  //--- make ROC plots
  TGraph* g_ROC_EB = new TGraph();
  
  std::map<float,int> ROC_nEvents_EB_prompt = ROC_nEvents_EB["old"]["prompt"];
  std::map<float,int> ROC_nEvents_EB_fake  = ROC_nEvents_EB["old"]["fake"];
  for(std::map<float,int>::const_iterator mapIt = ROC_nEvents_EB_prompt.begin(); mapIt != ROC_nEvents_EB_prompt.end(); ++mapIt)
  {
    float val = mapIt->first;
    g_ROC_EB -> SetPoint(g_ROC_EB->GetN(),1.*ROC_nEvents_EB_prompt[val]/ROC_nEvents_EB_prompt[-1.],1.*ROC_nEvents_EB_fake[val]/ROC_nEvents_EB_fake[-1.]);
  }

  g_ROC_EB -> Write("g_ROC_EB");

  
  for( auto MVA_label : MVALabels_EB )
  {
    TGraph* g_ROC_new_EB = new TGraph();
    
    std::map<float,int> ROC_new_nEvents_EB_prompt = ROC_nEvents_EB[MVA_label]["prompt"];
    std::map<float,int> ROC_new_nEvents_EB_fake  = ROC_nEvents_EB[MVA_label]["fake"];
    for(std::map<float,int>::const_iterator mapIt = ROC_new_nEvents_EB_prompt.begin(); mapIt != ROC_new_nEvents_EB_prompt.end(); ++mapIt)
    {
      float val = mapIt->first;
      g_ROC_new_EB -> SetPoint(g_ROC_new_EB->GetN(),1.*ROC_new_nEvents_EB_prompt[val]/ROC_new_nEvents_EB_prompt[-1.],1.*ROC_new_nEvents_EB_fake[val]/ROC_new_nEvents_EB_fake[-1.]);
    }
    
    g_ROC_new_EB -> Write(("g_ROC_"+MVA_label).c_str());
  }
  
  
  
  TGraph* g_ROC_EE = new TGraph();
  
  std::map<float,int> ROC_nEvents_EE_prompt = ROC_nEvents_EE["old"]["prompt"];
  std::map<float,int> ROC_nEvents_EE_fake  = ROC_nEvents_EE["old"]["fake"];
  for(std::map<float,int>::const_iterator mapIt = ROC_nEvents_EE_prompt.begin(); mapIt != ROC_nEvents_EE_prompt.end(); ++mapIt)
  {
    float val = mapIt->first;
    g_ROC_EE -> SetPoint(g_ROC_EE->GetN(),1.*ROC_nEvents_EE_prompt[val]/ROC_nEvents_EE_prompt[-1.],1.*ROC_nEvents_EE_fake[val]/ROC_nEvents_EE_fake[-1.]);
  }

  g_ROC_EE -> Write("g_ROC_EE");

  
  for( auto MVA_label : MVALabels_EE )
  {
    TGraph* g_ROC_new_EE = new TGraph();
    
    std::map<float,int> ROC_new_nEvents_EE_prompt = ROC_nEvents_EE[MVA_label]["prompt"];
    std::map<float,int> ROC_new_nEvents_EE_fake  = ROC_nEvents_EE[MVA_label]["fake"];
    for(std::map<float,int>::const_iterator mapIt = ROC_new_nEvents_EE_prompt.begin(); mapIt != ROC_new_nEvents_EE_prompt.end(); ++mapIt)
    {
      float val = mapIt->first;
      g_ROC_new_EE -> SetPoint(g_ROC_new_EE->GetN(),1.*ROC_new_nEvents_EE_prompt[val]/ROC_new_nEvents_EE_prompt[-1.],1.*ROC_new_nEvents_EE_fake[val]/ROC_new_nEvents_EE_fake[-1.]);
    }
    
    g_ROC_new_EE -> Write(("g_ROC_"+MVA_label).c_str());
  }
  
  
  

  
  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
  
}
