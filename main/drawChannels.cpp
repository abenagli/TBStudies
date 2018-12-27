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
    std::cout << ">>> drawChannels.exe::usage:   ./bin/drawChannels.exe   crysLayout [1=tile,2=barphi,3=barz,4=barzflat]" << std::endl;
    return -1;
  }
  int crysLayout = atoi(argv[1]);
  
  
  TFile* inFile;
  if( crysLayout == 1 ) inFile = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/data/tile_singleMuPt5.root",  "READ");
  if( crysLayout == 2 ) inFile = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/data/barphi_singleMuPt5.root","READ");
  if( crysLayout == 3 ) inFile = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/data/barz_singleMuPt5.root",  "READ");
  if( crysLayout == 4 ) inFile = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/data/barzflat_singleMuPtFlat.root",  "READ");
  
  TTree* tree = (TTree*)( inFile->Get("FTLDumpHits/hits_tree") );
  
  std::vector<float>* recHits_time = new std::vector<float>;
  std::vector<float>* recHits_energy = new std::vector<float>;
  std::vector<int>* recHits_module = new std::vector<int>;
  std::vector<int>* recHits_modType = new std::vector<int>;
  std::vector<int>* recHits_crystal = new std::vector<int>;
  std::vector<int>* recHits_ieta = new std::vector<int>;
  std::vector<int>* recHits_iphi = new std::vector<int>;
  
  tree -> SetBranchAddress("recHits_time",    &recHits_time);
  tree -> SetBranchAddress("recHits_energy",  &recHits_energy);
  tree -> SetBranchAddress("recHits_module",  &recHits_module);
  tree -> SetBranchAddress("recHits_modType", &recHits_modType);
  tree -> SetBranchAddress("recHits_crystal", &recHits_crystal);
  tree -> SetBranchAddress("recHits_ieta",    &recHits_ieta);
  tree -> SetBranchAddress("recHits_iphi",    &recHits_iphi);
  
  
  //--- open output file
  TFile* outFile;
  if( crysLayout == 1 ) outFile = TFile::Open(Form("/Users/abenagli/Work/TIMING/TBStudies/plots/tile_channels.root"),  "RECREATE");
  if( crysLayout == 2 ) outFile = TFile::Open(Form("/Users/abenagli/Work/TIMING/TBStudies/plots/barphi_channels.root"),"RECREATE");
  if( crysLayout == 3 ) outFile = TFile::Open(Form("/Users/abenagli/Work/TIMING/TBStudies/plots/barz_channels.root"),  "RECREATE");
  if( crysLayout == 4 ) outFile = TFile::Open(Form("/Users/abenagli/Work/TIMING/TBStudies/plots/barzflat_channels.root"),  "RECREATE");
  outFile -> cd();
  
  TProfile2D* p2_crystal = new TProfile2D("p2_crystal","",865,-0.5,864.5,577,-0.5,576.5);
  TProfile2D* p2_crystal_foldIEta = new TProfile2D("p2_crystal_foldIEta","",16,-0.5,15.5,577,-0.5,576.5);
  TProfile2D* p2_crystal_foldIEtaIPhi = new TProfile2D("p2_crystal_foldIEtaIPhi","",16,-0.5,15.5,64,-0.5,63.5);
  
  TProfile2D* p2_modType = new TProfile2D("p2_modType","",865,-0.5,864.5,577,-0.5,576.5);
  TProfile2D* p2_modType_foldIPhi = new TProfile2D("p2_modType_foldIPhi","",865,-0.5,864.5,1,-0.5,0.5);
  TProfile2D* p2_modType_foldIEtaIPhi = new TProfile2D("p2_modType_foldIEtaIPhi","",16,-0.5,15.5,1,-0.5,0.5);
  
  TProfile2D* p2_module = new TProfile2D("p2_module","",865,-0.5,864.5,577,-0.5,576.5);
  TProfile2D* p2_module_foldIPhi = new TProfile2D("p2_module_foldIPhi","",865,-0.5,864.5,1,-0.5,0.5);
  
  
  
  //--- loop over events
  int nEntries = tree->GetEntries();
  for(int entry = 0; entry < nEntries; ++entry)
  {
    tree -> GetEntry(entry);

    for(unsigned int jj = 0; jj < recHits_energy->size(); ++jj)
    {
      int ieta = recHits_ieta->at(jj);
      int iphi = recHits_iphi->at(jj);
      int crystal = recHits_crystal->at(jj);
      int modType = recHits_modType->at(jj);
      int module = recHits_module->at(jj);
      
      if( ieta < 0 ) continue;

      if( crysLayout == 1 )
      {
        p2_crystal -> Fill( ieta,iphi,crystal );
        p2_crystal_foldIEta -> Fill( int(ieta-1)%4,iphi,crystal );
        p2_crystal_foldIEtaIPhi -> Fill( int(ieta-1)%4,int(iphi-1)%16,crystal );
        
        p2_modType -> Fill( ieta,iphi,modType );
        p2_modType_foldIPhi -> Fill( ieta,int(iphi-1)%1,modType );
        p2_modType_foldIEtaIPhi -> Fill( int(ieta-1)/(18*4),int(iphi-1)%1,modType );
        
        p2_module -> Fill( ieta,iphi,module );
        p2_module_foldIPhi -> Fill( ieta,int(iphi-1)%1,module );
      }
      
      if( crysLayout == 2 )
      {
        p2_crystal -> Fill( ieta,iphi,crystal );
        p2_crystal_foldIEta -> Fill( int(ieta-1)%16,iphi,crystal );
        p2_crystal_foldIEtaIPhi -> Fill( int(ieta-1)%16,int(iphi-1)%4,crystal );
        
        p2_modType -> Fill( ieta,iphi,modType );
        p2_modType_foldIPhi -> Fill( ieta,int(iphi-1)%1,modType );
        p2_modType_foldIEtaIPhi -> Fill( int(ieta-1)/(18*16),int(iphi-1)%1,modType );
        
        p2_module -> Fill( ieta,iphi,module );
        p2_module_foldIPhi -> Fill( ieta,int(iphi-1)%1,module );        
      }
      
      if( crysLayout == 3 )
      {
        p2_crystal -> Fill( ieta,iphi,crystal );
        p2_crystal_foldIEta -> Fill( int(ieta-1)%1,iphi,crystal );
        p2_crystal_foldIEtaIPhi -> Fill( int(ieta-1)%1,int(iphi-1)%64,crystal );
        
        p2_modType -> Fill( ieta,iphi,modType );
        p2_modType_foldIPhi -> Fill( ieta,int(iphi-1)%1,modType );
        p2_modType_foldIEtaIPhi -> Fill( int(ieta-1)/(18),int(iphi-1)%1,modType );
        
        p2_module -> Fill( ieta,iphi,module );
        p2_module_foldIPhi -> Fill( ieta,int(iphi-1)%1,module );        
      }
      
      if( crysLayout == 4 )
      {
        p2_crystal -> Fill( ieta,iphi,crystal );
        p2_crystal_foldIEta -> Fill( int(ieta-1)%1,iphi,crystal );
        p2_crystal_foldIEtaIPhi -> Fill( int(ieta-1)%1,int(iphi-1)%64,crystal );
        
        p2_modType -> Fill( ieta,iphi,modType );
        p2_modType_foldIPhi -> Fill( ieta,int(iphi-1)%1,modType );
        p2_modType_foldIEtaIPhi -> Fill( int(ieta-1)/(14),int(iphi-1)%1,modType );
        
        p2_module -> Fill( ieta,iphi,module );
        p2_module_foldIPhi -> Fill( ieta,int(iphi-1)%1,module );        
      }
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
