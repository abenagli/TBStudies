#pragma link C++ class vector<vector<vector<int> > >+;
#pragma link C++ class vector<vector<vector<float> > >+;



void drawMappingPlots(const std::string& inFileName)
{
  std::string label = "default";
  TLatex* latexLabel = new TLatex(0.16,0.96,Form("bars"));
  latexLabel -> SetNDC();
  latexLabel -> SetTextFont(42);
  latexLabel -> SetTextSize(0.03);
  
  

  TChain* tree = new TChain("FTLDumpHits/DumpHits","FTLDumpHits/DumpHits");
  TFile* inFile = TFile::Open(inFileName.c_str(),"READ");
  TTree* tree = (TTree*)( inFile->Get("FTLDumpHits/DumpHits"));
  int nEntries = tree->GetEntries();
  
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
  tree -> SetBranchStatus("matchedRecHits_det",1);     tree -> SetBranchAddress("matchedRecHits_det",     &matchedRecHits_det);
  tree -> SetBranchStatus("matchedRecHits_energy",1);  tree -> SetBranchAddress("matchedRecHits_energy",  &matchedRecHits_energy);
  tree -> SetBranchStatus("matchedRecHits_runit",1);   tree -> SetBranchAddress("matchedRecHits_runit",   &matchedRecHits_runit);
  tree -> SetBranchStatus("matchedRecHits_modType",1); tree -> SetBranchAddress("matchedRecHits_modType", &matchedRecHits_modType);
  

  
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
      
      
      std::vector<float> recHitEs;
      float totEnergy = 0;
      for(unsigned int recHitIt = 0; recHitIt < (matchedRecHits_energy->at(trackIt)).size(); ++recHitIt)
      {
        if( (matchedRecHits_det->at(trackIt)).at(recHitIt) != 1 ) continue;
        float recHitE = (matchedRecHits_energy->at(trackIt)).at(recHitIt);
        if( recHitE < 1. ) continue;
        
        recHitEs.push_back( recHitE );
        totEnergy += recHitE;
      }
      if( totEnergy < 1. ) continue;
    }
  }
  
  

  /*
  TCanvas* c1 = new TCanvas("c1","c1",1400,500);
  c1 -> Divide(3,1);

  c1 -> cd(1);
  gPad -> SetLogy();
  TH1F* h_entry_x = new TH1F("h_entry_x","",30000,-3.,3.);
  t -> Draw("simHits_entry_local_x>>h_entry_x","simHits_det == 1","goff");
  h_entry_x -> SetTitle(";simHit entry/exit point x (cm);entries");
  h_entry_x -> Draw();
  TH1F* h_exit_x = new TH1F("h_exit_x","simHits_det == 1",30000,-3.,3.);
  t -> Draw("simHits_exit_local_x>>h_exit_x","","goff");
  h_exit_x -> SetLineColor(kRed);
  h_exit_x -> Draw("same");
  latexLabel -> Draw("same");
  */
}
