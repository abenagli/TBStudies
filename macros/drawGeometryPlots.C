template <class T>
void drawTH1F(TCanvas* c1,
              const std::string& histoName, const std::string& extra, std::map<std::string,TFile*>& inFiles,
              std::vector<std::string>& inputs1, std::vector<std::string>& inputs2, std::map<std::string,int>& colors, std::map<std::string,std::string>& labels, 
              const float& xMin, const float& yMin, const float& xMax, const float& yMax,
              const bool& normalize,
              const std::string& title)
{
  c1 = new TCanvas(Form("c1_%s%s",histoName.c_str(),extra.c_str()),Form("c1_%s%s",histoName.c_str(),extra.c_str()),1400,500);
  c1 -> Divide(2,1);

  c1 -> cd(1);
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(xMin,yMin,xMax,yMax) );
  hPad -> SetTitle(title.c_str());
  hPad -> Draw();
  
  TLegend* legend = new TLegend(0.35,0.70,0.75,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);

  for(auto label : inputs1)
  {
    T* histo = (T*)( inFiles[label]->Get(histoName.c_str()) );
    
    histo -> SetLineColor(colors[label]);
    histo -> SetMarkerColor(colors[label]);
    histo -> SetMarkerSize(0.7);
    if( normalize )
    {
      histo -> Scale( 1./histo->Integral() );
      histo -> SetLineWidth(3);
      histo -> SetFillColor(colors[label]);
      histo -> SetFillStyle(3001);
    }
    histo -> Draw("hist,same");
    if( !normalize ) histo -> Draw("PL,same");

    if( normalize )  legend -> AddEntry((T*)(histo),labels[label].c_str(),"F");
    if( !normalize ) legend -> AddEntry((T*)(histo),labels[label].c_str(),"PL");
  }
  
  legend -> Draw("same");
  
  
  c1 -> cd(2);
  
  hPad = (TH1F*)( gPad->DrawFrame(xMin,yMin,xMax,yMax) );
  hPad -> SetTitle(title.c_str());
  hPad -> Draw();
  
  legend = new TLegend(0.35,0.70,0.75,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);

  for(auto label : inputs2)
  {
    T* histo = (T*)( inFiles[label]->Get(histoName.c_str()) );
    
    histo -> SetLineColor(colors[label]);
    histo -> SetMarkerColor(colors[label]);
    histo -> SetMarkerSize(0.7);
    if( normalize )
    {
      histo -> Scale( 1./histo->Integral() );
      histo -> SetLineWidth(3);
      histo -> SetFillColor(colors[label]);
      histo -> SetFillStyle(3001);
    }
    histo -> Draw("hist,same");
    if( !normalize ) histo -> Draw("PL,same");

    if( normalize )  legend -> AddEntry((T*)(histo),labels[label].c_str(),"F");
    if( !normalize ) legend -> AddEntry((T*)(histo),labels[label].c_str(),"PL");
  }
  
  legend -> Draw("same");

  c1 -> Print(Form("plots/geometryPlots/c_%s%s.png",histoName.c_str(),extra.c_str()));
  c1 -> Print(Form("plots/geometryPlots/c_%s%s.pdf",histoName.c_str(),extra.c_str()));
}



void drawGeometryPlots()
{
  std::map<std::string,std::string> inFileNames;
  inFileNames["tile_10GeV"] = "plots/geometryPlots/geometryPlots_tile_singleMuPt10.root";
  inFileNames["tile_05GeV"] = "plots/geometryPlots/geometryPlots_tile_singleMuPt5.root";
  inFileNames["tile_03GeV"] = "plots/geometryPlots/geometryPlots_tile_singleMuPt3.root";
  inFileNames["tile_02GeV"] = "plots/geometryPlots/geometryPlots_tile_singleMuPt2.root";
  inFileNames["tile_01GeV"] = "plots/geometryPlots/geometryPlots_tile_singleMuPt1.root";
  inFileNames["bar_10GeV"] = "plots/geometryPlots/geometryPlots_bar_singleMuPt10.root";
  inFileNames["bar_05GeV"] = "plots/geometryPlots/geometryPlots_bar_singleMuPt5.root";
  inFileNames["bar_03GeV"] = "plots/geometryPlots/geometryPlots_bar_singleMuPt3.root";
  inFileNames["bar_02GeV"] = "plots/geometryPlots/geometryPlots_bar_singleMuPt2.root";
  inFileNames["bar_01GeV"] = "plots/geometryPlots/geometryPlots_bar_singleMuPt1.root";
  
  std::map<std::string,std::string> labels;
  labels["tile_10GeV"] = "tile geo. - 10 GeV p_{T} #mu gun";
  labels["tile_05GeV"] = "tile geo. -  5 GeV p_{T} #mu gun";
  labels["tile_03GeV"] = "tile geo. -  3 GeV p_{T} #mu gun";
  labels["tile_02GeV"] = "tile geo. -  2 GeV p_{T} #mu gun";
  labels["tile_01GeV"] = "tile geo. -  1 GeV p_{T} #mu gun";
  labels["bar_10GeV"] = "bar geo. - 10 GeV p_{T} #mu gun";
  labels["bar_05GeV"] = "bar geo. -  5 GeV p_{T} #mu gun";
  labels["bar_03GeV"] = "bar geo. -  3 GeV p_{T} #mu gun";
  labels["bar_02GeV"] = "bar geo. -  2 GeV p_{T} #mu gun";
  labels["bar_01GeV"] = "bar geo. -  1 GeV p_{T} #mu gun";
  
  std::map<std::string,int> colors;
  colors["tile_10GeV"] = kRed;
  colors["tile_05GeV"] = 97;
  colors["tile_03GeV"] = 94;
  colors["tile_02GeV"] = 94;
  colors["tile_01GeV"] = 91;
  colors["bar_10GeV"] = kBlue;
  colors["bar_05GeV"] = 65;
  colors["bar_03GeV"] = 68;
  colors["bar_02GeV"] = 68;
  colors["bar_01GeV"] = 71;
  
  std::map<std::string,TFile*> inFiles;
  for(auto mapIt : inFileNames)
  {
    inFiles[mapIt.first] = TFile::Open(mapIt.second.c_str(),"READ");
    std::cout << "opening file " << mapIt.second << std::endl;
  }
  
  
  
  TCanvas* c1;
  
  std::vector<std::string> inputs1;
  inputs1.push_back("tile_01GeV");
  inputs1.push_back("tile_03GeV");
  inputs1.push_back("tile_05GeV");
  inputs1.push_back("tile_10GeV");
  
  std::vector<std::string> inputs2;
  inputs2.push_back("bar_01GeV");
  inputs2.push_back("bar_03GeV");
  inputs2.push_back("bar_05GeV");
  inputs2.push_back("bar_10GeV");
  
  drawTH1F<TH1F>(c1,
                  "h1_recHits_n","",inFiles,
                  inputs1,inputs2,colors,labels,
                  0.,0.001,8.,1.,
                  true,
                  ";N_{BTL recHit} [MeV];event fracion");
  
  drawTH1F<TProfile>(c1,
                      "p1_recHits_n_vs_eta","",inFiles,
                      inputs1,inputs2,colors,labels,
                      0.,0.,1.5,5.,
                      false,
                      ";#eta_{#mu};#LT N_{BTL recHit} #GT [MeV]");
  
  drawTH1F<TProfile>(c1,
                      "p1_recHits_n_vs_phi","",inFiles,
                      inputs1,inputs2,colors,labels,
                      -3.15,0.,3.15,5.,
                      false,
                      ";#phi_{#mu};#LT N_{BTL recHit} #GT [MeV]");

  drawTH1F<TProfile>(c1,
                      "p1_recHits_n_vs_phi","_zoom",inFiles,
                      inputs1,inputs2,colors,labels,
                      0.,0.,1.,5.,
                      false,
                      ";#phi_{#mu};#LT N_{BTL recHit} #GT [MeV]");


  
  drawTH1F<TH1F>(c1,
                  "h1_recHits_energy","",inFiles,
                  inputs1,inputs2,colors,labels,
                  0.,0.001,10.,0.15,
                  true,
                  ";E_{BTL recHit} [MeV];event fracion");
  
  drawTH1F<TProfile>(c1,
                      "p1_recHits_energy_vs_eta","",inFiles,
                      inputs1,inputs2,colors,labels,
                      0.,0.,1.5,10.,
                      false,
                      ";#eta_{#mu};#LT E_{BTL recHit} #GT [MeV]");
  
  drawTH1F<TProfile>(c1,
                      "p1_recHits_energy_vs_phi","",inFiles,
                      inputs1,inputs2,colors,labels,
                      -3.15,0.,3.15,10.,
                      false,
                      ";#phi_{#mu};#LT E_{BTL recHit} #GT [MeV]");

  drawTH1F<TProfile>(c1,
                     "p1_recHits_energy_vs_phi","_zoom",inFiles,
                     inputs1,inputs2,colors,labels,
                     0.,0.,1.,15.,
                     false,
                     ";#phi_{#mu};#LT E_{BTL recHit} #GT [MeV]");
  
  
  
  drawTH1F<TH1F>(c1,
                  "h1_recHits_energySum","",inFiles,
                  inputs1,inputs2,colors,labels,
                  0.,0.001,15.,0.15,
                  true,
                  ";#Sigma E_{BTL recHit} [MeV];event fracion");
  
  drawTH1F<TProfile>(c1,
                      "p1_recHits_energySum_vs_eta","",inFiles,
                      inputs1,inputs2,colors,labels,
                      0.,0.,1.5,15.,
                      false,
                      ";#eta_{#mu};#LT #Sigma E_{BTL recHit} #GT [MeV]");
  
  drawTH1F<TProfile>(c1,
                      "p1_recHits_energySum_vs_phi","",inFiles,
                      inputs1,inputs2,colors,labels,
                      -3.15,0.,3.15,15.,
                      false,
                      ";#phi_{#mu};#LT #Sigma E_{BTL recHit} #GT [MeV]");

  drawTH1F<TProfile>(c1,
                     "p1_recHits_energySum_vs_phi","_zoom",inFiles,
                     inputs1,inputs2,colors,labels,
                     0.,0.,1.,15.,
                     false,
                     ";#phi_{#mu};#LT #Sigma E_{BTL recHit} #GT [MeV]");
  

  
  drawTH1F<TH1F>(c1,
                 "p1_eff_vs_eta","",inFiles,
                 inputs1,inputs2,colors,labels,
                 0.,0.,1.5,1.5,
                 false,
                 ";#eta_{#mu};efficiency");

  drawTH1F<TH1F>(c1,
                 "p1_eff_vs_phi","",inFiles,
                 inputs1,inputs2,colors,labels,
                 -3.15,0.,3.15,1.5,
                 false,
                 ";#phi_{#mu};efficiency");
  
  drawTH1F<TH1F>(c1,
                 "p1_eff_vs_phi","_zoom",inFiles,
                 inputs1,inputs2,colors,labels,
                 0.,0.,1.,1.5,
                 false,
                 ";#phi_{#mu};efficiency");
  
  
  
  drawTH1F<TH1F>(c1,
                 "p1_recHits_time_vs_eta","",inFiles,
                 inputs1,inputs2,colors,labels,
                 0.,0.,1.5,10.,
                 false,
                 ";#eta_{#mu};time");
}
