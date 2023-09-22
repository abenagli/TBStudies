void drawOscHisto()
{
  TCanvas* c1 = new TCanvas();
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,256,0.04) );
  hPad -> SetTitle(Form(";bin;counts"));
  hPad -> Draw();
  
  //--- legend
  TLegend* legend = new TLegend(0.40,0.70,0.75,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  
  std::vector<std::string> inFiles;
  inFiles.push_back("data/Na22_CPI_HPK50um_noGrease000.csv");
  inFiles.push_back("data/Na22_CPI_HPK50um_grease000.csv");
  inFiles.push_back("data/Na22_Hilger_HPK50um_noGrease000.csv");
  inFiles.push_back("data/Na22_Hilger_HPK50um_grease000.csv");
  
  std::vector<int> colors;
  colors.push_back(kRed);
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kBlue);
  
  std::vector<int> linestyles;
  linestyles.push_back(3);
  linestyles.push_back(1);
  linestyles.push_back(3);
  linestyles.push_back(1);
  
  std::vector<std::string> labels;
  labels.push_back("CPI 11.5#times11.5 mm^{2}- no grease");
  labels.push_back("CPI 11.5#times11.5 mm^{2}- grease");
  labels.push_back("Hilger 11#times11 mm^{2}- no grease");
  labels.push_back("Hilger 11#times11 mm^{2}- grease");
  
  int inFileIt = 0;
  for(auto inFile: inFiles)
  {
    TH1F* histo = new TH1F(Form("histo_%d",inFileIt),"",251,0,251);
    histo -> Sumw2();
    
    std::ifstream list(inFile.c_str(),std::ios::in);
    std::string line;
    
    int iBin = 1;
    while(1)
    {
      getline(list,line,'\n');
      if( !list.good() ) break;
      if( line.at(0) == '#' ) continue;
      
      std::stringstream ss(line);
      std::string token;
      int ii = 0;
      while( std::getline(ss,token,',') )
      {
        if( ii%2 == 1 )
          histo -> SetBinContent(iBin,atof(token.c_str()));
        histo -> SetBinError(iBin,sqrt(atof(token.c_str())));
        ++ii;
      }
      
      ++iBin;
    }
    
    histo -> SetLineWidth(2);
    histo -> SetLineStyle(linestyles.at(inFileIt));
    histo -> SetLineColor(colors.at(inFileIt));
    histo -> Scale(1./histo->Integral());
    histo -> Draw("hist,same");
    
    legend -> AddEntry(histo,labels.at(inFileIt).c_str(),"L");
    
    ++inFileIt;
  }
  
  legend -> Draw("same");
}


