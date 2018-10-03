


void drawTrends(const std::string& fileName)
{
  gStyle->SetPadRightMargin(0.50);

  TFile* inFile = TFile::Open(fileName.c_str(),"READ");
  TTree* tree = (TTree*)( inFile->Get("results") );
  
  std::vector<std::string> chs;
  std::map<std::string,std::string> labels;
  chs.push_back("CH1"); labels["CH1"] = "FBK 4x4 mm^{2} 25um - LYSO 11.5x11.5x3.75 mm^{3}";
  chs.push_back("CH2"); labels["CH2"] = "FBK 5x5 mm^{2} 20um - LYSO 11.5x11.5x3.75 mm^{3}";

  for(auto ch : chs)
  {
    TLatex* latexLabel = new TLatex(0.16,0.96,Form("%s",labels[ch].c_str()));
    latexLabel -> SetNDC();
    latexLabel -> SetTextFont(42);
    latexLabel -> SetTextSize(0.03);
    
    float t_Vbias;
    float t_NINOthr;
    float t_CTR_effSigma;
    float t_CTR_gausSigma;
    float t_CTR_gausSigmaErr;
    float t_CTR_ampCorr_effSigma;
    float t_CTR_ampCorr_gausSigma;
    float t_CTR_ampCorr_gausSigmaErr;
    float t_CTR_ampCorr_posCorr_effSigma;
    float t_CTR_ampCorr_posCorr_gausSigma;
    float t_CTR_ampCorr_posCorr_gausSigmaErr;
    float t_CTR_RTCorr_effSigma;
    float t_CTR_RTCorr_gausSigma;
    float t_CTR_RTCorr_gausSigmaErr;
    float t_CTR_RTCorr_posCorr_effSigma;
    float t_CTR_RTCorr_posCorr_gausSigma;
    float t_CTR_RTCorr_posCorr_gausSigmaErr;
    float t_CTR_distAmpCorr_effSigma;
    float t_CTR_distAmpCorr_gausSigma;
    float t_CTR_distAmpCorr_gausSigmaErr;
    float t_CTR_distRTCorr_effSigma;
    float t_CTR_distRTCorr_gausSigma;
    float t_CTR_distRTCorr_gausSigmaErr;

    tree -> SetBranchAddress(Form("Vbias_%s",ch.c_str()),&t_Vbias);
    tree -> SetBranchAddress(Form("NINOthr_%s",ch.c_str()),&t_NINOthr);
    tree -> SetBranchAddress(Form("CTR_effSigma_%s",ch.c_str()),&t_CTR_effSigma);
    tree -> SetBranchAddress(Form("CTR_gausSigma_%s",ch.c_str()),&t_CTR_gausSigma);
    tree -> SetBranchAddress(Form("CTR_gausSigmaErr_%s",ch.c_str()),&t_CTR_gausSigmaErr);
    tree -> SetBranchAddress(Form("CTR_ampCorr_effSigma_%s",ch.c_str()),&t_CTR_ampCorr_effSigma);
    tree -> SetBranchAddress(Form("CTR_ampCorr_gausSigma_%s",ch.c_str()),&t_CTR_ampCorr_gausSigma);
    tree -> SetBranchAddress(Form("CTR_ampCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_ampCorr_gausSigmaErr);
    tree -> SetBranchAddress(Form("CTR_ampCorr_posCorr_effSigma_%s",ch.c_str()),&t_CTR_ampCorr_posCorr_effSigma);
    tree -> SetBranchAddress(Form("CTR_ampCorr_posCorr_gausSigma_%s",ch.c_str()),&t_CTR_ampCorr_posCorr_gausSigma);
    tree -> SetBranchAddress(Form("CTR_ampCorr_posCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_ampCorr_posCorr_gausSigmaErr);
    tree -> SetBranchAddress(Form("CTR_RTCorr_effSigma_%s",ch.c_str()),&t_CTR_RTCorr_effSigma);
    tree -> SetBranchAddress(Form("CTR_RTCorr_gausSigma_%s",ch.c_str()),&t_CTR_RTCorr_gausSigma);
    tree -> SetBranchAddress(Form("CTR_RTCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_RTCorr_gausSigmaErr);
    tree -> SetBranchAddress(Form("CTR_RTCorr_posCorr_effSigma_%s",ch.c_str()),&t_CTR_RTCorr_posCorr_effSigma);
    tree -> SetBranchAddress(Form("CTR_RTCorr_posCorr_gausSigma_%s",ch.c_str()),&t_CTR_RTCorr_posCorr_gausSigma);
    tree -> SetBranchAddress(Form("CTR_RTCorr_posCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_RTCorr_posCorr_gausSigmaErr);
    tree -> SetBranchAddress(Form("CTR_distAmpCorr_effSigma_%s",ch.c_str()),&t_CTR_distAmpCorr_effSigma);
    tree -> SetBranchAddress(Form("CTR_distAmpCorr_gausSigma_%s",ch.c_str()),&t_CTR_distAmpCorr_gausSigma);
    tree -> SetBranchAddress(Form("CTR_distAmpCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_distAmpCorr_gausSigmaErr);
    tree -> SetBranchAddress(Form("CTR_distRTCorr_effSigma_%s",ch.c_str()),&t_CTR_distRTCorr_effSigma);
    tree -> SetBranchAddress(Form("CTR_distRTCorr_gausSigma_%s",ch.c_str()),&t_CTR_distRTCorr_gausSigma);
    tree -> SetBranchAddress(Form("CTR_distRTCorr_gausSigmaErr_%s",ch.c_str()),&t_CTR_distRTCorr_gausSigmaErr);
    
    std::map<float,std::vector<float> > vec_Vbias;
    std::map<float,std::vector<float> > vec_NINOthr;
    
    std::map<float,TGraphErrors*> g_vs_NINOthr_CTR_effSigma;
    std::map<float,TGraphErrors*> g_vs_NINOthr_CTR_ampCorr_effSigma;
    std::map<float,TGraphErrors*> g_vs_NINOthr_CTR_distAmpCorr_effSigma;
    
    std::map<float,TGraphErrors*> g_vs_NINOthr_CTR_gausSigma;
    std::map<float,TGraphErrors*> g_vs_NINOthr_CTR_ampCorr_gausSigma;
    std::map<float,TGraphErrors*> g_vs_NINOthr_CTR_distAmpCorr_gausSigma;
    
    std::map<float,TGraphErrors*> g_vs_Vbias_CTR_effSigma;
    std::map<float,TGraphErrors*> g_vs_Vbias_CTR_ampCorr_effSigma;
    std::map<float,TGraphErrors*> g_vs_Vbias_CTR_distAmpCorr_effSigma;
    
    std::map<float,TGraphErrors*> g_vs_Vbias_CTR_gausSigma;
    std::map<float,TGraphErrors*> g_vs_Vbias_CTR_ampCorr_gausSigma;
    std::map<float,TGraphErrors*> g_vs_Vbias_CTR_distAmpCorr_gausSigma;
    
    for(int entry = 0; entry < tree->GetEntries(); ++entry)
    {
      tree -> GetEntry(entry);
      
      vec_Vbias[t_NINOthr].push_back(t_Vbias);
      vec_NINOthr[t_Vbias].push_back(t_NINOthr);

      
      if( g_vs_NINOthr_CTR_effSigma[t_Vbias] == NULL ) g_vs_NINOthr_CTR_effSigma[t_Vbias] = new TGraphErrors();
      g_vs_NINOthr_CTR_effSigma[t_Vbias] -> SetPoint(g_vs_NINOthr_CTR_effSigma[t_Vbias]->GetN(),t_NINOthr,t_CTR_effSigma);
      std::cout << "Vbias: " << t_Vbias << "   NINOthr: " << t_NINOthr << std::endl;
      

      if( g_vs_NINOthr_CTR_ampCorr_effSigma[t_Vbias] == NULL ) g_vs_NINOthr_CTR_ampCorr_effSigma[t_Vbias] = new TGraphErrors();
      g_vs_NINOthr_CTR_ampCorr_effSigma[t_Vbias] -> SetPoint(g_vs_NINOthr_CTR_ampCorr_effSigma[t_Vbias]->GetN(),t_NINOthr,t_CTR_ampCorr_effSigma);
      
      if( g_vs_NINOthr_CTR_distAmpCorr_effSigma[t_Vbias] == NULL ) g_vs_NINOthr_CTR_distAmpCorr_effSigma[t_Vbias] = new TGraphErrors();
      g_vs_NINOthr_CTR_distAmpCorr_effSigma[t_Vbias] -> SetPoint(g_vs_NINOthr_CTR_distAmpCorr_effSigma[t_Vbias]->GetN(),t_NINOthr,t_CTR_distAmpCorr_effSigma);
      
      if( g_vs_NINOthr_CTR_gausSigma[t_Vbias] == NULL ) g_vs_NINOthr_CTR_gausSigma[t_Vbias] = new TGraphErrors();
      g_vs_NINOthr_CTR_gausSigma[t_Vbias] -> SetPoint(g_vs_NINOthr_CTR_gausSigma[t_Vbias]->GetN(),t_NINOthr,t_CTR_gausSigma);
      g_vs_NINOthr_CTR_gausSigma[t_Vbias] -> SetPointError(g_vs_NINOthr_CTR_gausSigma[t_Vbias]->GetN()-1,0.,t_CTR_gausSigmaErr);
      
      if( g_vs_NINOthr_CTR_ampCorr_gausSigma[t_Vbias] == NULL ) g_vs_NINOthr_CTR_ampCorr_gausSigma[t_Vbias] = new TGraphErrors();
      g_vs_NINOthr_CTR_ampCorr_gausSigma[t_Vbias] -> SetPoint(g_vs_NINOthr_CTR_ampCorr_gausSigma[t_Vbias]->GetN(),t_NINOthr,t_CTR_ampCorr_gausSigma);
      g_vs_NINOthr_CTR_ampCorr_gausSigma[t_Vbias] -> SetPointError(g_vs_NINOthr_CTR_ampCorr_gausSigma[t_Vbias]->GetN()-1,0.,t_CTR_ampCorr_gausSigmaErr);
      
      if( g_vs_NINOthr_CTR_distAmpCorr_gausSigma[t_Vbias] == NULL ) g_vs_NINOthr_CTR_distAmpCorr_gausSigma[t_Vbias] = new TGraphErrors();
      g_vs_NINOthr_CTR_distAmpCorr_gausSigma[t_Vbias] -> SetPoint(g_vs_NINOthr_CTR_distAmpCorr_gausSigma[t_Vbias]->GetN(),t_NINOthr,t_CTR_distAmpCorr_gausSigma);
      g_vs_NINOthr_CTR_distAmpCorr_gausSigma[t_Vbias] -> SetPointError(g_vs_NINOthr_CTR_distAmpCorr_gausSigma[t_Vbias]->GetN()-1,0.,t_CTR_distAmpCorr_gausSigmaErr);
      
      
      if( g_vs_Vbias_CTR_effSigma[t_NINOthr] == NULL ) g_vs_Vbias_CTR_effSigma[t_NINOthr] = new TGraphErrors();
      g_vs_Vbias_CTR_effSigma[t_NINOthr] -> SetPoint(g_vs_Vbias_CTR_effSigma[t_NINOthr]->GetN(),t_Vbias,t_CTR_effSigma);
      
      if( g_vs_Vbias_CTR_ampCorr_effSigma[t_NINOthr] == NULL ) g_vs_Vbias_CTR_ampCorr_effSigma[t_NINOthr] = new TGraphErrors();
      g_vs_Vbias_CTR_ampCorr_effSigma[t_NINOthr] -> SetPoint(g_vs_Vbias_CTR_ampCorr_effSigma[t_NINOthr]->GetN(),t_Vbias,t_CTR_ampCorr_effSigma);
      
      if( g_vs_Vbias_CTR_distAmpCorr_effSigma[t_NINOthr] == NULL ) g_vs_Vbias_CTR_distAmpCorr_effSigma[t_NINOthr] = new TGraphErrors();
      g_vs_Vbias_CTR_distAmpCorr_effSigma[t_NINOthr] -> SetPoint(g_vs_Vbias_CTR_distAmpCorr_effSigma[t_NINOthr]->GetN(),t_Vbias,t_CTR_distAmpCorr_effSigma);
      
      if( g_vs_Vbias_CTR_gausSigma[t_NINOthr] == NULL ) g_vs_Vbias_CTR_gausSigma[t_NINOthr] = new TGraphErrors();
      g_vs_Vbias_CTR_gausSigma[t_NINOthr] -> SetPoint(g_vs_Vbias_CTR_gausSigma[t_NINOthr]->GetN(),t_Vbias,t_CTR_gausSigma);
      g_vs_Vbias_CTR_gausSigma[t_NINOthr] -> SetPointError(g_vs_Vbias_CTR_gausSigma[t_NINOthr]->GetN()-1,0.,t_CTR_gausSigmaErr);
      
      if( g_vs_Vbias_CTR_ampCorr_gausSigma[t_NINOthr] == NULL ) g_vs_Vbias_CTR_ampCorr_gausSigma[t_NINOthr] = new TGraphErrors();
      g_vs_Vbias_CTR_ampCorr_gausSigma[t_NINOthr] -> SetPoint(g_vs_Vbias_CTR_ampCorr_gausSigma[t_NINOthr]->GetN(),t_Vbias,t_CTR_ampCorr_gausSigma);
      g_vs_Vbias_CTR_ampCorr_gausSigma[t_NINOthr] -> SetPointError(g_vs_Vbias_CTR_ampCorr_gausSigma[t_NINOthr]->GetN()-1,0.,t_CTR_ampCorr_gausSigmaErr);
      
      if( g_vs_Vbias_CTR_distAmpCorr_gausSigma[t_NINOthr] == NULL ) g_vs_Vbias_CTR_distAmpCorr_gausSigma[t_NINOthr] = new TGraphErrors();
      g_vs_Vbias_CTR_distAmpCorr_gausSigma[t_NINOthr] -> SetPoint(g_vs_Vbias_CTR_distAmpCorr_gausSigma[t_NINOthr]->GetN(),t_Vbias,t_CTR_distAmpCorr_gausSigma);
      g_vs_Vbias_CTR_distAmpCorr_gausSigma[t_NINOthr] -> SetPointError(g_vs_Vbias_CTR_distAmpCorr_gausSigma[t_NINOthr]->GetN()-1,0.,t_CTR_distAmpCorr_gausSigmaErr);
    }
    

    
    for(auto it : vec_NINOthr)
    {
      float Vbias = it.first;
      std::sort(vec_NINOthr[Vbias].begin(),vec_NINOthr[Vbias].end());
      
      TCanvas* c_vs_NINOthr = new TCanvas(Form("c_Vbias%.0f_vs_NINOthr_%s",Vbias,ch.c_str()),Form("c_Vbias%.0f_vs_NINOthr_%s",Vbias,ch.c_str()),1000,500);
      c_vs_NINOthr -> cd();
      
      TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,10.,1.5*vec_NINOthr[Vbias].at(vec_NINOthr[Vbias].size()-1),200.) );
      hPad -> SetTitle(";NINO thr. [mV];time resolution [ps]");
      hPad -> Draw();

      
      g_vs_NINOthr_CTR_effSigma[Vbias] -> SetMarkerStyle(20);
      g_vs_NINOthr_CTR_effSigma[Vbias] -> SetMarkerColor(kBlack);
      g_vs_NINOthr_CTR_effSigma[Vbias] -> SetLineColor(kBlack);
      g_vs_NINOthr_CTR_effSigma[Vbias] -> SetLineStyle(1);
      g_vs_NINOthr_CTR_effSigma[Vbias] -> Draw("P,same");
      
      g_vs_NINOthr_CTR_ampCorr_effSigma[Vbias] -> SetMarkerStyle(20);
      g_vs_NINOthr_CTR_ampCorr_effSigma[Vbias] -> SetMarkerColor(kRed);
      g_vs_NINOthr_CTR_ampCorr_effSigma[Vbias] -> SetLineColor(kRed);
      g_vs_NINOthr_CTR_ampCorr_effSigma[Vbias] -> SetLineStyle(1);
      g_vs_NINOthr_CTR_ampCorr_effSigma[Vbias] -> Draw("P,same");
      
      g_vs_NINOthr_CTR_distAmpCorr_effSigma[Vbias] -> SetMarkerStyle(20);
      g_vs_NINOthr_CTR_distAmpCorr_effSigma[Vbias] -> SetMarkerColor(kGreen+2);
      g_vs_NINOthr_CTR_distAmpCorr_effSigma[Vbias] -> SetLineColor(kGreen+2);
      g_vs_NINOthr_CTR_distAmpCorr_effSigma[Vbias] -> SetLineStyle(1);
      g_vs_NINOthr_CTR_distAmpCorr_effSigma[Vbias] -> Draw("P,same");

      
      g_vs_NINOthr_CTR_gausSigma[Vbias] -> SetMarkerStyle(24);
      g_vs_NINOthr_CTR_gausSigma[Vbias] -> SetMarkerColor(kBlack);
      g_vs_NINOthr_CTR_gausSigma[Vbias] -> SetLineColor(kBlack);
      g_vs_NINOthr_CTR_gausSigma[Vbias] -> SetLineStyle(7);
      g_vs_NINOthr_CTR_gausSigma[Vbias] -> Draw("P,same");
      
      g_vs_NINOthr_CTR_ampCorr_gausSigma[Vbias] -> SetMarkerStyle(24);
      g_vs_NINOthr_CTR_ampCorr_gausSigma[Vbias] -> SetMarkerColor(kRed);
      g_vs_NINOthr_CTR_ampCorr_gausSigma[Vbias] -> SetLineColor(kRed);
      g_vs_NINOthr_CTR_ampCorr_gausSigma[Vbias] -> SetLineStyle(7);
      g_vs_NINOthr_CTR_ampCorr_gausSigma[Vbias] -> Draw("P,same");
      
      g_vs_NINOthr_CTR_distAmpCorr_gausSigma[Vbias] -> SetMarkerStyle(24);
      g_vs_NINOthr_CTR_distAmpCorr_gausSigma[Vbias] -> SetMarkerColor(kGreen+2);
      g_vs_NINOthr_CTR_distAmpCorr_gausSigma[Vbias] -> SetLineColor(kGreen+2);
      g_vs_NINOthr_CTR_distAmpCorr_gausSigma[Vbias] -> SetLineStyle(7);
      g_vs_NINOthr_CTR_distAmpCorr_gausSigma[Vbias] -> Draw("P,same");

      latexLabel -> Draw("same");

      TLegend* legend = new TLegend(0.55,0.50,0.99,0.90);
      legend -> SetFillColor(0);
      legend -> SetFillStyle(1000);  
      legend -> SetTextFont(82);
      legend -> AddEntry(g_vs_NINOthr_CTR_effSigma[Vbias],            "                  raw - eff. #sigma","P");
      legend -> AddEntry(g_vs_NINOthr_CTR_ampCorr_effSigma[Vbias],    "      amp. walk corr. - eff. #sigma","P");
      legend -> AddEntry(g_vs_NINOthr_CTR_distAmpCorr_effSigma[Vbias],"amp. walk & pos. corr - eff. #sigma","P");
      legend -> AddEntry(g_vs_NINOthr_CTR_gausSigma[Vbias],            "                  raw - gaus #sigma","P");
      legend -> AddEntry(g_vs_NINOthr_CTR_ampCorr_gausSigma[Vbias],    "      amp. walk corr. - gaus #sigma","P");
      legend -> AddEntry(g_vs_NINOthr_CTR_distAmpCorr_gausSigma[Vbias],"amp. walk & pos. corr - gaus #sigma","P");
      legend -> Draw("same");

      TLatex* latexLabel2 = new TLatex(0.20,0.80,Form("V_{bias} = %.0f V",Vbias));
      latexLabel2 -> SetNDC();
      latexLabel2 -> SetTextFont(42);
      latexLabel2 -> SetTextSize(0.05);
      latexLabel2 -> Draw("same");
      
      c_vs_NINOthr -> Print(Form("c_Vbias%.0f_vs_NINOthr_%s.png",Vbias,ch.c_str()));
      c_vs_NINOthr -> Print(Form("c_Vbias%.0f_vs_NINOthr_%s.pdf",Vbias,ch.c_str()));
    }
    
    
    
    for(auto it : vec_Vbias)
    {
      float NINOthr = it.first;
      std::sort(vec_Vbias[NINOthr].begin(),vec_Vbias[NINOthr].end());
      
      TCanvas* c_vs_Vbias = new TCanvas(Form("c_NINOthr%.0f_vs_Vbias_%s",NINOthr,ch.c_str()),Form("c_NINOthr%.0f_vs_Vbias_%s",NINOthr,ch.c_str()),1000,500);
      c_vs_Vbias -> cd();
      
      TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,10.,1.5*vec_Vbias[NINOthr].at(vec_Vbias[NINOthr].size()-1),200.) );
      hPad -> SetTitle(";V_{bias} [V];time resolution [ps]");
      hPad -> Draw();
      
      
      g_vs_Vbias_CTR_effSigma[NINOthr] -> SetMarkerStyle(20);
      g_vs_Vbias_CTR_effSigma[NINOthr] -> SetMarkerColor(kBlack);
      g_vs_Vbias_CTR_effSigma[NINOthr] -> SetLineColor(kBlack);
      g_vs_Vbias_CTR_effSigma[NINOthr] -> SetLineStyle(1);
      g_vs_Vbias_CTR_effSigma[NINOthr] -> Draw("P,same");
      
      g_vs_Vbias_CTR_ampCorr_effSigma[NINOthr] -> SetMarkerStyle(20);
      g_vs_Vbias_CTR_ampCorr_effSigma[NINOthr] -> SetMarkerColor(kRed);
      g_vs_Vbias_CTR_ampCorr_effSigma[NINOthr] -> SetLineColor(kRed);
      g_vs_Vbias_CTR_ampCorr_effSigma[NINOthr] -> SetLineStyle(1);
      g_vs_Vbias_CTR_ampCorr_effSigma[NINOthr] -> Draw("P,same");
      
      g_vs_Vbias_CTR_distAmpCorr_effSigma[NINOthr] -> SetMarkerStyle(20);
      g_vs_Vbias_CTR_distAmpCorr_effSigma[NINOthr] -> SetMarkerColor(kGreen+2);
      g_vs_Vbias_CTR_distAmpCorr_effSigma[NINOthr] -> SetLineColor(kGreen+2);
      g_vs_Vbias_CTR_distAmpCorr_effSigma[NINOthr] -> SetLineStyle(1);
      g_vs_Vbias_CTR_distAmpCorr_effSigma[NINOthr] -> Draw("P,same");
      
      
      g_vs_Vbias_CTR_gausSigma[NINOthr] -> SetMarkerStyle(24);
      g_vs_Vbias_CTR_gausSigma[NINOthr] -> SetMarkerColor(kBlack);
      g_vs_Vbias_CTR_gausSigma[NINOthr] -> SetLineColor(kBlack);
      g_vs_Vbias_CTR_gausSigma[NINOthr] -> SetLineStyle(7);
      g_vs_Vbias_CTR_gausSigma[NINOthr] -> Draw("P,same");
      
      g_vs_Vbias_CTR_ampCorr_gausSigma[NINOthr] -> SetMarkerStyle(24);
      g_vs_Vbias_CTR_ampCorr_gausSigma[NINOthr] -> SetMarkerColor(kRed);
      g_vs_Vbias_CTR_ampCorr_gausSigma[NINOthr] -> SetLineColor(kRed);
      g_vs_Vbias_CTR_ampCorr_gausSigma[NINOthr] -> SetLineStyle(7);
      g_vs_Vbias_CTR_ampCorr_gausSigma[NINOthr] -> Draw("P,same");
      
      g_vs_Vbias_CTR_distAmpCorr_gausSigma[NINOthr] -> SetMarkerStyle(24);
      g_vs_Vbias_CTR_distAmpCorr_gausSigma[NINOthr] -> SetMarkerColor(kGreen+2);
      g_vs_Vbias_CTR_distAmpCorr_gausSigma[NINOthr] -> SetLineColor(kGreen+2);
      g_vs_Vbias_CTR_distAmpCorr_gausSigma[NINOthr] -> SetLineStyle(7);
      g_vs_Vbias_CTR_distAmpCorr_gausSigma[NINOthr] -> Draw("P,same");
      
      latexLabel -> Draw("same");
      
      TLegend* legend = new TLegend(0.55,0.50,0.99,0.90);
      legend -> SetFillColor(0);
      legend -> SetFillStyle(1000);  
      legend -> SetTextFont(82);
      legend -> AddEntry(g_vs_Vbias_CTR_effSigma[NINOthr],            "                  raw - eff. #sigma","P");
      legend -> AddEntry(g_vs_Vbias_CTR_ampCorr_effSigma[NINOthr],    "      amp. walk corr. - eff. #sigma","P");
      legend -> AddEntry(g_vs_Vbias_CTR_distAmpCorr_effSigma[NINOthr],"amp. walk & pos. corr - eff. #sigma","P");
      legend -> AddEntry(g_vs_Vbias_CTR_gausSigma[NINOthr],            "                  raw - gaus #sigma","P");
      legend -> AddEntry(g_vs_Vbias_CTR_ampCorr_gausSigma[NINOthr],    "      amp. walk corr. - gaus #sigma","P");
      legend -> AddEntry(g_vs_Vbias_CTR_distAmpCorr_gausSigma[NINOthr],"amp. walk & pos. corr - gaus #sigma","P");
      legend -> Draw("same");

      TLatex* latexLabel2 = new TLatex(0.20,0.80,Form("NINO thr. = %.0f mV",NINOthr));
      latexLabel2 -> SetNDC();
      latexLabel2 -> SetTextFont(42);
      latexLabel2 -> SetTextSize(0.05);
      latexLabel2 -> Draw("same");
      
      c_vs_Vbias -> Print(Form("c_NINOthr%.0f_vs_Vbias_%s.png",NINOthr,ch.c_str()));
      c_vs_Vbias -> Print(Form("c_NINOthr%.0f_vs_Vbias_%s.pdf",NINOthr,ch.c_str()));
    }
  }  
}
