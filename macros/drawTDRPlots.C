#include "/Users/abenagli/Work/TIMING/TBStudies/macros/CMS_lumi.h"

std::string plotDir = "/Users/abenagli/Work/TIMING/TBStudies/plots/TDR_v2";
TFile* inFile;
TFile* inFile1;
TFile* inFile2;
TFile* inFile3;
TFile* inFile4;
TFile* inFile5;
TLegend* legend;
TLegend* legend2;
TCanvas* c;
TEfficiency* e;
TEfficiency* e1;
TEfficiency* e2;
TEfficiency* e3;
TH1F* hPad;
TLatex* latexLabel;
TLatex* latexLabel2;
TLatex* latexLabel3;
TLatex* latex1;
TLatex* latex2;
TLatex* latex3;
TGraph* g;
TProfile* p1;
TProfile* p2;
TProfile* p3;
TProfile* p4;
TH1F* h1;
TH1F* h2;
TH1F* h3;

void drawEfficiency_vs_phi();
void drawEfficiency_vs_eta();
void drawEfficiency_vs_pt();
void drawOccupancy();
void drawEdep();
void drawTimeRes();


void DrawPadLabels(TPad* pad)
{
  writeExtraText = true;              // if extra text
  extraText  = "   Phase-2 Simulation";  // default extra text is "Preliminary"
  lumi_sqrtS = "";                    // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 0;             // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  
  // iPos=0  : out of frame
  // iPos=11 : top-left, left-aligned
  // iPos=22 : center, centered
  // iPos=33 : top-right, right-aligned
  // mode generally : iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
  int iPos = 0;
  
  CMS_lumi( pad,iPeriod,iPos );
  gPad -> Update();
}



void drawTDRPlots()
{
  drawEfficiency_vs_phi();
  // drawEfficiency_vs_eta();
  // drawEfficiency_vs_pt();
  // drawOccupancy();
  // drawEdep();
  //drawTimeRes();
}

void drawEfficiency_vs_phi()
{
  inFile1 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/SingleMuPtFlatPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");
  inFile2 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/SinglePiPtFlatPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");
  inFile3 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/MinBiasPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");

  c = new TCanvas("c_all_efficiency_vs_phi","c_all_efficiency_vs_phi");
  // gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(-0.0415,0.,0.4,1.25) );
  hPad -> SetTitle(";#phi at BTL;Efficiency");
  hPad -> Draw();
  DrawPadLabels(c);
  
  legend = new TLegend(0.30,0.80,0.80,0.93);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.05);

  for(int ii = 1; ii < 3; ++ii)
  {
    double x = 0.0415-0.006 + (2.*3.14159-4.*0.0415)/36. * ii;
    TLine* line = new TLine(x,hPad->GetMinimum(),x,1.);
    line -> SetLineWidth(3);
    line -> SetLineStyle(7);
    line -> Draw("same");

    for(int jj = 1; jj < 3; ++jj)
    {
      double y = x - (2.*3.14159-4.*0.0415)/36./3. * jj;
      std::cout << "y: " << y << std::endl;
      TLine* line = new TLine(y,hPad->GetMinimum(),y,1.);
      line -> SetLineWidth(2);
      line -> SetLineStyle(5);
      line -> Draw("same");      

      if( ii == 1 && jj == 1) legend -> AddEntry(line,"gap between SiPMs","L");
    }
    
    if( ii == 1 ) legend -> AddEntry(line,"gap between trays","L");
  }
  {
    double xMin = -(0.0415-0.006);
    double xMax = +(0.0415-0.006);
    TBox* box = new TBox(xMin,hPad->GetMinimum(),xMax,1.);
    box -> SetFillColor(kBlack);
    box -> SetFillStyle(3004);
    box -> Draw("same");
    legend -> AddEntry(box,"TST rail","F");
  }
  
  
  e1 = (TEfficiency*)( inFile1->Get(Form("p1_matchedRecHit_eff_vs_phiFold__pt00.8-10.0__Ethr1.0MeV")) );
  e2 = (TEfficiency*)( inFile2->Get(Form("p1_matchedRecHit_eff_vs_phiFold__pt00.8-10.0__Ethr1.0MeV")) );
  e3 = (TEfficiency*)( inFile3->Get(Form("p1_matchedRecHit_eff_vs_phiFold__pt00.8-10.0__Ethr1.0MeV")) );

  e1 -> SetMarkerSize(1.0);
  e1 -> SetMarkerColor(kBlue);
  e1 -> SetLineColor(kBlue);
  e1 -> Draw("same");  
  // e2 -> SetMarkerSize(1.0);
  // e2 -> SetMarkerColor(kRed);
  // e2 -> SetLineColor(kRed);
  // e2 -> Draw("same");
  // e3 -> SetMarkerSize(1.0);
  // e3 -> SetMarkerColor(kBlack);
  // e3 -> SetLineColor(kBlack);
  // e3 -> Draw("same");

  legend2 = new TLegend(0.15,0.89,0.40,0.93);
  legend2 -> SetFillColor(0);
  legend2 -> SetFillStyle(1000);  
  legend2 -> SetTextFont(42);
  legend2 -> SetTextSize(0.05);
  legend2 -> AddEntry(e1,"single muons","PL");
  // legend2 -> AddEntry(e2,"single pions","PL");
  // legend2 -> AddEntry(e3,"minimum bias","PL");
  
  legend -> Draw("same");
  legend2 -> Draw("same");
  
  c -> Print(Form("%s/c_all_efficiency_vs_phi.png",plotDir.c_str()));
  c -> Print(Form("%s/c_all_efficiency_vs_phi.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_all_efficiency_vs_phi.C",plotDir.c_str()));  
}



void drawEfficiency_vs_eta()
{
  inFile1 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/SingleMuPtFlatPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");
  inFile2 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/SinglePiPtFlatPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");
  inFile3 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/MinBiasPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");
  
  c = new TCanvas("c_all_efficiency_vs_eta","c_all_efficiency_vs_eta");
  // gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.6,1.5,1.1) );
  hPad -> SetTitle(";|#eta| at BTL;Efficiency");
  hPad -> Draw();
  DrawPadLabels(c);
  
  e1 = (TEfficiency*)( inFile1->Get(Form("p1_matchedRecHit_eff_vs_eta__pt00.8-10.0__Ethr1.0MeV")) );
  e2 = (TEfficiency*)( inFile2->Get(Form("p1_matchedRecHit_eff_vs_eta__pt00.8-10.0__Ethr1.0MeV")) );
  e3 = (TEfficiency*)( inFile3->Get(Form("p1_matchedRecHit_eff_vs_eta__pt00.8-10.0__Ethr1.0MeV")) );

  e1 -> SetMarkerSize(1.0);
  e1 -> SetMarkerColor(kBlue);
  e1 -> SetLineColor(kBlue);
  e1 -> Draw("same");  
  // e2 -> SetMarkerSize(1.0);
  // e2 -> SetMarkerColor(kRed);
  // e2 -> SetLineColor(kRed);
  // e2 -> Draw("same");
  // e3 -> SetMarkerSize(1.0);
  // e3 -> SetMarkerColor(kBlack);
  // e3 -> SetLineColor(kBlack);
  // e3 -> Draw("same");

  legend = new TLegend(0.30,0.80,0.85,0.93);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.05);

  TLine* line = new TLine(0.67,hPad->GetMinimum(),0.67,1.);
  line -> SetLineWidth(3);
  line -> SetLineStyle(7);
  line -> Draw("same");
  
  line = new TLine(1.15,hPad->GetMinimum(),1.15,1.);
  line -> SetLineWidth(3);
  line -> SetLineStyle(7);
  line -> Draw("same");
  
  legend -> AddEntry(line,"LYSO thickness change","L");
  legend -> Draw("same");
  
  legend2 = new TLegend(0.15,0.89,0.40,0.93);
  legend2 -> SetFillColor(0);
  legend2 -> SetFillStyle(1000);  
  legend2 -> SetTextFont(42);
  legend2 -> SetTextSize(0.05);
  legend2 -> AddEntry(e1,"single muons","PL");
  // legend2 -> AddEntry(e2,"single pions","PL");
  // legend2 -> AddEntry(e3,"minimum bias","PL");
  
//   legend2 -> Draw("same");
  
  c -> Print(Form("%s/c_all_efficiency_vs_eta.png",plotDir.c_str()));
  c -> Print(Form("%s/c_all_efficiency_vs_eta.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_all_efficiency_vs_eta.C",plotDir.c_str()));  
}



void drawEfficiency_vs_pt()
{
  inFile1 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/SingleMuPtFlatPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");
  inFile2 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/SinglePiPtFlatPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");
  inFile3 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/MinBiasPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");
  
  c = new TCanvas("c_all_efficiency_vs_pt","c_all_efficiency_vs_pt");
  // gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.8,0.6,10.,1.1) );
  hPad -> SetTitle(";p_{T} [GeV];Efficiency");
  hPad -> Draw();
  DrawPadLabels(c);
  
  e1 = (TEfficiency*)( inFile1->Get(Form("p1_matchedRecHit_eff_vs_pt__Ethr1.0MeV")) );
  e2 = (TEfficiency*)( inFile2->Get(Form("p1_matchedRecHit_eff_vs_pt__Ethr1.0MeV")) );
  e3 = (TEfficiency*)( inFile3->Get(Form("p1_matchedRecHit_eff_vs_pt__Ethr1.0MeV")) );

  e1 -> SetMarkerSize(1.0);
  e1 -> SetMarkerColor(kBlue);
  e1 -> SetLineColor(kBlue);
  e1 -> Draw("same");  
  e2 -> SetMarkerSize(1.0);
  e2 -> SetMarkerColor(kRed);
  e2 -> SetLineColor(kRed);
  e2 -> Draw("same");
  e3 -> SetMarkerSize(1.0);
  e3 -> SetMarkerColor(kBlack);
  e3 -> SetLineColor(kBlack);
  e3 -> Draw("same");
  
  legend2 = new TLegend(0.15,0.80,0.40,0.93);
  legend2 -> SetFillColor(0);
  legend2 -> SetFillStyle(1000);  
  legend2 -> SetTextFont(42);
  legend2 -> SetTextSize(0.05);
  legend2 -> AddEntry(e1,"single muons","PL");
  legend2 -> AddEntry(e2,"single pions","PL");
  legend2 -> AddEntry(e3,"minimum bias","PL");
  
  legend2 -> Draw("same");
  
  c -> Print(Form("%s/c_all_efficiency_vs_pt.png",plotDir.c_str()));
  c -> Print(Form("%s/c_all_efficiency_vs_pt.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_all_efficiency_vs_pt.C",plotDir.c_str()));  
}



void drawOccupancy()
{
  inFile = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/MinBiasPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");
  
  c = new TCanvas("c_simHits_occupancy","c_simHits_occupancy");
  gPad -> SetGridx();
  gPad -> SetGridy();
  gPad -> SetLogx();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.01,0.001,10.,0.30) );
  hPad -> SetTitle(";Readout energy threshold [MeV];PU200 crystal occupancy");
  hPad -> Draw();
  DrawPadLabels(c);
  
  int* colors = new int[6];
  colors[0] = kGreen; colors[1] = kGreen+1;
  colors[2] = kRed; colors[3] = kRed+1;
  colors[4] = kAzure+1; colors[5] = kAzure+2;

  legend = new TLegend(0.75,0.50,0.85,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextSize(0.05);
  legend -> SetTextFont(42);
  
  for(int iRU = 0; iRU < 6; ++iRU)
  {
    g = (TGraph*)( inFile->Get(Form("g_simHits_PU200_occ_vs_RU_energyCut_RU%d",iRU)) );
    g -> SetLineColor(colors[iRU]);
    g -> SetLineWidth(5);
    g -> SetLineStyle(1);
    g -> Draw("L,same");

    legend -> AddEntry(g,Form("RU%d",iRU),"L");
  }
  
  legend -> Draw("same");

  TLine* line_100keV = new TLine(0.1,hPad->GetMinimum(),0.1,hPad->GetMaximum());
  line_100keV -> SetLineWidth(2);
  line_100keV -> SetLineStyle(2);
  line_100keV -> Draw("same");
  
  TLine* line_1MeV = new TLine(1.,hPad->GetMinimum(),1.,hPad->GetMaximum());
  line_1MeV -> SetLineWidth(2);
  line_1MeV -> SetLineStyle(2);
  line_1MeV -> Draw("same");
  
  c -> Print(Form("%s/c_minBias_simHits_occupancy_PU200.png",plotDir.c_str()));
  c -> Print(Form("%s/c_minBias_simHits_occupancy_PU200.pdf",plotDir.c_str()));



  inFile1 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/SingleMuPtFlatPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");
  inFile2 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/SingleMuPtFlatPlots/barphiflat_PDE0.23_DCR56GHz/hitsPlots_barphiflat.root","READ");
  
  c = new TCanvas("c_eff_timeRes_vs_Ethr","c_eff_timeRes_vs_Ethr");
  gPad -> SetGridx();
  gPad -> SetGridy();
  gPad -> SetTicky(0);
  
  hPad = (TH1F*)( gPad->DrawFrame(0.1,0.,10.,1.) );
  hPad -> SetTitle(";Readout energy threshold [MeV];Efficiency");
  hPad -> Draw();
  DrawPadLabels(c);
  
  e1 = (TEfficiency*)( inFile1->Get("p1_matchedRecHit_eff_vs_Ethr") );
  e1 -> Draw("same");
  
  p1 = (TProfile*)( inFile2->Get("p1_matchedRecHit_timeRes_vs_Ethr_sumEnergy") );
  p1 -> SetMarkerColor(kRed);
  p1 -> SetMarkerStyle(23);
  p1 -> SetLineColor(kRed);
  float rightmax = 1.;
  // float scale = 1./p1->GetMaximum();
  float scale = 1./100;
  p1 -> Scale(scale);
  p1 -> Draw("same");

  p2 = (TProfile*)( inFile1->Get("p1_matchedRecHit_timeRes_vs_Ethr_sumEnergy") );
  p2 -> SetMarkerColor(kRed);
  p2 -> SetMarkerStyle(20);
  p2 -> SetLineColor(kRed);
  p2 -> Scale(scale);
  p2 -> Draw("same");
  
  TGaxis* axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),0.,100.,510,"+L");
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.05);
  axis->SetTitleFont(42);
  axis->SetTitleSize(0.05);
  axis->SetTitleColor(kRed);
  axis->SetTitle("#LT #sigma_{t} #GT [ps]");
  axis->Draw();
  
  gPad -> SetLogx();

  TLine* line = new TLine(1.,0.,1.,1.);
  line -> SetLineStyle(7);
  line -> Draw("same");

  legend2 = new TLegend(0.15,0.15,0.40,0.30);
  legend2 -> SetFillColor(0);
  legend2 -> SetFillStyle(1000);  
  legend2 -> SetTextFont(42);
  legend2 -> SetTextSize(0.05);
  legend2 -> AddEntry(p2,"#color[2]{#sigma_{t} at 0 fb^{-1}}","P");
  legend2 -> AddEntry(p1,"#color[2]{#sigma_{t} at 4000 fb^{-1}}","P");

  legend2 -> Draw("same");
  
  c -> Print(Form("%s/c_minBias_simHits_EthrOpt.png",plotDir.c_str()));
  c -> Print(Form("%s/c_minBias_simHits_EthrOpt.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_minBias_simHits_EthrOpt.C",plotDir.c_str()));
}




void drawEdep()
{
  inFile1 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd5/SingleMuPtFlatPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");
  inFile2 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd5/SinglePiPtFlatPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");
  inFile3 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd5/MinBiasPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat_40OOT.root","READ");
  
  c = new TCanvas("c_all_Edep","c_all_Edep");
  gPad -> SetGridx();
  gPad -> SetGridy();
  gPad -> SetLogy();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.00001,50,10.) );
  hPad -> SetTitle(";#Sigma E_{recHit}^{#DeltaR < 0.05} [MeV];Event fraction");
  hPad -> Draw();
  DrawPadLabels(c);
  
  h1 = (TH1F*)( inFile1->Get(Form("h1_matchedRecHit_totEnergy__pt00.8-10.0__Ethr1.0MeV")) );
  h2 = (TH1F*)( inFile2->Get(Form("h1_matchedRecHit_totEnergy__pt00.8-10.0__Ethr1.0MeV")) );
  h3 = (TH1F*)( inFile3->Get(Form("h1_matchedRecHit_totEnergy__pt00.8-10.0__Ethr1.0MeV")) );

  TH1F* h1_clone = (TH1F*)(h1->Clone());
  h1_clone -> GetXaxis() -> SetRangeUser(2.75,500.);
  h1 -> Scale(1./h1->Integral(h1->FindBin(2.),h1->FindBin(20.)));
  h1 -> Rebin(5);
  h1 -> SetLineColor(kBlue);
  h1 -> SetLineWidth(3);
  h1 -> Draw("hist,same");
  
  TH1F* h2_clone = (TH1F*)(h2->Clone());
  h2_clone -> GetXaxis() -> SetRangeUser(1.,500.);
  h2 -> Scale(1./h2->Integral(h2->FindBin(2.),h2->FindBin(20.)));
  h2 -> Rebin(5);
  h2 -> SetLineColor(kRed);
  h2 -> SetLineWidth(3);
  h2 -> Draw("hist,same");
  
  TH1F* h3_clone = (TH1F*)(h3->Clone());
  h3_clone -> GetXaxis() -> SetRangeUser(1.,500.);
  h3 -> Scale(1./h3->Integral(h3->FindBin(2.),h3->FindBin(20.)));
  h3 -> Rebin(5);
  h3 -> SetLineColor(kBlack);
  h3 -> SetLineWidth(3);
  h3 -> Draw("hist,same");  
  
  legend2 = new TLegend(0.15,0.80,0.55,0.93);
  legend2 -> SetFillColor(0);
  legend2 -> SetFillStyle(1000);  
  legend2 -> SetTextFont(42);
  legend2 -> SetTextSize(0.05);
  legend2 -> AddEntry(h1,Form("single muons"),"L");
  legend2 -> AddEntry(h2,Form("single pions"),"L");
  legend2 -> AddEntry(h3,Form("minimum bias"),"L");
  // legend2 -> AddEntry(h1,Form("single muons,  #LT E #GT = %.1f MeV",h1_clone->GetMean()),"L");
  // legend2 -> AddEntry(h2,Form("single pions,    #LT E #GT = %.1f MeV",h2_clone->GetMean()),"L");
  // legend2 -> AddEntry(h3,Form("minimum bias, #LT E #GT = %.1f MeV",h3_clone->GetMean()),"L");
  
  legend2 -> Draw("same");
  
  c -> Print(Form("%s/c_all_Edep.png",plotDir.c_str()));
  c -> Print(Form("%s/c_all_Edep.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_all_Edep.C",plotDir.c_str()));
  
  
  
  c = new TCanvas("c_all_Edep_vs_eta","c_all_Edep_vs_eta");
  gPad -> SetGridx();
  gPad -> SetGridy();
  // gPad -> SetLogy();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,1.5,16.) );
  hPad -> SetTitle(";|#eta| at BTL;#LT #Sigma E_{recHit}^{#DeltaR < 0.05} #GT [MeV]");
  hPad -> Draw();
  DrawPadLabels(c);
  
  p1 = (TProfile*)( inFile1->Get(Form("p1_matchedRecHit_totEnergy_vs_eta__pt00.8-10.0__Ethr1.0MeV")) );
  p2 = (TProfile*)( inFile2->Get(Form("p1_matchedRecHit_totEnergy_vs_eta__pt00.8-10.0__Ethr1.0MeV")) );
  p3 = (TProfile*)( inFile3->Get(Form("p1_matchedRecHit_totEnergy_vs_eta__pt00.8-10.0__Ethr1.0MeV")) );
  
  p1 -> SetLineColor(kBlue);
  p1 -> SetLineWidth(3);
  p1 -> SetMarkerSize(0);
  p1 -> DrawCopy("hist,same");
  p1 -> SetFillColor(kBlue);
  p1 -> SetFillStyle(3001);
  p1 -> Draw("E2,same");
  
  p2 -> SetLineColor(kRed);
  p2 -> SetLineWidth(3);
  p2 -> SetMarkerSize(0);
  p2 -> DrawCopy("hist,same");
  p2 -> SetFillColor(kRed);
  p2 -> SetFillStyle(3001);
  p2 -> Draw("E2,same");
  
  p3 -> SetLineColor(kBlack);
  p3 -> SetLineWidth(3);
  p3 -> SetMarkerSize(0);
  p3 -> DrawCopy("hist,same");
  p3 -> SetFillColor(kBlack);
  p3 -> SetFillStyle(3001);
  p3 -> Draw("E2,same");
  
  
  legend2 = new TLegend(0.15,0.80,0.55,0.93);
  legend2 -> SetFillColor(0);
  legend2 -> SetFillStyle(1000);  
  legend2 -> SetTextFont(42);
  legend2 -> SetTextSize(0.05);
  legend2 -> AddEntry(h1,"single muons","L");
  legend2 -> AddEntry(h2,"single pions","L");
  legend2 -> AddEntry(h3,"minimum bias","L");
  
  legend2 -> Draw("same");
  
  c -> Print(Form("%s/c_all_Edep_vs_eta.png",plotDir.c_str()));
  c -> Print(Form("%s/c_all_Edep_vs_eta.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_all_Edep_vs_eta.C",plotDir.c_str()));
  
  
  
  c = new TCanvas("c_singleMuPtFlat_maxEnergy_vs_eta","c_singleMuPtFlat_maxEnergy_vs_eta");
  gPad -> SetGridx();
  gPad -> SetGridy();
  gPad -> SetTicky(0);
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,1.5,16.) );
  hPad -> SetTitle(";|#eta| at BTL;#LT #Sigma E_{recHit}^{#DeltaR < 0.05} #GT [MeV]");
  hPad -> Draw();
  DrawPadLabels(c);
  
  legend = new TLegend(0.40,0.30-2*0.05,0.83,0.30);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.05);
  
  // latexLabel2 = new TLatex(0.15,0.23,Form("#splitline{E_{recHit}^{thr} = %.1f MeV}{p_{T} #in [0.8,%.1f] GeV}",1.,10.));
  // latexLabel2 -> SetNDC();
  // latexLabel2 -> SetTextFont(42);
  // latexLabel2 -> SetTextSize(0.05);
  // latexLabel2 -> Draw("same");
  
  p1 = (TProfile*)( inFile1->Get(Form("p1_matchedRecHit_totEnergy_vs_eta__pt00.8-%04.1f__Ethr%.1fMeV",10.,1.)) );
  p1 -> SetMarkerSize(1.0);
  p1 -> SetMarkerColor(kBlack);
  p1 -> SetLineColor(kBlack);
  p1 -> SetLineWidth(1);
  p1 -> Draw("same");
  legend -> AddEntry(p1,Form("total energy"),"PE");
  
  p1 = (TProfile*)( inFile1->Get(Form("p1_matchedRecHit_maxEnergy_vs_eta__pt00.8-%04.1f__Ethr%.1fMeV",10.,1.)) );
  p1 -> SetMarkerSize(1.0);
  p1 -> SetMarkerStyle(24);
  p1 -> SetMarkerColor(kBlack);
  p1 -> SetLineColor(kBlack);
  p1 -> Draw("same");
  legend -> AddEntry(p1,Form("seed recHit energy"),"PE");
  
  float rightmax = 1.;
  float scale = gPad->GetUymax()/rightmax;
  
  p1 = (TProfile*)( inFile1->Get(Form("p1_matchedRecHit_maxOverTotEnergy_vs_eta__pt00.8-%04.1f__Ethr%.1fMeV",10.,1.)) );
  p1 -> SetLineColor(kBlack);
  p1 -> SetLineWidth(3);
  p1->SetLineColor(kRed);
  p1->Scale(scale);
  p1->Draw("hist,same");
  
  //draw an axis on the right side
  TGaxis* axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),0,rightmax,510,"+L");
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.05);
  axis->SetTitleFont(42);
  axis->SetTitleSize(0.06);
  axis->SetTitleColor(kRed);
  axis->SetTitle("E_{recHit}^{seed} / #Sigma E_{recHit}");
  axis->Draw();

  latexLabel3 = new TLatex(0.65,0.90,Form("single muons"));
  latexLabel3 -> SetNDC();
  latexLabel3 -> SetTextFont(42);
  latexLabel3 -> SetTextSize(0.05);
//   latexLabel3 -> Draw("same");
  
  legend2 = new TLegend(0.30,0.88,0.85,0.93);
  legend2 -> SetFillColor(0);
  legend2 -> SetFillStyle(1000);  
  legend2 -> SetTextFont(42);
  legend2 -> SetTextSize(0.05);
  
  TLine* line = new TLine(0.67,hPad->GetMinimum(),0.67,16.);
  line -> SetLineWidth(2);
  line -> SetLineStyle(7);
  line -> Draw("same");
  
  line = new TLine(1.15,hPad->GetMinimum(),1.15,16.);
  line -> SetLineWidth(2);
  line -> SetLineStyle(7);
  line -> Draw("same");
  
  legend2 -> AddEntry(line,"LYSO thickness change","L");
  legend2 -> Draw("same");
  
//   legend -> Draw("same");
  
  c -> Print(Form("%s/c_singleMuPtFlat_maxEnergy_vs_eta.png",plotDir.c_str()));
  c -> Print(Form("%s/c_singleMuPtFlat_maxEnergy_vs_eta.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_singleMuPtFlat_maxEnergy_vs_eta.C",plotDir.c_str()));
  

  
  c = new TCanvas("c_minBias_maxEnergy_vs_eta","c_minBias_maxEnergy_vs_eta");
  gPad -> SetGridx();
  gPad -> SetGridy();
  gPad -> SetTicky(0);
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,1.5,16.) );
  hPad -> SetTitle(";|#eta| at BTL;#LT #Sigma E_{recHit}^{#DeltaR < 0.05} #GT [MeV]");
  hPad -> Draw();
  DrawPadLabels(c);
  
  // latexLabel2 = new TLatex(0.15,0.23,Form("#splitline{E_{recHit}^{thr} = %.1f MeV}{p_{T} #in [0.8,%.1f] GeV}",1.,10.));
  // latexLabel2 -> SetNDC();
  // latexLabel2 -> SetTextFont(42);
  // latexLabel2 -> SetTextSize(0.05);
  // latexLabel2 -> Draw("same");
  
  p3 = (TProfile*)( inFile3->Get(Form("p1_matchedRecHit_totEnergy_vs_eta__pt00.8-%04.1f__Ethr%.1fMeV",10.,1.)) );
  p3 -> SetMarkerSize(1.0);
  p3 -> SetMarkerColor(kBlack);
  p3 -> SetLineColor(kBlack);
  p3 -> SetLineWidth(1);
  p3 -> Draw("same");
  
  p3 = (TProfile*)( inFile3->Get(Form("p1_matchedRecHit_maxEnergy_vs_eta__pt00.8-%04.1f__Ethr%.1fMeV",10.,1.)) );
  p3 -> SetMarkerSize(1.0);
  p3 -> SetMarkerStyle(24);
  p3 -> SetMarkerColor(kBlack);
  p3 -> SetLineColor(kBlack);
  p3 -> Draw("same");
  
  rightmax = 1.;
  scale = gPad->GetUymax()/rightmax;
  
  p3 = (TProfile*)( inFile3->Get(Form("p1_matchedRecHit_maxOverTotEnergy_vs_eta__pt00.8-%04.1f__Ethr%.1fMeV",10.,1.)) );
  p3 -> SetLineColor(kBlack);
  p3 -> SetLineWidth(3);
  p3->SetLineColor(kRed);
  p3->Scale(scale);
  p3->Draw("hist,same");
  
  //draw an axis on the right side
  axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),0,rightmax,510,"+L");
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.05);
  axis->SetTitleFont(42);
  axis->SetTitleSize(0.06);
  axis->SetTitleColor(kRed);
  axis->SetTitle("E_{recHit}^{seed} / #Sigma E_{recHit}");
  axis->Draw();

  latexLabel3 = new TLatex(0.55,0.90,Form("minimum bias events"));
  latexLabel3 -> SetNDC();
  latexLabel3 -> SetTextFont(42);
  latexLabel3 -> SetTextSize(0.05);
//   latexLabel3 -> Draw("same");

  line = new TLine(0.67,hPad->GetMinimum(),0.67,16.);
  line -> SetLineWidth(2);
  line -> SetLineStyle(7);
  line -> Draw("same");
  
  line = new TLine(1.15,hPad->GetMinimum(),1.15,16.);
  line -> SetLineWidth(2);
  line -> SetLineStyle(7);
  line -> Draw("same");

  legend2 -> Draw("same");
  legend -> Draw("same");
  
  c -> Print(Form("%s/c_minBias_maxEnergy_vs_eta.png",plotDir.c_str()));
  c -> Print(Form("%s/c_minBias_maxEnergy_vs_eta.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_minBias_maxEnergy_vs_eta.C",plotDir.c_str()));  
}




void drawTimeRes()
{
  //inFile1 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/MinBiasPlots/barphiflat_PDE0.38_DCR0GHz/hitsPlots_barphiflat.root","READ");
  inFile1 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR_v2/CMSSW_10_4_0_mtd5/MinBiasPlots/barphiflat_PDE0.38_DCR0GHz/hitsPlots_barphiflat.root","READ");
  inFile2 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR_v2/CMSSW_10_4_0_mtd5/MinBiasPlots/barphiflat_PDE0.33_DCR33GHz/hitsPlots_barphiflat.root","READ");
  inFile3 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR_v2/CMSSW_10_4_0_mtd5/MinBiasPlots/barphiflat_PDE0.24_DCR55GHz/hitsPlots_barphiflat.root","READ");
  inFile4 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR_v2/CMSSW_10_4_0_mtd5/MinBiasPlots/barphiflat_PDE0.21_DCR63GHz/hitsPlots_barphiflat.root","READ");
  // inFile1 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/MinBiasPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat.root","READ");
  // inFile2 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd5/MinBiasPlots/barphiflat_PDE0.39_DCR0GHz/hitsPlots_barphiflat_40OOT.root","READ");
  // inFile3 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR_v2/CMSSW_10_4_0_mtd5/MinBiasPlots/barphiflat_PDE0.38_DCR0GHz/hitsPlots_barphiflat.root","READ");
  inFile5 = TFile::Open("/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/MinBiasPlots/barphiflat_PDE0.38_DCR0GHz/hitsPlots_barphiflat.root","READ");
  
  c = new TCanvas("c_MinBias_timeRes","c_MinBias_timeRes");
  
  gPad -> SetGridx();
  gPad -> SetGridy();

  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,1.5,120.) );
  hPad -> SetTitle(";|#eta| at BTL;#sigma_{t} [ps]");
  hPad -> Draw();
  DrawPadLabels(c);
    
  legend = new TLegend(0.20,0.92-5*0.05,0.85,0.92);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.05);
  
  p1 = (TProfile*)( inFile1->Get(Form("p1_matchedRecHit_timeRes_vs_eta__sumEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",0.8,10.,1.)) );
  p1 -> SetLineColor(kBlue);
  p1 -> SetLineWidth(3);
  p1 -> Draw("hist,same");
  legend -> AddEntry(p1,Form("      0 fb^{-1}, parametric model"),"L");
  
  p2 = (TProfile*)( inFile2->Get(Form("p1_matchedRecHit_timeRes_vs_eta__sumEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",0.8,10.,1.)) );
  p2 -> SetLineColor(kMagenta);
  p2 -> SetLineWidth(3);
  p2 -> Draw("hist,same");
  legend -> AddEntry(p2,Form("1000 fb^{-1}, parametric model"),"L");

  p3 = (TProfile*)( inFile3->Get(Form("p1_matchedRecHit_timeRes_vs_eta__sumEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",0.8,10.,1.)) );
  p3 -> SetLineColor(kOrange+1);
  p3 -> SetLineWidth(3);
  p3 -> Draw("hist,same");
  legend -> AddEntry(p3,Form("3000 fb^{-1}, parametric model"),"L");
  
  p4 = (TProfile*)( inFile4->Get(Form("p1_matchedRecHit_timeRes_vs_eta__sumEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",0.8,10.,1.)) );
  p4 -> SetLineColor(kGreen+1);
  p4 -> SetLineWidth(3);
  p4 -> Draw("hist,same");
  legend -> AddEntry(p4,Form("4000 fb^{-1}, parametric model"),"L");
  
  legend -> Draw("same");

  TGraphErrors* graph = new TGraphErrors();
  std::vector<float> etaRanges;
  etaRanges.push_back(0.);
  etaRanges.push_back(0.125);
  etaRanges.push_back(0.250);
  etaRanges.push_back(0.375);
  etaRanges.push_back(0.500);
  etaRanges.push_back(0.625);
  etaRanges.push_back(0.750);
  etaRanges.push_back(0.875);
  etaRanges.push_back(1.000);
  etaRanges.push_back(1.125);
  etaRanges.push_back(1.250);
  etaRanges.push_back(1.375);
  etaRanges.push_back(1.500);

  for(unsigned int jj = 0; jj < etaRanges.size()-1; ++jj)
  {
    float etaMin = etaRanges.at(jj);
    float etaMax = etaRanges.at(jj+1);
    h1 = (TH1F*)( inFile5->Get(Form("h1_matchedRecHit_timeSum_energyWeighted__eta%.2f-%.2f__Ethr%.1fMeV",etaMin,etaMax,1.)) ); 

    TF1* f_gaus = new TF1("f_gaus","gaus(0)",-25.,25.);
    h1->Fit(f_gaus,"QNRS+","",h1->GetMean()-0.5*h1->GetRMS(),h1->GetMean()+h1->GetRMS());
    
    graph -> SetPoint(jj,0.5*(etaMin+etaMax),1000.*f_gaus->GetParameter(2));
    graph -> SetPointError(jj,0.5*(etaMax-etaMin),1000.*f_gaus->GetParError(2));
  }

  graph -> SetMarkerStyle(20);
  graph -> Draw("PL,same");

  legend -> AddEntry(graph,Form("full simulation"),"PE");
  
  c -> Print(Form("%s/c_MinBias_timeRes_vs_eta.png",plotDir.c_str()));
  c -> Print(Form("%s/c_MinBias_timeRes_vs_eta.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_MinBias_timeRes_vs_eta.C",plotDir.c_str()));
}
