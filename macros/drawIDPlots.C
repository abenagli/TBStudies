#include "/Users/abenagli/Work/TIMING/TBStudies/macros/CMS_lumi.h"

//std::string inFileName = "/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd5/IDPlots/IDPlots_eleID_evalMVA.root";
//std::string plotDir = "/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd5/IDPlots/";
//std::string plotLabel = "eleID_evalMVA";

std::string inFileName = "/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/IDPlots/IDPlots_eleID.root";
std::string plotDir = "/Users/abenagli/Work/TIMING/TBStudies/plots/TDR/CMSSW_10_4_0_mtd3/IDPlots/";
std::string plotLabel = "eleID";

TFile* inFile;
TLegend* legend;
TLegend* legend2;
TCanvas* c;
TH1F* hPad;
TLatex* latexLabel1;
TLatex* latexLabel2;
TGraph* g1;
TGraph* g2;
TGraph* g3;
TGraph* g4;
TGraph* g5;
TGraph* g6;
TProfile* p1;
TProfile* p2;
TH1F* h1;
TH1F* h2;
TH1F* h3;



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



void drawHistogram1D(const std::string& label,
                     const std::string& title, const int& rebin,
                     const float& xMin, const float& yMin, const float& xMax, const float& yMax,
                     TLatex* latex = NULL)
{
  inFile = TFile::Open(inFileName.c_str(),"READ");
  
  c = new TCanvas(Form("c_%s",label.c_str()),Form("c_%s",label.c_str()));
//   gPad -> SetGridx();
//   gPad -> SetGridy();
  gPad -> SetLogy();
  
  hPad = (TH1F*)( gPad->DrawFrame(xMin,yMin,xMax,yMax) );
  hPad -> SetTitle(title.c_str());
  hPad -> Draw();
  DrawPadLabels(c);
  
//   h1 = (TH1F*)( inFile->Get(Form("h1_promptNoPU_%s",label.c_str())) );
//   h1 -> Rebin(rebin);
//   h1 -> SetLineColor(kRed-2);
//   h1 -> SetLineWidth(2);
//   h1 -> SetFillColor(kRed-2);
//   h1 -> SetFillStyle(3003);
//   h1 -> Scale(1./h1->Integral());
//   h1 -> Draw("hist,same");
  
  h2 = (TH1F*)( inFile->Get(Form("h1_ele_%s",label.c_str())) );
  h2 -> Rebin(rebin);
  h2 -> SetLineColor(kRed);
  h2 -> SetLineWidth(4);
  h2 -> Scale(1./h2->Integral());
  h2 -> Draw("hist,same");

  h3 = (TH1F*)( inFile->Get(Form("h1_pi_%s",label.c_str())) );
  h3 -> Rebin(rebin);
  h3 -> SetLineColor(kBlue);
  h3 -> SetLineWidth(4);
  h3 -> Scale(1./h3->Integral());
  h3 -> Draw("hist,same");
  
  legend2 = new TLegend(0.50,0.75,0.85,0.83);
  legend2 -> SetFillColor(0);
  legend2 -> SetFillStyle(1000);  
  legend2 -> SetTextFont(42);
  legend2 -> SetTextSize(0.05);
//   legend2 -> AddEntry(h1,"prompt electrons (no PU)","L");
//   legend2 -> AddEntry(h2,"prompt electrons (PU200)","L");
//   legend2 -> AddEntry(h3,"fake electrons (PU200)",  "L");
  legend2 -> AddEntry(h2,"single electrons","L");
  legend2 -> AddEntry(h3,"single pions",  "L");
  
  hPad -> GetYaxis() -> SetTitle(Form("Event fraction / %.1e MeV",h2->GetBinWidth(1)));
  legend2 -> Draw("same");

  if( latex ) latex -> Draw("same");
  
  c -> Print(Form("%s/c_%s_%s.png",plotDir.c_str(),plotLabel.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_%s.pdf",plotDir.c_str(),plotLabel.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_%s.C",  plotDir.c_str(),plotLabel.c_str(),label.c_str()));
}




void drawROC(const std::string& label,
             const std::string& title,
             const float& xMin, const float& yMin, const float& xMax, const float& yMax,
             TLatex* latex = NULL)
{
  inFile = TFile::Open(inFileName.c_str(),"READ");
  
  c = new TCanvas(Form("c_%s",label.c_str()),Form("c_%s",label.c_str()));
//   gPad -> SetGridx();
//   gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(xMin,yMin,xMax,yMax) );
  hPad -> SetTitle(title.c_str());
  hPad -> Draw();
  DrawPadLabels(c);
  
  g1 = (TGraph*)( inFile->Get(Form("g_ROC_%s",label.c_str())) );
  g1 -> SetLineColor(kRed);
  g1 -> SetLineWidth(3);
  g1 -> Draw("L,same");

  g2 = (TGraph*)( inFile->Get(Form("g_ROC2_%s",label.c_str())) );
  g2 -> SetLineColor(kBlue);
  g2 -> SetLineWidth(3);
  g2 -> Draw("L,same");
  
//   g2 = (TGraph*)( inFile->Get(Form("g_ROC_eleID_MTDOnly_%s",label.c_str())) );
//   g2 -> SetLineColor(kBlue);
//   g2 -> SetLineWidth(3);
//   g2 -> Draw("L,same");
  
  legend2 = new TLegend(0.17,0.75,0.40,0.83);
  legend2 -> SetFillColor(0);
  legend2 -> SetFillStyle(1000);  
  legend2 -> SetTextFont(42);
  legend2 -> SetTextSize(0.05);
  legend2 -> AddEntry(g1,"old MVA","L");
  legend2 -> AddEntry(g2,"new MVA (old BDT output + MTD)","L");
  
  g3 = (TGraph*)( inFile->Get(Form("g_ROC_eleID_eleOnlyNoIso_%s",label.c_str())) );
  if( g3 )
  {
    g3 -> SetLineColor(kGreen);
    g3 -> SetLineWidth(2);
    g3 -> Draw("L,same");
    legend2 -> AddEntry(g3,"new MVA (ele. vars, w/o isolation)","L");
  }

  g4 = (TGraph*)( inFile->Get(Form("g_ROC_eleID_eleOnlyIso_%s",label.c_str())) );
  if( g4 )
  {
    g4 -> SetLineColor(kMagenta);
    g4 -> SetLineWidth(2);
    g4 -> Draw("L,same");
    legend2 -> AddEntry(g4,"new MVA (ele. vars, w/ isolation)","L");
  }
  
  g5 = (TGraph*)( inFile->Get(Form("g_ROC_eleID_eleMTDNoIso_%s",label.c_str())) );
  if( g5 )
  {
    g5 -> SetLineColor(kGreen);
    g5 -> SetLineWidth(2);
    g5 -> SetLineStyle(2);
    g5 -> Draw("L,same");
    legend2 -> AddEntry(g5,"new MVA (ele. vars + MTD, w/o isolation)","L");
  }
  
  g6 = (TGraph*)( inFile->Get(Form("g_ROC_eleID_eleMTDIso_%s",label.c_str())) );
  if( g6 )
  {
    g6 -> SetLineColor(kMagenta);
    g6 -> SetLineWidth(2);
    g6 -> SetLineStyle(2);
    g6 -> Draw("L,same");
    legend2 -> AddEntry(g6,"new MVA (ele. vars + MTD, w/ isolation)","L");
  }
  
  legend2 -> Draw("same");

  if( latex ) latex -> Draw("same");
  
  c -> Print(Form("%s/c_%s_ROC_%s.png",plotDir.c_str(),plotLabel.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_ROC_%s.pdf",plotDir.c_str(),plotLabel.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_ROC_%s.C",  plotDir.c_str(),plotLabel.c_str(),label.c_str()));
}



void drawIDPlots()
{
  latexLabel1 = new TLatex(0.55,0.85,Form("Barrel"));
  latexLabel1 -> SetNDC();
  latexLabel1 -> SetTextFont(42);
  latexLabel1 -> SetTextSize(0.05);
  
  latexLabel2 = new TLatex(0.55,0.85,Form("Endcap"));
  latexLabel2 -> SetNDC();
  latexLabel2 -> SetTextFont(42);
  latexLabel2 -> SetTextSize(0.05);
  
// //   drawHistogram1D("electrons_t_BTL",               ";t_{0} [ns];",                            1,-1., 0.00001,1.,  1.,  latexLabel1);
// //   drawHistogram1D("electrons_pathLength_BTL",      ";track length [mm];",                     1, 0., 0.00001,500.,1.,  latexLabel1);
//   drawHistogram1D("electrons_mva_BTL",             ";mva;",                                   5,-1., 0.00001,1.,  1.,  latexLabel1);
//   drawHistogram1D("electrons_trackRelIso_BTL",     ";I_{track}/p_{T};",                       5, 0., 0.00001,0.5, 1., latexLabel1);
//   drawHistogram1D("matchedRecHit_energySum_BTL",   ";#Sigma E_{recHit} [MeV];",               5, 0., 0.0002, 50., 0.5, latexLabel1);
   drawHistogram1D("matchedRecHit_energySeed_BTL",  ";E_{recHit}^{seed} [MeV];",               4, 0., 0.001,  20., 3.,  latexLabel1);
//   drawHistogram1D("matchedRecHit_energyRatio_BTL", ";E_{recHit}^{seed} / #Sigma E_{recHit};", 2,-1., 0.001,  1.1,  3., latexLabel1);
//   drawHistogram1D("matchedRecHit_n_BTL",           ";N_{recHits};",                           1,-0.5,0.001,  19.5, 3., latexLabel1);
//   drawHistogram1D("matchedRecHit_sieie_BTL",       ";#sigma_{i#etai#eta};",                   2,-1., 0.001,  10., 3.,  latexLabel1);
//   drawHistogram1D("matchedRecHit_sipip_BTL",       ";#sigma_{i#phii#phi};",                   2,-1., 0.0001, 20., 3.,  latexLabel1);
  
// //   drawHistogram1D("electrons_t_ETL",               ";t_{0} [ns];",                            1,-1., 0.00001,1.,  1.,  latexLabel1);
// //   drawHistogram1D("electrons_pathLength_ETL",      ";track length [mm];",                     1, 0., 0.00001,500.,1.,  latexLabel1);
//   drawHistogram1D("electrons_mva_ETL",             ";mva;",                                   5,-1., 0.00001,1.,  1.,  latexLabel2);  
//   drawHistogram1D("electrons_trackRelIso_ETL",     ";I_{track}/p_{T};",                       5, 0., 0.00001,0.5, 1., latexLabel2);
//   drawHistogram1D("matchedRecHit_energySum_ETL",   ";#Sigma E_{recHit} [MeV];",               4, 0., 0.0002, 2.,  3., latexLabel2);
//   drawHistogram1D("matchedRecHit_energySeed_ETL",  ";E_{recHit}^{seed} [MeV];",               4, 0., 0.001,  2.,  3.,  latexLabel2);
//   drawHistogram1D("matchedRecHit_energyRatio_ETL", ";E_{recHit}^{seed} / #Sigma E_{recHit};", 2,-1., 0.001,  1.1,  3., latexLabel2);
//   drawHistogram1D("matchedRecHit_n_ETL",           ";N_{recHits};",                           1,-0.5, 0.001,  19.5, 3.,latexLabel2);
//   drawHistogram1D("matchedRecHit_sieie_ETL",       ";#sigma_{i#etai#eta};",                   2,-1., 0.001,  0.1, 3.,  latexLabel2);
//   drawHistogram1D("matchedRecHit_sipip_ETL",       ";#sigma_{i#phii#phi};",                   2,-1., 0.0001, 0.1, 3.,  latexLabel2);
  
  latexLabel1 = new TLatex(0.20,0.85,Form("Barrel"));
  latexLabel1 -> SetNDC();
  latexLabel1 -> SetTextFont(42);
  latexLabel1 -> SetTextSize(0.05);
  
  latexLabel2 = new TLatex(0.20,0.85,Form("Endcap"));
  latexLabel2 -> SetNDC();
  latexLabel2 -> SetTextFont(42);
  latexLabel2 -> SetTextSize(0.05);
  
  drawROC("EB",";Electron efficiency;Pion efficiency", 0.85, 0., 1.005, 1.01, latexLabel1);
//   drawROC("EB",";#epsilon_{prompt};#epsilon_{fake}", 0.9, 0., 1.01, 0.31, latexLabel1);
//   drawROC("EE",";#epsilon_{prompt};#epsilon_{fake}", 0.9, 0., 1.01, 0.03, latexLabel2);
}
