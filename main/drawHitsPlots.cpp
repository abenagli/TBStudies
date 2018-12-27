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
#include "TGaxis.h"



int main(int argc, char** argv)
{
  setTDRStyle();

  if( argc < 2 )
  {
    std::cout << ">>> studyHits::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  
  
  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);

  int nRUs = opts.GetOpt<int>("Input.nRUs");
  
  float PDE = opts.GetOpt<float>("Options.PDE");
  float DCR = opts.GetOpt<float>("Options.DCR");
  
  std::string label = opts.GetOpt<std::string>("Input.label");
  std::string plotDir = opts.GetOpt<std::string>("Input.plotDir");
  plotDir += std::string(Form("%s_PDE%.2f_DCR%.0fGHz",label.c_str(),PDE,DCR));
  
  std::vector<std::string> particleLabel = opts.GetOpt<std::vector<std::string> >("Input.particleLabel");
  TFile* inFile = TFile::Open(Form("%s/hitsPlots_%s.root",plotDir.c_str(),label.c_str()),"READ");
  

  
  std::string particleLabel2 = "";
  for(auto jj : particleLabel) particleLabel2 += jj + " ";
  TLatex* latexLabel = NULL;
  if( label == "tile" )
    latexLabel = new TLatex(0.16,0.96,Form("11.5#times11.5 mm^{2} tiles -- %s",particleLabel2.c_str()));
  if( label == "barphi" )
    latexLabel = new TLatex(0.16,0.96,Form("3#times50 mm^{2} bars along #phi -- %s",particleLabel2.c_str()));
  if( label == "barz" )
    latexLabel = new TLatex(0.16,0.96,Form("3#times50 mm^{2} bars along z -- %s",particleLabel2.c_str()));
  if( label == "barzflat" )
    latexLabel = new TLatex(0.16,0.96,Form("3#times56 mm^{2} flat bars along z -- %s",particleLabel2.c_str()));
  latexLabel -> SetNDC();
  latexLabel -> SetTextFont(42);
  latexLabel -> SetTextSize(0.03);

  TCanvas* c;
  TH1F* hPad;
  int* colors;
  TLegend* legend;

  TH1F* h1;
  TEfficiency* e;
  TProfile* p;
  TGraph* g;

  TLatex* latexLabel2;
  
  std::vector<float> ptRanges = opts.GetOpt<std::vector<float> >("Options.ptRanges");
  std::vector<float> etaRanges = opts.GetOpt<std::vector<float> >("Options.etaRanges");
  std::vector<float> EthrVals = opts.GetOpt<std::vector<float> >("Options.EthrVals");
  
  float etaMin = opts.GetOpt<float>("Options.etaMin");
  float etaMax = opts.GetOpt<float>("Options.etaMax");
  
  float phiMin = opts.GetOpt<float>("Options.phiMin");
  float phiMax = opts.GetOpt<float>("Options.phiMax");
  
  
  

  
  c = new TCanvas("c_tracks","c_trakcs",2500,1200);
  c -> Divide(2,1);

  c -> cd(1);
  gPad -> SetLogy();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.00001,10.,0.1) );
  hPad -> SetTitle(";track p_{T} [GeV];event fraction");
  hPad -> Draw();
  
  h1 = (TH1F*)( inFile->Get(Form("h1_tracks_pt")) );
  h1 -> Scale(1./h1->Integral());
  h1 -> SetLineColor(kRed);
  h1 -> SetFillColor(kRed);
  h1 -> SetFillStyle(3003);
  h1 -> SetLineWidth(2);
  h1 -> SetLineStyle(1);
  h1 -> Draw("hist,same");
  
  latexLabel2 = new TLatex(0.50,0.80,Form("#LT p_{T} #GT = %.1f GeV",h1->GetMean()));
  latexLabel2 -> SetNDC();
  latexLabel2 -> SetTextFont(82);
  latexLabel2 -> SetTextSize(0.04);
  latexLabel2 -> Draw("same");
  
  c -> cd(2);
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(-4.0,0.001,4.0,0.05) );
  hPad -> SetTitle(";track #eta;event fraction");
  hPad -> Draw();
  
  h1 = (TH1F*)( inFile->Get(Form("h1_tracks_eta")) );
  h1 -> Scale(1./h1->Integral());
  h1 -> SetLineColor(kRed);
  h1 -> SetFillColor(kRed);
  h1 -> SetFillStyle(3003);
  h1 -> SetLineWidth(2);
  h1 -> SetLineStyle(1);
  h1 -> Draw("hist,same");
  
  c -> Print(Form("%s/c_%s_tracks.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_tracks.pdf",plotDir.c_str(),label.c_str()));
  delete c;
  
  
  
  c = new TCanvas("c_simHits_occupancy","c_simHits_occupancy",1400,1200);
  gPad -> SetLogx();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.01,0.001,10.,0.30) );
  hPad -> SetTitle(";energy threshold [MeV];PU200 occupancy");
  hPad -> Draw();
  
  colors = new int[9];
  colors[0] = kGreen; colors[1] = kGreen+1; colors[2] = kGreen+2;
  colors[3] = kRed; colors[4] = kRed+1; colors[5] = kRed+2;
  colors[6] = kAzure+1; colors[7] = kAzure+2; colors[8] = kAzure+3;

  legend = new TLegend(0.81,0.50,0.99,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);
  
  for(int iRU = 0; iRU < nRUs; ++iRU)
  {
    g = (TGraph*)( inFile->Get(Form("g_simHits_PU200_occ_vs_energyCut_RU%d",iRU)) );
    g -> SetLineColor(colors[iRU]);
    g -> SetLineWidth(5);
    g -> SetLineStyle(1);
    g -> Draw("L,same");

    legend -> AddEntry(g,Form("RU%d",iRU),"L");
  }
  
  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  c -> Print(Form("%s/c_%s_simHits_occupancy_PU200.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_simHits_occupancy_PU200.pdf",plotDir.c_str(),label.c_str()));
  delete c;

  
  c = new TCanvas("c_recHits_occupancy","c_recHits_occupancy",1400,1200);
  gPad -> SetLogx();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.01,0.001,10.,0.30) );
  hPad -> SetTitle(";energy threshold [MeV];PU200 occupancy");
  hPad -> Draw();
  
  colors = new int[9];
  colors[0] = kGreen; colors[1] = kGreen+1; colors[2] = kGreen+2;
  colors[3] = kRed; colors[4] = kRed+1; colors[5] = kRed+2;
  colors[6] = kAzure+1; colors[7] = kAzure+2; colors[8] = kAzure+3;

  legend = new TLegend(0.81,0.50,0.99,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);
  
  for(int iRU = 0; iRU < nRUs; ++iRU)
  {
    g = (TGraph*)( inFile->Get(Form("g_recHits_PU200_occ_vs_energyCut_RU%d",iRU)) );
    g -> SetLineColor(colors[iRU]);
    g -> SetLineWidth(5);
    g -> SetLineStyle(1);
    g -> Draw("L,same");

    legend -> AddEntry(g,Form("RU%d",iRU),"L");
  }
  
  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  c -> Print(Form("%s/c_%s_recHits_occupancy_PU200.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_recHits_occupancy_PU200.pdf",plotDir.c_str(),label.c_str()));
  delete c;  
  
  
  
  

  c = new TCanvas(Form("c_trackerRes"),Form("c_trackerRes"),2500,1200);
  c -> Divide(2,1);

  c -> cd(1);
  gPad -> SetGridx();
  gPad -> SetGridy();
  gPad -> SetLogy();
  
  hPad = (TH1F*)( gPad->DrawFrame(etaMin,0.05,etaMax,10.) );
  hPad -> SetTitle(";|#eta|;#sigma_{R#DeltaPhi} [mm]");
  hPad -> Draw();
  
  legend = new TLegend(0.60,0.16,0.83,0.16+(ptRanges.size()-1)*0.04);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);    
  legend -> SetTextSize(0.02);
  
  for(unsigned int jj = 0; jj < ptRanges.size()-1; ++jj)
  {
    float ptRangesMin = ptRanges.at(jj);
    float ptRangesMax = ptRanges.at(jj+1);

    TGraphErrors* g = new TGraphErrors();
    
    for(unsigned int jj = 0; jj < etaRanges.size()-1; ++jj)
    {
      float etaRangesMin = etaRanges.at(jj);
      float etaRangesMax = etaRanges.at(jj+1);    
      
      h1 = (TH1F*)( inFile->Get(Form("h1_matchedSimHit_track_RDphi__pt%04.1f-%04.1f__eta%02.1f-%02.1f",ptRangesMin,ptRangesMax,etaRangesMin,etaRangesMax)) );
      TF1* func = (TF1*)( h1->GetFunction("gaus") );
      g -> SetPoint(g->GetN(),0.5*(etaRangesMin+etaRangesMax),10.*func->GetParameter(2));
      g -> SetPointError(g->GetN()-1,0.,10.*func->GetParError(2));
    }
    
    g -> SetMarkerSize(1.0);
    g -> SetMarkerColor(51+int(48/(ptRanges.size()-2))*jj);
    g -> SetLineColor(51+int(48/(ptRanges.size()-2))*jj);
    g -> SetMarkerSize(1.);
    g -> SetLineWidth(1);
    g -> Draw("PL,same");
    
    legend -> AddEntry(g,Form("p_{T} #in [%.1f,%.1f] GeV",ptRangesMin,ptRangesMax),"PL");
  }

  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  
  c -> cd(2);
  gPad -> SetGridx();
  gPad -> SetGridy();
  gPad -> SetLogy();
  
  hPad = (TH1F*)( gPad->DrawFrame(etaMin,0.05,etaMax,10.) );
  hPad -> SetTitle(";|#eta|;#sigma_{#Deltaz} [mm]");
  hPad -> Draw();
  
  legend = new TLegend(0.60,0.16,0.83,0.16+(ptRanges.size()-1)*0.04);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);    
  legend -> SetTextSize(0.02);
  
  for(unsigned int jj = 0; jj < ptRanges.size()-1; ++jj)
  {
    float ptRangesMin = ptRanges.at(jj);
    float ptRangesMax = ptRanges.at(jj+1);

    TGraphErrors* g = new TGraphErrors();
    
    for(unsigned int jj = 0; jj < etaRanges.size()-1; ++jj)
    {
      float etaRangesMin = etaRanges.at(jj);
      float etaRangesMax = etaRanges.at(jj+1);    
      
      h1 = (TH1F*)( inFile->Get(Form("h1_matchedSimHit_track_Dz__pt%04.1f-%04.1f__eta%02.1f-%02.1f",ptRangesMin,ptRangesMax,etaRangesMin,etaRangesMax)) );
      TF1* func = (TF1*)( h1->GetFunction("gaus") );
      g -> SetPoint(g->GetN(),0.5*(etaRangesMin+etaRangesMax),10.*func->GetParameter(2));
      g -> SetPointError(g->GetN()-1,0.,10.*func->GetParError(2));
    }
    
    g -> SetMarkerSize(1.0);
    g -> SetMarkerColor(51+int(48/(ptRanges.size()-2))*jj);
    g -> SetLineColor(51+int(48/(ptRanges.size()-2))*jj);
    g -> SetMarkerSize(1.);
    g -> SetLineWidth(1);
    g -> Draw("PL,same");
    
    legend -> AddEntry(g,Form("p_{T} #in [%.1f,%.1f] GeV",ptRangesMin,ptRangesMax),"PL");
  }

  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  c -> Print(Form("%s/c_%s_trackerRes.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_trackerRes.pdf",plotDir.c_str(),label.c_str()));
  delete c;
  
  
  
  
  
  
  c = new TCanvas("c_maxEnergy_vs_eta","c_maxEnergy_vs_eta",1400,1200);
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(etaMin,0.,etaMax,12.) );
  hPad -> SetTitle(";|#eta|;#Sigma E_{recHit} [MeV]");
  hPad -> Draw();
  
  legend = new TLegend(0.50,0.16,0.83,0.16+2*0.04);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);
  legend -> SetTextSize(0.02);
  
  latexLabel2 = new TLatex(0.18,0.20,Form("#splitline{E_{thr} = %.1f MeV}{p_{T} #in [0.8,%.1f] GeV }",1.,ptRanges.at(ptRanges.size()-1)));
  latexLabel2 -> SetNDC();
  latexLabel2 -> SetTextFont(82);
  latexLabel2 -> SetTextSize(0.02);
  latexLabel2 -> Draw("same");

  p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_totEnergy_vs_eta__pt00.8-%04.1f__Ethr%.1fMeV",ptRanges.at(ptRanges.size()-1),1.)) );
  p -> SetMarkerSize(1.0);
  p -> SetMarkerColor(kBlack);
  p -> SetLineColor(kBlack);
  p -> Draw("same");
  legend -> AddEntry(p,Form("total energy"),"PL");
  
  p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_maxEnergy_vs_eta__pt00.8-%04.1f__Ethr%.1fMeV",ptRanges.at(ptRanges.size()-1),1.)) );
  p -> SetMarkerSize(1.0);
  p -> SetMarkerStyle(24);
  p -> SetMarkerColor(kBlack);
  p -> SetLineColor(kBlack);
  p -> Draw("same");
  legend -> AddEntry(p,Form("max. recHit energy"),"PL");
  
  float rightmax = 1.;
  float scale = gPad->GetUymax()/rightmax;
  
  p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_maxOverTotEnergy_vs_eta__pt00.8-%04.1f__Ethr%.1fMeV",ptRanges.at(ptRanges.size()-1),1.)) );
  p -> SetLineColor(kBlack);
  p -> SetLineWidth(2);
  p->SetLineColor(kRed);
  p->Scale(scale);
  p->Draw("hist,same");
  
  //draw an axis on the right side
  TGaxis* axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),0,rightmax,510,"+L");
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  axis->SetLabelFont(42);
  axis->SetLabelSize(0.04);
  axis->SetTitleFont(42);
  axis->SetTitleSize(0.06);
  axis->SetTitleColor(kRed);
  axis->SetTitle("E_{recHit}^{max} / #Sigma E_{recHit}");
  axis->Draw();
  
  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  c -> Print(Form("%s/c_%s_maxEnergy_vs_eta.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_maxEnergy_vs_eta.pdf",plotDir.c_str(),label.c_str()));
  delete c;

  
  for(unsigned int jj = 0; jj < ptRanges.size()-1; ++jj)
  {
    float ptMin = ptRanges.at(jj);
    float ptMax = ptRanges.at(jj+1);
    
    c = new TCanvas(Form("c_maxEnergy_vs_eta__pt%04.1f-%04.1f",ptMin,ptMax),Form("c_maxEnergy_vs_eta__pt%04.1f-%04.1f",ptMin,ptMax),1400,1200);
    gPad -> SetGridx();
    gPad -> SetGridy();
    
    hPad = (TH1F*)( gPad->DrawFrame(etaMin,0.,etaMax,12.) );
    hPad -> SetTitle(";|#eta|;#Sigma E_{recHit} [MeV]");
    hPad -> Draw();
    
    legend = new TLegend(0.50,0.16,0.83,0.16+2*0.04);
    legend -> SetFillColor(0);
    legend -> SetFillStyle(1000);  
    legend -> SetTextFont(82);
    legend -> SetTextSize(0.02);
    
    latexLabel2 = new TLatex(0.18,0.20,Form("#splitline{E_{thr} = %.1f MeV}{p_{T} #in [%.1f,%.1f] GeV }",1.,ptMin,ptMax));
    latexLabel2 -> SetNDC();
    latexLabel2 -> SetTextFont(82);
    latexLabel2 -> SetTextSize(0.02);
    latexLabel2 -> Draw("same");

    p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_totEnergy_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,1.)) );
    p -> SetMarkerSize(1.0);
    p -> SetMarkerColor(kBlack);
    p -> SetLineColor(kBlack);
    p -> Draw("same");
    legend -> AddEntry(p,Form("total energy"),"PL");
    
    p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_maxEnergy_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,1.)) );
    p -> SetMarkerSize(1.0);
    p -> SetMarkerStyle(24);
    p -> SetMarkerColor(kBlack);
    p -> SetLineColor(kBlack);
    p -> Draw("same");
    legend -> AddEntry(p,Form("max. recHit energy"),"PL");
    
    float rightmax = 1.;
    float scale = gPad->GetUymax()/rightmax;
    
    p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_maxOverTotEnergy_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,1.)) );
    p -> SetLineColor(kBlack);
    p -> SetLineWidth(2);
    p->SetLineColor(kRed);
    p->Scale(scale);
    p->Draw("hist,same");
    
    //draw an axis on the right side
    TGaxis* axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),0,rightmax,510,"+L");
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->SetLabelFont(42);
    axis->SetLabelSize(0.04);
    axis->SetTitleFont(42);
    axis->SetTitleSize(0.06);
    axis->SetTitleColor(kRed);
    axis->SetTitle("E_{recHit}^{max} / #Sigma E_{recHit}");
    axis->Draw();
    
    legend -> Draw("same");
    latexLabel -> Draw("same");
    
    c -> Print(Form("%s/c_%s_maxEnergy_vs_eta__pt%04.1f-%04.1f.png",plotDir.c_str(),label.c_str(),ptMin,ptMax));
    c -> Print(Form("%s/c_%s_maxEnergy_vs_eta__pt%04.1f-%04.1f.pdf",plotDir.c_str(),label.c_str(),ptMin,ptMax));
    delete c;    
  }
  
  
  
  
  
  
  c = new TCanvas("c_timeRes_vs_eta","c_timeRes_vs_eta",1400,1200);
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(etaMin,0.,etaMax,120.) );
  hPad -> SetTitle(";|#eta|;#sigma_{t} [ps]");
  hPad -> Draw();
  
  legend = new TLegend(0.50,0.87-2*0.05,0.83,0.87);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(82);
  legend -> SetTextSize(0.04);
  
  latexLabel2 = new TLatex(0.17,0.21,Form("#splitline{E_{thr} = %.1f MeV}{#splitline{p_{T} #in [0.8,%.1f] GeV}{PDE = %.0f%%, DCR = %.0fGHz}}",1.,ptRanges.at(ptRanges.size()-1),100.*PDE,DCR));
  latexLabel2 -> SetNDC();
  latexLabel2 -> SetTextFont(82);
  latexLabel2 -> SetTextSize(0.04);
  latexLabel2 -> Draw("same");
  
  p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_timeRes_vs_eta__totEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",0.8,ptRanges.at(ptRanges.size()-1),1.)) );
  p -> SetMarkerSize(1.5);
  p -> SetMarkerColor(kRed);
  p -> SetLineColor(kRed);
  p -> Draw("same");
  legend -> AddEntry(p,Form("single crystal"),"PL");
  
  // p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_timeRes_vs_eta__maxEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",0.8,ptRanges.at(ptRanges.size()-1),,1.)) );
  // p -> SetMarkerSize(1.0);
  // p -> SetMarkerStyle(24);
  // p -> SetMarkerColor(kBlack);
  // p -> SetLineColor(kBlack);
  // p -> Draw("same");
  // legend -> AddEntry(p,Form("max. recHit energy"),"PL");
  
  p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_timeRes_vs_eta__sumEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",0.8,ptRanges.at(ptRanges.size()-1),1.)) );
  p -> SetMarkerSize(1.5);
  p -> SetMarkerStyle(22);
  p -> SetMarkerColor(kBlue);
  p -> SetLineColor(kBlue);
  p -> Draw("same");
  legend -> AddEntry(p,Form("sum of recHits"),"PL");
  
  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  c -> Print(Form("%s/c_%s_timeRes_vs_eta.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_timeRes_vs_eta.pdf",plotDir.c_str(),label.c_str()));
  delete c;
  
  
  for(unsigned int jj = 0; jj < ptRanges.size()-1; ++jj)
  {
    float ptMin = ptRanges.at(jj);
    float ptMax = ptRanges.at(jj+1);
    
    c = new TCanvas(Form("c_timeRes_vs_eta__pt%04.1f-%04.1f",ptMin,ptMax),Form("c_timeRes_vs_eta__pt%04.1f-%04.1f",ptMin,ptMax),1400,1200);
    gPad -> SetGridx();
    gPad -> SetGridy();
    
    hPad = (TH1F*)( gPad->DrawFrame(etaMin,0.,etaMax,120.) );
    hPad -> SetTitle(";|#eta|;#sigma_{t} [ps]");
    hPad -> Draw();
    
    legend = new TLegend(0.50,0.87-2*0.05,0.83,0.87);
    legend -> SetFillColor(0);
    legend -> SetFillStyle(1000);
    legend -> SetTextFont(82);
    legend -> SetTextSize(0.04);
    
    latexLabel2 = new TLatex(0.17,0.21,Form("#splitline{E_{thr} = %.1f MeV}{#splitline{p_{T} #in [%.1f,%.1f] GeV}{PDE = %.0f%%, DCR = %.0fGHz}}",1.,ptMin,ptMax,100.*PDE,DCR));
    latexLabel2 -> SetNDC();
    latexLabel2 -> SetTextFont(82);
    latexLabel2 -> SetTextSize(0.04);
    latexLabel2 -> Draw("same");
    
    p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_timeRes_vs_eta__totEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,1.)) );
    p -> SetMarkerSize(1.5);
    p -> SetMarkerColor(kRed);
    p -> SetLineColor(kRed);
    p -> Draw("same");
    legend -> AddEntry(p,Form("single crystal"),"PL");
    
    // p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_timeRes_vs_eta__maxEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,1.)) );
    // p -> SetMarkerSize(1.0);
    // p -> SetMarkerStyle(24);
    // p -> SetMarkerColor(kBlack);
    // p -> SetLineColor(kBlack);
    // p -> Draw("same");
    // legend -> AddEntry(p,Form("max. recHit energy"),"PL");
    
    p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_timeRes_vs_eta__sumEnergy__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,1.)) );
    p -> SetMarkerSize(1.5);
    p -> SetMarkerStyle(22);
    p -> SetMarkerColor(kBlue);
    p -> SetLineColor(kBlue);
    p -> Draw("same");
    legend -> AddEntry(p,Form("sum of recHits"),"PL");
    
    legend -> Draw("same");
    latexLabel -> Draw("same");
    
    c -> Print(Form("%s/c_%s_timeRes_vs_eta__pt%04.1f-%04.1f.png",plotDir.c_str(),label.c_str(),ptMin,ptMax));
    c -> Print(Form("%s/c_%s_timeRes_vs_eta__pt%04.1f-%04.1f.pdf",plotDir.c_str(),label.c_str(),ptMin,ptMax));
    delete c;    
  }
  
  
  
  
  
  
  // for(unsigned int jj = 0; jj < ptRanges.size()-1; ++jj)
  // {
  //   float ptMin = ptRanges.at(jj);
  //   float ptMax = ptRanges.at(jj+1);
    
  //   c = new TCanvas(Form("c_maxEnergy_vs_local_x__pt%04.1f-%04.1f",ptMin,ptMax),Form("c_maxEnergy_vs_local_x__pt%04.1f-%04.1f",ptMin,ptMax),1400,1200);
  //   gPad -> SetGridx();
  //   gPad -> SetGridy();

  //   if( label == "tile" )     hPad = (TH1F*)( gPad->DrawFrame(-1., 0.,1., 12.) );
  //   if( label == "barphi" )   hPad = (TH1F*)( gPad->DrawFrame(-5., 0.,5., 12.) );
  //   if( label == "barz" )     hPad = (TH1F*)( gPad->DrawFrame(-0.2,0.,0.2,12.) );
  //   if( label == "barzflat" ) hPad = (TH1F*)( gPad->DrawFrame(-0.2,0.,0.2,12.) );
  //   hPad -> SetTitle(";local x (#phi direction) [cm];#Sigma E_{simHit} [MeV]");
  //   hPad -> Draw();
    
  //   legend = new TLegend(0.50,0.16,0.83,0.16+2*0.04);
  //   legend -> SetFillColor(0);
  //   legend -> SetFillStyle(1000);  
  //   legend -> SetTextFont(82);
  //   legend -> SetTextSize(0.02);
    
  //   latexLabel2 = new TLatex(0.18,0.20,Form("#splitline{E_{thr} = %.1f MeV}{p_{T} #in [%.1f,%.1f] GeV }",1.,ptMin,ptMax));
  //   latexLabel2 -> SetNDC();
  //   latexLabel2 -> SetTextFont(82);
  //   latexLabel2 -> SetTextSize(0.02);
  //   latexLabel2 -> Draw("same");
    
  //   p = (TProfile*)( inFile->Get(Form("p1_matchedSimHit_totEnergy_vs_local_x__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,1.)) );
  //   p -> SetMarkerSize(1.0);
  //   p -> SetMarkerColor(kBlack);
  //   p -> SetLineColor(kBlack);
  //   p -> Draw("same");
  //   legend -> AddEntry(p,Form("total energy"),"PL");
    
  //   p = (TProfile*)( inFile->Get(Form("p1_matchedSimHit_maxEnergy_vs_local_x__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,1.)) );
  //   p -> SetMarkerSize(1.0);
  //   p -> SetMarkerStyle(24);
  //   p -> SetMarkerColor(kBlack);
  //   p -> SetLineColor(kBlack);
  //   p -> Draw("same");
  //   legend -> AddEntry(p,Form("max. simHit energy"),"PL");
    
  //   float rightmax = 1.;
  //   float scale = gPad->GetUymax()/rightmax;
    
  //   p = (TProfile*)( inFile->Get(Form("p1_matchedSimHit_maxOverTotEnergy_vs_local_x__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,1.)) );
  //   p -> SetLineColor(kBlack);
  //   p -> SetLineWidth(2);
  //   p->SetLineColor(kRed);
  //   p->Scale(scale);
  //   p->Draw("hist,same");
    
  //   //draw an axis on the right side
  //   TGaxis* axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),0,rightmax,510,"+L");
  //   axis->SetLineColor(kRed);
  //   axis->SetLabelColor(kRed);
  //   axis->SetLabelFont(42);
  //   axis->SetLabelSize(0.04);
  //   axis->SetTitleFont(42);
  //   axis->SetTitleSize(0.06);
  //   axis->SetTitleColor(kRed);
  //   axis->SetTitle("E_{simHit}^{max} / #Sigma E_{simHit}");
  //   axis->Draw();
    
  //   legend -> Draw("same");
  //   latexLabel -> Draw("same");
    
  //   c -> Print(Form("%s/c_%s_maxEnergy_vs_local_x__pt%04.1f-%04.1f.png",plotDir.c_str(),label.c_str(),ptMin,ptMax));
  //   c -> Print(Form("%s/c_%s_maxEnergy_vs_local_x__pt%04.1f-%04.1f.pdf",plotDir.c_str(),label.c_str(),ptMin,ptMax));
  //   delete c;    
  // }
  
  // for(unsigned int jj = 0; jj < ptRanges.size()-1; ++jj)
  // {
  //   float ptMin = ptRanges.at(jj);
  //   float ptMax = ptRanges.at(jj+1);
    
  //   c = new TCanvas(Form("c_maxEnergy_vs_local_y__pt%04.1f-%04.1f",ptMin,ptMax),Form("c_maxEnergy_vs_local_y__pt%04.1f-%04.1f",ptMin,ptMax),1400,1200);
  //   gPad -> SetGridx();
  //   gPad -> SetGridy();

  //   if( label == "tile" )  hPad = (TH1F*)( gPad->DrawFrame(-1., 0.,1., 12.) );
  //   if( label == "barphi") hPad = (TH1F*)( gPad->DrawFrame(-0.2,0.,0.2,12.) );
  //   if( label == "barz" )  hPad = (TH1F*)( gPad->DrawFrame(-5., 0.,5., 12.) );
  //   if( label == "barzflat" )  hPad = (TH1F*)( gPad->DrawFrame(-5., 0.,5., 12.) );
  //   hPad -> SetTitle(";local y (#eta direction) [cm];#Sigma E_{simHit} [MeV]");
  //   hPad -> Draw();
    
  //   legend = new TLegend(0.50,0.16,0.83,0.16+2*0.04);
  //   legend -> SetFillColor(0);
  //   legend -> SetFillStyle(1000);  
  //   legend -> SetTextFont(82);
  //   legend -> SetTextSize(0.02);
    
  //   latexLabel2 = new TLatex(0.18,0.20,Form("#splitline{E_{thr} = %.1f MeV}{p_{T} #in [%.1f,%.1f] GeV }",1.,ptMin,ptMax));
  //   latexLabel2 -> SetNDC();
  //   latexLabel2 -> SetTextFont(82);
  //   latexLabel2 -> SetTextSize(0.02);
  //   latexLabel2 -> Draw("same");
    
  //   p = (TProfile*)( inFile->Get(Form("p1_matchedSimHit_totEnergy_vs_local_y__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,1.)) );
  //   p -> SetMarkerSize(1.0);
  //   p -> SetMarkerColor(kBlack);
  //   p -> SetLineColor(kBlack);
  //   p -> Draw("same");
  //   legend -> AddEntry(p,Form("total energy"),"PL");
    
  //   p = (TProfile*)( inFile->Get(Form("p1_matchedSimHit_maxEnergy_vs_local_y__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,1.)) );
  //   p -> SetMarkerSize(1.0);
  //   p -> SetMarkerStyle(24);
  //   p -> SetMarkerColor(kBlack);
  //   p -> SetLineColor(kBlack);
  //   p -> Draw("same");
  //   legend -> AddEntry(p,Form("max. simHit energy"),"PL");
    
  //   float rightmax = 1.;
  //   float scale = gPad->GetUymax()/rightmax;
    
  //   p = (TProfile*)( inFile->Get(Form("p1_matchedSimHit_maxOverTotEnergy_vs_local_y__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,1.)) );
  //   p -> SetLineColor(kBlack);
  //   p -> SetLineWidth(2);
  //   p->SetLineColor(kRed);
  //   p->Scale(scale);
  //   p->Draw("hist,same");
    
  //   //draw an axis on the right side
  //   TGaxis* axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),gPad->GetUymax(),0,rightmax,510,"+L");
  //   axis->SetLineColor(kRed);
  //   axis->SetLabelColor(kRed);
  //   axis->SetLabelFont(42);
  //   axis->SetLabelSize(0.04);
  //   axis->SetTitleFont(42);
  //   axis->SetTitleSize(0.06);
  //   axis->SetTitleColor(kRed);
  //   axis->SetTitle("E_{simHit}^{max} / #Sigma E_{simHit}");
  //   axis->Draw();
    
  //   legend -> Draw("same");
  //   latexLabel -> Draw("same");
    
  //   c -> Print(Form("%s/c_%s_maxEnergy_vs_local_y__pt%04.1f-%04.1f.png",plotDir.c_str(),label.c_str(),ptMin,ptMax));
  //   c -> Print(Form("%s/c_%s_maxEnergy_vs_local_y__pt%04.1f-%04.1f.pdf",plotDir.c_str(),label.c_str(),ptMin,ptMax));
  //   delete c;    
  // }
  
  
  
  
  
  
  c = new TCanvas("c_recHits_energy","c_recHits_energy",1400,1200);
  gPad -> SetLogy();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.0001,20.,1.) );
  hPad -> SetTitle(";recHit energy [MeV];event fraction");
  hPad -> Draw();
  
  legend = new TLegend(0.31,0.80,0.80,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);
  legend -> SetTextSize(0.02);
  
  h1 = (TH1F*)( inFile->Get("h1_recHit_energy") );
  h1 -> Scale(1./h1->Integral());
  h1 -> SetLineColor(kBlack);
  h1 -> SetLineWidth(3);
  h1 -> Draw("hist,same");
  legend -> AddEntry(h1,Form("all recHits"),"L");

  h1 = (TH1F*)( inFile->Get("h1_matchedRecHit_energy") );
  h1 -> Scale(1./h1->Integral());
  h1 -> SetLineColor(kRed);
  h1 -> SetLineWidth(3);
  h1 -> Draw("hist,same");
  legend -> AddEntry(h1,Form("track-matched recHits"),"L");
  
  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  c -> Print(Form("%s/c_%s_recHits_energy.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_recHits_energy.pdf",plotDir.c_str(),label.c_str()));
  delete c;
  

  
  c = new TCanvas("c_matchedRecHits_energySum","c_matchedRecHits_energySum",1400,1200);
  gPad -> SetLogy();
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.0001,20.,1.) );
  hPad -> SetTitle(";#Sigma E_{recHit} [MeV];event fraction");
  hPad -> Draw();
  
  legend = new TLegend(0.31,0.80,0.80,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);
  legend -> SetTextSize(0.02);
  
  h1 = (TH1F*)( inFile->Get("h1_matchedRecHit_energySum") );
  h1 -> Scale(1./h1->Integral());
  h1 -> SetLineColor(kBlack);
  h1 -> SetLineWidth(3);
  h1 -> Draw("hist,same");
  legend -> AddEntry(h1,Form("track-matched recHits"),"L");

  h1 = (TH1F*)( inFile->Get("h1_matchedRecHit_energySumCorr") );
  h1 -> Scale(1./h1->Integral());
  h1 -> SetLineColor(kRed);
  h1 -> SetLineWidth(3);
  h1 -> Draw("hist,same");
  legend -> AddEntry(h1,Form("track-matched recHits (slant thicnkess corr.)"),"L");
  
  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  c -> Print(Form("%s/c_%s_matchedRecHits_energySum.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_matchedRecHits_energySum.pdf",plotDir.c_str(),label.c_str()));
  delete c;
  
  
  
  c = new TCanvas("c_matchedRecHits_time_vs_eta","c_matchedRecHits_time_vs_eta",2500,1200);
  c -> Divide(2,1);
  
  c -> cd(1);

  h1 = (TH1F*)( inFile->Get("h1_matchedRecHit_time") );
  hPad = (TH1F*)( gPad->DrawFrame(0.,1.,20.,10.*h1->GetMaximum()) );
  hPad -> SetTitle(";time [ns];entries");
  hPad -> Draw();
  
  gPad -> SetGridx();
  gPad -> SetGridy();
  gPad -> SetLogy();

  legend = new TLegend(0.31,0.80,0.80,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);
  legend -> SetTextSize(0.02);
  
  h1 = (TH1F*)( inFile->Get("h1_recHit_time") );
  h1 -> SetLineColor(kBlack);
  h1 -> SetMarkerColor(kBlack);
  h1 -> Draw("hist,same");
  legend -> AddEntry(h1,Form("all recHits"),"L");
  
  h1 = (TH1F*)( inFile->Get("h1_matchedRecHit_time") );
  h1 -> SetLineColor(kRed);
  h1 -> SetMarkerColor(kRed);
  h1 -> Draw("hist,same");
  legend -> AddEntry(h1,Form("track-matched recHits"),"L");

  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  c -> cd(2);
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(etaMin,0.,etaMax,10.) );
  hPad -> SetTitle(";|#eta|;time [ns]");
  hPad -> Draw();
  
  p = (TProfile*)( inFile->Get("p1_matchedRecHit_time_vs_eta") );
  p -> SetLineColor(kRed);
  p -> SetMarkerColor(kRed);
  p -> Draw("same");
  
  latexLabel -> Draw("same");
  
  c -> Print(Form("%s/c_%s_matchedRecHits_time.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_matchedRecHits_time.pdf",plotDir.c_str(),label.c_str()));
  delete c;





  
  c = new TCanvas("c_efficiency_vs_eta","c_efficiency_vs_eta",1400,1200);
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(etaMin,0.,etaMax,1.1) );
  hPad -> SetTitle(";|#eta|;efficiency");
  hPad -> Draw();

  legend = new TLegend(0.60,0.16,0.83,0.16+EthrVals.size()*0.04);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);
  legend -> SetTextSize(0.02);
  
  int EthrIt = 0;
  for(auto Ethr : EthrVals)
  {
    e = (TEfficiency*)( inFile->Get(Form("p1_matchedRecHit_eff_vs_eta__pt00.8-%04.1f__Ethr%.1fMeV",ptRanges.at(ptRanges.size()-1),Ethr)) );
    
    e -> SetMarkerSize(1.0);
    e -> SetMarkerColor(51+int(48/(EthrVals.size()-1))*EthrIt);
    e -> SetLineColor(51+int(48/(EthrVals.size()-1))*EthrIt);
    e -> Draw("same");

    legend -> AddEntry(e,Form("E_{thr} = %.1f MeV",Ethr),"PL");
    ++EthrIt;
  }
  
  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  c -> Print(Form("%s/c_%s_efficiency_vs_eta.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_efficiency_vs_eta.pdf",plotDir.c_str(),label.c_str()));
  delete c;
  


  EthrIt = 0;
  for(auto Ethr : EthrVals)
  {
    c = new TCanvas(Form("c_efficiency_vs_eta_Ethr%.1fMeV",Ethr),Form("c_efficiency_vs_eta_Ethr%.1fMeV",Ethr),1400,1200);
    gPad -> SetGridx();
    gPad -> SetGridy();
    
    hPad = (TH1F*)( gPad->DrawFrame(etaMin,0.,etaMax,1.1) );
    hPad -> SetTitle(";|#eta|;efficiency");
    hPad -> Draw();

    latexLabel2 = new TLatex(0.18,0.90,Form("E_{thr} = %.1f MeV",Ethr));
    latexLabel2 -> SetNDC();
    latexLabel2 -> SetTextFont(82);
    latexLabel2 -> SetTextSize(0.02);
    latexLabel2 -> Draw("same");
    
    legend = new TLegend(0.60,0.16,0.83,0.16+(ptRanges.size()-1)*0.04);
    legend -> SetFillColor(0);
    legend -> SetFillStyle(1000);  
    legend -> SetTextFont(82);    
    legend -> SetTextSize(0.02);
    
    for(unsigned int jj = 0; jj < ptRanges.size()-1; ++jj)
    {
      float ptMin = ptRanges.at(jj);
      float ptMax = ptRanges.at(jj+1);

      e = (TEfficiency*)( inFile->Get(Form("p1_matchedRecHit_eff_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,Ethr)) );
      
      e -> SetMarkerSize(1.0);
      e -> SetMarkerColor(51+int(48/(ptRanges.size()-2))*jj);
      e -> SetLineColor(51+int(48/(ptRanges.size()-2))*jj);
      e -> Draw("same");
      
      legend -> AddEntry(e,Form("p_{T} #in [%.1f,%.1f] GeV",ptMin,ptMax),"PL");
      ++EthrIt;      
    }
    
    legend -> Draw("same");
    latexLabel -> Draw("same");
    
    c -> Print(Form("%s/c_%s_efficiency_vs_eta_Ethr%.1fMeV.png",plotDir.c_str(),label.c_str(),Ethr));
    c -> Print(Form("%s/c_%s_efficiency_vs_eta_Ethr%.1fMeV.pdf",plotDir.c_str(),label.c_str(),Ethr));
    delete c;
  }
  
  
  
  c = new TCanvas("c_totEnergy_vs_eta","c_totEnergy_vs_eta",1400,1200);
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(etaMin,0.,etaMax,12.) );
  hPad -> SetTitle(";|#eta|;#Sigma E_{recHit} [MeV]");
  hPad -> Draw();
  
  legend = new TLegend(0.60,0.16,0.83,0.16+EthrVals.size()*0.04);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);
  legend -> SetTextSize(0.02);
  
  EthrIt = 0;
  for(auto Ethr : EthrVals)
  {
    p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_totEnergy_vs_eta__pt00.8-%04.1f__Ethr%.1fMeV",ptRanges.at(ptRanges.size()-1),Ethr)) );
    
    p -> SetMarkerSize(1.0);
    p -> SetMarkerColor(51+int(48/(EthrVals.size()-1))*EthrIt);
    p -> SetLineColor(51+int(48/(EthrVals.size()-1))*EthrIt);
    p -> Draw("same");
    
    legend -> AddEntry(p,Form("E_{thr} = %.1f MeV",Ethr),"PL");
    ++EthrIt;
  }
  
  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  c -> Print(Form("%s/c_%s_totEnergy_vs_eta.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_totEnergy_vs_eta.pdf",plotDir.c_str(),label.c_str()));
  delete c;
  
  
  
  EthrIt = 0;
  for(auto Ethr : EthrVals)
  {
    c = new TCanvas(Form("c_totEnergy_vs_eta_Ethr%.1fMeV",Ethr),Form("c_totEnergy_vs_eta_Ethr%.1fMeV",Ethr),1400,1200);
    gPad -> SetGridx();
    gPad -> SetGridy();
    
    hPad = (TH1F*)( gPad->DrawFrame(etaMin,0.,etaMax,12.) );
    hPad -> SetTitle(";|#eta|;#Sigma E_{recHit} [MeV]");
    hPad -> Draw();
    
    latexLabel2 = new TLatex(0.18,0.90,Form("E_{thr} = %.1f MeV",Ethr));
    latexLabel2 -> SetNDC();
    latexLabel2 -> SetTextFont(82);
    latexLabel2 -> SetTextSize(0.02);
    latexLabel2 -> Draw("same");
    
    legend = new TLegend(0.60,0.16,0.83,0.16+(ptRanges.size()-1)*0.04);
    legend -> SetFillColor(0);
    legend -> SetFillStyle(1000);  
    legend -> SetTextFont(82);    
    legend -> SetTextSize(0.02);
    
    for(unsigned int jj = 0; jj < ptRanges.size()-1; ++jj)
    {
      float ptMin = ptRanges.at(jj);
      float ptMax = ptRanges.at(jj+1);
      
      p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_totEnergy_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,Ethr)) );
      
      p -> SetMarkerSize(1.0);
      p -> SetMarkerColor(51+int(48/(ptRanges.size()-2))*jj);
      p -> SetLineColor(51+int(48/(ptRanges.size()-2))*jj);
      p -> Draw("same");

      legend -> AddEntry(p,Form("p_{T} #in [%.1f,%.1f] GeV",ptMin,ptMax),"PL");
      ++EthrIt;
    }
    
    legend -> Draw("same");
    latexLabel -> Draw("same");
    
    c -> Print(Form("%s/c_%s_totEnergy_vs_eta_Ethr%.1fMeV.png",plotDir.c_str(),label.c_str(),Ethr));
    c -> Print(Form("%s/c_%s_totEnergy_vs_eta_Ethr%.1fMeV.pdf",plotDir.c_str(),label.c_str(),Ethr));
    delete c;
  }
  
  
  
  c = new TCanvas("c_n_vs_eta","c_n_vs_eta",1400,1200);
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(etaMin,0.,etaMax,3.) );
  hPad -> SetTitle(";|#eta|;N_{recHit}");
  hPad -> Draw();
  
  legend = new TLegend(0.60,0.90-EthrVals.size()*0.04,0.83,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);
  legend -> SetTextSize(0.02);
  
  EthrIt = 0;
  for(auto Ethr : EthrVals)
  {
    p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_n_vs_eta__pt00.8-%04.1f__Ethr%.1fMeV",ptRanges.at(ptRanges.size()-1),Ethr)) );
    
    p -> SetMarkerSize(1.0);
    p -> SetMarkerColor(51+int(48/(EthrVals.size()-1))*EthrIt);
    p -> SetLineColor(51+int(48/(EthrVals.size()-1))*EthrIt);
    p -> Draw("same");
    
    legend -> AddEntry(p,Form("E_{thr} = %.1f MeV",Ethr),"PL");
    ++EthrIt;
  }
  
  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  c -> Print(Form("%s/c_%s_n_vs_eta.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_n_vs_eta.pdf",plotDir.c_str(),label.c_str()));
  delete c;
  

  
  EthrIt = 0;
  for(auto Ethr : EthrVals)
  {
    c = new TCanvas(Form("c_n_vs_eta_Ethr%.1fMeV",Ethr),Form("c_n_vs_eta_Ethr%.1fMeV",Ethr),1400,1200);
    gPad -> SetGridx();
    gPad -> SetGridy();
    
    hPad = (TH1F*)( gPad->DrawFrame(etaMin,0.,etaMax,4.) );
    hPad -> SetTitle(";|#eta|;N_{recHit}");
    hPad -> Draw();
    
    latexLabel2 = new TLatex(0.18,0.90,Form("E_{thr} = %.1f MeV",Ethr));
    latexLabel2 -> SetNDC();
    latexLabel2 -> SetTextFont(82);
    latexLabel2 -> SetTextSize(0.02);
    latexLabel2 -> Draw("same");
    
    legend = new TLegend(0.60,0.90-(ptRanges.size()-1)*0.04,0.83,0.90);
    legend -> SetFillColor(0);
    legend -> SetFillStyle(1000);  
    legend -> SetTextFont(82);    
    legend -> SetTextSize(0.02);
    
    for(unsigned int jj = 0; jj < ptRanges.size()-1; ++jj)
    {
      float ptMin = ptRanges.at(jj);
      float ptMax = ptRanges.at(jj+1);
      
      p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_n_vs_eta__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,Ethr)) );
      
      p -> SetMarkerSize(1.0);
      p -> SetMarkerColor(51+int(48/(ptRanges.size()-2))*jj);
      p -> SetLineColor(51+int(48/(ptRanges.size()-2))*jj);
      p -> Draw("same");

      legend -> AddEntry(p,Form("p_{T} #in [%.1f,%.1f] GeV",ptMin,ptMax),"PL");
      ++EthrIt;
    }
    
    legend -> Draw("same");
    latexLabel -> Draw("same");
    
    c -> Print(Form("%s/c_%s_n_vs_eta_Ethr%.1fMeV.png",plotDir.c_str(),label.c_str(),Ethr));
    c -> Print(Form("%s/c_%s_n_vs_eta_Ethr%.1fMeV.pdf",plotDir.c_str(),label.c_str(),Ethr));
    delete c;
  }
  
  
  
  
  
  
  c = new TCanvas("c_efficiency_vs_phi","c_efficiency_vs_phi",1400,1200);
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(phiMin,0.,phiMax,1.1) );
  hPad -> SetTitle(";|#phi|;efficiency");
  hPad -> Draw();
  
  legend = new TLegend(0.60,0.16,0.83,0.16+EthrVals.size()*0.04);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);
  legend -> SetTextSize(0.02);
  
  EthrIt = 0;
  for(auto Ethr : EthrVals)
  {
    e = (TEfficiency*)( inFile->Get(Form("p1_matchedRecHit_eff_vs_phi__pt00.8-%04.1f__Ethr%.1fMeV",ptRanges.at(ptRanges.size()-1),Ethr)) );

    e -> SetMarkerSize(1.0);
    e -> SetMarkerColor(51+int(48/(EthrVals.size()-1))*EthrIt);
    e -> SetLineColor(51+int(48/(EthrVals.size()-1))*EthrIt);
    e -> Draw("same");

    legend -> AddEntry(e,Form("E_{thr} = %.1f MeV",Ethr),"PL");
    ++EthrIt;
  }
  
  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  c -> Print(Form("%s/c_%s_efficiency_vs_phi.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_efficiency_vs_phi.pdf",plotDir.c_str(),label.c_str()));
  delete c;


  EthrIt = 0;
  for(auto Ethr : EthrVals)
  {
    c = new TCanvas(Form("c_efficiency_vs_phi_Ethr%.1fMeV",Ethr),Form("c_efficiency_vs_phi_Ethr%.1fMeV",Ethr),1400,1200);
    gPad -> SetGridx();
    gPad -> SetGridy();
    
    hPad = (TH1F*)( gPad->DrawFrame(phiMin,0.,phiMax,1.1) );
    hPad -> SetTitle(";|#phi|;efficiency");
    hPad -> Draw();
    
    latexLabel2 = new TLatex(0.18,0.90,Form("E_{thr} = %.1f MeV",Ethr));
    latexLabel2 -> SetNDC();
    latexLabel2 -> SetTextFont(82);
    latexLabel2 -> SetTextSize(0.02);
    latexLabel2 -> Draw("same");
    
    legend = new TLegend(0.60,0.16,0.83,0.16+(ptRanges.size()-1)*0.04);
    legend -> SetFillColor(0);
    legend -> SetFillStyle(1000);  
    legend -> SetTextFont(82);    
    legend -> SetTextSize(0.02);    
    
    for(unsigned int jj = 0; jj < ptRanges.size()-1; ++jj)
    {
      float ptMin = ptRanges.at(jj);
      float ptMax = ptRanges.at(jj+1);

      e = (TEfficiency*)( inFile->Get(Form("p1_matchedRecHit_eff_vs_phi__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,Ethr)) );
      
      e -> SetMarkerSize(1.0);
      e -> SetMarkerColor(51+int(48/(ptRanges.size()-2))*jj);
      e -> SetLineColor(51+int(48/(ptRanges.size()-2))*jj);
      e -> Draw("same");

      legend -> AddEntry(e,Form("p_{T} #in [%.1f,%.1f] GeV",ptMin,ptMax),"PL");
      ++EthrIt;      
    }
    
    legend -> Draw("same");
    latexLabel -> Draw("same");
    
    c -> Print(Form("%s/c_%s_efficiency_vs_phi_Ethr%.1fMeV.png",plotDir.c_str(),label.c_str(),Ethr));
    c -> Print(Form("%s/c_%s_efficiency_vs_phi_Ethr%.1fMeV.pdf",plotDir.c_str(),label.c_str(),Ethr));
    delete c;
  }
  
  



  
  c = new TCanvas("c_totEnergy_vs_phi","c_totEnergy_vs_phi",1400,1200);
  gPad -> SetGridx();
  gPad -> SetGridy();
  
  hPad = (TH1F*)( gPad->DrawFrame(phiMin,0.,phiMax,12.) );
  hPad -> SetTitle(";|#phi|;#Sigma E_{recHit} [MeV]");
  hPad -> Draw();
  
  legend = new TLegend(0.60,0.16,0.83,0.16+EthrVals.size()*0.04);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(82);
  legend -> SetTextSize(0.02);
  
  EthrIt = 0;
  for(auto Ethr : EthrVals)
  {
    p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_totEnergy_vs_phi__pt00.8-%04.1f__Ethr%.1fMeV",ptRanges.at(ptRanges.size()-1),Ethr)) );
    
    p -> SetMarkerSize(1.0);
    p -> SetMarkerColor(51+int(48/(EthrVals.size()-1))*EthrIt);
    p -> SetLineColor(51+int(48/(EthrVals.size()-1))*EthrIt);
    p -> Draw("same");
    
    legend -> AddEntry(p,Form("E_{thr} = %.1f MeV",Ethr),"PL");
    ++EthrIt;
  }
  
  legend -> Draw("same");
  latexLabel -> Draw("same");
  
  c -> Print(Form("%s/c_%s_totEnergy_vs_phi.png",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s_totEnergy_vs_phi.pdf",plotDir.c_str(),label.c_str()));
  delete c;
  
  
  
  EthrIt = 0;
  for(auto Ethr : EthrVals)
  {
    c = new TCanvas(Form("c_totEnergy_vs_phi_Ethr%.1fMeV",Ethr),Form("c_totEnergy_vs_phi_Ethr%.1fMeV",Ethr),1400,1200);
    gPad -> SetGridx();
    gPad -> SetGridy();
    
    hPad = (TH1F*)( gPad->DrawFrame(phiMin,0.,phiMax,12.) );
    hPad -> SetTitle(";|#phi|;#Sigma E_{recHit} [MeV]");
    hPad -> Draw();
    
    latexLabel2 = new TLatex(0.18,0.90,Form("E_{thr} = %.1f MeV",Ethr));
    latexLabel2 -> SetNDC();
    latexLabel2 -> SetTextFont(82);
    latexLabel2 -> SetTextSize(0.02);
    latexLabel2 -> Draw("same");
    
    legend = new TLegend(0.60,0.16,0.83,0.16+(ptRanges.size()-1)*0.04);
    legend -> SetFillColor(0);
    legend -> SetFillStyle(1000);  
    legend -> SetTextFont(82);    
    legend -> SetTextSize(0.02);
    
    for(unsigned int jj = 0; jj < ptRanges.size()-1; ++jj)
    {
      float ptMin = ptRanges.at(jj);
      float ptMax = ptRanges.at(jj+1);
      
      p = (TProfile*)( inFile->Get(Form("p1_matchedRecHit_totEnergy_vs_phi__pt%04.1f-%04.1f__Ethr%.1fMeV",ptMin,ptMax,Ethr)) );
      
      p -> SetMarkerSize(1.0);
      p -> SetMarkerColor(51+int(48/(ptRanges.size()-2))*jj);
      p -> SetLineColor(51+int(48/(ptRanges.size()-2))*jj);
      p -> Draw("same");

      legend -> AddEntry(p,Form("p_{T} #in [%.1f,%.1f] GeV",ptMin,ptMax),"PL");
      ++EthrIt;
    }
    
    legend -> Draw("same");
    latexLabel -> Draw("same");
    
    c -> Print(Form("%s/c_%s_totEnergy_vs_phi_Ethr%.1fMeV.png",plotDir.c_str(),label.c_str(),Ethr));
    c -> Print(Form("%s/c_%s_totEnergy_vs_phi_Ethr%.1fMeV.pdf",plotDir.c_str(),label.c_str(),Ethr));
    delete c;
  }
  
  
  
  
  
  
}
