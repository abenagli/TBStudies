
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
#include "THStack.h"
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
#include "TLorentzVector.h"

#define PI 3.14159265359

float ptMin = 20.;



void DrawPlot(std::map<int,TH1F*> h1, std::map<int,TH1F*> h2,
              const std::string& lab1, const std::string& lab2,
              const std::string& label, const std::string& title, const std::string& text,
              const std::string plotDir = "./",
              const bool& logy = false,
              const float& norm1 = 1., const float& norm2 = 1.)
{
  TCanvas* c = new TCanvas();
  c -> cd();
  
  TPad* pad1 = new TPad("pad1","pad1",0,0.3,1,1.0);
  pad1->SetBottomMargin(0.03);
  pad1->SetGridx();
  pad1->Draw();
  pad1->cd();
  if( logy ) gPad -> SetLogy();
  
  THStack* hs1 = new THStack("hs1","");
  THStack* hs2 = new THStack("hs2","");

  for(int it = 0; it < 3; ++it)
  {
    h2[it] -> SetTitle(title.c_str());
    h2[it] -> SetLineColor(19-2*it);
    h2[it] -> SetFillColor(19-2*it);
    h2[it] -> SetLineWidth(2);
    h2[it]->GetXaxis()->SetTitleFont(43);
    h2[it]->GetYaxis()->SetTitleFont(43);
    h2[it]->GetXaxis()->SetTitleSize(50);
    h2[it]->GetYaxis()->SetTitleSize(22);
    h2[it]->GetXaxis()->SetTitleOffset(2.);
    h2[it]->GetYaxis()->SetTitleOffset(1.85);
    h2[it]->GetXaxis()->SetLabelFont(43);
    h2[it]->GetYaxis()->SetLabelFont(43);
    h2[it]->GetXaxis()->SetLabelSize(17);
    h2[it]->GetYaxis()->SetLabelSize(17);
    hs2->Add(h2[it]);

    h1[it] -> SetLineColor(kRed);
    h1[it] -> SetFillStyle(0);
    h1[it] -> SetLineWidth(2);
    hs1->Add(h1[it]);
  }

  TH1F* h_sum1 = (TH1F*)( h1[0]->Clone("h_sum1") );
  h_sum1 -> Add(h1[1]);
  h_sum1 -> Add(h1[2]);
  
  TH1F* h_sum2 = (TH1F*)( h2[0]->Clone("h_sum2") );
  h_sum2 -> Add(h2[1]);
  h_sum2 -> Add(h2[2]);
  

  TH1F* hist;
  TList* hs_histos1 = hs1->GetHists();
  TIter next1(hs_histos1);
  TList* hs_histos2 = hs2->GetHists();
  TIter next2(hs_histos2);
  if( norm1 != 1. )
  {
    while( (hist =(TH1F*)next1()) )
      hist -> Scale(1./norm1);
    h_sum1 -> Scale(1./norm1);
  }
  if( norm2 != 1. )
  {
    while( (hist =(TH1F*)next2()) )
      hist -> Scale(1./norm2);
    h_sum2 -> Scale(1./norm2);
  }
  
  
  float yMin = 999999.;
  float yMax = -999999.;
  for(int bin = 1; bin <= h_sum1->GetNbinsX(); ++bin)
  {
    float content1 = h_sum1->GetBinContent(bin);
    if( content1 < yMin && content1 > 0 ) yMin = content1;
    if( content1 > yMax ) yMax = content1;
  }
  for(int bin = 1; bin <= h_sum2->GetNbinsX(); ++bin)
  {
    float content2 = h_sum2->GetBinContent(bin);
    if( content2 < yMin && content2 > 0 ) yMin = content2;
    if( content2 > yMax ) yMax = content2;
  }

  if( !logy ) { h_sum1 -> SetMinimum(0.); hs1 -> SetMinimum(0.); }
  if( !logy ) { h_sum2 -> SetMinimum(0.); hs2 -> SetMinimum(0.); }
  if( logy )  { h_sum1 -> SetMinimum(yMin); hs1 -> SetMinimum(yMin); }
  if( logy )  { h_sum2 -> SetMinimum(yMin); hs2 -> SetMinimum(yMin); }
  h_sum1 -> SetMaximum(yMax); hs1 -> SetMaximum(yMax);
  h_sum2 -> SetMaximum(yMax); hs2 -> SetMaximum(yMax);

  hs2 -> SetTitle(title.c_str());
  hs2 -> Draw("hist");
  // hs1 -> Draw("hist,same");
  h_sum1 -> Draw("hist,same");
  
  TLegend* legend = new TLegend(0.60,0.90-0.08*5,0.99,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  legend -> AddEntry(h2[0],Form("%s - prompt jet",lab2.c_str()),"F");
  legend -> AddEntry(h2[1],Form("%s - fake jet",lab2.c_str()),"F");
  legend -> AddEntry(h2[2],Form("%s - undecided",lab2.c_str()),"F");
  legend -> AddEntry(h1[0],Form("%s",lab1.c_str()),"L");
  legend -> Draw("same");
  
  hs2->GetXaxis()->SetLabelSize(0.);
  hs2->GetXaxis()->SetTitleSize(0.);

  TLatex* latex = new TLatex(0.80,0.96,text.c_str());
  latex -> SetNDC();
  latex -> SetTextFont(42);
  latex -> SetTextSize(0.05);
  latex ->SetTextAlign(31);
  latex -> Draw("same");
  

  c -> cd();
  
  TPad* pad2 = new TPad("pad2","pad2", 0,0.,1,0.3);
  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.4);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  
  TH1F* h_ratio  = (TH1F*)( h_sum2->Clone(Form("%s_ratio",h_sum2->GetName())) );
  h_ratio -> Divide(h_sum1);
  h_ratio -> GetYaxis() -> SetRangeUser(0.,3.);
  h_ratio -> SetMarkerColor(kBlack);
  h_ratio -> SetLineColor(kBlack);
  h_ratio -> SetMarkerSize(0.7);
  h_ratio -> SetLineWidth(1);
  h_ratio -> SetTitle(Form(";%s;%s / %s",h_sum2->GetXaxis()->GetTitle(),lab2.c_str(),lab1.c_str()));
  h_ratio -> Draw();

  h_ratio->GetXaxis()->SetTitleFont(43);
  h_ratio->GetYaxis()->SetTitleFont(43);
  h_ratio->GetXaxis()->SetTitleSize(30);
  h_ratio->GetYaxis()->SetTitleSize(18);
  h_ratio->GetXaxis()->SetTitleOffset(3.);
  h_ratio->GetYaxis()->SetTitleOffset(1.85);
  h_ratio->GetXaxis()->SetLabelFont(43);
  h_ratio->GetYaxis()->SetLabelFont(43);
  h_ratio->GetXaxis()->SetLabelSize(15);
  h_ratio->GetYaxis()->SetLabelSize(15);
  h_ratio->GetYaxis()->SetNdivisions(110);
  TF1* f_line1 = new TF1("f_line1","1.",-999999.,999999.);
  f_line1 -> SetLineColor(kBlack);
  f_line1 -> SetLineStyle(2);
  f_line1 -> Draw("same");
  
  c -> Print(Form("%s/c_%s.pdf",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s.png",plotDir.c_str(),label.c_str()));
  delete c;
}



void DrawPlot(TH1F* h1, TH1F* h2,
              const std::string& lab1, const std::string& lab2,
              const std::string& label, const std::string& title, const std::string& text,
              const std::string plotDir = "./",
              const bool& logy = false,
              const float& norm1 = 1., const float& norm2 = 1.)
{
  TCanvas* c = new TCanvas();
  c -> cd();
  
  h1 -> SetLineColor(kRed);
  h1 -> SetLineWidth(2);
  h2 -> SetFillColor(19);
  
  TPad* pad1 = new TPad("pad1","pad1",0,0.3,1,1.0);
  pad1->SetBottomMargin(0.03);
  pad1->SetGridx();
  pad1->Draw();
  pad1->cd();
  if( logy ) gPad -> SetLogy();
  
  float yMin = 999999.;
  float yMax = -999999.;
  for(int bin = 1; bin <= h1->GetNbinsX(); ++bin)
  {
    float content1 = h1->GetBinContent(bin);
    if( content1 < yMin && content1 > 0 ) yMin = content1;
    if( content1 > yMax ) yMax = content1;
  }
  for(int bin = 1; bin <= h2->GetNbinsX(); ++bin)
  {
    float content2 = h2->GetBinContent(bin);
    if( content2 < yMin && content2 > 0 ) yMin = content2;
    if( content2 > yMax ) yMax = content2;
  }

  if( !logy ) { h1 -> SetMinimum(0.); }
  if( !logy ) { h2 -> SetMinimum(0.); }
  if( logy )  { h1 -> SetMinimum(yMin); }
  if( logy )  { h2 -> SetMinimum(yMin); }
  h1 -> SetMaximum(yMax);
  h2 -> SetMaximum(yMax);

  h2 -> SetTitle(title.c_str());
  h2 -> Draw("hist");
  h1 -> Draw("hist,same");
  
  TLegend* legend = new TLegend(0.60,0.90-0.08*5,0.99,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);  
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  legend -> AddEntry(h2,Form("%s - prompt jet",lab2.c_str()),"F");
  legend -> AddEntry(h1,Form("%s - prompt jet",lab1.c_str()),"L");
  legend -> Draw("same");
  
  h2->GetXaxis()->SetLabelSize(0.);
  h2->GetXaxis()->SetTitleSize(0.);

  TLatex* latex = new TLatex(0.80,0.96,text.c_str());
  latex -> SetNDC();
  latex -> SetTextFont(42);
  latex -> SetTextSize(0.05);
  latex ->SetTextAlign(31);
  latex -> Draw("same");
  

  c -> cd();
  
  TPad* pad2 = new TPad("pad2","pad2", 0,0.,1,0.3);
  pad2->SetTopMargin(0.03);
  pad2->SetBottomMargin(0.4);
  pad2->SetGridx();
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  
  TH1F* h_ratio  = (TH1F*)( h2->Clone(Form("%s_ratio",h2->GetName())) );
  h_ratio -> Divide(h1);
  h_ratio -> GetYaxis() -> SetRangeUser(0.,3.);
  h_ratio -> SetMarkerColor(kBlack);
  h_ratio -> SetLineColor(kBlack);
  h_ratio -> SetMarkerSize(0.7);
  h_ratio -> SetLineWidth(1);
  h_ratio -> SetTitle(Form(";%s;%s/%s",h2->GetXaxis()->GetTitle(),lab2.c_str(),lab1.c_str()));
  h_ratio -> Draw();

  h_ratio->GetXaxis()->SetTitleFont(43);
  h_ratio->GetYaxis()->SetTitleFont(43);
  h_ratio->GetXaxis()->SetTitleSize(30);
  h_ratio->GetYaxis()->SetTitleSize(18);
  h_ratio->GetXaxis()->SetTitleOffset(3.);
  h_ratio->GetYaxis()->SetTitleOffset(1.85);
  h_ratio->GetXaxis()->SetLabelFont(43);
  h_ratio->GetYaxis()->SetLabelFont(43);
  h_ratio->GetXaxis()->SetLabelSize(15);
  h_ratio->GetYaxis()->SetLabelSize(15);
  h_ratio->GetYaxis()->SetNdivisions(110);
  TF1* f_line1 = new TF1("f_line1","1.",-999999.,999999.);
  f_line1 -> SetLineColor(kBlack);
  f_line1 -> SetLineStyle(2);
  f_line1 -> Draw("same");
  
  c -> Print(Form("%s/c_%s.pdf",plotDir.c_str(),label.c_str()));
  c -> Print(Form("%s/c_%s.png",plotDir.c_str(),label.c_str()));
  delete c;
}






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
  
  
  
  std::vector<float> ptRanges;
  ptRanges.push_back(20.);
  ptRanges.push_back(30.);
  ptRanges.push_back(35.);
  ptRanges.push_back(40.);
  ptRanges.push_back(45.);
  ptRanges.push_back(50.);
  ptRanges.push_back(60.);
  ptRanges.push_back(70.);
  ptRanges.push_back(80.);
  ptRanges.push_back(100.);
  ptRanges.push_back(150.);
  ptRanges.push_back(200.);
  int nBins_pt = ptRanges.size() - 1;
  double* xAxis_pt = new double[ptRanges.size()];
  for(unsigned int jj = 0; jj < ptRanges.size(); ++jj)
    xAxis_pt[jj] = ptRanges.at(jj);
  
  
  
  //-----------------
  // define variables
  int genLeptons_n;
  std::vector<int>*   genLeptons_pdgId = new std::vector<int>;
  std::vector<float>* genLeptons_pt    = new std::vector<float>;
  std::vector<float>* genLeptons_eta   = new std::vector<float>;
  std::vector<float>* genLeptons_phi   = new std::vector<float>;
  std::vector<float>* genLeptons_energy= new std::vector<float>;
  std::vector<float>* genLeptons_vtx_x = new std::vector<float>;
  std::vector<float>* genLeptons_vtx_y = new std::vector<float>;
  std::vector<float>* genLeptons_vtx_z = new std::vector<float>;
  
  float genVtx_x;
  float genVtx_y;
  float genVtx_z;
  float genVtx_t;
  
  int vtxs_n;
  std::vector<float>* vtxs_x = new std::vector<float>;
  std::vector<float>* vtxs_y = new std::vector<float>;
  std::vector<float>* vtxs_z = new std::vector<float>;
  std::vector<float>* vtxs_t = new std::vector<float>;
  
  int jets_n;
  std::vector<float>* jets_pt     = new std::vector<float>;
  std::vector<float>* jets_eta    = new std::vector<float>;
  std::vector<float>* jets_phi    = new std::vector<float>;
  std::vector<float>* jets_mass   = new std::vector<float>;
  std::vector<int>*   jets_isPU   = new std::vector<int>;
  
  std::vector<float>* jets_matchedGenJet_pt = new std::vector<float>;
  std::vector<float>* jets_matchedGenJet_eta = new std::vector<float>;
  std::vector<float>* jets_matchedGenJet_phi = new std::vector<float>;
  std::vector<float>* jets_matchedGenJet_energy = new std::vector<float>;
  
  std::vector<float>* jets_myMatchedGenJet_pt = new std::vector<float>;
  std::vector<float>* jets_myMatchedGenJet_eta = new std::vector<float>;
  std::vector<float>* jets_myMatchedGenJet_phi = new std::vector<float>;
  std::vector<float>* jets_myMatchedGenJet_energy = new std::vector<float>;
  std::vector<std::vector<float> >* jets_myMatchedGenJet_constituents_pt = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* jets_myMatchedGenJet_constituents_eta = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* jets_myMatchedGenJet_constituents_phi = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* jets_myMatchedGenJet_constituents_mass = new std::vector<std::vector<float> >;
  std::vector<std::vector<int> >* jets_myMatchedGenJet_constituents_charge = new std::vector<std::vector<int> >;
  std::vector<std::vector<int> >* jets_myMatchedGenJet_constituents_pdgId = new std::vector<std::vector<int> >;
  
  std::vector<std::vector<float> >* jets_constituents_pt    = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* jets_constituents_eta   = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* jets_constituents_phi   = new std::vector<std::vector<float> >;
  std::vector<std::vector<float> >* jets_constituents_mass  = new std::vector<std::vector<float> >;
  std::vector<std::vector<int> >*   jets_constituents_pdgId = new std::vector<std::vector<int> >;
  std::vector<std::vector<float> >* jets_constituents_puppiWeight = new std::vector<std::vector<float> >;
  
  //----------
  // get trees
  std::vector<std::string> sampleLabels;
  std::map<std::string,std::string> labels;
  std::vector<std::string> uniqueLabels;
  std::vector<float> xsecs;
  std::map<std::string,TFile*> inFiles_bkg;
  std::map<std::string,TChain*> trees_bkg;
  
  std::vector<std::string> inputFiles = opts.GetOpt<std::vector<std::string> >("Input.inputFiles");
  for(unsigned int fileIt = 0; fileIt < inputFiles.size()/4; ++fileIt)
  {
    std::string sampleLabel = inputFiles.at(0+fileIt*4);
    sampleLabels.push_back(sampleLabel);
    
    std::string label = inputFiles.at(1+fileIt*4);
    labels[sampleLabel] = label;
    uniqueLabels.push_back(label);
    
    float xsec = atof( inputFiles.at(3+fileIt*4).c_str() );
    xsecs.push_back(xsec);
    
    trees_bkg[sampleLabel] = new TChain("FTLDumpJets/jet_tree");
    
    trees_bkg[sampleLabel] -> SetBranchAddress("genLeptons_n",    &genLeptons_n);
    trees_bkg[sampleLabel] -> SetBranchAddress("genLeptons_pt",   &genLeptons_pt);
    trees_bkg[sampleLabel] -> SetBranchAddress("genLeptons_eta",  &genLeptons_eta);
    trees_bkg[sampleLabel] -> SetBranchAddress("genLeptons_phi",  &genLeptons_phi);
    trees_bkg[sampleLabel] -> SetBranchAddress("genLeptons_energy",&genLeptons_energy);
    trees_bkg[sampleLabel] -> SetBranchAddress("genLeptons_pdgId",&genLeptons_pdgId);
    trees_bkg[sampleLabel] -> SetBranchAddress("genLeptons_vtx_x",&genLeptons_vtx_x);
    trees_bkg[sampleLabel] -> SetBranchAddress("genLeptons_vtx_y",&genLeptons_vtx_y);
    trees_bkg[sampleLabel] -> SetBranchAddress("genLeptons_vtx_z",&genLeptons_vtx_z);
    
    trees_bkg[sampleLabel] -> SetBranchAddress("genVtx_x",&genVtx_x);
    trees_bkg[sampleLabel] -> SetBranchAddress("genVtx_y",&genVtx_y);
    trees_bkg[sampleLabel] -> SetBranchAddress("genVtx_z",&genVtx_z);
    trees_bkg[sampleLabel] -> SetBranchAddress("genVtx_t",&genVtx_t);
    
    trees_bkg[sampleLabel] -> SetBranchAddress("vtxs_n",&vtxs_n);
    trees_bkg[sampleLabel] -> SetBranchAddress("vtxs_x",&vtxs_x);
    trees_bkg[sampleLabel] -> SetBranchAddress("vtxs_y",&vtxs_y);
    trees_bkg[sampleLabel] -> SetBranchAddress("vtxs_z",&vtxs_z);
    trees_bkg[sampleLabel] -> SetBranchAddress("vtxs_t",&vtxs_t);
    
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_n",     &jets_n);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_pt",    &jets_pt);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_eta",   &jets_eta);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_phi",   &jets_phi);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_mass",  &jets_mass);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_isPU",  &jets_isPU);
    
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_matchedGenJet_pt", &jets_matchedGenJet_pt);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_matchedGenJet_eta", &jets_matchedGenJet_eta);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_matchedGenJet_phi", &jets_matchedGenJet_phi);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_matchedGenJet_energy", &jets_matchedGenJet_energy);
    
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_myMatchedGenJet_pt", &jets_myMatchedGenJet_pt);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_myMatchedGenJet_eta", &jets_myMatchedGenJet_eta);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_myMatchedGenJet_phi", &jets_myMatchedGenJet_phi);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_myMatchedGenJet_energy", &jets_myMatchedGenJet_energy);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_myMatchedGenJet_constituents_pt", &jets_myMatchedGenJet_constituents_pt);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_myMatchedGenJet_constituents_eta", &jets_myMatchedGenJet_constituents_eta);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_myMatchedGenJet_constituents_phi", &jets_myMatchedGenJet_constituents_phi);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_myMatchedGenJet_constituents_mass", &jets_myMatchedGenJet_constituents_mass);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_myMatchedGenJet_constituents_charge", &jets_myMatchedGenJet_constituents_charge);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_myMatchedGenJet_constituents_pdgId", &jets_myMatchedGenJet_constituents_pdgId);
    
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_constituents_pt",    &jets_constituents_pt);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_constituents_eta",   &jets_constituents_eta);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_constituents_phi",   &jets_constituents_phi);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_constituents_mass",  &jets_constituents_mass);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_constituents_pdgId", &jets_constituents_pdgId);
    trees_bkg[sampleLabel] -> SetBranchAddress("jets_constituents_puppiWeight", &jets_constituents_puppiWeight);
    
    trees_bkg[sampleLabel] -> Add( Form("%s/*.root",inputFiles.at(2+fileIt*4).c_str()) );
  }
  
  std::vector<std::string>::iterator iter;
  iter = std::unique(uniqueLabels.begin(),uniqueLabels.end());
  uniqueLabels.resize( std::distance(uniqueLabels.begin(),iter) );  
  
  
  
  //------------------
  // define histograms
  TFile* outFile = TFile::Open("studyJets.root","RECREATE");
  
  std::map<std::string,std::map<int,TH1F*> > h_nJets;
  std::map<std::string,std::map<int,TH1F*> > h_pt;
  std::map<std::string,std::map<int,TH1F*> > h_eta;
  
  std::map<std::string,TH1F*> h_ERes_all;
  std::map<std::string,TH1F*> h_ERes_charged_all;
  std::map<std::string,TH1F*> h_ERes_em_all;
  std::map<std::string,TH1F*> h_ERes_neutral_all;
  std::map<std::string,TH1F*> h_ERes_lepton_all;
  std::map<std::string,std::map<std::pair<float,float>,TH1F*> > h_ERes_pt;
  std::map<std::string,std::map<std::pair<float,float>,TH1F*> > h_ERes_charged_pt;
  std::map<std::string,std::map<std::pair<float,float>,TH1F*> > h_ERes_em_pt;
  std::map<std::string,std::map<std::pair<float,float>,TH1F*> > h_ERes_neutral_pt;
  std::map<std::string,std::map<std::pair<float,float>,TH1F*> > h_ERes_lepton_pt;
  
  std::map<std::string,TProfile*> p_chargedPt_gen_prompt;
  std::map<std::string,TProfile*> p_emPt_gen_prompt;
  std::map<std::string,TProfile*> p_neutralPt_gen_prompt;
  std::map<std::string,TProfile*> p_leptonPt_gen_prompt;
  std::map<std::string,TProfile*> p_chargedPt_raw_prompt;
  std::map<std::string,TProfile*> p_emPt_raw_prompt;
  std::map<std::string,TProfile*> p_neutralPt_raw_prompt;
  std::map<std::string,TProfile*> p_leptonPt_raw_prompt;
  std::map<std::string,TProfile*> p_chargedPt_weighted_prompt;
  std::map<std::string,TProfile*> p_emPt_weighted_prompt;
  std::map<std::string,TProfile*> p_neutralPt_weighted_prompt;
  std::map<std::string,TProfile*> p_leptonPt_weighted_prompt;
  
  std::map<std::string,TProfile*> p_chargedPtFrac_gen_prompt;
  std::map<std::string,TProfile*> p_emPtFrac_gen_prompt;
  std::map<std::string,TProfile*> p_neutralPtFrac_gen_prompt;
  std::map<std::string,TProfile*> p_leptonPtFrac_gen_prompt;
  std::map<std::string,TProfile*> p_chargedPtFrac_raw_prompt;
  std::map<std::string,TProfile*> p_emPtFrac_raw_prompt;
  std::map<std::string,TProfile*> p_neutralPtFrac_raw_prompt;
  std::map<std::string,TProfile*> p_leptonPtFrac_raw_prompt;
  std::map<std::string,TProfile*> p_chargedPtFrac_weighted_prompt;
  std::map<std::string,TProfile*> p_emPtFrac_weighted_prompt;
  std::map<std::string,TProfile*> p_neutralPtFrac_weighted_prompt;
  std::map<std::string,TProfile*> p_leptonPtFrac_weighted_prompt;
  
  for(auto uniqueLabel : uniqueLabels)
  {
    for(int it = 0; it < 4; ++it)
    {
      std::string PULabel = "prompt";
      if( it == 1 ) PULabel = "fake";
      if( it == 2 ) PULabel = "undecided";
      if( it == 3 ) PULabel = "all";
      
      h_nJets[uniqueLabel][it] = new TH1F(Form("h_nJets_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",10,-0.5,9.5);
      h_nJets[uniqueLabel][it] -> Sumw2();
      
      h_pt[uniqueLabel][it] = new TH1F(Form("h_pt_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
      h_pt[uniqueLabel][it] -> Sumw2();
      
      h_eta[uniqueLabel][it] = new TH1F(Form("h_eta_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",50,-5.,5.);
      h_eta[uniqueLabel][it] -> Sumw2();
      
      if(it == 0 )
      {
        h_ERes_all[uniqueLabel] = new TH1F(Form("h_ERes_all_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",250,0.,5.);
        h_ERes_all[uniqueLabel] -> Sumw2();
        h_ERes_charged_all[uniqueLabel] = new TH1F(Form("h_ERes_charged_all_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",250,0.,2.);
        h_ERes_charged_all[uniqueLabel] -> Sumw2();
        h_ERes_em_all[uniqueLabel] = new TH1F(Form("h_ERes_em_all_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",250,0.,2.);
        h_ERes_em_all[uniqueLabel] -> Sumw2();
        h_ERes_neutral_all[uniqueLabel] = new TH1F(Form("h_ERes_neutral_all_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",250,0.,2.);
        h_ERes_neutral_all[uniqueLabel] -> Sumw2();
        h_ERes_lepton_all[uniqueLabel] = new TH1F(Form("h_ERes_lepton_all_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",250,0.,2.);
        h_ERes_lepton_all[uniqueLabel] -> Sumw2();
        
        for(unsigned int ptIt = 0; ptIt < ptRanges.size()-1; ++ptIt)
        {
          float ptRangeMin = ptRanges.at(ptIt);
          float ptRangeMax = ptRanges.at(ptIt+1);
          
          h_ERes_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)] = new TH1F(Form("h_ERes_pt%.0f-%.0f__%s_%s",ptRangeMin,ptRangeMax,uniqueLabel.c_str(),PULabel.c_str()),"",250,0.,5.);
          h_ERes_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)] -> Sumw2();
          h_ERes_charged_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)] = new TH1F(Form("h_ERes_charged_pt%.0f-%.0f__%s_%s",ptRangeMin,ptRangeMax,uniqueLabel.c_str(),PULabel.c_str()),"",250,0.,5.);
          h_ERes_charged_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)] -> Sumw2();
          h_ERes_em_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)] = new TH1F(Form("h_ERes_em_pt%.0f-%.0f__%s_%s",ptRangeMin,ptRangeMax,uniqueLabel.c_str(),PULabel.c_str()),"",250,0.,5.);
          h_ERes_em_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)] -> Sumw2();
          h_ERes_neutral_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)] = new TH1F(Form("h_ERes_neutral_pt%.0f-%.0f__%s_%s",ptRangeMin,ptRangeMax,uniqueLabel.c_str(),PULabel.c_str()),"",250,0.,5.);
          h_ERes_neutral_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)] -> Sumw2();
          h_ERes_lepton_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)] = new TH1F(Form("h_ERes_lepton_pt%.0f-%.0f__%s_%s",ptRangeMin,ptRangeMax,uniqueLabel.c_str(),PULabel.c_str()),"",250,0.,5.);
          h_ERes_lepton_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)] -> Sumw2();
        }
        
        p_chargedPt_gen_prompt[uniqueLabel] = new TProfile(Form("p_chargedPt_gen_%s_%s",      uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_emPt_gen_prompt[uniqueLabel] = new TProfile(Form("p_emPt_gen_%s_%s",      uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_neutralPt_gen_prompt[uniqueLabel] = new TProfile(Form("p_neutralPt_gen_%s_%s",      uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_leptonPt_gen_prompt[uniqueLabel] = new TProfile(Form("p_leptonPt_gen_%s_%s",      uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_chargedPt_raw_prompt[uniqueLabel] = new TProfile(Form("p_chargedPt_raw_%s_%s",      uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_emPt_raw_prompt[uniqueLabel] = new TProfile(Form("p_emPt_raw_%s_%s",      uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_neutralPt_raw_prompt[uniqueLabel] = new TProfile(Form("p_neutralPt_raw_%s_%s",      uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_leptonPt_raw_prompt[uniqueLabel] = new TProfile(Form("p_leptonPt_raw_%s_%s",      uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_chargedPt_weighted_prompt[uniqueLabel] = new TProfile(Form("p_chargedPt_weighted_%s_%s",      uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_emPt_weighted_prompt[uniqueLabel] = new TProfile(Form("p_emPt_weighted_%s_%s",      uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_neutralPt_weighted_prompt[uniqueLabel] = new TProfile(Form("p_neutralPt_weighted_%s_%s",      uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_leptonPt_weighted_prompt[uniqueLabel] = new TProfile(Form("p_leptonPt_weighted_%s_%s",      uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        
        p_chargedPtFrac_gen_prompt[uniqueLabel] = new TProfile(Form("p_chargedPtFrac_gen_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_emPtFrac_gen_prompt[uniqueLabel] = new TProfile(Form("p_emPtFrac_gen_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_neutralPtFrac_gen_prompt[uniqueLabel] = new TProfile(Form("p_neutralPtFrac_gen_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_leptonPtFrac_gen_prompt[uniqueLabel] = new TProfile(Form("p_leptonPtFrac_gen_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_chargedPtFrac_raw_prompt[uniqueLabel] = new TProfile(Form("p_chargedPtFrac_raw_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_emPtFrac_raw_prompt[uniqueLabel] = new TProfile(Form("p_emPtFrac_raw_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_neutralPtFrac_raw_prompt[uniqueLabel] = new TProfile(Form("p_neutralPtFrac_raw_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_leptonPtFrac_raw_prompt[uniqueLabel] = new TProfile(Form("p_leptonPtFrac_raw_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_chargedPtFrac_weighted_prompt[uniqueLabel] = new TProfile(Form("p_chargedPtFrac_weighted_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_emPtFrac_weighted_prompt[uniqueLabel] = new TProfile(Form("p_emPtFrac_weighted_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_neutralPtFrac_weighted_prompt[uniqueLabel] = new TProfile(Form("p_neutralPtFrac_weighted_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
        p_leptonPtFrac_weighted_prompt[uniqueLabel] = new TProfile(Form("p_leptonPtFrac_weighted_%s_%s",uniqueLabel.c_str(),PULabel.c_str()),"",nBins_pt,xAxis_pt);
      }
    }
  }
  
  
  
  //-----------------
  // loop over events

  float etaMin = opts.GetOpt<float>("Input.etaMin");
  float etaMax = opts.GetOpt<float>("Input.etaMax");
  
  // NoPU / PU200 samples
  for(unsigned int sampleLabelIt = 0; sampleLabelIt < sampleLabels.size(); ++sampleLabelIt)
  {
    std::string sampleLabel = sampleLabels.at(sampleLabelIt);
    std::string label = labels[sampleLabel];
     
    int nEntries = trees_bkg[sampleLabel]->GetEntries();
    // nEntries = 10000;
    float weight = xsecs.at(sampleLabelIt) / nEntries * 1.;
    
    std::cout << ">>> processing " << sampleLabel << " (" << sampleLabelIt << "/" << sampleLabels.size() << ")   -   " << label << "   -   " << nEntries << " entries" << std::endl;
    for(int entry = 0; entry < nEntries; ++entry)
    {
      if( entry%10000 == 0) std::cout << ">>> " << sampleLabel << "   reading entry " << entry << " / " << nEntries << "\r" << std::flush;
      trees_bkg[sampleLabel] -> GetEntry(entry);
      
      if( genLeptons_n != 2 ) continue;
      TLorentzVector lep1;
      lep1.SetPtEtaPhiE(genLeptons_pt->at(0),genLeptons_eta->at(0),genLeptons_phi->at(0),genLeptons_energy->at(0));
      TLorentzVector lep2;
      lep2.SetPtEtaPhiE(genLeptons_pt->at(1),genLeptons_eta->at(1),genLeptons_phi->at(1),genLeptons_energy->at(1));
      if( fabs((lep1+lep2).M()-91.1876) > 15. ) continue;
      if( (fabs(genLeptons_eta->at(0)) > 2.4) || (fabs(genLeptons_eta->at(1)) > 2.4) ) continue;
      
      if( fabs(vtxs_z->at(0)-genVtx_z) > 0.1 ) continue;
      // // if( fabs(vtxs_z->at(0)-genVtx_z) < 0.1 ) continue; // reversed selection
      // if( (vtxs_t->at(0) != 0) && ( fabs(vtxs_t->at(0)-genVtx_t) > 0.090) ) continue;
      
      int nJets_prompt = 0;
      int nJets_fake = 0;
      int nJets_und = 0;
      int nJets_all = 0;
      int itJet = 0;
      
      for(unsigned int jetIt = 0; jetIt < jets_pt->size(); ++jetIt)
      {
        if( jets_pt->at(jetIt) < ptMin ) continue;
        if( jets_pt->at(jetIt) > ptRanges.at(ptRanges.size()-1) ) continue;
        if( fabs(jets_eta->at(jetIt))  < etaMin ) continue;
        if( fabs(jets_eta->at(jetIt)) >= etaMax ) continue;
        
        TLorentzVector jet_p4;
        jet_p4.SetPtEtaPhiM(jets_pt->at(jetIt),jets_eta->at(jetIt),jets_phi->at(jetIt),jets_mass->at(jetIt));
        
        h_pt[label][3] -> Fill( jets_pt->at(jetIt),weight );
        h_eta[label][3] -> Fill( jets_eta->at(jetIt),weight );
        ++nJets_all;
          
        if( jets_isPU->at(jetIt) == 0 )
        {
          h_pt[label][0] -> Fill( jets_pt->at(jetIt),weight );
          h_eta[label][0] -> Fill( jets_eta->at(jetIt),weight );
          ++nJets_prompt;
          
          if( jets_myMatchedGenJet_pt->at(jetIt) > 0. && jets_matchedGenJet_pt->at(jetIt) > 0. )
          {
            TLorentzVector genJet_p4;
            genJet_p4.SetPtEtaPhiE(jets_matchedGenJet_pt->at(jetIt),jets_matchedGenJet_eta->at(jetIt),jets_matchedGenJet_phi->at(jetIt),jets_matchedGenJet_energy->at(jetIt));
            TLorentzVector myGenJet_p4;
            myGenJet_p4.SetPtEtaPhiE(jets_myMatchedGenJet_pt->at(jetIt),jets_myMatchedGenJet_eta->at(jetIt),jets_myMatchedGenJet_phi->at(jetIt),jets_myMatchedGenJet_energy->at(jetIt));
            
            if( myGenJet_p4.Pt() != genJet_p4.Pt() ) continue;
            
            
            float ptRangeMin = ptRanges.at(0);
            float ptRangeMax = ptRanges.at(1);
            for(unsigned int ptIt = 0; ptIt < ptRanges.size()-1; ++ptIt)
              if( genJet_p4.Pt() >= ptRanges.at(ptIt) &&
                  genJet_p4.Pt() < ptRanges.at(ptIt+1) )
              {
                ptRangeMin = ptRanges.at(ptIt);
                ptRangeMax = ptRanges.at(ptIt+1);
              }
            
            h_ERes_all[label] -> Fill( jet_p4.E()/genJet_p4.E(),weight );
            h_ERes_pt[label][std::make_pair(ptRangeMin,ptRangeMax)] -> Fill( jet_p4.E()/genJet_p4.E(),weight );
            
            float chargedE_gen = 0.;
            float emE_gen = 0.;
            float neutralE_gen = 0.;
            float leptonE_gen = 0.;
            float chargedPt_gen = 0.;
            float emPt_gen = 0.;
            float neutralPt_gen = 0.;
            float leptonPt_gen = 0.;
            for(unsigned int jj = 0; jj < jets_myMatchedGenJet_constituents_pt->at(jetIt).size(); ++jj)
            {
              TLorentzVector jetConst_gen_p4;
              jetConst_gen_p4.SetPtEtaPhiM(jets_myMatchedGenJet_constituents_pt->at(jetIt)[jj],
                                           jets_myMatchedGenJet_constituents_eta->at(jetIt)[jj],
                                           jets_myMatchedGenJet_constituents_phi->at(jetIt)[jj],
                                           jets_myMatchedGenJet_constituents_mass->at(jetIt)[jj]);
              
              if( fabs(jets_myMatchedGenJet_constituents_pdgId->at(jetIt)[jj]) == 22 )
                emE_gen += jetConst_gen_p4.E();
              else if( fabs(jets_myMatchedGenJet_constituents_pdgId->at(jetIt)[jj]) == 11 ||
                       fabs(jets_myMatchedGenJet_constituents_pdgId->at(jetIt)[jj]) == 13 )
                leptonE_gen += jetConst_gen_p4.E();
              else if( fabs(jets_myMatchedGenJet_constituents_charge->at(jetIt)[jj]) > 0 )
                chargedE_gen += jetConst_gen_p4.E();
              else
                neutralE_gen += jetConst_gen_p4.E();
              
              if( fabs(jets_myMatchedGenJet_constituents_pdgId->at(jetIt)[jj]) == 22 )
                emPt_gen += jetConst_gen_p4.Pt();
              else if( fabs(jets_myMatchedGenJet_constituents_pdgId->at(jetIt)[jj]) == 11 ||
                       fabs(jets_myMatchedGenJet_constituents_pdgId->at(jetIt)[jj]) == 13 )
                leptonPt_gen += jetConst_gen_p4.Pt();
              else if( fabs(jets_myMatchedGenJet_constituents_charge->at(jetIt)[jj]) > 0 )
                chargedPt_gen += jetConst_gen_p4.Pt();
              else
                neutralPt_gen += jetConst_gen_p4.Pt();
            }
            
            float chargedE_raw = 0.;
            float emE_raw = 0.;
            float neutralE_raw = 0.;
            float leptonE_raw = 0.;
            float chargedE_weighted = 0.;
            float emE_weighted = 0.;
            float neutralE_weighted = 0.;
            float leptonE_weighted = 0.;
            
            float chargedPt_raw = 0.;
            float emPt_raw = 0.;
            float neutralPt_raw = 0.;
            float leptonPt_raw = 0.;
            float chargedPt_weighted = 0.;
            float emPt_weighted = 0.;
            float neutralPt_weighted = 0.;
            float leptonPt_weighted = 0.;
            for(unsigned int jj = 0; jj < jets_constituents_pt->at(jetIt).size(); ++jj)
            {
              TLorentzVector jetConst_raw_p4;
              jetConst_raw_p4.SetPtEtaPhiM(jets_constituents_pt->at(jetIt)[jj],
                                           jets_constituents_eta->at(jetIt)[jj],
                                           jets_constituents_phi->at(jetIt)[jj],
                                           jets_constituents_mass->at(jetIt)[jj]);
              TLorentzVector jetConst_weighted_p4;
              jetConst_weighted_p4.SetPtEtaPhiM(jets_constituents_puppiWeight->at(jetIt)[jj]*jets_constituents_pt->at(jetIt)[jj],
                                                jets_constituents_puppiWeight->at(jetIt)[jj]*jets_constituents_eta->at(jetIt)[jj],
                                                jets_constituents_puppiWeight->at(jetIt)[jj]*jets_constituents_phi->at(jetIt)[jj],
                                                jets_constituents_puppiWeight->at(jetIt)[jj]*jets_constituents_mass->at(jetIt)[jj]);
              
              if(      fabs(jets_constituents_pdgId->at(jetIt)[jj]) == 211 ) chargedE_raw += jetConst_raw_p4.E();
              else if( fabs(jets_constituents_pdgId->at(jetIt)[jj]) == 22  ) emE_raw += jetConst_raw_p4.E();
              else if( fabs(jets_constituents_pdgId->at(jetIt)[jj]) == 130 ) neutralE_raw += jetConst_raw_p4.E();
              else if( fabs(jets_constituents_pdgId->at(jetIt)[jj]) <= 16 )  leptonE_raw += jetConst_raw_p4.E();
              if(      fabs(jets_constituents_pdgId->at(jetIt)[jj]) == 211 ) chargedE_weighted += jetConst_weighted_p4.E();
              else if( fabs(jets_constituents_pdgId->at(jetIt)[jj]) == 22  ) emE_weighted += jetConst_weighted_p4.E();
              else if( fabs(jets_constituents_pdgId->at(jetIt)[jj]) == 130 ) neutralE_weighted += jetConst_weighted_p4.E();
              else if( fabs(jets_constituents_pdgId->at(jetIt)[jj]) <= 16 )  leptonE_weighted += jetConst_weighted_p4.E();
              
              if(      fabs(jets_constituents_pdgId->at(jetIt)[jj]) == 211 ) chargedPt_raw += jetConst_raw_p4.Pt();
              else if( fabs(jets_constituents_pdgId->at(jetIt)[jj]) == 22  ) emPt_raw += jetConst_raw_p4.Pt();
              else if( fabs(jets_constituents_pdgId->at(jetIt)[jj]) == 130 ) neutralPt_raw += jetConst_raw_p4.Pt();
              else if( fabs(jets_constituents_pdgId->at(jetIt)[jj]) <= 16 )  leptonPt_raw += jetConst_raw_p4.Pt();
              if(      fabs(jets_constituents_pdgId->at(jetIt)[jj]) == 211 ) chargedPt_weighted += jetConst_weighted_p4.Pt();
              else if( fabs(jets_constituents_pdgId->at(jetIt)[jj]) == 22  ) emPt_weighted += jetConst_weighted_p4.Pt();
              else if( fabs(jets_constituents_pdgId->at(jetIt)[jj]) == 130 ) neutralPt_weighted += jetConst_weighted_p4.Pt();
              else if( fabs(jets_constituents_pdgId->at(jetIt)[jj]) <= 16 )  leptonPt_weighted += jetConst_weighted_p4.Pt();
            }
            
            /*
            if( chargedE_weighted/genJet_p4.E() < 0.05 )
            {
              std::cout << "\n\n\n" << std::endl;
              std::cout << "genLepton1:   pt = " << genLeptons_pt->at(0) << "   eta = " << genLeptons_eta->at(0) << "   phi = " << genLeptons_phi->at(0) << "   pdgId = " << genLeptons_pdgId->at(0) << std::endl;
              std::cout << "genLepton2:   pt = " << genLeptons_pt->at(1) << "   eta = " << genLeptons_eta->at(1) << "   phi = " << genLeptons_phi->at(1) << "   pdgId = " << genLeptons_pdgId->at(1) << std::endl;
              std::cout << "jet:   pt = " << jet_p4.Pt() << "   eta = " << jet_p4.Eta() << "   phi = " << jet_p4.Phi() << std::endl;
              for(unsigned int jj = 0; jj < jets_constituents_pt->at(jetIt).size(); ++jj)                                                                                                                                                  
              {                                                                                                                                                                                                                            
                TLorentzVector jetConst_raw_p4;                                                                                                                                                                                            
                jetConst_raw_p4.SetPtEtaPhiM(jets_constituents_pt->at(jetIt)[jj],                                                                                                                                                          
                                             jets_constituents_eta->at(jetIt)[jj],                                                                                                                                                         
                                             jets_constituents_phi->at(jetIt)[jj],                                                                                                                                                         
                                             jets_constituents_mass->at(jetIt)[jj]);                                                                                                                                                       
                std::cout << ">>> const " << jj << ":   pt =  " << jets_constituents_pt->at(jetIt)[jj] << "   eta = " << jets_constituents_eta->at(jetIt)[jj] << "   phi = " << jets_constituents_phi->at(jetIt)[jj] << "   pdgId: " << jets_constituents_pdgId->at(jetIt)[jj] << "   puppiWeight: " << jets_constituents_puppiWeight->at(jetIt)[jj] << std::endl;
              }
              std::cout << ">>> charged E = " << chargedE_weighted << std::endl;
              std::cout << ">>>      em E = " << emE_weighted << std::endl;
              std::cout << ">>> neutral E = " << neutralE_weighted << std::endl;
              std::cout << ">>>  lepton E = " << leptonE_weighted << std::endl;
            }
            */
            
            h_ERes_charged_all[label] -> Fill( chargedE_weighted/chargedE_gen,weight );
            h_ERes_em_all[label] -> Fill( emE_weighted/emE_gen,weight );
            h_ERes_neutral_all[label] -> Fill( neutralE_weighted/neutralE_gen,weight );
            h_ERes_lepton_all[label] -> Fill( leptonE_weighted/leptonE_gen,weight );
            
            if( chargedE_weighted > 0.01 )
              h_ERes_charged_pt[label][std::make_pair(ptRangeMin,ptRangeMax)] -> Fill( chargedE_weighted/chargedE_gen,weight );
            if( emE_weighted > 0.01 )
              h_ERes_em_pt[label][std::make_pair(ptRangeMin,ptRangeMax)] -> Fill( emE_weighted/emE_gen,weight );
            if( neutralE_weighted > 0.01 )
              h_ERes_neutral_pt[label][std::make_pair(ptRangeMin,ptRangeMax)] -> Fill( neutralE_weighted/neutralE_gen,weight );
            if( leptonE_weighted > 0.01 )
              h_ERes_lepton_pt[label][std::make_pair(ptRangeMin,ptRangeMax)] -> Fill( leptonE_weighted/leptonE_gen,weight );
            
            p_chargedPt_gen_prompt[label] -> Fill( genJet_p4.Pt(),chargedPt_gen );
            p_emPt_gen_prompt[label] -> Fill( genJet_p4.Pt(),emPt_gen );
            p_neutralPt_gen_prompt[label] -> Fill( genJet_p4.Pt(),neutralPt_gen );
            p_leptonPt_gen_prompt[label] -> Fill( genJet_p4.Pt(),leptonPt_gen );
            p_chargedPt_raw_prompt[label] -> Fill( genJet_p4.Pt(),chargedPt_raw );
            p_emPt_raw_prompt[label] -> Fill( genJet_p4.Pt(),emPt_raw );
            p_neutralPt_raw_prompt[label] -> Fill( genJet_p4.Pt(),neutralPt_raw );
            p_leptonPt_raw_prompt[label] -> Fill( genJet_p4.Pt(),leptonPt_raw );
            p_chargedPt_weighted_prompt[label] -> Fill( genJet_p4.Pt(),chargedPt_weighted );
            p_emPt_weighted_prompt[label] -> Fill( genJet_p4.Pt(),emPt_weighted );
            p_neutralPt_weighted_prompt[label] -> Fill( genJet_p4.Pt(),neutralPt_weighted );
            p_leptonPt_weighted_prompt[label] -> Fill( genJet_p4.Pt(),leptonPt_weighted );
            
            p_chargedPtFrac_gen_prompt[label] -> Fill( genJet_p4.Pt(),chargedPt_gen/genJet_p4.Pt() );
            p_emPtFrac_gen_prompt[label] -> Fill( genJet_p4.Pt(),emPt_gen/genJet_p4.Pt() );
            p_neutralPtFrac_gen_prompt[label] -> Fill( genJet_p4.Pt(),neutralPt_gen/genJet_p4.Pt() );
            p_leptonPtFrac_gen_prompt[label] -> Fill( genJet_p4.Pt(),leptonPt_gen/genJet_p4.Pt() );
            p_chargedPtFrac_raw_prompt[label] -> Fill( genJet_p4.Pt(),chargedPt_raw/genJet_p4.Pt() );
            p_emPtFrac_raw_prompt[label] -> Fill( genJet_p4.Pt(),emPt_raw/genJet_p4.Pt() );
            p_neutralPtFrac_raw_prompt[label] -> Fill( genJet_p4.Pt(),neutralPt_raw/genJet_p4.Pt() );
            p_leptonPtFrac_raw_prompt[label] -> Fill( genJet_p4.Pt(),leptonPt_raw/genJet_p4.Pt() );
            p_chargedPtFrac_weighted_prompt[label] -> Fill( genJet_p4.Pt(),chargedPt_weighted/genJet_p4.Pt() );
            p_emPtFrac_weighted_prompt[label] -> Fill( genJet_p4.Pt(),emPt_weighted/genJet_p4.Pt() );
            p_neutralPtFrac_weighted_prompt[label] -> Fill( genJet_p4.Pt(),neutralPt_weighted/genJet_p4.Pt() );
            p_leptonPtFrac_weighted_prompt[label] -> Fill( genJet_p4.Pt(),leptonPt_weighted/genJet_p4.Pt() );
            
            // if( leptonPt_weighted/neutralPt_weighted > 0.3 )
            // {
            //   std::cout << std::endl;
            //   std::cout << std::endl;
            //   std::cout << std::endl;
            //   std::cout << "genLepton1:   pt = " << genLeptons_pt->at(0) << "   eta = " << genLeptons_eta->at(0) << "   phi = " << genLeptons_phi->at(0) << "   pdgId = " << genLeptons_pdgId->at(0) << std::endl;
            //   std::cout << "genLepton2:   pt = " << genLeptons_pt->at(1) << "   eta = " << genLeptons_eta->at(1) << "   phi = " << genLeptons_phi->at(1) << "   pdgId = " << genLeptons_pdgId->at(1) << std::endl;
            //   std::cout << std::endl;
            //   std::cout << "jet:   pt = " << jet_p4.Pt() << "   eta = " << jet_p4.Eta() << "   phi: " << jet_p4.Phi() << std::endl;
            //   std::cout << std::endl;
            //   for(unsigned int jj = 0; jj < jets_constituents_pt->at(jetIt).size(); ++jj)
            //   {
            //     TLorentzVector jetConst_raw_p4;
            //     jetConst_raw_p4.SetPtEtaPhiM(jets_constituents_pt->at(jetIt)[jj],
            //                                  jets_constituents_eta->at(jetIt)[jj],
            //                                  jets_constituents_phi->at(jetIt)[jj],
            //                                  jets_constituents_mass->at(jetIt)[jj]);
            //     std::cout << ">>> const " << jj << ":   pt =  " << jets_constituents_pt->at(jetIt)[jj] << "   eta = " << jets_constituents_eta->at(jetIt)[jj] << "   phi = " << jets_constituents_phi->at(jetIt)[jj] << "   pdgId: " << jets_constituents_pdgId->at(jetIt)[jj] << std::endl;
            //   }              
            // }
          }
        }
        else if( jets_isPU->at(jetIt) == 1 )
        {
          h_pt[label][1] -> Fill( jets_pt->at(jetIt),weight );
          h_eta[label][1] -> Fill( jets_eta->at(jetIt),weight );
          ++nJets_fake;
        }
        else
        {
          h_pt[label][2] -> Fill( jets_pt->at(jetIt),weight );
          h_eta[label][2] -> Fill( jets_eta->at(jetIt),weight );
          ++nJets_und;
        }
        
        ++itJet;
      }
      if( nJets_prompt > 0 ) h_nJets[label][0] -> Fill( nJets_prompt,weight );
      if( nJets_fake > 0 )   h_nJets[label][1] -> Fill( nJets_fake,weight );
      if( nJets_und > 0 )    h_nJets[label][2] -> Fill( nJets_und,weight );
      if( nJets_all > 0 )    h_nJets[label][3] -> Fill( nJets_und,weight );
    }
    std::cout << std::endl;
    
  } // NoPU / PU200
  
  
  
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir");
  DrawPlot(h_nJets[uniqueLabels.at(0)],h_nJets[uniqueLabels.at(1)],uniqueLabels.at(0),uniqueLabels.at(1),"nJets",";N_{jets};events / 1 pb^{-1}",        "",plotDir);
  DrawPlot(h_pt[uniqueLabels.at(0)],   h_pt[uniqueLabels.at(1)],   uniqueLabels.at(0),uniqueLabels.at(1),"pt",   ";jet p_{T} [GeV];events / 1 pb^{-1}", "",plotDir);
  DrawPlot(h_eta[uniqueLabels.at(0)],  h_eta[uniqueLabels.at(1)],  uniqueLabels.at(0),uniqueLabels.at(1),"eta",   ";jet #eta;events / 1 pb^{-1}",       "",plotDir);  
  
  DrawPlot(h_ERes_all[uniqueLabels.at(0)],h_ERes_all[uniqueLabels.at(1)],uniqueLabels.at(0),uniqueLabels.at(1),"ERes",   ";jet energy / genJet energy;events / 1 pb^{-1}", "",plotDir);
  DrawPlot(h_ERes_charged_all[uniqueLabels.at(0)],h_ERes_charged_all[uniqueLabels.at(1)],uniqueLabels.at(0),uniqueLabels.at(1),"ERes_charged",   ";jet charged energy / genJet energy;events / 1 pb^{-1}", "",plotDir);
  DrawPlot(h_ERes_em_all[uniqueLabels.at(0)],h_ERes_em_all[uniqueLabels.at(1)],uniqueLabels.at(0),uniqueLabels.at(1),"ERes_em",   ";jet energy / genJet em energy;events / 1 pb^{-1}", "",plotDir);
  DrawPlot(h_ERes_neutral_all[uniqueLabels.at(0)],h_ERes_neutral_all[uniqueLabels.at(1)],uniqueLabels.at(0),uniqueLabels.at(1),"ERes_neutral",   ";jet neutral energy / genJet energy;events / 1 pb^{-1}", "",plotDir,true);
  DrawPlot(h_ERes_lepton_all[uniqueLabels.at(0)],h_ERes_lepton_all[uniqueLabels.at(1)],uniqueLabels.at(0),uniqueLabels.at(1),"ERes_lepton",   ";jet lepton energy / genJet energy;events / 1 pb^{-1}", "",plotDir,true);
  
  
  float* vals = new float[6];
  TH1F* histo;
  
  
  //--- resolution plots
  std::map<std::string,TGraphErrors*> g_ERes_vs_pt;
  std::map<std::string,TGraphErrors*> g_ERes_charged_vs_pt;
  std::map<std::string,TGraphErrors*> g_ERes_em_vs_pt;
  std::map<std::string,TGraphErrors*> g_ERes_neutral_vs_pt;
  std::map<std::string,TGraphErrors*> g_ERes_lepton_vs_pt;
  
  for(unsigned int uniqueLabelIt = 0; uniqueLabelIt < uniqueLabels.size(); ++uniqueLabelIt)
  {
    std::string uniqueLabel = uniqueLabels.at(uniqueLabelIt);
    g_ERes_vs_pt[uniqueLabel] = new TGraphErrors();
    g_ERes_charged_vs_pt[uniqueLabel] = new TGraphErrors();
    g_ERes_em_vs_pt[uniqueLabel] = new TGraphErrors();
    g_ERes_neutral_vs_pt[uniqueLabel] = new TGraphErrors();
    g_ERes_lepton_vs_pt[uniqueLabel] = new TGraphErrors();
    
    for(unsigned int ptIt = 0; ptIt < ptRanges.size()-1; ++ptIt)
    {
      float ptRangeMin = ptRanges.at(ptIt);
      float ptRangeMax = ptRanges.at(ptIt+1);
      
      histo = h_ERes_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)];
      FindSmallestInterval(vals,histo,0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      g_ERes_vs_pt[uniqueLabel] -> SetPoint(ptIt,0.5*(ptRangeMin+ptRangeMax),sigma);
      g_ERes_vs_pt[uniqueLabel] -> SetPointError(ptIt,0.5*(ptRangeMax-ptRangeMin),0.01);
      
      histo = h_ERes_charged_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)];
      FindSmallestInterval(vals,histo,0.68);
      mean = vals[0];
      min = vals[4];
      max = vals[5];
      delta = max-min;
      sigma = 0.5*delta;
      g_ERes_charged_vs_pt[uniqueLabel] -> SetPoint(ptIt,0.5*(ptRangeMin+ptRangeMax),sigma/mean);
      g_ERes_charged_vs_pt[uniqueLabel] -> SetPointError(ptIt,0.5*(ptRangeMax-ptRangeMin),0.01);
      
      histo = h_ERes_em_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)];
      FindSmallestInterval(vals,histo,0.68);
      mean = vals[0];
      min = vals[4];
      max = vals[5];
      delta = max-min;
      sigma = 0.5*delta;
      g_ERes_em_vs_pt[uniqueLabel] -> SetPoint(ptIt,0.5*(ptRangeMin+ptRangeMax),sigma/mean);
      g_ERes_em_vs_pt[uniqueLabel] -> SetPointError(ptIt,0.5*(ptRangeMax-ptRangeMin),0.01);
      
      histo = h_ERes_neutral_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)];
      FindSmallestInterval(vals,histo,0.68);
      mean = vals[0];
      min = vals[4];
      max = vals[5];
      delta = max-min;
      sigma = 0.5*delta;
      g_ERes_neutral_vs_pt[uniqueLabel] -> SetPoint(ptIt,0.5*(ptRangeMin+ptRangeMax),sigma/mean);
      g_ERes_neutral_vs_pt[uniqueLabel] -> SetPointError(ptIt,0.5*(ptRangeMax-ptRangeMin),0.01);
      
      histo = h_ERes_lepton_pt[uniqueLabel][std::make_pair(ptRangeMin,ptRangeMax)];
      FindSmallestInterval(vals,histo,0.68);
      mean = vals[0];
      min = vals[4];
      max = vals[5];
      delta = max-min;
      sigma = 0.5*delta;
      g_ERes_lepton_vs_pt[uniqueLabel] -> SetPoint(ptIt,0.5*(ptRangeMin+ptRangeMax),sigma/mean);
      g_ERes_lepton_vs_pt[uniqueLabel] -> SetPointError(ptIt,0.5*(ptRangeMax-ptRangeMin),0.01);
    }
  }
  
  TCanvas* c = new TCanvas();
  c -> cd();
  
  TH1F* hPad = (TH1F*)( gPad->DrawFrame(0.,0.,200.,1.) );
  hPad -> SetTitle("; p_{T} [GeV];#sigma_{E} / E");
  hPad -> Draw();
  gPad -> SetGridy();
  
  std::vector<int> colors;
  colors.push_back(kRed);
  colors.push_back(16);
  
  TLegend* legend = new TLegend(0.52,0.90-0.08*4,0.99,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  
  for(unsigned int uniqueLabelIt = 0; uniqueLabelIt < uniqueLabels.size(); ++uniqueLabelIt)
  {
    std::string uniqueLabel = uniqueLabels.at(uniqueLabelIt);
    
    g_ERes_vs_pt[uniqueLabel] -> SetLineColor(colors.at(uniqueLabelIt));
    g_ERes_vs_pt[uniqueLabel] -> SetMarkerColor(colors.at(uniqueLabelIt));
    g_ERes_vs_pt[uniqueLabel] -> SetMarkerStyle(20);
    g_ERes_vs_pt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(g_ERes_vs_pt[uniqueLabel],Form("%s - %.1f < |#eta| < %.1f",uniqueLabel.c_str(),etaMin,etaMax),"PL");
  }
  
  legend -> Draw("same");
  
  c -> Print(Form("%s/c_ERes_vs_pt.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_ERes_vs_pt.png",plotDir.c_str()));
  delete c;
  
  
  
  c = new TCanvas();
  c -> cd();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,200.,1.) );
  hPad -> SetTitle("; p_{T} [GeV];#sigma_{E} / E");
  hPad -> Draw();
  gPad -> SetGridy();
  
  for(unsigned int uniqueLabelIt = 0; uniqueLabelIt < uniqueLabels.size(); ++uniqueLabelIt)
  {
    std::string uniqueLabel = uniqueLabels.at(uniqueLabelIt);
    
    g_ERes_charged_vs_pt[uniqueLabel] -> SetLineColor(colors.at(uniqueLabelIt));
    g_ERes_charged_vs_pt[uniqueLabel] -> SetMarkerColor(colors.at(uniqueLabelIt));
    g_ERes_charged_vs_pt[uniqueLabel] -> SetMarkerStyle(20);
    g_ERes_charged_vs_pt[uniqueLabel] -> Draw("PL,same");
  }
  
  legend -> Draw("same");
  
  c -> Print(Form("%s/c_ERes_charged_vs_pt.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_ERes_charged_vs_pt.png",plotDir.c_str()));
  delete c;
  
  
  
  c = new TCanvas();
  c -> cd();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,200.,1.) );
  hPad -> SetTitle("; p_{T} [GeV];#sigma_{E} / E");
  hPad -> Draw();
  gPad -> SetGridy();
  
  for(unsigned int uniqueLabelIt = 0; uniqueLabelIt < uniqueLabels.size(); ++uniqueLabelIt)
  {
    std::string uniqueLabel = uniqueLabels.at(uniqueLabelIt);
    
    g_ERes_em_vs_pt[uniqueLabel] -> SetLineColor(colors.at(uniqueLabelIt));
    g_ERes_em_vs_pt[uniqueLabel] -> SetMarkerColor(colors.at(uniqueLabelIt));
    g_ERes_em_vs_pt[uniqueLabel] -> SetMarkerStyle(20);
    g_ERes_em_vs_pt[uniqueLabel] -> Draw("PL,same");
  }
  
  legend -> Draw("same");
  
  c -> Print(Form("%s/c_ERes_em_vs_pt.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_ERes_em_vs_pt.png",plotDir.c_str()));
  delete c;
  
  
  
  c = new TCanvas();
  c -> cd();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,200.,1.) );
  hPad -> SetTitle("; p_{T} [GeV];#sigma_{E} / E");
  hPad -> Draw();
  gPad -> SetGridy();
  
  for(unsigned int uniqueLabelIt = 0; uniqueLabelIt < uniqueLabels.size(); ++uniqueLabelIt)
  {
    std::string uniqueLabel = uniqueLabels.at(uniqueLabelIt);
    
    g_ERes_neutral_vs_pt[uniqueLabel] -> SetLineColor(colors.at(uniqueLabelIt));
    g_ERes_neutral_vs_pt[uniqueLabel] -> SetMarkerColor(colors.at(uniqueLabelIt));
    g_ERes_neutral_vs_pt[uniqueLabel] -> SetMarkerStyle(20);
    g_ERes_neutral_vs_pt[uniqueLabel] -> Draw("PL,same");
  }
  
  legend -> Draw("same");
  
  c -> Print(Form("%s/c_ERes_neutral_vs_pt.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_ERes_neutral_vs_pt.png",plotDir.c_str()));
  delete c;
  
  
  
  c = new TCanvas();
  c -> cd();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,200.,1.) );
  hPad -> SetTitle("; p_{T} [GeV];#sigma_{E} / E");
  hPad -> Draw();
  gPad -> SetGridy();
  
  for(unsigned int uniqueLabelIt = 0; uniqueLabelIt < uniqueLabels.size(); ++uniqueLabelIt)
  {
    std::string uniqueLabel = uniqueLabels.at(uniqueLabelIt);
    
    g_ERes_lepton_vs_pt[uniqueLabel] -> SetLineColor(colors.at(uniqueLabelIt));
    g_ERes_lepton_vs_pt[uniqueLabel] -> SetMarkerColor(colors.at(uniqueLabelIt));
    g_ERes_lepton_vs_pt[uniqueLabel] -> SetMarkerStyle(20);
    g_ERes_lepton_vs_pt[uniqueLabel] -> Draw("PL,same");
  }
  
  legend -> Draw("same");
  
  c -> Print(Form("%s/c_ERes_lepton_vs_pt.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_ERes_lepton_vs_pt.png",plotDir.c_str()));
  delete c;
  
  
  
  
  
  
  //--- ptSum plots
  c = new TCanvas();
  c -> cd();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,200.,1.2) );
  hPad -> SetTitle("; p_{T} [GeV];#Sigma_{comp.}p_{T} / p_{T}_{genJet}");
  hPad -> Draw();
  gPad -> SetGridy();
  
  legend = new TLegend(0.52,0.95-0.04*8,0.99,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  
  for(unsigned int uniqueLabelIt = 0; uniqueLabelIt < uniqueLabels.size(); ++uniqueLabelIt)
  {
    std::string uniqueLabel = uniqueLabels.at(uniqueLabelIt);
    
    p_chargedPtFrac_raw_prompt[uniqueLabel] -> SetLineColor(kRed);
    p_chargedPtFrac_raw_prompt[uniqueLabel] -> SetMarkerColor(kRed);
    p_chargedPtFrac_raw_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_chargedPtFrac_raw_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_chargedPtFrac_raw_prompt[uniqueLabel],Form("%s - charged component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_chargedPtFrac_raw_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_emPtFrac_raw_prompt[uniqueLabel] -> SetLineColor(kOrange);
    p_emPtFrac_raw_prompt[uniqueLabel] -> SetMarkerColor(kOrange);
    p_emPtFrac_raw_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_emPtFrac_raw_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_emPtFrac_raw_prompt[uniqueLabel],Form("%s - em component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_emPtFrac_raw_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_neutralPtFrac_raw_prompt[uniqueLabel] -> SetLineColor(kBlue);
    p_neutralPtFrac_raw_prompt[uniqueLabel] -> SetMarkerColor(kBlue);
    p_neutralPtFrac_raw_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_neutralPtFrac_raw_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_neutralPtFrac_raw_prompt[uniqueLabel],Form("%s - neutral component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_neutralPtFrac_raw_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_leptonPtFrac_raw_prompt[uniqueLabel] -> SetLineColor(kGreen+2);
    p_leptonPtFrac_raw_prompt[uniqueLabel] -> SetMarkerColor(kGreen+2);
    p_leptonPtFrac_raw_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_leptonPtFrac_raw_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_leptonPtFrac_raw_prompt[uniqueLabel],Form("%s - lepton component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_leptonPtFrac_raw_prompt[uniqueLabel] -> SetMarkerStyle(24);
  }
  
  legend -> Draw("same");
  
  c -> Print(Form("%s/c_ptFrac_raw_vs_pt.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_ptFrac_raw_vs_pt.png",plotDir.c_str()));
  delete c;
  
  
  
  c = new TCanvas();
  c -> cd();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,200.,150.) );
  hPad -> SetTitle("; p_{T} [GeV];#Sigma_{comp.}p_{T} [GeV]");
  hPad -> Draw();
  gPad -> SetGridy();
  
  legend = new TLegend(0.52,0.95-0.04*8,0.99,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  
  for(unsigned int uniqueLabelIt = 0; uniqueLabelIt < uniqueLabels.size(); ++uniqueLabelIt)
  {
    std::string uniqueLabel = uniqueLabels.at(uniqueLabelIt);
    
    p_chargedPt_raw_prompt[uniqueLabel] -> SetLineColor(kRed);
    p_chargedPt_raw_prompt[uniqueLabel] -> SetMarkerColor(kRed);
    p_chargedPt_raw_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_chargedPt_raw_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_chargedPt_raw_prompt[uniqueLabel],Form("%s - charged component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_chargedPt_raw_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_emPt_raw_prompt[uniqueLabel] -> SetLineColor(kOrange);
    p_emPt_raw_prompt[uniqueLabel] -> SetMarkerColor(kOrange);
    p_emPt_raw_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_emPt_raw_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_emPt_raw_prompt[uniqueLabel],Form("%s - em component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_emPt_raw_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_neutralPt_raw_prompt[uniqueLabel] -> SetLineColor(kBlue);
    p_neutralPt_raw_prompt[uniqueLabel] -> SetMarkerColor(kBlue);
    p_neutralPt_raw_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_neutralPt_raw_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_neutralPt_raw_prompt[uniqueLabel],Form("%s - neutral component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_neutralPt_raw_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_leptonPt_raw_prompt[uniqueLabel] -> SetLineColor(kGreen+2);
    p_leptonPt_raw_prompt[uniqueLabel] -> SetMarkerColor(kGreen+2);
    p_leptonPt_raw_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_leptonPt_raw_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_leptonPt_raw_prompt[uniqueLabel],Form("%s - lepton component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_leptonPt_raw_prompt[uniqueLabel] -> SetMarkerStyle(24);
  }
  
  legend -> Draw("same");
  
  c -> Print(Form("%s/c_ptSum_raw_vs_pt.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_ptSum_raw_vs_pt.png",plotDir.c_str()));
  delete c;
  
  
  
  c = new TCanvas();
  c -> cd();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,200.,1.2) );
  hPad -> SetTitle("; p_{T} [GeV];#Sigma_{comp.}p_{T} / p_{T}_{genJet}");
  hPad -> Draw();
  gPad -> SetGridy();
  
  legend = new TLegend(0.52,0.95-0.04*8,0.99,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  
  for(unsigned int uniqueLabelIt = 0; uniqueLabelIt < uniqueLabels.size(); ++uniqueLabelIt)
  {
    std::string uniqueLabel = uniqueLabels.at(uniqueLabelIt);
    
    p_chargedPtFrac_weighted_prompt[uniqueLabel] -> SetLineColor(kRed);
    p_chargedPtFrac_weighted_prompt[uniqueLabel] -> SetMarkerColor(kRed);
    p_chargedPtFrac_weighted_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_chargedPtFrac_weighted_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_chargedPtFrac_weighted_prompt[uniqueLabel],Form("%s - charged component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_chargedPtFrac_weighted_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_emPtFrac_weighted_prompt[uniqueLabel] -> SetLineColor(kOrange);
    p_emPtFrac_weighted_prompt[uniqueLabel] -> SetMarkerColor(kOrange);
    p_emPtFrac_weighted_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_emPtFrac_weighted_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_emPtFrac_weighted_prompt[uniqueLabel],Form("%s - em component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_emPtFrac_weighted_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_neutralPtFrac_weighted_prompt[uniqueLabel] -> SetLineColor(kBlue);
    p_neutralPtFrac_weighted_prompt[uniqueLabel] -> SetMarkerColor(kBlue);
    p_neutralPtFrac_weighted_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_neutralPtFrac_weighted_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_neutralPtFrac_weighted_prompt[uniqueLabel],Form("%s - neutral component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_neutralPtFrac_weighted_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_leptonPtFrac_weighted_prompt[uniqueLabel] -> SetLineColor(kGreen+2);
    p_leptonPtFrac_weighted_prompt[uniqueLabel] -> SetMarkerColor(kGreen+2);
    p_leptonPtFrac_weighted_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_leptonPtFrac_weighted_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_leptonPtFrac_weighted_prompt[uniqueLabel],Form("%s - lepton component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_leptonPtFrac_weighted_prompt[uniqueLabel] -> SetMarkerStyle(24);
  }
  
  legend -> Draw("same");
  
  c -> Print(Form("%s/c_ptFrac_weighted_vs_pt.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_ptFrac_weighted_vs_pt.png",plotDir.c_str()));
  delete c;
  
  
  
  c = new TCanvas();
  c -> cd();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,200.,150.) );
  hPad -> SetTitle("; p_{T} [GeV];#Sigma_{comp.}p_{T} [GeV]");
  hPad -> Draw();
  gPad -> SetGridy();
  
  legend = new TLegend(0.52,0.95-0.04*8,0.99,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  
  for(unsigned int uniqueLabelIt = 0; uniqueLabelIt < uniqueLabels.size(); ++uniqueLabelIt)
  {
    std::string uniqueLabel = uniqueLabels.at(uniqueLabelIt);
    
    p_chargedPt_weighted_prompt[uniqueLabel] -> SetLineColor(kRed);
    p_chargedPt_weighted_prompt[uniqueLabel] -> SetMarkerColor(kRed);
    p_chargedPt_weighted_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_chargedPt_weighted_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_chargedPt_weighted_prompt[uniqueLabel],Form("%s - charged component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_chargedPt_weighted_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_emPt_weighted_prompt[uniqueLabel] -> SetLineColor(kOrange);
    p_emPt_weighted_prompt[uniqueLabel] -> SetMarkerColor(kOrange);
    p_emPt_weighted_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_emPt_weighted_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_emPt_weighted_prompt[uniqueLabel],Form("%s - em component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_emPt_weighted_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_neutralPt_weighted_prompt[uniqueLabel] -> SetLineColor(kBlue);
    p_neutralPt_weighted_prompt[uniqueLabel] -> SetMarkerColor(kBlue);
    p_neutralPt_weighted_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_neutralPt_weighted_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_neutralPt_weighted_prompt[uniqueLabel],Form("%s - neutral component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_neutralPt_weighted_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_leptonPt_weighted_prompt[uniqueLabel] -> SetLineColor(kGreen+2);
    p_leptonPt_weighted_prompt[uniqueLabel] -> SetMarkerColor(kGreen+2);
    p_leptonPt_weighted_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_leptonPt_weighted_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_leptonPt_weighted_prompt[uniqueLabel],Form("%s - lepton component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_leptonPt_weighted_prompt[uniqueLabel] -> SetMarkerStyle(24);
  }
  
  legend -> Draw("same");
  
  c -> Print(Form("%s/c_ptSum_weighted_vs_pt.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_ptSum_weighted_vs_pt.png",plotDir.c_str()));
  delete c;
  
  
  
  c = new TCanvas();
  c -> cd();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,200.,1.2) );
  hPad -> SetTitle("; p_{T} [GeV];#Sigma_{comp.}p_{T} / p_{T}_{genJet}");
  hPad -> Draw();
  gPad -> SetGridy();
  
  legend = new TLegend(0.52,0.95-0.04*8,0.99,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  
  for(unsigned int uniqueLabelIt = 0; uniqueLabelIt < uniqueLabels.size(); ++uniqueLabelIt)
  {
    std::string uniqueLabel = uniqueLabels.at(uniqueLabelIt);
    
    p_chargedPtFrac_gen_prompt[uniqueLabel] -> SetLineColor(kRed);
    p_chargedPtFrac_gen_prompt[uniqueLabel] -> SetMarkerColor(kRed);
    p_chargedPtFrac_gen_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_chargedPtFrac_gen_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_chargedPtFrac_gen_prompt[uniqueLabel],Form("%s - charged component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_chargedPtFrac_gen_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_emPtFrac_gen_prompt[uniqueLabel] -> SetLineColor(kOrange);
    p_emPtFrac_gen_prompt[uniqueLabel] -> SetMarkerColor(kOrange);
    p_emPtFrac_gen_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_emPtFrac_gen_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_emPtFrac_gen_prompt[uniqueLabel],Form("%s - em component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_emPtFrac_gen_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_neutralPtFrac_gen_prompt[uniqueLabel] -> SetLineColor(kBlue);
    p_neutralPtFrac_gen_prompt[uniqueLabel] -> SetMarkerColor(kBlue);
    p_neutralPtFrac_gen_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_neutralPtFrac_gen_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_neutralPtFrac_gen_prompt[uniqueLabel],Form("%s - neutral component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_neutralPtFrac_gen_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_leptonPtFrac_gen_prompt[uniqueLabel] -> SetLineColor(kGreen+2);
    p_leptonPtFrac_gen_prompt[uniqueLabel] -> SetMarkerColor(kGreen+2);
    p_leptonPtFrac_gen_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_leptonPtFrac_gen_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_leptonPtFrac_gen_prompt[uniqueLabel],Form("%s - lepton component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_leptonPtFrac_gen_prompt[uniqueLabel] -> SetMarkerStyle(24);
  }
  
  legend -> Draw("same");
  
  c -> Print(Form("%s/c_ptFrac_gen_vs_pt.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_ptFrac_gen_vs_pt.png",plotDir.c_str()));
  delete c;
  
  
  
  c = new TCanvas();
  c -> cd();
  
  hPad = (TH1F*)( gPad->DrawFrame(0.,0.,200.,150.) );
  hPad -> SetTitle("; p_{T} [GeV];#Sigma_{comp.}p_{T} [GeV]");
  hPad -> Draw();
  gPad -> SetGridy();
  
  legend = new TLegend(0.52,0.95-0.04*8,0.99,0.90);
  legend -> SetFillColor(0);
  legend -> SetFillStyle(1000);
  legend -> SetTextFont(42);
  legend -> SetTextSize(0.03);
  
  for(unsigned int uniqueLabelIt = 0; uniqueLabelIt < uniqueLabels.size(); ++uniqueLabelIt)
  {
    std::string uniqueLabel = uniqueLabels.at(uniqueLabelIt);
    
    p_chargedPt_gen_prompt[uniqueLabel] -> SetLineColor(kRed);
    p_chargedPt_gen_prompt[uniqueLabel] -> SetMarkerColor(kRed);
    p_chargedPt_gen_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_chargedPt_gen_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_chargedPt_gen_prompt[uniqueLabel],Form("%s - charged component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_chargedPt_gen_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_emPt_gen_prompt[uniqueLabel] -> SetLineColor(kOrange);
    p_emPt_gen_prompt[uniqueLabel] -> SetMarkerColor(kOrange);
    p_emPt_gen_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_emPt_gen_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_emPt_gen_prompt[uniqueLabel],Form("%s - em component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_emPt_gen_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_neutralPt_gen_prompt[uniqueLabel] -> SetLineColor(kBlue);
    p_neutralPt_gen_prompt[uniqueLabel] -> SetMarkerColor(kBlue);
    p_neutralPt_gen_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_neutralPt_gen_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_neutralPt_gen_prompt[uniqueLabel],Form("%s - neutral component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_neutralPt_gen_prompt[uniqueLabel] -> SetMarkerStyle(24);
    
    p_leptonPt_gen_prompt[uniqueLabel] -> SetLineColor(kGreen+2);
    p_leptonPt_gen_prompt[uniqueLabel] -> SetMarkerColor(kGreen+2);
    p_leptonPt_gen_prompt[uniqueLabel] -> SetMarkerStyle(20);
    p_leptonPt_gen_prompt[uniqueLabel] -> Draw("PL,same");
    legend -> AddEntry(p_leptonPt_gen_prompt[uniqueLabel],Form("%s - lepton component",uniqueLabel.c_str()),"P");
    if( uniqueLabelIt == 1 ) p_leptonPt_gen_prompt[uniqueLabel] -> SetMarkerStyle(24);
  }
  
  legend -> Draw("same");
  
  c -> Print(Form("%s/c_ptSum_gen_vs_pt.pdf",plotDir.c_str()));
  c -> Print(Form("%s/c_ptSum_gen_vs_pt.png",plotDir.c_str()));
  delete c;  
  
  
  
  outFile -> Write();
}
