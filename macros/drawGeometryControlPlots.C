void drawGeometryControlPlots(const std::string& inFileName)
{
  TFile* inFile = TFile::Open(inFileName.c_str(),"READ");
  TTree* t = (TTree*)( inFile->Get("FTLDumpHits/hits_tree"));
  
  

  // std::string label = "tile";
  // TLatex* latexLabel = new TLatex(0.16,0.96,Form("11.5#times11.5#timesd mm^{3} tiles"));
  // std::string label = "barphi";
  // TLatex* latexLabel = new TLatex(0.16,0.96,Form("3#times50#timesd mm^{3} bars along #phi"));
  // std::string label = "barz";
  // TLatex* latexLabel = new TLatex(0.16,0.96,Form("3#times50#timesd mm^{3} bars along z"));
  std::string label = "barzflat";
  TLatex* latexLabel = new TLatex(0.16,0.96,Form("3#times56#timesd mm^{3} flat bars along z"));
  latexLabel -> SetNDC();
  latexLabel -> SetTextFont(42);
  latexLabel -> SetTextSize(0.03);
  
  
  
  TCanvas* c1 = new TCanvas("c1","c1",1400,500);
  c1 -> Divide(3,1);

  c1 -> cd(1);
  gPad -> SetLogy();
  TH1F* h_entry_x = new TH1F("h_entry_x","",30000,-3.,3.);
  t -> Draw("simHits_entry_local_x>>h_entry_x","","goff");
  h_entry_x -> SetTitle(";simHit entry/exit point x (cm);entries");
  h_entry_x -> Draw();
  TH1F* h_exit_x = new TH1F("h_exit_x","",30000,-3.,3.);
  t -> Draw("simHits_exit_local_x>>h_exit_x","","goff");
  h_exit_x -> SetLineColor(kRed);
  h_exit_x -> Draw("same");
  latexLabel -> Draw("same");
  
  c1 -> cd(2);
  gPad -> SetLogy();
  TH1F* h_entry_y = new TH1F("h_entry_y","",30000,-3.,3.);
  t -> Draw("simHits_entry_local_y>>h_entry_y","","goff");
  h_entry_y -> SetTitle(";simHit entry/exit point y (cm);entries");
  h_entry_y -> Draw();
  TH1F* h_exit_y = new TH1F("h_exit_y","",30000,-3.,3.);
  t -> Draw("simHits_exit_local_y>>h_exit_y","","goff");
  h_exit_y -> SetLineColor(kRed);
  h_exit_y -> Draw("same");
  latexLabel -> Draw("same");
  
  c1 -> cd(3);
  gPad -> SetLogy();
  TH1F* h_entry_z = new TH1F("h_entry_z","",10000,-0.5,0.5);
  t -> Draw("simHits_entry_local_z>>h_entry_z","","goff");
  h_entry_z -> SetTitle(";simHit entry/exit point z (cm);entries");
  h_entry_z -> Draw();
  TH1F* h_exit_z = new TH1F("h_exit_z","",10000,-0.5,0.5);
  t -> Draw("simHits_exit_local_z>>h_exit_z","","goff");
  h_exit_z -> SetLineColor(kRed);
  h_exit_z -> Draw("same");
  latexLabel -> Draw("same");
  
  c1 -> Print(Form("c_%s_simHit_entryPoint.png",label.c_str()));
  
  

  TCanvas* c2 = new TCanvas("c2","c2",1000,500);
  c2 -> Divide(2,1);
  
  c2 -> cd(1);
  gPad -> SetLogy();
  TH1F* h_Dz = new TH1F("h_Dz","",500,-5.,5.);
  t -> Draw("matchedSimHits_track_Dz>>h_Dz","track_pt > 0.7 && track_isHighPurity == 1","goff");
  h_Dz -> SetTitle(";simHit-track #Deltaz (cm);entries");
  h_Dz -> Draw();
  h_Dz -> Fit("gaus","QLS+","",-1.,1.);
  TF1* f_gaus_Dz = (TF1*)( h_Dz->GetFunction("gaus") );
  TLatex* latex1 = new TLatex(0.20,0.80,Form("#sigma = %.2f mm",10.*f_gaus_Dz->GetParameter(2)));
  latex1 -> SetNDC();
  latex1 -> SetTextFont(42);
  latex1 -> SetTextColor(kRed);
  latex1 -> SetTextSize(0.03);
  latex1 -> Draw("same");
  latexLabel -> Draw("same");
  
  c2 -> cd(2);
  gPad -> SetLogy();
  TH1F* h_RDphi = new TH1F("h_RDphi","",500,-1.,1.);
  t -> Draw("matchedSimHits_track_RDphi>>h_RDphi","","goff");
  h_RDphi -> SetTitle(";simHit-track R#Delta#phi (cm);entries");
  h_RDphi -> Draw();
  h_RDphi -> Fit("gaus","QLS+","",-0.05,0.05);
  TF1* f_gaus_RDphi = (TF1*)( h_RDphi->GetFunction("gaus") );
  TLatex* latex2 = new TLatex(0.20,0.80,Form("#sigma = %.2f mm",10.*f_gaus_RDphi->GetParameter(2)));
  latex2 -> SetNDC();
  latex2 -> SetTextFont(42);
  latex2 -> SetTextColor(kRed);
  latex2 -> SetTextSize(0.03);
  latex2 -> Draw("same");
  latexLabel -> Draw("same");
  
  c2 -> Print(Form("c_%s_simHit_track.png",label.c_str()));
  
  
  
  c2 = new TCanvas("c2bis","c2bis",1000,500);
  c2 -> Divide(2,1);
  
  c2 -> cd(1);
  gPad -> SetLogy();
  h_Dz = new TH1F("h_Dz_bis","",500,-5.,5.);
  t -> Draw("matchedRecHits_track_Dz>>h_Dz_bis","track_pt > 0.7 && track_isHighPurity == 1","goff");
  h_Dz -> SetTitle(";recHit-track #Deltaz (cm);entries");
  h_Dz -> Draw();
  h_Dz -> Fit("gaus","QLS+","",-0.4,0.4);
  latexLabel -> Draw("same");
  
  c2 -> cd(2);
  gPad -> SetLogy();
  h_RDphi = new TH1F("h_RDphi_bis","",500,-5.,5.);
  t -> Draw("matchedRecHits_track_RDphi>>h_RDphi_bis","track_pt > 0.7 && track_isHighPurity == 1","goff");
  h_RDphi -> SetTitle(";recHit-track R#Delta#phi (cm);entries");
  h_RDphi -> Draw();
  latexLabel -> Draw("same");
  
  c2 -> Print(Form("c_%s_recHit_track.png",label.c_str()));


  
  TCanvas* c3 = new TCanvas("c3","c3");
  
  TH2F* h2_R_vs_iphi = new TH2F("h2_R_vs_iphi","",64,-0.5,63.5,1000,117.,118);
  t -> Draw("matchedSimHits_entry_global_R:((matchedSimHits_iphi-1)%64)>>h2_R_vs_iphi","","COLZ");
  // t -> Draw("matchedSimHits_entry_global_R:((matchedSimHits_iphi-1)%64)>>h2_R_vs_iphi","matchedSimHits_modType == 1","COLZ");
  // TH2F* h2_R_vs_iphi = new TH2F("h2_R_vs_iphi","",4,-0.5,3.5,200,117,118);
  // t -> Draw("recHits_global_R:((recHits_iphi-1)%16)>>h2_R_vs_iphi","","COLZ");
  h2_R_vs_iphi -> SetTitle(";simHit local iphi;simHit global R (cm)");
  c3 -> Print("c_simHit_global_R_vs_local_iphi.png");
  latexLabel -> Draw("same");

  c3 -> Print(Form("c_%s_simHit_R_vs_iphi.png",label.c_str()));
}
