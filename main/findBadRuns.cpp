#include "interface/TrackTree.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <map>
#include <string>
#include <cstdlib>
#include <stdlib.h>
#include <utility>

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TProfile2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TApplication.h"
#include "TFormula.h"

#include "TStyle.h"
#include "TSystem.h"
#include "TKey.h"

#include "TString.h"
#include "TTree.h"
#include "TBranch.h"

#include "TSpline.h"
#include "TCanvas.h"
#include "TObject.h"



int main(int argc, char** argv)
{
  if( argc < 3 )
  {
    std::cout << ">>> drawCTRSingles::usage:   " << argv[0] << " configFile.cfg runNumber" << std::endl;
    return -1;
  }
  
  
  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int iRun = atoi(argv[2]);
  
  int iRunMax = iRun;
  if( argc > 3 ) iRunMax = atoi(argv[3]);
  
  
  
  //--- get options
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir");
  int drawPlots = opts.GetOpt<int>("Input.drawPlots");
  int useFNALReco = opts.GetOpt<int>("Input.useFNALReco");
  
  
  
  // define channels
  const int NCH = 6;
  
  int ampch_id [NCH];
  if( useFNALReco )
  {
    ampch_id[0] = 10;
    ampch_id[1] = 13;
    ampch_id[2] = 11;
    ampch_id[3] = 14;
    ampch_id[4] = 12;
    ampch_id[5] = 15;
  }
  else
  {
    ampch_id[0] = 3;
    ampch_id[1] = 4;
    ampch_id[2] = 5;
    ampch_id[3] = 6;
    ampch_id[4] = 7;
    ampch_id[5] = 8;  
  }
  
  std::string chname[NCH/2];
  chname[0] = "top left";
  chname[1] = "mid left";
  chname[2] = "btm left";
  
  std::map<int,TLatex*> latex;
  for(int iCh = 0; iCh < NCH/2; ++iCh)
  {
    latex[iCh] = new TLatex(0.16,0.96,Form("%s",chname[iCh].c_str()));
    latex[iCh] -> SetNDC();
    latex[iCh] -> SetTextFont(42);
    latex[iCh] -> SetTextSize(0.03);
  }
  
  double cut_ampMin[NCH/2];
  cut_ampMin[0] = 8.;
  cut_ampMin[1] = 8.;
  cut_ampMin[2] = 8.;
  
  double cut_ampMax[NCH/2];
  cut_ampMax[0] = 1000.;
  cut_ampMax[1] = 1000.;
  cut_ampMax[2] = 1000.;
  
  
  
  {
    std::cout << "processing run " << iRun << std::endl;
    
    
    
    // define beamspot
    float BSX = 0.;
    float BSY = 0.;
    float centerX[NCH/2];
    float centerY[NCH/2];
    if( iRun >= 6369 && iRun <= 6401 )
    {
      centerX[0] = 9.; centerY[0] = 30.;
      centerX[1] = 9.; centerY[1] = 26.;
      centerX[2] = 9.; centerY[2] = 22.;
      
      BSX = 6.;
      BSY = 1.5;

      cut_ampMin[0] = 100.;
      cut_ampMin[1] = 100.;
      cut_ampMin[2] = 100.;
    }
    if( iRun >= 6402 && iRun <= 6440 )
    {
      centerX[0] = 9.; centerY[0] = 30.;
      centerX[1] = 9.; centerY[1] = 26.;
      centerX[2] = 9.; centerY[2] = 22.;
      
      BSX = 6.;
      BSY = 1.5;
      
      cut_ampMin[0] = 65.;
      cut_ampMin[1] = 65.;
      cut_ampMin[2] = 65.;
    }
    if( iRun >= 6442 && iRun <= 6485 )
    {
      centerX[0] = 9.; centerY[0] = 30.;
      centerX[1] = 9.; centerY[1] = 26.;
      centerX[2] = 9.; centerY[2] = 22.;
      
      BSX = 6.;
      BSY = 1.5;
      
      cut_ampMin[0] = 40.;
      cut_ampMin[1] = 40.;
      cut_ampMin[2] = 40.;
    }
    if( iRun >= 6486 && iRun <= 6520 )
    {
      centerX[0] = 9.; centerY[0] = 30.;
      centerX[1] = 9.; centerY[1] = 26.;
      centerX[2] = 9.; centerY[2] = 22.;
      
      BSX = 6.;
      BSY = 1.5;
      
      cut_ampMin[0] = 25.;
      cut_ampMin[1] = 25.;
      cut_ampMin[2] = 25.;
    }
    if( iRun >= 6521 && iRun <= 6522 )
    {
      centerX[0] = 9.; centerY[0] = 30.;
      centerX[1] = 9.; centerY[1] = 26.;
      centerX[2] = 9.; centerY[2] = 22.;
      
      BSX = 6.;
      BSY = 1.5;
      
      cut_ampMin[0] = 13.;
      cut_ampMin[1] = 13.;
      cut_ampMin[2] = 13.;
    }
    if( iRun >= 6523 && iRun <= 6710 )
    {
      centerX[0] = 9.; centerY[0] = 30.;
      centerX[1] = 9.; centerY[1] = 26.;
      centerX[2] = 9.; centerY[2] = 22.;
      
      BSX = 6.;
      BSY = 1.5;
      
      cut_ampMin[0] = 100.;
      cut_ampMin[1] = 100.;
      cut_ampMin[2] = 100.;
    }
    if( iRun >= 6754 && iRun <= 6866 )
    {
      centerX[0] = 12.; centerY[0] = 29.5;
      centerX[1] = 12.; centerY[1] = 25.5;
      centerX[2] = 12.; centerY[2] = 22.;
      
      BSX = 15.;
      BSY = 1.5;
      
      cut_ampMin[0] = 90.;
      cut_ampMin[1] = 90.;
      cut_ampMin[2] = 90.;
    }
    if( iRun >= 6867 && iRun <= 6913 )
    {
      centerX[0] = 14.; centerY[0] = 29.5;
      centerX[1] = 14.; centerY[1] = 25.5;
      centerX[2] = 14.; centerY[2] = 22.;
      
      BSX = 12.;
      BSY = 1.5;
      
      cut_ampMin[0] = 90.;
      cut_ampMin[1] = 90.;
      cut_ampMin[2] = 90.;
    }
    if( iRun >= 6914 && iRun <= 6951 )
    {
      centerX[0] = 14.; centerY[0] = 29.5;
      centerX[1] = 14.; centerY[1] = 25.5;
      centerX[2] = 14.; centerY[2] = 22.;
      
      BSX = 12.;
      BSY = 1.5;
      
      cut_ampMin[0] = 60.;
      cut_ampMin[1] = 60.;
      cut_ampMin[2] = 60.;
    }
    if( iRun >= 6952 && iRun <= 7011 )
    {
      centerX[0] = 14.; centerY[0] = 29.5;
      centerX[1] = 14.; centerY[1] = 25.5;
      centerX[2] = 14.; centerY[2] = 22.;
      
      BSX = 12.;
      BSY = 1.5;
      
      cut_ampMin[0] = 35.;
      cut_ampMin[1] = 35.;
      cut_ampMin[2] = 35.;
    }
    if( iRun >= 7012 && iRun <= 7058 )
    {
      centerX[0] = 14.; centerY[0] = 29.5;
      centerX[1] = 14.; centerY[1] = 25.5;
      centerX[2] = 14.; centerY[2] = 22.;
      
      BSX = 12.;
      BSY = 1.5;
      
      cut_ampMin[0] = 20.;
      cut_ampMin[1] = 20.;
      cut_ampMin[2] = 20.;
    }
    if( iRun >= 7059 && iRun <= 7102 )
    {
      centerX[0] = 14.; centerY[0] = 29.5;
      centerX[1] = 14.; centerY[1] = 25.5;
      centerX[2] = 14.; centerY[2] = 22.;
      
      BSX = 12.;
      BSY = 1.5;
      
      cut_ampMin[0] = 10.;
      cut_ampMin[1] = 10.;
      cut_ampMin[2] = 10.;
    }
    if( iRun >= 7111 && iRun <= 7135 )
    {
      centerX[0] = 14.; centerY[0] = 29.5;
      centerX[1] = 14.; centerY[1] = 25.5;
      centerX[2] = 14.; centerY[2] = 22.;
      
      BSX = 12.;
      BSY = 1.5;
      
      cut_ampMin[0] = 8.;
      cut_ampMin[1] = 8.;
      cut_ampMin[2] = 8.;
    }
    if( iRun >= 7136 && iRun <= 7208 )
    {
      centerX[0] = 14.; centerY[0] = 29.5;
      centerX[1] = 14.; centerY[1] = 25.5;
      centerX[2] = 14.; centerY[2] = 22.;
      
      BSX = 12.;
      BSY = 1.5;
      
      cut_ampMin[0] = 8.;
      cut_ampMin[1] = 8.;
      cut_ampMin[2] = 8.;
    }
    if( iRun >= 7211 && iRun <= 7335 )
    {
      centerX[0] = 14.; centerY[0] = 29.5;
      centerX[1] = 14.; centerY[1] = 25.5;
      centerX[2] = 14.; centerY[2] = 22.;
      
      BSX = 12.;
      BSY = 1.5;
      
      cut_ampMin[0] = 8.;
      cut_ampMin[1] = 8.;
      cut_ampMin[2] = 8.;
    }
    if( iRun >= 7336 && iRun <= 7403 )
    {
      centerX[0] = 17.; centerY[0] = 29.5;
      centerX[1] = 17.; centerY[1] = 25.5;
      centerX[2] = 17.; centerY[2] = 22.;
      
      BSX = 12.;
      BSY = 1.5;
      
      cut_ampMin[0] = 100.;
      cut_ampMin[1] = 100.;
      cut_ampMin[2] = 100.;
    }
    if( iRun >= 7404 && iRun <= 7497 )
    {
      centerX[0] = 17.; centerY[0] = 29.5;
      centerX[1] = 17.; centerY[1] = 25.5;
      centerX[2] = 17.; centerY[2] = 22.;
      
      BSX = 12.;
      BSY = 1.5;
      
      cut_ampMin[0] = 100.;
      cut_ampMin[1] = 100.;
      cut_ampMin[2] = 100.;
    }
    if( iRun >= 7498 && iRun <= 7516 )
    {
      centerX[0] = 14.; centerY[0] = 29.5;
      centerX[1] = 14.; centerY[1] = 25.5;
      centerX[2] = 14.; centerY[2] = 22.;
      
      BSX = 13.;
      BSY = 1.5;
      
      cut_ampMin[0] = 90.;
      cut_ampMin[1] = 90.;
      cut_ampMin[2] = 90.;
    }
    if( iRun >= 7518 && iRun <= 7553 )
    {
      centerX[0] = 14.; centerY[0] = 29.5;
      centerX[1] = 14.; centerY[1] = 25.5;
      centerX[2] = 14.; centerY[2] = 22.;
      
      BSX = 13.;
      BSY = 1.5;
      
      cut_ampMin[0] = 20.;
      cut_ampMin[1] = 20.;
      cut_ampMin[2] = 20.;
    }
    if( iRun >= 7555 && iRun <= 7615 )
    {
      centerX[0] = 17.; centerY[0] = 29.5;
      centerX[1] = 17.; centerY[1] = 25.5;
      centerX[2] = 17.; centerY[2] = 22.;
      
      BSX = 12.;
      BSY = 1.5;
      
      cut_ampMin[0] = 25.;
      cut_ampMin[1] = 25.;
      cut_ampMin[2] = 25.;
    }
    if( iRun >= 7629 && iRun <= 7638 )
    {
      centerX[0] = 14.; centerY[0] = 29.5;
      centerX[1] = 14.; centerY[1] = 25.5;
      centerX[2] = 14.; centerY[2] = 22.;
      
      BSX = 13.;
      BSY = 1.5;
      
      cut_ampMin[0] = 20.;
      cut_ampMin[1] = 20.;
      cut_ampMin[2] = 20.;
    }
    if( iRun >= 7639 && iRun <= 7679 )
    {
      centerX[0] = 10.; centerY[0] = 28.5;
      centerX[1] = 10.; centerY[1] = 24.5;
      centerX[2] = 10.; centerY[2] = 21.;
      
      BSX = 14.;
      BSY = 1.5;
      
      cut_ampMin[0] = 90.;
      cut_ampMin[1] = 90.;
      cut_ampMin[2] = 90.;
    }
    if( iRun >= 7682 && iRun <= 7800 )
    {
      centerX[0] = 14.; centerY[0] = 28.5;
      centerX[1] = 14.; centerY[1] = 24.5;
      centerX[2] = 14.; centerY[2] = 21.;
      
      BSX = 13.;
      BSY = 1.5;
      
      cut_ampMin[0] = 120.;
      cut_ampMin[1] = 120.;
      cut_ampMin[2] = 120.;
    }
    if( iRun >= 7802 && iRun <= 7845 )
    {
      centerX[0] = 14.; centerY[0] = 28.5;
      centerX[1] = 14.; centerY[1] = 24.5;
      centerX[2] = 14.; centerY[2] = 21.;
      
      BSX = 13.;
      BSY = 1.5;
      
      cut_ampMin[0] = 25.;
      cut_ampMin[1] = 25.;
      cut_ampMin[2] = 25.;
    }
    if( iRun >= 7850 && iRun <= 7877 )
    {
      centerX[0] = 10.; centerY[0] = 28.5;
      centerX[1] = 10.; centerY[1] = 24.5;
      centerX[2] = 10.; centerY[2] = 21.;
      
      BSX = 14.;
      BSY = 1.5;
      
      cut_ampMin[0] = 90.;
      cut_ampMin[1] = 90.;
      cut_ampMin[2] = 90.;
    }
    if( iRun >= 7881 && iRun <= 7953 )
    {
      centerX[0] = 16.; centerY[0] = 28.5;
      centerX[1] = 16.; centerY[1] = 24.5;
      centerX[2] = 16.; centerY[2] = 21.;
      
      BSX = 11.;
      BSY = 1.5;
      
      cut_ampMin[0] = 150.;
      cut_ampMin[1] = 150.;
      cut_ampMin[2] = 150.;
    }
    if( iRun >= 8029 && iRun <= 8062 )
    {
      centerX[0] = 18.; centerY[0] = 29.;
      centerX[1] = 18.; centerY[1] = 25.;
      centerX[2] = 18.; centerY[2] = 21.5;
      
      BSX = 10.;
      BSY = 1.5;
      
      cut_ampMin[0] = 150.;
      cut_ampMin[1] = 150.;
      cut_ampMin[2] = 150.;
    }
    if( iRun >= 8064 && iRun <= 8088 )
    {
      centerX[0] = 18.; centerY[0] = 29.;
      centerX[1] = 18.; centerY[1] = 25.;
      centerX[2] = 18.; centerY[2] = 21.5;
      
      BSX = 10.;
      BSY = 1.5;
      
      cut_ampMin[0] = 40.;
      cut_ampMin[1] = 40.;
      cut_ampMin[2] = 40.;
    }
    if( iRun >= 8113 && iRun <= 8164 )
    {
      centerX[0] = 18.; centerY[0] = 29.;
      centerX[1] = 18.; centerY[1] = 25.;
      centerX[2] = 18.; centerY[2] = 21.5;
      
      BSX = 10.;
      BSY = 1.5;
      
      cut_ampMin[0] = 40.;
      cut_ampMin[1] = 40.;
      cut_ampMin[2] = 40.;
    }
    if( iRun >= 8173 && iRun <= 8261 )
    {
      centerX[0] = 12.; centerY[0] = 29.;
      centerX[1] = 12.; centerY[1] = 25.5;
      centerX[2] = 12.; centerY[2] = 22.;
      
      BSX = 9.;
      BSY = 1.5;
      
      cut_ampMin[0] = 200.;
      cut_ampMin[1] = 200.;
      cut_ampMin[2] = 200.;
    }
    if( iRun >= 8262 && iRun <= 8298 )
    {
      centerX[0] = 12.; centerY[0] = 29.;
      centerX[1] = 12.; centerY[1] = 25.5;
      centerX[2] = 12.; centerY[2] = 22.;
      
      BSX = 9.;
      BSY = 1.5;
      
      cut_ampMin[0] = 60.;
      cut_ampMin[1] = 60.;
      cut_ampMin[2] = 60.;
    }
    if( iRun >= 8299 && iRun <= 8334 )
    {
      centerX[0] = 12.; centerY[0] = 29.;
      centerX[1] = 12.; centerY[1] = 25.5;
      centerX[2] = 12.; centerY[2] = 22.;
      
      BSX = 7.;
      BSY = 1.5;
      
      cut_ampMin[0] = 80.;
      cut_ampMin[1] = 80.;
      cut_ampMin[2] = 80.;
    }
    if( iRun >= 8335 && iRun <= 8443 )
    {
      centerX[0] = 12.; centerY[0] = 29.;
      centerX[1] = 12.; centerY[1] = 25.5;
      centerX[2] = 12.; centerY[2] = 22.;
      
      BSX = 7.;
      BSY = 1.5;
      
      cut_ampMin[0] = 250.;
      cut_ampMin[1] = 250.;
      cut_ampMin[2] = 250.;
    }
    if( iRun >= 8461 && iRun <= 8555 )
    {
      centerX[0] = 14.; centerY[0] = 28.5;
      centerX[1] = 14.; centerY[1] = 25.;
      centerX[2] = 14.; centerY[2] = 21.5;
      
      BSX = 13.;
      BSY = 1.5;
      
      cut_ampMin[0] = 110.;
      cut_ampMin[1] = 110.;
      cut_ampMin[2] = 110.;
    }
    if( iRun >= 8560 && iRun <= 8627 )
    {
      centerX[0] = 15.; centerY[0] = 28.5;
      centerX[1] = 15.; centerY[1] = 25.;
      centerX[2] = 15.; centerY[2] = 21.5;
      
      BSX = 11.;
      BSY = 1.5;
      
      cut_ampMin[0] = 150.;
      cut_ampMin[1] = 150.;
      cut_ampMin[2] = 150.;
    }
    if( iRun >= 8629 && iRun <= 8922 )
    {
      centerX[0] = 10.; centerY[0] = 27.5;
      centerX[1] = 10.; centerY[1] = 24.5;
      centerX[2] = 10.; centerY[2] = 22.5;
      
      BSX = 14.;
      BSY = 1.5;
      
      cut_ampMin[0] = 10.;
      cut_ampMin[1] = 10.;
      cut_ampMin[2] = 10.;
    }
    if( iRun >= 8943 && iRun <= 9068 )
    {
      centerX[0] = 10.; centerY[0] = 28.5;
      centerX[1] = 10.; centerY[1] = 25.;
      centerX[2] = 10.; centerY[2] = 21.5;
      
      BSX = 14.;
      BSY = 1.5;
      
      cut_ampMin[0] = 50.;
      cut_ampMin[1] = 50.;
      cut_ampMin[2] = 50.;
    }
    if( iRun >= 9069 && iRun <= 9362 )
    {
      centerX[0] = 10.; centerY[0] = 26.5;
      centerX[1] = 10.; centerY[1] = 25.;
      centerX[2] = 10.; centerY[2] = 23.;
      
      BSX = 14.;
      BSY = 1.5;
      
      cut_ampMin[0] = 10.;
      cut_ampMin[1] = 20.;
      cut_ampMin[2] = 10.;
    }
    if( iRun >= 9363 && iRun <= 9651 )
    {
      centerX[0] = 10.; centerY[0] = 28.;
      centerX[1] = 10.; centerY[1] = 25.;
      centerX[2] = 10.; centerY[2] = 22.;
      
      BSX = 14.;
      BSY = 1.5;
      
      cut_ampMin[0] = 10.;
      cut_ampMin[1] = 20.;
      cut_ampMin[2] = 10.;
    }
    if( iRun >= 9652 && iRun <= 10147 )
    {
      centerX[0] = 10.; centerY[0] = 28.5;
      centerX[1] = 10.; centerY[1] = 25.;
      centerX[2] = 10.; centerY[2] = 21.5;
      
      BSX = 14.;
      BSY = 1.5;
      
      cut_ampMin[0] = 50.;
      cut_ampMin[1] = 50.;
      cut_ampMin[2] = 50.;
    }
    if( iRun >= 10161 && iRun <= 10205 )
    {
      centerX[0] = 15.; centerY[0] = 28.5;
      centerX[1] = 15.; centerY[1] = 25.;
      centerX[2] = 15.; centerY[2] = 21.5;
      
      BSX = 11.;
      BSY = 1.5;
      
      cut_ampMin[0] = 80.;
      cut_ampMin[1] = 80.;
      cut_ampMin[2] = 80.;
    }
    if( iRun >= 10411 && iRun <= 10546 )
    {
      centerX[0] = 9.; centerY[0] = 30.5;
      centerX[1] = 9.; centerY[1] = 27.;
      centerX[2] = 9.; centerY[2] = 23.5;
      
      BSX = 11.;
      BSY = 1.5;
      
      cut_ampMin[0] = 80.;
      cut_ampMin[1] = 80.;
      cut_ampMin[2] = 80.;
    }
    if( iRun >= 10644 && iRun <= 10839 )
    {
      centerX[0] = 13.; centerY[0] = 34.;
      centerX[1] = 13.; centerY[1] = 30.5;
      centerX[2] = 13.; centerY[2] = 27.;
      
      BSX = 11.;
      BSY = 1.5;
      
      cut_ampMin[0] = 100.;
      cut_ampMin[1] = 100.;
      cut_ampMin[2] = 100.;
    }
    // else if( iRun < 7498 )
    // {
    //   centerX[0] = 15.; centerY[0] = 29.;
    //   centerX[1] = 15.; centerY[1] = 25.5;
    //   centerX[2] = 15.; centerY[2] = 22.;
    // }
    // else if( iRun < 7553 )
    // {
    //   centerX[0] = 10.; centerY[0] = 29.;
    //   centerX[1] = 10.; centerY[1] = 25.5;
    //   centerX[2] = 10.; centerY[2] = 22.;
    // }
    // else if( iRun < 7616 )
    // {
    //   centerX[0] = 15.; centerY[0] = 29.;
    //   centerX[1] = 15.; centerY[1] = 25.5;
    //   centerX[2] = 15.; centerY[2] = 22.;
    // }
    // else if( iRun < 7616 )
    // {
    //   centerX[0] = 15.; centerY[0] = 29.;
    //   centerX[1] = 15.; centerY[1] = 25.5;
    //   centerX[2] = 15.; centerY[2] = 22.;
    // }
    // else if( iRun < 10410 )
    // {
    //   centerX[0] = 0.; centerY[0] = 29.;
    //   centerX[1] = 0.; centerY[1] = 25.5;
    //   centerX[2] = 0.; centerY[2] = 22.;
    // }
    // else if( iRun < 10569 )
    // {
    //   centerX[0] = 2.; centerY[0] = 31.;
    //   centerX[1] = 2.; centerY[1] = 27.;
    //   centerX[2] = 2.; centerY[2] = 24.;
    // }
    // else if( iRun < 10640 )
    // {
    //   centerX[0] = 7.; centerY[0] = 31.;
    //   centerX[1] = 7.; centerY[1] = 27.;
    //   centerX[2] = 7.; centerY[2] = 24.;
    // }
    
    
    // get tree
    TChain* myTree;
    
    // myTree -> Add(Form("%s/*7069*.root",folder.c_str()));
    // myTree -> Add(Form("%s/*7070*.root",folder.c_str()));
    // myTree -> Add(Form("%s/*7071*.root",folder.c_str()));
    // myTree -> Add(Form("%s/*7072*.root",folder.c_str()));
    // myTree -> Add(Form("%s/*7073*.root",folder.c_str()));
    // myTree -> Add(Form("%s/*7074*.root",folder.c_str()));
    // myTree -> Add(Form("%s/*7075*.root",folder.c_str()));
    // myTree -> Add(Form("%s/*7076*.root",folder.c_str()));
    // myTree -> Add(Form("%s/*7077*.root",folder.c_str()));
    // myTree -> Add(Form("%s/*7078*.root",folder.c_str()));
    // myTree -> Add(Form("%s/*7079*.root",folder.c_str()));
    // myTree -> Add(Form("%s/*7080*.root",folder.c_str()));
    // myTree -> Add(Form("%s/*7081*.root",folder.c_str()));
    
    if( useFNALReco )
    {
      myTree = new TChain("pulse","pulse");
      // std::cout << ">>> opening file " << Form("%s/RawDataSaver0CMSVMETiming_Run%d_0_Raw.root",inputDir.c_str(),iRun) << std::endl;
      // TFile* inFile = TFile::Open(Form("%s/RawDataSaver0CMSVMETiming_Run%d_0_Raw.root",inputDir.c_str(),iRun));
      // if( !inFile ) return -1;
      // myTree = (TTree*)( inFile->Get("pulse") );
      // if( myTree == 0 ) return -1;
      // std::cout << ">>> got tree" << std::endl;
    }
    else
    {
      myTree = new TChain("h4","h4");
      // std::cout << ">>> opening file " << Form("%s/%d.root",inputDir.c_str(),iRun) << std::endl;
      // TFile* inFile = TFile::Open(Form("%s/%d.root",inputDir.c_str(),iRun));
      // if( !inFile ) return -1;
      // myTree = (TTree*)( inFile->Get("h4") );
      // if( myTree == 0 ) return -1;
      // std::cout << ">>> got tree" << std::endl;
      for(int jj = iRun; jj <= iRunMax; ++jj)
        myTree -> Add(Form("%s/%d.root",inputDir.c_str(),jj));
    }
    
    float  amp[100];
    float x_dut[4];
    float y_dut[4];
    int nTracks;
    std::vector<TrackPar>* tracks = new std::vector<TrackPar>;    
    
    if( useFNALReco )
    {
      myTree -> SetBranchStatus("*",0);
      myTree -> SetBranchStatus("amp",   1); myTree -> SetBranchAddress("amp",     amp);
      myTree -> SetBranchStatus("x_dut", 1); myTree -> SetBranchAddress("x_dut", x_dut);
      myTree -> SetBranchStatus("y_dut", 1); myTree -> SetBranchAddress("y_dut", y_dut);
    }
    else
    {
      myTree -> SetBranchStatus("time",0);
      myTree -> SetBranchStatus("time_chi2",0);
      myTree -> SetBranchStatus("time_error",0);
      myTree -> SetBranchStatus("time_slope",0);
      myTree -> SetBranchStatus("pedestal",0);
      myTree -> SetBranchStatus("b_charge",0);
      myTree -> SetBranchStatus("b_slope",0);
      myTree -> SetBranchStatus("amp_max",   1); myTree -> SetBranchAddress("amp_max",        amp);
      myTree -> SetBranchStatus("n_tracks",  1); myTree -> SetBranchAddress("n_tracks",  &nTracks);
      myTree -> SetBranchStatus("fitResult", 1); myTree -> SetBranchAddress("fitResult",  &tracks);
    }
    
    
    
    // define histograms
    float minX = -10, minY =   0;
    float maxX = +40, maxY = +50;
    
    TH2F hBeamXY("hBeamXY", "hBeamXY", 100, minX, maxX, 100, minY, maxY);
    
    std::map<int,TH1F> h_amp;
    std::map<int,TProfile2D> p2_XY_amp;
    std::map<int,TProfile2D> p2_XY_eff;
    for(int iCh = 0; iCh<NCH/2; ++iCh)
    {
      h_amp[iCh] = TH1F(Form("h_amp_%s",chname[iCh].c_str()),"",1000,0.,250.);
      
      p2_XY_amp[iCh] = TProfile2D(Form("p2_XY_amp_%s",chname[iCh].c_str()),"", 100,minX,maxX, 100,minY,maxY);
      p2_XY_eff[iCh] = TProfile2D(Form("p2_XY_eff_%s",chname[iCh].c_str()),"", 100,minX,maxX, 100,minY,maxY);
    }
    
    
    
    // loop over entries
    int counter = 0;  
    int efficiency = 0;

    int nEntries = myTree->GetEntries();  
    if( nEntries == 0 ) return -1;
    std::cout << ">>> read " << nEntries << " entries" << std::endl;
    for(int entry = 0; entry < nEntries; ++entry)
    {
      myTree -> GetEntry(entry);    
      
      float myX;
      float myY;
      if( useFNALReco )
      {
        myX = x_dut[0];
        myY = y_dut[0];
      }
      else
      {
        if( tracks->size() >= 1 )
        {
          myX = (tracks->at(0)).x();
          myY = (tracks->at(0)).y();
        }
      }
      
      if (myX == -999 || myY == -999) continue;
      
      hBeamXY.Fill(myX, myY);
      
      for(int iCh = 0; iCh < NCH/2; ++iCh)
      {
        float amp1 = amp[ampch_id[iCh*2]];
        float amp2 = amp[ampch_id[iCh*2+1]];
        if( !useFNALReco )
        {
          amp1 = amp1/4096.*1000.;
          amp2 = amp2/4096.*1000.;
        }
        
        h_amp[iCh].Fill( 0.5*(amp1+amp2) );    
        p2_XY_amp[iCh].Fill(myX, myY, 0.5*(amp1+amp2) );
        
        //efficiency plots
        if( 0.5*(amp1+amp2) > cut_ampMin[iCh] )
          p2_XY_eff[iCh].Fill(myX, myY, 1);
        else
          p2_XY_eff[iCh].Fill(myX, myY, 0 );
        
        if( fabs(myX-centerX[iCh]) > BSX ) continue;
        if( fabs(myY-centerY[iCh]) > BSY ) continue;
        
        //efficiency plots
        if( 0.5*(amp1+amp2) > cut_ampMin[iCh] )
        {
          p2_XY_eff[iCh].Fill(myX, myY, 1);
          ++counter;
          ++efficiency;
        }
        else
        {
          p2_XY_eff[iCh].Fill(myX, myY, 0 );
          ++counter;
        }
      }
    }
    
    if( counter == 0 ) return -1;
    
    
    std::string plotDir;
    if( iRunMax == iRun ) plotDir = std::string(Form("/afs/cern.ch/user/a/abenagli/www/TIMING/TBatFNALApr2019/goodRuns/all/%d",iRun));
    else                  plotDir = std::string(Form("/afs/cern.ch/user/a/abenagli/www/TIMING/TBatFNALApr2019/goodRuns/%d-%d",iRun,iRunMax));
    system(Form("mkdir -p %s",plotDir.c_str()));
    
    std::string outFileName(Form("%s/efficiency.txt",plotDir.c_str()));
    FILE* outFile = fopen(outFileName.c_str(),"w");
    fprintf(outFile, "%d %f \n",iRun,1.*efficiency/counter);
    std::cout << ">>> Run: " << iRun << " --> Events read: " << counter << " --> selectedEntries: " << efficiency << " :: EFFICIENCY = " << 1.*efficiency/counter << std::endl;
    fclose(outFile);
    
    
    if( drawPlots )
    {
      system(Form("cp /afs/cern.ch/user/a/abenagli/www/TIMING/TBatFNALApr2019/goodRuns/index.php %s",plotDir.c_str()));
      
      TCanvas* cXYAmps = new TCanvas("cXYAmps","cXYAmps",800,1600);
      cXYAmps -> Divide(1,NCH/2);
      int c_id = 1;
      for (int iCh = 0; iCh < NCH/2; ++iCh)
      {
        TLine * x_min = new TLine(cut_ampMin[iCh], minY, cut_ampMin[iCh], maxY);
        TLine * x_max = new TLine(cut_ampMax[iCh], minY, cut_ampMax[iCh], maxY);
        x_min->SetLineColor(kRed);
        x_min->SetLineWidth(1);
        x_max->SetLineColor(kRed);
        x_max->SetLineWidth(1);
        
        cXYAmps -> cd(c_id);
        h_amp[iCh].SetFillColor(kBlue);
        h_amp[iCh].Draw();
        gPad -> SetLogy();
        c_id++;
        
        x_min -> Draw("same");
        x_max -> Draw("same");
        
        latex[iCh] -> Draw("same");
      }
      cXYAmps -> Print(Form("%s/cXYAmps.png",plotDir.c_str()));
      delete cXYAmps;
      
      TCanvas* cXYEffs = new TCanvas("cXYEffs","cXYEffs",800,1600);
      cXYEffs -> Divide(1,NCH/2);
      c_id = 1;
      for (int iCh = 0; iCh < NCH/2; ++iCh)
      {
        TLine * x_bs_min = new TLine(centerX[iCh]-BSX, minY, centerX[iCh]-BSX, maxY);
        TLine * x_bs_max = new TLine(centerX[iCh]+BSX, minY, centerX[iCh]+BSX, maxY);
        TLine * y_bs_min = new TLine(minX, centerY[iCh]-BSY, maxX, centerY[iCh]-BSY);
        TLine * y_bs_max = new TLine(minX, centerY[iCh]+BSY, maxX, centerY[iCh]+BSY);
        
        x_bs_min->SetLineColor(kRed);
        x_bs_min->SetLineWidth(1);
        x_bs_max->SetLineColor(kRed);
        x_bs_max->SetLineWidth(1);
        y_bs_min->SetLineColor(kRed);
        y_bs_min->SetLineWidth(1);
        y_bs_max->SetLineColor(kRed);
        y_bs_max->SetLineWidth(1);
        
        cXYEffs -> cd(c_id);
        p2_XY_eff[iCh].Draw("COLZ");
        p2_XY_eff[iCh].SetStats(0);
        p2_XY_eff[iCh].GetXaxis()->SetTitle("X [mm]");
        p2_XY_eff[iCh].GetYaxis()->SetTitle("Y [mm]");
        gPad -> SetLogz();
        c_id++;
        
        x_bs_min->Draw("same");
        x_bs_max->Draw("same");
        y_bs_min->Draw("same");
        y_bs_max->Draw("same");                
        
        latex[iCh] -> Draw("same");
      }
      cXYEffs -> Print(Form("%s/cXYEffs.png",plotDir.c_str()));
      delete cXYEffs;
      
      
      // TCanvas* cXYEffs = new TCanvas("cXYEffs","cXYEffs",1500,800);
      // cXYEffs -> Divide(2,2);
      // c_id = 1;
      // for (int iCh = 1; iCh < NCH; iCh++)
      // {
      //   cXYEffs -> cd(c_id);
      //   p2_XY_eff[iCh]->Draw("COLZ");
      //   p2_XY_eff[iCh]->SetStats(0);
      //   p2_XY_eff[iCh] -> GetXaxis()->SetTitle("X [mm]");
      //   p2_XY_eff[iCh] -> GetYaxis()->SetTitle("Y [mm]");
      //   gPad -> SetLogz();
      //   c_id++;
      
      //   x_bs_min->Draw("same");
      //   x_bs_max->Draw("same");
      //   y_bs_min->Draw("same");
      //   y_bs_max->Draw("same");                
      // }
    }
  }
}
