#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <string>
#include <map>

#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TPad.h"
#include "TLatex.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLine.h"


/*** find effective sigma ***/
void FindSmallestInterval(float* ret, TH1F* histo, const float& fraction);

/*** draw time resolution plot ***/
struct CTRResult
{
  float effSigma;
  float gausSigma;
  float gausSigmaErr;
};

CTRResult drawCTRPlot(TH1F* histo, const std::string& label, const int& rebin, const int& isMCP0, const int& isMCP1, const float& MCPIntrinsic,
                      const std::string& title, TLatex* latexLabel,
                      TH1F* histo_center, const std::string& center_label, TH1F* histo_border, const std::string& border_label);

/*** print canvas ***/
void PrintCanvas(TCanvas* c, CfgManager& opts, const std::string& ch, const std::string& plotDir, const std::string& label);


/*** draw 1D histogram ***/
void DrawHistogram(CfgManager& opts,
                   const std::string& ch, TH1F* histo,
                   const std::string& title,
                   const float& xMin, const float& xMax, const int& rebin, const bool& logy,
                   TLatex* latexLabel, std::vector<TLine*>* lines = NULL,
                   const bool& landauFit = false, const float& fitXMin = 0., const float& fitXMax = 0.);

void DrawHistogram(CfgManager& opts,
                   const std::string& ch, std::vector<TH1F*>& histos,
                   const std::string& title,
                   const float& xMin, const float& xMax, const int& rebin, const bool& logy,
                   TLatex* latexLabel, std::vector<TLine*>* lines = NULL,
                   const bool& landauFit = false, const float& fitXMin = 0., const float& fitXMax = 0.);

/*** draw 2D histogram ***/
void DrawHistogram2D(CfgManager& opts,
                     const std::string& ch, TH2F* histo2,
                     const std::string& title,
                     const float& xMin, const float& xMax, const float& yMin, const float& yMax, const float& zMin, const float& zMax,
                     const bool& logx, const bool& logy, const std::string& drawOpt,
                     TLatex* latexLabel);

/*** draw 1D profile ***/
void DrawProfile(CfgManager& opts,
                 const std::string& ch, TProfile* prof,
                 const std::string& title,
                 const float& xMin, const float& xMax, const float& yMin, const float& yMax,
                 const int& color, const std::string& drawOpt,
                 TLatex* latexLabel = NULL, TF1* func = NULL, const float& x0 = -999);

void DrawProfile(CfgManager& opts,
                 const std::string& ch, std::vector<TProfile*> profs,
                 const std::string& title,
                 const float& xMin, const float& xMax, const float& yMin, const float& yMax,
                 const std::vector<int>& colors, const std::string& drawOpt,
                 TLatex* latexLabel = NULL);

/*** draw 2D profile ***/
void DrawProfile2D(CfgManager& opts,
                   const std::string& ch, TProfile2D* prof2,
                   const std::string& title,
                   const float& xMin, const float& xMax, const float& yMin, const float& yMax, const float& zMin, const float& zMax,
                   TLatex* latexLabel = NULL, std::vector<TLine*>* lines = NULL);

/*** draw 1D graph ***/
void DrawGraph(CfgManager& opts,
               const std::string& ch, TGraph* graph,
               const std::string& title,
               const float& xMin, const float& xMax, const float& yMin, const float& yMax,
               const int& color, const std::string& drawOpt,
               TLatex* latexLabel = NULL);

/*** draw 1D graph ***/
void DrawGraph(CfgManager& opts,
               const std::string& ch, std::vector<TGraph*> graphs,
               const std::string& title,
               const float& xMin, const float& xMax, const float& yMin, const float& yMax,
               const std::vector<int>& colors, const std::string& drawOpt,
               TLatex* latexLabel = NULL);
