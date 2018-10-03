#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>

#include "TTree.h"



/*** tree variables ***/
struct TreeVars
{
  int event;
  int spill;
  int run;
  
  float* time;
  float* amp_max;
  int* nFibresOnX;
  int* nFibresOnY;
  float* hodoX;
  float* hodoY;
  float* wireX;
  float* wireY;
  
  std::map<std::string,int> channelIds;
  std::map<std::string,int> timeMethodIds;
  std::map<std::string,float> VbiasVals;
  std::map<std::string,float> NINOthrVals;
  
  float beamX;
  float beamY;
};

struct AnalysisVars
{
  float amp;
  float time;
  float time2;
  float rt;
  float dist;
  
  float timeRef;
};

/*** initialize tree variables ***/
void InitTreeVars(TTree* tree, TreeVars& treeVars, CfgManager& opts);

/*** reconstruct beam position ***/
void ReconstructHodoPosition(TreeVars& tv, CfgManager& opts, const std::string& ch, const bool& tracking, const int& plane);

/*** default event selection ***/
bool AcceptEvent(AnalysisVars& av, TreeVars& tv, CfgManager& opts, const std::string& ch, const int& nFibresMin, const int& nFibresMax, const float& CTRMin, const float& CTRMax);

