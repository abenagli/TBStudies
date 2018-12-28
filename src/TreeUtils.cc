#include "interface/TreeUtils.h"



void InitTreeVars(TTree* tree, TreeVars& treeVars, CfgManager& opts)
{
  treeVars.time    = new float[1000];
  treeVars.amp_max = new float[1000];
  treeVars.nFibresOnX = new int[2];
  treeVars.nFibresOnY = new int[2];
  treeVars.hodoX = new float[2];
  treeVars.hodoY = new float[2];
  treeVars.wireX = new float[2];
  treeVars.wireY = new float[2];
  
  treeVars.tracks = new std::vector<TrackPar>;

  tree -> SetBranchAddress("event",&treeVars.event);
  tree -> SetBranchAddress("spill",&treeVars.spill);
  tree -> SetBranchAddress("run",  &treeVars.run);
  
  tree -> SetBranchAddress("time",   treeVars.time);
  tree -> SetBranchAddress("amp_max",treeVars.amp_max);
  // tree -> SetBranchAddress("hodo.nFibresOnX",treeVars.nFibresOnX);
  // tree -> SetBranchAddress("hodo.nFibresOnY",treeVars.nFibresOnY);
  // tree -> SetBranchAddress("hodo.X",treeVars.hodoX);
  // tree -> SetBranchAddress("hodo.Y",treeVars.hodoY);
  // tree -> SetBranchAddress("wire.X",wireX);
  // tree -> SetBranchAddress("wire.Y",wireY);
  tree -> SetBranchAddress("n_tracks",  &treeVars.nTracks);
  tree -> SetBranchAddress("fitResult", &treeVars.tracks);
  
  if( opts.OptExist("Input.tableX0") ) tree -> SetBranchAddress("tableX",&treeVars.tableX);
  if( opts.OptExist("Input.tableY0") ) tree -> SetBranchAddress("tableY",&treeVars.tableY);
  
  std::vector<std::string> channels = opts.GetOpt<std::vector<std::string> >("Channels.channels");  
  for(auto ch : channels)
  {
    std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
    if( ampCh != "NULL" ) tree -> SetBranchAddress(ampCh.c_str(), &treeVars.channelIds[ampCh]);
    
    std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
    if( timeCh != "NULL" ) tree -> SetBranchAddress(timeCh.c_str(), &treeVars.channelIds[timeCh]);
    
    std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));    
    for(auto timeMethod: timeMethods)
      if( timeMethod != "NULL" ) tree -> SetBranchAddress(timeMethod.c_str(), &treeVars.timeMethodIds[timeMethod]);
    
    std::string Vbias = opts.OptExist(Form("%s.Vbias", ch.c_str())) ? opts.GetOpt<std::string>(Form("%s.Vbias", ch.c_str())) : "NULL";
    if( Vbias != "NULL" ) tree -> SetBranchAddress(Vbias.c_str(), &treeVars.VbiasVals[Vbias]);
    else treeVars.VbiasVals[Vbias] = -1.;
    
    std::string NINOthr = opts.OptExist(Form("%s.NINOthr", ch.c_str())) ? opts.GetOpt<std::string>(Form("%s.NINOthr", ch.c_str())) : "NULL";
    if( NINOthr != "NULL" ) tree -> SetBranchAddress(NINOthr.c_str(), &treeVars.NINOthrVals[NINOthr]);
    else treeVars.NINOthrVals[NINOthr] = -1.;
  }
}



void InitTreeVarsFNAL(TTree* tree, TreeVars& treeVars, CfgManager& opts)
{
  treeVars.time    = new float[468];
  treeVars.amp_max = new float[36];
  treeVars.nFibresOnX = new int[1];
  treeVars.nFibresOnY = treeVars.nFibresOnX;
  treeVars.hodoX = new float[4];
  treeVars.hodoY = new float[4];

  tree -> SetBranchStatus("*",0);
  
  tree -> SetBranchStatus("i_evt",1);   tree -> SetBranchAddress("i_evt",&treeVars.event);
  
  tree -> SetBranchStatus("gaus_mean",1);  tree -> SetBranchAddress("gaus_mean",&treeVars.time[0]);
  tree -> SetBranchStatus("IL_10mV",1);    tree -> SetBranchAddress("IL_10mV",  &treeVars.time[36]);
  tree -> SetBranchStatus("IL_20mV",1);    tree -> SetBranchAddress("IL_20mV",  &treeVars.time[72]);
  tree -> SetBranchStatus("IL_30mV",1);    tree -> SetBranchAddress("IL_30mV",  &treeVars.time[108]);
  tree -> SetBranchStatus("IL_50mV",1);    tree -> SetBranchAddress("IL_50mV",  &treeVars.time[144]);
  tree -> SetBranchStatus("IL_70mV",1);    tree -> SetBranchAddress("IL_70mV",  &treeVars.time[180]);
  tree -> SetBranchStatus("IL_90mV",1);    tree -> SetBranchAddress("IL_90mV",  &treeVars.time[216]);
  tree -> SetBranchStatus("IL_100mV",1);   tree -> SetBranchAddress("IL_100mV", &treeVars.time[252]);
  tree -> SetBranchStatus("IL_120mV",1);   tree -> SetBranchAddress("IL_120mV", &treeVars.time[288]);
  tree -> SetBranchStatus("IL_140mV",1);   tree -> SetBranchAddress("IL_140mV", &treeVars.time[324]);
  tree -> SetBranchStatus("IL_160mV",1);   tree -> SetBranchAddress("IL_160mV", &treeVars.time[360]);
  tree -> SetBranchStatus("IL_180mV",1);   tree -> SetBranchAddress("IL_180mV", &treeVars.time[396]);
  tree -> SetBranchStatus("IL_200mV",1);   tree -> SetBranchAddress("IL_200mV", &treeVars.time[432]);
  tree -> SetBranchStatus("amp",1);        tree -> SetBranchAddress("amp",      treeVars.amp_max);
  tree -> SetBranchStatus("ntracks",1);    tree -> SetBranchAddress("ntracks",  &treeVars.nFibresOnX[0]);
  tree -> SetBranchStatus("x_dut",1);      tree -> SetBranchAddress("x_dut",    treeVars.hodoX);
  tree -> SetBranchStatus("y_dut",1);      tree -> SetBranchAddress("y_dut",    treeVars.hodoY);
  
  if( opts.OptExist("Input.tableX0") ) tree -> SetBranchAddress("tableX",&treeVars.tableX);
  if( opts.OptExist("Input.tableY0") ) tree -> SetBranchAddress("tableY",&treeVars.tableY);

  for(int jj = 0; jj < 36; ++jj)
  {
    treeVars.channelIds[Form("%d",jj)] = jj;
  }
  for(int jj = 0; jj < 468; ++jj)
  {
    treeVars.timeMethodIds[Form("%d",jj)] = jj;
  }
}



void ReconstructHodoPosition(TreeVars& tv, CfgManager& opts, const std::string& ch, const bool& tracking, const int& plane)
{
  float hodo1X_z = opts.GetOpt<float>("Input.hodo1X_z");
  float hodo1Y_z = opts.GetOpt<float>("Input.hodo1X_z");
  float hodo2X_z = opts.GetOpt<float>("Input.hodo2Y_z");
  float hodo2Y_z = opts.GetOpt<float>("Input.hodo2Y_z");
  
  if( tracking )
  {
    float mX = (tv.hodoX[1]-tv.hodoX[0]-3.50) / (hodo2X_z-hodo1X_z);
    float mY = (tv.hodoY[1]-tv.hodoY[0]-0.12) / (hodo2Y_z-hodo1Y_z);
    float z = opts.GetOpt<float>(Form("%s.zPos",ch.c_str()));
    tv.beamX = tv.hodoX[0] + mX * (z-hodo1X_z);
    tv.beamY = tv.hodoY[0] + mY * (z-hodo1Y_z);
  }
  else
  {
    // tv.beamX = tv.hodoX[plane];
    // tv.beamY = tv.hodoY[plane];
    if( tv.tracks->size() >= 1 )
    {
      tv.beamX = (tv.tracks->at(0)).x();
      tv.beamY = (tv.tracks->at(0)).y();
    }
  }
  
  if( tv.run == 12426 ) tv.tableX = 42.;
  if( opts.OptExist("Input.tableX0") ) tv.beamX = (tv.tableX-opts.GetOpt<float>("Input.tableX0")) - tv.beamX;
  if( opts.OptExist("Input.tableY0") ) tv.beamY = (tv.tableY-opts.GetOpt<float>("Input.tableY0")) + tv.beamY;
}



bool AcceptEvent(AnalysisVars& av,
                 TreeVars& tv, CfgManager& opts,
                 const std::string& ch,
                 const int& nFibresMin, const int& nFibresMax,
                 const float& CTRMin, const float& CTRMax)
{
  std::string Vbias = opts.GetOpt<std::string>(Form("%s.Vbias", ch.c_str()));
  if( Vbias != "NULL" )
  {
    float VbiasMin = opts.GetOpt<float>(Form("%s.VbiasMin", ch.c_str()));
    float VbiasMax = opts.GetOpt<float>(Form("%s.VbiasMax", ch.c_str()));
    
    if( tv.VbiasVals[Vbias] < VbiasMin || tv.VbiasVals[Vbias] > VbiasMax) return false;
  }
  
  std::string NINOthr = opts.GetOpt<std::string>(Form("%s.NINOthr", ch.c_str()));
  if( NINOthr != "NULL" )
  {
    float NINOthrMin = opts.GetOpt<float>(Form("%s.NINOthrMin", ch.c_str()));
    float NINOthrMax = opts.GetOpt<float>(Form("%s.NINOthrMax", ch.c_str()));
    
    if( tv.NINOthrVals[NINOthr] < NINOthrMin || tv.NINOthrVals[NINOthr] > NINOthrMax) return false;
  }
  
  
  if( tv.nTracks != 1 ) return false;
  
  // if( (tv.nFibresOnX[0] < nFibresMin || tv.nFibresOnX[0] > nFibresMax) ) return false;
  // if( (tv.nFibresOnX[1] < nFibresMin || tv.nFibresOnX[1] > nFibresMax) ) return false;
  // if( (tv.nFibresOnY[0] < nFibresMin || tv.nFibresOnY[0] > nFibresMax) ) return false;
  // if( (tv.nFibresOnY[1] < nFibresMin || tv.nFibresOnY[1] > nFibresMax) ) return false;
  
  float xtalXMin = opts.GetOpt<float>(Form("%s.xtalXMin",ch.c_str()));
  float xtalXMax = opts.GetOpt<float>(Form("%s.xtalXMax",ch.c_str()));
  float xtalXCen = opts.GetOpt<float>(Form("%s.xtalXCen",ch.c_str()));
  float xtalYMin = opts.GetOpt<float>(Form("%s.xtalYMin",ch.c_str())); 
  float xtalYMax = opts.GetOpt<float>(Form("%s.xtalYMax",ch.c_str()));
  float xtalYCen = opts.GetOpt<float>(Form("%s.xtalYCen",ch.c_str()));
  if( (tv.beamX < xtalXMin) || (tv.beamX >= xtalXMax) ) return false;
  if( (tv.beamY < xtalYMin) || (tv.beamY >= xtalYMax) ) return false;
  
  std::string ampCh  = opts.GetOpt<std::string>(Form("%s.ampCh", ch.c_str()));
  std::string timeCh = opts.GetOpt<std::string>(Form("%s.timeCh",ch.c_str()));
  std::string trgCh  = opts.GetOpt<std::string>(Form("%s.trgCh",ch.c_str()));
  std::vector<std::string> timeMethods = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",ch.c_str()));
  float ampMin = opts.GetOpt<float>(Form("%s.ampMin",ch.c_str())); 
  float ampMax = opts.GetOpt<float>(Form("%s.ampMax",ch.c_str()));
  float timeMin = opts.GetOpt<float>(Form("%s.timeMin",ch.c_str())); 
  float timeMax = opts.GetOpt<float>(Form("%s.timeMax",ch.c_str()));
  float totMin = opts.GetOpt<float>(Form("%s.totMin",ch.c_str())); 
  float totMax = opts.GetOpt<float>(Form("%s.totMax",ch.c_str()));
  float rtMin = opts.GetOpt<float>(Form("%s.rtMin",ch.c_str()));
  float rtMax = opts.GetOpt<float>(Form("%s.rtMax",ch.c_str()));
  
  if( ampCh == "NULL" ) return false;
  float amp = tv.amp_max[tv.channelIds[ampCh]] / 4096.;
  if( amp < ampMin || amp >= ampMax ) return false;
  if( isnan(amp) ) return false;
  
  if( timeCh == "NULL" ) return false;
  float timTrg = trgCh != "NULL" ? tv.time[tv.channelIds[trgCh]+tv.timeMethodIds["LED"]] : 0;
  float tim = tv.time[tv.channelIds[timeCh]+tv.timeMethodIds[timeMethods.at(0)]] - timTrg;
  float tim2 = timeMethods.at(1) != "NULL" ? tv.time[tv.channelIds[timeCh]+tv.timeMethodIds[timeMethods.at(1)]] : 0.;
  float tim50 = tv.time[tv.channelIds[ampCh]+tv.timeMethodIds[timeMethods.at(2)]];
  float tim1000 = tv.time[tv.channelIds[ampCh]+tv.timeMethodIds[timeMethods.at(3)]];
  if( isinf(tim) ) return false;
  if( isnan(tim) ) return false;
  if( tim < timeMin || tim >= timeMax ) return false;
  if( (tim2 > tim) && ((tim2-tim) < totMin) ) return false;
  if( (tim2 > tim) && ((tim2-tim) >= totMax) ) return false;
  if( ((tim1000-tim50) > 0.) && ((tim1000-tim50) < rtMin) ) return false;
  if( ((tim1000-tim50) > 0.) && ((tim1000-tim50) >= rtMax) ) return false;
  
  std::string refCh  = opts.GetOpt<std::string>(Form("%s.refCh", ch.c_str()));
  std::string ampChRef  = opts.GetOpt<std::string>(Form("%s.ampCh", refCh.c_str()));
  std::vector<std::string> timeMethodsRef = opts.GetOpt<std::vector<std::string> >(Form("%s.timeMethods",refCh.c_str()));
  std::string timeChRef = opts.GetOpt<std::string>(Form("%s.timeCh",refCh.c_str()));
  std::string trgChRef  = opts.GetOpt<std::string>(Form("%s.trgCh",refCh.c_str()));
  float ampMinRef = opts.GetOpt<float>(Form("%s.ampMin",refCh.c_str())); 
  float ampMaxRef = opts.GetOpt<float>(Form("%s.ampMax",refCh.c_str()));
  float timeMinRef = opts.GetOpt<float>(Form("%s.timeMin",refCh.c_str())); 
  float timeMaxRef = opts.GetOpt<float>(Form("%s.timeMax",refCh.c_str()));
  float totMinRef = opts.GetOpt<float>(Form("%s.totMin",refCh.c_str())); 
  float totMaxRef = opts.GetOpt<float>(Form("%s.totMax",refCh.c_str()));
  float rtMinRef = opts.GetOpt<float>(Form("%s.rtMin",refCh.c_str()));
  float rtMaxRef = opts.GetOpt<float>(Form("%s.rtMax",refCh.c_str()));
    
  float ampRef = tv.amp_max[tv.channelIds[ampChRef]] / 4096.;
  if( ampRef < ampMinRef || ampRef >= ampMaxRef ) return false;
  if( isnan(ampRef) ) return false;
  
  float timTrgRef = trgChRef != "NULL" ? tv.time[tv.channelIds[trgChRef]+tv.timeMethodIds["LED"]] : 0;
  float timRef = tv.time[tv.channelIds[timeChRef]+tv.timeMethodIds[timeMethodsRef.at(0)]] - timTrgRef;
  float tim2Ref = timeMethodsRef.at(1) != "NULL" ? tv.time[tv.channelIds[timeChRef]+tv.timeMethodIds[timeMethodsRef.at(1)]] : 0.;
  float tim50Ref = timeMethodsRef.at(2) != "NULL" ? tv.time[tv.channelIds[ampChRef]+tv.timeMethodIds[timeMethodsRef.at(2)]] : 0.;
  float tim1000Ref = timeMethodsRef.at(3) != "NULL" ? tv.time[tv.channelIds[ampChRef]+tv.timeMethodIds[timeMethodsRef.at(3)]] : 0.;
  if( isinf(timRef) ) return false;
  if( isnan(timRef) ) return false;
  if( timRef < timeMinRef || timRef >= timeMaxRef ) return false;
  if( (tim2Ref > timRef) && ((tim2Ref-timRef) < totMinRef) ) return false;
  if( (tim2Ref > timRef) && ((tim2Ref-timRef) >= totMaxRef) ) return false;
  if( ((tim1000Ref-tim50Ref) > 0.) && ((tim1000Ref-tim50Ref) < rtMinRef) ) return false;
  if( ((tim1000Ref-tim50Ref) > 0.) && ((tim1000Ref-tim50Ref) >= rtMaxRef) ) return false;
  
  if( (tim-timRef) < CTRMin ) return false;
  if( (tim-timRef) >= CTRMax ) return false;
  
  av.amp = amp;
  av.time = tim;
  av.time2 = tim2;
  av.rt = tim1000 - tim50;
  av.dist = sqrt( pow(tv.beamX-xtalXCen,2) + pow(tv.beamY-xtalYCen,2) );
  av.timeRef = timRef;
  
  return true;
}
