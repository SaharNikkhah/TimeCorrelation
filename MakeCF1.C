
//  file: MakeCF
// Making correlation of photon arrival time going over all the samples  method 2                                            
//
//  Programmer:  Sahar Nikkhah nikkhah.1@osu.edu
//  Revision history:    23/04/2021 
//
//
  
#include "TFile.h"
#include "TNtuple.h"
#include <TROOT.h>
#include <TStyle.h>
#include "stdio.h"
#include "TText.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include <fstream>
#include "TObjString.h"
#include <iostream>
#include <string>
#include <ostream>
#include "TF1.h"



//========== function prototypes===================
void ReadData(std::ifstream* data, int Nsamples, float* wf1, float* wf2);// reading data and put them in a w1 and w2
void Correlation(float dt, float* wf1, float* wf2, int Nsamples,  TH1F* CF);//creating correlation
void AutoCorrelation(float dt, float* wf1,  int Nsamples,  TH1F* AutoCF);//creating autocorrelation
void Normalize(TH1F* CF);// normalizing data
//void RemovePhaseSpace(TH1D* CF, float taumax, int Nsamples); we do not need to remove  phase space because we are doing slide window for doingg correlation
TF1* FitCentralPeak(TH1F* cf, bool isAutoCorrelation=false);// fitting the result

void MakeCF1(int runnum=702134142,
	    int maxEvent=-1,         
	   TString dirName="/home/sahar/Documents/Physics6810/Finalproject/TimeCorrelation/method2")
  
{




  // =============================getting info from files=========================
  TFile* tf = new TFile(Form("%s/Run%d.root",dirName.Data(),runnum),"READ");
  TTree* Settings;
  tf->GetObject("ScopeSettings",Settings);
  float TimeBase,SampleRate;
  int TimeStamp,Ntrigs,Nsamples;
  TObjString* TrgSource = new TObjString();
  TObjString* Comment = new TObjString();
  Settings->GetBranch("TimeBase")->SetAddress(&TimeBase);
  Settings->GetBranch("TimeStamp")->SetAddress(&TimeStamp);// size of each event
  Settings->GetBranch("SampleRate")->SetAddress(&SampleRate);// speed of data taking
  Settings->GetBranch("Ntriggers")->SetAddress(&Ntrigs);//number of events in each run
  Settings->GetBranch("Nsamples")->SetAddress(&Nsamples);// number of sample in each event
  Settings->GetBranch("TrgSource")->SetAddress(&TrgSource);
  Settings->GetBranch("UserComment")->SetAddress(&Comment);
  Settings->GetEvent();
  float dt = 1.0/SampleRate;

  tf->Close();
  // ====================================================================
  

  
  float wf1[500000],wf2[500000];//these will contain the waveforms data in them
  std::ifstream* data = new std::ifstream(Form("%s/Run%d.data",dirName.Data(),runnum), std::ios::binary);
  float taumax = dt*((float)Nsamples);// dt is the time difference in each sample in our case is around 10^-9


  
  //============================== Creating Histograms====================
  TFile* CFfile = new TFile(Form("Run%dCF.root",TimeStamp),"RECREATE"); //creating a new file to save the cf file in it
  TH1F* CF = new TH1F("CF",Form("CF Run%d",TimeStamp),1001,-3e-6+0.5*dt, 3e-6-0.5*dt);// I am using this time range because of my sliding window size in correlation function
  TH1F* AutoCF1 = new TH1F("AutoCF1",Form("AutoCF Ch1 Run%d",TimeStamp),1001,-3e-6+0.5*dt, 3e-6-0.5*dt);
  TH1F* AutoCF2 = new TH1F("AutoCF2",Form("AutoCF Ch2 Run%d",TimeStamp),1001,-3e-6+0.5*dt, 3e-6-0.5*dt);




  
  
  //==========================Correlation loop=============================
  int evMax = (maxEvent>0)?maxEvent:Ntrigs;// looping over all the events
  cout<<evMax;
  for (int iev=0; iev<evMax; iev++){

  

    ReadData(data,Nsamples,wf1,wf2);// reading data
    Correlation(dt, wf1,  wf2, Nsamples,CF);
    AutoCorrelation(dt, wf1,  Nsamples,  AutoCF1);
    AutoCorrelation(dt, wf2,  Nsamples,  AutoCF2);
    cout<<Nsamples;
  
    if (iev%10==0){//to see the progress in the code
      cout<<"event number is; " <<iev<<endl;
    } 
 
  }//end of events loop
  
  // ==========="normalization region"=========================                                                     
  Normalize(CF);
  Normalize(AutoCF1);
  Normalize(AutoCF2);

  
  //==================---Fitting===============================
  TF1* result;
  result = FitCentralPeak(CF,false);
  result = FitCentralPeak(AutoCF1,true);
  result = FitCentralPeak(AutoCF2,true);


  //===================== now drawing=========================

  TCanvas* tc = new TCanvas("CorrFctn",Form("CF for Run%d",TimeStamp));
  tc->Divide(3,1);
  tc->Draw();
  tc->cd(1);
  gPad->SetGridx();
  CF->DrawCopy();
  tc->cd(2);
  gPad->SetGridx();
  AutoCF1->DrawCopy();
  tc->cd(3);
  gPad->SetGridx();
  AutoCF2->DrawCopy();

  
  //cf file write and close
  CFfile->Write();
  CFfile->Close();
}//=======================end of main function========================






//===========================Functions===============================

void ReadData(std::ifstream* data, int Nsamples, float* wf1, float* wf2){
  char buff1[500000];
  char buff2[500000];
  data->read(buff1,Nsamples);
  data->read(buff2,Nsamples);
  for (int ipt=0; ipt<Nsamples; ipt++){
    buff1[ipt]+=128;               
    buff2[ipt]+=128;
    wf1[ipt] = (float)buff1[ipt];
    wf2[ipt] = (float)buff2[ipt];   
  }
}

//====================================================================

void Correlation(float dt, float* wf1, float* wf2, int Nsamples,TH1F* CF){
  
  float hiVal1 = -999.0;// this is just to make sure hivalue is less than our signals(signals are negative)
  for (int i=0; i<Nsamples; i++) if (wf1[i]>hiVal1) hiVal1=wf1[i];// to determine the real hivalue
  float baseline1 = hiVal1-6;//removing the noise see the waveform
  float hiVal2 = -999.0;
  for (int i=0; i<Nsamples; i++) if (wf2[i]>hiVal2) hiVal2=wf2[i];
  float baseline2 = hiVal2-6;//-1 one is noise, check the waveforms
  cout<<"===="<<baseline1<< " and" << baseline2<<endl;
  // looping over samples in sliding window 
  for (int i=3000 ; i<Nsamples-3000; i++){
    if ((baseline1-wf1[i])<0) wf1[i]=baseline1;
    for (int j=i-3000; j<i+3000; j++){
      if ((baseline2-wf2[j])<0) wf2[j]=baseline2;
      float weight=(baseline1-wf1[i])*(baseline2-wf2[j]);// I think this is needed, the size of bin matters not where it is. see waveforms
      CF-> Fill(dt*(float)(i-j),weight );
    }
  }
}
//==================================================================

void AutoCorrelation(float dt, float* wf1,  int Nsamples,TH1F* AutoCF){
      
  float hiVal1 = -999.0;
  for (int i=0; i<Nsamples; i++) if (wf1[i]>hiVal1) hiVal1=wf1[i];
  float baseline1 = hiVal1-6;
  
  for (int i=3000 ; i<Nsamples-3000; i++){
    if ((baseline1-wf1[i])<0) wf1[i]=baseline1;
    for (int j=i-3000; j<i+3000; j++){
      if (j==i) continue;// the only difference between the correlation and auto is here
      if ((baseline1-wf1[j])<0) wf1[j]=baseline1;
      float weight=(baseline1-wf1[i])*(baseline1-wf1[j]);// I think this is needed, the size of bin matters not where it is. see waveforms
      AutoCF-> Fill(dt*(float)(i-j),weight );
   }
  }
}
//============================================================

void Normalize(TH1F* CF){
  double ave=0.0;
  int nbinAve=0;
  for (int i=1; i<CF->GetXaxis()->GetNbins()/4; i++){
    ave += CF->GetBinContent(i);
    nbinAve++;
  }
  for (int i=3*CF->GetXaxis()->GetNbins()/4; i<=CF->GetXaxis()->GetNbins(); i++){
    ave += CF->GetBinContent(i);
    nbinAve++;
  }
  ave /= (double)nbinAve;
  CF->Scale(1.0/ave);
}
//===========================================================

TF1* FitCentralPeak(TH1F* cf, bool isAutoCorrelation){
  TF1* mg2 = new TF1("limG2","1.0+fabs([0])*exp(-x*x/(2.0*[1]*[1]))",-1.e-6,1.e-6);//guassian
  mg2->SetParameter(0,0.1);
  mg2->SetParameter(1,1.0e-6);
  if (isAutoCorrelation){
    // just one side of it for auto  since it's  symmetrized
    // also, deadtime effect  means exclude bin closest to zero.
    double bw = cf->GetXaxis()->GetBinWidth(1);
    cf->Fit(mg2,"","",bw*1.5,1.e-6);
  }
  else{
    cf->Fit(mg2,"R");
  }
  return mg2;
}
//======================End of code=======================
