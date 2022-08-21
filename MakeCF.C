//  file: MakeCF
// Making correlation of photon arrival time using the position of photons this is mehtod 1                                                
//
//  Programmer:  Sahar Nikkhah nikkhah.1@osu.edu
//  Revision history:    23/04/2021
//====================================================
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
#include <iostream>
#include <fstream>

//============ function prototypes==================


void ReadData(std::ifstream* data, int Nsamples, float* wf1, float* wf2);//this function reads the data and put them into arrays
int ReturnPeakPositions(float* wf, int Nsamples, float dt, int* pos);  // note that number of peaks is the return value
void Correlation(TH1D* cf, int Nphoton1,int Nphoton2, float dt, int* pos1, int* pos2);// making correlation and auto correlation
void RemovePhaseSpace(TH1D* CF, float taumax, int Nsamples);// this removes a triangle that is made by correlating over -t to t. this is for the fact that smaller numbers will repeat more than larger ones, you can remove it to see the result. 
void Normalize(TH1D* cf); // normalazing 
TF1* FitCentralPeak(TH1D* cf, bool isAutoCorrelation=false);// fiting the the plot at center (Around tau=0)



//=====================Start of main function===========================
void MakeCF(int runnum=702134142,// name of the run
	    int evMax=300,         // number of events
	    TString dirName="/home/sahar/Documents/Physics6810/Finalproject/TimeCorrelation/Method1")//please input your own directory
{

  
  //================getting info from the file ttree and print them out======
  TFile* tf = new TFile(Form("%s/Run%d.root",dirName.Data(),runnum),"READ");// reading the file
  TTree* Settings;
  tf->GetObject("ScopeSettings",Settings);//opening the tree and put the info in settings
  float TimeBase,SampleRate;
  int TimeStamp,Ntrigs,Nsamples;
  TObjString* TrgSource = new TObjString();
  TObjString* Comment = new TObjString();
  Settings->GetBranch("TimeBase")->SetAddress(&TimeBase);
  Settings->GetBranch("TimeStamp")->SetAddress(&TimeStamp);// the size of each event
  Settings->GetBranch("SampleRate")->SetAddress(&SampleRate);// it is the  speed of data taking
  Settings->GetBranch("Ntriggers")->SetAddress(&Ntrigs);// Number of events
  Settings->GetBranch("Nsamples")->SetAddress(&Nsamples);// number of samples
  Settings->GetEvent();
  float dt = 1.0/SampleRate;// delta t of each sample to next one
  tf->Close();
  // ============================Defining variables and histograms==========
  
  float wf1[300000],wf2[300000];
  std::ifstream* data = new std::ifstream(Form("%s/Run%d.data",dirName.Data(),runnum), std::ios::binary);
  float taumax = dt*((float)Nsamples);
  cout << "Taumax is " << taumax << endl;
  TFile* CFfile = new TFile(Form("Run%dCF.root",TimeStamp),"RECREATE");// creating a CF file with the same name just has CF at the end
  int Nbins = 500;// I played with this and found this number works fine
  double TauLo = -3e-6-0.5*dt;// just to make this the same range as method 2
  double TauHi =  3e-6-0.5*dt;//same
  TH1D* CF = new TH1D("CF",Form("CF Run%d",TimeStamp),Nbins,TauLo,TauHi);// we are going to fill these 
  TH1D* AutoCF1 = new TH1D("AutoCF1",Form("AutoCF Ch1 Run%d",TimeStamp),Nbins,TauLo,TauHi);
  TH1D* AutoCF2 = new TH1D("AutoCF2",Form("AutoCF Ch2 Run%d",TimeStamp),Nbins,TauLo,TauHi);
  

  //============ looping over events  starts=============================  
  for (int iev=0; iev<evMax; iev++){
    
    cout << "Looking at event " << iev+1 << " out of " << evMax << endl;
    ReadData(data,Nsamples,wf1,wf2);

    int pos1[50000],pos2[50000];      // these are the positions of the "photons"
    int Nphotons1 = ReturnPeakPositions(wf1,Nsamples,dt,pos1);
    int Nphotons2 = ReturnPeakPositions(wf2,Nsamples,dt,pos2);
    cout << "Channel 1 saw " << Nphotons1 << " photons, and Channel 2 saw " << Nphotons2 << endl;//shows the number of photon in each event



    //=================correlation====================================
    for (int i=0; i<Nphotons1; i++){
      for (int j=0; j<Nphotons2; j++){
	float tau = dt*((float)(pos1[i]-pos2[j]));
	CF->Fill(tau);

      }
    }//end of correlation loop
   

    //====Auto correlation functions are all possible pairs in one phototube
 
    int Nphot;               // this is set then reset for detectors 1, 2
    TH1D *Auto;    // these are set then reset for detectors 1, 2
    int* pos;                // this is set then reset for detectors 1, 2
    Nphot = Nphotons1;
    Auto = AutoCF1;
    pos = pos1;
    for (int i=0; i<Nphot; i++){
      for (int j=0; j<Nphot; j++){
	if (j==i) continue;// skipping the same photon
	double tau = dt*((double)(pos[i]-pos[j]));
	Auto->Fill(tau);
      }
    }
    Nphot = Nphotons2;
    Auto = AutoCF2;
   
    pos = pos2;
    for (int i=0; i<Nphot; i++){
      for (int j=0; j<Nphot; j++){
	if (j==i) continue;
	double tau = dt*((double)(pos[i]-pos[j]));
	Auto->Fill(tau);
	
      }
    }//end of autocorrelation loop
  }//end of events loop
 

  //================removing phase space====================
  RemovePhaseSpace(CF,taumax,Nsamples);
  RemovePhaseSpace(AutoCF1,taumax,Nsamples);
  RemovePhaseSpace(AutoCF2,taumax,Nsamples);

  //===============Normalizing=============================
  Normalize(CF);
  Normalize(AutoCF1);
  Normalize(AutoCF2);
  
  //===============Fitting and Results=====================
  TF1* result;
  result = FitCentralPeak(CF,false);
  double PH12  = result->GetParameter(0);
  double dPH12 = result->GetParError(0);
  result = FitCentralPeak(AutoCF1,true);
  double PH11  = result->GetParameter(0);
  double dPH11 = result->GetParError(0);
  result = FitCentralPeak(AutoCF2,true);
  double PH22  = result->GetParameter(0);
  double dPH22 = result->GetParError(0);

  //======================Draw==========================
  TCanvas* tc = new TCanvas("CorrFctn",Form("CF for Run%d",TimeStamp),1000,1000);
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
  
  //====CfFile write and close===========================
  CFfile->Write();
  CFfile->Close();
}// end of main function


//============================Functions==================================

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
// ====================================================================

int ReturnPeakPositions(float* wf, int Nsamples, float dt, int* pos){
  float threshold    = 7;        // from the waveform
  float deadtime     =15e-9;    //this is the width explained it power point 
  int   nDeadBuckets = deadtime/dt;// number of samples that should be skipped (explained in the power point)
  float hiVal = -999.0;
  for (int i=0; i<Nsamples; i++) if (wf[i]>hiVal) hiVal=wf[i];
  float baseline = hiVal-1.0;//getiing rid of noise see waveform
  cout << "Baseline is " << baseline << endl;
  
  //  now find the peaks!
  int numPeaks =0;
  for (int i=0; i<Nsamples; i++){
    if (wf[i]<baseline-threshold){
      pos[numPeaks++]=i;
      i += nDeadBuckets;
    }
  }
  return numPeaks;
}

//==================================================================
void RemovePhaseSpace(TH1D* CF, float taumax, int Nsamples){
  // this just divides out the "triangle" phasespace
  TH1D* triangle = (TH1D*)CF->Clone();
  for (int ibin=1; ibin<triangle->GetNbinsX()+1; ibin++){
    float taubin = triangle->GetBinCenter(ibin);
    float triangleVal = Nsamples*(1.0-fabs(taubin)/taumax);
    triangle->SetBinContent(ibin,triangleVal);
    triangle->SetBinError(ibin,0.0);
  }
  CF->Divide(triangle);
  delete triangle;   // deleting it
}   

//=====================================================

void Normalize(TH1D* cf){
  // assumes that left-most and right-most quarters of the data are                      
  //  in the "normalization region"                                                      
  double ave=0.0;
  int nbinAve=0;
  for (int i=1; i<cf->GetXaxis()->GetNbins()/4; i++){
    ave += cf->GetBinContent(i);
    nbinAve++;
  }
  for (int i=3*cf->GetXaxis()->GetNbins()/4; i<=cf->GetXaxis()->GetNbins(); i++){
    ave += cf->GetBinContent(i);
    nbinAve++;
  }
  ave /= (double)nbinAve;
  cf->Scale(1.0/ave);
 
}
//=========================================================

TF1* FitCentralPeak(TH1D* cf, bool isAutoCorrelation){
  TF1* mg2 = new TF1("limG2","1.0+fabs([0])*exp(-x*x/(2.0*[1]*[1]))",-1.e-6,1.e-6);
  mg2->SetParameter(0,0.1);
  mg2->SetParameter(1,1.0e-6);
  if (isAutoCorrelation){
    // for autocorrelation, must only fit *either* \tau>0 *or* \tau<0, since it's trivially symmetrized
    // also, deadtime effect ("track merging") means exclude bin closest to zero.
    double bw = cf->GetXaxis()->GetBinWidth(1);
    cf->Fit(mg2,"","",bw*1.5,1.e-6);
  }
  else{
    cf->Fit(mg2,"R");
  }
  return mg2;
}
//========================End of the code=================


