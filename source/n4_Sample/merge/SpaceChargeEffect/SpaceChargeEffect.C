#include <iostream>
#include<fstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TFile.h>
#include "Plotting.hh"
#include "TF1.h"

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"
#include "Random.hh"
#include <cmath>
#include <TProfile.h>


using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {


  randomEngine.Seed(123456);
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

   // Histograms
  TH1::StatOverflows(kTRUE);
  gRandom->SetSeed(257); //233 //245 //189  //257
 
  double alpha0, eta, vx, vy, vdrift;
  ifstream input_stream1("town50.txt");
  ifstream input_stream2("vdrift50.txt");

  input_stream1 >> alpha0 >> eta;                                  
  input_stream2 >> vx >> vy >> vdrift;   

  input_stream1.close();
  input_stream2.close();

  TFile* cluster = new TFile("cluster.root");
  
  //get cluster num
  TH1F* den = (TH1F*)cluster-> Get("hNc");
  //get cluster size
  TH1F* size = (TH1F*)cluster->Get("hnc");

/////////////////////////////////////////////////////////////////////////////////////////
  double ncluster1,ncluster2,ncluster_rnd;
  int ncluster0,ncluster;
  
  ncluster1 = den->GetRandom();
  ncluster0 = int(ncluster1);
  ncluster2 = ncluster1 - ncluster0;
  ncluster_rnd = gRandom->Rndm();
  if(ncluster_rnd<ncluster2){
    ncluster = ncluster0 + 1;
  }else{
    ncluster = ncluster0;
  }
   
 // cout<<"ncluster = "<<ncluster<<endl;

  //write xcluster， sizcluster 
  double x0cluster[ncluster];
  double xcluster[ncluster] = {0.};
  double siz0cluster[ncluster];
  double sizcluster[ncluster] = {0};
  double tcluster[ncluster];

 
  for (int i = 0; i < ncluster; i++ ){

    //get cluster position
    const double gap = 0.2;
    x0cluster[i] = (gRandom->Rndm())*gap;
  //  cout<< "x0cluster[" << i << "] = "<< x0cluster[i] <<endl;
    //get cluster time
    tcluster[i] = x0cluster[i] / SpeedOfLight;

    //get cluster size
    double siz1,siz2,siz_rnd;
    int siz0,cl_size;
    siz1 = size->GetRandom();
    siz0 = int(siz1);
    siz2 = siz1 - siz0;
    siz_rnd = gRandom -> Rndm();
    if(siz_rnd<siz2) {
      cl_size = siz0 + 1;
    }else{
      cl_size = siz0;
    } 
    siz0cluster[i] = cl_size;
  //  cout<< "siz0cluster[" << i << "] = "<< siz0cluster[i] << endl;
  //  cout << i << endl;
    
  }
 // cout << "ncluster = " << ncluster << endl;
 
  


/////////////////////////////////////////////////////////////////////////////////
   
  const double tStep = 20 * 0.001;
  const double dStep = tStep * vdrift;
  const int nSteps = 6000;
  long long int N_t[nSteps] = {0};
  
  double alpha[nSteps+1] = {alpha0}; //  alpha[0] =alpha0;
  double kratio;
  double Astep;
  double s_rnd;   
  
  double rnddis ;

  int N_amp;
  double cltmean,cltsigma,N_ava;
  double stepsigma ;
  double t_aval[nSteps] = {0.};
  int nesatu = 5e7;

  for (int iStep = 0; iStep < nSteps; iStep++){
 //   cout << "iStep is " << iStep << endl ; 
    t_aval[iStep] = iStep * tStep;
//    cout << "t_aval[iStep] = " << t_aval[iStep] << endl;

 //   calculate relate Townsend coefficients
    kratio = eta / alpha[iStep];
    Astep = exp( (alpha[iStep] - eta)*dStep );
    rnddis = kratio *( Astep - 1 )/( Astep - kratio );
    stepsigma = sqrt(Astep*(Astep-1)*(1+kratio)/(1-kratio));

    for(int icluster = 0; icluster < ncluster; icluster++){
      if ( tcluster[icluster] > t_aval[iStep] + tStep) {
        sizcluster[icluster] = 0;

      }
      else if(  tcluster[icluster] >= t_aval[iStep] && tcluster[icluster] <= t_aval[iStep] + tStep){
        sizcluster[icluster] = siz0cluster[icluster];
        xcluster[icluster] = x0cluster[icluster];


      }
      else if( tcluster[icluster] < t_aval[iStep]){ 


        if(xcluster[icluster] >=0.2) {
          sizcluster[icluster] = 0; 
 //         cout << "iStep = "<< iStep <<", border xcluster = " << xcluster[icluster] << endl; 
          continue;
        }
        int N_this = 0;
        if(sizcluster[icluster] <  500){

          for(int ielectron = 0; ielectron < sizcluster[icluster]; ielectron++){
            s_rnd = gRandom -> Rndm();
            if(s_rnd < rnddis){ 
              N_this += 0;
            }
            else{
              N_amp = int(log((Astep-kratio)*(1-s_rnd)/(1-kratio)/Astep)/
                          log(1-(1-kratio)/(Astep-kratio)));
              N_this += (N_amp + 1);
            }
          }
         }
        else{
          N_ava = sizcluster[icluster];
          cltmean = N_ava * Astep;
          cltsigma = sqrt(N_ava*1.)*stepsigma;
          N_amp = int(gRandom->Gaus(cltmean,cltsigma));
          N_this = N_amp;
        }
        sizcluster[icluster] = N_this;
        xcluster[icluster] += dStep;
      }
      N_t[iStep] += sizcluster[icluster]; 
    }
    if(iStep < 2000){
    cout << "N_t[iStep] is " << N_t[iStep]<< endl;
    cout << "alpha =" << alpha[iStep] << endl;}
  //  alpha[iStep+1] = alpha0 * nesatu/(nesatu+N_t[iStep]);
      alpha[iStep + 1] = alpha0;
  }

 
  
  const double gap = 0.2, thickness_insulatingfilm = 0.03, permitivity_insulatingfilm = 3.3;
  const double thickness_bakelite = 0.18, permitivity_bakelite = 5; 

  double Ew = 1/(gap + 2.* thickness_bakelite/permitivity_bakelite +
              2.*thickness_insulatingfilm/permitivity_insulatingfilm);

  double ElementaryCharge = 1.6e-19;
  double i_t[nSteps] = {0.};
  const double R = 1000.; //Unit [ohm]
  const double C = 10.;  //Unit[pF]
  double v_t[nSteps] = {0.};
  for(int i_result = 0; i_result < nSteps; i_result++){
 //current calculate
    i_t[i_result] = -N_t[i_result] * ElementaryCharge * Ew * vdrift*1e9 *1e3;  //m_A
 //voltage calculate
    v_t[i_result] = i_t[i_result]*tStep*1e-9*exp((i_result*tStep)/(R*C*1e-3));
    if(i_result > 0 ) {
      v_t[i_result] += v_t[i_result -1];}
  }
  for(int i_result = 0; i_result < nSteps; i_result++){
    v_t[i_result] *=exp((-i_result*tStep)/(R*C*1e-3))/(C*1e-12);                    //m_V
   // cout << "v_t[" << i_result <<  "] = " << v_t[i_result] << endl;
  }
    //plot
  TGraph *gr1 = new TGraph (nSteps, t_aval, i_t);
  TCanvas *c1 = new TCanvas("c1","Graph Draw Current");
  TGraph *gr2 = new TGraph (nSteps, t_aval, v_t);
  TCanvas *c2 = new TCanvas("c2","Graph Draw Voltage");
  
  c1->cd();
  gr1->Draw("AC");
  gr1->SetTitle("signal_current");
  gr1->GetXaxis()->SetTitle("time[ms]");
  gr1->GetYaxis()->SetTitle("signal_mA]");

  c2->cd();
  gr2->Draw("AC");
  gr2->SetTitle("signal_voltage");
  gr2->GetXaxis()->SetTitle("time[ms]");
  gr2->GetYaxis()->SetTitle("signal_mV]");
 
  app.Run(kTRUE);

}
