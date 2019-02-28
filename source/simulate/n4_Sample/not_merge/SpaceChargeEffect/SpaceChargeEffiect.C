#include <iostream>
#include<fstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"
#include "Random.hh"
#include <cmath>
#include <TGraph.h>


using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {

  randomEngine.Seed(123456);
  TApplication app("app", &argc, argv);
  gRandom->SetSeed(5643);
 
  double alpha0, eta, vx, vy, vdrift;
  ifstream input_stream1("town60.txt");
  ifstream input_stream2("vdrift60.txt");
  ifstream input_stream3("cluster.txt");

  input_stream1 >> alpha0 >> eta;                                  
  input_stream2 >> vx >> vy >> vdrift;   

 
  int ncluster; 
  input_stream3 >> ncluster; 
   
  double x0cluster[ncluster];
  double xcluster[ncluster] = {0.};
  double siz0cluster[ncluster];
  double sizcluster[ncluster] = {0};
  double tcluster[ncluster];

  for (int icluster = 0; icluster < ncluster; icluster++) {
    input_stream3 >> x0cluster[icluster] >> siz0cluster[icluster];
    tcluster[icluster] = x0cluster[icluster] / SpeedOfLight;
  }
 
  input_stream1.close();
  input_stream2.close();
  input_stream3.close();

  const double tStep = 20 * 0.001;
  const double dStep = tStep * vdrift;
  const int nSteps = 6000;
  int N_t[nSteps] = {0};
  
  int alpha[nSteps+1];   alpha[0] =alpha0;
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
    cout << "iStep is " << iStep << endl ; 
    t_aval[iStep] = iStep * tStep;
    cout << "t_aval[iStep] = " << t_aval[iStep] << endl;

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
          cout << "iStep = "<< iStep <<", border xcluster = " << xcluster[icluster] << endl; 
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
    cout << "N_t[iStep] is " << N_t[iStep]<< endl;
    alpha[iStep+1] = alpha0 * nesatu/(nesatu+N_t[iStep]);
  }
  

  double i_t[nSteps] = {0.};
  
  const double gap = 0.2, thickness_insulatingfilm = 0.03, permitivity_insulatingfilm = 3.3;
  const double thickness_bakelite = 0.18, permitivity_bakelite = 5; 

 // double Ew = 1/(gap+2.*thickness_bakelite/permitivity_bakelite+
 //       2.*thickness_insulatingfilm/permitivity_insulatingfilm);
  double Ew = 1/(gap + 2.* thickness_bakelite/permitivity_bakelite +
              2.*thickness_insulatingfilm/permitivity_insulatingfilm)*10.;  //*10?

  double ElementaryCharge = 1.6e-19;
  for(int i_result = 0; i_result < nSteps; i_result++){
    i_t[i_result] = -N_t[i_result] * ElementaryCharge * Ew * vdrift;  
  }
  //plot
  TGraph *gr1 = new TGraph (nSteps, t_aval, i_t);
  TCanvas *c1 = new TCanvas("c1","Graph Draw Current");
 
  c1->cd();
  gr1->Draw("AC");
 // gr1->SetTitle("signal_current");
  //gr1->GetXaxis()->SetTitle("time[ns]");
  //gr1->GetYaxis()->SetTitle("signalmiu_A]");
  
  app.Run(kTRUE);

}
