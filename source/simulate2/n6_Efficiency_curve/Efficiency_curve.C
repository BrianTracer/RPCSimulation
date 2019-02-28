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

#include <math.h>


using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {


  randomEngine.Seed(123456);
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

   // Histograms
  TH1::StatOverflows(kTRUE);
  gRandom->SetSeed(123425); //233 //245 //189  //257
  int xmax =45;
 
  double alpha[15], eta[15],lor[15], vx[15], vy[15], vdrift[15];
  ifstream input_stream1("town.txt");         //10-150KeV/cm
  ifstream input_stream2("vdrift.txt");       //10-150KeV/cm
  for(int i = 0; i < 15; i++){
    input_stream1 >> alpha[i] >> eta[i] >> lor[i];                                  
    cout << "alpha = "<< alpha[i] <<", eta = "<< eta[i] << endl;
    input_stream2 >> vx[i] >> vy[i] >> vdrift[i];   
    cout << "vx = "<< vx[i] << ", vy = "<< vy[i] << ", vdrift = " << vdrift[i] << endl;
  }
  input_stream1.close();
  input_stream2.close();
  
  float gap = 0.1;
  float bakelite = 0.18;
  float permitivity_bakelite = 5;
  int num_vol = 1;
  double voltage[num_vol];
  double Ef[num_vol];
  double alpha_Ef[num_vol];
  double beta_Ef[num_vol];
  double vdrift_Ef[num_vol];
  for (int i = 0; i<num_vol; i++){
  //  voltage[i] = 10000 + i*100;
  //  Ef[i] = voltage[i] /gap / 1000; // KeV
  //  Ef[i] = voltage[i]/(gap + 2*bakelite / permitivity_bakelite) / 1000;
    Ef[i] = 40; ///////////////////////////////////////////////////////////////////////////////////////////////////////////
    int tenPlace = floor(Ef[i]) / 10;
    float unitPlace = Ef[i] - tenPlace * 10;
    cout << "tenPlace = " << tenPlace << endl;
    cout << "unitPlace = " << unitPlace << endl;
    alpha_Ef[i] = alpha[tenPlace -1 ] + (float)(unitPlace / 10) * (alpha[tenPlace] - alpha[tenPlace - 1]);
    beta_Ef[i] = eta[tenPlace -1 ] + (float)(unitPlace / 10)*(eta[tenPlace] - eta[tenPlace -1]);
    vdrift_Ef[i] = vdrift[tenPlace -1 ] + (float)(unitPlace / 10)*(vdrift[tenPlace] - vdrift[tenPlace -1]);
    cout << "Ef[" << i << "] =" << Ef[i] << ", alpha = " << alpha_Ef[i] << ", beta = " << beta_Ef[i] << ", vdrift = "<< vdrift_Ef[i]<< endl;

  }
  
  TFile* cluster = new TFile("cluster1mm.root");
  
  //get cluster num
  TH1F* den = (TH1F*)cluster-> Get("hNc");
  //get cluster size
  TH1F* size = (TH1F*)cluster->Get("hnc");
//------------------------------------------------------------------------------------------------------------------

for(int i_vol = 0; i_vol < num_vol; i_vol++){
 /*   
  alpha_Ef[i_vol] = 61.306;
  beta_Ef[i_vol] = 24.9491;
  vdrift_Ef[i_vol] = 0.00924346;
 */ 
  int nevents = 2000;
  double T[nevents] = {0.}, TOT[nevents] = {0.};
  int save[nevents] ={0};
 for (int ievent = 0; ievent < nevents ;ievent++){
  cout << "ievent = " << ievent << endl;
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
    const double gap = 0.1;
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
  const double dStep = tStep * vdrift_Ef[i_vol];
  const int nSteps = 6000;
  long long int N_t[nSteps] = {0};
  
  double alpha[nSteps+1] = {alpha_Ef[i_vol]}; //  alpha[0] =alpha0;
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
    kratio = beta_Ef[i_vol] / alpha[iStep];
    Astep = exp( (alpha[iStep] - beta_Ef[i_vol])*dStep );
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


        if(xcluster[icluster] >=0.1) {
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
 //   if(iStep < 2000){
 //   cout << "N_t[iStep] is " << N_t[iStep]<< endl;
 //   cout << "alpha =" << alpha[iStep] << endl;
 //   }
      alpha[iStep+1] = alpha_Ef[i_vol] * nesatu/(nesatu+N_t[iStep]);  // space charge effect
  //    alpha[iStep + 1] = alpha_Ef[i_vol];
  }

 
  
  const double gap = 0.1, thickness_insulatingfilm = 0.03, permitivity_insulatingfilm = 3.3;
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
    i_t[i_result] = -N_t[i_result] * ElementaryCharge * Ew * vdrift_Ef[i_vol]*1e9 *1e6;  //miu_A
 //voltage calculate
    v_t[i_result] = i_t[i_result]*tStep*1e-9*exp((i_result*tStep)/(R*C*1e-3));
    if(i_result > 0 ) {
      v_t[i_result] += v_t[i_result -1];}
  }
  for(int i_result = 0; i_result < nSteps; i_result++){
    v_t[i_result] *=exp((-i_result*tStep)/(R*C*1e-3))/(C*1e-12);               //miu_V
   // cout << "v_t[" << i_result <<  "] = " << v_t[i_result] << endl;
  }
/////////////////////////////////////////////////////////////////////////////////////////
  double threshold_voltage = -0.5;
  for(int i_result = 0; i_result < nSteps; i_result++){
    if(i_result >0 && v_t[i_result] < threshold_voltage && v_t[i_result -1] > threshold_voltage){
      T[ievent] = (i_result - 1) * tStep + (v_t[i_result -1 ] - threshold_voltage)/(v_t[i_result -1]-v_t[i_result]) * tStep;
      save[ievent] = 1;
      break;
    }
  }
//  hT->Fill(ievent);

  for(int i_result = nSteps; i_result > 0; i_result--){
    if(i_result < nSteps && v_t[i_result] < threshold_voltage && v_t[i_result +1] > threshold_voltage){
      TOT[ievent] = (i_result - 1) * tStep + (v_t[i_result -1 ] - threshold_voltage)/(v_t[i_result -1]-v_t[i_result]) * tStep;
      TOT[ievent] -= T[ievent];
      break;
    }
  } 
 // cout << "ievent = " << ievent << ", T = " << T[ievent] << ", TOT = " << TOT[ievent] << endl;   
 }
  cluster->Close();

  TH1F* ht = new TH1F("ht", "before TOT",100, 0., 25.);
  TProfile *hprof = new TProfile("hprof","Profile of T versus TOT",100,0,xmax,0,25);  // 0,70

  TCanvas *c0 = new TCanvas("c0","ht");
  
  double T1[nevents] = {0.};
  double TOT1[nevents] = {0.};
  int nevents_cut = 0;
  for ( Int_t i=0; i<nevents; i++) {
    if (save[i] == 1){
      T1[nevents_cut] = T[i];
      TOT1[nevents_cut] = TOT[i];
      nevents_cut++;
      ht->Fill(T[i]);
      hprof->Fill(TOT[i],T[i],1);
    }
  }
  c0->cd();
  ht-> Draw();
  ht->GetXaxis()->SetTitle("Time [ns]");

 
  TGraph *gr1 = new TGraph (nevents_cut, TOT1, T1);
  TCanvas *c1 = new TCanvas("c1","Graph T-TOT");
  c1->cd();
  gr1->Draw("A*");
  gr1->SetMarkerSize(0.2);
  gr1->SetTitle("T-TOT");
  gr1->GetXaxis()->SetTitle("TOT[ns]");
  gr1->GetYaxis()->SetTitle("time[ns]");
  

   // Create a canvas giving the coordinates and the size
  TCanvas *c2 = new TCanvas("c2", "Profile example");
  c2->cd();
  hprof->Draw();
  hprof->GetXaxis()->SetTitle("TOT[ns]");
  hprof->GetYaxis()->SetTitle("time[ns]");

  TF1 *g1 = new TF1("m1","pol4",0,xmax);
  hprof->Fit(g1);
  double Pol4[4] = {0.};
  g1->GetParameters(&Pol4[0]);
  for (int i = 0; i<=4;i++){
 //   cout << Pol4[i] << endl;
  }
   
  //T-TOT correct
  double T_correct[nevents_cut] = {0.};
 
  TH1F *hc = new TH1F("hc", "after TOT",100, -10., 10.);
  for ( Int_t i=0; i<nevents_cut; i++) {
    T_correct[i] = Pol4[0] + Pol4[1]*TOT1[i] + Pol4[2]*pow(TOT1[i],2) + Pol4[3]*pow(TOT1[i],3) +Pol4[4]*pow(TOT1[i],4); 
    hc->Fill(T_correct[i] - T1[i]);
  }
  TCanvas *c3 = new TCanvas("c3", "TOT correct");
  c3->cd();
  hc->Draw();
  hc->GetXaxis()->SetTitle("Time_correction - Time [ns]");
//  hc->Fit("gaus");

  TF1 *g2 = new TF1("m2","gaus",-10,10);
  hc->Fit(g2);
  double par[3] = {0.};
  g2->GetParameters(&par[0]);
  double Time_resolution;
  Time_resolution = par[2];
 // cout << nevents_cut << ", "<< nevents<< endl;
  double efficiency = (double)nevents_cut/nevents;
  cout<< "i_vol = " << i_vol << endl;
  cout<< "efficiency = " << efficiency << endl;;
  cout<< "Time_resolution = " << Time_resolution << endl;

}
//---------------------------------------------------------------------------------------------------------------------------------------
  app.Run(kTRUE);

}
