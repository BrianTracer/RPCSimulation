#include <iostream>
#include<fstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TGraph.h>

#include "/home/joshua/garfield/Include/MediumMagboltz.hh"
#include "/home/joshua/garfield/Include/FundamentalConstants.hh"

using namespace Garfield;
using namespace std;

int TownsendCof() {

  Float_t e[15],alpha[15],beta[15],lor[15],eta[15];
  
  ifstream input_stream("town.txt");
  if(input_stream.is_open()){
    int i=0;
    while(1){
    e[i]=10.e3+i*10.e3;
    input_stream >> alpha[i] >> beta[i] >> lor[i];
    eta[i]=alpha[i]-beta[i];
    if(!input_stream.good()) break;                                   
    cout<<i<<" : "<<e[i]<<"    "<<alpha[i]<<"    "<<eta[i]<<endl;
    i++;   
     }
  }
  input_stream.close();
  
  Float_t vx[15],vy[15],vdrift[15];
  ifstream input_stream2("vdrift.txt");
  if(input_stream2.is_open()){
    int i=0;
    while(1){
    input_stream2 >> vx[i] >> vy[i] >> vdrift[i];
    if(!input_stream2.good()) break;                                   
    cout<<i<<" : "<<vx[i]<<"    "<<vy[i]<<"    "<<vdrift[i]<<endl;
    i++;   
     }
  }
  input_stream2.close();


  TGraph *gr1 = new TGraph (15, e, alpha);
  TGraph *gr2 = new TGraph (15, e, beta);
  TGraph *gr3 = new TGraph (15, &e[0], &eta[0]);
  TGraph *gr4 = new TGraph (15, e, vdrift);

  TCanvas *c1 = new TCanvas("c1","Graph Draw Alpha");
  TCanvas *c2 = new TCanvas("c2","Graph Draw Beta");
  TCanvas *c3 = new TCanvas("c3","Graph Draw Eta");
  TCanvas *c4 = new TCanvas("c4","Graph Draw Drift");

  c1->cd();
  gr1->Draw("AC*");
  gr1->SetTitle("gasfile_alpha");
  gr1->GetXaxis()->SetTitle("E[V/cm]");

  gr1->GetXaxis()->SetTitle("E[V/cm]");
  gr1->GetYaxis()->SetTitle("alpha[cm^-1]");
  c1->SaveAs("alpha.pdf");

  c2->cd();
  c2->SetLogy();
  gr2->Draw("AC*");
  gr2->SetTitle("gasfile_beta");
  gr2->GetXaxis()->SetTitle("E[V/cm]");
  gr2->GetYaxis()->SetTitle("beta[cm^-1]");
  c2->SaveAs("beta.pdf");

  c3->cd();
  gr3->Draw("AC*");
  gr3->SetTitle("gasfile_eta");
  gr3->GetXaxis()->SetTitle("E[V/cm]");
  gr3->GetYaxis()->SetTitle("eta[cm^-1]");
  c3->SaveAs("eta.pdf");

  c4->cd();
  gr4->Draw("AC*");
  gr4->SetTitle("gasfile_vdrift");
  gr4->GetXaxis()->SetTitle("E[V/cm]");
  gr4->GetYaxis()->SetTitle("vdrift[cm/ns]");
  c4->SaveAs("drift.pdf");

  return 0;
}
