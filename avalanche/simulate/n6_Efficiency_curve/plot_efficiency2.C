#include <iostream>
#include <fstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TGraphErrors.h>

using namespace std;

int plot_efficiency2() {

  float gap = 0.2;
  Float_t Voltage[20],Ef[20],efficiency[20],time_resolution[20],time_error[20],voltage_error[20] = {0};
  
  ifstream input_stream("correct2mm.txt");
  int num_vol = 0;
  while (!input_stream.eof()) {
    input_stream >> Ef[num_vol] >> efficiency[num_vol] >> time_resolution[num_vol] >> time_error[num_vol];
    //Voltage[num_vol] = Ef[num_vol]*(gap)*1000;
    Voltage[num_vol] = Ef[num_vol]*(gap+2*0.18/10)*1000;
    cout << num_vol << ", " << Ef[num_vol] << ", " << efficiency[num_vol] << ", " << time_resolution[num_vol] << ", " << time_error[num_vol] << endl; 
    num_vol++;  
  }
  input_stream.close();

  TCanvas *c1 = new TCanvas("c1","Graph Draw efficiency");

  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.49,0.95,0.99);
  TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.,0.95,0.49);

  pad1->SetGrid();
  pad2->SetGrid();
  pad1->Draw();
  pad2->Draw();
  
  
  TGraph *gr1 = new TGraph (num_vol-1, Voltage, efficiency);
//  TGraph *gr2 = new TGraph (num_vol-1, Voltage, time_resolution);
  TGraphErrors *gr2 = new TGraphErrors(num_vol-1,Voltage,time_resolution,voltage_error,time_error);


 // TCanvas *c2 = new TCanvas("c2","Graph Draw time_resolution");

  pad1->cd();
  gr1->Draw("ACP");
  gr1->SetMarkerColor(4);
  gr1->SetMarkerStyle(21);
  gr1->SetMarkerSize(0.8);
  gr1->SetTitle("Efficiency(2mm gas gap)");
  //gr1->GetXaxis()->SetTitle("Voltage[V]");
  gr1->GetYaxis()->SetTitle("Efficiency");
  //gr1->SaveAs("efficiency.pdf");

  pad2->cd();
  gr2->Draw("AP");
  gr2->SetMarkerColor(4);
  gr2->SetMarkerStyle(21);
  gr2->SetMarkerSize(0.7);
  gr2->SetTitle("Time_resolution(2mm gas gap)");
  gr2->GetXaxis()->SetTitle("Voltage[V]");
  gr2->GetYaxis()->SetTitle("Time_resolution[ns]");

  c1->Update();
  c1->SaveAs("2mm Efficiency and Time_resolution2.pdf");
  return 0;
}
