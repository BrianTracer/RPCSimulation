#include <iostream>
#include <fstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TGraph.h>

using namespace std;

int plot_efficiency() {

  float gap = 0.1;
  Float_t Voltage[25],Ef[25],efficiency[25],time_resolution[25];
  
  ifstream input_stream("correct.txt");
  int num_vol = 0;
  while (!input_stream.eof()) {
    input_stream >> Ef[num_vol] >> efficiency[num_vol] >> time_resolution[num_vol];
    Voltage[num_vol] = Ef[num_vol]*(gap)*1000;
    cout << num_vol << ", " << Ef[num_vol] << ", " << efficiency[num_vol] << ", " << time_resolution[num_vol] << endl; 
    num_vol++;  
  }
  input_stream.close();

  TCanvas *c1 = new TCanvas("c1","Graph Draw efficiency");
  TPad *pad1 = new TPad("pad1","This is pad1",0.05,0.49,0.95,0.99);
  TPad *pad2 = new TPad("pad2","This is pad2",0.05,0.,0.95,0.49);

  pad1->Draw();
  pad2->Draw();
  
  
  TGraph *gr1 = new TGraph (num_vol-1, Voltage, efficiency);
  TGraph *gr2 = new TGraph (num_vol-1, Voltage, time_resolution);


 // TCanvas *c2 = new TCanvas("c2","Graph Draw time_resolution");

  pad1->cd();
  gr1->Draw("AC*");
  gr1->SetTitle("Efficiency(1mm gas gap)");
  //gr1->GetXaxis()->SetTitle("Voltage[V]");
  gr1->GetYaxis()->SetTitle("Efficiency");
  //c1->SaveAs("efficiency.pdf");

  pad2->cd();
  gr2->Draw("A*");
  gr2->SetTitle("Time_resolution(1mm gas gap)");
  gr2->GetXaxis()->SetTitle("Voltage[V]");
  gr2->GetYaxis()->SetTitle("Time_resolution[ns]");
  c1->SaveAs("Efficiency and Time_resolution.pdf");

  return 0;
}
