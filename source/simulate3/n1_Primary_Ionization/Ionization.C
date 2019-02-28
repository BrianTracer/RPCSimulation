#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TFile.h>

#include "MediumMagboltz.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "ComponentConstant.hh"
#include "Sensor.hh"
#include "TrackHeed.hh"
#include "Plotting.hh"
#include "Random.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  randomEngine.Seed(123456);
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  // Histograms
  TH1::StatOverflows(kTRUE);
 
  TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons",
                              120, 0, 120);
  TH1F* hEdep = new TH1F("hEdep", "Energy Loss",
                         100, 0., 6.);

  TH1F* hxc = new TH1F("hxc", "track", 100, -0.1, 1.6);
  TH1F* hyc = new TH1F("hyc", "Y-axis", 300, -10., 10.);
  TH1F* hzc = new TH1F("hzc", "Z-axis", 300, -10., 10.);
  TH1F* htc = new TH1F("htc", "t-axis", 100, -0.001, 0.008);
  TH1F* hnc = new TH1F("hnc", "Cluster Size", 200,0., 500.);
  TH1F* hNc = new TH1F("hNc", "Number of Clusters", 50, 0., 50.);
  TH1F* hxe = new TH1F("hxe", "track of electrons", 200, -0.1, 1.6);
  TH1F* hye = new TH1F("hye", "Y-axis of electrons", 300, -10., 10.);
  TH1F* hze = new TH1F("hze", "Z-axis of electrons", 300, -10., 10.);
  // Make a medium
  const double temperature = ZeroCelsius + 20;
  const double pressure = 1.0 * AtmosphericPressure;  

  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition("c2f4h2", 94.7, "ic4h10",5.,"sf6",0.3);
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);

  // Detector geometry
  // Gap [cm]	include gas and box
  const double width = 0.15;  //Unit[cm]                 //11111111111111111111111111111111111111111111111111111111
  SolidBox* box = new SolidBox(width / 2., 0., 0., width / 2., 10., 10.);
  GeometrySimple* geo = new GeometrySimple();
  geo->AddSolid(box, gas);

  // Make a component    include geo and electricfield     
  ComponentConstant* comp = new ComponentConstant();
  comp->SetGeometry(geo); 	
  comp->SetElectricField(0., 0., 0.);

  // Make a sensor     
  Sensor* sensor = new Sensor();
  sensor->AddComponent(comp);

  // Track class    include Sensor, Particle and Momentum
  TrackHeed* track = new TrackHeed();
  track->SetSensor(sensor);
  track->SetParticle("muon");
  track->SetMomentum(120.e9);

  const int nEvents = 100000;
  track->EnableDebugging();
  for (int i = 0; i < nEvents; ++i) {
    if (i == 1) track->DisableDebugging();
    if (i % 1000 == 0) std::cout << i << "/" << nEvents << "\n";
    // Initial position and direction 
    double x0 = 0., y0 = 0., z0 = 0., t0 = 0.;
    double dx0 = 1., dy0 = 0., dz0 = 0.; 
    track->NewTrack(x0, y0, z0, t0, dx0, dy0, dz0);
    // Cluster coordinates
    double xc = 0., yc = 0., zc = 0., tc = 0.;
    // Number of electrons produced in a collision
    int nc = 0;
    // Energy loss in a collision
    double ec = 0.;
    // Dummy variable (not used at present)
    double extra = 0.;
    // Total energy loss along the track
    double esum = 0.;
    // Total number of electrons produced along the track
    int nsum = 0;
    // The number of cycles
    int Nc = 0; 
    // Loop over the clusters.
    
    double xe,ye,ze,te,ee,dxe,dye,dze;
    while (track->GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
      esum += ec;
      nsum += nc;
      Nc++;
      hxc->Fill(10*xc);
      hyc->Fill(yc);
      hzc->Fill(zc);
      htc->Fill(tc);
      hnc->Fill(nc);

      for (int ie =0; ie<nc;ie++ ){
        track ->GetElectron(ie,xe,ye,ze,te,ee,dxe,dye,dze);
        hxe->Fill(10*xe);
        hye->Fill(ye);
        hze->Fill(ze);
      }
    }
    hNc->Fill(Nc);
    hElectrons->Fill(nsum);
    hEdep->Fill(esum * 1.e-3);
  }

  TFile *fcluster = new TFile("cluster1s5mm.root","RECREATE");
  //空间分布
  TCanvas* c1 = new TCanvas();
  hxc->GetXaxis()->SetTitle("Cluster Num along the track[mm]");
  hxc->Draw();
  hxc->Write();
  c1->SaveAs("Cluster Num along the track.pdf");

  TCanvas* c2 = new TCanvas();
  hyc->GetXaxis()->SetTitle("Cluster Num on the Y-axis[cm]");
  hyc->Draw();
  hyc->Write();
  c2->SaveAs("Cluster Num on the Y-axis.pdf");

  TCanvas* c3 = new TCanvas();
  hzc->GetXaxis()->SetTitle("Cluster Num on the Z-axis[cm]");
  hzc->Draw();
  hzc->Write();
  c3->SaveAs("Cluster Num on the Z-axis.pdf");
  //团簇时间分布
  TCanvas* c4 = new TCanvas();
  htc->GetXaxis()->SetTitle("Cluter Num on the t-axis[ns]");
  htc->Draw();
  htc->Write();
  c4->SaveAs("Cluter Num on the t-axis.pdf");
  //能损分布
  TCanvas* c5 = new TCanvas();
  hEdep->GetXaxis()->SetTitle("energy loss [keV]");
  hEdep->Draw();
  hEdep->Write();;
  c5->SaveAs("edep.pdf");
  //团簇的尺寸分布
  TCanvas* c6 = new TCanvas();
  hnc->GetXaxis()->SetTitle("Cluster Size");
  c6->SetLogy();
  hnc->Draw();
  hnc->Write();
  c6->SaveAs("Cluster Size.pdf");
  //团簇的数目
  TCanvas* c7 = new TCanvas();
  hNc->GetXaxis()->SetTitle("number of clusters");
  hNc->Draw();
  hNc->Write();
  c7->SaveAs("number of clusters.pdf");

  fcluster->Close();

  TFile *felectron = new TFile("electron.root","RECREATE");
  //电子数目分布（每个事例在探测器中电离产生的总电子数目的分布情况）
  TCanvas* c8 = new TCanvas();
  hElectrons->GetXaxis()->SetTitle("number of electrons"); 
  hElectrons->Draw();
  hElectrons->Write();
  c8->SaveAs("ne.pdf");
  //电子空间分布
  TCanvas* c9 = new TCanvas();
  hxe->GetXaxis()->SetTitle("Electron Num along the track[mm]");
  hxe->Draw();
  hxe->Write();
  c9->SaveAs("Electron Num along the track.pdf");

  TCanvas* c10 = new TCanvas();
  hye->GetXaxis()->SetTitle("Electron Num on the Y-axis[cm]");
  hye->Draw();
  hye->Write();
  c10->SaveAs("Electron Num on the Y-axis.pdf");

  TCanvas* c11 = new TCanvas();
  hze->GetXaxis()->SetTitle("Electron Num on the Z-axis[cm]");
  hze->Draw();
  hze->Write();
  c11->SaveAs("Electron Num on the Z-axis.pdf");

  felectron->Close();

  app.Run(kTRUE); 

}
