#include <iostream>
#include<fstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TGraph.h>

#include "MediumMagboltz.hh"
#include "MediumGas.hh"
#include "FundamentalConstants.hh"
#include "ComponentConstant.hh"
#include "Sensor.hh"
#include "AvalancheMC.hh"
#include "AvalancheMicroscopic.hh"
#include "ViewSignal.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "TrackHeed.hh"
#include "Plotting.hh"
#include "Random.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {

 // randomEngine.Seed(123456);
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  const double temperature = ZeroCelsius + 20;
  const double pressure = 1.0 * AtmosphericPressure;  

  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition("c2f4h2", 97, "ic4h10",3);
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->LoadGasFile("c2f4h2_97_ic4h10_3.gas");
 
   const double width = 0.2;  //Unit[cm]
  SolidBox* box = new SolidBox(width / 2., 0., 0., width / 2., 10., 10.);
  GeometrySimple* geo = new GeometrySimple();
  geo->AddSolid(box, gas);



  const double gap = 0.2, thickness_insulatingfilm = 0.03, permitivity_insulatingfilm = 3.3;
  const double thickness_bakelite = 0.18, permitivity_bakelite = 5; 

  double Ew = 1/(gap+2.*thickness_bakelite/permitivity_bakelite+
        2.*thickness_insulatingfilm/permitivity_insulatingfilm);

   // Make a component    include geo and electricfield     
  ComponentConstant* comp = new ComponentConstant();
  comp->SetGeometry(geo); 	
  comp->SetElectricField(40.7e3,0.,0.);
  comp->SetWeightingField(Ew, 0., 0., "s");
  
  // Make a sensor     
  Sensor* sensor = new Sensor();
  sensor->AddComponent(comp);
  sensor->AddElectrode(comp,"s");

  const double tStart = -5;
  const double tStop  = 50;
  const unsigned int nSteps = 1050;
  const double tStep = int((tStop - tStart)/nSteps);
 // sensor->SetTimeWindow(tStart, tStep, nSteps);

  AvalancheMC* aval = new AvalancheMC();
  aval->SetSensor(sensor);
//  aval->SetTimeSteps(0.005);
//  aval->SetCollisionSteps(100);
  aval->SetDistanceSteps(0.002);
  aval->EnableSignalCalculation();
  
  // Track class    include Sensor, Particle and Momentum
  TrackHeed* track = new TrackHeed();
  track->SetSensor(sensor);
  track->SetParticle("muon");
  track->SetMomentum(120.e9);
  
  sensor->SetTimeWindow(-5,0.2, 600);
  const int nEvents = 1;
  track->EnableDebugging();
  for (int i = 0; i < nEvents; ++i) {
    if (i == 1) track->DisableDebugging();
    if (i % 1 == 0) std::cout << i << "/" << nEvents << "\n";
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
   
    double xe,ye,ze,te,ee,dxe,dye,dze;
    while (track->GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
      for (int ie =0; ie<nc;ie++ ){
        track ->GetElectron(ie,xe,ye,ze,te,ee,dxe,dye,dze);
        aval->AvalancheElectron(xe, ye, ze, te);
      }
    }

  }
 
  //double signalbin = 10000;
  //sensor->GetSignal("s",signalbin); 
  ViewSignal* signalView = new ViewSignal();
  signalView->SetSensor(sensor);
  signalView->PlotSignal("s");
  
  app.Run(kTRUE); 

  return 0;}
