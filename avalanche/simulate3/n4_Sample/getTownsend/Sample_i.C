#include <iostream>
#include<fstream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

   TApplication app("app", &argc, argv);
 
  const double pressure = 1 * AtmosphericPressure;
  const double temperature = 293.15;
 
  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("C2F4H2", 94.7, "iC4H10", 5., "SF6", 0.3);
 
  double e=40.7e3, b=0, btheta=90;
  bool verbose=true;
  double vx, vy, vz, dl, dt, alpha, eta, lor, vxerr, vyerr, vzerr, dlerr, dterr, alphaerr, etaerr, lorerr,
          alphatof;
  const int ncoll = 10;
 
  gas-> RunMagboltz(e,b,btheta,
              ncoll, verbose,
              vx, vy, vz, dl,  dt,
              alpha, eta, lor,
              vxerr, vyerr, vzerr,
              dlerr, dterr,
              alphaerr, etaerr, lorerr,
              alphatof);
 
  std::cout<<"vx = "<< vx <<", vy = "<< vy <<", vz = "<< vz <<"/n";
  std::cout<<"alpha = "<< alpha << ", eta = " << eta << std::endl;

  std::ofstream output_stream;
  std::ofstream output_stream2;

  output_stream.open("town.txt"); 
  output_stream2.open("vdrift.txt");

  output_stream << alpha <<" "<<eta << std::endl;
  output_stream2 << vx <<" "<< vy <<" "<< vz <<" "<<std::endl;
  
  output_stream.close();
  output_stream2.close();


  app.Run(kTRUE);

}
