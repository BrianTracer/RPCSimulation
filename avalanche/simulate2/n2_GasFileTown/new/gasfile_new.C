#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
 
  const double temperature = ZeroCelsius + 20;
  const double pressure = 1.0 * AtmosphericPressure;  
 
  // Setup the gas.
  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetTemperature(temperature);
  gas->SetPressure(pressure);
  gas->SetComposition("C2F4H2", 94.7, "iC4H10", 5., "SF6", 0.3);
 
  // Set the field range to be covered by the gas table. 
  const int nFields = 15;
  const double emin = 10.e3;
  const double emax = 150.e3;
  // Flag to request logarithmic spacing.
  const bool useLog = false;
  gas->SetFieldGrid(emin, emax, nFields, useLog); 

  const int ncoll = 10;
  // Switch on debugging to print the Magboltz output.
  gas->EnableDebugging();
  // Run Magboltz to generate the gas table.
  gas->GenerateGasTable(ncoll);
  gas->DisableDebugging();
  // Save the table. 
  gas->SetComposition("C2F4H2", 94.7, "iC4H10", 5., "SF6", 0.3);

  app.Run(kTRUE);

}
