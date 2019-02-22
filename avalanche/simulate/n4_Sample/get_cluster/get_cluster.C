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
#include "fstream"
#include "RandomEngine.hh"

using namespace Garfield;
using namespace std;

int main(int argc, char * argv[]) {

 // randomEngine.Seed(222226);
  gRandom->SetSeed(237886);
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  ofstream output_stream;
  output_stream.open("cluster_new.txt"); 

  TFile *cluster = new TFile("cluster.root");
  
  //get cluster num
  TH1F *den = (TH1F*)cluster-> Get("hNc");
  //get cluster size
  TH1F *size = (TH1F*)cluster->Get("hnc");

  double ncluster1,ncluster2,ncluster_rnd;
  int ncluster0,ncluster;
  
  const int nEvents = 1;
  int Ncluster = 0;
  for (int iEvent = 0; iEvent < nEvents; iEvent++){
    ncluster1 = den->GetRandom();
    ncluster0 = int(ncluster1);
    ncluster2 = ncluster1 - ncluster0;
    ncluster_rnd = gRandom->Rndm();
    if(ncluster_rnd<ncluster2){
      ncluster = ncluster0 + 1;
    }else{
      ncluster = ncluster0;
    }
    Ncluster += ncluster ;
  }
 
  output_stream << Ncluster << endl;
  cout<<"Ncluster = "<<Ncluster<<endl;


  //write xcluster， sizcluster 
  double xcluster[Ncluster] = {0};
  double sizcluster[Ncluster] = {0};


 
  for (int i = 0; i < Ncluster; i++ ){

    //get cluster position
    const double gap = 0.2;
    xcluster[i] = (gRandom->Rndm())*gap;
    cout<< "xcluster[" << i << "] = "<< xcluster[i] <<endl;

 
  
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
    sizcluster[i] = cl_size;
    cout<< "sizcluster[" << i << "] = "<< sizcluster[i] << endl;
    
    output_stream << xcluster[i] <<" "<< sizcluster[i] << endl;
    cout << i << endl;
    
  }
  cout << "Ncluster = " << Ncluster << endl;
 
  cluster->Close();
  output_stream.close();
  

  
  app.Run(kTRUE); 

}
