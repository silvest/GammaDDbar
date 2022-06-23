#include <iostream>
#include <cstdio>
#include <map>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMatrixF.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TDatime.h>
#include <TF1.h>
#include <string>
#include <vector>
#include <math.h>
#include <TString.h>
#include <TComplex.h>

using namespace std;

void rescalehistos(TString file){

  TString outfile = file +"_rescaled.root";
  TString infile = file + ".root";
  
  TFile *input = new TFile(infile);

  double scalex = 1000.;
  double scaley = 1000.;
  
  TH2D* phiM12_vs_x12_in = (TH2D*) input->Get("phiM12_vs_x12");
  int nbinsx = phiM12_vs_x12_in->GetXaxis()->GetNbins();
  int nbinsy = phiM12_vs_x12_in->GetYaxis()->GetNbins();
  double xlow = scalex*phiM12_vs_x12_in->GetXaxis()->GetXmin();
  double xhi = scalex*phiM12_vs_x12_in->GetXaxis()->GetXmax();
  double ylow = scaley*phiM12_vs_x12_in->GetYaxis()->GetXmin();
  double yhi = scaley*phiM12_vs_x12_in->GetYaxis()->GetXmax();
  
  TH2D * phiM12_vs_x12_r = new TH2D("phiM12_vs_x12_r","phiM12_vs_x12_r",nbinsx,xlow,xhi,nbinsy,ylow,yhi);
  for (int i=1;i<=nbinsx;i++)
    for (int j=1;j<=nbinsy;j++)
      phiM12_vs_x12_r->Fill(scalex * phiM12_vs_x12_in->GetXaxis()->GetBinCenter(i),scaley * phiM12_vs_x12_in->GetYaxis()->GetBinCenter(j),phiM12_vs_x12_in->GetBinContent(i,j));

TH2D* phiG12_vs_y12_in = (TH2D*) input->Get("phiG12_vs_y12");
nbinsx = phiG12_vs_y12_in->GetXaxis()->GetNbins();
nbinsy = phiG12_vs_y12_in->GetYaxis()->GetNbins();
xlow = scalex*phiG12_vs_y12_in->GetXaxis()->GetXmin();
xhi = scalex*phiG12_vs_y12_in->GetXaxis()->GetXmax();
ylow = scaley*phiG12_vs_y12_in->GetYaxis()->GetXmin();
yhi = scaley*phiG12_vs_y12_in->GetYaxis()->GetXmax();

  
TH2D * phiG12_vs_y12_r = new TH2D("phiG12_vs_y12_r","phiG12_vs_y12_r",nbinsx,xlow,xhi,nbinsy,ylow,yhi);
for (int i=1;i<=nbinsx;i++)
  for (int j=1;j<=nbinsy;j++)
     phiG12_vs_y12_r->Fill(scalex * phiG12_vs_y12_in->GetXaxis()->GetBinCenter(i),scaley * phiG12_vs_y12_in->GetYaxis()->GetBinCenter(j),phiG12_vs_y12_in->GetBinContent(i,j));

  TH2D* phi_vs_qopm1_in = (TH2D*) input->Get("phi_vs_qopm1");
  nbinsx = phi_vs_qopm1_in->GetXaxis()->GetNbins();
  nbinsy = phi_vs_qopm1_in->GetYaxis()->GetNbins();
  xlow = scalex*phi_vs_qopm1_in->GetXaxis()->GetXmin();
  xhi = scalex*phi_vs_qopm1_in->GetXaxis()->GetXmax();
  ylow = scaley*phi_vs_qopm1_in->GetYaxis()->GetXmin();
  yhi = scaley*phi_vs_qopm1_in->GetYaxis()->GetXmax();

  
  TH2D * phi_vs_qopm1_r = new TH2D("phi_vs_qopm1_r","phi_vs_qopm1_r",nbinsx,xlow,xhi,nbinsy,ylow,yhi);
  for (int i=1;i<=nbinsx;i++)
    for (int j=1;j<=nbinsy;j++)
      phi_vs_qopm1_r->Fill(scalex * phi_vs_qopm1_in->GetXaxis()->GetBinCenter(i),scaley * phi_vs_qopm1_in->GetYaxis()->GetBinCenter(j),phi_vs_qopm1_in->GetBinContent(i,j));


  // TH1D* qopm1_in = (TH1D*) input->Get("qopm1");
  // int nbins = qopm1_in->GetXaxis()->GetNbins();
  // TH1D * qopm1 = new TH1D("qopm1","qopm1",nbins,100.*qopm1_in->GetXaxis()->GetBinLowEdge(qopm1_in->GetXaxis()->GetFirst()),100.*qopm1_in->GetXaxis()->GetBinUpEdge(qopm1_in->GetXaxis()->GetLast()));
  // for (int i=1;i<=nbins;i++0)
  //   qopm1->Fill(100. * qopm1_in->GetXaxis()->GetBinCenter(i),qopm1_in->getBinContent(i));

  TFile *output = new TFile(outfile,"RECREATE");
  output->cd();
  
  phi_vs_qopm1_r->Write();
  phiG12_vs_y12_r->Write();
  phiM12_vs_x12_r->Write();
  input->Close();
  output->Close();
  
}
