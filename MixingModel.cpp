#include "MixingModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCGaussianPrior.h>

#include <TRandom3.h>
using namespace std;

// ---------------------------------------------------------

MixingModel::MixingModel(int combination, bool phig12_i, bool noagamma_i, bool nokspipi_i, bool nok3pi_i, bool nokp_i, bool rad_i, bool epsK_i, bool ckmcorr_i, bool combgamma_all_i, bool combgamma_delta_i) : BCModel(), histos(obs)
{
  TH1::SetDefaultBufferSize(1000000);
  phig12 = phig12_i; // turn on the phi_Gamma measurement
  noagamma = noagamma_i; // turn off the Agamma measurement
  nokspipi = nokspipi_i; // unused in this code
  nok3pi = nok3pi_i; // unused in this code
  nokp = nokp_i; // unused in this code
  rad = rad_i; // If rad is true print the root histograms in rad, otherwise they are defined in degrees
  epsK = epsK_i; // If true takes into account the correction due to epsK
  ckmcorr = ckmcorr_i; // If true insert the Rcp parameter
  combgamma_all = combgamma_all_i; // add the LHCb comb results
  combgamma_delta = combgamma_delta_i; // add the LHCb comb results in D-> K pi sector
  comb = combination; // set the combination
  r2d = (rad ? 1. : 180. / M_PI);
  d2r = M_PI / 180.;

  //------------------------------------------ Inserting the measurements --------------------------------------------------------------------------

  if(comb == 0){
    Add_time_integrated_meas(); //B decay chain time integrated
    Add_time_dependent_Bmeas(); // time dependent B decay and mixing
    Add_time_dependent_Dmeas(); // time dependent D decay and mixing
    Add_other_meas(); // Other useful inputs
  } // LHCb Combination
  else if(comb == 1){
    Add_old_meas();  //mesures of the old code
  } // Old Combination
  else if(comb == 2){
    Add_time_integrated_meas(); //B decay chain time integrated
    Add_time_dependent_Bmeas(); // time dependent B decay and mixing
    Add_time_dependent_Dmeas(); // time dependent D decay and mixing
    Add_other_meas(); // Other useful inputs
    Add_old_meas(); // old code measurements
  } // Global combination LHCb + Old

  //------------------------------------------- Defining the Histograms --------------------------------------------------------------

  DefineHistograms();

  //------------------------------------------- Defining the parameters ---------------------------------------------------------------------------

  DefineParameters();

};

// ---------------------------------------------------------
MixingModel::~MixingModel()
{  // default destructor
};
// ---------------------------------------------------------

// ---------------------------------------------------------
void MixingModel::Add_time_integrated_meas(){

  vector<dato> CorrData;
  TMatrixDSym Corr, Corr2;

  //------------------------------------------------- Time Integrated B Decay Chain-------------------------------------------------------------------------

  //1. PDF: glwads-dh-hh-dmix (UID0)
  //Observables 8:
  CorrData.clear();
  CorrData.push_back(dato(0.136,0.009,0.001)); //acp_dk
  CorrData.push_back(dato(-0.008,0.002,0.002)); //acp_dpi
  CorrData.push_back(dato(-0.011,0.003,0.002)); //afav_dk
  CorrData.push_back(dato(0.95,0.009,0.01)); //rcp
  CorrData.push_back(dato(0.0095,0.0005,0.0003)); // rm_dk
  CorrData.push_back(dato(0.00415,0.00008,0.00004)); // rm_dpi
  CorrData.push_back(dato(0.02520,0.0008,0.0004)); // rp_dk
  CorrData.push_back(dato(0.00320,0.00007,0.00004)); // rp_dpi
  //Correlation Matrix (stat):
  Corr.ResizeTo(8,8);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.;
  Corr(0,2) = 0.02;
  Corr(0,3) = -0.02;
  Corr(0,4) = 0.;
  Corr(0,5) = 0.;
  Corr(0,6) = 0.;
  Corr(0,7) = 0.;
  Corr(1,1) = 1.;
  Corr(1,2) = 0.08;
  Corr(1,3) = 0.;
  Corr(1,4) = 0.01;
  Corr(1,5) = 0.01;
  Corr(1,6) = 0.;
  Corr(1,7) = -0.01;
  Corr(2,2) = 1.;
  Corr(2,3) = 0.;
  Corr(2,4) = 0.;
  Corr(2,5) = 0.01;
  Corr(2,6) = 0.;
  Corr(2,7) = -0.01;
  Corr(3,3) = 1.;
  Corr(3,4) = 0.03;
  Corr(3,5) = 0.;
  Corr(3,6) = 0.04;
  Corr(3,7) = 0.;
  Corr(4,4) = 1.;
  Corr(4,5) = -0.03;
  Corr(4,6) = 0.05;
  Corr(4,7) = 0.02;
  Corr(5,5) = 1.;
  Corr(5,6) = 0.02;
  Corr(5,7) = 0.08;
  Corr(6,6) = 1.;
  Corr(6,7) = -0.04;
  Corr(7,7) = 1.;
  //Correlation Matrix (syst)
  Corr2.ResizeTo(8,8);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.13;
  Corr2(0,2) = 0.07;
  Corr2(0,3) = -0.22;
  Corr2(0,4) = -0.03;
  Corr2(0,5) = 0.;
  Corr2(0,6) = 0.02;
  Corr2(0,7) = -0.07;
  Corr2(1,1) = 1.;
  Corr2(1,2) = -0.74;
  Corr2(1,3) = 0.18;
  Corr2(1,4) = 0.04;
  Corr2(1,5) = 0.31;
  Corr2(1,6) = -0.07;
  Corr2(1,7) = -0.21;
  Corr2(2,2) = 1.;
  Corr2(2,3) = -0.02;
  Corr2(2,4) = 0.05;
  Corr2(2,5) = -0.24;
  Corr2(2,6) = 0.11;
  Corr2(2,7) = 0.22;
  Corr2(3,3) = 1.;
  Corr2(3,4) = 0.09;
  Corr2(3,5) = 0.11;
  Corr2(3,6) = 0.17;
  Corr2(3,7) = 0.14;
  Corr2(4,4) = 1.;
  Corr2(4,5) = 0.16;
  Corr2(4,6) = 0.93;
  Corr2(4,7) = 0.13;
  Corr2(5,5) = 1.;
  Corr2(5,6) = 0.21;
  Corr2(5,7) = 0.84;
  Corr2(6,6) = 1.;
  Corr2(6,7) = 0.26;
  Corr2(7,7) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID0", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //2. PDF: glwads-dh-h3pi-dmix (UID1)
  //Observables 8:
  CorrData.clear();
  CorrData.push_back(dato(-0.313,0.102,0.038)); //aads_dk_k3pi
  CorrData.push_back(dato(0.023,0.048,0.005)); //aads_dpi_k3pi
  CorrData.push_back(dato(0.1,0.034,0.018)); //acp_dk_4pi
  CorrData.push_back(dato(-0.0041,0.0079,0.0024)); //acp_dpi_4pi
  CorrData.push_back(dato(0.00,0.012,0.002)); // afav_dk_k3pi
  CorrData.push_back(dato(0.014,0.0015,0.0006)); // rads_dk_k3pi
  CorrData.push_back(dato(0.00377,0.00018,0.00006)); // rads_dpi_k3pi
  CorrData.push_back(dato(0.975,0.037,0.019)); // rcp_4pi
  //Correlation Matrix (stat):
  Corr.ResizeTo(8,8);
  Corr(0,0) = 1.;
  Corr(0,1) = -0.06;
  Corr(0,2) = 0.0;
  Corr(0,3) = 0.01;
  Corr(0,4) = 0.01;
  Corr(0,5) = 0.08;
  Corr(0,6) = 0.;
  Corr(0,7) = 0.;
  Corr(1,1) = 1.;
  Corr(1,2) = 0.01;
  Corr(1,3) = 0.02;
  Corr(1,4) = 0.02;
  Corr(1,5) = 0.01;
  Corr(1,6) = -0.02;
  Corr(1,7) = 0.;
  Corr(2,2) = 1.;
  Corr(2,3) = -0.02;
  Corr(2,4) = 0.02;
  Corr(2,5) = 0.;
  Corr(2,6) = 0.;
  Corr(2,7) = -0.02;
  Corr(3,3) = 1.;
  Corr(3,4) = 0.1;
  Corr(3,5) = 0.;
  Corr(3,6) = 0.;
  Corr(3,7) = 0.;
  Corr(4,4) = 1.;
  Corr(4,5) = 0.;
  Corr(4,6) = 0.;
  Corr(4,7) = 0.;
  Corr(5,5) = 1.;
  Corr(5,6) = -0.05;
  Corr(5,7) = 0.04;
  Corr(6,6) = 1.;
  Corr(6,7) = 0.;
  Corr(7,7) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(8,8);
  Corr2(0,0) = 1.;
  Corr2(0,1) = -0.09;
  Corr2(0,2) = -0.04;
  Corr2(0,3) = 0.02;
  Corr2(0,4) = 0.02;
  Corr2(0,5) = 0.87;
  Corr2(0,6) = 0.05;
  Corr2(0,7) = -0.04;
  Corr2(1,1) = 1.;
  Corr2(1,2) = -0.34;
  Corr2(1,3) = 0.43;
  Corr2(1,4) = 0.05;
  Corr2(1,5) = 0.10;
  Corr2(1,6) = 0.46;
  Corr2(1,7) = -0.04;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.31;
  Corr2(2,4) = 0.09;
  Corr2(2,5) = 0.03;
  Corr2(2,6) = -0.35;
  Corr2(2,7) = 0.07;
  Corr2(3,3) = 1.;
  Corr2(3,4) = 0.32;
  Corr2(3,5) = 0.01;
  Corr2(3,6) = 0.24;
  Corr2(3,7) = -0.07;
  Corr2(4,4) = 1.;
  Corr2(4,5) = 0.02;
  Corr2(4,6) = -0.02;
  Corr2(4,7) = 0.02;
  Corr2(5,5) = 1.;
  Corr2(5,6) = 0.14;
  Corr2(5,7) = 0.04;
  Corr2(6,6) = 1.;
  Corr2(6,7) = -0.06;
  Corr2(7,7) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID1", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //3. PDF: glwads-dh-hhpipi0-dmix (UID2)
  //Observables 11:
  CorrData.clear();
  CorrData.push_back(dato(-0.2,0.27,0.04)); //aads_dk_kpipipi0
  CorrData.push_back(dato(0.438,0.19,0.011)); //aads_dpi_kpipipi0
  CorrData.push_back(dato(0.3,0.2,0.02)); //acp_dk_kkpi0
  CorrData.push_back(dato(0.054,0.091,0.011)); //acp_dk_pipipi0
  CorrData.push_back(dato(-0.03,0.04,0.005)); // acp_dpi_kkpi0
  CorrData.push_back(dato(-0.016,0.02,0.004)); // acp_dpi_pipipi0
  CorrData.push_back(dato(0.01,0.026,0.005)); // afav_dk_kpipi0
  CorrData.push_back(dato(0.014,0.0047,0.0021)); // rads_dk_kpipi0
  CorrData.push_back(dato(0.00235,0.00049,0.00006)); //rads_dpi_kpipi0
  CorrData.push_back(dato(0.95,0.22,0.05)); //rcp_kkpi0
  CorrData.push_back(dato(0.98,0.11,0.05)); //rcp_pipipi0
  //Correlation Matrix (stat):
  Corr.ResizeTo(11,11);
  Corr(0,0) = 1.;
  Corr(0,1) = -0.04;
  Corr(0,2) = 0.0;
  Corr(0,3) = 0.0;
  Corr(0,4) = 0.0;
  Corr(0,5) = 0.01;
  Corr(0,6) = 0.01;
  Corr(0,7) = 0.13;
  Corr(0,8) = 0.;
  Corr(0,9) = 0.;
  Corr(0,10) = 0.;
  Corr(1,1) = 1.;
  Corr(1,2) = 0.0;
  Corr(1,3) = 0.0;
  Corr(1,4) = 0.0;
  Corr(1,5) = 0.01;
  Corr(1,6) = 0.01;
  Corr(1,7) = -0.01;
  Corr(1,8) = -0.34;
  Corr(1,9) = 0.;
  Corr(1,10) = 0.;
  Corr(2,2) = 1.;
  Corr(2,3) = 0.;
  Corr(2,4) = -0.04;
  Corr(2,5) = 0.01;
  Corr(2,6) = 0.01;
  Corr(2,7) = 0.0;
  Corr(2,8) = 0.;
  Corr(2,9) = -0.2;
  Corr(2,10) = -0.01;
  Corr(3,3) = 1.;
  Corr(3,4) = 0.01;
  Corr(3,5) = -0.04;
  Corr(3,6) = 0.02;
  Corr(3,7) = 0.;
  Corr(3,8) = 0.;
  Corr(3,9) = 0.;
  Corr(3,10) = -0.04;
  Corr(4,4) = 1.;
  Corr(4,5) = 0.04;
  Corr(4,6) = 0.04;
  Corr(4,7) = 0.;
  Corr(4,8) = 0.;
  Corr(4,9) = 0.;
  Corr(4,10) = 0.;
  Corr(5,5) = 1.;
  Corr(5,6) = 0.08;
  Corr(5,7) = 0.;
  Corr(5,8) = 0.;
  Corr(5,9) = 0.;
  Corr(5,10) = 0.;
  Corr(6,6) = 1.;
  Corr(6,7) = 0.;
  Corr(6,8) = 0.;
  Corr(6,9) = 0.;
  Corr(6,10) = 0.;
  Corr(7,7) = 1.;
  Corr(7,8) = 0.03;
  Corr(7,9) = 0.;
  Corr(7,10) = 0.01;
  Corr(8,8) = 1.;
  Corr(8,9) = 0.;
  Corr(8,10) = 0.;
  Corr(9,9) = 1.;
  Corr(9,10) = 0.02;
  Corr(10,10) = 1.;


  //Correlation Matrix (syst)
  Corr2.ResizeTo(11,11);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.03;
  Corr2(0,2) = 0.07;
  Corr2(0,3) = 0.07;
  Corr2(0,4) = 0.18;
  Corr2(0,5) = 0.17;
  Corr2(0,6) = -0.16;
  Corr2(0,7) = 0.81;
  Corr2(0,8) = 0.32;
  Corr2(0,9) = 0.02;
  Corr2(0,10) = 0.13;
  Corr2(1,1) = 1.;
  Corr2(1,2) = 0.28;
  Corr2(1,3) = 0.31;
  Corr2(1,4) = 0.67;
  Corr2(1,5) = 0.68;
  Corr2(1,6) = -0.63;
  Corr2(1,7) = -0.18;
  Corr2(1,8) = -0.49;
  Corr2(1,9) = 0.;
  Corr2(1,10) = -0.04;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.77;
  Corr2(2,4) = 0.07;
  Corr2(2,5) = 0.05;
  Corr2(2,6) = 0.05;
  Corr2(2,7) = 0.08;
  Corr2(2,8) = -0.08;
  Corr2(2,9) = -0.33;
  Corr2(2,10) = -0.18;
  Corr2(3,3) = 1.;
  Corr2(3,4) = 0.05;
  Corr2(3,5) = 0.02;
  Corr2(3,6) = -0.06;
  Corr2(3,7) = 0.13;
  Corr2(3,8) = -0.11;
  Corr2(3,9) = -0.14;
  Corr2(3,10) = -0.25;
  Corr2(4,4) = 1.;
  Corr2(4,5) = 0.88;
  Corr2(4,6) = -0.82;
  Corr2(4,7) = -0.04;
  Corr2(4,8) = 0.02;
  Corr2(4,9) = -0.04;
  Corr2(4,10) = 0.02;
  Corr2(5,5) = 1.;
  Corr2(5,6) = -0.87;
  Corr2(5,7) = -0.03;
  Corr2(5,8) = 0.;
  Corr2(5,9) = 0.;
  Corr2(5,10) = 0.01;
  Corr2(6,6) = 1.;
  Corr2(6,7) = -0.05;
  Corr2(6,8) = 0.06;
  Corr2(6,9) = 0.04;
  Corr2(6,10) = 0.;
  Corr2(7,7) = 1.;
  Corr2(7,8) = 0.33;
  Corr2(7,9) = -0.03;
  Corr2(7,10) = -0.02;
  Corr2(8,8) = 1.;
  Corr2(8,9) = 0.02;
  Corr2(8,10) = -0.02;
  Corr2(9,9) = 1.;
  Corr2(9,10) = 0.38;
  Corr2(10,10) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID2", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  //4. PDF: ggsz-dh (UID3)
  //Observables 6:
  CorrData.clear();
  CorrData.push_back(dato(0.0568,0.0096,0.0031)); //xm_dk
  CorrData.push_back(dato(0.06550,0.01140,0.004300)); //ym_dk
  CorrData.push_back(dato(-0.09300,0.0098,0.003)); //xp_dk
  CorrData.push_back(dato(-0.0125,0.0123,0.0038)); //yp_dk
  CorrData.push_back(dato(-0.05470,0.0199,0.0035)); // xi_x_dpi
  CorrData.push_back(dato(0.00710,0.02330,0.0057)); // xi_y_dpi
  //Correlation Matrix (stat):
  Corr.ResizeTo(6,6);
  Corr(0,0) = 1.;
  Corr(0,1) = -0.12;
  Corr(0,2) = -0.01;
  Corr(0,3) = 0.02;
  Corr(0,4) = 0.04;
  Corr(0,5) = -0.16;
  Corr(1,1) = 1.;
  Corr(1,2) = -0.01;
  Corr(1,3) = -0.01;
  Corr(1,4) = 0.1;
  Corr(1,5) = 0.04;
  Corr(2,2) = 1.;
  Corr(2,3) = 0.1;
  Corr(2,4) = -0.11;
  Corr(2,5) = 0.03;
  Corr(3,3) = 1.;
  Corr(3,4) = -0.07;
  Corr(3,5) = -0.15;
  Corr(4,4) = 1.;
  Corr(4,5) = 0.15;
  Corr(5,5) = 1.;
  //Correlation Matrix (syst)
  Corr2.ResizeTo(6,6);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.30;
  Corr2(0,2) = 0.16;
  Corr2(0,3) = 0.58;
  Corr2(0,4) = 0.27;
  Corr2(0,5) = 0.23;
  Corr2(1,1) = 1.;
  Corr2(1,2) = 0.44;
  Corr2(1,3) = 0.22;
  Corr2(1,4) = 0.18;
  Corr2(1,5) = 0.17;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.45;
  Corr2(2,4) = 0.41;
  Corr2(2,5) = 0.31;
  Corr2(3,3) = 1.;
  Corr2(3,4) = 0.35;
  Corr2(3,5) = 0.24;
  Corr2(4,4) = 1.;
  Corr2(4,5) = 0.50;
  Corr2(5,5) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID3", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //5. PDF: glwads-dkdpi-kskpi-dmix (UID4)
  //Observables 7:
  CorrData.clear();
  CorrData.push_back(dato(-0.02,0.011,0.003)); //afav_dpi_kskpi
  CorrData.push_back(dato(0.007,0.017,0.003)); //asup_dpi_kskpi
  CorrData.push_back(dato(0.084,0.049,0.008)); //afav_dk_kskpi
  CorrData.push_back(dato(0.021,0.094,0.017)); //asup_dk_kskpi
  CorrData.push_back(dato(2.585,0.057,0.019)); //rfavsup_dpi_kskpi
  CorrData.push_back(dato(0.079,0.004,0.002)); // rfav_dkdpi_kskpi
  CorrData.push_back(dato(0.0620,0.006,0.003)); // rsup_dkdpi_kskpi

  //Correlation Matrix (stat):
  Corr.ResizeTo(7,7);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.;
  Corr(0,2) = -0.05;
  Corr(0,3) = 0.;
  Corr(0,4) = 0.;
  Corr(0,5) = -0.01;
  Corr(0,6) = 0.;
  Corr(1,1) = 1.;
  Corr(1,2) = 0.;
  Corr(1,3) = -0.05;
  Corr(1,4) = 0.;
  Corr(1,5) = 0.;
  Corr(1,6) = 0.;
  Corr(2,2) = 1.;
  Corr(2,3) = 0.;
  Corr(2,4) = 0.;
  Corr(2,5) = 0.;
  Corr(2,6) = 0.;
  Corr(3,3) = 1.;
  Corr(3,4) = 0.;
  Corr(3,5) = 0.;
  Corr(3,6) = -0.02;
  Corr(4,4) = 1.;
  Corr(4,5) = -0.11;
  Corr(4,6) = 0.15;
  Corr(5,5) = 1.;
  Corr(5,6) = 0.06;
  Corr(6,6) = 1.;
  //Correlation Matrix (syst)
  Corr2.ResizeTo(7,7);
  Corr2(0,0) = 1.;
  Corr2(0,1) = -0.88;
  Corr2(0,2) = 0.74;
  Corr2(0,3) = 0.05;
  Corr2(0,4) = 0.;
  Corr2(0,5) = 0.;
  Corr2(0,6) = 0.;
  Corr2(1,1) = 1.;
  Corr2(1,2) = -0.73;
  Corr2(1,3) = -0.08;
  Corr2(1,4) = 0.;
  Corr2(1,5) = 0.;
  Corr2(1,6) = 0.;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.63;
  Corr2(2,4) = 0.;
  Corr2(2,5) = -0.17;
  Corr2(2,6) = 0.;
  Corr2(3,3) = 1.;
  Corr2(3,4) = 0.;
  Corr2(3,5) = 0.;
  Corr2(3,6) = -0.1;
  Corr2(4,4) = 1.;
  Corr2(4,5) = 0.;
  Corr2(4,6) = 0.;
  Corr2(5,5) = 1.;
  Corr2(5,6) = 0.25;
  Corr2(6,6) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID4", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  //6. PDF: glwads-dsth-hh-dmix (UID5)
  //Observables 16:
  CorrData.clear();
  CorrData.push_back(dato(0.123,0.054,0.031)); //acp_dstk_dg
  CorrData.push_back(dato(-0.115,0.019,0.009)); //acp_dstk_dp
  CorrData.push_back(dato(-0.004,0.014,0.003)); //afav_dstk_dg
  CorrData.push_back(dato(0.02,0.007,0.003)); //afav_dstk_dp
  CorrData.push_back(dato(0.952,0.062,0.065));// rcp_dg
  CorrData.push_back(dato(1.051,0.022,0.028));// rcp_dp
  CorrData.push_back(dato(0.01170,0.02150,0.0313));// rm_dstk_dg
  CorrData.push_back(dato(0.0202,0.0035,0.0023));// rm_dstk_dp
  CorrData.push_back(dato(0.0292,0.0214,0.0312));//rp_dstk_dg
  CorrData.push_back(dato(0.0033,0.0035,0.0022));//rp_dstk_dp
  CorrData.push_back(dato(0.,0.014,0.006));//acp_dstpi_dg
  CorrData.push_back(dato(0.013,0.007,0.003));//acp_dstpi_dp
  CorrData.push_back(dato(0.00472,0.00092,0.00118));//rm_dstpi_dg
  CorrData.push_back(dato(0.00405,0.00056,0.00059));//rm_dstpi_dp
  CorrData.push_back(dato(0.00403,0.00091,0.00114));//rp_dstpi_dg
  CorrData.push_back(dato(0.00536,0.00056,0.00058));//rp_dstpi_dp

  //Correlation Matrix (stat):
  Corr.ResizeTo(16,16);
  Corr(0,0) = 1.;
  Corr(0,1) = -0.61;
  Corr(0,2) = 0.0;
  Corr(0,3) = 0.0;
  Corr(0,4) = -0.15;
  Corr(0,5) = 0.07;
  Corr(0,6) = -0.03;
  Corr(0,7) = 0.01;
  Corr(0,8) = -0.03;
  Corr(0,9) = 0.01;
  Corr(0,10) = -0.01;
  Corr(0,11) = -0.02;
  Corr(0,12) = -0.01;
  Corr(0,13) = 0.;
  Corr(0,14) = -0.02;
  Corr(0,15) = 0.;
  Corr(1,1) = 1.;
  Corr(1,2) = 0.0;
  Corr(1,3) = 0.01;
  Corr(1,4) = -0.05;
  Corr(1,5) = 0.08;
  Corr(1,6) = 0.0;
  Corr(1,7) = 0.;
  Corr(1,8) = 0.;
  Corr(1,9) = 0.;
  Corr(1,10) = 0.03;
  Corr(1,11) = 0.05;
  Corr(1,12) = 0.;
  Corr(1,13) = 0.;
  Corr(1,14) = 0.;
  Corr(1,15) = 0.;
  Corr(2,2) = 1.;
  Corr(2,3) = -0.59;
  Corr(2,4) = 0.;
  Corr(2,5) = 0.;
  Corr(2,6) = 0.0;
  Corr(2,7) = 0.0;
  Corr(2,8) = 0.;
  Corr(2,9) = 0.;
  Corr(2,10) = 0.;
  Corr(2,11) = 0.01;
  Corr(2,12) = 0.;
  Corr(2,13) = 0.;
  Corr(2,14) = 0.;
  Corr(2,15) = 0.;
  Corr(3,3) = 1.;
  Corr(3,4) = 0.;
  Corr(3,5) = 0.;
  Corr(3,6) = 0.;
  Corr(3,7) = 0.;
  Corr(3,8) = 0.;
  Corr(3,9) = 0.;
  Corr(3,10) = 0.01;
  Corr(3,11) = 0.02;
  Corr(3,12) = 0.;
  Corr(3,13) = 0.01;
  Corr(3,14) = 0.;
  Corr(3,15) = -0.01;
  Corr(4,4) = 1.;
  Corr(4,5) = -0.44;
  Corr(4,6) = 0.23;
  Corr(4,7) = -0.08;
  Corr(4,8) = 0.23;
  Corr(4,9) = -0.08;
  Corr(4,10) = 0.;
  Corr(4,11) = 0.;
  Corr(4,12) = 0.10;
  Corr(4,13) = 0.;
  Corr(4,14) = 0.1;
  Corr(4,15) = 0.;
  Corr(5,5) = 1.;
  Corr(5,6) = -0.04;
  Corr(5,7) = 0.03;
  Corr(5,8) = -0.04;
  Corr(5,9) = 0.02;
  Corr(5,10) = 0.;
  Corr(5,11) = 0.;
  Corr(5,12) = -0.02;
  Corr(5,13) = 0.;
  Corr(5,14) = -0.02;
  Corr(5,15) = 0.;
  Corr(6,6) = 1.;
  Corr(6,7) = -0.59;
  Corr(6,8) = 0.79;
  Corr(6,9) = -0.27;
  Corr(6,10) = 0.;
  Corr(6,11) = 0.;
  Corr(6,12) = 0.3;
  Corr(6,13) = -0.03;
  Corr(6,14) = 0.33;
  Corr(6,15) = -0.02;
  Corr(7,7) = 1.;
  Corr(7,8) = -0.27;
  Corr(7,9) = 0.1;
  Corr(7,10) = 0.;
  Corr(7,11) = 0.;
  Corr(7,12) = -0.07;
  Corr(7,13) = 0.05;
  Corr(7,14) = -0.1;
  Corr(7,15) = 0.05;
  Corr(8,8) = 1.;
  Corr(8,9) = -0.6;
  Corr(8,10) = 0.;
  Corr(8,11) = 0.;
  Corr(8,12) = 0.32;
  Corr(8,13) = -0.01;
  Corr(8,14) = 0.3;
  Corr(8,15) = -0.03;
  Corr(9,9) = 1.;
  Corr(9,10) = 0.;
  Corr(9,11) = 0.;
  Corr(9,12) = -0.09;
  Corr(9,13) = 0.04;
  Corr(9,14) = -0.07;
  Corr(9,15) = 0.05;
  Corr(10,10) = 1.;
  Corr(10,11) = 0.05;
  Corr(10,12) = 0.;
  Corr(10,13) = 0.;
  Corr(10,14) = 0.;
  Corr(10,15) = 0.;
  Corr(11,11) = 1.;
  Corr(11,12) = 0.;
  Corr(11,13) = 0.;
  Corr(11,14) = 0.;
  Corr(11,15) = -0.01;
  Corr(12,12) = 1.;
  Corr(12,13) = -0.11;
  Corr(12,14) = 0.33;
  Corr(12,15) = 0.19;
  Corr(13,13) = 1.;
  Corr(13,14) = 0.19;
  Corr(13,15) = 0.58;
  Corr(14,14) = 1.;
  Corr(14,15) = -0.11;
  Corr(15,15) = 1.;


  //Correlation Matrix (syst)
  Corr2.ResizeTo(16,16);
  Corr2(0,0) = 1.;
  Corr2(0,1) = -0.52;
  Corr2(0,2) = 0.50;
  Corr2(0,3) = -0.39;
  Corr2(0,4) = -0.31;
  Corr2(0,5) = -0.71;
  Corr2(0,6) = -0.15;
  Corr2(0,7) = 0.10;
  Corr2(0,8) = -0.05;
  Corr2(0,9) = 0.12;
  Corr2(0,10) = 0.74;
  Corr2(0,11) = -0.02;
  Corr2(0,12) = -0.15;
  Corr2(0,13) = -0.14;
  Corr2(0,14) = -0.05;
  Corr2(0,15) = -0.17;
  Corr2(1,1) = 1.;
  Corr2(1,2) = -0.44;
  Corr2(1,3) = -0.08;
  Corr2(1,4) = 0.33;
  Corr2(1,5) = 0.62;
  Corr2(1,6) = 0.19;
  Corr2(1,7) = -0.19;
  Corr2(1,8) = 0.16;
  Corr2(1,9) = -0.08;
  Corr2(1,10) = -0.49;
  Corr2(1,11) = 0.;
  Corr2(1,12) = 0.18;
  Corr2(1,13) = 0.01;
  Corr2(1,14) = 0.1;
  Corr2(1,15) = 0.03;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.23;
  Corr2(2,4) = 0.02;
  Corr2(2,5) = -0.35;
  Corr2(2,6) = -0.04;
  Corr2(2,7) = 0.03;
  Corr2(2,8) = 0.04;
  Corr2(2,9) = 0.05;
  Corr2(2,10) = 0.57;
  Corr2(2,11) = -0.38;
  Corr2(2,12) = -0.15;
  Corr2(2,13) = -0.19;
  Corr2(2,14) = 0.04;
  Corr2(2,15) = -0.09;
  Corr2(3,3) = 1.;
  Corr2(3,4) = -0.19;
  Corr2(3,5) = 0.14;
  Corr2(3,6) = -0.08;
  Corr2(3,7) = 0.06;
  Corr2(3,8) = -0.12;
  Corr2(3,9) = -0.01;
  Corr2(3,10) = -0.54;
  Corr2(3,11) = -0.64;
  Corr2(3,12) = -0.05;
  Corr2(3,13) = 0.04;
  Corr2(3,14) = 0.02;
  Corr2(3,15) = 0.24;
  Corr2(4,4) = 1.;
  Corr2(4,5) = 0.49;
  Corr2(4,6) = 0.36;
  Corr2(4,7) = -0.15;
  Corr2(4,8) = 0.41;
  Corr2(4,9) = -0.09;
  Corr2(4,10) = 0.22;
  Corr2(4,11) = 0.16;
  Corr2(4,12) = 0.16;
  Corr2(4,13) = -0.06;
  Corr2(4,14) = 0.19;
  Corr2(4,15) = -0.07;
  Corr2(5,5) = 1.;
  Corr2(5,6) = 0.14;
  Corr2(5,7) = -0.04;
  Corr2(5,8) = 0.10;
  Corr2(5,9) = -0.05;
  Corr2(5,10) = -0.35;
  Corr2(5,11) = 0.07;
  Corr2(5,12) = 0.15;
  Corr2(5,13) = 0.08;
  Corr2(5,14) = 0.09;
  Corr2(5,15) = 0.1;
  Corr2(6,6) = 1.;
  Corr2(6,7) = -0.55;
  Corr2(6,8) = 0.98;
  Corr2(6,9) = -0.51;
  Corr2(6,10) = 0.;
  Corr2(6,11) = 0.07;
  Corr2(6,12) = 0.37;
  Corr2(6,13) = -0.04;
  Corr2(6,14) = 0.4;
  Corr2(6,15) = -0.04;
  Corr2(7,7) = 1.;
  Corr2(7,8) = -0.53;
  Corr2(7,9) = 0.94;
  Corr2(7,10) = 0.04;
  Corr2(7,11) = -0.01;
  Corr2(7,12) = -0.13;
  Corr2(7,13) = 0.08;
  Corr2(7,14) = -0.15;
  Corr2(7,15) = 0.06;
  Corr2(8,8) = 1.;
  Corr2(8,9) = -0.48;
  Corr2(8,10) = 0.12;
  Corr2(8,11) = 0.04;
  Corr2(8,12) = 0.37;
  Corr2(8,13) = -0.06;
  Corr2(8,14) = 0.41;
  Corr2(8,15) = -0.06;
  Corr2(9,9) = 1.;
  Corr2(9,10) = 0.06;
  Corr2(9,11) = -0.03;
  Corr2(9,12) = -0.14;
  Corr2(9,13) = 0.02;
  Corr2(9,14) = -0.13;
  Corr2(9,15) = 0.04;
  Corr2(10,10) = 1.;
  Corr2(10,11) = 0.26;
  Corr2(10,12) = -0.1;
  Corr2(10,13) = -0.18;
  Corr2(10,14) = 0.02;
  Corr2(10,15) = -0.25;
  Corr2(11,11) = 1.;
  Corr2(11,12) = 0.09;
  Corr2(11,13) = 0.08;
  Corr2(11,14) = -0.07;
  Corr2(11,15) = -0.14;
  Corr2(12,12) = 1.;
  Corr2(12,13) = 0.53;
  Corr2(12,14) = 0.79;
  Corr2(12,15) = 0.16;
  Corr2(13,13) = 1.;
  Corr2(13,14) = 0.22;
  Corr2(13,15) = 0.60;
  Corr2(14,14) = 1.;
  Corr2(14,15) = 0.41;
  Corr2(15,15) = 1.;


  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID5", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //7. PDF: glwads-dkst-hh-h3pi-dmix (UID6)
  //Observables 12:
  CorrData.clear();
  CorrData.push_back(dato(-0.004,0.023,0.008));//afav_dkst_kpi
  CorrData.push_back(dato(0.06,0.07,0.01));//acp_dkst_kk
  CorrData.push_back(dato(0.15,0.13,0.01));//acp_dkst_pipi
  CorrData.push_back(dato(1.22,0.09,0.02));//rcp_dkst_kk
  CorrData.push_back(dato(1.08,0.14,0.03));//rcp_dkst_pipi
  CorrData.push_back(dato(0.02,0.006,0.001));// rp_dkst_kpi
  CorrData.push_back(dato(0.002,0.004,0.001));// rm_dkst_kpi
  CorrData.push_back(dato(-0.013,0.031,0.009));// afav_dkst_k3pi
  CorrData.push_back(dato(0.02,0.11,0.01));//acp_dkst_pipipipi
  CorrData.push_back(dato(1.08,0.13,0.03));//rcp_dkst_pipipipi
  CorrData.push_back(dato(0.016,0.007,0.003));// rp_dkst_k3pi
  CorrData.push_back(dato(0.006,0.006,0.004));//rm_dkst_k3pi

  //Correlation Matrix (stat):
  Corr.ResizeTo(12,12);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.;
  Corr(0,2) = 0.0;
  Corr(0,3) = 0.0;
  Corr(0,4) = 0.0;
  Corr(0,5) = 0.08;
  Corr(0,6) = -0.01;
  Corr(0,7) = 0.;
  Corr(0,8) = 0.;
  Corr(0,9) = 0.;
  Corr(0,10) = 0.;
  Corr(0,11) = 0.;
  Corr(1,1) = 1.;
  Corr(1,2) = 0.0;
  Corr(1,3) = 0.0;
  Corr(1,4) = 0.0;
  Corr(1,5) = 0.0;
  Corr(1,6) = 0.0;
  Corr(1,7) = 0.;
  Corr(1,8) = 0.;
  Corr(1,9) = 0.;
  Corr(1,10) = 0.;
  Corr(1,11) = 0.;
  Corr(2,2) = 1.;
  Corr(2,3) = 0.;
  Corr(2,4) = -0.02;
  Corr(2,5) = 0.0;
  Corr(2,6) = 0.0;
  Corr(2,7) = 0.0;
  Corr(2,8) = 0.;
  Corr(2,9) = 0.;
  Corr(2,10) = 0.0;
  Corr(2,11) = 0.;
  Corr(3,3) = 1.;
  Corr(3,4) = 0.05;
  Corr(3,5) = 0.02;
  Corr(3,6) = -0.01;
  Corr(3,7) = 0.;
  Corr(3,8) = 0.;
  Corr(3,9) = 0.;
  Corr(3,10) = 0.;
  Corr(3,11) = 0.;
  Corr(4,4) = 1.;
  Corr(4,5) = 0.03;
  Corr(4,6) = 0.02;
  Corr(4,7) = 0.;
  Corr(4,8) = 0.;
  Corr(4,9) = 0.;
  Corr(4,10) = 0.;
  Corr(4,11) = 0.;
  Corr(5,5) = 1.;
  Corr(5,6) = 0.02;
  Corr(5,7) = 0.;
  Corr(5,8) = 0.;
  Corr(5,9) = 0.;
  Corr(5,10) = 0.;
  Corr(5,11) = 0.;
  Corr(6,6) = 1.;
  Corr(6,7) = 0.;
  Corr(6,8) = 0.;
  Corr(6,9) = 0.;
  Corr(6,10) = 0.;
  Corr(6,11) = 0.;
  Corr(7,7) = 1.;
  Corr(7,8) = 0.0;
  Corr(7,9) = 0.;
  Corr(7,10) = 0.07;
  Corr(7,11) = -0.03;
  Corr(8,8) = 1.;
  Corr(8,9) = 0.01;
  Corr(8,10) = 0.;
  Corr(8,11) = 0.;
  Corr(9,9) = 1.;
  Corr(9,10) = 0.04;
  Corr(9,11) = 0.04;
  Corr(10,10) = 1.;
  Corr(10,11) = 0.03;
  Corr(11,11) = 1.;


  //Correlation Matrix (syst)
  Corr2.ResizeTo(12,12);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.82;
  Corr2(0,2) = 0.72;
  Corr2(0,3) = 0.0;
  Corr2(0,4) = 0.;
  Corr2(0,5) = 0.01;
  Corr2(0,6) = -0.02;
  Corr2(0,7) = 0.94;
  Corr2(0,8) = 0.84;
  Corr2(0,9) = 0.;
  Corr2(0,10) = -0.01;
  Corr2(0,11) = 0.;
  Corr2(1,1) = 1.;
  Corr2(1,2) = 0.65;
  Corr2(1,3) = -0.04;
  Corr2(1,4) = 0.02;
  Corr2(1,5) = 0.01;
  Corr2(1,6) = -0.02;
  Corr2(1,7) = 0.83;
  Corr2(1,8) = 0.77;
  Corr2(1,9) = 0.;
  Corr2(1,10) = 0.;
  Corr2(1,11) = 0.;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.;
  Corr2(2,4) = -0.03;
  Corr2(2,5) = 0.0;
  Corr2(2,6) = -0.02;
  Corr2(2,7) = 0.72;
  Corr2(2,8) = 0.68;
  Corr2(2,9) = 0.;
  Corr2(2,10) = 0.;
  Corr2(2,11) = 0.01;
  Corr2(3,3) = 1.;
  Corr2(3,4) = 0.0;
  Corr2(3,5) = 0.05;
  Corr2(3,6) = 0.03;
  Corr2(3,7) = -0.01;
  Corr2(3,8) = 0.;
  Corr2(3,9) = -0.01;
  Corr2(3,10) = -0.01;
  Corr2(3,11) = -0.01;
  Corr2(4,4) = 1.;
  Corr2(4,5) = 0.06;
  Corr2(4,6) = 0.08;
  Corr2(4,7) = -0.01;
  Corr2(4,8) = 0.0;
  Corr2(4,9) = -0.01;
  Corr2(4,10) = -0.02;
  Corr2(4,11) = 0.01;
  Corr2(5,5) = 1.;
  Corr2(5,6) = 0.08;
  Corr2(5,7) = -0.01;
  Corr2(5,8) = 0.;
  Corr2(5,9) = 0.;
  Corr2(5,10) = -0.01;
  Corr(5,11) = -0.01;
  Corr2(6,6) = 1.;
  Corr2(6,7) = -0.01;
  Corr2(6,8) = -0.01;
  Corr2(6,9) = -0.01;
  Corr2(6,10) = 0.01;
  Corr2(6,11) = 0.03;
  Corr2(7,7) = 1.;
  Corr2(7,8) = 0.84;
  Corr2(7,9) = 0.;
  Corr2(7,10) = -0.01;
  Corr2(7,11) = -0.02;
  Corr2(8,8) = 1.;
  Corr2(8,9) = 0.03;
  Corr2(8,10) = 0.01;
  Corr2(8,11) = 0.;
  Corr2(9,9) = 1.;
  Corr2(9,10) = 0.01;
  Corr2(9,11) = -0.01;
  Corr2(10,10) = 1.;
  Corr2(10,11) = 0.05;
  Corr2(11,11) = 1.;


  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID6", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //8. PDF: glwads-dkst-hh-h3pi-dmix (UID7)
  //Observables 12:
  CorrData.clear();
  CorrData.push_back(dato(-0.05,0.1,0.01));//acp_dkstz_kk
  CorrData.push_back(dato(-0.18,0.14,0.01));//acp_dkst_pipi
  CorrData.push_back(dato(0.92,0.1,0.02));//rcp_dkstz_kk
  CorrData.push_back(dato(1.32,0.19,0.03));//rcp_dkst_pipi
  CorrData.push_back(dato(-0.03,0.15,0.01));//acp_dkstz_4pi
  CorrData.push_back(dato(1.01,0.16,0.04));// rcp_dkstz_4pi
  CorrData.push_back(dato(0.064,0.021,0.002));// rp_dkstz_kpi
  CorrData.push_back(dato(0.095,0.021,0.003));// rm_dkstz_kpi
  CorrData.push_back(dato(0.074,0.026,0.002));// rp_dkstz_k3pi
  CorrData.push_back(dato(0.072,0.025,0.003));// rm_dkstz_k3pi
  CorrData.push_back(dato(0.047,0.027,0.01));// afav_dkstz_kpi
  CorrData.push_back(dato(0.037,0.032,0.01));//afav_dkstz_k3pi

  //Correlation Matrix (stat):
  Corr.ResizeTo(12,12);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.;
  Corr(0,2) = 0.03;
  Corr(0,3) = 0.0;
  Corr(0,4) = 0.0;
  Corr(0,5) = 0.;
  Corr(0,6) = 0.;
  Corr(0,7) = 0.;
  Corr(0,8) = 0.;
  Corr(0,9) = 0.;
  Corr(0,10) = 0.;
  Corr(0,11) = 0.;
  Corr(1,1) = 1.;
  Corr(1,2) = 0.0;
  Corr(1,3) = 0.07;
  Corr(1,4) = 0.0;
  Corr(1,5) = 0.0;
  Corr(1,6) = 0.0;
  Corr(1,7) = 0.;
  Corr(1,8) = 0.;
  Corr(1,9) = 0.;
  Corr(1,10) = 0.;
  Corr(1,11) = 0.;
  Corr(2,2) = 1.;
  Corr(2,3) = 0.05;
  Corr(2,4) = 0.;
  Corr(2,5) = 0.04;
  Corr(2,6) = 0.02;
  Corr(2,7) = 0.03;
  Corr(2,8) = 0.01;
  Corr(2,9) = 0.01;
  Corr(2,10) = 0.0;
  Corr(2,11) = 0.;
  Corr(3,3) = 1.;
  Corr(3,4) = 0.;
  Corr(3,5) = 0.03;
  Corr(3,6) = 0.01;
  Corr(3,7) = 0.02;
  Corr(3,8) = 0.01;
  Corr(3,9) = 0.01;
  Corr(3,10) = 0.;
  Corr(3,11) = 0.;
  Corr(4,4) = 1.;
  Corr(4,5) = 0.01;
  Corr(4,6) = 0.;
  Corr(4,7) = 0.;
  Corr(4,8) = 0.;
  Corr(4,9) = 0.;
  Corr(4,10) = 0.;
  Corr(4,11) = 0.;
  Corr(5,5) = 1.;
  Corr(5,6) = 0.01;
  Corr(5,7) = 0.01;
  Corr(5,8) = 0.02;
  Corr(5,9) = 0.02;
  Corr(5,10) = 0.;
  Corr(5,11) = 0.;
  Corr(6,6) = 1.;
  Corr(6,7) = 0.05;
  Corr(6,8) = 0.01;
  Corr(6,9) = 0.01;
  Corr(6,10) = 0.09;
  Corr(6,11) = 0.;
  Corr(7,7) = 1.;
  Corr(7,8) = 0.01;
  Corr(7,9) = 0.01;
  Corr(7,10) = -0.12;
  Corr(7,11) = 0.;
  Corr(8,8) = 1.;
  Corr(8,9) = 0.06;
  Corr(8,10) = 0.;
  Corr(8,11) = 0.09;
  Corr(9,9) = 1.;
  Corr(9,10) = 0.;
  Corr(9,11) = -0.09;
  Corr(10,10) = 1.;
  Corr(10,11) = 0.;
  Corr(11,11) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(12,12);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.72;
  Corr2(0,2) = 0.31;
  Corr2(0,3) = -0.65;
  Corr2(0,4) = 0.31;
  Corr2(0,5) = -0.36;
  Corr2(0,6) = -0.15;
  Corr2(0,7) = -0.63;
  Corr2(0,8) = -0.67;
  Corr2(0,9) = -0.66;
  Corr2(0,10) = -0.51;
  Corr2(0,11) = -0.46;
  Corr2(1,1) = 1.;
  Corr2(1,2) = 0.39;
  Corr2(1,3) = -0.77;
  Corr2(1,4) = 0.19;
  Corr2(1,5) = -0.33;
  Corr2(1,6) = -0.16;
  Corr2(1,7) = -0.77;
  Corr2(1,8) = -0.79;
  Corr2(1,9) = -0.75;
  Corr2(1,10) = -0.65;
  Corr2(1,11) = -0.6;
  Corr2(2,2) = 1.;
  Corr2(2,3) = -0.45;
  Corr2(2,4) = 0.09;
  Corr2(2,5) = -0.15;
  Corr2(2,6) = 0.;
  Corr2(2,7) = -0.49;
  Corr2(2,8) = -0.56;
  Corr2(2,9) = -0.54;
  Corr2(2,10) = -0.46;
  Corr2(2,11) = -0.47;
  Corr2(3,3) = 1.;
  Corr2(3,4) = -0.15;
  Corr2(3,5) = 0.32;
  Corr2(3,6) = 0.1;
  Corr2(3,7) = 0.81;
  Corr2(3,8) = 0.81;
  Corr2(3,9) = 0.81;
  Corr2(3,10) = 0.71;
  Corr2(3,11) = 0.68;
  Corr2(4,4) = 1.;
  Corr2(4,5) = -0.1;
  Corr2(4,6) = -0.02;
  Corr2(4,7) = -0.16;
  Corr2(4,8) = -0.14;
  Corr2(4,9) = -0.2;
  Corr2(4,10) = 0.01;
  Corr2(4,11) = 0.1;
  Corr2(5,5) = 1.;
  Corr2(5,6) = 0.06;
  Corr2(5,7) = 0.33;
  Corr2(5,8) = 0.35;
  Corr2(5,9) = 0.37;
  Corr2(5,10) = 0.22;
  Corr(5,11) = 0.22;
  Corr2(6,6) = 1.;
  Corr2(6,7) = 0.1;
  Corr2(6,8) = 0.12;
  Corr2(6,9) = 0.12;
  Corr2(6,10) = 0.08;
  Corr2(6,11) = 0.06;
  Corr2(7,7) = 1.;
  Corr2(7,8) = 0.83;
  Corr2(7,9) = 0.81;
  Corr2(7,10) = 0.72;
  Corr2(7,11) = 0.69;
  Corr2(8,8) = 1.;
  Corr2(8,9) = 0.85;
  Corr2(8,10) = 0.75;
  Corr2(8,11) = 0.75;
  Corr2(9,9) = 1.;
  Corr2(9,10) = 0.72;
  Corr2(9,11) = 0.69;
  Corr2(10,10) = 1.;
  Corr2(10,11) = 0.72;
  Corr2(11,11) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID7", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //9. PDF: ggsz-dkstz-lhcb-md (UID8)
  //Observables 4:
  CorrData.clear();
  CorrData.push_back(dato(-0.15,0.14,0.03162));//xm_dkstz
  CorrData.push_back(dato(0.25,0.15,0.06083));//ym_dkstz
  CorrData.push_back(dato(0.05,0.24,0.04123));//xp_dkstz
  CorrData.push_back(dato(-0.65,0.235,0.08062));//yp_dkstz

  //Correlation Matrix (stat):
  Corr.ResizeTo(4,4);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.14;
  Corr(0,2) = 0.;
  Corr(0,3) = 0.;
  Corr(1,1) = 1.;
  Corr(1,2) = 0.;
  Corr(1,3) = 0.;
  Corr(2,2) = 1.;
  Corr(2,3) = 0.14;
  Corr(3,3) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(4,4);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.;
  Corr2(0,2) = 0.;
  Corr2(0,3) = 0.;
  Corr2(1,1) = 1.;
  Corr2(1,2) = 0.;
  Corr2(1,3) = 0.;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.;
  Corr2(3,3) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID8", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //10. PDF: glwads-dhpipi-hh-dmix (UID9)
  //Observables 11:
  CorrData.clear();
  CorrData.push_back(dato(1.0400,0.064,0.0));//rcp_dkpipi
  CorrData.push_back(dato(0.013,0.019,0.013));//afav_dkpipi_kpi
  CorrData.push_back(dato(-0.002,0.003,0.011));//afav_dpipipi_kpi
  CorrData.push_back(dato(-0.045,0.064,0.011));//acp_dkpipi_kk
  CorrData.push_back(dato(-0.054,0.101,0.011));// acp_dkpipi_pipi
  CorrData.push_back(dato(-0.019,0.011,0.01));//acp_dpipipi_kk
  CorrData.push_back(dato(-0.013,0.016,0.01));// acp_dpipipi_pipi
  CorrData.push_back(dato(0.0107,0.006,0.00110));// rp_dkpipi
  CorrData.push_back(dato(0.00530,0.0045,0.0006));// rm_dkpipi
  CorrData.push_back(dato(0.00432,0.00053,0.00021));// rp_dpipipi
  CorrData.push_back(dato(0.00421,0.00053,0.00021));// rm_dpipipi

  //Correlation Matrix (stat):
  Corr.ResizeTo(11,11);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.;
  Corr(0,2) = 0.0;
  Corr(0,3) = 0.0;
  Corr(0,4) = 0.0;
  Corr(0,5) = 0.0;
  Corr(0,6) = 0.0;
  Corr(0,7) = 0.;
  Corr(0,8) = 0.;
  Corr(0,9) = 0.;
  Corr(0,10) = 0.;
  Corr(1,1) = 1.;
  Corr(1,2) = 0.0;
  Corr(1,3) = 0.0;
  Corr(1,4) = 0.0;
  Corr(1,5) = 0.0;
  Corr(1,6) = 0.0;
  Corr(1,7) = 0.;
  Corr(1,8) = 0.;
  Corr(1,9) = 0.;
  Corr(1,10) = 0.;
  Corr(2,2) = 1.;
  Corr(2,3) = 0.;
  Corr(2,4) = 0.;
  Corr(2,5) = 0.0;
  Corr(2,6) = 0.0;
  Corr(2,7) = 0.0;
  Corr(2,8) = 0.;
  Corr(2,9) = 0.;
  Corr(2,10) = 0.;
  Corr(3,3) = 1.;
  Corr(3,4) = 0.2;
  Corr(3,5) = 0.;
  Corr(3,6) = 0.;
  Corr(3,7) = 0.;
  Corr(3,8) = 0.;
  Corr(3,9) = 0.;
  Corr(3,10) = 0.;
  Corr(4,4) = 1.;
  Corr(4,5) = 0.0;
  Corr(4,6) = 0.0;
  Corr(4,7) = 0.;
  Corr(4,8) = 0.;
  Corr(4,9) = 0.;
  Corr(4,10) = 0.;
  Corr(5,5) = 1.;
  Corr(5,6) = 0.08;
  Corr(5,7) = 0.;
  Corr(5,8) = 0.;
  Corr(5,9) = 0.;
  Corr(5,10) = 0.;
  Corr(6,6) = 1.;
  Corr(6,7) = 0.;
  Corr(6,8) = 0.;
  Corr(6,9) = 0.;
  Corr(6,10) = 0.;
  Corr(7,7) = 1.;
  Corr(7,8) = 0.0;
  Corr(7,9) = 0.;
  Corr(7,10) = 0.0;
  Corr(8,8) = 1.;
  Corr(8,9) = 0.;
  Corr(8,10) = 0.;
  Corr(9,9) = 1.;
  Corr(9,10) = 0.0;
  Corr(10,10) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(11,11);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.;
  Corr2(0,2) = 0.;
  Corr2(0,3) = 0.;
  Corr2(0,4) = 0.;
  Corr2(0,5) = 0.;
  Corr2(0,6) = 0.;
  Corr2(0,7) = 0.;
  Corr2(0,8) = 0.;
  Corr2(0,9) = 0.;
  Corr2(0,10) = 0.;
  Corr2(1,1) = 1.;
  Corr2(1,2) = 0.;
  Corr2(1,3) = 0.;
  Corr2(1,4) = 0.;
  Corr2(1,5) = 0.;
  Corr2(1,6) = 0.;
  Corr2(1,7) = 0.;
  Corr2(1,8) = 0.;
  Corr2(1,9) = 0.;
  Corr2(1,10) = 0.;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.;
  Corr2(2,4) = 0.;
  Corr2(2,5) = 0.;
  Corr2(2,6) = 0.;
  Corr2(2,7) = 0.;
  Corr2(2,8) = 0.;
  Corr2(2,9) = 0.;
  Corr2(2,10) = 0.;
  Corr2(3,3) = 1.;
  Corr2(3,4) = 0.;
  Corr2(3,5) = 0.;
  Corr2(3,6) = 0.;
  Corr2(3,7) = 0.;
  Corr2(3,8) = 0.;
  Corr2(3,9) = 0.;
  Corr2(3,10) = 0.;
  Corr2(4,4) = 1.;
  Corr2(4,5) = 0.;
  Corr2(4,6) = 0.;
  Corr2(4,7) = 0.;
  Corr2(4,8) = 0.;
  Corr2(4,9) = 0.;
  Corr2(4,10) = 0.;
  Corr2(5,5) = 1.;
  Corr2(5,6) = 0.;
  Corr2(5,7) = 0.;
  Corr2(5,8) = 0.;
  Corr2(5,9) = 0.;
  Corr2(5,10) = 0.;
  Corr2(6,6) = 1.;
  Corr2(6,7) = 0.;
  Corr2(6,8) = 0.;
  Corr2(6,9) = 0.;
  Corr2(6,10) = 0.;
  Corr2(7,7) = 1.;
  Corr2(7,8) = 0.;
  Corr2(7,9) = 0.;
  Corr2(7,10) = 0.;
  Corr2(8,8) = 1.;
  Corr2(8,9) = 0.;
  Corr2(8,10) = 0.;
  Corr2(9,9) = 1.;
  Corr2(9,10) = 0.;
  Corr2(10,10) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID9", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //-------------------------------------------  Other useful inputs for the parameters of the Time Integrated B Decay Chain -------------------------------------------------------------------------

  //20. PDF: dk3pi_dkpipi0_constraints (UID19)
  //Observables 6:
  CorrData.clear();
  CorrData.push_back(dato(0.44,0.095,0.));//kD_k3pi
  CorrData.push_back(dato(2.80998, 0.40143 , 0.));//dD_k3pi
  //CorrData.push_back(dato(2*M_PI - 2.80998, 0.40143 , 0.));//dD_k3pi
  //NB: I used 2pi-x, where x is the value written in gammacharm_lhc document because I have the phase convention dD = -dD
  CorrData.push_back(dato(0.79,0.04,0.));// kD_kpipi0
  CorrData.push_back(dato(3.42085, 0.19199,0.));// dD_kpipi0
  //CorrData.push_back(dato(2 * M_PI - 3.42085, 0.19199,0.));// dD_kpipi0
  //NB: I used 2pi-x, where x is the value written in gammacharm_lhc document because I have the phase convention dD = -dD
  CorrData.push_back(dato(0.055,0.0007,0.));// rD_k3pi
  CorrData.push_back(dato(0.0441,0.0011,0.));// rD_kpipi0

  //Correlation Matrix (stat):
  Corr.ResizeTo(6,6);
  Corr(0,0) = 1.;
  Corr(0,1) = -0.75;
  Corr(0,2) = 0.;
  Corr(0,3) = -0.07;
  Corr(0,4) = 0.52;
  Corr(0,5) = -0.06;
  Corr(1,1) = 1.;
  Corr(1,2) = 0.03;
  Corr(1,3) = 0.17;
  Corr(1,4) = -0.42;
  Corr(1,5) = 0.01;
  Corr(2,2) = 1.;
  Corr(2,3) = 0.19;
  Corr(2,4) = -0.01;
  Corr(2,5) = -0.01;
  Corr(3,3) = 1.;
  Corr(3,4) = -0.02;
  Corr(3,5) = 0.25;
  Corr(4,4) = 1.;
  Corr(4,5) = -0.12;
  Corr(5,5) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(6,6);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.;
  Corr2(0,2) = 0.;
  Corr2(0,3) = 0.;
  Corr2(0,4) = 0.;
  Corr2(0,5) = 0.;
  Corr2(1,1) = 1.;
  Corr2(1,2) = 0.;
  Corr2(1,3) = 0.;
  Corr2(1,4) = 0.;
  Corr2(1,5) = 0.;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.;
  Corr2(2,4) = 0.;
  Corr2(2,5) = 0.;
  Corr2(3,3) = 1.;
  Corr2(3,4) = 0.;
  Corr2(3,5) = 0.;
  Corr2(4,4) = 1.;
  Corr2(4,5) = 0.;
  Corr2(5,5) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID19", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //21. PDF: d4pi_dmixing_cleo (UID20)

  //Observables 1:
  meas.insert(pair<string,dato>("UID20", dato(0.737,0.028,0.))); //Fpipipipi

  //22. PDF: CleoDhhpi0Diluition (UID21)

  //Observables 2:
  CorrData.clear();
  CorrData.push_back(dato(0.973,0.017,0.));//F_pipipi0
  CorrData.push_back(dato(0.732,0.055,0.));//F_kkpi0

  //Correlation Matrix (stat):
  Corr.ResizeTo(2,2);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.;
  Corr(1,1) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(2,2);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.;
  Corr2(1,1) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID21", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //25. PDF: dkstcoherence (UID24)

  //Observables 1:
  meas.insert(pair<string,dato>("UID24", dato(0.95,0.06,0.)));//k_dkst

  //26. PDF: dkstzcoherence (UID25)

  //Observables 1:
  meas.insert(pair<string,dato>("UID25", dato(0.958,0.0075,0.024)));//k_dkstz

}
// ---------------------------------------------------------

// ---------------------------------------------------------
void MixingModel::Add_time_dependent_Bmeas(){

  vector<dato> CorrData;
  TMatrixDSym Corr, Corr2;

  //-------------------------------------- Time Dependent B0 decay-------------------------------------------------------------------------

  //11. PDF: dsk (UID10)

  //Observables 5:
  CorrData.clear();
  CorrData.push_back(dato(0.73,0.14,0.05));//c_dsk
  CorrData.push_back(dato(0.39,0.28,0.15));// d_dsk
  CorrData.push_back(dato(0.31,0.28,0.15));//db_dsk
  CorrData.push_back(dato(-0.52,0.2,0.07));//s_dsk
  CorrData.push_back(dato(-0.49,0.2,0.07));// sb_dsk

  //Correlation Matrix (stat):
  Corr.ResizeTo(5,5);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.09;
  Corr(0,2) = 0.08;
  Corr(0,3) = 0.01;
  Corr(0,4) = -0.06;
  Corr(1,1) = 1.;
  Corr(1,2) = 0.51;
  Corr(1,3) = -0.08;
  Corr(1,4) = 0.;
  Corr(2,2) = 1.;
  Corr(2,3) = -0.04;
  Corr(2,4) = 0.;
  Corr(3,3) = 1.;
  Corr(3,4) = 0.;
  Corr(4,4) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(5,5);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.05;
  Corr2(0,2) = 0.03;
  Corr2(0,3) = 0.03;
  Corr2(0,4) = -0.01;
  Corr2(1,1) = 1.;
  Corr2(1,2) = 0.42;
  Corr2(1,3) = 0.02;
  Corr2(1,4) = 0.02;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.03;
  Corr2(2,4) = 0.03;
  Corr2(3,3) = 1.;
  Corr2(3,4) = 0.01;
  Corr2(4,4) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID10", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //12. PDF: dskpipi (UID11)

  //Observables 5:
  CorrData.clear();
  CorrData.push_back(dato(0.631,0.096,0.032));//c_dskpipi
  CorrData.push_back(dato(-0.334,0.232,0.097));// d_dskpipi
  CorrData.push_back(dato(-0.695,0.215,0.081));//db_dskpipi
  CorrData.push_back(dato(-0.424,0.135,0.033));//s_dskpipi
  CorrData.push_back(dato(-0.463,0.134,0.031));// sb_dskpipi

  //Correlation Matrix (stat):
  Corr.ResizeTo(5,5);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.14;
  Corr(0,2) = 0.17;
  Corr(0,3) = 0.04;
  Corr(0,4) = 0.04;
  Corr(1,1) = 1.;
  Corr(1,2) = 0.50;
  Corr(1,3) = -0.08;
  Corr(1,4) = -0.05;
  Corr(2,2) = 1.;
  Corr(2,3) = -0.04;
  Corr(2,4) = -0.11;
  Corr(3,3) = 1.;
  Corr(3,4) = 0.01;
  Corr(4,4) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(5,5);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.04;
  Corr2(0,2) = 0.04;
  Corr2(0,3) = -0.14;
  Corr2(0,4) = -0.01;
  Corr2(1,1) = 1.;
  Corr2(1,2) = 0.21;
  Corr2(1,3) = -0.03;
  Corr2(1,4) = -0.02;
  Corr2(2,2) = 1.;
  Corr2(2,3) = -0.02;
  Corr2(2,4) = -0.04;
  Corr2(3,3) = 1.;
  Corr2(3,4) = -0.28;
  Corr2(4,4) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID11", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //13. PDF: dmpi (UID12)

  //Observables 2:
  CorrData.clear();
  CorrData.push_back(dato(0.058,0.02,0.011));//s_dmpi
  CorrData.push_back(dato(0.038,0.02,0.007));// sb_dmpi
  //changed the sign of the observables with respect to the document

  //Correlation Matrix (stat):
  Corr.ResizeTo(2,2);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.6;
  Corr(1,1) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(2,2);
  Corr2(0,0) = 1.;
  Corr2(0,1) = -0.41;
  Corr2(1,1) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID12", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //-------------------------------------  Other useful inputs for the Time Dependent B0 Decay and mixing parameters -------------------------------------------------------------------------

  //27. PDF: phis (UID26)

  //Observables 1:
  meas.insert(pair<string,dato>("UID26", dato(-0.05,0.019,0.)));//phis = - 2 beta_s

  //28. PDF: beta (UID27)

  //Observables 1:
  meas.insert(pair<string,dato>("UID27", dato(0.38746,0.0122,0.)));// beta

}
// ---------------------------------------------------------

// ---------------------------------------------------------
void MixingModel::Add_time_dependent_Dmeas(){

  vector<dato> CorrData;
  TMatrixDSym Corr, Corr2;

  //-------------------------------------- Time Dependent D0 decay-------------------------------------------------------------------------

  //15. PDF: charm-kskpi (UID14)

  //Observables 4:
  CorrData.clear();
  CorrData.push_back(dato(0.0027,0.0016,0.0004));//xcp
  CorrData.push_back(dato(0.00740,0.0036,0.0011));//ycp
  CorrData.push_back(dato(-0.00053,0.0007,0.00022));// dx
  CorrData.push_back(dato(0.0006,0.0016,0.0003));//dy

  //Correlation Matrix (stat):
  Corr.ResizeTo(4,4);
  Corr(0,0) = 1.;
  Corr(0,1) = -0.17;
  Corr(0,2) = 0.04;
  Corr(0,3) = -0.02;
  Corr(1,1) = 1.;
  Corr(1,2) = -0.03;
  Corr(1,3) = 0.01;
  Corr(2,2) = 1.;
  Corr(2,3) = -0.13;
  Corr(3,3) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(4,4);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.15;
  Corr2(0,2) = 0.01;
  Corr2(0,3) = -0.02;
  Corr2(1,1) = 1.;
  Corr2(1,2) = -0.05;
  Corr2(1,3) = -0.03;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.14;
  Corr2(3,3) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID14", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //16. PDF: charm-kskpi (UID15)

  //Observables 4:
  CorrData.clear();
  CorrData.push_back(dato(0.00397,0.00046,0.00029));//xcp
  CorrData.push_back(dato(0.00459,0.00120,0.00085));//ycp
  CorrData.push_back(dato(-0.00027,0.00018,0.00001));// dx
  CorrData.push_back(dato(0.00020,0.00036,0.00013));//dy

  //Correlation Matrix (stat):
  Corr.ResizeTo(4,4);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.11;
  Corr(0,2) = -0.02;
  Corr(0,3) = -0.01;
  Corr(1,1) = 1.;
  Corr(1,2) = -0.03;
  Corr(1,3) = -0.05;
  Corr(2,2) = 1.;
  Corr(2,3) = 0.08;
  Corr(3,3) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(4,4);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.13;
  Corr2(0,2) = 0.01;
  Corr2(0,3) = 0.01;
  Corr2(1,1) = 1.;
  Corr2(1,2) = -0.02;
  Corr2(1,3) = 0.01;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.31;
  Corr2(3,3) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID15", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //17. PDF: charm-kpi (UID16)
  //Observables 6:
  CorrData.clear();
  CorrData.push_back(dato(0.00345,0.00004,0.00002));//Rdp
  CorrData.push_back(dato(0.00501,0.00064,0.00038));//yp
  CorrData.push_back(dato(0.00006,0.00003,0.00002));// xpsq
  CorrData.push_back(dato(0.00345,0.00004,0.00002));// Rdm
  CorrData.push_back(dato(0.00554,0.00064,0.00038));// ym
  CorrData.push_back(dato(0.00002,0.00003,0.00002));// xmsq

  //Correlation Matrix (stat):
  Corr.ResizeTo(6,6);
  Corr(0,0) = 1.;
  Corr(0,1) = -0.94;
  Corr(0,2) = 0.84;
  Corr(0,3) = -0.01;
  Corr(0,4) = 0.;
  Corr(0,5) = 0.;
  Corr(1,1) = 1.;
  Corr(1,2) = -0.96;
  Corr(1,3) = 0.;
  Corr(1,4) = 0.;
  Corr(1,5) = 0.;
  Corr(2,2) = 1.;
  Corr(2,3) = 0.;
  Corr(2,4) = 0.;
  Corr(2,5) = 0.;
  Corr(3,3) = 1.;
  Corr(3,4) = -0.94;
  Corr(3,5) = 0.85;
  Corr(4,4) = 1.;
  Corr(4,5) = -0.96;
  Corr(5,5) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(6,6);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.;
  Corr2(0,2) = 0.;
  Corr2(0,3) = 0.;
  Corr2(0,4) = 0.;
  Corr2(0,5) = 0.;
  Corr2(1,1) = 1.;
  Corr2(1,2) = 0.;
  Corr2(1,3) = 0.;
  Corr2(1,4) = 0.;
  Corr2(1,5) = 0.;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.;
  Corr2(2,4) = 0.;
  Corr2(2,5) = 0.;
  Corr2(3,3) = 1.;
  Corr2(3,4) = 0.;
  Corr2(3,5) = 0.;
  Corr2(4,4) = 1.;
  Corr2(4,5) = 0.;
  Corr2(5,5) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID16", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //18. PDF: charm-deltaacp-diff (UID17)

  //Observables 2:
  CorrData.clear();
  CorrData.push_back(dato(0.115,0.002,0.));// dtotau
  CorrData.push_back(dato(-0.00154,0.00029,0.));//dacp

  //Correlation Matrix (stat):
  Corr.ResizeTo(2,2);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.;
  Corr(1,1) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(2,2);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.;
  Corr2(1,1) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID17", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //29. PDF: Charm_ycp_minus_ycp_rs_lhcb-r1 (UID28)

  //Observables 1:
  if(comb == 0){
    meas.insert(pair<string,dato>("UID28", dato(0.0057,0.0013,0.0009)));// ycp //old measurement
  }
  else{
    meas.insert(pair<string,dato>("UID28", dato(0.00696,0.00026,0.00013)));// ycp  // new measurement https://arxiv.org/abs/2202.09106
  } //global combination

  //30. PDF: charm-dy-rs (UID29)

  //Observables 1:
  meas.insert(pair<string,dato>("UID29", dato(-0.0001,0.00011,0.00003)));// DY

  //31. PDF: charm-kpi (UID30)
  //Observables 6:
  CorrData.clear();
  CorrData.push_back(dato(0.00338,0.00015,0.00006));//Rdp
  CorrData.push_back(dato(0.00581,0.00525,0.00032));//yp
  CorrData.push_back(dato(-0.00002,0.00045,0.00003));// xpsq
  CorrData.push_back(dato(0.0036,0.00015,0.00007));// Rdm
  CorrData.push_back(dato(0.00332,0.00521,0.0004));// ym
  CorrData.push_back(dato(0.00008,0.00043,0.00004));// xmsq

  //Correlation Matrix (stat):
  Corr.ResizeTo(6,6);
  Corr(0,0) = 1.;
  Corr(0,1) = -0.92;
  Corr(0,2) = 0.82;
  Corr(0,3) = 0.;
  Corr(0,4) = 0.;
  Corr(0,5) = 0.;
  Corr(1,1) = 1.;
  Corr(1,2) = -0.96;
  Corr(1,3) = 0.;
  Corr(1,4) = 0.;
  Corr(1,5) = 0.;
  Corr(2,2) = 1.;
  Corr(2,3) = 0.;
  Corr(2,4) = 0.;
  Corr(2,5) = 0.;
  Corr(3,3) = 1.;
  Corr(3,4) = -0.92;
  Corr(3,5) = 0.81;
  Corr(4,4) = 1.;
  Corr(4,5) = -0.96;
  Corr(5,5) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(6,6);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.;
  Corr2(0,2) = 0.;
  Corr2(0,3) = 0.;
  Corr2(0,4) = 0.;
  Corr2(0,5) = 0.;
  Corr2(1,1) = 1.;
  Corr2(1,2) = 0.;
  Corr2(1,3) = 0.;
  Corr2(1,4) = 0.;
  Corr2(1,5) = 0.;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.;
  Corr2(2,4) = 0.;
  Corr2(2,5) = 0.;
  Corr2(3,3) = 1.;
  Corr2(3,4) = 0.;
  Corr2(3,5) = 0.;
  Corr2(4,4) = 1.;
  Corr2(4,5) = 0.;
  Corr2(5,5) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID30", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  if(comb == 2){

  // https://arxiv.org/abs/2208.09402

  //First measurements
  //Observables 2:
  CorrData.clear();
  CorrData.push_back(dato(0.132 ,0.011, 0.007)); //A_kpi
  CorrData.push_back(dato( 0.130, 0.012, 0.008)); // A_kpi ^ pipipi0

  //Correlation Matrix (stat):
  Corr.ResizeTo(2,2);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.38;
  Corr(1,1) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(2,2);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.16;
  Corr2(1,1) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("BESIII_Adk", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //Second measurements
  //Observables 2:
  CorrData.clear();
  CorrData.push_back(dato( -0.0562 , 0.0081, 0.0051)); //rD_kpi cos()
  CorrData.push_back(dato( -0.011, 0.012, 0.0076)); // rD_kpi sin()

  //Correlation Matrix (stat):
  Corr.ResizeTo(2,2);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.02;
  Corr(1,1) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("BESIII_rDkpi_polar", CorrelatedGaussianObservables(CorrData, Corr)));

  //https://arxiv.org/pdf/2208.06512.pdf
  //Observables 4:
  CorrData.clear();
  CorrData.push_back(dato(0.0040,0.00040,0.0002));//xcp
  CorrData.push_back(dato(0.00550,0.00120,0.0006));//ycp
  CorrData.push_back(dato(-0.0003,0.0002,0.0000));// dx
  CorrData.push_back(dato(0.0003,0.0003,0.0001));//dy

  //Correlation Matrix (stat):
  Corr.ResizeTo(4,4);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.12;
  Corr(0,2) = -0.02;
  Corr(0,3) = -0.02;
  Corr(1,1) = 1.;
  Corr(1,2) = -0.01;
  Corr(1,3) = -0.06;
  Corr(2,2) = 1.;
  Corr(2,3) = 0.07;
  Corr(3,3) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(4,4);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.08;
  Corr2(0,2) = 0.;
  Corr2(0,3) = -0.01;
  Corr2(1,1) = 1.;
  Corr2(1,2) = -0.02;
  Corr2(1,3) = -0.04;
  Corr2(2,2) = 1.;
  Corr2(2,3) = 0.33;
  Corr2(3,3) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("LHCb_kspp_Au2022", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  }

}
// ---------------------------------------------------------

// ---------------------------------------------------------
void MixingModel::Add_other_meas(){

  vector<dato> CorrData;
  TMatrixDSym Corr, Corr2;

  //-------------------------------------------- Other Observables-------------------------------------------------------------------------

  //14. PDF: charm-kspipi-nocpv (UID13)

  //2 Observables:
  CorrData.clear();
  CorrData.push_back(dato(-0.0086,0.0053,0.0017));//x
  CorrData.push_back(dato(0.0003,0.0046,0.0013));//y

  //Correlation Matrix (stat):
  Corr.ResizeTo(2,2);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.37;
  Corr(1,1) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(2,2);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.;
  Corr2(1,1) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID13", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //19. PDF: charm-k3pi (UID18)

  //Observables 1:
  meas.insert(pair<string,dato>("UID18", dato(0.00005,0.00002,0.))); //k3pi // informations about x,y

  //23. PDF: dkskpiRWS (UID22)

  //Observables 1:
  meas.insert(pair<string,dato>("UID22", dato(0.37,0.003,0.012))); //RD_kskpi

  //24. PDF: dkskpi (UID23)

  //Observables 3:
  CorrData.clear();
  CorrData.push_back(dato(0.356,0.03471,0.)); //RD_kskpi
  CorrData.push_back(dato( - 16.6 * d2r, - 18.4 * d2r,0.));//dD_kskpi // here again the convention dD = - dD
  //CorrData.push_back(dato( (16.6 / 180)* M_PI, (18.4/180.)*M_PI,0.));//dD_kskpi // here again the convention dD = - dD
  CorrData.push_back(dato(0.94,0.12,0.));// kD_kskpi

  //Correlation Matrix (stat):
  Corr.ResizeTo(3,3);
  Corr(0,0) = 1.;
  Corr(0,1) = 0.;
  Corr(0,2) = 0.;
  Corr(1,1) = 1.;
  Corr(1,2) = -0.5;
  Corr(2,2) = 1.;

  //Correlation Matrix (syst)
  Corr2.ResizeTo(3,3);
  Corr2(0,0) = 1.;
  Corr2(0,1) = 0.;
  Corr2(0,2) = 0.;
  Corr2(1,1) = 1.;
  Corr2(1,2) = 0.;
  Corr2(2,2) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID23", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  if(comb == 2){
   // https://arxiv.org/pdf/2208.10098.pdf
   meas.insert(pair<string,dato>("Fpipipipi_BESIII", dato(0.735, 0.015, 0.005))); //F_pipipipi
  }

}
// ---------------------------------------------------------

// ---------------------------------------------------------
void MixingModel::Add_old_meas(){

  tauD = .410; //ps

  epsI = 2.228 * sin(43.52 / 180. * M_PI);

  /*r2d = (rad ? 1. : 180. / M_PI);
  d2r = M_PI / 180.;*/

  vector<dato> CorrData;
  TMatrixDSym Corr, Corr2;

  // UTfit summer 2018 web page A^2 lambda^4 eta
  RCP_in = 0.826 * 0.826 * 0.225 * 0.225 * 0.225 * 0.225 * 0.357;

  // Yellow report current
  RCP_err = (2. * 1.5e-2 + 4. * 0.12e-2 + 3.e-2) * RCP_in;

  // Experimental data from HFLAV Moriond 2019 + LHCb AGamma 1911.01114 + Belle yCP 1912.10912 (irrilevante)

  // 1st Block
  if (combgamma_all){
    meas.insert(pair<string, dato>("ycp", dato(0.835e-2, 0.155e-2))); //without LHCb 2019
  }else{
    meas.insert(pair<string, dato>("ycp", dato(0.719e-2, 0.113e-2))); //HFLAV 2020
  }

  // 2nd Block
  // Because this measurement is just present in UID29 of the LHCb Combo
  if(comb != 2){
    if (!noagamma && !combgamma_all) {
      if (phig12){
        meas.insert(pair<string, dato>("AGamma", dato(1.9e-4, 1.3e-4, 0.4e-4))); //LHCb 2021
      }else{
        meas.insert(pair<string, dato>("AGamma", dato(1.0e-4, 1.1e-4, 0.3e-4))); //LHCb 2021
      }
    }
  }

  //      meas.insert(pair<string, dato>("AGamma", dato(-0.0309e-2, 0.0204e-2))); //HFLAV 2020

  // 3rd Block
  CorrData.clear();

  CorrData.push_back(dato(0.53e-2, .19e-2, .06e-2, .07e-2)); // x_exp_belle_kpp
  CorrData.push_back(dato(0.28e-2, .15e-2, .05e-2, .05e-2)); // y_exp_belle_kpp
  CorrData.push_back(dato(0.91, .16, .05, .06)); // qop_exp_belle_kpp
  CorrData.push_back(dato(-6. * d2r, 11. * d2r, 3. * d2r, 4. * d2r)); // phi_exp_belle_kpp
  //Check HFAG results for all CPV
  //        CorrData.push_back(dato(0.56e-2, .19e-2, .09e-2)); // x_exp_belle_kpp
  //        CorrData.push_back(dato(0.25e-2, .16e-2, .07e-2)); // y_exp_belle_kpp
  //        CorrData.push_back(dato(0.83, .19, .07)); // qop_exp_belle_kpp
  //        CorrData.push_back(dato(-13. / 180. * M_PI, 12. / 180. * M_PI, 4. / 180. * M_PI)); // phi_exp_belle_kpp

  Corr.ResizeTo(4, 4);

  Corr(0, 0) = 1.;
  Corr(0, 1) = .054;
  Corr(0, 2) = -0.074;
  Corr(0, 3) = -0.031;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.034;
  Corr(1, 3) = -0.019;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.044;
  Corr(3, 3) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpp_belle",
  CorrelatedGaussianObservables(CorrData, Corr)));

  //5th Block
  CorrData.clear();

  CorrData.push_back(dato(0.16e-2, 0.23e-2, 0.12e-2, 0.08e-2)); // x_exp_babar
  CorrData.push_back(dato(0.57e-2, 0.20e-2, 0.13e-2, 0.07e-2)); // y_exp_babar

  Corr.ResizeTo(2, 2);

  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.0615;
  Corr(1, 1) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kppkk",
  CorrelatedGaussianObservables(CorrData, Corr)));

  // 6th Block
  meas.insert(pair<string, dato>("RM", dato(0.013e-2, 0.0269e-2))); //HFLAV

  // 8th Block p1
  CorrData.clear();

  CorrData.push_back(dato(2.48e-2, 0.59e-2, 0.39e-2)); // xpp_kpp_babar_exp
  CorrData.push_back(dato(-.07e-2, .65e-2, .5e-2)); // ypp_kpp_babar_exp

  Corr.ResizeTo(2, 2);

  Corr(0, 0) = 1.;
  Corr(0, 1) = -0.69;
  Corr(1, 1) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpp_babar_plus",
  CorrelatedGaussianObservables(CorrData, Corr)));

  // 8th Block p2
  CorrData.clear();

  CorrData.push_back(dato(3.5e-2, 0.78e-2, 0.65e-2)); // xpm_kpp_babar_exp
  CorrData.push_back(dato(-0.82e-2, .68e-2, .41e-2)); // ypm_kpp_babar_exp

  Corr.ResizeTo(2, 2);

  Corr(0, 0) = 1.;
  Corr(0, 1) = -0.66;
  Corr(1, 1) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpp_babar_minus",
  CorrelatedGaussianObservables(CorrData, Corr)));

  // 9th Block
  CorrData.clear();

  CorrData.push_back(dato(0.533e-2, 0.107e-2, .045e-2)); // RD_exp_cleoc
  CorrData.push_back(dato(0.06e-2, 0.23e-2, 0.11e-2)); // xsq_exp_cleoc
  CorrData.push_back(dato(4.2e-2, 2.e-2, 1.e-2)); // y_exp_cleoc
  CorrData.push_back(dato(0.84, 0.20, .06)); // cd_exp_cleoc
  CorrData.push_back(dato(-0.01, 0.41, .04)); // sd_exp_cleoc

  Corr.ResizeTo(5, 5);

  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.;
  Corr(0, 2) = 0.;
  Corr(0, 3) = -0.42;
  Corr(0, 4) = 0.01;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -0.73;
  Corr(1, 3) = 0.39;
  Corr(1, 4) = 0.02;
  Corr(2, 2) = 1.;
  Corr(2, 3) = -0.53;
  Corr(2, 4) = -0.03;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.04;
  Corr(4, 4) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("cleoc",
  CorrelatedGaussianObservables(CorrData, Corr)));

  // 10th Block
  CorrData.clear();

  CorrData.push_back(dato(0.303e-2, .0189e-2)); // RD_exp_babar_kp
  CorrData.push_back(dato(-0.024e-2, .052e-2)); // xp2plus_exp_babar_kp
  CorrData.push_back(dato(0.98e-2, .78e-2)); // ypplus_exp_babar_kp

  Corr.ResizeTo(3, 3);

  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.77;
  Corr(0, 2) = -0.87;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -0.94;
  Corr(2, 2) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpi_babar_plus",
  CorrelatedGaussianObservables(CorrData, Corr)));

  // 11th Block
  CorrData.clear();

  CorrData.push_back(dato(-2.1e-2, 5.4e-2)); // AD_exp_babar_kp
  CorrData.push_back(dato(-0.020e-2, .050e-2)); // xp2minus_exp_babar_kp
  CorrData.push_back(dato(0.96e-2, .75e-2)); // ypminus_exp_babar_kp

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpi_babar_minus",
  CorrelatedGaussianObservables(CorrData, Corr)));

  // 12th Block
  CorrData.clear();

  CorrData.push_back(dato(0.364e-2, .018e-2)); // RD_exp_belle_kp
  CorrData.push_back(dato(0.032e-2, .037e-2)); // xp2plus_exp_belle_kp
  CorrData.push_back(dato(-0.12e-2, .58e-2)); // ypplus_exp_belle_kp

  Corr.ResizeTo(3, 3);

  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.655;
  Corr(0, 2) = -0.834;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -0.909;
  Corr(2, 2) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpi_belle_plus",
  CorrelatedGaussianObservables(CorrData, Corr)));

  //13th Block
  CorrData.clear();

  CorrData.push_back(dato(2.3e-2, 4.7e-2)); // AD_exp_belle_kp
  CorrData.push_back(dato(0.006e-2, .034e-2)); // xp2minus_exp_belle_kp
  CorrData.push_back(dato(0.20e-2, .54e-2)); // ypminus_exp_belle_kp

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpi_belle_minus",
  CorrelatedGaussianObservables(CorrData, Corr)));


  // 14th Block
  CorrData.clear();

  CorrData.push_back(dato(3.51e-3, 0.35e-3)); // RD_exp_cdf_kp
  CorrData.push_back(dato(4.3e-3, 4.3e-3)); // yprime_exp_cdf_kp
  CorrData.push_back(dato(0.08e-3, 0.18e-3)); // xprimesq_exp_cfg_kp

  Corr.ResizeTo(3, 3);

  Corr(0, 0) = 1.;
  Corr(0, 1) = -0.967;
  Corr(0, 2) = 0.9;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -0.975;
  Corr(2, 2) = 1.;

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpi_cdf",
  CorrelatedGaussianObservables(CorrData, Corr)));

  // because they are already present in UID16; UID18; UID14 OF LHCb combination
  if(comb != 2){
    if (!combgamma_all) {

      //15th Block
      CorrData.clear();

      CorrData.push_back(dato(3.454e-3, 0.028e-3, 0.014e-3)); // RD_exp_lhcb_kp
      CorrData.push_back(dato(5.01e-3, 0.48e-3, 0.29e-3)); // ypp_exp_lhcb_kp
      CorrData.push_back(dato(0.61e-4, 0.26e-4, 0.16e-4)); // xpp2_exp_lhcb_kp
      CorrData.push_back(dato(5.54e-3, 0.48e-3, 0.29e-3)); // ypm_exp_lhcb_kp
      CorrData.push_back(dato(0.16e-4, 0.26e-4, 0.16e-4)); // xpm2_exp_lhcb_kp

      Corr.ResizeTo(5, 5);

      Corr(0, 0) = 1.;
      Corr(0, 1) = -0.883;
      Corr(0, 2) = 0.745;
      Corr(0, 3) = -0.883;
      Corr(0, 4) = 0.749;
      Corr(1, 1) = 1.;
      Corr(1, 2) = -0.944;
      Corr(1, 3) = 0.758;
      Corr(1, 4) = -0.644;
      Corr(2, 2) = 1.;
      Corr(2, 3) = -0.642;
      Corr(2, 4) = 0.545;
      Corr(3, 3) = 1.;
      Corr(3, 4) = -0.946;
      Corr(4, 4) = 1.;

      corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpi_lhcb",
      CorrelatedGaussianObservables(CorrData, Corr)));

      // 7th Block
      meas.insert(pair<string, dato>("RMKppp", dato(0.0096e-2, 0.0036e-2))); //LHCb Kppp multiplied by 2

      //4th Block
      CorrData.clear();

      CorrData.push_back(dato(3.97e-3, .46e-3, .29e-3)); // XCP_LHCb_kspp
      CorrData.push_back(dato(4.59e-3, 1.20e-3, .85e-3)); // YCP_LHCb_kspp
      CorrData.push_back(dato(-.27e-3, .18e-3, .01e-3)); // DX_LHCb_kspp
      CorrData.push_back(dato(.20e-3, .36e-3, .13e-3)); // DY_LHCb_kspp

      Corr.ResizeTo(4, 4);

      Corr(0, 0) = 1.;
      Corr(0, 1) = .11;
      Corr(0, 2) = -.02;
      Corr(0, 3) = -.01;
      Corr(1, 1) = 1.;
      Corr(1, 2) = -.01;
      Corr(1, 3) = -.05;
      Corr(2, 2) = 1.;
      Corr(2, 3) = .08;
      Corr(3, 3) = 1.;

      Corr2.ResizeTo(4, 4);

      Corr2(0, 0) = 1.;
      Corr2(0, 1) = .13;
      Corr2(0, 2) = .01;
      Corr2(0, 3) = 0.01;
      Corr2(1, 1) = 1.;
      Corr2(1, 2) = -.02;
      Corr2(1, 3) = .01;
      Corr2(2, 2) = 1.;
      Corr2(2, 3) = .31;
      Corr2(3, 3) = 1.;

      corrmeas.insert(pair<string, CorrelatedGaussianObservables>("LHCb_kspp",
      CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

    } else {

      //LHCb combo results
      CorrData.clear();

      CorrData.push_back(dato(4.00e-3, .52e-3)); // X_LHCb_combo
      CorrData.push_back(dato(6.31e-3, .32e-3)); // Y_LHCb_combo
      CorrData.push_back(dato(.997, .016)); // qop_LHCb_combo
      CorrData.push_back(dato(-2.4, 1.2)); // phi_LHCb_combo (degrees)
      CorrData.push_back(dato(.05867, .00015)); // rDkp_LHCb_combo
      CorrData.push_back(dato(10., 4.2)); // deltaDKpi_LHCb_combo (degrees)


      Corr.ResizeTo(6, 6);

      Corr(0, 0) = 1.;
      Corr(0, 1) = .013;
      Corr(0, 2) = -.129;
      Corr(0, 3) = .083;
      Corr(0, 4) = .282;
      Corr(0, 5) = .029;
      Corr(1, 1) = 1.;
      Corr(1, 2) = -.018;
      Corr(1, 3) = .04;
      Corr(1, 4) = -.095;
      Corr(1, 5) = .891;
      Corr(2, 2) = 1.;
      Corr(2, 3) = .554;
      Corr(2, 4) = -.07;
      Corr(2, 5) = -.02;
      Corr(3, 3) = 1.;
      Corr(3, 4) = .015;
      Corr(3, 5) = .055;
      Corr(4, 4) = 1.;
      Corr(4, 5) = .295;
      Corr(5, 5) = 1.;


      corrmeas.insert(pair<string, CorrelatedGaussianObservables>("LHCb_combo",
      CorrelatedGaussianObservables(CorrData, Corr)));
    }


    if (combgamma_delta) {

      //LHCb D -> Kpi results
      CorrData.clear();

      CorrData.push_back(dato(.05867, .00015)); // rDkp_LHCb_combo
      CorrData.push_back(dato(10., 4.2)); // deltaDKpi_LHCb_combo (degrees)

      Corr.ResizeTo(2, 2);

      Corr(0, 0) = 1.;
      Corr(0, 1) = 0.295;
      Corr(1, 1) = 1.;

      corrmeas.insert(pair<string, CorrelatedGaussianObservables>("LHCb_combo_delta",
      CorrelatedGaussianObservables(CorrData, Corr)));

    }
  }


}
// ---------------------------------------------------------


// ---------------------------------------------------------
void MixingModel::DefineParameters()
{
  //------------------------------- Adding model parameters for all the combinations ------------------------------------------------------------------------

  if(comb == 0){

    //------------------------ Time Integrated B decay chain parameters -----------------------------------------------
    //Global parameters
    AddParameter("g", 0.7  ,  1.5, "#gamma"); //All combined interval [1, 1.6]
    AddParameter("x12", 0.5e-3, 8.e-3,"x_{12}"); // ACI [2e-3, 6e-3]
    AddParameter("y12", 3e-3, 9e-3, "y_{12}"); // ACI [4e-3, 1e-2]

    //1. PDF: glwads-dh-hh-dmix (UID0)
    //Parametri:9
    AddParameter("r_dk", 0.07, 0.13, "r_{DK}"); // LHCb 0.0984 [0.0958,0.1011]
    AddParameter("r_dpi",0.,0.012,"r_{D#pi}"); // LHCb 0.00480 [0.00424, 0.00550]
    AddParameter("rD_kpi", 0.055, 0.065, "r^{D}_{K#pi}"); // LHCb 0.05867 [0.05852, 0.05882]
    AddParameter("d_dk", 1.5,  3., "#delta_{DK}"); // LHCb  2.23 [2.153, 2.2968] //ACI [2.1, 2.6]
    AddParameter("d_dpi", 3., 2*M_PI, "#delta_{D#pi}"); //LHCb [4.76, 5.30] //ACI [-1.5, 1.]
    AddParameter("dD_kpi", 2.2, 3.8, "#delta^{D}_{K#pi}"); //LHCb [120, 230]G //ACI[4.,5.]

    //2. PDF: glwads-dh-h3pi-dmix (UID1)
    //Parametri:11
    AddParameter("rD_k3pi", 0.048, 0.065, "r^{D}_{K3#pi}"); // 0.0556 [0.05502,0.05619]
    AddParameter("dD_k3pi", 1.8, 5.5, "#delta^{D}_{K3#pi}"); // 3.595 [3.87,3.33]
    AddParameter("kD_k3pi", 0., 1., "#kappa^{D}_{K3#pi}"); // 0.480
    AddParameter("F_pipipipi", 0.5, 1., "F_{#pi#pi#pi#pi}"); // 0.737

    //3. PDF: glwads-dh-hhpi0-dmix (UID2)
    //Parametri:12
    AddParameter("rD_kpipi0", 0.03, 0.06, "r^{D}_{K#pi#pi^{0}}"); //0.044
    AddParameter("dD_kpipi0", 1.5, 4.5, "#delta^{D}_{K#pi#pi^{0}}"); // [3.28,2.86]
    AddParameter("kD_kpipi0", 0.4, 1., "#kappa^{D}_{K#pi#pi^{0}}"); // 0.79
    AddParameter("F_pipipi0", 0.6, 1., "F_{#pi#pi#pi^{0}}"); // 0.973
    AddParameter("F_kkpi0", 0.2, 1., "F_{KK#pi^0}"); // 0.732

    //4. PDF: ggsz-dh (UID3)
    //Parametri: 5

    //5. PDF: glwads-dkdpi-kskpi-dmix (UID4)
    //Parametri: 11
    AddParameter("rD_kskpi", 0.55, 0.7, "r^{D}_{K^{0}_S K #pi}"); //0.6150
    AddParameter("dD_kskpi", -2., 2., "#delta^{D}_{K^{0}_S K #pi}"); // [6.23, 5.75]
    AddParameter("kD_kskpi", 0.3, 1., "#kappa_{K^{0}_S K #pi}"); // 0.82
    AddParameter("RBRdkdpi", 0.05, 0.1, "#frac{B^{-} ---> D^0 K^{-}}{B^{-} ---> D^0 #pi^{-}}"); //

    //6. PDF: glwads-dsth-hh-dmix (UID5)
    //Parametri: 9
    AddParameter("r_dstk", 0., 0.3, "r_{D^{*}K}"); // 0.099
    AddParameter("d_dstk", -4., 1., "#delta_{D^{*}K}"); // 0.873
    AddParameter("r_dstpi", 0., 0.06, "r_{D^{*} #pi}"); //0.0095 ACI [0., 0.03]
    AddParameter("d_dstpi", -1., 2*M_PI -1., "#delta_{D^{*} #pi}"); // 2.426 ACI [0., 4.5]

    //7. PDF: glwads-dkst-hh-h3pi-dmix-newvars (UID6)
    //Parametri:12
    AddParameter("r_dkst", 0., 0.3, "r_{DK^{*}}"); // 0.106
    AddParameter("d_dkst", -1.5, 4.5, "#delta_{DK^{*}}"); // 0.61 //ACI [0.,4.]
    AddParameter("k_dkst", 0.5, 1., "#kappa_{DK^{*}}"); //

    //8. PDF: glwads-dkstz-hh-h3pi-dmix (UID7)
    //Parametri:12
    AddParameter("r_dkstz", 0., 0.5, "r_{D K^{0*}}");  //0.25
    AddParameter("d_dkstz", 2., 5.5, "#delta_{D K^{0*}}"); //2.884
    AddParameter("k_dkstz", 0.75, 1., "#kappa_{D K^{0*}}"); //

    //9. PDF: ggsz_dkstz_lhcb_md (UID8)
    //Parametri:4

    //10. PDF: glwads-dhpipi-hh-dmix (UID9)
    //Parametri:11
    AddParameter("r_dkpipi", 0., 0.2, "r_{DK #pi #pi}"); //0.079
    AddParameter("d_dkpipi", -M_PI -1.5, M_PI -1.5, "#delta_{DK #pi #pi}"); // //ACI [1., 2*M_PI]
    AddParameter("k_dkpipi", 0., 1., "#kappa_{DK #pi #pi}"); //
    AddParameter("r_dpipipi",0., 0.2, "r_{D #pi #pi #pi}"); //0.067
    AddParameter("d_dpipipi", -M_PI - 1.5, M_PI - 1.5, "#delta_{D #pi #pi #pi}"); //
    AddParameter("k_dpipipi", 0., 1., "#kappa_{D #pi #pi #pi}"); //

    //---------------------------------- Parametri Time Dependent B decay  ---------------------------------------------------

    //11. PDF: dsk (UID10)
    //Parametri:4
    AddParameter("l_dsk", 0., 1.2, "#lambda_{D^{-}_{s} K}"); // [0.218,0.408]
    AddParameter("d_dsk", -2.5, 2.5, "#delta_{D^{-}_{s} K}"); //[4.5, 2*M_PI ]
    AddParameter("phis", -0.3, 0.2, "#phi_s");

    //12. PDF: dskpipi (UID11)
    //Parametri:5
    AddParameter("l_dskpipi", 0., 1.2, "#lambda_{D^{-}_{s} K #pi #pi}");
    AddParameter("d_dskpipi", -2., 1.5, "#delta_{D_s K #pi #pi}");
    AddParameter("k_dskpipi", 0., 1., "#kappa_{D_s K #pi #pi}");

    //13. PDF: dmpi (UID12)
    //Parametri:4
    AddParameter("l_dmpi", 0., 0.4, "#lambda_{D^{-} #pi}");
    AddParameter("d_dmpi", -M_PI, M_PI, "#delta_{D^{-} #pi}");
    AddParameter("beta", 0.2, 0.5, "#beta");

    //---------------------------------- Parametri Time Dependent D decay  ---------------------------------------------------

    //14. PDF: charm-kspipi-nocpv (UID13)
    //Parametri:2

    //15. PDF: charm-kspipi (UID14)
    //Parametri:4
    AddParameter("PhiM12",-0.2,0.3,"#phi_{M}");
    AddParameter("PhiG12", -0.2, 0.3, "#phi_{#Gamma}");


    //16. PDF: charm-kspipi (UID15)
    //Parametri:4

    //17. PDF: charm-kpi (UID16)
    //Parametri:7
    AddParameter("AD", -4e-2, 3e-2, "A_D");

    //18. PDF: charm-deltaacp-diff (UID17)
    //Parametri: 6
    AddParameter("Delta_totau", 0.1, 0.13, "<#Delta t>"); //0.115
    AddParameter("DeltaAcp", -3.5e-3, 0., "#Delta a_{CP}"); //-0.00152

    SetPriorConstantAll();
  }
  else if(comb == 1){

    AddParameter("x12", 0., 2.e-2, "x_{12}");
    AddParameter("y12", 0., 5.e-2, "y_{12}");
    AddParameter("PhiM12", -0.2, 0.2, "#phi_{M}");
    AddParameter("dD_kpi", 0., 2*M_PI, "#delta^D_{K #pi}"); //dD_kpi
    AddParameter("dD_kpipi0", 0., 2*M_PI, "#delta^D_{K #pi #pi^0}"); // dD_kpipi0
    AddParameter("rD_kpi", 0., 0.1, "r^D_{K #pi}"); // rD_kpi
    if (phig12)
    AddParameter("PhiG12", -0.2, 0.2, "#phi_{#Gamma}");
    else
    AddParameter("PhiG12", 0., 0.);

    if (ckmcorr)
    AddParameter("RCP", RCP_in - 5. * RCP_err, RCP_in + 5. * RCP_err);

    SetPriorConstantAll();
    if (ckmcorr)
    GetParameter("RCP").SetPrior(new BCGaussianPrior(RCP_in, RCP_err));

  }
  else if(comb == 2){
    //---------------------------  Time Integrated B decay chain Parameters  -----------------------------------------------
    //Global parameters
    AddParameter("g", 0.7  ,  1.5, "#gamma"); //All combined interval [1, 1.6]
    AddParameter("x12", 0.5e-3, 8.e-3,"x_{12}"); // ACI [2e-3, 6e-3]
    AddParameter("y12", 3e-3, 9e-3, "y_{12}"); // ACI [4e-3, 1e-2]

    //1. PDF: glwads-dh-hh-dmix (UID0)
    //Parametri:9
    AddParameter("r_dk", 0.07, 0.13, "r_{DK}"); // LHCb 0.0984 [0.0958,0.1011]
    AddParameter("r_dpi",0.,0.012,"r_{D#pi}"); // LHCb 0.00480 [0.00424, 0.00550]
    AddParameter("rD_kpi", 0.055, 0.065, "r^{D}_{K#pi}"); // LHCb 0.05867 [0.05852, 0.05882]
    AddParameter("d_dk", 1.5,  3., "#delta_{DK}"); // LHCb  2.23 [2.153, 2.2968] //ACI [2.1, 2.6]
    AddParameter("d_dpi", 3., 2*M_PI, "#delta_{D#pi}"); //LHCb [4.76, 5.30] //ACI [-1.5, 1.]
    AddParameter("dD_kpi", 2.2, 3.8, "#delta^{D}_{K#pi}"); //LHCb [120, 230]G //ACI[4.,5.]

    //2. PDF: glwads-dh-h3pi-dmix (UID1)
    //Parametri:11
    AddParameter("rD_k3pi", 0.048, 0.065, "r^{D}_{K3#pi}"); // 0.0556 [0.05502,0.05619]
    AddParameter("dD_k3pi", 1.8, 5.5, "#delta^{D}_{K3#pi}"); // 3.595 [3.87,3.33]
    AddParameter("kD_k3pi", 0., 1., "#kappa^{D}_{K3#pi}"); // 0.480
    AddParameter("F_pipipipi", 0.5, 1., "F_{#pi#pi#pi#pi}"); // 0.737

    //3. PDF: glwads-dh-hhpi0-dmix (UID2)
    //Parametri:12
    AddParameter("rD_kpipi0", 0.03, 0.06, "r^{D}_{K#pi#pi^{0}}"); //0.044
    AddParameter("dD_kpipi0", 1.5, 4.5, "#delta^{D}_{K#pi#pi^{0}}"); // [3.28,2.86]
    AddParameter("kD_kpipi0", 0.4, 1., "#kappa^{D}_{K#pi#pi^{0}}"); // 0.79
    AddParameter("F_pipipi0", 0.6, 1., "F_{#pi#pi#pi^{0}}"); // 0.973
    AddParameter("F_kkpi0", 0.2, 1., "F_{KK#pi^0}"); // 0.732

    //4. PDF: ggsz-dh (UID3)
    //Parametri: 5

    //5. PDF: glwads-dkdpi-kskpi-dmix (UID4)
    //Parametri: 11
    AddParameter("rD_kskpi", 0.55, 0.7, "r^{D}_{K^{0}_S K #pi}"); //0.6150
    AddParameter("dD_kskpi", -2., 2., "#delta^{D}_{K^{0}_S K #pi}"); // [6.23, 5.75]
    AddParameter("kD_kskpi", 0.3, 1., "#kappa_{K^{0}_S K #pi}"); // 0.82
    AddParameter("RBRdkdpi", 0.05, 0.1, "#frac{B^{-} ---> D^0 K^{-}}{B^{-} ---> D^0 #pi^{-}}"); //

    //6. PDF: glwads-dsth-hh-dmix (UID5)
    //Parametri: 9
    AddParameter("r_dstk", 0., 0.3, "r_{D^{*}K}"); // 0.099
    AddParameter("d_dstk", -4., 1., "#delta_{D^{*}K}"); // 0.873
    AddParameter("r_dstpi", 0., 0.06, "r_{D^{*} #pi}"); //0.0095 ACI [0., 0.03]
    AddParameter("d_dstpi", -1., 2*M_PI -1., "#delta_{D^{*} #pi}"); // 2.426 ACI [0., 4.5]

    //7. PDF: glwads-dkst-hh-h3pi-dmix-newvars (UID6)
    //Parametri:12
    AddParameter("r_dkst", 0., 0.3, "r_{DK^{*}}"); // 0.106
    AddParameter("d_dkst", -1.5, 4.5, "#delta_{DK^{*}}"); // 0.61 //ACI [0.,4.]
    AddParameter("k_dkst", 0.5, 1., "#kappa_{DK^{*}}"); //

    //8. PDF: glwads-dkstz-hh-h3pi-dmix (UID7)
    //Parametri:12
    AddParameter("r_dkstz", 0., 0.5, "r_{D K^{0*}}");  //0.25
    AddParameter("d_dkstz", 2., 5.5, "#delta_{D K^{0*}}"); //2.884
    AddParameter("k_dkstz", 0.75, 1., "#kappa_{D K^{0*}}"); //

    //9. PDF: ggsz_dkstz_lhcb_md (UID8)
    //Parametri:4

    //10. PDF: glwads-dhpipi-hh-dmix (UID9)
    //Parametri:11
    AddParameter("r_dkpipi", 0., 0.2, "r_{DK #pi #pi}"); //0.079
    AddParameter("d_dkpipi", -M_PI -1.5 , M_PI - 1.5, "#delta_{DK #pi #pi}"); // //ACI [1., 2*M_PI]
    AddParameter("k_dkpipi", 0., 1., "#kappa_{DK #pi #pi}"); //
    AddParameter("r_dpipipi",0., 0.2, "r_{D #pi #pi #pi}"); //0.067
    AddParameter("d_dpipipi", -M_PI - 1.5, M_PI -1.5, "#delta_{D #pi #pi #pi}"); //
    AddParameter("k_dpipipi", 0., 1., "#kappa_{D #pi #pi #pi}"); //

    //---------------------------------- Parametri Time Dependent B decay  ---------------------------------------------------

    //11. PDF: dsk (UID10)
    //Parametri:4
    AddParameter("l_dsk", 0., 1.2, "#lambda_{D^{-}_{s} K}"); // [0.218,0.408]
    AddParameter("d_dsk", -2.5, 2.5, "#delta_{D^{-}_{s} K}"); //[4.5, 2*M_PI ]
    AddParameter("phis", -0.3, 0.2, "#phi_s");

    //12. PDF: dskpipi (UID11)
    //Parametri:5
    AddParameter("l_dskpipi", 0., 1.2, "#lambda_{D^{-}_{s} K #pi #pi}");
    AddParameter("d_dskpipi", -2., 1.5, "#delta_{D_s K #pi #pi}");
    AddParameter("k_dskpipi", 0., 1., "#kappa_{D_s K #pi #pi}");

    //13. PDF: dmpi (UID12)
    //Parametri:4
    AddParameter("l_dmpi", 0., 0.4, "#lambda_{D^{-} #pi}");
    AddParameter("d_dmpi", -M_PI, M_PI, "#delta_{D^{-} #pi}");
    AddParameter("beta", 0.2, 0.5, "#beta");

    //---------------------------------- Parametri Time Dependent D decay  ---------------------------------------------------

    //14. PDF: charm-kspipi-nocpv (UID13)
    //Parametri:2

    //15. PDF: charm-kspipi (UID14)
    //Parametri:4
    AddParameter("PhiM12",-0.2,0.3,"#phi_{M}");
    if (phig12){
      AddParameter("PhiG12", -0.2, 0.2, "#phi_{#Gamma}");
    }
    else{
      AddParameter("phiG12", 0., 0.);
    }

    //16. PDF: charm-kspipi (UID15)
    //Parametri:4

    //17. PDF: charm-kpi (UID16)
    //Parametri:7
    AddParameter("AD", -4e-2, 3e-2, "A_D");

    //18. PDF: charm-deltaacp-diff (UID17)
    //Parametri: 6
    AddParameter("Delta_totau", 0.1, 0.13, "<#Delta t>"); //0.115
    AddParameter("DeltaAcp", -3.5e-3, 0., "#Delta a_{CP}"); //-0.00152

    //19. PDF: charm-k3pi (UID18)
    //Parametri: 2

    //20. PDF: dk3pi_dkpipi0_constraints (UID19)
    //Parametri: 6 kD_k3pi, dD_k3pi, kD_kpipi0, dD_kpipi0, rD_k3pi, rD_kpipi0

    //21. PDF: d4pi_dmixing_cleo (UID20)
    //Parametri:1 F_pipipipi

    //22. PDF: CleoDhhpi0Dilution (UID21)
    //Parametri:2  F_pipipi0, F_kkpi0

    //23. PDF: dkskpiRWS (UID22)
    //Parametri:5  rD_kskpi, dD_kskpi, kD_kskpi, xD, yD

    //24. PDF: dkskpi (UID23)
    //Parametri 5:  rD_kskpi, dD_kskpi, kD_kskpi, xD, yD

    //25. PDF: dkstcoherence (UID24)
    //Parametri 1: k_dkst

    //26. PDF: dkstzcoherence (UID25)
    //Parametri 1: k_dkstz

    //27. PDF: phis (UID26)
    //Parametri 1: phis

    //28. PDF: beta (UID27)
    //Parametri 1: beta

    //29. PDF: Charm_ycp_minus_ycp_rs_lhcb-r1 (UID28)
    //Parametri 6: rD_kpi, dD_kpi, xD, yD, qopD, phiD

    //30. PDF: charm-dy-rs (UID29)
    //Parametri 4:  xD, yD, qopD, phiD

    //31. PDF: charm-kpi (UID30
    //Parametri 7: rD_kpi, dD_kpi, xD, yD, qopD, phiD, AD

    //#Parameters = 52

    if (ckmcorr){
      AddParameter("RCP", RCP_in - 5. * RCP_err, RCP_in + 5. * RCP_err);
    }
    SetPriorConstantAll();
    if (ckmcorr){
      GetParameter("RCP").SetPrior(new BCGaussianPrior(RCP_in, RCP_err));
    }
  }

}
// ---------------------------------------------------------

// ---------------------------------------------------------
void MixingModel::DefineHistograms(){

  if(comb == 0){
    //1D Histograms of the parameters

    //Time integrated decay chain
    histos.createH1D("g", 200, 1., -1.);
    histos.createH1D("x12", 200, 1., -1.);
    histos.createH1D("y12", 200, 1., -1.);
    histos.createH1D("x", 200, 1., -1.);
    histos.createH1D("y", 200, 1., -1.);
    histos.createH1D("r_dk", 200, 1., -1.);
    histos.createH1D("r_dpi", 200, 1., -1.);
    histos.createH1D("rD_kpi", 200, 1., -1.);
    histos.createH1D("d_dk", 200, 1., -1.);
    histos.createH1D("d_dpi", 200, 1., -1.);
    histos.createH1D("dD_kpi", 200, 1., -1.);
    histos.createH1D("dD_kpi_LHCb", 200, 1., -1.);
    histos.createH1D("rD_k3pi", 200, 1., -1.);
    histos.createH1D("dD_k3pi", 200, 1., -1.);
    histos.createH1D("dD_k3pi_LHCb", 200, 1., -1.);
    histos.createH1D("kD_k3pi", 200, 1., -1.);
    histos.createH1D("F_pipipipi", 200, 1., -1.);
    histos.createH1D("rD_kpipi0", 200, 1., -1.);
    histos.createH1D("dD_kpipi0", 200, 1., -1.);
    histos.createH1D("dD_kpipi0_LHCb", 200, 1., -1.);
    histos.createH1D("kD_kpipi0", 200, 1., -1.);
    histos.createH1D("F_pipipi0", 200, 1., -1.);
    histos.createH1D("F_kkpi0", 200, 1., -1.);
    histos.createH1D("rD_kskpi", 200, 1., -1.);
    histos.createH1D("dD_kskpi", 200, 1., -1.);
    histos.createH1D("dD_kskpi_LHCb", 200, 1., -1.);
    histos.createH1D("kD_kskpi", 200, 1., -1.);
    histos.createH1D("RBRdkdpi", 200, 1., -1.);
    histos.createH1D("r_dstk", 200, 1., -1.);
    histos.createH1D("d_dstk", 200, 1., -1.);
    histos.createH1D("r_dstpi", 200, 1., -1.);
    histos.createH1D("d_dstpi", 200, 1., -1.);
    histos.createH1D("r_dkst", 200, 1., -1.);
    histos.createH1D("d_dkst", 200, 1., -1.);
    histos.createH1D("k_dkst", 200, 1., -1.);
    histos.createH1D("r_dkstz", 200, 1., -1.);
    histos.createH1D("d_dkstz", 200, 1., -1.);
    histos.createH1D("k_dkstz", 200, 1., -1.);
    histos.createH1D("r_dkpipi", 200, 1., -1.);
    histos.createH1D("d_dkpipi", 200, 1., -1.);
    histos.createH1D("k_dkpipi", 200, 1., -1.);
    histos.createH1D("r_dpipipi", 200, 1., -1.);
    histos.createH1D("d_dpipipi", 200, 1., -1.);
    histos.createH1D("k_dpipipi", 200, 1., -1.);

    //Time dependent B decay and mixing
    histos.createH1D("l_dsk", 200, 1., -1.);
    histos.createH1D("d_dsk", 200, 1., -1.);
    histos.createH1D("phis", 200, 1., -1.);
    histos.createH1D("beta_s", 200, 1., -1.);

    histos.createH1D("l_dskpipi", 200, 1., -1.);
    histos.createH1D("d_dskpipi", 200, 1., -1.);
    histos.createH1D("k_dskpipi", 200, 1., -1.);

    histos.createH1D("l_dmpi", 200, 1., -1.);
    histos.createH1D("d_dmpi", 200, 1., -1.);
    histos.createH1D("beta", 200, 1., -1.);

    //Time dependent D decay and mixing
    histos.createH1D("PhiM12", 200, 1., -1.);
    histos.createH1D("PhiG12", 200, 1., -1.);
    histos.createH1D("AD", 200, 1., -1.);
    histos.createH1D("Delta_totau", 200, 1., -1.);
    histos.createH1D("DeltaAcp", 200, 1., -1.);
    histos.createH1D("phi", 200, 1., -1.);
    histos.createH1D("qopm1", 200, 1., -1.);
    histos.createH1D("qop", 200, 1., -1.);
    histos.createH1D("phi12", 200, 1., -1.);
    histos.createH1D("delta", 200, 1., -1.);

    // 2D Histograms

    //Time integrated decay chain
    histos.createH2D("r_dk", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dk", "r_dk", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpi", "r_dpi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpi", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dk", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpi", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkstz", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkstz", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpi", "rD_kpi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y12", "x12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "x", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y12", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x12", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "x", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "y", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "qopm1", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x", "qopm1", 200, 1., -1., 200, 1., -1);

    //Time dependent D mixing and Decay
    histos.createH2D("y12", "qopm1", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x12", "qopm1", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "qopm1", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "x12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "y12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "x12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "y12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "x", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "y", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiG12", "x12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiG12", "y12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiG12", "x", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiG12", "y", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "PhiG12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpi", "dD_kpi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_k3pi", "dD_k3pi", 200, 1., -1., 200, 1., -1);

    histos.createH2D("g", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x12", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y12", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dk", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dk", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_k3pi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_k3pi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_k3pi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_pipipipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_pipipi0", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpipi0", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpipi0", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_kpipi0", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_kkpi0", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kskpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kskpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_kskpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("RBRdkdpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dstk", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dstk", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dstpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dstpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkst", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkst", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkst", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkstz", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkstz", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkstz", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkpipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkpipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkpipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpipipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpipipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dpipipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dsk", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dsk", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phis", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("beta_s", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dskpipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dskpipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dskpipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dmpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dmpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("beta", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiG12", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("AD", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("Delta_totau", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("DeltaAcp", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi12", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("delta", "qop", 200, 1., -1., 200, 1., -1);

    histos.createH2D("g", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x12", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y12", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dk", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dk", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_k3pi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_k3pi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_k3pi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_pipipipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_pipipi0", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpipi0", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpipi0", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_kpipi0", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_kkpi0", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kskpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kskpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_kskpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("RBRdkdpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dstk", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dstk", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dstpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dstpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkst", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkst", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkst", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkstz", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkstz", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkstz", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkpipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkpipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkpipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpipipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpipipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dpipipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dsk", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dsk", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phis", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("beta_s", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dskpipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dskpipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dskpipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dmpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dmpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("beta", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiG12", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("AD", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("Delta_totau", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("DeltaAcp", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi12", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("delta", "phi", 200, 1., -1., 200, 1., -1);

    histos.createH2D("g", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x12", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y12", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dk", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dk", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_k3pi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_k3pi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_k3pi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_pipipipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_pipipi0", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpipi0", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpipi0", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_kpipi0", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_kkpi0", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kskpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kskpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_kskpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("RBRdkdpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dstk", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dstk", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dstpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dstpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkst", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkst", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkst", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkstz", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkstz", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkstz", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkpipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkpipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkpipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpipipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpipipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dpipipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dsk", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dsk", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phis", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("beta_s", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dskpipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dskpipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dskpipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dmpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dmpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("beta", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiG12", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("AD", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("Delta_totau", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("DeltaAcp", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("delta", "phi12", 200, 1., -1., 200, 1., -1);


  }
  else if(comb == 1){

    if (phig12) {
      histos.createH2D("phi", "phiG12", 200, 1., -1., 200, 1., -1);
      histos.createH2D("phiM12", "phiG12", 200, 1., -1., 200, 1., -1);
      histos.createH2D("phiG12", "y12", 200, 1., -1., 200, 1., -1);
      histos.createH2D("phipphig12", "phimphig12", 200, 1., -1., 200, 1., -1);
      histos.createH1D("phipphig12", 200, 1., -1.);
      histos.createH1D("phimphig12", 200, 1., -1.);
      histos.createH1D("phiG12", 200, 1., -1.);
    }

    histos.createH2D("phi", "qopm1", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "x", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "y", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phiM12", "x12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "x", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y12", "x12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "qopm1", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x", "qopm1", 200, 1., -1., 200, 1., -1);
    histos.createH1D("qopm1", 200, 1., -1.);
    histos.createH1D("delta", 200, 1., -1.);
    histos.createH1D("x", 200, 1., -1.);
    histos.createH1D("y", 200, 1., -1.);
    histos.createH1D("x12", 200, 1., -1.);
    histos.createH1D("y12", 200, 1., -1.);
    histos.createH1D("phi", 200, 1., -1.);
    histos.createH1D("phiM12", 200, 1., -1.);
    histos.createH1D("phi12", 200, 1., -1.);
    histos.createH1D("RD", 200, 1., -1.);
    histos.createH1D("Rm", 200, 1., -1.);
    histos.createH1D("Am", 200, 1., -1.);
    histos.createH1D("AGamma", 200, 1., -1.);
    histos.createH1D("ycp", 200, 1., -1.);
    histos.createH1D("xcp", 200, 1., -1.);
    histos.createH1D("dx", 200, 1., -1.);
    histos.createH1D("ypp", 200, 1., -1.);
    histos.createH1D("ypm", 200, 1., -1.);
    histos.createH1D("xpp", 200, 1., -1.);
    histos.createH1D("xpm", 200, 1., -1.);
    histos.createH1D("xpp2", 200, 1., -1.);
    histos.createH1D("xpm2", 200, 1., -1.);
    histos.createH1D("ypp_kpp", 200, 1., -1.);
    histos.createH1D("ypm_kpp", 200, 1., -1.);
    histos.createH1D("xpp_kpp", 200, 1., -1.);
    histos.createH1D("xpm_kpp", 200, 1., -1.);
    histos.createH1D("delta_kp", 200, 1., -1.);
    histos.createH1D("delta_kpp", 200, 1., -1.);
    histos.createH1D("dD_kpi", 200, 1., -1.);
    histos.createH1D("dD_kpipi0", 200, 1., -1.);
    if (ckmcorr) {
      histos.createH1D("phiKsLHCb", 200, 1., -1.);
      histos.createH1D("RCP", 200, 1., -1.);
    }


  }
  else if(comb == 2){

    if (phig12) {
        histos.createH2D("phi", "PhiG12", 200, 1., -1., 200, 1., -1);
        histos.createH2D("PhiM12", "PhiG12", 200, 1., -1., 200, 1., -1);
        histos.createH2D("PhiG12", "y12", 200, 1., -1., 200, 1., -1);
        histos.createH2D("PhiG12", "y", 200, 1., -1., 200, 1., -1);
        histos.createH2D("PhiG12", "x", 200, 1., -1., 200, 1., -1);
        histos.createH2D("PhiG12", "x12", 200, 1., -1., 200, 1., -1);
        histos.createH2D("phipphig12", "phimphig12", 200, 1., -1., 200, 1., -1);
        histos.createH1D("phipphig12", 200, 1., -1.);
        histos.createH1D("phimphig12", 200, 1., -1.);
        histos.createH1D("PhiG12", 200, 1., -1.);
        histos.createH2D("PhiG12", "qop", 200, 1., -1., 200, 1., -1);
        histos.createH2D("PhiG12", "phi", 200, 1., -1., 200, 1., -1);
        histos.createH2D("PhiG12", "phi12", 200, 1., -1., 200, 1., -1);

    }

    //1D Histograms of the parameters

    //Time integrated decay chain
    histos.createH1D("g", 200, 1., -1.);
    histos.createH1D("x12", 200, 1., -1.);
    histos.createH1D("y12", 200, 1., -1.);
    histos.createH1D("x", 200, 1., -1.);
    histos.createH1D("y", 200, 1., -1.);
    histos.createH1D("r_dk", 200, 1., -1.);
    histos.createH1D("r_dpi", 200, 1., -1.);
    histos.createH1D("rD_kpi", 200, 1., -1.);
    histos.createH1D("d_dk", 200, 1., -1.);
    histos.createH1D("d_dpi", 200, 1., -1.);
    histos.createH1D("dD_kpi", 200, 1., -1.);
    histos.createH1D("dD_kpi_LHCb", 200, 1., -1.);
    histos.createH1D("rD_k3pi", 200, 1., -1.);
    histos.createH1D("dD_k3pi", 200, 1., -1.);
    histos.createH1D("dD_k3pi_LHCb", 200, 1., -1.);
    histos.createH1D("kD_k3pi", 200, 1., -1.);
    histos.createH1D("F_pipipipi", 200, 1., -1.);
    histos.createH1D("rD_kpipi0", 200, 1., -1.);
    histos.createH1D("dD_kpipi0", 200, 1., -1.);
    histos.createH1D("dD_kpipi0_LHCb", 200, 1., -1.);
    histos.createH1D("kD_kpipi0", 200, 1., -1.);
    histos.createH1D("F_pipipi0", 200, 1., -1.);
    histos.createH1D("F_kkpi0", 200, 1., -1.);
    histos.createH1D("rD_kskpi", 200, 1., -1.);
    histos.createH1D("dD_kskpi", 200, 1., -1.);
    histos.createH1D("dD_kskpi_LHCb", 200, 1., -1.);
    histos.createH1D("kD_kskpi", 200, 1., -1.);
    histos.createH1D("RBRdkdpi", 200, 1., -1.);
    histos.createH1D("r_dstk", 200, 1., -1.);
    histos.createH1D("d_dstk", 200, 1., -1.);
    histos.createH1D("r_dstpi", 200, 1., -1.);
    histos.createH1D("d_dstpi", 200, 1., -1.);
    histos.createH1D("r_dkst", 200, 1., -1.);
    histos.createH1D("d_dkst", 200, 1., -1.);
    histos.createH1D("k_dkst", 200, 1., -1.);
    histos.createH1D("r_dkstz", 200, 1., -1.);
    histos.createH1D("d_dkstz", 200, 1., -1.);
    histos.createH1D("k_dkstz", 200, 1., -1.);
    histos.createH1D("r_dkpipi", 200, 1., -1.);
    histos.createH1D("d_dkpipi", 200, 1., -1.);
    histos.createH1D("k_dkpipi", 200, 1., -1.);
    histos.createH1D("r_dpipipi", 200, 1., -1.);
    histos.createH1D("d_dpipipi", 200, 1., -1.);
    histos.createH1D("k_dpipipi", 200, 1., -1.);

    //Time dependent B decay and mixing
    histos.createH1D("l_dsk", 200, 1., -1.);
    histos.createH1D("d_dsk", 200, 1., -1.);
    histos.createH1D("phis", 200, 1., -1.);
    histos.createH1D("beta_s", 200, 1., -1.);

    histos.createH1D("l_dskpipi", 200, 1., -1.);
    histos.createH1D("d_dskpipi", 200, 1., -1.);
    histos.createH1D("k_dskpipi", 200, 1., -1.);

    histos.createH1D("l_dmpi", 200, 1., -1.);
    histos.createH1D("d_dmpi", 200, 1., -1.);
    histos.createH1D("beta", 200, 1., -1.);

    //Time dependent D decay and mixing
    histos.createH1D("PhiM12", 200, 1., -1.);
    histos.createH1D("AD", 200, 1., -1.);
    histos.createH1D("Delta_totau", 200, 1., -1.);
    histos.createH1D("DeltaAcp", 200, 1., -1.);
    histos.createH1D("phi", 200, 1., -1.);
    histos.createH1D("qopm1", 200, 1., -1.);
    histos.createH1D("qop", 200, 1., -1.);
    histos.createH1D("phi12", 200, 1., -1.);
    histos.createH1D("delta", 200, 1., -1.);
    histos.createH1D("delta_kp", 200, 1., -1.);
    histos.createH1D("delta_kpp", 200, 1., -1.);

    // 2D Histograms

    //Time integrated decay chain
    histos.createH2D("r_dk", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpi", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dk", "r_dk", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpi", "r_dpi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dk", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpi", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkstz", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkstz", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpi", "rD_kpi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y12", "x12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "x", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y12", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x12", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "g", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x", "g", 200, 1., -1., 200, 1., -1);


    //Time dependent D mixing and Decay
    histos.createH2D("y12", "qopm1", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x12", "qopm1", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "qopm1", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x", "qopm1", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "qopm1", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "x12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "y12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "x", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "y", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "y12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "y", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "x12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "x", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpi", "dD_kpi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_k3pi", "dD_k3pi", 200, 1., -1., 200, 1., -1);

    if (ckmcorr) {
        histos.createH1D("phiKsLHCb", 200, 1., -1.);
        histos.createH1D("RCP", 200, 1., -1.);
    }

    histos.createH2D("g", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x12", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y12", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dk", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dk", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_k3pi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_k3pi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_k3pi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_pipipipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_pipipi0", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpipi0", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpipi0", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_kpipi0", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_kkpi0", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kskpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kskpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_kskpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("RBRdkdpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dstk", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dstk", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dstpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dstpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkst", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkst", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkst", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkstz", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkstz", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkstz", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkpipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkpipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkpipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpipipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpipipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dpipipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dsk", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dsk", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phis", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("beta_s", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dskpipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dskpipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dskpipi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dmpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dmpi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("beta", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("AD", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("Delta_totau", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("DeltaAcp", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi12", "qop", 200, 1., -1., 200, 1., -1);
    histos.createH2D("delta", "qop", 200, 1., -1., 200, 1., -1);

    histos.createH2D("g", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x12", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y12", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dk", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dk", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_k3pi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_k3pi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_k3pi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_pipipipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_pipipi0", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpipi0", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpipi0", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_kpipi0", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_kkpi0", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kskpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kskpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_kskpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("RBRdkdpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dstk", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dstk", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dstpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dstpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkst", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkst", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkst", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkstz", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkstz", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkstz", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkpipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkpipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkpipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpipipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpipipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dpipipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dsk", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dsk", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phis", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("beta_s", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dskpipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dskpipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dskpipi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dmpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dmpi", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("beta", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("AD", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("Delta_totau", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("DeltaAcp", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phi12", "phi", 200, 1., -1., 200, 1., -1);
    histos.createH2D("delta", "phi", 200, 1., -1., 200, 1., -1);

    histos.createH2D("g", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x12", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("x", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y12", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("y", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dk", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dk", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_k3pi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_k3pi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_k3pi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_pipipipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_pipipi0", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kpipi0", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kpipi0", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_kpipi0", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("F_kkpi0", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("rD_kskpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("dD_kskpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("kD_kskpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("RBRdkdpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dstk", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dstk", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dstpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dstpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkst", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkst", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkst", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkstz", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkstz", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkstz", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dkpipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dkpipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dkpipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("r_dpipipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dpipipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dpipipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dsk", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dsk", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("phis", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("beta_s", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dskpipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dskpipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("k_dskpipi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("l_dmpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("d_dmpi", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("beta", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("PhiM12", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("AD", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("Delta_totau", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("DeltaAcp", "phi12", 200, 1., -1., 200, 1., -1);
    histos.createH2D("delta", "phi12", 200, 1., -1., 200, 1., -1);



  }



}
// ---------------------------------------------------------


double MixingModel::Acp(double rB, double delta_B, double kB , double F_D , double alpha ){
  return (2 * rB * kB * sin(g) * sin(delta_B) * ( (2 * F_D -1)  - alpha * y12) ) /
  ( (1+ rB*rB) * (1 - alpha * y12 * (2 * F_D - 1 ) )  + 2 * rB * kB * cos(g) * cos(delta_B) * ( (2 * F_D -1) - alpha * y12  )  );
}

// ---------------------------------------------------------

double MixingModel::Afav(double rB, double rD, double delta_B, double delta_D, double kB , double kD , double alpha){
  return (2 * rD * rB * kD * kB * sin(g) * sin(delta_B + delta_D) - alpha * rB * kB * sin(g) * (x12 * cos(delta_B) * (1 - rD *rD) + y12 * sin(delta_B) * (1 + rD * rD) ) ) /
  (1 + rD*rD*rB*rB + 2*rB*rD*kB*kD*cos(g)*cos(delta_B + delta_D)
  - alpha * y12 * ( rD * kD *(1 + rB * rB) * cos(delta_D) + rB * kB* (1 + rD * rD) *cos(g) * cos(delta_B) )
  + alpha * x12 * (rB * kB * (1 - rD * rD) * cos(g)* sin(delta_B) - rD * kD * (1 - rB*rB) * sin(- delta_D) ) );
}

// ---------------------------------------------------------

double MixingModel::Rcp_h(double rBCP, double delta_BCP, double rBCF, double rD, double delta_BCF, double delta_D, double kB, double kD, double F_D, double alpha){
  return ( (1 + rBCP * rBCP) * (1 - (2 * F_D - 1)*alpha*y12)  + 2 * kB * rBCP * cos(g) * cos(delta_BCP) * ( (2 * F_D - 1) - alpha * y12 ) ) /
  (1 + rD*rD * rBCF*rBCF + 2 * rBCF * rD * kB * kD * cos(g) * cos(delta_BCF + delta_D)
  - alpha * y12 * ( rD*kD*(1+rBCF*rBCF)*cos(delta_D) + kB *rBCF*(1 + rD*rD)*cos(g)*cos(delta_BCF) )
  + alpha * x12 * ( kB *rBCF*(1 - rD*rD)*cos(g)*sin(delta_BCF) - rD*kD*(1 - rBCF*rBCF)*sin(- delta_D) )   );
}

// ---------------------------------------------------------

double MixingModel::Rm(double rB, double rD, double delta_B, double delta_D, double kB, double kD, double alpha){
  return (rD*rD + rB*rB + 2*kB*kD*rB*rD*cos(delta_B - g - delta_D)
  - alpha*y12*(rD*kD*(1+rB*rB)*cos(delta_D) + rB*kB*(1+rD*rD)*cos(delta_B - g))
  - alpha*x12*(rD*kD*(1-rB*rB)*sin(delta_D) + rB*kB*(1 - rD*rD)*sin(delta_B - g)) )/
  (1 + rD*rD*rB*rB + 2*rD*rB*kD*kB*cos(delta_B - g + delta_D)
  - alpha*y12*(rD*kD*(1+rB*rB)*cos(delta_D) + rB*kB*(1+rD*rD)*cos(delta_B - g))
  + alpha*x12*(rD*kD*(1 - rB*rB)*sin(delta_D) + rB*kB*(1 - rD*rD)*sin(delta_B - g)) );
}

// ---------------------------------------------------------

double MixingModel::Rp(double rB, double rD, double delta_B, double delta_D, double kB, double kD, double alpha){
  return (rD*rD + rB*rB + 2*kB*kD*rB*rD*cos(delta_B + g - delta_D)
  - alpha*y12*(rD*kD*(1+rB*rB)*cos(delta_D) + rB*kB*(1+rD*rD)*cos(delta_B + g))
  - alpha*x12*(rD*kD*(1-rB*rB)*sin(delta_D) + rB*kB*(1 - rD*rD)*sin(delta_B + g)) )/
  (1 + rD*rD*rB*rB + 2*rD*rB*kD*kB*cos(delta_B + g + delta_D)
  - alpha*y12*(rD*kD*(1+rB*rB)*cos(delta_D) + rB*kB*(1+rD*rD)*cos(delta_B + g))
  + alpha*x12*(rD*kD*(1 - rB*rB)*sin(delta_D) + rB*kB*(1 - rD*rD)*sin(delta_B + g)) );
}

// ---------------------------------------------------------

double MixingModel::Asup(double rB, double rD, double delta_B, double delta_D, double kB, double kD , double alpha){
  return (2*rB*rD*kB*kD*sin(g)*sin(delta_B - delta_D) - alpha*rB*kB*sin(g)*( y12*(1+rD*rD)*sin(delta_B) - x12*(1-rD*rD)*cos(delta_B)) )/
  (rB*rB + rD*rD + 2*rB*rD*kB*kD*cos(g)*cos(delta_B - delta_D)
  - alpha*y12*(rB*kB*(1+rD*rD)*cos(g)*cos(delta_B) + rD*kD*(1+rB*rB)*cos(delta_D))
  - alpha*x12*(rB*kB*(1-rD*rD)*cos(g)*sin(delta_B) + rD*kD*(rB*rB-1)*sin(-delta_D)) );
}

// ---------------------------------------------------------

double MixingModel::Rads(double rB, double rD, double delta_B, double delta_D, double kB, double kD, double alpha){
  return (rD*rD + rB*rB + 2*rB*rD*kB*kD*cos(g)*cos(delta_B - delta_D)
  - alpha*y12*(rD*kD*(1+rB*rB)*cos(delta_D) + kB*rB*(1+rD*rD)*cos(g)*cos(delta_B))
  - alpha*x12*(rD*kD*(1-rB*rB)*sin(delta_D) + rB*kB*(1-rD*rD)*cos(g)*sin(delta_B)) )/
  (1+rD*rD*rB*rB + 2*rB*rD*kB*kD*cos(g)*cos(delta_B + delta_D)
  - alpha*y12*(rD*kD*(1+rB*rB)*cos(delta_D) + rB*kB*(1+rD*rD)*cos(g)*cos(delta_B))
  + alpha*x12*(rD*kD*(1-rB*rB)*sin(delta_D) + rB*kB*(1-rD*rD)*cos(g)*sin(delta_B)));
}

// ---------------------------------------------------------

double MixingModel::Rfav(double rB1, double rB2, double rD, double delta_B1, double delta_B2, double delta_D, double BR, double kD, double alpha){
  return BR * (
    (1 + rB1*rB1*rD*rD + 2*kD*rD*rB1*cos(g)*cos(delta_B1 + delta_D)
    - alpha*y12*(rD*kD*(1+rB1*rB1)*cos(delta_D) + rB1*(1+rD*rD)*cos(g)*cos(delta_B1))
    + alpha*x12*(rD*kD*(1-rB1*rB1)*sin(delta_D) + rB1*(1 - rD*rD)*cos(g)*sin(delta_B1) ) ) /
    (1 + rB2*rB2*rD*rD + 2*kD*rD*rB2*cos(g)*cos(delta_B2 + delta_D)
    - alpha*y12*(rD*kD*(1+rB2*rB2)*cos(delta_D) + rB2*(1+rD*rD)*cos(g)*cos(delta_B2))
    + alpha*x12*(rD*kD*(1-rB2*rB2)*sin(delta_D) + rB2*(1 - rD*rD)*cos(g)*sin(delta_B2) ) )
  );
}

// ---------------------------------------------------------

double MixingModel::Rsup(double rB1, double rB2, double rD, double delta_B1, double delta_B2, double delta_D, double BR, double kD, double alpha){
  return BR * (
    (rD*rD + rB1*rB1 + 2*rD*kD*rB1*cos(g)*cos(delta_B1 - delta_D)
    - alpha*y12*(rD*kD*(1+rB1*rB1)*cos(delta_D) + rB1*(1+rD*rD)*cos(g)*cos(delta_B1) )
    - alpha*x12*(rD*kD*(1 - rB1*rB1)*sin(delta_D) + rB1*(1 - rD*rD)*cos(g)*sin(delta_B1) ) )/
    (rD*rD + rB2*rB2 + 2*rD*kD*rB2*cos(g)*cos(delta_B2 - delta_D)
    - alpha*y12*(rD*kD*(1+rB2*rB2)*cos(delta_D) + rB2*(1+rD*rD)*cos(g)*cos(delta_B2))
    - alpha*x12*(rD*kD*(1 - rB2*rB2)*sin(delta_D) + rB2*(1 - rD*rD)*cos(g)*sin(delta_B2) ) )
  );
}


// --------------------------------------------------------- //phi = phiD di LHCb

double MixingModel::y_plus(double delta_D){
  return qop * (   sin(phi) * (x * cos(delta_D) - y * sin(delta_D) ) - cos(phi) * (x * sin(delta_D) + y*cos(delta_D))  );
}

// ---------------------------------------------------------

double MixingModel::y_minus(double delta_D){
  return ( 1./qop) * ( - sin(phi) * (x * cos(delta_D) - y * sin(delta_D) ) - cos(phi) * (x * sin(delta_D) + y*cos(delta_D))  );
}

// ---------------------------------------------------------

double MixingModel::x_plus(double delta_D){
  return ( qop) * ( - cos(phi) * (x * cos(delta_D) - y * sin(delta_D) ) - sin(phi) * (x * sin(delta_D) + y*cos(delta_D))  );
}

// ---------------------------------------------------------

double MixingModel::x_minus(double delta_D){
  return ( -1./qop) * (  cos(phi) * (x * cos(delta_D) - y * sin(delta_D) ) - sin(phi) * (x * sin(delta_D) + y*cos(delta_D))  );
}

// ---------------------------------------------------------
double MixingModel::Calculate_time_integrated_observables(){

  double ll1;
  ll1 = 0.;

  //------------------------------------------------------ Calculating time integrated Observables ----------------------------------------------------

  // 1. PDF: glwads-dh-hh-dmix (UID0)
  // 8 Observables
  acp_dk_uid0 = Acp(r_dk, d_dk, 1., 1., 2 * 0.523);
  acp_dpi_uid0 = Acp(r_dpi, d_dpi,1., 1.,  2 * 0.523);
  afav_dk_uid0 = Afav(r_dk, rD_kpi, d_dk, dD_kpi, 1.,1.,2 * 0.523);
  rcp_uid0 = Rcp_h(r_dk, d_dk, r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 1., 2 * 0.523  ) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_kpi, d_dpi, dD_kpi,1., 1., 1., 2 * 0.523  ) ;
  rm_dk_uid0 = Rm(r_dk, rD_kpi, d_dk, dD_kpi, 1.,1., 2 * 0.523);
  rm_dpi_uid0 = Rm(r_dpi, rD_kpi, d_dpi, dD_kpi, 1.,1., 2*0.523);
  rp_dk_uid0 = Rp(r_dk,rD_kpi,d_dk,dD_kpi,1.,1.,2*0.523);
  rp_dpi_uid0 = Rp(r_dpi,rD_kpi,d_dpi,dD_kpi,1.,1.,2*0.523);

  // 2. PDF: glwads-dh-h3pi-dmix (UID1)
  // 8 Observables
  aads_dk_k3pi_uid1 = Asup(r_dk,rD_k3pi,d_dk,dD_k3pi,1.,kD_k3pi, 2*0.570);
  aads_dpi_k3pi_uid1 = Asup(r_dpi, rD_k3pi, d_dpi, dD_k3pi, 1., kD_k3pi, 2*0.570);
  acp_dk_4pi_uid1 = Acp(r_dk,d_dk,1.,F_pipipipi, 2*0.570);
  acp_dpi_4pi_uid1 = Acp(r_dpi, d_dpi, 1., F_pipipipi, 2*0.570);
  afav_dk_k3pi_uid1 = Afav(r_dk,rD_k3pi,d_dk,dD_k3pi,1.,kD_k3pi,2*0.570);
  rads_dk_k3pi_uid1 = Rads(r_dk, rD_k3pi, d_dk, dD_k3pi, 1., kD_k3pi, 2*0.570);
  rads_dpi_k3pi_uid1 = Rads(r_dpi, rD_k3pi, d_dpi, dD_k3pi, 1., kD_k3pi, 2*0.570);
  rcp_4pi_uid1 = Rcp_h(r_dk, d_dk, r_dk, rD_k3pi, d_dk, dD_k3pi, 1., kD_k3pi, F_pipipipi, 2*0.57) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_k3pi, d_dpi, dD_k3pi, 1., kD_k3pi, F_pipipipi, 2*0.57);

  // 3. PDF: glwads-dh-hhpi0-dmix (UID2)
  // 11 Observables
  aads_dk_kpipi0_uid2 = Asup(r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, 2*0.592);
  aads_dpi_kpipi0_uid2 = Asup(r_dpi, rD_kpipi0, d_dpi, dD_kpipi0, 1., kD_kpipi0, 2*0.592);
  acp_dk_kkpi0_uid2 = Acp(r_dk, d_dk, 1., F_kkpi0, 2*0.592);
  acp_dk_pipipi0_uid2 = Acp(r_dk,d_dk, 1., F_pipipi0, 2*0.592);
  acp_dpi_kkpi0_uid2 = Acp(r_dpi, d_dpi, 1., F_kkpi0, 2*0.592);
  acp_dpi_pipipi0_uid2 = Acp(r_dpi, d_dpi, 1., F_pipipi0, 2 * 0.592);
  afav_dk_kpipi0_uid2 = Afav(r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, 2*0.592);
  rads_dk_kpipi0_uid2 = Rads(r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, 2*0.592);
  rads_dpi_kpipi0_uid2 = Rads(r_dpi, rD_kpipi0, d_dpi, dD_kpipi0, 1., kD_kpipi0, 2*0.592 );
  rcp_kkpi0_uid2 = Rcp_h(r_dk, d_dk, r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, F_kkpi0, 2*0.592) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_kpipi0, d_dpi, dD_kpipi0, 1., kD_kpipi0,  F_kkpi0, 2*0.592);
  rcp_pipipi0_uid2 = Rcp_h(r_dk, d_dk, r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, F_pipipi0, 2*0.592) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_kpipi0, d_dpi, dD_kpipi0, 1., kD_kpipi0,  F_pipipi0, 2*0.592);

  //4. PDF: ggsz-dh (UID3)
  // 4 Observables
  xm_dk_uid3 = r_dk * cos(d_dk - g);
  ym_dk_uid3 = r_dk * sin(d_dk - g);
  xp_dk_uid3 = r_dk * cos(d_dk + g);
  yp_dk_uid3 = r_dk * sin(d_dk + g);
  xi_x_dpi_uid3 = (r_dpi/r_dk) * cos(d_dpi - d_dk);
  xi_y_dpi_uid3 = (r_dpi/r_dk) * sin(d_dpi - d_dk);

  //5. PDF: glwads-dkdpi-kskpi-dmix (UID4)
  // 7 Observables
  afav_dpi_kskpi_uid4 = Afav(r_dpi, rD_kskpi, d_dpi, dD_kskpi, 1., kD_kskpi, 1.);
  asup_dpi_kskpi_uid4 = Asup(r_dpi, rD_kskpi, d_dpi, dD_kskpi, 1., kD_kskpi, 1.);
  afav_dk_kskpi_uid4 = Afav(r_dk, rD_kskpi, d_dk, dD_kskpi, 1., kD_kskpi, 1.);
  asup_dk_kskpi_uid4 = Asup(r_dk, rD_kskpi, d_dk, dD_kskpi, 1., kD_kskpi, 1.);
  rfavsup_dpi_kskpi_uid4 = 1./Rads(r_dpi, rD_kskpi, d_dpi, dD_kskpi, 1., kD_kskpi, 1.);
  rfav_dkdpi_kskpi_uid4 = Rfav(r_dk, r_dpi, rD_kskpi, d_dk, d_dpi, dD_kskpi, RBRdkdpi, kD_kskpi, 1.);
  rsup_dkdpi_kskpi_uid4 = Rsup(r_dk, r_dpi, rD_kskpi, d_dk, d_dpi, dD_kskpi, RBRdkdpi, kD_kskpi, 1.);

  // 6. PDF: glwads-dsth-hh-dmix (UID5)
  // 16 Observables
  acp_dstk_dg_uid5 = Acp(r_dstk, d_dstk + M_PI, 1., 1., 2*0.523 );
  acp_dstk_dp_uid5 = Acp(r_dstk, d_dstk, 1., 1., 2*0.523);
  afav_dstk_dg_uid5 = Afav(r_dstk, rD_kpi, d_dstk + M_PI,dD_kpi, 1., 1., 2 * 0.523);
  afav_dstk_dp_uid5 = Afav(r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 2 * 0.523);
  rcp_dg_uid5 = Rcp_h(r_dstk, d_dstk + M_PI, r_dstk, rD_kpi, d_dstk + M_PI, dD_kpi, 1., 1., 1., 2 * 0.523) / Rcp_h(r_dstpi, d_dstpi + M_PI, r_dstpi, rD_kpi, d_dstpi + M_PI, dD_kpi, 1., 1., 1., 2 * 0.523);
  rcp_dp_uid5 = Rcp_h(r_dstk, d_dstk, r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 1., 2 * 0.523) / Rcp_h(r_dstpi, d_dstpi, r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 1., 2 * 0.523);
  rm_dstk_dg_uid5 = Rm(r_dstk, rD_kpi, d_dstk + M_PI, dD_kpi, 1.,1., 2*0.523);
  rm_dstk_dp_uid5 = Rm(r_dstk, rD_kpi, d_dstk, dD_kpi, 1.,1., 2*0.523);
  rp_dstk_dg_uid5 = Rp(r_dstk, rD_kpi, d_dstk + M_PI, dD_kpi, 1.,1., 2*0.523);
  rp_dstk_dp_uid5 = Rp(r_dstk, rD_kpi, d_dstk, dD_kpi, 1.,1., 2*0.523);
  acp_dstpi_dg_uid5 = Acp(r_dstpi, d_dstpi + M_PI, 1., 1., 2*0.523 );
  acp_dstpi_dp_uid5 = Acp(r_dstpi, d_dstpi, 1., 1., 2*0.523);
  rm_dstpi_dg_uid5 = Rm(r_dstpi, rD_kpi, d_dstpi + M_PI, dD_kpi, 1.,1., 2*0.523);
  rm_dstpi_dp_uid5 = Rm(r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1.,1., 2*0.523);
  rp_dstpi_dg_uid5 = Rp(r_dstpi, rD_kpi, d_dstpi + M_PI, dD_kpi, 1.,1., 2*0.523);
  rp_dstpi_dp_uid5 = Rp(r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1.,1., 2*0.523);

  // 7. PDF: glwads-dkst-hh-h3pi-dmix-newvars (UID6)
  // 12 Observables
  afav_dkst_kpi_uid6 = Afav(r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 2 * 0.594 );
  acp_dkst_kk_uid6 = Acp(r_dkst, d_dkst, k_dkst, 1., 2*0.594);
  acp_dkst_pipi_uid6 = Acp(r_dkst, d_dkst, k_dkst, 1., 2*0.594);
  rcp_dkst_kk_uid6 = Rcp_h(r_dkst, d_dkst, r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 1., 2*0.594);
  rcp_dkst_pipi_uid6 = Rcp_h(r_dkst, d_dkst, r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 1., 2*0.594);
  rp_dkst_kpi_uid6 = Rp(r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 2*0.594);
  rm_dkst_kpi_uid6 = Rm(r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 2*0.594);
  afav_dkst_k3pi_uid6 = Afav(r_dkst, rD_k3pi, d_dkst, dD_k3pi, k_dkst, kD_k3pi, 2 * 0.594);
  acp_dkst_pipipipi_uid6 = Acp(r_dkst, d_dkst, k_dkst, F_pipipipi, 2 * 0.594);
  rcp_dkst_pipipipi_uid6 = Rcp_h(r_dkst, d_dkst, r_dkst, rD_k3pi, d_dkst, dD_k3pi, k_dkst, kD_k3pi, F_pipipipi, 2*0.594);
  rp_dkst_k3pi_uid6 = Rp(r_dkst, rD_k3pi, d_dkst, dD_k3pi, k_dkst, kD_k3pi, 2 * 0.594);
  rm_dkst_k3pi_uid6 = Rm(r_dkst, rD_k3pi, d_dkst, dD_k3pi, k_dkst, kD_k3pi, 2 * 0.594);

  // 8. PDF: glwads-dkstz-hh-h3pi-dmix (UID7)
  //12 Observables
  acp_dkstz_kk_uid7 = Acp(r_dkstz, d_dkstz, k_dkstz, 1., 2*0.6);
  acp_dkstz_pipi_uid7 = Acp(r_dkstz, d_dkstz, k_dkstz, 1., 2 * 0.6);
  rcp_dkstz_kk_uid7 = Rcp_h(r_dkstz, d_dkstz, r_dkstz, rD_kpi, d_dkstz, dD_kpi, k_dkstz, 1., 1.,  2 * 0.6);
  rcp_dkstz_pipi_uid7 = Rcp_h(r_dkstz, d_dkstz, r_dkstz, rD_kpi, d_dkstz, dD_kpi, k_dkstz, 1.,1.,  2 * 0.6);
  acp_dkstz_4pi_uid7 = Acp(r_dkstz, d_dkstz, k_dkstz, F_pipipipi, 2 * 0.6);
  rcp_dkstz_4pi_uid7 = Rcp_h(r_dkstz, d_dkstz, r_dkstz, rD_k3pi, d_dkstz, dD_k3pi, k_dkstz, kD_k3pi, F_pipipipi, 2 * 0.6 );
  rp_dkstz_kpi_uid7 = Rp(r_dkstz, rD_kpi, d_dkstz, dD_kpi, k_dkstz, 1., 2 * 0.6);
  rm_dkstz_kpi_uid7 = Rm(r_dkstz, rD_kpi, d_dkstz, dD_kpi, k_dkstz, 1., 2 * 0.6);
  rp_dkstz_k3pi_uid7 = Rp(r_dkstz, rD_k3pi, d_dkstz, dD_k3pi, k_dkstz, kD_k3pi, 2 * 0.6);
  rm_dkstz_k3pi_uid7 = Rm(r_dkstz, rD_k3pi, d_dkstz, dD_k3pi, k_dkstz, kD_k3pi, 2 * 0.6);
  afav_dkstz_kpi_uid7 = Afav(r_dkstz, rD_kpi, d_dkstz, dD_kpi, k_dkstz, 1., 2 * 0.6);
  afav_dkstz_k3pi_uid7 = Afav(r_dkstz, rD_k3pi, d_dkstz, dD_k3pi, k_dkstz, kD_k3pi, 2 * 0.6);

  //9. PDF: ggsz_dkstz_lhcb_md (UID8)
  // 4 Observables
  xm_dkstz_uid8 = k_dkstz * r_dkstz * cos(d_dkstz - g);
  ym_dkstz_uid8 = k_dkstz * r_dkstz * sin(d_dkstz - g);
  xp_dkstz_uid8 = k_dkstz * r_dkstz * cos(d_dkstz + g);
  yp_dkstz_uid8 = k_dkstz * r_dkstz * sin(d_dkstz + g);

  //10. PDF: glwads-dhpipi-hh-dmix (UID9)
  // 11 Observables
  rcp_dkpipi_uid9 = Rcp_h(r_dkpipi, d_dkpipi, r_dkpipi, rD_kpi, d_dkpipi, dD_kpi, k_dkpipi, 1., 1.,  2 * 0.6) / Rcp_h(r_dpipipi, d_dpipipi, r_dpipipi, rD_kpi, d_dpipipi, dD_kpi, k_dpipipi, 1., 1., 2 * 0.6);
  afav_dkpipi_kpi_uid9 = Afav(r_dkpipi, rD_kpi, d_dkpipi, dD_kpi, k_dkpipi, 1., 2 * 0.6);
  afav_dpipipi_kpi_uid9 = Afav(r_dpipipi, rD_kpi, d_dpipipi, dD_kpi, k_dpipipi, 1., 2*0.6);
  acp_dkpipi_kk_uid9 = Acp(r_dkpipi, d_dkpipi, k_dkpipi, 1., 2 * 0.6);
  acp_dkpipi_pipi_uid9 = Acp(r_dkpipi, d_dkpipi, k_dkpipi, 1., 2*0.6);
  acp_dpipipi_kk_uid9 = Acp(r_dpipipi, d_dpipipi, k_dpipipi, 1., 2*0.6);
  acp_dpipipi_pipi_uid9 = Acp(r_dpipipi, d_dpipipi, k_dpipipi, 1., 2*0.6);
  rp_dkpipi_uid9 = Rp(r_dkpipi, rD_kpi, d_dkpipi, dD_kpi, k_dkpipi, 1., 2 * 0.6);
  rm_dkpipi_uid9 = Rm(r_dkpipi, rD_kpi, d_dkpipi, dD_kpi, k_dkpipi, 1., 2*0.6);
  rp_dpipipi_uid9 = Rp(r_dpipipi, rD_kpi, d_dpipipi, dD_kpi, k_dpipipi, 1., 2*0.6);
  rm_dpipipi_uid9 = Rm(r_dpipipi, rD_kpi, d_dpipipi, dD_kpi, k_dpipipi, 1., 2*0.6);

  //------------------------------------------------------ Calculating the contirbutions of other useful inputs ----------------------------------------------------

  //20. PDF: dk3pi_dkpipi0_constraints (UID19)
  // 6 Observables
  kD_k3pi_uid19 = kD_k3pi;
  dD_k3pi_uid19 = 2*M_PI - dD_k3pi; //because of the phase sign convention
  kD_kpipi0_uid19 = kD_kpipi0;
  dD_kpipi0_uid19 = 2*M_PI - dD_kpipi0; //because of the phase sign convention
  rD_k3pi_uid19 = rD_k3pi;
  rD_kpipi0_uid19 = rD_kpipi0;

  //21. PDF: d4pi_dmixing_cleo (UID20)
  // 1 Osservabile
  F_pipipipi_uid20 = F_pipipipi;

  //22. PDF: CleoDhhpi0Dilution (UID21)
  F_pipipi0_uid21 = F_pipipi0;
  F_kkpi0_uid21 = F_kkpi0;

  //25. PDF: dkstcoherence (UID24)
  // 1 Osservabile
  k_dkst_uid24 = k_dkst;

  //26. PDF: dkstzcoherence (UID25)
  // 1 Osservabile
  k_dkstz_uid25 = k_dkstz;

  //------------------------------------------------------ Calculating the contribution to the LogLikelihood ----------------------------------------------------

  TVectorD corr(8);

  //1. PDF: glwads-dh-hh-dmix (UID0)
  //Observables 8:
  corr(0) = acp_dk_uid0;
  corr(1) = acp_dpi_uid0;
  corr(2) = afav_dk_uid0;
  corr(3) = rcp_uid0;
  corr(4) = rm_dk_uid0;
  corr(5) = rm_dpi_uid0;
  corr(6) = rp_dk_uid0;
  corr(7) = rp_dpi_uid0;

  ll1 += corrmeas.at("UID0").logweight(corr);

  //2. PDF: glwads-dh-h3pi-dmix (UID1)
  //Observables 8:
  corr.ResizeTo(8);
  corr(0) = aads_dk_k3pi_uid1;
  corr(1) = aads_dpi_k3pi_uid1;
  corr(2) = acp_dk_4pi_uid1;
  corr(3) = acp_dpi_4pi_uid1;
  corr(4) = afav_dk_k3pi_uid1;
  corr(5) = rads_dk_k3pi_uid1;
  corr(6) = rads_dpi_k3pi_uid1;
  corr(7) = rcp_4pi_uid1;

  ll1 += corrmeas.at("UID1").logweight(corr);

  //3. PDF: glwads-dh-hhpipi0-dmix (UID2)
  //Observables 11:
  corr.ResizeTo(11);
  corr(0) = aads_dk_kpipi0_uid2;
  corr(1) = aads_dpi_kpipi0_uid2;
  corr(2) = acp_dk_kkpi0_uid2;
  corr(3) = acp_dk_pipipi0_uid2;
  corr(4) = acp_dpi_kkpi0_uid2;
  corr(5) = acp_dpi_pipipi0_uid2;
  corr(6) = afav_dk_kpipi0_uid2;
  corr(7) = rads_dk_kpipi0_uid2;
  corr(8) = rads_dpi_kpipi0_uid2;
  corr(9) = rcp_kkpi0_uid2;
  corr(10) = rcp_pipipi0_uid2;

  ll1 += corrmeas.at("UID2").logweight(corr);

  //4. PDF: ggsz-dh (UID3)
  //Observables 6:
  corr.ResizeTo(6);
  corr(0) = xm_dk_uid3;
  corr(1) = ym_dk_uid3;
  corr(2) = xp_dk_uid3;
  corr(3) = yp_dk_uid3;
  corr(4) = xi_x_dpi_uid3;
  corr(5) = xi_y_dpi_uid3;

  ll1 += corrmeas.at("UID3").logweight(corr);

  //5. PDF: glwads-dkdpi-kskpi-dmix (UID4)
  //Observables 7:
  corr.ResizeTo(7);
  corr(0) = afav_dpi_kskpi_uid4;
  corr(1) = asup_dpi_kskpi_uid4;
  corr(2) = afav_dk_kskpi_uid4;
  corr(3) = asup_dk_kskpi_uid4;
  corr(4) = rfavsup_dpi_kskpi_uid4;
  corr(5) = rfav_dkdpi_kskpi_uid4;
  corr(6) = rsup_dkdpi_kskpi_uid4;

  ll1 += corrmeas.at("UID4").logweight(corr);

  //6. PDF: glwads-dsth-hh-dmix (UID5)
  //Observables 16:
  corr.ResizeTo(16);
  corr(0) = acp_dstk_dg_uid5;
  corr(1) = acp_dstk_dp_uid5;
  corr(2) = afav_dstk_dg_uid5;
  corr(3) = afav_dstk_dp_uid5;
  corr(4) = rcp_dg_uid5;
  corr(5) = rcp_dp_uid5;
  corr(6) = rm_dstk_dg_uid5;
  corr(7) = rm_dstk_dp_uid5;
  corr(8) = rp_dstk_dg_uid5;
  corr(9) = rp_dstk_dp_uid5;
  corr(10) = acp_dstpi_dg_uid5;
  corr(11) = acp_dstpi_dp_uid5;
  corr(12) = rm_dstpi_dg_uid5;
  corr(13) = rm_dstpi_dp_uid5;
  corr(14) = rp_dstpi_dg_uid5;
  corr(15) = rp_dstpi_dp_uid5;

  ll1 += corrmeas.at("UID5").logweight(corr);

  //7. PDF: glwads-dkst-hh-h3pi-dmix (UID6)
  //Observables 12:
  corr.ResizeTo(12);
  corr(0) = afav_dkst_kpi_uid6;
  corr(1) = acp_dkst_kk_uid6;
  corr(2) = acp_dkst_pipi_uid6;
  corr(3) = rcp_dkst_kk_uid6;
  corr(4) = rcp_dkst_pipi_uid6;
  corr(5) = rp_dkst_kpi_uid6;
  corr(6) = rm_dkst_kpi_uid6;
  corr(7) = afav_dkst_k3pi_uid6;
  corr(8) = acp_dkst_pipipipi_uid6;
  corr(9) = rcp_dkst_pipipipi_uid6;
  corr(10) = rp_dkst_k3pi_uid6;
  corr(11) = rm_dkst_k3pi_uid6;

  ll1 += corrmeas.at("UID6").logweight(corr);


  //8. PDF: glwads-dkst-hh-h3pi-dmix (UID7)
  //Observables 12:
  corr.ResizeTo(12);
  corr(0) = acp_dkstz_kk_uid7;
  corr(1) = acp_dkstz_pipi_uid7;
  corr(2) = rcp_dkstz_kk_uid7;
  corr(3) = rcp_dkstz_pipi_uid7;
  corr(4) = acp_dkstz_4pi_uid7;
  corr(5) = rcp_dkstz_4pi_uid7;
  corr(6) = rp_dkstz_kpi_uid7;
  corr(7) = rm_dkstz_kpi_uid7;
  corr(8) = rp_dkstz_k3pi_uid7;
  corr(9) = rm_dkstz_k3pi_uid7;
  corr(10) = afav_dkstz_kpi_uid7;
  corr(11) = afav_dkstz_k3pi_uid7;

  ll1 += corrmeas.at("UID7").logweight(corr);


  //9. PDF: ggsz-dkstz-lhcb-md (UID8)
  //Observables 4:
  corr.ResizeTo(4);
  corr(0) = xm_dkstz_uid8;
  corr(1) = ym_dkstz_uid8;
  corr(2) = xp_dkstz_uid8;
  corr(3) = yp_dkstz_uid8;

  ll1 += corrmeas.at("UID8").logweight(corr);

  //10. PDF: glwads-dhpipi-hh-dmix (UID9)
  //Observables 11:
  corr.ResizeTo(11);
  corr(0) = rcp_dkpipi_uid9;
  corr(1) = afav_dkpipi_kpi_uid9;
  corr(2) = afav_dpipipi_kpi_uid9;
  corr(3) = acp_dkpipi_kk_uid9;
  corr(4) = acp_dkpipi_pipi_uid9;
  corr(5) = acp_dpipipi_kk_uid9;
  corr(6) = acp_dpipipi_pipi_uid9;
  corr(7) = rp_dkpipi_uid9;
  corr(8) = rm_dkpipi_uid9;
  corr(9) = rp_dpipipi_uid9;
  corr(10) = rm_dpipipi_uid9;

  ll1 += corrmeas.at("UID9").logweight(corr);

  //20. PDF: dk3pi_dkpipi0_constraints (UID19)
  //Observables 6:
  corr.ResizeTo(6);
  corr(0) = kD_k3pi_uid19;
  corr(1) = 2*M_PI - dD_k3pi_uid19;
  corr(2) = kD_kpipi0_uid19;
  corr(3) = 2*M_PI - dD_kpipi0_uid19;
  corr(4) = rD_k3pi_uid19;
  corr(5) = rD_kpipi0_uid19;

  ll1 += corrmeas.at("UID19").logweight(corr);

  //21. PDF: d4pi_dmixing_cleo (UID20)
  //Observables 1:

  ll1 += meas.at("UID20").logweight(F_pipipipi_uid20);

  //22. PDF: CleoDhhpi0Diluition (UID21)
  //Observables 2:
  corr.ResizeTo(2);
  corr(0) = F_pipipi0;
  corr(1) = F_kkpi0;

  ll1 += corrmeas.at("UID21").logweight(corr);

  //25. PDF: dkstcoherence (UID24)
  //Observables 1:

  ll1 += meas.at("UID24").logweight(k_dkst_uid24);

  //26. PDF: dkstzcoherence (UID25)
  //Observables 1:

  ll1 += meas.at("UID25").logweight(k_dkstz_uid25);

  return ll1;

}
// ---------------------------------------------------------

// ---------------------------------------------------------
double MixingModel::Calculate_time_dependent_Bobservables(){

  double ll2;
  ll2 = 0.;

  //----------------------------------------------- Calculating time dependent B observables -------------------------------------------------------------------------

  //11. PDF: dsk (UID10)
  // 5 Observables
  c_dsk_uid10 = (1 - l_dsk * l_dsk)/(1 + l_dsk * l_dsk);
  d_dsk_uid10 = - (2*l_dsk*cos(d_dsk - (g+phis))) /(1 + l_dsk*l_dsk); //phis = - 2 beta_s
  db_dsk_uid10 = - (2*l_dsk*cos(d_dsk + (g + phis))) / (1 + l_dsk*l_dsk);
  s_dsk_uid10 = (2*l_dsk*sin(d_dsk - (g+phis)))/(1 + l_dsk*l_dsk);
  sb_dsk_uid10 = - (2*l_dsk*sin(d_dsk + (g+phis))) / (1+l_dsk*l_dsk);

  //12. PDF: dskpipi (UID11)
  // 5 Observables
  c_dskpipi_uid11 = (1 - l_dskpipi * l_dskpipi) / (1 + l_dskpipi * l_dskpipi);
  d_dskpipi_uid11 =  - (2*k_dskpipi*l_dskpipi*cos(d_dskpipi - (g+phis))) /(1 + l_dskpipi*l_dskpipi);
  db_dskpipi_uid11 = - (2*k_dskpipi*l_dskpipi*cos(d_dskpipi + (g + phis))) / (1 + l_dskpipi*l_dskpipi);
  s_dskpipi_uid11 = (2*k_dskpipi*l_dskpipi*sin(d_dskpipi - (g+phis)))/(1 + l_dskpipi*l_dskpipi);
  sb_dskpipi_uid11 = - (2*k_dskpipi*l_dskpipi*sin(d_dskpipi + (g+phis))) / (1+l_dskpipi*l_dskpipi);

  //13. PDF: dmpi (UID12)
  // 2 Observables
  s_dmpi_uid12 = - (2 * l_dmpi * sin(d_dmpi - (2 *beta + g))) / (1 + l_dmpi*l_dmpi);
  sb_dmpi_uid12 = (2 * l_dmpi * sin(d_dmpi + (2 *beta + g))) / (1 + l_dmpi*l_dmpi);

  //27. PDF: phis (UID26)
  // 1 Osservabile
  phis_uid26 = phis;

  //28. PDF: beta (UID27)
  // 1 Osservabile
  beta_uid27 = beta;

  //----------------------------------------------- Contribution to the LogLikelihood -------------------------------------------------------------------------

  TVectorD corr(8);

  //11. PDF: dsk (UID10)
  //Observables 5:
  corr.ResizeTo(5);
  corr(0) = c_dsk_uid10;
  corr(1) = d_dsk_uid10;
  corr(2) = db_dsk_uid10;
  corr(3) = s_dsk_uid10;
  corr(4) = sb_dsk_uid10;

  ll2 += corrmeas.at("UID10").logweight(corr);

  //12. PDF: dskpipi (UID11)
  //Observables 5:
  corr.ResizeTo(5);
  corr(0) = c_dskpipi_uid11;
  corr(1) = d_dskpipi_uid11;
  corr(2) = db_dskpipi_uid11;
  corr(3) = s_dskpipi_uid11;
  corr(4) = sb_dskpipi_uid11;

  ll2 += corrmeas.at("UID11").logweight(corr);

  //13. PDF: dmpi (UID12)
  //Observables 2:
  corr.ResizeTo(2);
  corr(0) = s_dmpi_uid12;
  corr(1) = sb_dmpi_uid12;

  ll2 += corrmeas.at("UID12").logweight(corr);

  //27. PDF: phis (UID26)
  //Observables 1:

  ll2 += meas.at("UID26").logweight(phis_uid26);

  //28. PDF: beta (UID27)
  //Observables 1:

  ll2 += meas.at("UID27").logweight(beta_uid27);


  return ll2;

}
// ---------------------------------------------------------

// ---------------------------------------------------------
double MixingModel::Calculate_time_dependent_Dobservables(){

  double ll3;
  ll3 = 0.;

  //------------------------------------------------------  Time Dependent D observables  ----------------------------------------------------

  //15. PDF: charm-kspipi (UID14)
  // 4 Observables
  xcp_uid14 = xcp;
  ycp_uid14 = ycp;
  dx_uid14 = dx;
  dy_uid14 = dy;

  //16. PDF: charm-kspipi (UID15)
  // 4 Observables
  xcp_uid15 = xcp;
  ycp_uid15 = ycp;
  dx_uid15 = dx;
  dy_uid15 = dy;

  //17. PDF: charm-kpi (UID16)
  // 6 Observables
  Rdp_uid16 = rD_kpi * rD_kpi * (1. + AD);
  yp_uid16 = y_plus(dD_kpi);
  xpsq_uid16 = x_plus(dD_kpi) * x_plus(dD_kpi);
  Rdm_uid16 = rD_kpi * rD_kpi * (1. - AD);
  ym_uid16 = y_minus(dD_kpi);
  xmsq_uid16 = x_minus(dD_kpi) * x_minus(dD_kpi);

  //18. PDF: charm-deltaacp-diff (UID17)
  // 2 Observables
  dtotau_uid17 = Delta_totau;
  dacp_uid17 = DeltaAcp - 0.5*Delta_totau*(y*cos(phi)*(qop - 1./qop) - x*sin(phi)*(qop + 1./qop));

  //29. PDF: Charm_ycp_minus_ycp_rs_lhcb-r1 (UID28)
  // 1 Osservabile
  ycp_uid28 = ycp + 0.5 * rD_kpi * ( cos(phi)*(-y*cos(dD_kpi) + x*sin(dD_kpi))*(qop + 1./qop) + sin(phi)*(y*sin(dD_kpi)+x*cos(dD_kpi))*(qop - 1./qop) );

  //30. PDF: charm-dy-rs (UID29)
  // 1 Osservabile
  DY_uid29 = 0.5*( -y*cos(phi)*(qop - 1./qop) + x*sin(phi)*(qop + 1./qop) );

  //31. PDF: charm-kpi (UID30)
  // 6 Observables
  Rdp_uid30 = rD_kpi * rD_kpi * (1 + AD);
  yp_uid30 = y_plus(dD_kpi);
  xpsq_uid30 = x_plus(dD_kpi) * x_plus(dD_kpi);
  Rdm_uid30 = rD_kpi * rD_kpi * (1-AD);
  ym_uid30 = y_minus(dD_kpi);
  xmsq_uid30 = x_minus(dD_kpi) * x_minus(dD_kpi);

  if(comb == 2){

  // https://arxiv.org/abs/2208.09402
  Akpi_BESIII = (- 2 * rD_kpi * cos(dD_kpi) + y) / (1 + rD_kpi * rD_kpi);
  Akpi_kpipi0_BESIII = ( F_pipipi0 * (- 2 * rD_kpi * cos(dD_kpi) + y )  ) / ( 1 + rD_kpi * rD_kpi + (1 - F_pipipi0)  * ( 2 * rD_kpi * cos(dD_kpi) + y  ) );

  xi_x_BESIII = rD_kpi * cos(dD_kpi);
  xi_y_BESIII = - rD_kpi * sin(dD_kpi); // Because of the different phase convention dD = (-dD)BESIII

  //https://arxiv.org/pdf/2208.06512.pdf
  //they are the same of xcp, ycp, dx, dy

  }

  //----------------------------------------------- Contribution to the LogLikelihood -------------------------------------------------------------------------
  TVectorD corr(8);

  //15. PDF: charm-kskpi (UID14)
  //Observables 4:
  corr.ResizeTo(4);
  corr(0) = xcp_uid14;
  corr(1) = ycp_uid14;
  corr(2) = dx_uid14;
  corr(3) = dy_uid14;

  ll3 += corrmeas.at("UID14").logweight(corr);

  //16. PDF: charm-kskpi (UID15)
  //Observables 4:
  corr.ResizeTo(4);
  corr(0) = xcp_uid15;
  corr(1) = ycp_uid15;
  corr(2) = dx_uid15;
  corr(3) = dy_uid15;

  ll3 += corrmeas.at("UID15").logweight(corr);

  //17. PDF: charm-kpi (UID16)
  //Observables 6:
  corr.ResizeTo(6);
  corr(0) = Rdp_uid16;
  corr(1) = yp_uid16;
  corr(2) = xpsq_uid16;
  corr(3) = Rdm_uid16;
  corr(4) = ym_uid16;
  corr(5) = xmsq_uid16;

  ll3 += corrmeas.at("UID16").logweight(corr);


  //18. PDF: charm-deltaacp-diff (UID17)
  //Observables 2:
  corr.ResizeTo(2);
  corr(0) = dtotau_uid17;
  corr(1) = dacp_uid17;

  ll3 += corrmeas.at("UID17").logweight(corr);

  //29. PDF: Charm_ycp_minus_ycp_rs_lhcb-r1 (UID28)
  //Observables 1:
  ll3 += meas.at("UID28").logweight(ycp_uid28);

  //30. PDF: charm-dy-rs (UID29)
  //Observables 1:
  ll3 += meas.at("UID29").logweight(DY_uid29);

  //31. PDF: charm-kpi (UID30)
  //Observables 6:
  corr.ResizeTo(6);
  corr(0) = Rdp_uid30;
  corr(1) = yp_uid30;
  corr(2) = xpsq_uid30;
  corr(3) = Rdm_uid30;
  corr(4) = ym_uid30;
  corr(5) = xmsq_uid30;

  ll3 += corrmeas.at("UID30").logweight(corr);

  if(comb == 2){

  // https://arxiv.org/abs/2208.09402
  corr.ResizeTo(2);
  corr(0) = Akpi_BESIII;
  corr(1) = Akpi_kpipi0_BESIII;
  ll3 += corrmeas.at("BESIII_Adk").logweight(corr);

  corr.ResizeTo(2);
  corr(0) = xi_x_BESIII;
  corr(1) = xi_y_BESIII;
  ll3 += corrmeas.at("BESIII_rDkpi_polar").logweight(corr);

  //https://arxiv.org/pdf/2208.06512.pdf
  //Observables 4:
  corr.ResizeTo(4);
  corr(0) = xcp;
  corr(1) = ycp;
  corr(2) = dx;
  corr(3) = dy;

  ll3 += corrmeas.at("LHCb_kspp_Au2022").logweight(corr);

  }

  return ll3;

}
// ---------------------------------------------------------

// ---------------------------------------------------------
double MixingModel::Calculate_other_observables(){

  double ll4;
  ll4 = 0.;

  //----------------------------------------------- Other Observables-------------------------------------------------------------------------
  //14. PDF: charm-kspipi-nocpv (UID13)
  // 2 Observables
  x_uid13 = x;
  y_uid13 = y;

  //19. PDF: charm-k3pi (UID18)
  // 1 Osservabile
  k3pi_uid18 = 0.25 * (x*x + y*y);

  //23. PDF: dkskpiRWS (UID22)
  // 1 Osservabile
  RD_kskpi_uid22 = (rD_kskpi * rD_kskpi - rD_kskpi * kD_kskpi * (y * cos(dD_kskpi) + x*sin(dD_kskpi))) / (1. - rD_kskpi*kD_kskpi*(y*cos(dD_kskpi) - x*sin(dD_kskpi)));
  //Relative Sign Problem

  //24. PDF: dkskpi (UID23)
  // 3 Observables
  RD_kskpi_uid23 = (rD_kskpi * rD_kskpi - rD_kskpi * kD_kskpi * (y * cos(dD_kskpi) + x*sin(dD_kskpi)) ) / (1. - rD_kskpi*kD_kskpi*(y*cos(dD_kskpi) - x*sin(dD_kskpi)) );
  dD_kskpi_uid23 = - dD_kskpi; // change sign of the phase convention
  kD_kskpi_uid23 = kD_kskpi;

  if(comb == 2){
    F_pipipipi_BESIII = F_pipipipi;
  }

  //----------------------------------------------- Contribution to the LogLikelihood -------------------------------------------------------------------------
  TVectorD corr(8);

  //14. PDF: charm-kspipi-nocpv (UID13)
  //Observables 2:
  corr.ResizeTo(2);
  corr(0) = x;
  corr(1) = y;

  ll4 += corrmeas.at("UID13").logweight(corr);

  //19. PDF: charm-k3pi (UID18)
  //Observables 1:
  ll4 += meas.at("UID18").logweight(k3pi_uid18);

  //23. PDF: dkskpiRWS (UID22)
  //Observables 1:
  ll4 += meas.at("UID22").logweight(RD_kskpi_uid22);

  //24. PDF: dkskpi (UID23)
  //Observables 3:
  corr.ResizeTo(3);
  corr(0) = RD_kskpi_uid23;
  corr(1) = dD_kskpi_uid23;
  corr(2) = kD_kskpi_uid23;

  ll4 += corrmeas.at("UID23").logweight(corr);

  if(comb == 2){
    ll4 += meas.at("Fpipipipi_BESIII").logweight(F_pipipipi_BESIII);
  }

  return ll4;

}
// ---------------------------------------------------------

// ---------------------------------------------------------
double MixingModel::Calculate_old_observables(){

  double llo;
  llo = 0.;

  delta_kpi =  - dD_kpi + M_PI; //because of the phase convention (see LHCb_combo_delta)
  delta_kpipi =  - dD_kpipi0 + M_PI;

  Rd = rD_kpi * rD_kpi;

  // 2nd Block
  Ag = (qop - 1. / qop) / 2. * y * cos(phi)- (qop + 1. / qop) / 2. * x * sin(phi);

  //6th Block
  rm =  (x * x + y * y) / 2.;

  // 8th Block P1
  double xpkpp = x * cos(delta_kpipi) + y * sin(delta_kpipi);
  double ypkpp = -x * sin(delta_kpipi) + y * cos(delta_kpipi);

  yp_kpp_plus = qop * (ypkpp * cos(phi) - xpkpp * sin(phi));
  xp_kpp_plus = qop * (xpkpp * cos(phi) + ypkpp * sin(phi));

  // 8th Block P2
  yp_kpp_minus = 1. / qop * (ypkpp * cos(phi) + xpkpp * sin(phi));
  xp_kpp_minus = (xpkpp * cos(phi) - ypkpp * sin(phi)) / qop;

  xp = x * cos(delta_kpi) + y * sin(delta_kpi);
  yp = -x * sin(delta_kpi) + y * cos(delta_kpi);

  // 10th and 12th Block
  yp_plus = qop * (yp * cos(phi) - xp * sin(phi));
  xp_plus = qop * (xp * cos(phi) + yp * sin(phi));
  xp_plus_sq = xp_plus * xp_plus;

  // 11th and 13th Block
  yp_minus = 1. / qop * (yp * cos(phi) + xp * sin(phi));
  xp_minus = (xp * cos(phi) - yp * sin(phi)) / qop;
  xp_minus_sq = xp_minus * xp_minus;

  //no direct CPV in tree decays
  // in the global combination It is allowed
  if(comb != 2){
    AD = 0.;
  }

  Rdkp = Rd * (1. + AD);

  TVectorD corr(4);

  // 1st Block
  llo += meas.at("ycp").logweight(ycp);

  // 2nd Block
  if(comb != 2){
    if (!noagamma && !combgamma_all){
      llo += meas.at("AGamma").logweight(Ag);
    }
  }

  // 3rd Block
  corr(0) = x;
  corr(1) = y;
  corr(2) = qop;
  corr(3) = phiKsBfact;

  llo += corrmeas.at("kpp_belle").logweight(corr);


  //Babar D0 -> Ks pi pi and D0 -> Ks K K
  // 5th Block
  corr.ResizeTo(2);
  corr(0) = x;
  corr(1) = y;

  llo += corrmeas.at("kppkk").logweight(corr);

  // 6th Block
  llo += meas.at("RM").logweight(rm);


  // 8th Block P1
  corr.ResizeTo(2);

  corr(0) =  yp_kpp_plus;
  corr(1) =  xp_kpp_plus;

  llo += corrmeas.at("kpp_babar_plus").logweight(corr);

  // 8th Block P2
  corr(0) =  yp_kpp_minus;
  corr(1) =  xp_kpp_minus;

  llo += corrmeas.at("kpp_babar_minus").logweight(corr);


  //CLEOc Psi results
  // 9th Block
  corr.ResizeTo(5);

  corr(0) = Rd;
  corr(1) = x*x;
  corr(2) = y;
  corr(3) = cos(delta_kpi);
  corr(4) = sin(delta_kpi);

  llo += corrmeas.at("cleoc").logweight(corr);

  // Babar Kpi & Belle Kpi
  //10th and 12th Block
  corr.ResizeTo(3);

  corr(0) = Rd;
  corr(1) = xp_plus_sq;
  corr(2) = yp_plus;

  llo += corrmeas.at("kpi_babar_plus").logweight(corr);
  llo += corrmeas.at("kpi_belle_plus").logweight(corr);

  // 11th and 13th Block
  corr(0) = AD;
  corr(1) = xp_minus_sq;
  corr(2) = yp_minus;

  llo += corrmeas.at("kpi_babar_minus").logweight(corr);
  llo += corrmeas.at("kpi_belle_minus").logweight(corr);

  // 14th Block
  corr(0) = Rd;
  corr(1) = (yp_plus + yp_minus) / 2.;
  corr(2) = (xp_plus_sq + xp_minus_sq + yp_plus * yp_plus + yp_minus * yp_minus) / 2. - corr(1) * corr(1);

  llo += corrmeas.at("kpi_cdf").logweight(corr);

  if(comb != 2){

    if (!combgamma_all) {

      //15th Block
      corr.ResizeTo(5);

      corr(0) = Rd;
      corr(1) = yp_plus;
      corr(2) = xp_plus_sq;
      corr(3) = yp_minus;
      corr(4) = xp_minus_sq;

      llo += corrmeas.at("kpi_lhcb").logweight(corr);

      // 7th Block
      llo += meas.at("RMKppp").logweight(rm);

      //4th Block
      corr.ResizeTo(4);

      corr(0) = xcpphiKsLHCb;
      corr(1) = ycpksppLHCb;
      corr(2) = dxphiKsLHCb;
      corr(3) = AgksppLHCb;

      llo += corrmeas.at("LHCb_kspp").logweight(corr);

    } else {

      // LHCb Combination Resultss
      corr.ResizeTo(6);

      corr(0) = x;
      corr(1) = y;
      corr(2) = qop;
      corr(3) = phi / d2r;
      corr(4) = sqrt(Rd);
      corr(5) = delta_kpi / d2r;

      llo += corrmeas.at("LHCb_combo").logweight(corr);
    }
    if (combgamma_delta) {

      // DKpi LHCb combination results
      corr.ResizeTo(2);

      corr(0) = sqrt(Rd);
      corr(1) = delta_kpi / d2r;

      llo += corrmeas.at("LHCb_combo_delta").logweight(corr);

    }
  }


  return llo;
}
// ---------------------------------------------------------

double MixingModel::LogLikelihood(const std::vector<double> &parameters)
{

  double ll = 0.;

  if(comb == 0){
    //--------------------------   assignment of the parameters variables to the extracted values --------------------------------------------

    g = parameters[0];
    x12 = parameters[1];
    y12 = parameters[2];

    // 1. PDF: glwads-dh-hh-dmix (UID0)
    r_dk = parameters[3];
    r_dpi = parameters[4];
    rD_kpi = parameters[5];
    d_dk = parameters[6];
    d_dpi = parameters[7];
    dD_kpi = parameters[8];

    //  2. PDF: glwads-dh-h3pi-dmix (UID1)
    rD_k3pi = parameters[9];
    dD_k3pi = parameters[10];
    kD_k3pi = parameters[11];
    F_pipipipi = parameters[12];

    // 3. PDF: glwads-dh-hhpi0-dmix (UID2)
    rD_kpipi0 = parameters[13];
    dD_kpipi0 = parameters[14];
    kD_kpipi0 = parameters[15];
    F_pipipi0 = parameters[16];
    F_kkpi0 = parameters[17];

    //5. PDF: glwads-dkdpi-kskpi-dmix (UID4)
    rD_kskpi = parameters[18];
    dD_kskpi = parameters[19];
    kD_kskpi = parameters[20];
    RBRdkdpi = parameters[21];

    // 6. PDF: glwads-dsth-hh-dmix (UID5)
    r_dstk = parameters[22];
    d_dstk = parameters[23];
    r_dstpi = parameters[24];
    d_dstpi = parameters[25];

    // 7. PDF: glwads-dkst-hh-h3pi-dmix-newvars (UID6)
    r_dkst = parameters[26];
    d_dkst = parameters[27];
    k_dkst = parameters[28];

    //  8. PDF: glwads-dkstz-hh-h3pi-dmix (UID7)
    r_dkstz = parameters[29];
    d_dkstz = parameters[30];
    k_dkstz = parameters[31];

    //10. PDF: glwads-dhpipi-hh-dmix (UID9)
    r_dkpipi = parameters[32];
    d_dkpipi = parameters[33];
    k_dkpipi = parameters[34];
    r_dpipipi = parameters[35];
    d_dpipipi = parameters[36];
    k_dpipipi = parameters[37];

    //11. PDF: dsk (UID10)
    l_dsk = parameters[38];
    d_dsk = parameters[39];
    phis = parameters[40];

    //12. PDF: dskpipi (UID11)
    l_dskpipi = parameters[41];
    d_dskpipi = parameters[42];
    k_dskpipi = parameters[43];

    //13. PDF: dmpi (UID12)
    l_dmpi = parameters[44];
    d_dmpi = parameters[45];
    beta = parameters[46];

    //15. PDF: charm-kspipi (UID14)
    PhiM12 = parameters[47];
    PhiG12 = parameters[48];

    //17. PDF: charm-kpi (UID16)
    AD = parameters[49];

    //18. PDF: charm-deltaacp-diff (UID17)
    Delta_totau = parameters[50];
    DeltaAcp = parameters[51];

    //General parameters
    phi12 = remainder(-PhiG12 + PhiM12, 2. * M_PI);
    x = x12;
    y = y12;
    qop = 1. + x12 * y12 * sin(phi12) / (x12 * x12 + y12 * y12);
    phi = atan(- (x12 * x12 * sin(2. * PhiM12) + y12 * y12 * sin(2. * PhiG12)) /
    (x12 * x12 * cos(2. * PhiM12) + y12 * y12 * cos(2. * PhiG12)) ) / 2.;
    d = (1. - qop * qop) / (1. + qop * qop);

    xcp = 0.5 * ( x * cos(phi) * (qop + 1./qop) + y * sin(phi) * (qop - 1./qop)  );
    ycp = 0.5 * ( y * cos(phi) * (qop + 1./qop) - x * sin(phi) * (qop - 1./qop));
    dx = 0.5 * ( x*cos(phi)*(qop - 1./qop) + y*sin(phi)*(qop + 1./qop) );
    dy = 0.5 * ( y*cos(phi)*(qop - 1./qop) - x*sin(phi)*(qop + 1./qop)  );

    //------------------------------------------ Calculating the observables and their contributions to the Log-Likelihood ---------------------------------

    //---------------------------------------------------- Time integrated observables ----------------------------------------------------
    ll += Calculate_time_integrated_observables();
    //------------------------------------------------------  Time Dependent B observables  ----------------------------------------------------
    ll += Calculate_time_dependent_Bobservables();
    //------------------------------------------------------  Time Dependent D observables  ----------------------------------------------------
    ll += Calculate_time_dependent_Dobservables();
    //----------------------------------------------- Other Observables-------------------------------------------------------------------------
    ll += Calculate_other_observables();

    //---------------------------------- Fill the Histograms ---------------------

    // Time Integrated B decay chain
    obs["g"] = g * r2d;
    obs["x12"] = x12;
    obs["y12"] = y12;
    obs["r_dk"] = r_dk;
    obs["r_dpi"] = r_dpi;
    obs["rD_kpi"] = rD_kpi;
    obs["d_dk"] = d_dk * r2d;
    obs["d_dpi"] = d_dpi * r2d;
    obs["dD_kpi"] = dD_kpi * r2d;
    obs["dD_kpi_LHCb"] = (2*M_PI - dD_kpi) * r2d;
    obs["rD_k3pi"] = rD_k3pi;
    obs["dD_k3pi"] = dD_k3pi * r2d;
    obs["dD_k3pi_LHCb"] = (2*M_PI - dD_k3pi) * r2d;
    obs["kD_k3pi"] = kD_k3pi;
    obs["F_pipipipi"] = F_pipipipi;
    obs["rD_kpipi0"] = rD_kpipi0;
    obs["dD_kpipi0"] = dD_kpipi0 * r2d;
    obs["dD_kpipi0_LHCb"] = (2*M_PI - dD_kpipi0) * r2d;
    obs["kD_kpipi0"] = kD_kpipi0;
    obs["F_pipipi0"] = F_pipipi0;
    obs["F_kkpi0"] = F_kkpi0;
    obs["rD_kskpi"] = rD_kskpi;
    obs["dD_kskpi"] = dD_kskpi * r2d;
    obs["dD_kskpi_LHCb"] = -dD_kskpi * r2d;
    obs["kD_kskpi"] = kD_kskpi;
    obs["RBRdkdpi"] = RBRdkdpi;
    obs["r_dstk"] = r_dstk;
    obs["d_dstk"] = d_dstk * r2d;
    obs["r_dstpi"] = r_dstpi;
    obs["d_dstpi"] = d_dstpi * r2d;
    obs["r_dkst"] = r_dkst;
    obs["d_dkst"] = d_dkst * r2d;
    obs["k_dkst"] = k_dkst;
    obs["r_dkstz"] = r_dkstz;
    obs["d_dkstz"] = d_dkstz * r2d;
    obs["k_dkstz"] = k_dkstz;
    obs["r_dkpipi"] = r_dkpipi;
    obs["d_dkpipi"] = d_dkpipi * r2d;
    obs["k_dkpipi"] = k_dkpipi;
    obs["r_dpipipi"] = r_dpipipi;
    obs["d_dpipipi"] = d_dpipipi * r2d;
    obs["k_dpipipi"] = k_dpipipi;

    // Time dependent B decay chain

    obs["l_dsk"] = l_dsk;
    obs["d_dsk"] = d_dsk;
    obs["phis"] = phis * r2d;
    obs["beta_s"] = - 0.5 * phis;
    obs["l_dskpipi"] = l_dskpipi;
    obs["d_dskpipi"] = d_dskpipi;
    obs["k_dskpipi"] = k_dskpipi;
    obs["l_dmpi"] = l_dmpi;
    obs["d_dmpi"] = d_dmpi * r2d;
    obs["beta"] = beta * r2d;

    //Time dependent D decay and mixing
    obs["PhiM12"] = PhiM12 * r2d;
    obs["PhiG12"] = PhiG12 * r2d;
    obs["phipphig12"] = remainder(PhiG12 + phi, 2. * M_PI) * r2d;
    obs["phimphig12"] = remainder(-PhiG12 + phi, 2. * M_PI) * r2d;
    obs["AD"] = AD;
    obs["Delta_totau"] = Delta_totau;
    obs["DeltaAcp"] = DeltaAcp;
    obs["qopm1"] = qop -1;
    obs["qop"] = qop;
    obs["phi"] = phi * r2d;
    obs["phi12"] = phi12 * r2d;
    obs["delta"] = d;
    obs["x"] = x;
    obs["y"] = y;

  }
  else if(comb == 1){

    x12 = parameters[0];
    y12 = parameters[1];
    PhiM12 = parameters[2];
    dD_kpi = parameters[3];
    dD_kpipi0 = parameters[4];
    rD_kpi = parameters[5];
    PhiG12 = parameters[6];
    if (ckmcorr){
      RCP = parameters[7];
    }

    //General parameters
    phi12 = remainder(-PhiG12 + PhiM12, 2. * M_PI);
    x = x12;
    y = y12;
    qop = 1. + x12 * y12 * sin(phi12) / (x12 * x12 + y12 * y12);
    phi = atan(- (x12 * x12 * sin(2. * PhiM12) + y12 * y12 * sin(2. * PhiG12)) /
    (x12 * x12 * cos(2. * PhiM12) + y12 * y12 * cos(2. * PhiG12)) ) / 2.;
    d = (1. - qop * qop) / (1. + qop * qop);

    phiKsBfact = phi;
    phiKsLHCb = phi;
    if (epsK){
    phiKsBfact += -2 * epsI;
  }
    if (ckmcorr) {
      //    std::cout << RCP << std::endl;
      phiKsBfact -= RCP;
      phiKsLHCb -= RCP;
    }

    ycpksppLHCb = (qop + 1. / qop) / 2. * y * cos(phiKsLHCb) - (qop - 1. / qop) / 2. * x * sin(phiKsLHCb);
    AgksppLHCb = (qop - 1. / qop) / 2. * y * cos(phiKsLHCb)- (qop + 1. / qop) / 2. * x * sin(phiKsLHCb);
    xcpphiKsLHCb = 0.5 * ( x * cos(phiKsLHCb) * (qop + 1./qop) + y * sin(phiKsLHCb) * (qop - 1./qop)  );
    dxphiKsLHCb = 0.5 * ( x*cos(phiKsLHCb)*(qop - 1./qop) + y*sin(phiKsLHCb)*(qop + 1./qop) );
    ycpksppLHCb = (qop + 1. / qop) / 2. * y * cos(phiKsLHCb) - (qop - 1. / qop) / 2. * x * sin(phiKsLHCb);
    am = (qop * qop - 1. / (qop * qop)) / (qop * qop + 1. / (qop * qop));

    xcp = 0.5 * ( x * cos(phi) * (qop + 1./qop) + y * sin(phi) * (qop - 1./qop)  );
    ycp = 0.5 * ( y * cos(phi) * (qop + 1./qop) - x * sin(phi) * (qop - 1./qop));
    dx = 0.5 * ( x*cos(phi)*(qop - 1./qop) + y*sin(phi)*(qop + 1./qop) );
    dy = 0.5 * ( y*cos(phi)*(qop - 1./qop) - x*sin(phi)*(qop + 1./qop)  );


    //----------------------------------Contributions to the LL-------------------------------------------------------------------------
    ll += Calculate_old_observables();

    //----------------------------------Fill the Histograms-------------------------------------------------------------------------
    obs["y"] = y;
    obs["x"] = x;
    obs["dD_kpi"] = dD_kpi * r2d;
    obs["dD_kpipi0"] = dD_kpipi0 * r2d;
    obs["delta_kp"] = delta_kpi * r2d;
    obs["delta_kpp"] = delta_kpipi * r2d;
    obs["phi"] = phi * r2d;
    obs["phiG12"] = PhiG12 * r2d;
    obs["phipphig12"] = remainder(PhiG12 + phi, 2. * M_PI) * r2d;
    obs["phimphig12"] = remainder(-PhiG12 + phi, 2. * M_PI) * r2d;
    obs["phiM12"] = PhiM12 * r2d;
    obs["phi12"] = phi12 * r2d;
    obs["x12"] = x12;
    obs["y12"] = y12;
    obs["qopm1"] = qop - 1.;
    obs["delta"] = d;
    obs["RD"] = Rd;
    obs["AD"] = AD;
    obs["RDkp"] = Rdkp;
    obs["rDkp"] = sqrt(Rdkp);
    obs["Rm"] = rm;
    obs["Am"] = am;
    obs["AGamma"] = Ag;
    obs["ycp"] = ycp;
    obs["dx"] = dx;
    obs["dxphiKsLHCb"] = dxphiKsLHCb;
    obs["xcp"] = xcp;
    obs["xcpphiKsLHCb"] = xcpphiKsLHCb;
    obs["ypp"] = yp_plus;
    obs["ypm"] = yp_minus;
    obs["xpp"] = xp_plus;
    obs["xpm"] = xp_minus;
    obs["xpp2"] = xp_plus_sq;
    obs["xpm2"] = xp_minus_sq;
    obs["ypp_kpp"] = yp_kpp_plus;
    obs["ypm_kpp"] = yp_kpp_minus;
    obs["xpp_kpp"] = xp_kpp_plus;
    obs["xpm_kpp"] = xp_kpp_minus;
    obs["phiKsLHCb"] = phiKsLHCb * r2d;
    obs["RCP"] = RCP;

  }
  else if(comb == 2){

    //Variabili globali
    g = parameters[0];
    x12 = parameters[1];
    y12 = parameters[2];

    // 1. PDF: glwads-dh-hh-dmix (UID0)
    r_dk = parameters[3];
    r_dpi = parameters[4];
    rD_kpi = parameters[5];
    d_dk = parameters[6];
    d_dpi = parameters[7];
    dD_kpi = parameters[8];

    //  2. PDF: glwads-dh-h3pi-dmix (UID1)
    rD_k3pi = parameters[9];
    dD_k3pi = parameters[10];
    kD_k3pi = parameters[11];
    F_pipipipi = parameters[12];

    // 3. PDF: glwads-dh-hhpi0-dmix (UID2)
    rD_kpipi0 = parameters[13];
    dD_kpipi0 = parameters[14];
    kD_kpipi0 = parameters[15];
    F_pipipi0 = parameters[16];
    F_kkpi0 = parameters[17];

    //5. PDF: glwads-dkdpi-kskpi-dmix (UID4)
    rD_kskpi = parameters[18];
    dD_kskpi = parameters[19];
    kD_kskpi = parameters[20];
    RBRdkdpi = parameters[21];

    // 6. PDF: glwads-dsth-hh-dmix (UID5)
    r_dstk = parameters[22];
    d_dstk = parameters[23];
    r_dstpi = parameters[24];
    d_dstpi = parameters[25];

    // 7. PDF: glwads-dkst-hh-h3pi-dmix-newvars (UID6)
    r_dkst = parameters[26];
    d_dkst = parameters[27];
    k_dkst = parameters[28];

    //  8. PDF: glwads-dkstz-hh-h3pi-dmix (UID7)
    r_dkstz = parameters[29];
    d_dkstz = parameters[30];
    k_dkstz = parameters[31];

    //10. PDF: glwads-dhpipi-hh-dmix (UID9)
    r_dkpipi = parameters[32];
    d_dkpipi = parameters[33];
    k_dkpipi = parameters[34];
    r_dpipipi = parameters[35];
    d_dpipipi = parameters[36];
    k_dpipipi = parameters[37];

    //11. PDF: dsk (UID10)
    l_dsk = parameters[38];
    d_dsk = parameters[39];
    phis = parameters[40];

    //12. PDF: dskpipi (UID11)
    l_dskpipi = parameters[41];
    d_dskpipi = parameters[42];
    k_dskpipi = parameters[43];

    //13. PDF: dmpi (UID12)
    l_dmpi = parameters[44];
    d_dmpi = parameters[45];
    beta = parameters[46];

    //15. PDF: charm-kspipi (UID14)
    PhiM12 = parameters[47];
    PhiG12 = parameters[48];

    //17. PDF: charm-kpi (UID16)
    AD = parameters[49];

    //18. PDF: charm-deltaacp-diff (UID17)
    Delta_totau = parameters[50];
    DeltaAcp = parameters[51];

    if (ckmcorr){
      RCP = parameters[52];
    }

    //General parameters
    phi12 = remainder(-PhiG12 + PhiM12, 2. * M_PI);
    x = x12;
    y = y12;
    qop = 1. + x12 * y12 * sin(phi12) / (x12 * x12 + y12 * y12);
    phi = atan(- (x12 * x12 * sin(2. * PhiM12) + y12 * y12 * sin(2. * PhiG12)) /
    (x12 * x12 * cos(2. * PhiM12) + y12 * y12 * cos(2. * PhiG12)) ) / 2.;
    d = (1. - qop * qop) / (1. + qop * qop);

    xcp = 0.5 * ( x * cos(phi) * (qop + 1./qop) + y * sin(phi) * (qop - 1./qop)  );
    ycp = 0.5 * ( y * cos(phi) * (qop + 1./qop) - x * sin(phi) * (qop - 1./qop));
    dx = 0.5 * ( x*cos(phi)*(qop - 1./qop) + y*sin(phi)*(qop + 1./qop) );
    dy = 0.5 * ( y*cos(phi)*(qop - 1./qop) - x*sin(phi)*(qop + 1./qop)  );

    phiKsBfact = phi;
    phiKsLHCb = phi;
    if (epsK){
    phiKsBfact += -2 * epsI;
  }
    if (ckmcorr) {
      //    std::cout << RCP << std::endl;
      phiKsBfact -= RCP;
      phiKsLHCb -= RCP;
    }

    ycpksppLHCb = (qop + 1. / qop) / 2. * y * cos(phiKsLHCb) - (qop - 1. / qop) / 2. * x * sin(phiKsLHCb);
    AgksppLHCb = (qop - 1. / qop) / 2. * y * cos(phiKsLHCb)- (qop + 1. / qop) / 2. * x * sin(phiKsLHCb);
    xcpphiKsLHCb = 0.5 * ( x * cos(phiKsLHCb) * (qop + 1./qop) + y * sin(phiKsLHCb) * (qop - 1./qop)  );
    dxphiKsLHCb = 0.5 * ( x*cos(phiKsLHCb)*(qop - 1./qop) + y*sin(phiKsLHCb)*(qop + 1./qop) );
    ycpksppLHCb = (qop + 1. / qop) / 2. * y * cos(phiKsLHCb) - (qop - 1. / qop) / 2. * x * sin(phiKsLHCb);
    am = (qop * qop - 1. / (qop * qop)) / (qop * qop + 1. / (qop * qop));

    //---------------------------------------------------- Time integrated observables ----------------------------------------------------
    ll += Calculate_time_integrated_observables();
    //------------------------------------------------------  Time Dependent B observables  ----------------------------------------------------
    ll += Calculate_time_dependent_Bobservables();
    //------------------------------------------------------  Time Dependent D observables  ----------------------------------------------------
    ll += Calculate_time_dependent_Dobservables();
    //----------------------------------------------- Other Observables-------------------------------------------------------------------------
    ll += Calculate_other_observables();
    //----------------------------------------------- Old Observables-------------------------------------------------------------------------
    ll += Calculate_old_observables();

    //---------------------------------- Fill the Histograms ---------------------

    // Time Integrated B decay chain
    obs["g"] = g * r2d;
    obs["x12"] = x12;
    obs["y12"] = y12;
    obs["r_dk"] = r_dk;
    obs["r_dpi"] = r_dpi;
    obs["rD_kpi"] = rD_kpi;
    obs["d_dk"] = d_dk * r2d;
    obs["d_dpi"] = d_dpi * r2d;
    obs["dD_kpi"] = dD_kpi * r2d;
    obs["dD_kpi_LHCb"] = (2*M_PI - dD_kpi) * r2d;
    obs["rD_k3pi"] = rD_k3pi;
    obs["dD_k3pi"] = dD_k3pi * r2d;
    obs["dD_k3pi_LHCb"] = (2*M_PI - dD_k3pi) * r2d;
    obs["kD_k3pi"] = kD_k3pi;
    obs["F_pipipipi"] = F_pipipipi;
    obs["rD_kpipi0"] = rD_kpipi0;
    obs["dD_kpipi0"] = dD_kpipi0 * r2d;
    obs["dD_kpipi0_LHCb"] = (2*M_PI - dD_kpipi0) * r2d;
    obs["kD_kpipi0"] = kD_kpipi0;
    obs["F_pipipi0"] = F_pipipi0;
    obs["F_kkpi0"] = F_kkpi0;
    obs["rD_kskpi"] = rD_kskpi;
    obs["dD_kskpi"] = dD_kskpi * r2d;
    obs["dD_kskpi_LHCb"] = -dD_kskpi * r2d;
    obs["kD_kskpi"] = kD_kskpi;
    obs["RBRdkdpi"] = RBRdkdpi;
    obs["r_dstk"] = r_dstk;
    obs["d_dstk"] = d_dstk * r2d;
    obs["r_dstpi"] = r_dstpi;
    obs["d_dstpi"] = d_dstpi * r2d;
    obs["r_dkst"] = r_dkst;
    obs["d_dkst"] = d_dkst * r2d;
    obs["k_dkst"] = k_dkst;
    obs["r_dkstz"] = r_dkstz;
    obs["d_dkstz"] = d_dkstz * r2d;
    obs["k_dkstz"] = k_dkstz;
    obs["r_dkpipi"] = r_dkpipi;
    obs["d_dkpipi"] = d_dkpipi * r2d;
    obs["k_dkpipi"] = k_dkpipi;
    obs["r_dpipipi"] = r_dpipipi;
    obs["d_dpipipi"] = d_dpipipi * r2d;
    obs["k_dpipipi"] = k_dpipipi;

    // Time dependent B decay chain

    obs["l_dsk"] = l_dsk;
    obs["d_dsk"] = d_dsk;
    obs["phis"] = phis * r2d;
    obs["beta_s"] = - 0.5 * phis;
    obs["l_dskpipi"] = l_dskpipi;
    obs["d_dskpipi"] = d_dskpipi;
    obs["k_dskpipi"] = k_dskpipi;
    obs["l_dmpi"] = l_dmpi;
    obs["d_dmpi"] = d_dmpi * r2d;
    obs["beta"] = beta * r2d;

    //Time dependent D decay and mixing
    obs["PhiM12"] = PhiM12 * r2d;
    obs["PhiG12"] = PhiG12 * r2d;
    obs["phipphig12"] = remainder(PhiG12 + phi, 2. * M_PI) * r2d;
    obs["phimphig12"] = remainder(-PhiG12 + phi, 2. * M_PI) * r2d;
    obs["AD"] = AD;
    obs["Delta_totau"] = Delta_totau;
    obs["DeltaAcp"] = DeltaAcp;
    obs["qopm1"] = qop -1;
    obs["qop"] = qop;
    obs["phi"] = phi * r2d;
    obs["phi12"] = phi12 * r2d;
    obs["delta"] = d;
    obs["delta_kp"] = delta_kpi * r2d;
    obs["delta_kpp"] = delta_kpipi * r2d;

    obs["phiKsLHCb"] = phiKsLHCb * r2d;
    obs["RCP"] = RCP;
    obs["x"] = x;
    obs["y"] = y;

  }

  return ll;

}
// ---------------------------------------------------------

void MixingModel::MCMCUserIterationInterface()
{
  std::vector<double> pars;
  for (unsigned int i = 0; i < fMCMCNChains; ++i) {
    //        std::cout << "interface: chain " << i << ", ";
    pars = fMCMCStates.at(i).parameters;
    //    std::cout << pars.size() << std::endl;
    LogLikelihood(pars);
    histos.fillh1d();
    histos.fillh2d();
  }
}

void MixingModel::PrintHistogram()
{
  histos.write();
}
