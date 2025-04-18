#include "MixingModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCGaussianPrior.h>

#include <TRandom3.h>
using namespace std;

// ---------------------------------------------------------

MixingModel::MixingModel(std::vector<string> nParam, int combination) : BCModel(), histos(obs)
{

  //------------------------------------------ Setting the auxiliary variables --------------------------------------------------------------------------

  TH1::SetDefaultBufferSize(1000000);
  comb = combination; // set the combination type
  r2d = 180. / M_PI;  // to go form radiants to degrees
  d2r = M_PI / 180.;  // degrees to radiants
  tau = 4.1e-1;       // ps D lifetime
  nVarab = nParam;    // copy the names of the variables of interest to print the histograms

  //------------------------------------------ Inserting the measurements --------------------------------------------------------------------------

  if (comb == 0)
  {                             // Charged B measurements
    Add_ChargedB_meas();        // Charged B decay chains
    Add_time_dependent_Dmeas(); // time dependent D decay and mixing
    Add_other_meas();           // Other useful inputs
    Add_old_meas();             // old code measurements
  }
  else if (comb == 1)
  {                             // Neutral Bd measurements
    Add_NeutralBd_meas();       // time dependent Bd decay and mixing
    Add_time_dependent_Dmeas(); // time dependent D decay and mixing
    Add_other_meas();           // Other useful inputs
    Add_old_meas();             // old code measurements
  }
  else if (comb == 2)
  {                             // Neutral Bs measurements
    Add_NeutralBs_meas();       // time dependent Bs decay and mixing
    Add_time_dependent_Dmeas(); // time dependent D decay and mixing
    Add_other_meas();           // Other useful inputs
    Add_old_meas();             // old code measurements
  }
  else if (comb == 3)
  {                             // All modes
    Add_ChargedB_meas();        // Charged B decay chains
    Add_NeutralBd_meas();       // time dependent Bd decay and mixing
    Add_NeutralBs_meas();       // time dependent Bs decay and mixing
    Add_time_dependent_Dmeas(); // time dependent D decay and mixing
    Add_other_meas();           // Other useful inputs
    Add_old_meas();             // old code measurements
  }
  else if(comb == 4){ // Only charm modes
    Add_time_dependent_Dmeas(); // time dependent D decay and mixing
    Add_other_meas();           // Other useful inputs
    Add_old_meas();             // old code measurements
  }

  //------------------------------------------- Defining the Histograms --------------------------------------------------------------

  DefineHistograms();

  //------------------------------------------- Defining the parameters ---------------------------------------------------------------------------

  DefineParameters();

};

// ---------------------------------------------------------
MixingModel::~MixingModel() { // default destructor
};
// ---------------------------------------------------------

void MixingModel::SymmetrizeUpperTriangularMatrix(TMatrixDSym& Corr) // Function to symmetrize an upper triangular matrix
{
  for (int i = 0; i < Corr.GetNrows(); ++i)
  {
    for (int j = i + 1; j < Corr.GetNcols(); ++j)
    {
      Corr(j, i) = Corr(i, j);
    }
  }
}

// ---------------------------------------------------------
void MixingModel::Add_ChargedB_meas()
{

  vector<dato> CorrData;
  TMatrixDSym Corr, Corr2, Corr3;

  //-------------------------------------------------  Babar measurements  -------------------------------------------------------------------------

  // B -> DK, D -> KK, D -> pipi normalized to D -> Kpi
  // https://journals.aps.org/prd/pdf/10.1103/PhysRevD.82.072004
  // 4 Observables :
  CorrData.clear();
  CorrData.push_back(dato(0.25, 0.06, 0.02));   // Acp+
  CorrData.push_back(dato(-0.09, 0.07, 0.02));   // Acp-
  CorrData.push_back(dato(1.18, 0.09, 0.05));   // RCP+
  CorrData.push_back(dato(1.07, 0.08, 0.04));   // RCP-

  // Correlation Matrix (stat):
  Corr.ResizeTo(4, 4);
  Corr(0,0)=1.;
  Corr(0,1)=0.;
  Corr(0,2)=-0.08;
  Corr(0,3)=0.;
  Corr(1,1)=1.;
  Corr(1,2)=0.;
  Corr(1,3)=0.04;
  Corr(2,2)=1.;
  Corr(2,3)=0.09;
  Corr(3,3)=1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(4, 4);
  Corr2(0,0)=1.;
  Corr2(0,1)=0.56;
  Corr2(0,2)=-0.06;
  Corr2(0,3)=0.;
  Corr2(1,1)=1.;
  Corr2(1,2)=0.;
  Corr2(1,3)=0.;
  Corr2(2,2)=1.;
  Corr2(2,3)=0.12;
  Corr2(3,3)=1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Babar_PRD82_072004", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  //  https://arxiv.org/pdf/0807.2408
  // B -> DstarK
  // D -> KK, pipi + fcp-
  // 4 Observables :
  meas.insert(pair<string, dato>("Babar_0807.2408_ACPp", dato(-0.11, 0.09, 0.01))); // ACP+
  meas.insert(pair<string, dato>("Babar_0807.2408_ACPm", dato(0.06, 0.1, 0.02))); // ACP-
  meas.insert(pair<string, dato>("Babar_0807.2408_RCPp", dato(1.31, 0.13, 0.03))); // RCP+
  meas.insert(pair<string, dato>("Babar_0807.2408_RCPm", dato(1.10, 0.12, 0.04))); // RCP-


  // https://journals.aps.org/prd/pdf/10.1103/PhysRevD.80.092001
  // B -> DKstar
  // D -> KK, pipi, fcp-
  meas.insert(pair<string, dato>("Babar_PRD80_092001_ACPp", dato(0.09, 0.13, 0.06))); // ACP+
  meas.insert(pair<string, dato>("Babar_PRD80_092001_ACPm", dato(-0.23, 0.21, 0.07))); // ACP-
  meas.insert(pair<string, dato>("Babar_PRD80_092001_RCPp", dato(2.17, 0.35, 0.09))); // RCP+
  meas.insert(pair<string, dato>("Babar_PRD80_092001_RCPm", dato(1.03, 0.27, 0.13))); // RCP-


  // https://arxiv.org/pdf/hep-ex/0703037
  // B -> DK
  // GLW D -> pi+pi-pi0
  meas.insert(pair<string, dato>("Babar_0703037", dato(-0.02, 0.15, 0.03))); // ACP
  CorrData.clear();
  CorrData.push_back(dato(0.75, 0.11, 0.04));   // rho+
  CorrData.push_back(dato( 147 / 180. * M_PI, 23. / 180. * M_PI, 1./180 * M_PI));   // theta+
  CorrData.push_back(dato(0.72, 0.11, 0.04));   // rho-
  CorrData.push_back(dato(173. / 180. * M_PI, 42. / 180. * M_PI, 2. / 180. * M_PI));   // theta-

  // Correlation Matrix (stat):
  Corr.ResizeTo(4, 4);
  Corr(0,0)=1.;
  Corr(0,1)=0.;
  Corr(0,2)=0.;
  Corr(0,3)=0.;
  Corr(1,1)=1.;
  Corr(1,2)=0.;
  Corr(1,3)=0.;
  Corr(2,2)=1.;
  Corr(2,3)=0.;
  Corr(3,3)=1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(4, 4);
  Corr2(0,0)=1.;
  Corr2(0,1)=0.;
  Corr2(0,2)=0.14;
  Corr2(0,3)=0.;
  Corr2(1,1)=1.;
  Corr2(1,2)=0.;
  Corr2(1,3)=0.;
  Corr2(2,2)=1.;
  Corr2(2,3)=0.;
  Corr2(3,3)=1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Babar_0703037_rhotheta", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // https://arxiv.org/pdf/1006.4241
  // B -> DK
  // D -> Kpi
  meas.insert(pair<string, dato>("Babar_1006.4241_RDK", dato(1.1e-2, 0.5e-2, 0.2e-2))); // R_DK_Kpi
  meas.insert(pair<string, dato>("Babar_1006.4241_ADK", dato(-0.88, 0.47, 0.14))); // A_DK_Kpi
  // B -> [Dpi0]_Dstar K
  // D -> Kpi
  meas.insert(pair<string, dato>("Babar_1006.4241_RDstarKpi0", dato(1.8e-2, 0.9e-2, 0.4e-2))); // R_[Dpi0]_DstarK_Kpi
  meas.insert(pair<string, dato>("Babar_1006.4241_ADstarKpi0", dato(0.77, 0.35, 0.12))); // A_[Dpi0]_DstarK_Kpi
  // B -> [Dg]_Dstar K
  // D -> Kpi
  meas.insert(pair<string, dato>("Babar_1006.4241_RDstarKg", dato(1.3e-2, 1.4e-2, 0.8e-2))); // R_[Dg]_DstarK_Kpi
  meas.insert(pair<string, dato>("Babar_1006.4241_ADstarKg", dato(0.28, 0.94, 0.33))); // A_[Dg]_DstarK_Kpi
  // B -> Dpi
  // D -> Kpi
  meas.insert(pair<string, dato>("Babar_1006.4241_RDpi", dato(3.3e-3, 0.6e-3, 0.4e-3))); // R_Dpi_Kpi
  meas.insert(pair<string, dato>("Babar_1006.4241_ADpi", dato(0.03, 0.17, 0.04))); // A_Dpi_Kpi
  // B -> [Dpi0]_Dstarpi
  // D -> Kpi
  meas.insert(pair<string, dato>("Babar_1006.4241_RDstarpipi0", dato(3.2e-3, 0.9e-3, 0.8e-3))); // R_[Dpi0]_Dstarpi_Kpi
  meas.insert(pair<string, dato>("Babar_1006.4241_ADstarpipi0", dato(-0.09, 0.27, 0.05))); // A_[Dpi0]_Dstarpi_Kpi
  // B -> [Dpig]_Dstarpi
  // D -> Kpi
  meas.insert(pair<string, dato>("Babar_1006.4241_RDstarpig", dato(2.7e-3, 1.4e-3, 2.2e-3))); // R_[Dg]_Dstarpi_Kpi
  meas.insert(pair<string, dato>("Babar_1006.4241_ADstarpig", dato(-0.65, 0.55, 0.22))); // A_[Dg]_Dstarpi_Kpi


  // https://arxiv.org/pdf/1104.4472
  // B -> DK
  // D -> Kpipi0
  meas.insert(pair<string, dato>("Babar_1104.4472_Rp_DK_Kpipi0", dato(4.5e-3, 11e-3, 2.5e-3))); // Rp_DK_Kpipi0
  meas.insert(pair<string, dato>("Babar_1104.4472_Rm_DK_Kpipi0", dato(12e-3, 11e-3, 3e-3))); // Rm_DK_Kpipi0


  //https://arxiv.org/pdf/0909.3981
  // B -> DK
  // D -> Kpi
  meas.insert(pair<string, dato>("Babar_0909.3981_RADS", dato(0.066, 0.031, 0.010))); // RADS_DKstar_Kpi
  meas.insert(pair<string, dato>("Babar_0909.3981_Asup", dato(-0.34, 0.43, 0.16))); // Asup_DKstar_Kpi


  // https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.105.121801
  // B -> D(star)K(star) BPGGSZ
  // D -> K0spipi + D -> K0sKK
  CorrData.clear();
  CorrData.push_back(dato(6.e-2, 3.9e-2, 0.7e-2, 0.6e-2));   // x_DK-
  CorrData.push_back(dato(6.2e-2, 4.5e-2, 0.4e-2, 0.6e-2));   // y_DK-
  CorrData.push_back(dato(-10.3e-2, 3.7e-2, 0.6e-2, 0.7e-2));   // x_DK+
  CorrData.push_back(dato(-2.1e-2, 4.8e-2, 0.4e-2, 0.9e-2));   // y_DK+
  CorrData.push_back(dato(-10.4e-2, 5.1e-2, 1.9e-2, 0.2e-2));   // x_DstarK-
  CorrData.push_back(dato(-5.2e-2, 6.3e-2, 0.9e-2, 0.7e-2));   // y_DstarK-
  CorrData.push_back(dato(14.7e-2, 5.3e-2, 1.7e-2, 0.3e-2));   // x_DstarK+
  CorrData.push_back(dato(-3.2e-2, 7.7e-2, 0.8e-2, 0.6e-2));   // y_DstarK+
  CorrData.push_back(dato(7.5e-2, 9.6e-2, 2.9e-2, 0.7e-2));   // x_DKstar-
  CorrData.push_back(dato(12.7e-2, 9.5e-2, 2.7e-2, 0.6e-2));   // y_DKstar-
  CorrData.push_back(dato(-15.1e-2, 8.3e-2, 2.9e-2, 0.6e-2));   // x_DKstar+
  CorrData.push_back(dato(4.5e-2, 10.6e-2, 3.6e-2, 0.8e-2));   // y_DKstar+

  // Correlation Matrix (stat):
  Corr.ResizeTo(12, 12);
  Corr(0,0)=1.;
  Corr(0,1)=-0.02;
  Corr(0,2)=0.;
  Corr(0,3)=0.;
  Corr(0,4)=0.;
  Corr(0,5)=0.;
  Corr(0,6)=0.;
  Corr(0,7)=0.;
  Corr(0,8)=0.;
  Corr(0,9)=0.;
  Corr(0,10)=0.;
  Corr(0,11)=0.;
  Corr(1,1)=1.;
  Corr(1,2)=0.;
  Corr(1,3)=0.;
  Corr(1,4)=0.;
  Corr(1,5)=0.;
  Corr(1,6)=0.;
  Corr(1,7)=0.;
  Corr(1,8)=0.;
  Corr(1,9)=0.;
  Corr(1,10)=0.;
  Corr(1,11)=0.;
  Corr(2,2)=1.;
  Corr(2,3)=0.06;
  Corr(2,4)=0.;
  Corr(2,5)=0.;
  Corr(2,6)=0.;
  Corr(2,7)=0.;
  Corr(2,8)=0.;
  Corr(2,9)=0.;
  Corr(2,10)=0.;
  Corr(2,11)=0.;
  Corr(3,3)=1.;
  Corr(3,4)=0.;
  Corr(3,5)=0.;
  Corr(3,6)=0.;
  Corr(3,7)=0.;
  Corr(3,8)=0.;
  Corr(3,9)=0.;
  Corr(3,10)=0.;
  Corr(3,11)=0.;
  Corr(4,4)=1.;
  Corr(4,5)=-0.01;
  Corr(4,6)=0.;
  Corr(4,7)=0.;
  Corr(4,8)=0.;
  Corr(4,9)=0.;
  Corr(4,10)=0.;
  Corr(4,11)=0.;
  Corr(5,5)=1.;
  Corr(5,6)=0.;
  Corr(5,7)=0.;
  Corr(5,8)=0.;
  Corr(5,9)=0.;
  Corr(5,10)=0.;
  Corr(5,11)=0.;
  Corr(6,6)=1.;
  Corr(6,7)=-0.03;
  Corr(6,8)=0.;
  Corr(6,9)=0.;
  Corr(6,10)=0.;
  Corr(6,11)=0.;
  Corr(7,7)=1.;
  Corr(7,8)=0.;
  Corr(7,9)=0.;
  Corr(7,10)=0.;
  Corr(7,11)=0.;
  Corr(8,8)=1.;
  Corr(8,9)=-0.2;
  Corr(8,10)=0.;
  Corr(8,11)=0.;
  Corr(9,9)=1.;
  Corr(9,10)=0.;
  Corr(9,11)=0.;
  Corr(10,10)=1.;
  Corr(10,11)=0.1;
  Corr(11,11)=1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(12, 12);
  Corr2(0,0)=1.;
  Corr2(0,1)=-0.23;
  Corr2(0,2)=-0.02;
  Corr2(0,3)=-0.05;
  Corr2(0,4)=-0.32;
  Corr2(0,5)=-0.15;
  Corr2(0,6)=-0.12;
  Corr2(0,7)=0.12;
  Corr2(0,8)=-0.19;
  Corr2(0,9)=0.03;
  Corr2(0,10)=-0.20;
  Corr2(0,11)=-0.01;
  Corr2(1,1)=1.;
  Corr2(1,2)=0.18;
  Corr2(1,3)=-0.05;
  Corr2(1,4)=-0.1;
  Corr2(1,5)=0.13;
  Corr2(1,6)=-0.09;
  Corr2(1,7)=0.17;
  Corr2(1,8)=-0.08;
  Corr2(1,9)=-0.06;
  Corr2(1,10)=0.;
  Corr2(1,11)=-0.06;
  Corr2(2,2)=1.;
  Corr2(2,3)=0.;
  Corr2(2,4)=-0.13;
  Corr2(2,5)=-0.02;
  Corr2(2,6)=-0.09;
  Corr2(2,7)=0.16;
  Corr2(2,8)=-0.15;
  Corr2(2,9)=-0.03;
  Corr2(2,10)=-0.07;
  Corr2(2,11)=-0.03;
  Corr2(3,3)=1.;
  Corr2(3,4)=0.05;
  Corr2(3,5)=-0.05;
  Corr2(3,6)=0.05;
  Corr2(3,7)=-0.09;
  Corr2(3,8)=-0.07;
  Corr2(3,9)=-0.03;
  Corr2(3,10)=0.;
  Corr2(3,11)=0.;
  Corr2(4,4)=1.;
  Corr2(4,5)=-0.02;
  Corr2(4,6)=0.44;
  Corr2(4,7)=-0.5;
  Corr2(4,8)=0.49;
  Corr2(4,9)=-0.08;
  Corr2(4,10)=0.33;
  Corr2(4,11)=-0.03;
  Corr2(5,5)=1.;
  Corr2(5,6)=-0.6;
  Corr2(5,7)=0.19;
  Corr2(5,8)=0.;
  Corr2(5,9)=-0.05;
  Corr2(5,10)=0.;
  Corr2(5,11)=-0.04;
  Corr2(6,6)=1.;
  Corr2(6,7)=-0.44;
  Corr2(6,8)=0.25;
  Corr2(6,9)=-0.02;
  Corr2(6,10)=0.16;
  Corr2(6,11)=-0.03;
  Corr2(7,7)=1.;
  Corr2(7,8)=-0.25;
  Corr2(7,9)=0.05;
  Corr2(7,10)=-0.14;
  Corr2(7,11)=0.02;
  Corr2(8,8)=1.;
  Corr2(8,9)=0.63;
  Corr2(8,10)=0.73;
  Corr2(8,11)=0.66;
  Corr2(9,9)=1.;
  Corr2(9,10)=0.65;
  Corr2(9,11)=0.97;
  Corr2(10,10)=1.;
  Corr2(10,11)=0.74;
  Corr2(11,11)=1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  // Correlation Matrix (Model)
  Corr3.ResizeTo(12, 12);
  Corr3(0,0)=1.;
  Corr3(0,1)=0.09;
  Corr3(0,2)=0.83;
  Corr3(0,3)=0.25;
  Corr3(0,4)=0.28;
  Corr3(0,5)=0.06;
  Corr3(0,6)=0.24;
  Corr3(0,7)=0.;
  Corr3(0,8)=0.86;
  Corr3(0,9)=0.11;
  Corr3(0,10)=0.77;
  Corr3(0,11)=0.30;
  Corr3(1,1)=1.;
  Corr3(1,2)=-0.31;
  Corr3(1,3)=-0.56;
  Corr3(1,4)=-0.64;
  Corr3(1,5)=-0.74;
  Corr3(1,6)=0.80;
  Corr3(1,7)=0.51;
  Corr3(1,8)=-0.23;
  Corr3(1,9)=0.85;
  Corr3(1,10)=-0.38;
  Corr3(1,11)=-0.51;
  Corr3(2,2)=1.;
  Corr3(2,3)=0.57;
  Corr3(2,4)=0.71;
  Corr3(2,5)=0.44;
  Corr3(2,6)=-0.15;
  Corr3(2,7)=-0.15;
  Corr3(2,8)=0.91;
  Corr3(2,9)=-0.20;
  Corr3(2,10)=0.96;
  Corr3(2,11)=0.56;
  Corr3(3,3)=1.;
  Corr3(3,4)=0.61;
  Corr3(3,5)=0.93;
  Corr3(3,6)=-0.74;
  Corr3(3,7)=-0.74;
  Corr3(3,8)=0.54;
  Corr3(3,9)=-0.18;
  Corr3(3,10)=0.61;
  Corr3(3,11)=0.95;
  Corr3(4,4)=1.;
  Corr3(4,5)=0.64;
  Corr3(4,6)=-0.5;
  Corr3(4,7)=-0.16;
  Corr3(4,8)=0.52;
  Corr3(4,9)=-0.54;
  Corr3(4,10)=0.73;
  Corr3(4,11)=0.56;
  Corr3(5,5)=1.;
  Corr3(5,6)=-0.86;
  Corr3(5,7)=-0.76;
  Corr3(5,8)=0.39;
  Corr3(5,9)=-0.41;
  Corr3(5,10)=0.50;
  Corr3(5,11)=0.89;
  Corr3(6,6)=1.;
  Corr3(6,7)=0.74;
  Corr3(6,8)=-0.12;
  Corr3(6,9)=0.52;
  Corr3(6,10)=-0.23;
  Corr3(6,11)=-0.70;
  Corr3(7,7)=1.;
  Corr3(7,8)=-0.22;
  Corr3(7,9)=0.16;
  Corr3(7,10)=-0.18;
  Corr3(7,11)=-0.77;
  Corr3(8,8)=1.;
  Corr3(8,9)=-0.11;
  Corr3(8,10)=0.91;
  Corr3(8,11)=0.53;
  Corr3(9,9)=1.;
  Corr3(9,10)=-0.26;
  Corr3(9,11)=-0.14;
  Corr3(10,10)=1.;
  Corr3(10,11)=0.58;
  Corr3(11,11)=1.;
  SymmetrizeUpperTriangularMatrix(Corr3);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Babar_PRL105.121801", CorrelatedGaussianObservables(CorrData, Corr, Corr2, Corr3)));

  //-------------------------------------------------------------------------------------------------------------------------

  //-------------------------------------------------  CDF measurements  -------------------------------------------------------------------------

  // https://journals.aps.org/prd/pdf/10.1103/PhysRevD.81.031105
  // B -> DK
  // D -> KK, D-> pipi
  meas.insert(pair<string, dato>("CDF_PRD81_031105_ACPp", dato(0.39, 0.17, 0.04))); // Acp+
  meas.insert(pair<string, dato>("CDF_PRD81_031105_RCPp", dato(1.30, 0.24, 0.12)));  //RCP+

  // https://arxiv.org/pdf/1108.5765
  // B -> DK
  // D -> Kpi
  meas.insert(pair<string, dato>("CDF_1108.5765_RDK", dato(22.e-3, 8.6e-3, 2.6e-3))); // R_DK_Kpi
  meas.insert(pair<string, dato>("CDF_1108.5765_ADK", dato(-0.82, 0.44, 0.09))); // A_DK_Kpi
  // B -> Dpi
  // D -> Kpi
  meas.insert(pair<string, dato>("CDF_1108.5765_RDpi", dato(2.8e-3, 0.7e-3, 0.4e-3))); // R_Dpi_Kpi
  meas.insert(pair<string, dato>("CDF_1108.5765_ADpi", dato(0.13, 0.25, 0.02))); // A_Dpi_Kpi

  //-------------------------------------------------------------------------------------------------------------------------


  //-------------------------------------------------  Bpm -> Dhpm  -------------------------------------------------------------------------

  // GLW: D -> KK, pipi; ADS: D -> Kpi
  // https://arxiv.org/pdf/2012.09903.pdf
  // Observables 8:
  CorrData.clear();
  CorrData.push_back(dato(0.136, 0.009, 0.001));       // acp_dk
  CorrData.push_back(dato(-0.008, 0.002, 0.002));      // acp_dpi
  CorrData.push_back(dato(-0.011, 0.003, 0.002));      // afav_dk
  CorrData.push_back(dato(0.95, 0.009, 0.01));         // rcp
  CorrData.push_back(dato(0.0095, 0.0005, 0.0003));    // rm_dk
  CorrData.push_back(dato(0.00415, 0.00008, 0.00004)); // rm_dpi
  CorrData.push_back(dato(0.02520, 0.0008, 0.0004));   // rp_dk
  CorrData.push_back(dato(0.00320, 0.00007, 0.00004)); // rp_dpi

  // Correlation Matrix (stat):
  Corr.ResizeTo(8, 8);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.;
  Corr(0, 2) = 0.02;
  Corr(0, 3) = -0.02;
  Corr(0, 4) = 0.;
  Corr(0, 5) = 0.;
  Corr(0, 6) = 0.;
  Corr(0, 7) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.08;
  Corr(1, 3) = 0.;
  Corr(1, 4) = 0.01;
  Corr(1, 5) = 0.01;
  Corr(1, 6) = 0.;
  Corr(1, 7) = -0.01;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.;
  Corr(2, 4) = 0.;
  Corr(2, 5) = 0.01;
  Corr(2, 6) = 0.;
  Corr(2, 7) = -0.01;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.03;
  Corr(3, 5) = 0.;
  Corr(3, 6) = 0.04;
  Corr(3, 7) = 0.;
  Corr(4, 4) = 1.;
  Corr(4, 5) = -0.03;
  Corr(4, 6) = 0.05;
  Corr(4, 7) = 0.02;
  Corr(5, 5) = 1.;
  Corr(5, 6) = 0.02;
  Corr(5, 7) = 0.08;
  Corr(6, 6) = 1.;
  Corr(6, 7) = -0.04;
  Corr(7, 7) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);
  // Correlation Matrix (syst)
  Corr2.ResizeTo(8, 8);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.13;
  Corr2(0, 2) = 0.07;
  Corr2(0, 3) = -0.22;
  Corr2(0, 4) = -0.03;
  Corr2(0, 5) = 0.;
  Corr2(0, 6) = 0.02;
  Corr2(0, 7) = -0.07;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = -0.74;
  Corr2(1, 3) = 0.18;
  Corr2(1, 4) = 0.04;
  Corr2(1, 5) = 0.31;
  Corr2(1, 6) = -0.07;
  Corr2(1, 7) = -0.21;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = -0.02;
  Corr2(2, 4) = 0.05;
  Corr2(2, 5) = -0.24;
  Corr2(2, 6) = 0.11;
  Corr2(2, 7) = 0.22;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = 0.09;
  Corr2(3, 5) = 0.11;
  Corr2(3, 6) = 0.17;
  Corr2(3, 7) = 0.14;
  Corr2(4, 4) = 1.;
  Corr2(4, 5) = 0.16;
  Corr2(4, 6) = 0.93;
  Corr2(4, 7) = 0.13;
  Corr2(5, 5) = 1.;
  Corr2(5, 6) = 0.21;
  Corr2(5, 7) = 0.84;
  Corr2(6, 6) = 1.;
  Corr2(6, 7) = 0.26;
  Corr2(7, 7) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID0", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // GLW: D -> KKpipi, D -> 4pi
  // https://arxiv.org/abs/2301.10328
  // 6 Observables
  CorrData.clear();
  CorrData.push_back(dato(0.095, 0.023, 0.002));     // acp_dk_kkpipi
  CorrData.push_back(dato(-0.009, 0.006, 0.001));    // acp_dpi_kkpipi
  CorrData.push_back(dato(0.061, 0.013, 0.002));     // acp_dk_4pi
  CorrData.push_back(dato(-0.0082, 0.0031, 0.0007)); // acp_dpi_4pi
  CorrData.push_back(dato(0.974, 0.024, 0.015));     // rcp_kkpipi
  CorrData.push_back(dato(0.978, 0.014, 0.010));     // rcp_4pi
  // Correlation Matrix (stat):
  Corr.ResizeTo(6, 6);
  Corr(0, 0) = 1.;
  Corr(0, 1) = -0.025;
  Corr(0, 2) = 0.00;
  Corr(0, 3) = 0.;
  Corr(0, 4) = 0.015;
  Corr(0, 5) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.;
  Corr(1, 3) = 0.;
  Corr(1, 4) = 0.002;
  Corr(1, 5) = 0.;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.028;
  Corr(2, 4) = 0.002;
  Corr(2, 5) = 0.016;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.;
  Corr(3, 5) = 0.002;
  Corr(4, 4) = 1.;
  Corr(4, 5) = 0.068;
  Corr(5, 5) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(6, 6);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.345;
  Corr2(0, 2) = 0.758;
  Corr2(0, 3) = 0.343;
  Corr2(0, 4) = -0.589;
  Corr2(0, 5) = -0.056;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.423;
  Corr2(1, 3) = 0.999;
  Corr2(1, 4) = 0.002;
  Corr2(1, 5) = 0.007;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.420;
  Corr2(2, 4) = -0.033;
  Corr2(2, 5) = -0.246;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = 0.009;
  Corr2(3, 5) = 0.014;
  Corr2(4, 4) = 1.;
  Corr2(4, 5) = 0.236;
  Corr2(5, 5) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);


  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("GLW_2301.10328", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  // GLW: D -> pipipi0, KKpi0; D -> Kpipi0
  //https://arxiv.org/pdf/2112.10617
  // Observables 11:
  CorrData.clear();
  CorrData.push_back(dato(1.021, 0.079, 0.005));       // rcp_kkpi0
  CorrData.push_back(dato(0.902, 0.041, 0.004));       // rcp_pipipi0
  CorrData.push_back(dato(-0.024, 0.013, 0.002));      // afav_dk_kpipi0
  CorrData.push_back(dato(0.067, 0.073, 0.003));       // acp_dk_kkpi0
  CorrData.push_back(dato(0.109, 0.043, 0.003));       // acp_dk_pipipi0
  CorrData.push_back(dato(-0.001, 0.019, 0.002));      // acp_dpi_kkpi0
  CorrData.push_back(dato(0.001, 0.010, 0.002));       // acp_dpi_pipipi0
  CorrData.push_back(dato(0.0179, 0.0024, 0.0003));    // R_k+
  CorrData.push_back(dato(0.0085, 0.0020, 0.0004));    // R_k-
  CorrData.push_back(dato(0.00188, 0.00027, 0.00005)); // R_pi+
  CorrData.push_back(dato(0.00227, 0.00028, 0.00004)); // R_pi-
  // Correlation Matrix (stat):
  Corr.ResizeTo(11, 11);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.05;
  Corr(0, 2) = 0.0;
  Corr(0, 3) = -0.01;
  Corr(0, 4) = 0.0;
  Corr(0, 5) = 0.0;
  Corr(0, 6) = 0.0;
  Corr(0, 7) = 0.02;
  Corr(0, 8) = 0.01;
  Corr(0, 9) = 0.;
  Corr(0, 10) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.0;
  Corr(1, 3) = 0.0;
  Corr(1, 4) = -0.05;
  Corr(1, 5) = 0.0;
  Corr(1, 6) = 0.0;
  Corr(1, 7) = 0.02;
  Corr(1, 8) = 0.01;
  Corr(1, 9) = 0.;
  Corr(1, 10) = 0.;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.01;
  Corr(2, 4) = 0.02;
  Corr(2, 5) = 0.04;
  Corr(2, 6) = 0.08;
  Corr(2, 7) = -0.01;
  Corr(2, 8) = 0.;
  Corr(2, 9) = -0.01;
  Corr(2, 10) = 0.01;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.0;
  Corr(3, 5) = -0.01;
  Corr(3, 6) = 0.01;
  Corr(3, 7) = 0.;
  Corr(3, 8) = 0.;
  Corr(3, 9) = 0.;
  Corr(3, 10) = 0.0;
  Corr(4, 4) = 1.;
  Corr(4, 5) = 0.01;
  Corr(4, 6) = -0.01;
  Corr(4, 7) = 0.;
  Corr(4, 8) = 0.;
  Corr(4, 9) = 0.;
  Corr(4, 10) = 0.;
  Corr(5, 5) = 1.;
  Corr(5, 6) = 0.04;
  Corr(5, 7) = 0.;
  Corr(5, 8) = 0.;
  Corr(5, 9) = 0.;
  Corr(5, 10) = 0.;
  Corr(6, 6) = 1.;
  Corr(6, 7) = -0.01;
  Corr(6, 8) = 0.;
  Corr(6, 9) = -0.01;
  Corr(6, 10) = 0.01;
  Corr(7, 7) = 1.;
  Corr(7, 8) = 0.09;
  Corr(7, 9) = -0.01;
  Corr(7, 10) = 0.0;
  Corr(8, 8) = 1.;
  Corr(8, 9) = 0.;
  Corr(8, 10) = -0.01;
  Corr(9, 9) = 1.;
  Corr(9, 10) = 0.01;
  Corr(10, 10) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(11, 11);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.83;
  Corr2(0, 2) = -0.29;
  Corr2(0, 3) = 0.09;
  Corr2(0, 4) = 0.15;
  Corr2(0, 5) = -0.23;
  Corr2(0, 6) = -0.25;
  Corr2(0, 7) = -0.1;
  Corr2(0, 8) = -0.05;
  Corr2(0, 9) = -0.15;
  Corr2(0, 10) = -0.18;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = -0.39;
  Corr2(1, 3) = -0.07;
  Corr2(1, 4) = 0.07;
  Corr2(1, 5) = -0.26;
  Corr2(1, 6) = -0.28;
  Corr2(1, 7) = -0.22;
  Corr2(1, 8) = -0.08;
  Corr2(1, 9) = -0.15;
  Corr2(1, 10) = -0.14;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = -0.06;
  Corr2(2, 4) = -0.15;
  Corr2(2, 5) = 0.91;
  Corr2(2, 6) = 0.92;
  Corr2(2, 7) = -0.04;
  Corr2(2, 8) = -0.12;
  Corr2(2, 9) = -0.09;
  Corr2(2, 10) = -0.06;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = 0.89;
  Corr2(3, 5) = -0.2;
  Corr2(3, 6) = -0.19;
  Corr2(3, 7) = 0.94;
  Corr2(3, 8) = 0.90;
  Corr2(3, 9) = 0.66;
  Corr2(3, 10) = 0.57;
  Corr2(4, 4) = 1.;
  Corr2(4, 5) = -0.23;
  Corr2(4, 6) = -0.23;
  Corr2(4, 7) = 0.88;
  Corr2(4, 8) = 0.84;
  Corr2(4, 9) = 0.54;
  Corr2(4, 10) = 0.44;
  Corr2(5, 5) = 1.;
  Corr2(5, 6) = 1.;
  Corr2(5, 7) = -0.18;
  Corr2(5, 8) = -0.17;
  Corr2(5, 9) = -0.05;
  Corr2(5, 10) = 0.0;
  Corr2(6, 6) = 1.;
  Corr2(6, 7) = -0.17;
  Corr2(6, 8) = -0.16;
  Corr2(6, 9) = -0.06;
  Corr2(6, 10) = -0.02;
  Corr2(7, 7) = 1.;
  Corr2(7, 8) = 0.92;
  Corr2(7, 9) = .73;
  Corr2(7, 10) = .63;
  Corr2(8, 8) = 1.;
  Corr2(8, 9) = .75;
  Corr2(8, 10) = .71;
  Corr2(9, 9) = 1.;
  Corr2(9, 10) = .97;
  Corr2(10, 10) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("2112.10617", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // ADS: D -> K0sKpi
  // https://arxiv.org/pdf/2002.08858
  // Observables 7:
  CorrData.clear();
  CorrData.push_back(dato(-0.02, 0.011, 0.003));  // afav_dpi_kskpi
  CorrData.push_back(dato(0.007, 0.017, 0.003));  // asup_dpi_kskpi
  CorrData.push_back(dato(0.084, 0.049, 0.008));  // afav_dk_kskpi
  CorrData.push_back(dato(0.021, 0.094, 0.017));  // asup_dk_kskpi
  CorrData.push_back(dato(2.585, 0.057, 0.019));  // rfavsup_dpi_kskpi
  CorrData.push_back(dato(0.079, 0.004, 0.002));  // rfav_dkdpi_kskpi
  CorrData.push_back(dato(0.0620, 0.006, 0.003)); // rsup_dkdpi_kskpi

  // Correlation Matrix (stat):
  Corr.ResizeTo(7, 7);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.;
  Corr(0, 2) = -0.05;
  Corr(0, 3) = 0.;
  Corr(0, 4) = 0.;
  Corr(0, 5) = -0.01;
  Corr(0, 6) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.;
  Corr(1, 3) = -0.05;
  Corr(1, 4) = 0.;
  Corr(1, 5) = 0.;
  Corr(1, 6) = 0.;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.;
  Corr(2, 4) = 0.;
  Corr(2, 5) = 0.;
  Corr(2, 6) = 0.;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.;
  Corr(3, 5) = 0.;
  Corr(3, 6) = -0.02;
  Corr(4, 4) = 1.;
  Corr(4, 5) = -0.11;
  Corr(4, 6) = 0.15;
  Corr(5, 5) = 1.;
  Corr(5, 6) = 0.06;
  Corr(6, 6) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);
  // Correlation Matrix (syst)
  Corr2.ResizeTo(7, 7);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = -0.88;
  Corr2(0, 2) = 0.74;
  Corr2(0, 3) = 0.05;
  Corr2(0, 4) = 0.;
  Corr2(0, 5) = 0.;
  Corr2(0, 6) = 0.;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = -0.73;
  Corr2(1, 3) = -0.08;
  Corr2(1, 4) = 0.;
  Corr2(1, 5) = 0.;
  Corr2(1, 6) = 0.;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.63;
  Corr2(2, 4) = 0.;
  Corr2(2, 5) = -0.17;
  Corr2(2, 6) = 0.;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = 0.;
  Corr2(3, 5) = 0.;
  Corr2(3, 6) = -0.1;
  Corr2(4, 4) = 1.;
  Corr2(4, 5) = 0.;
  Corr2(4, 6) = 0.;
  Corr2(5, 5) = 1.;
  Corr2(5, 6) = 0.25;
  Corr2(6, 6) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID4", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // GLW: D -> KK, D -> K0spi0 (Belle)
  // https://arxiv.org/abs/2308.05048
  // Observables 4:
  CorrData.clear();
  CorrData.push_back(dato(-0.167, 0.057, 0.006)); // acp_dk_k0pi0
  CorrData.push_back(dato(1.151, 0.074, 0.019));  // Rcp_dk_k0pi0
  CorrData.push_back(dato(0.125, 0.058, 0.014));  // acp_dk_kk
  CorrData.push_back(dato(1.164, 0.081, 0.036));  // Rcp_dk_kk

  // Correlation Matrix (stat):
  Corr.ResizeTo(4, 4);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.056;
  Corr(0, 2) = 0.;
  Corr(0, 3) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.;
  Corr(1, 3) = -0.081;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.060;
  Corr(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (stat):
  Corr2.ResizeTo(4, 4);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = -0.490;
  Corr2(0, 2) = 0.542;
  Corr2(0, 3) = 0.005;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = -0.128;
  Corr2(1, 3) = -0.063;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.342;
  Corr2(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("2308.05048", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // ADS: D -> Kpipi0 (Belle)
  // https://arxiv.org/pdf/1310.1741
  // Observables 4:
  CorrData.clear();
  CorrData.push_back(dato(0.0198, 0.0062, 0.0024));    // RADS_dk_kpipi0
  CorrData.push_back(dato(0.41, 0.307, 0.05));         // asup_dk_kpipi0
  CorrData.push_back(dato(0.00189, 0.00054, 0.00024)); // RADS_dpi_kpipi0
  CorrData.push_back(dato(0.16, 0.27, 0.04));          // asup_dpi_kpipi0

  // Correlation Matrix (stat):
  Corr.ResizeTo(4, 4);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.;
  Corr(0, 2) = 0.;
  Corr(0, 3) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.;
  Corr(1, 3) = 0.;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.;
  Corr(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (stat):
  Corr2.ResizeTo(4, 4);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.;
  Corr2(0, 2) = 0.;
  Corr2(0, 3) = 0.;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.;
  Corr2(1, 3) = 0.;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.;
  Corr2(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Belle_PRD88_2013", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // ADS: D -> Kpi (Belle)
  // https://arxiv.org/abs/1103.5951
  // Observables 4:
  CorrData.clear();
  CorrData.push_back(dato(0.0163, 0.0042, 0.001));     // RADS_dk_kpi
  CorrData.push_back(dato(-0.39, 0.27, 0.04));         // asup_dk_kpi
  CorrData.push_back(dato(0.00328, 0.00037, 0.00015)); // RADS_dpi_kpi
  CorrData.push_back(dato(-0.04, 0.11, 0.02));         // asup_dpi_kpi

  // Correlation Matrix (stat):
  Corr.ResizeTo(4, 4);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.242;
  Corr(0, 2) = 0.;
  Corr(0, 3) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.;
  Corr(1, 3) = 0.;
  Corr(2, 2) = 1.;
  Corr(2, 3) = -0.032;
  Corr(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (stat):
  Corr2.ResizeTo(4, 4);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.;
  Corr2(0, 2) = 0.;
  Corr2(0, 3) = 0.;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.;
  Corr2(1, 3) = 0.;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.;
  Corr2(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Belle_PRL106_2011", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // K^*+- region fit: ADS: D -> K0sKpi (Belle)
  // 2306.02940
  // Observables 7:
  CorrData.clear();
  CorrData.push_back(dato(0.055, 0.119, 0.020)); // afav_dK_kskpi
  CorrData.push_back(dato(0.231, 0.184, 0.014)); // asup_dK_kskpi
  CorrData.push_back(dato(0.046, 0.029, 0.016)); // afav_dpi_kskpi
  CorrData.push_back(dato(0.009, 0.046, 0.009)); // asup_dpi_kskpi
  CorrData.push_back(dato(0.093, 0.012, 0.005)); // rfav_dkdpi_kskpi
  CorrData.push_back(dato(0.103, 0.020, 0.006)); // rsup_dkdpi_kskpi
  CorrData.push_back(dato(2.412, 0.132, 0.019)); // rfavsup_dpi_kskpi

  // Correlation Matrix (stat):
  Corr.ResizeTo(7, 7);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.003;
  Corr(0, 2) = -0.012;
  Corr(0, 3) = 0.001;
  Corr(0, 4) = -0.052;
  Corr(0, 5) = -0.013;
  Corr(0, 6) = 0.002;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.001;
  Corr(1, 3) = -0.011;
  Corr(1, 4) = -0.004;
  Corr(1, 5) = -0.034;
  Corr(1, 6) = 0.002;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.001;
  Corr(2, 4) = 0.002;
  Corr(2, 5) = -0.004;
  Corr(2, 6) = -0.011;
  Corr(3, 3) = 1.;
  Corr(3, 4) = -0.002;
  Corr(3, 5) = -0.002;
  Corr(3, 6) = 0.014;
  Corr(4, 4) = 1.;
  Corr(4, 5) = 0.034;
  Corr(4, 6) = -0.133;
  Corr(5, 5) = 1.;
  Corr(5, 6) = 0.208;
  Corr(6, 6) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);
  // Correlation Matrix (syst)
  Corr2.ResizeTo(7, 7);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.195;
  Corr2(0, 2) = 0.047;
  Corr2(0, 3) = 0.013;
  Corr2(0, 4) = 0.120;
  Corr2(0, 5) = -0.053;
  Corr2(0, 6) = 0.191;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.038;
  Corr2(1, 3) = 0.004;
  Corr2(1, 4) = 0.344;
  Corr2(1, 5) = 0.21;
  Corr2(1, 6) = 0.007;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.023;
  Corr2(2, 4) = -0.004;
  Corr2(2, 5) = -0.037;
  Corr2(2, 6) = 0.018;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = -0.017;
  Corr2(3, 5) = -0.024;
  Corr2(3, 6) = 0.006;
  Corr2(4, 4) = 1.;
  Corr2(4, 5) = 0.915;
  Corr2(4, 6) = 0.015;
  Corr2(5, 5) = 1.;
  Corr2(5, 6) = -0.097;
  Corr2(6, 6) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("2306.02940", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // GGSZ: D -> K3pi
  // https://arxiv.org/pdf/2209.03692
  // Observables 6:
  CorrData.clear();
  CorrData.push_back(dato(-0.0917, 0.0074, 0.));       // xp_dk
  CorrData.push_back(dato(0.01670, 0.0175, 0.));       // xm_dk
  CorrData.push_back(dato(-0.0154, 0.0234, 0.));       // yp_dk
  CorrData.push_back(dato(0.0929, 0.0047, 0.));        // ym_dk
  CorrData.push_back(dato(-0.0471, 0.0122, 0.));       // xi_x_dpi
  CorrData.push_back(dato(0.0024, 0.0178, 0.));        // xi_y_dpi

  // Correlation Matrix (stat):
  Corr.ResizeTo(6, 6);
  Corr(0, 0) = 1.;
  Corr(0, 1) = -.357;
  Corr(0, 2) = -0.740;
  Corr(0, 3) = -0.310;
  Corr(0, 4) = 0.047;
  Corr(0, 5) = 0.298;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.274;
  Corr(1, 3) = -0.352;
  Corr(1, 4) = 0.114;
  Corr(1, 5) = -0.112;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.174;
  Corr(2, 4) = -0.13;
  Corr(2, 5) = -0.321;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.38;
  Corr(3, 5) = -0.139;
  Corr(4, 4) = 1.;
  Corr(4, 5) = 0.095;
  Corr(5, 5) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(6, 6);
  Corr2.UnitMatrix();

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("2209.03692", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // GGSZ D -> K0spipipi0 (Belle)
  // https://arxiv.org/pdf/1908.09499
  // Observables 8:
  CorrData.clear();
  CorrData.push_back(dato(0.095, 0.121, 0.029));  // xm_dk
  CorrData.push_back(dato(0.354, 0.170, 0.045));  // ym_dk
  CorrData.push_back(dato(-0.03, 0.121, 0.026));  // xp_dk
  CorrData.push_back(dato(0.22, 0.376, 0.079));   // yp_dk
  CorrData.push_back(dato(-0.014, 0.021, 0.021)); // xm_dpi
  CorrData.push_back(dato(-0.033, 0.059, 0.023)); // ym_dpi
  CorrData.push_back(dato(0.039, 0.024, 0.020));  // xp_dpi
  CorrData.push_back(dato(-0.196, 0.069, 0.048)); // yp_dpi

  // Correlation Matrix (stat):
  Corr.ResizeTo(8, 8);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.486;
  Corr(0, 2) = 0.172;
  Corr(0, 3) = -0.231;
  Corr(0, 4) = 0.;
  Corr(0, 5) = 0.;
  Corr(0, 6) = 0.;
  Corr(0, 7) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -0.127;
  Corr(1, 3) = 0.179;
  Corr(1, 4) = 0.;
  Corr(1, 5) = 0.;
  Corr(1, 6) = 0.;
  Corr(1, 7) = 0.;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.365;
  Corr(2, 4) = 0.;
  Corr(2, 5) = 0.;
  Corr(2, 6) = 0.;
  Corr(2, 7) = 0.;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.;
  Corr(3, 5) = 0.;
  Corr(3, 6) = 0.;
  Corr(3, 7) = 0.;
  Corr(4, 4) = 1.;
  Corr(4, 5) = -0.364;
  Corr(4, 6) = 0.314;
  Corr(4, 7) = 0.050;
  Corr(5, 5) = 1.;
  Corr(5, 6) = 0.347;
  Corr(5, 7) = 0.055;
  Corr(6, 6) = 1.;
  Corr(6, 7) = -0.032;
  Corr(7, 7) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (stat):
  Corr2.ResizeTo(8, 8);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.;
  Corr2(0, 2) = 0.;
  Corr2(0, 3) = 0.;
  Corr2(0, 4) = 0.;
  Corr2(0, 5) = 0.;
  Corr2(0, 6) = 0.;
  Corr2(0, 7) = 0.;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.;
  Corr2(1, 3) = 0.;
  Corr2(1, 4) = 0.;
  Corr2(1, 5) = 0.;
  Corr2(1, 6) = 0.;
  Corr2(1, 7) = 0.;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.;
  Corr2(2, 4) = 0.;
  Corr2(2, 5) = 0.;
  Corr2(2, 6) = 0.;
  Corr2(2, 7) = 0.;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = 0.;
  Corr2(3, 5) = 0.;
  Corr2(3, 6) = 0.;
  Corr2(3, 7) = 0.;
  Corr2(4, 4) = 1.;
  Corr2(4, 5) = 0.;
  Corr2(4, 6) = 0.;
  Corr2(4, 7) = 0.;
  Corr2(5, 5) = 1.;
  Corr2(5, 6) = 0.;
  Corr2(5, 7) = 0.;
  Corr2(6, 6) = 1.;
  Corr2(6, 7) = 0.;
  Corr2(7, 7) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("1908.09499", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // GGSZ: D -> K0spipi, D -> K0sKK (Belle)
  // https://arxiv.org/abs/2110.12125
  // Observables 6:
  CorrData.clear();
  CorrData.push_back(dato(0.0924, 0.0327, 0.0017, 0.0023));  // xm_dk
  CorrData.push_back(dato(0.10, 0.0420, 0.0023, 0.0067));    // ym_dk
  CorrData.push_back(dato(-0.1128, 0.0315, 0.0018, 0.0022)); // xp_dk
  CorrData.push_back(dato(-0.0455, 0.0420, 0.0011, 0.0055)); // yp_dk
  CorrData.push_back(dato(-0.1109, 0.0475, 0.0051, 0.0073)); // xi_x_dpi
  CorrData.push_back(dato(-0.079, 0.0544, 0.0019, 0.0082));  // xi_y_dpi

  // Correlation Matrix (stat):
  Corr.ResizeTo(6, 6);
  Corr(0, 0) = 1.;
  Corr(0, 1) = -0.204;
  Corr(0, 2) = -0.051;
  Corr(0, 3) = 0.063;
  Corr(0, 4) = 0.365;
  Corr(0, 5) = -0.151;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.014;
  Corr(1, 3) = -0.051;
  Corr(1, 4) = -0.090;
  Corr(1, 5) = 0.404;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.152;
  Corr(2, 4) = -0.330;
  Corr(2, 5) = -0.057;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.026;
  Corr(3, 5) = -0.391;
  Corr(4, 4) = 1.;
  Corr(4, 5) = 0.080;
  Corr(5, 5) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);
  // Correlation Matrix (syst)
  Corr2.ResizeTo(6, 6);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.777;
  Corr2(0, 2) = 0.483;
  Corr2(0, 3) = 0.268;
  Corr2(0, 4) = 0.839;
  Corr2(0, 5) = 0.698;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.411;
  Corr2(1, 3) = .504;
  Corr2(1, 4) = .802;
  Corr2(1, 5) = .797;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.680;
  Corr2(2, 4) = 0.766;
  Corr2(2, 5) = .377;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = 0.480;
  Corr2(3, 5) = .303;
  Corr2(4, 4) = 1.;
  Corr2(4, 5) = 0.638;
  Corr2(5, 5) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);
  // Correlation Matrix (syst)
  Corr3.ResizeTo(6, 6);
  Corr3(0, 0) = 1.;
  Corr3(0, 1) = -.113;
  Corr3(0, 2) = .069;
  Corr3(0, 3) = .406;
  Corr3(0, 4) = -.016;
  Corr3(0, 5) = .114;
  Corr3(1, 1) = 1.;
  Corr3(1, 2) = 0.038;
  Corr3(1, 3) = -.196;
  Corr3(1, 4) = -.692;
  Corr3(1, 5) = -.106;
  Corr3(2, 2) = 1.;
  Corr3(2, 3) = 0.412;
  Corr3(2, 4) = -.226;
  Corr3(2, 5) = -0.469;
  Corr3(3, 3) = 1.;
  Corr3(3, 4) = 0.180;
  Corr3(3, 5) = -0.069;
  Corr3(4, 4) = 1.;
  Corr3(4, 5) = 0.622;
  Corr3(5, 5) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr3);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("2110.12125", CorrelatedGaussianObservables(CorrData, Corr, Corr2, Corr3)));


  // GGSZ: D -> K0sKK, D -> K0spipi
  // https://arxiv.org/abs/2301.10328
  // 6 Observables
  CorrData.clear();
  CorrData.push_back(dato(0.079, 0.029, 0.0057));  // xm_dk
  CorrData.push_back(dato(-0.033, 0.034, 0.036));  // ym_dk
  CorrData.push_back(dato(-0.125, 0.025, 0.017));  // xp_dk
  CorrData.push_back(dato(-0.042, 0.031, 0.013));  // yp_dk
  CorrData.push_back(dato(-0.031, 0.035, 0.0071)); // xi_x_dpi
  CorrData.push_back(dato(-0.017, 0.047, 0.012));  // xi_y_dpi
  // Correlation Matrix (stat):
  Corr.ResizeTo(6, 6);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.032;
  Corr(0, 2) = 0.008;
  Corr(0, 3) = -0.010;
  Corr(0, 4) = 0.034;
  Corr(0, 5) = 0.102;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.017;
  Corr(1, 3) = 0.;
  Corr(1, 4) = -0.091;
  Corr(1, 5) = 0.08;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.007;
  Corr(2, 4) = -0.100;
  Corr(2, 5) = 0.051;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.012;
  Corr(3, 5) = -0.097;
  Corr(4, 4) = 1.;
  Corr(4, 5) = 0.014;
  Corr(5, 5) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);
  // Correlation Matrix (syst)
  Corr2.ResizeTo(6, 6);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = -0.678;
  Corr2(0, 2) = 0.751;
  Corr2(0, 3) = 0.736;
  Corr2(0, 4) = -0.048;
  Corr2(0, 5) = -0.650;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = -0.973;
  Corr2(1, 3) = -0.961;
  Corr2(1, 4) = -0.200;
  Corr2(1, 5) = 0.898;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.971;
  Corr2(2, 4) = 0.166;
  Corr2(2, 5) = -0.862;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = 0.065;
  Corr2(3, 5) = -0.913;
  Corr2(4, 4) = 1.;
  Corr2(4, 5) = -0.057;
  Corr2(5, 5) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("GGSZ_2301.10328", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));



  // If combiining chargedB modes only
  // If combining all the modes see the neutralBd section since they are all correlated
  // https://arxiv.org/pdf/2010.08483 + https://arxiv.org/abs/2310.04277 +  https://arxiv.org/abs/2311.10434 + LHCB-PAPER-2024-023
  if (comb == 0){

    // GGSZ LHCb only ChargedB
    // Observables 22

    // D -> K0spipi, D -> K0sKK
    // https://arxiv.org/pdf/2010.08483
    // Observables 6:
    CorrData.clear();
    CorrData.push_back(dato(0.0568, 0.0096, 0.0020, 0.0023));     // xm_dk
    CorrData.push_back(dato(0.06550, 0.01140, 0.0025, 0.0035)); // ym_dk
    CorrData.push_back(dato(-0.09300, 0.0098, 0.0024, 0.0018));    // xp_dk
    CorrData.push_back(dato(-0.0125, 0.0123, 0.0026, 0.0028));    // yp_dk
    CorrData.push_back(dato(-0.05470, 0.0199, 0.0032, 0.0014));   // xi_x_dpi
    CorrData.push_back(dato(0.00710, 0.02330, 0.0054, 0.0018));   // xi_y_dpi

    //----------------------------------------------------------------------------------------------------------------------------------

    //-------------------------------------------------  Bpm -> D*hpm -------------------------------------------------------------------------

    // GGSZ: D -> K0spipi, D -> K0sKK
    // https://arxiv.org/abs/2310.04277 LHCb
    // Observables 6:
    CorrData.push_back(dato(11.42e-2, 3.16e-2, 1.26e-2, 0.41e-2));  // xp_dstk
    CorrData.push_back(dato(-8.91e-2, 3.55e-2, 2.04e-2, 0.23e-2));  // xm_dstk
    CorrData.push_back(dato(3.60e-2, 4.41e-2, 2.12e-2, 0.30e-2));   // yp_dstk
    CorrData.push_back(dato(-16.75e-2, 3.98e-2, 1.48e-2, 0.64e-2)); // ym_dstk
    CorrData.push_back(dato(0.51e-2, 5.00e-2, 2.66e-2, 0.93e-2));   // Re xi_dstpi
    CorrData.push_back(dato(7.92e-2, 5.04e-2, 3.78e-2, 0.83e-2));   // Im xi_dstpi

    // GGSZ: D -> K0spipi, D -> K0sKK
    // https://arxiv.org/abs/2311.10434 LHCb
    // Observables 6:
    CorrData.push_back(dato(-6.3e-2, 2.9e-2, 1.1e-2, 0.6e-2)); // xm_dstk
    CorrData.push_back(dato(-4.8e-2, 5.7e-2, 1.4e-2, 1.5e-2)); // ym_dstk
    CorrData.push_back(dato(6.0e-2, 2.6e-2, 0.9e-2, 0.2e-2));  // xp_dstk
    CorrData.push_back(dato(5.4e-2, 2.9e-2, 0.9e-2, 0.4e-2));  // yp_dstk
    CorrData.push_back(dato(11.5e-2, 9.4e-2, 3.3e-2, 2.3e-2)); // Re xi_dstpi
    CorrData.push_back(dato(-0.9e-2, 9.7e-2, 2.5e-2, 2.1e-2)); // Im xi_dstpi

    //----------------------------------------------------------------------------------------------------------------------------------

    //-------------------------------------------------  Bpm -> DK*pm -------------------------------------------------------------------------

    // GGSZ: D -> K0spipi, D -> K0sKK
    // LHCB-PAPER-2024-023
    // Observables 4:
    CorrData.push_back(dato(0.135, 0.056, 0.019, 0.009)); // xm_dkst
    CorrData.push_back(dato(-0.170, 0.108, 0.013, 0.039)); // ym_dkst
    CorrData.push_back(dato(0.003, 0.052, 0.018, 0.004)); // xp_dkst
    CorrData.push_back(dato(0.054, 0.061, 0.009, 0.019)); // yp_dkst


    // Correlation Matrix (stat):
    Corr.ResizeTo(22, 22);
    Corr = 0.;

    Corr(0, 0) = 1.;
    Corr(0, 1) = -0.125;
    Corr(0, 2) = -0.013;
    Corr(0, 3) = 0.019;
    Corr(0, 4) = 0.037;
    Corr(0, 5) = -0.161;
    Corr(1, 1) = 1.;
    Corr(1, 2) = -0.011;
    Corr(1, 3) = -0.010;
    Corr(1, 4) = 0.097;
    Corr(1, 5) = 0.041;
    Corr(2, 2) = 1.;
    Corr(2, 3) = 0.105;
    Corr(2, 4) = -0.108;
    Corr(2, 5) = 0.032;
    Corr(3, 3) = 1.;
    Corr(3, 4) = -0.07;
    Corr(3, 5) = -0.147;
    Corr(4, 4) = 1.;
    Corr(4, 5) = 0.15;
    Corr(5, 5) = 1.;

    Corr(6, 6) = 1.;
    Corr(6, 7) = -0.05;
    Corr(6, 8) = -0.08;
    Corr(6, 9) = -0.08;
    Corr(6, 10) = -0.16;
    Corr(6, 11) = -0.14;
    Corr(7, 7) = 1.;
    Corr(7, 8) = 0.03;
    Corr(7, 9) = -0.08;
    Corr(7, 10) = 0.25;
    Corr(7, 11) = 0.06;
    Corr(8, 8) = 1.;
    Corr(8, 9) = -0.09;
    Corr(8, 10) = 0.;
    Corr(8, 11) = -0.19;
    Corr(9, 9) = 1.;
    Corr(9, 10) = 0.;
    Corr(9, 11) = 0.33;
    Corr(10, 10) = 1.;
    Corr(10, 11) = 0.18;
    Corr(11, 11) = 1.;

    Corr(12, 12) = 1.;
    Corr(12, 13) = 0.420;
    Corr(12, 14) = 0.158;
    Corr(12, 15) = 0.105;
    Corr(12, 16) = 0.445;
    Corr(12, 17) = 0.422;
    Corr(13, 13) = 1.;
    Corr(13, 14) = 0.115;
    Corr(13, 15) = 0.232;
    Corr(13, 16) = 0.765;
    Corr(13, 17) = 0.631;
    Corr(14, 14) = 1.;
    Corr(14, 15) = -0.095;
    Corr(14, 16) = 0.012;
    Corr(14, 17) = 0.409;
    Corr(15, 15) = 1.;
    Corr(15, 16) = 0.263;
    Corr(15, 17) = 0.112;
    Corr(16, 16) = 1.;
    Corr(16, 17) = 0.597;
    Corr(17, 17) = 1.;

    Corr(18,18) = 1.;
    Corr(18,19) = 0.448;
    Corr(18,20) = 0.;
    Corr(18,21) = 0.;
    Corr(19,19) = 1.;
    Corr(19,20) = 0.;
    Corr(19,21) = 0.;
    Corr(20,20) = 1.;
    Corr(20,21) = -0.144;
    Corr(21,21) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr);

    // Correlation Matrix (syst)
    Corr2.ResizeTo(22, 22);
    Corr2 = 0.;

    Corr2(0, 0) = 1.;
    Corr2(0, 1) = 0.864;
    Corr2(0, 2) = 0.734;
    Corr2(0, 3) = 0.897;
    Corr2(0, 4) = 0.349;
    Corr2(0, 5) = 0.318;
    Corr2(1, 1) = 1.;
    Corr2(1, 2) = 0.874;
    Corr2(1, 3) = 0.903;
    Corr2(1, 4) = 0.408;
    Corr2(1, 5) = 0.362;
    Corr2(2, 2) = 1.;
    Corr2(2, 3) = 0.771;
    Corr2(2, 4) = 0.563;
    Corr2(2, 5) = 0.447;
    Corr2(3, 3) = 1.;
    Corr2(3, 4) = 0.507;
    Corr2(3, 5) = 0.451;
    Corr2(4, 4) = 1.;
    Corr2(4, 5) = 0.484;
    Corr2(5, 5) = 1.;

    Corr2(6, 6) = 1.;
    Corr2(6, 7) = -0.01;
    Corr2(6, 8) = 0.;
    Corr2(6, 9) = -0.01;
    Corr2(6, 10) = -0.01;
    Corr2(6, 11) = 0.06;
    Corr2(7, 7) = 1.;
    Corr2(7, 8) = 0.02;
    Corr2(7, 9) = 0.01;
    Corr2(7, 10) = -0.02;
    Corr2(7, 11) = -0.05;
    Corr2(8, 8) = 1.;
    Corr2(8, 9) = -0.01;
    Corr2(8, 10) = -0.09;
    Corr2(8, 11) = -0.02;
    Corr2(9, 9) = 1.;
    Corr2(9, 10) = 0.01;
    Corr2(9, 11) = -0.1;
    Corr2(10, 10) = 1.;
    Corr2(10, 11) = 0.02;
    Corr2(11, 11) = 1.;

    Corr2(12, 12) = 1.;
    Corr2(12, 13) = 0.630;
    Corr2(12, 14) = -0.241;
    Corr2(12, 15) = -0.016;
    Corr2(12, 16) = 0.602;
    Corr2(12, 17) = 0.083;
    Corr2(13, 13) = 1.;
    Corr2(13, 14) = 0.008;
    Corr2(13, 15) = 0.154;
    Corr2(13, 16) = 0.735;
    Corr2(13, 17) = 0.230;
    Corr2(14, 14) = 1.;
    Corr2(14, 15) = 0.515;
    Corr2(14, 16) = 0.232;
    Corr2(14, 17) = 0.618;
    Corr2(15, 15) = 1.;
    Corr2(15, 16) = 0.237;
    Corr2(15, 17) = 0.112;
    Corr2(16, 16) = 1.;
    Corr2(16, 17) = -0.201;
    Corr2(17, 17) = 1.;

    Corr2(18,18) = 1.;
    Corr2(18,19) = 0.319;
    Corr2(18,20) = 0.030;
    Corr2(18,21) = -0.158;
    Corr2(19,19) = 1.;
    Corr2(19,20) = 0.022;
    Corr2(19,21) = -0.478;
    Corr2(20,20) = 1.;
    Corr2(20,21) = -0.055;
    Corr2(21,21) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr2);

    // Correlation Matrix (syst)
    Corr3.ResizeTo(22, 22);
    Corr3 = 0.;


    Corr3(0, 0) = 1.;
    Corr3(0, 1) = -0.047;
    Corr3(0, 2) = -0.490;
    Corr3(0, 3) = 0.322;
    Corr3(0, 4) = 0.189;
    Corr3(0, 5) = 0.144;
    Corr3(0, 6) = 0.27;
    Corr3(0, 7) = 0.07;
    Corr3(0, 8) = 0.;
    Corr3(0, 9) = 0.21;
    Corr3(0, 10) = -0.23;
    Corr3(0, 11) = -0.05;
    Corr3(0, 12) = 0.08;
    Corr3(0, 13) = 0.09;
    Corr3(0, 14) = 0.29;
    Corr3(0, 15) = 0.38;
    Corr3(0, 16) = 0.21;
    Corr3(0, 17) = 0.06;
    Corr3(0, 18) = -0.27;
    Corr3(0, 19) = -0.07;
    Corr3(0, 20) = 0.20;
    Corr3(0, 21) = -0.04;
    Corr3(1, 1) = 1.;
    Corr3(1, 2) = 0.059;
    Corr3(1, 3) = -0.237;
    Corr3(1, 4) = -0.13;
    Corr3(1, 5) = -0.117;
    Corr3(1, 6) = 0.09;
    Corr3(1, 7) = 0.04;
    Corr3(1, 8) = 0.17;
    Corr3(1, 9) = - 0.05;
    Corr3(1, 10) = 0.07;
    Corr3(1, 11) = 0.11;
    Corr3(1, 12) = 0.12;
    Corr3(1, 13) = 0.07;
    Corr3(1, 14) = 0.22;
    Corr3(1, 15) = 0.14;
    Corr3(1, 16) = 0.24;
    Corr3(1, 17) = 0.35;
    Corr3(1, 18) = -0.53;
    Corr3(1, 19) = -0.65;
    Corr3(1, 20) = -0.24;
    Corr3(1, 21) = 0.46;
    Corr3(2, 2) = 1.;
    Corr3(2, 3) = 0.061;
    Corr3(2, 4) = 0.004;
    Corr3(2, 5) = -0.139;
    Corr3(2, 6) = -0.14;
    Corr3(2, 7) = 0.04;
    Corr3(2, 8) = -0.16;
    Corr3(2, 9) = -0.07;
    Corr3(2, 10) = 0.09;
    Corr3(2, 11) = 0.05;
    Corr3(2, 12) = -0.16;
    Corr3(2, 13) = -0.17;
    Corr3(2, 14) = -0.01;
    Corr3(2, 15) = -0.06;
    Corr3(2, 16) = -0.20;
    Corr3(2, 17) = -0.07;
    Corr3(2, 18) = -0.14;
    Corr3(2, 19) = 0.;
    Corr3(2, 20) = -0.16;
    Corr3(2, 21) = 0.01;
    Corr3(3, 3) = 1.;
    Corr3(3, 4) = 0.14;
    Corr3(3, 5) = -0.199;
    Corr3(3, 6) = 0.1;
    Corr3(3, 7) = -0.03;
    Corr3(3, 8) = -0.2;
    Corr3(3, 9) = 0.09;
    Corr3(3, 10) = -0.09;
    Corr3(3, 11) = 0.;
    Corr3(3, 12) = -0.05;
    Corr3(3, 13) = -0.05;
    Corr3(3, 14) = 0.41;
    Corr3(3, 15) = -0.17;
    Corr3(3, 16) = -0.08;
    Corr3(3, 17) = 0.;
    Corr3(3, 18) = -0.42;
    Corr3(3, 19) = -0.25;
    Corr3(3, 20) = 0.27;
    Corr3(3, 21) = -0.08;
    Corr3(4, 4) = 1.;
    Corr3(4, 5) = 0.638;
    Corr3(4, 6) = 0.11;
    Corr3(4, 7) = 0.09;
    Corr3(4, 8) = 0.02;
    Corr3(4, 9) = 0.05;
    Corr3(4, 10) = -0.1;
    Corr3(4, 11) = 0.06;
    Corr3(4, 12) = -0.47;
    Corr3(4, 13) = -0.50;
    Corr3(4, 14) = -0.23;
    Corr3(4, 15) = 0.54;
    Corr3(4, 16) = -0.41;
    Corr3(4, 17) = -0.50;
    Corr3(4, 18) = 0.01;
    Corr3(4, 19) = 0.27;
    Corr3(4, 20) = -0.08;
    Corr3(4, 21) = -0.28;
    Corr3(5, 5) = 1.;
    Corr3(5, 6) = 0.08;
    Corr3(5, 7) = 0.04;
    Corr3(5, 8) = 0.16;
    Corr3(5, 9) = 0.03;
    Corr3(5, 10) = -0.05;
    Corr3(5, 11) = 0.02;
    Corr3(5, 12) = -0.37;
    Corr3(5, 13) = -0.54;
    Corr3(5, 14) = -0.13;
    Corr3(5, 15) = 0.56;
    Corr3(5, 16) = -0.42;
    Corr3(5, 17) = -0.54;
    Corr3(5, 18) = 0.30;
    Corr3(5, 19) = 0.53;
    Corr3(5, 20) = 0.32;
    Corr3(5, 21) = -0.53;


    Corr3(6, 6) = 1.;
    Corr3(6, 7) = 0.26;
    Corr3(6, 8) = 0.44;
    Corr3(6, 9) = 0.07;
    Corr3(6, 10) = -0.59;
    Corr3(6, 11) = 0.24;
    Corr3(6, 12) = 0.03;
    Corr3(6, 13) = 0.01;
    Corr3(6, 14) = 0.11;
    Corr3(6, 15) = 0.16;
    Corr3(6, 16) = 0.07;
    Corr3(6, 17) = 0.04;
    Corr3(6, 18) = -0.03;
    Corr3(6, 19) = -0.05;
    Corr3(6, 20) = 0.01;
    Corr3(6, 21) = 0.04;
    Corr3(7, 7) = 1.;
    Corr3(7, 8) = -0.21;
    Corr3(7, 9) = -0.58;
    Corr3(7, 10) = -0.53;
    Corr3(7, 11) = -0.39;
    Corr3(7, 12) = -0.02;
    Corr3(7, 13) = -0.03;
    Corr3(7, 14) = -0.04;
    Corr3(7, 15) = 0.1;
    Corr3(7, 16) = -0.02;
    Corr3(7, 17) = -0.04;
    Corr3(7, 18) = -0.05;
    Corr3(7, 19) = -0.05;
    Corr3(7, 20) = 0.01;
    Corr3(7, 21) = 0.02;
    Corr3(8, 8) = 1.;
    Corr3(8, 9) = -0.04;
    Corr3(8, 10) = 0.08;
    Corr3(8, 11) = 0.44;
    Corr3(8, 12) = 0.19;
    Corr3(8, 13) = 0.26;
    Corr3(8, 14) = -0.20;
    Corr3(8, 15) = 0.14;
    Corr3(8, 16) = 0.26;
    Corr3(8, 17) = 0.17;
    Corr3(8, 18) = -0.03;
    Corr3(8, 19) = -0.05;
    Corr3(8, 20) = 0.03;
    Corr3(8, 21) = 0.02;
    Corr3(9, 9) = 1.;
    Corr3(9, 10) = -0.24;
    Corr3(9, 11) = -0.03;
    Corr3(9, 12) = 0.04;
    Corr3(9, 13) = 0.05;
    Corr3(9, 14) = 0.05;
    Corr3(9, 15) = 0.1;
    Corr3(9, 16) = 0.06;
    Corr3(9, 17) = 0.03;
    Corr3(9, 18) = 0.02;
    Corr3(9, 19) = 0.02;
    Corr3(9, 20) = -0.02;
    Corr3(9, 21) = -0.01;
    Corr3(10, 10) = 1.;
    Corr3(10, 11) = 0.57;
    Corr3(10, 12) = 0.04;
    Corr3(10, 13) = 0.02;
    Corr3(10, 14) = 0.01;
    Corr3(10, 15) = -0.16;
    Corr3(10, 16) = 0.;
    Corr3(10, 17) = 0.05;
    Corr3(10, 18) = 0.06;
    Corr3(10, 19) = 0.06;
    Corr3(10, 20) = 0.;
    Corr3(10, 21) = -0.02;
    Corr3(11, 11) = 1.;
    Corr3(11, 12) = 0.01;
    Corr3(11, 13) = 0.01;
    Corr3(11, 14) = 0.07;
    Corr3(11, 15) = -0.04;
    Corr3(11, 16) = 0.02;
    Corr3(11, 17) = 0.06;
    Corr3(11, 18) = 0.05;
    Corr3(11, 19) = 0.03;
    Corr3(11, 20) = -0.01;
    Corr3(11, 21) = 0.01;

    Corr3(12, 12) = 1.;
    Corr3(12, 13) = 0.886;
    Corr3(12, 14) = -0.117;
    Corr3(12, 15) = -0.063;
    Corr3(12, 16) = 0.916;
    Corr3(12, 17) = 0.904;
    Corr3(12, 18) = -0.04;
    Corr3(12, 19) = -0.01;
    Corr3(12, 20) = 0.03;
    Corr3(12, 21) = 0.04;
    Corr3(13, 13) = 1.;
    Corr3(13, 14) = -0.124;
    Corr3(13, 15) = -0.115;
    Corr3(13, 16) = 0.935;
    Corr3(13, 17) = 0.889;
    Corr3(13, 18) = -0.03;
    Corr3(13, 19) = 0.;
    Corr3(13, 20) = 0.02;
    Corr3(13, 21) = 0.01;
    Corr3(14, 14) = 1.;
    Corr3(14, 15) = -0.280;
    Corr3(14, 16) = -0.130;
    Corr3(14, 17) = 0.063;
    Corr3(14, 18) = 0.01;
    Corr3(14, 19) = 0.04;
    Corr3(14, 20) = -0.04;
    Corr3(14, 21) = 0.;
    Corr3(15, 15) = 1.;
    Corr3(15, 16) = 0.085;
    Corr3(15, 17) = -0.077;
    Corr3(15, 18) = -0.02;
    Corr3(15, 19) = -0.05;
    Corr3(15, 20) = 0.01;
    Corr3(15, 21) = 0.01;
    Corr3(16, 16) = 1.;
    Corr3(16, 17) = 0.941;
    Corr3(16, 18) = -0.04;
    Corr3(16, 19) = -0.02;
    Corr3(16, 20) = 0.01;
    Corr3(16, 21) = 0.02;
    Corr3(17, 17) = 1.;
    Corr3(17, 18) = -0.04;
    Corr3(17, 19) = -0.01;
    Corr3(17, 20) = 0.;
    Corr3(17, 21) = 0.03;

    Corr3(18,18) = 1.;
    Corr3(18,19) = 0.71;
    Corr3(18,20) = -0.01;
    Corr3(18,21) = -0.02;
    Corr3(19,19) = 1.;
    Corr3(19,20) = -0.05;
    Corr3(19,21) = -0.29;
    Corr3(20,20) = 1.;
    Corr3(20,21) = -0.15;
    Corr3(21,21) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr3);


    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("GGSZ_LHCb_Cb", CorrelatedGaussianObservables(CorrData, Corr, Corr2, Corr3)));
  } // end of if comb == 0

  //----------------------------------------------------------------------------------------------------------------------------------

  //-------------------------------------------------  Bpm -> D*hpm -------------------------------------------------------------------------


  // GLW: D -> KK, D -> pipi; ADS: D -> Kpi
  // https://arxiv.org/abs/2012.09903
  // Observables 18:
  CorrData.clear();
  CorrData.push_back(dato(0.123, 0.054, 0.031));       // acp_dstk_dg
  CorrData.push_back(dato(-0.115, 0.019, 0.009));      // acp_dstk_dp
  CorrData.push_back(dato(-0.004, 0.014, 0.003));      // afav_dstk_dg
  CorrData.push_back(dato(0.02, 0.007, 0.003));        // afav_dstk_dp
  CorrData.push_back(dato(0.952, 0.062, 0.065));       // rcp_dg
  CorrData.push_back(dato(1.051, 0.022, 0.028));       // rcp_dp
  CorrData.push_back(dato(0.01170, 0.02150, 0.0313));  // rm_dstk_dg
  CorrData.push_back(dato(0.0202, 0.0035, 0.0023));    // rm_dstk_dp
  CorrData.push_back(dato(0.0292, 0.0214, 0.0312));    // rp_dstk_dg
  CorrData.push_back(dato(0.0033, 0.0035, 0.0022));    // rp_dstk_dp
  CorrData.push_back(dato(0., 0.014, 0.006));          // acp_dstpi_dg
  CorrData.push_back(dato(0.013, 0.007, 0.003));       // acp_dstpi_dp
  CorrData.push_back(dato(0.00472, 0.00092, 0.00118)); // rm_dstpi_dg
  CorrData.push_back(dato(0.00405, 0.00056, 0.00059)); // rm_dstpi_dp
  CorrData.push_back(dato(0.00403, 0.00091, 0.00114)); // rp_dstpi_dg
  CorrData.push_back(dato(0.00536, 0.00056, 0.00058)); // rp_dstpi_dp
  CorrData.push_back(dato(-0.004, 0.004, 0.001));      // afav_dstpi_dg
  CorrData.push_back(dato(0.001, 0.002, 0.001));        // afav_dstpi_dp

  // Correlation Matrix (stat):
  Corr.ResizeTo(18, 18);
  Corr(0, 0) = 1.;
  Corr(0, 1) = -0.61;
  Corr(0, 2) = 0.0;
  Corr(0, 3) = 0.0;
  Corr(0, 4) = -0.15;
  Corr(0, 5) = 0.07;
  Corr(0, 6) = -0.03;
  Corr(0, 7) = 0.01;
  Corr(0, 8) = -0.03;
  Corr(0, 9) = 0.01;
  Corr(0, 10) = -0.01;
  Corr(0, 11) = -0.02;
  Corr(0, 12) = -0.01;
  Corr(0, 13) = 0.;
  Corr(0, 14) = -0.02;
  Corr(0, 15) = 0.;
  Corr(0, 16) = 0.01;
  Corr(0, 17) = 0.01;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.0;
  Corr(1, 3) = 0.01;
  Corr(1, 4) = -0.05;
  Corr(1, 5) = 0.08;
  Corr(1, 6) = 0.05;
  Corr(1, 7) = 0.;
  Corr(1, 8) = 0.;
  Corr(1, 9) = 0.;
  Corr(1, 10) = 0.03;
  Corr(1, 11) = 0.05;
  Corr(1, 12) = 0.;
  Corr(1, 13) = 0.;
  Corr(1, 14) = 0.;
  Corr(1, 15) = 0.;
  Corr(1, 16) = 0.02;
  Corr(1, 17) = 0.03;
  Corr(2, 2) = 1.;
  Corr(2, 3) = -0.59;
  Corr(2, 4) = 0.;
  Corr(2, 5) = 0.;
  Corr(2, 6) = 0.0;
  Corr(2, 7) = 0.0;
  Corr(2, 8) = 0.;
  Corr(2, 9) = 0.;
  Corr(2, 10) = 0.;
  Corr(2, 11) = 0.01;
  Corr(2, 12) = 0.;
  Corr(2, 13) = 0.;
  Corr(2, 14) = 0.;
  Corr(2, 15) = 0.;
  Corr(2, 16) = 0.;
  Corr(2, 17) = 0.01;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.;
  Corr(3, 5) = 0.;
  Corr(3, 6) = -0.02;
  Corr(3, 7) = 0.;
  Corr(3, 8) = 0.;
  Corr(3, 9) = 0.;
  Corr(3, 10) = 0.01;
  Corr(3, 11) = 0.02;
  Corr(3, 12) = 0.;
  Corr(3, 13) = 0.01;
  Corr(3, 14) = 0.;
  Corr(3, 15) = -0.01;
  Corr(3, 16) = 0.05;
  Corr(3, 17) = 0.08;
  Corr(4, 4) = 1.;
  Corr(4, 5) = -0.44;
  Corr(4, 6) = 0.23;
  Corr(4, 7) = -0.08;
  Corr(4, 8) = 0.23;
  Corr(4, 9) = -0.08;
  Corr(4, 10) = 0.;
  Corr(4, 11) = 0.;
  Corr(4, 12) = 0.10;
  Corr(4, 13) = 0.;
  Corr(4, 14) = 0.1;
  Corr(4, 15) = 0.;
  Corr(4, 16) = 0.;
  Corr(4, 17) = 0.;
  Corr(5, 5) = 1.;
  Corr(5, 6) = -0.04;
  Corr(5, 7) = 0.03;
  Corr(5, 8) = -0.04;
  Corr(5, 9) = 0.02;
  Corr(5, 10) = 0.;
  Corr(5, 11) = 0.;
  Corr(5, 12) = -0.02;
  Corr(5, 13) = 0.;
  Corr(5, 14) = -0.02;
  Corr(5, 15) = 0.;
  Corr(5, 16) = 0.;
  Corr(5, 17) = 0.;
  Corr(6, 6) = 1.;
  Corr(6, 7) = -0.59;
  Corr(6, 8) = 0.79;
  Corr(6, 9) = -0.27;
  Corr(6, 10) = 0.;
  Corr(6, 11) = 0.;
  Corr(6, 12) = 0.3;
  Corr(6, 13) = -0.03;
  Corr(6, 14) = 0.33;
  Corr(6, 15) = -0.02;
  Corr(6, 16) = 0.01;
  Corr(6, 17) = 0.01;
  Corr(7, 7) = 1.;
  Corr(7, 8) = -0.27;
  Corr(7, 9) = 0.1;
  Corr(7, 10) = 0.;
  Corr(7, 11) = 0.;
  Corr(7, 12) = -0.07;
  Corr(7, 13) = 0.05;
  Corr(7, 14) = -0.1;
  Corr(7, 15) = 0.05;
  Corr(7, 16) = 0.0;
  Corr(7, 17) = 0.0;
  Corr(8, 8) = 1.;
  Corr(8, 9) = -0.6;
  Corr(8, 10) = 0.;
  Corr(8, 11) = 0.;
  Corr(8, 12) = 0.32;
  Corr(8, 13) = -0.01;
  Corr(8, 14) = 0.3;
  Corr(8, 15) = -0.03;
  Corr(8, 16) = -0.01;
  Corr(8, 17) = -0.01;
  Corr(9, 9) = 1.;
  Corr(9, 10) = 0.;
  Corr(9, 11) = 0.;
  Corr(9, 12) = -0.09;
  Corr(9, 13) = 0.04;
  Corr(9, 14) = -0.07;
  Corr(9, 15) = 0.05;
  Corr(9, 16) = 0.0;
  Corr(9, 17) = 0.0;
  Corr(10, 10) = 1.;
  Corr(10, 11) = 0.05;
  Corr(10, 12) = 0.;
  Corr(10, 13) = 0.;
  Corr(10, 14) = 0.;
  Corr(10, 15) = 0.;
  Corr(10, 16) = 0.02;
  Corr(10, 17) = 0.03;
  Corr(11, 11) = 1.;
  Corr(11, 12) = 0.;
  Corr(11, 13) = 0.;
  Corr(11, 14) = 0.;
  Corr(11, 15) = -0.01;
  Corr(11, 16) = 0.05;
  Corr(11, 17) = 0.08;
  Corr(12, 12) = 1.;
  Corr(12, 13) = -0.11;
  Corr(12, 14) = 0.33;
  Corr(12, 15) = 0.19;
  Corr(12, 16) = 0.01;
  Corr(12, 17) = 0.01;
  Corr(13, 13) = 1.;
  Corr(13, 14) = 0.19;
  Corr(13, 15) = 0.58;
  Corr(13, 16) = 0.01;
  Corr(13, 17) = 0.02;
  Corr(14, 14) = 1.;
  Corr(14, 15) = -0.11;
  Corr(14, 16) = -0.01;
  Corr(14, 17) = -0.01;
  Corr(15, 15) = 1.;
  Corr(15, 16) = -0.01;
  Corr(15, 17) = -0.02;
  Corr(16, 16) = 1.;
  Corr(16, 17) = -0.28;
  Corr(17,17) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(18, 18);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = -0.52;
  Corr2(0, 2) = 0.50;
  Corr2(0, 3) = -0.39;
  Corr2(0, 4) = -0.31;
  Corr2(0, 5) = -0.71;
  Corr2(0, 6) = -0.15;
  Corr2(0, 7) = 0.10;
  Corr2(0, 8) = -0.05;
  Corr2(0, 9) = 0.12;
  Corr2(0, 10) = 0.74;
  Corr2(0, 11) = -0.02;
  Corr2(0, 12) = -0.16;
  Corr2(0, 13) = -0.14;
  Corr2(0, 14) = -0.05;
  Corr2(0, 15) = -0.17;
  Corr2(0, 16) = -0.70;
  Corr2(0, 17) = .05;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = -0.44;
  Corr2(1, 3) = -0.08;
  Corr2(1, 4) = 0.33;
  Corr2(1, 5) = 0.62;
  Corr2(1, 6) = 0.19;
  Corr2(1, 7) = -0.19;
  Corr2(1, 8) = 0.16;
  Corr2(1, 9) = -0.08;
  Corr2(1, 10) = -0.49;
  Corr2(1, 11) = 0.;
  Corr2(1, 12) = 0.18;
  Corr2(1, 13) = 0.01;
  Corr2(1, 14) = 0.1;
  Corr2(1, 15) = 0.03;
  Corr2(1, 16) = .53;
  Corr2(1, 17) = -.19;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.23;
  Corr2(2, 4) = 0.02;
  Corr2(2, 5) = -0.35;
  Corr2(2, 6) = -0.04;
  Corr2(2, 7) = 0.03;
  Corr2(2, 8) = 0.04;
  Corr2(2, 9) = 0.05;
  Corr2(2, 10) = 0.57;
  Corr2(2, 11) = -0.38;
  Corr2(2, 12) = -0.15;
  Corr2(2, 13) = -0.19;
  Corr2(2, 14) = 0.04;
  Corr2(2, 15) = -0.09;
  Corr2(2, 16) = -.38;
  Corr2(2, 17) = .32;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = -0.19;
  Corr2(3, 5) = 0.15;
  Corr2(3, 6) = -0.08;
  Corr2(3, 7) = 0.06;
  Corr2(3, 8) = -0.12;
  Corr2(3, 9) = -0.01;
  Corr2(3, 10) = -0.54;
  Corr2(3, 11) = -0.64;
  Corr2(3, 12) = -0.05;
  Corr2(3, 13) = 0.04;
  Corr2(3, 14) = 0.02;
  Corr2(3, 15) = 0.24;
  Corr2(3, 16) = .42;
  Corr2(3, 17) = .28;
  Corr2(4, 4) = 1.;
  Corr2(4, 5) = 0.49;
  Corr2(4, 6) = 0.36;
  Corr2(4, 7) = -0.15;
  Corr2(4, 8) = 0.41;
  Corr2(4, 9) = -0.09;
  Corr2(4, 10) = 0.22;
  Corr2(4, 11) = 0.16;
  Corr2(4, 12) = 0.16;
  Corr2(4, 13) = -0.06;
  Corr2(4, 14) = 0.19;
  Corr2(4, 15) = -0.07;
  Corr2(4, 16) = .07;
  Corr2(4, 17) = 0.;
  Corr2(5, 5) = 1.;
  Corr2(5, 6) = 0.14;
  Corr2(5, 7) = -0.04;
  Corr2(5, 8) = 0.10;
  Corr2(5, 9) = -0.05;
  Corr2(5, 10) = -0.35;
  Corr2(5, 11) = 0.07;
  Corr2(5, 12) = 0.15;
  Corr2(5, 13) = 0.08;
  Corr2(5, 14) = 0.09;
  Corr2(5, 15) = 0.1;
  Corr2(5, 16) = .42;
  Corr2(5, 17) = -0.11;
  Corr2(6, 6) = 1.;
  Corr2(6, 7) = -0.55;
  Corr2(6, 8) = 0.98;
  Corr2(6, 9) = -0.51;
  Corr2(6, 10) = 0.;
  Corr2(6, 11) = 0.07;
  Corr2(6, 12) = 0.37;
  Corr2(6, 13) = -0.04;
  Corr2(6, 14) = 0.4;
  Corr2(6, 15) = -0.04;
  Corr2(6, 16) = 0.1;
  Corr2(6, 17) = .0;
  Corr2(7, 7) = 1.;
  Corr2(7, 8) = -0.53;
  Corr2(7, 9) = 0.94;
  Corr2(7, 10) = 0.04;
  Corr2(7, 11) = -0.01;
  Corr2(7, 12) = -0.13;
  Corr2(7, 13) = 0.08;
  Corr2(7, 14) = -0.15;
  Corr2(7, 15) = 0.06;
  Corr2(7, 16) = -0.1;
  Corr2(7, 17) = 0.;
  Corr2(8, 8) = 1.;
  Corr2(8, 9) = -0.48;
  Corr2(8, 10) = 0.12;
  Corr2(8, 11) = 0.04;
  Corr2(8, 12) = 0.37;
  Corr2(8, 13) = -0.06;
  Corr2(8, 14) = 0.41;
  Corr2(8, 15) = -0.06;
  Corr2(8, 16) = 0.03;
  Corr2(8, 17) = 0.02;
  Corr2(9, 9) = 1.;
  Corr2(9, 10) = 0.06;
  Corr2(9, 11) = -0.03;
  Corr2(9, 12) = -0.14;
  Corr2(9, 13) = 0.02;
  Corr2(9, 14) = -0.13;
  Corr2(9, 15) = 0.04;
  Corr2(9, 16) = -0.1;
  Corr2(9, 17) = 0.0;
  Corr2(10, 10) = 1.;
  Corr2(10, 11) = 0.26;
  Corr2(10, 12) = -0.1;
  Corr2(10, 13) = -0.18;
  Corr2(10, 14) = 0.02;
  Corr2(10, 15) = -0.25;
  Corr2(10, 16) = -0.67;
  Corr2(10, 17) = 0.17;
  Corr2(11, 11) = 1.;
  Corr2(11, 12) = 0.09;
  Corr2(11, 13) = 0.08;
  Corr2(11, 14) = -0.07;
  Corr2(11, 15) = -0.14;
  Corr2(11, 16) = -0.02;
  Corr2(11, 17) = 0.06;
  Corr2(12, 12) = 1.;
  Corr2(12, 13) = 0.53;
  Corr2(12, 14) = 0.79;
  Corr2(12, 15) = 0.16;
  Corr2(12, 16) = 0.13;
  Corr2(12, 17) = -0.04;
  Corr2(13, 13) = 1.;
  Corr2(13, 14) = 0.22;
  Corr2(13, 15) = 0.60;
  Corr2(13, 16) = 0.09;
  Corr2(13, 17) = -0.05;
  Corr2(14, 14) = 1.;
  Corr2(14, 15) = 0.41;
  Corr2(14, 16) = 0.06;
  Corr2(14, 17) = .02;
  Corr2(15, 15) = 1.;
  Corr2(15, 16) = 0.17;
  Corr2(15, 17) = 0.;
  Corr2(16, 16) = 1.;
  Corr2(16, 17) = 0.43;
  Corr2(17, 17) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID5", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // GLW: D -> KK, D -> pipi, D -> K0spi0, .... (Belle)
  // https://arxiv.org/pdf/hep-ex/0601032
  // Observables 8:
  CorrData.clear();
  CorrData.push_back(dato(0.13, 0.30, 0.08));  // acp_dstk_dp_CPm
  CorrData.push_back(dato(1.15, 0.31, 0.12));  // rcp_dstk_dp_CPm
  CorrData.push_back(dato(-0.20, 0.22, 0.04)); // acp_dstk_dp_CPp
  CorrData.push_back(dato(1.41, 0.25, 0.06));  // rcp_dstk_dp_CPp

  // Correlation Matrix (stat):
  Corr.ResizeTo(4, 4);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.;
  Corr(0, 2) = 0.;
  Corr(0, 3) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.;
  Corr(1, 3) = 0.;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.;
  Corr(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (stat):
  Corr2.ResizeTo(4, 4);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.;
  Corr2(0, 2) = 0.;
  Corr2(0, 3) = 0.;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.;
  Corr2(1, 3) = 0.;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.;
  Corr2(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Belle_PRD73_2006", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // GGSZ: D -> K0spipi (Belle)
  // https://arxiv.org/abs/1003.3360
  // Observables 8:
  CorrData.clear();
  CorrData.push_back(dato(0.024, 0.140, 0.018));  // xm_dstk_dp
  CorrData.push_back(dato(-0.243, 0.137, 0.022)); // ym_dstk_dp
  CorrData.push_back(dato(0.133, 0.083, 0.018));  // xp_dstk_dp
  CorrData.push_back(dato(0.130, 0.120, 0.022));  // yp_dstk_dp
  CorrData.push_back(dato(0.144, 0.208, 0.025));  // xm_dstk_dg
  CorrData.push_back(dato(0.196, 0.215, 0.037));  // ym_dstk_dg
  CorrData.push_back(dato(-0.006, 0.147, 0.025)); // xp_dstk_dg
  CorrData.push_back(dato(-0.190, 0.177, 0.037)); // yp_dstk_dg

  // Correlation Matrix (stat):
  Corr.ResizeTo(8, 8);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.440;
  Corr(0, 2) = 0.;
  Corr(0, 3) = 0.;
  Corr(0, 4) = 0.;
  Corr(0, 5) = 0.;
  Corr(0, 6) = 0.;
  Corr(0, 7) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.;
  Corr(1, 3) = 0.;
  Corr(1, 4) = 0.;
  Corr(1, 5) = 0.;
  Corr(1, 6) = 0.;
  Corr(1, 7) = 0.;
  Corr(2, 2) = 1.;
  Corr(2, 3) = -0.101;
  Corr(2, 4) = 0.;
  Corr(2, 5) = 0.;
  Corr(2, 6) = 0.;
  Corr(2, 7) = 0.;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.;
  Corr(3, 5) = 0.;
  Corr(3, 6) = 0.;
  Corr(3, 7) = 0.;
  Corr(4, 4) = 1.;
  Corr(4, 5) = -0.207;
  Corr(4, 6) = 0.;
  Corr(4, 7) = 0.;
  Corr(5, 5) = 1.;
  Corr(5, 6) = 0.;
  Corr(5, 7) = 0.;
  Corr(6, 6) = 1.;
  Corr(6, 7) = 0.08;
  Corr(7, 7) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (stat):
  Corr2.ResizeTo(8, 8);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.;
  Corr2(0, 2) = 0.;
  Corr2(0, 3) = 0.;
  Corr2(0, 4) = 0.;
  Corr2(0, 5) = 0.;
  Corr2(0, 6) = 0.;
  Corr2(0, 7) = 0.;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.;
  Corr2(1, 3) = 0.;
  Corr2(1, 4) = 0.;
  Corr2(1, 5) = 0.;
  Corr2(1, 6) = 0.;
  Corr2(1, 7) = 0.;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.;
  Corr2(2, 4) = 0.;
  Corr2(2, 5) = 0.;
  Corr2(2, 6) = 0.;
  Corr2(2, 7) = 0.;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = 0.;
  Corr2(3, 5) = 0.;
  Corr2(3, 6) = 0.;
  Corr2(3, 7) = 0.;
  Corr2(4, 4) = 1.;
  Corr2(4, 5) = 0.;
  Corr2(4, 6) = 0.;
  Corr2(4, 7) = 0.;
  Corr2(5, 5) = 1.;
  Corr2(5, 6) = 0.;
  Corr2(5, 7) = 0.;
  Corr2(6, 6) = 1.;
  Corr2(6, 7) = 0.;
  Corr2(7, 7) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("Belle_PRD81_2010", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //----------------------------------------------------------------------------------------------------------------------------------

  //-------------------------------------------------  Bpm -> DKstpm -------------------------------------------------------------------------

  // GLW: D -> KK, D -> pipi, D -> 4pi; ADS: D -> Kpi, D -> K3pi
  // LHCb-PAPER-2024-023
  // Observables 12:
  CorrData.clear();
  CorrData.push_back(dato(-0.024, 0.014, 0.002));  // afav_dkst_kpi
  CorrData.push_back(dato(0.14, 0.04, 0.001));     // acp_dkst_kk
  CorrData.push_back(dato(0.31, 0.09, 0.05));     // acp_dkst_pipi
  CorrData.push_back(dato(-0.73, 0.16, 0.03));     // asup_dkst_kpi
  CorrData.push_back(dato(1.10, 0.05, 0.01));     // rcp_dkst_kk
  CorrData.push_back(dato(0.96, 0.09, 0.05));     // rcp_dkst_pipi
  CorrData.push_back(dato(0.0098, 0.0019, 0.0002));   // rsup_dkst_kpi
  CorrData.push_back(dato(-0.024, 0.018, 0.002)); // afav_dkst_k3pi
  CorrData.push_back(dato(0.04, 0.06, 0.03));    // acp_dkst_pipipipi
  CorrData.push_back(dato(-.19, 0.22, 0.01));     // asup_dkst_k3pi
  CorrData.push_back(dato(1.05, 0.07, 0.05));  // rcp_dkst_pipipipi
  CorrData.push_back(dato(0.0118, 0.0026, 0.0001));  // rsup_dkst_k3pi

  // Correlation Matrix (stat):
  Corr.ResizeTo(12, 12);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.;
  Corr(0, 2) = 0.0;
  Corr(0, 3) = 0.0;
  Corr(0, 4) = 0.0;
  Corr(0, 5) = 0.0;
  Corr(0, 6) = 0.0;
  Corr(0, 7) = 0.;
  Corr(0, 8) = 0.;
  Corr(0, 9) = 0.;
  Corr(0, 10) = 0.;
  Corr(0, 11) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.0;
  Corr(1, 3) = 0.0;
  Corr(1, 4) = -0.019;
  Corr(1, 5) = 0.0;
  Corr(1, 6) = 0.0;
  Corr(1, 7) = 0.;
  Corr(1, 8) = 0.;
  Corr(1, 9) = 0.;
  Corr(1, 10) = 0.;
  Corr(1, 11) = 0.;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.;
  Corr(2, 4) = 0.0;
  Corr(2, 5) = -0.085;
  Corr(2, 6) = 0.0;
  Corr(2, 7) = 0.0;
  Corr(2, 8) = 0.;
  Corr(2, 9) = 0.;
  Corr(2, 10) = 0.0;
  Corr(2, 11) = 0.;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.0;
  Corr(3, 5) = 0.0;
  Corr(3, 6) = 0.183;
  Corr(3, 7) = 0.;
  Corr(3, 8) = 0.;
  Corr(3, 9) = 0.;
  Corr(3, 10) = 0.;
  Corr(3, 11) = 0.;
  Corr(4, 4) = 1.;
  Corr(4, 5) = 0.055;
  Corr(4, 6) = 0.028;
  Corr(4, 7) = 0.;
  Corr(4, 8) = 0.;
  Corr(4, 9) = 0.;
  Corr(4, 10) = 0.;
  Corr(4, 11) = 0.;
  Corr(5, 5) = 1.;
  Corr(5, 6) = 0.02;
  Corr(5, 7) = 0.;
  Corr(5, 8) = 0.;
  Corr(5, 9) = 0.;
  Corr(5, 10) = 0.;
  Corr(5, 11) = 0.;
  Corr(6, 6) = 1.;
  Corr(6, 7) = 0.;
  Corr(6, 8) = 0.;
  Corr(6, 9) = 0.;
  Corr(6, 10) = 0.;
  Corr(6, 11) = 0.;
  Corr(7, 7) = 1.;
  Corr(7, 8) = 0.0;
  Corr(7, 9) = 0.;
  Corr(7, 10) = 0.0;
  Corr(7, 11) = 0.;
  Corr(8, 8) = 1.;
  Corr(8, 9) = 0.0;
  Corr(8, 10) = -0.22;
  Corr(8, 11) = 0.;
  Corr(9, 9) = 1.;
  Corr(9, 10) = 0.0;
  Corr(9, 11) = 0.04;
  Corr(10, 10) = 1.;
  Corr(10, 11) = 0.029;
  Corr(11, 11) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(12, 12);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.314;
  Corr2(0, 2) = 0.0;
  Corr2(0, 3) = -0.014;
  Corr2(0, 4) = 0.;
  Corr2(0, 5) = 0.0;
  Corr2(0, 6) = -0.016;
  Corr2(0, 7) = 0.918;
  Corr2(0, 8) = 0.01;
  Corr2(0, 9) = -0.169;
  Corr2(0, 10) = 0.0;
  Corr2(0, 11) = 0.;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.;
  Corr2(1, 3) = -0.017;
  Corr2(1, 4) = -0.039;
  Corr2(1, 5) = 0.0;
  Corr2(1, 6) = -0.011;
  Corr2(1, 7) = 0.319;
  Corr2(1, 8) = 0.032;
  Corr2(1, 9) = 0.155;
  Corr2(1, 10) = -0.019;
  Corr2(1, 11) = 0.048;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.;
  Corr2(2, 4) = 0.0;
  Corr2(2, 5) = -0.292;
  Corr2(2, 6) = 0.0;
  Corr2(2, 7) = 0.;
  Corr2(2, 8) = 0.;
  Corr2(2, 9) = 0.;
  Corr2(2, 10) = 0.;
  Corr2(2, 11) = 0.0;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = 0.0;
  Corr2(3, 5) = 0.0;
  Corr2(3, 6) = 0.220;
  Corr2(3, 7) = -0.017;
  Corr2(3, 8) = 0.;
  Corr2(3, 9) = 0.086;
  Corr2(3, 10) = 0.0;
  Corr2(3, 11) = 0.023;
  Corr2(4, 4) = 1.;
  Corr2(4, 5) = 0.051;
  Corr2(4, 6) = 0.0;
  Corr2(4, 7) = 0.0;
  Corr2(4, 8) = 0.0;
  Corr2(4, 9) = 0.0;
  Corr2(4, 10) = 0.0;
  Corr2(4, 11) = 0.0;
  Corr2(5, 5) = 1.;
  Corr2(5, 6) = 0.0;
  Corr2(5, 7) = 0.0;
  Corr2(5, 8) = 0.;
  Corr2(5, 9) = 0.;
  Corr2(5, 10) = 0.0;
  Corr(5, 11) = 0.0;
  Corr2(6, 6) = 1.;
  Corr2(6, 7) = 0.0;
  Corr2(6, 8) = 0.0;
  Corr2(6, 9) = -0.013;
  Corr2(6, 10) = 0.0;
  Corr2(6, 11) = 0.063;
  Corr2(7, 7) = 1.;
  Corr2(7, 8) = 0.013;
  Corr2(7, 9) = -0.173;
  Corr2(7, 10) = 0.0;
  Corr2(7, 11) = 0.013;
  Corr2(8, 8) = 1.;
  Corr2(8, 9) = 0.0;
  Corr2(8, 10) = -0.605;
  Corr2(8, 11) = 0.025;
  Corr2(9, 9) = 1.;
  Corr2(9, 10) = 0.0;
  Corr2(9, 11) = 0.345;
  Corr2(10, 10) = 1.;
  Corr2(10, 11) = 0.0;
  Corr2(11, 11) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("LHCB-PAPER-2024-023-GLWADS", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  // Coherence factor
  // https://arxiv.org/pdf/1709.05855
  meas.insert(pair<string, dato>("UID24", dato(0.95, 0.06, 0.))); // k_dkst


  // Bpm -> DK^*pm
  // https://arxiv.org/pdf/hep-ex/0604054 (Belle)
  // Observables 4:
  CorrData.clear();
  CorrData.push_back(dato(-0.10, 0.172, 0.006, 0.088));  // xp_dkst
  CorrData.push_back(dato(0., 0.16, 0.013, 0.095)); // yp_dkst
  CorrData.push_back(dato(-0.807, 0.272, 0.029, 0.097)); // xm_dkst
  CorrData.push_back(dato(-0.228, 0.387,0.046, 0.086 )); // ym_dkst

  // Correlation Matrix (stat):
  Corr.ResizeTo(4, 4);
  Corr(0, 0) = 1.;
  Corr(0, 1) = -0.201;
  Corr(0, 2) = 0.;
  Corr(0, 3) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.;
  Corr(1, 3) = 0.;
  Corr(2, 2) = 1.;
  Corr(2, 3) = -0.338;
  Corr(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(4, 4);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.;
  Corr2(0, 2) = 0.;
  Corr2(0, 3) = 0.;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.;
  Corr2(1, 3) = 0.;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.;
  Corr2(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  // Correlation Matrix (ext)
  Corr3.ResizeTo(4, 4);
  Corr3(0, 0) = 1.;
  Corr3(0, 1) = -0.103;
  Corr3(0, 2) = -0.386;
  Corr3(0, 3) = -0.741;
  Corr3(1, 1) = 1.;
  Corr3(1, 2) = 0.953;
  Corr3(1, 3) = -0.121;
  Corr3(2, 2) = 1.;
  Corr3(2, 3) = 0.044;
  Corr3(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr3);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("0604054", CorrelatedGaussianObservables(CorrData, Corr, Corr2, Corr3)));


  //----------------------------------------------------------------------------------------------------------------------------------

  //-------------------------------------------------  Bpm -> Dhpipi -------------------------------------------------------------------------

  // GLW: D -> KK, D -> pipi; ADS: D -> Kpi
  // https://arxiv.org/pdf/1505.07044
  // Observables 11:
  CorrData.clear();
  CorrData.push_back(dato(1.0400, 0.064, 0.0));        // rcp_dkpipi
  CorrData.push_back(dato(0.013, 0.019, 0.013));       // afav_dkpipi_kpi
  CorrData.push_back(dato(-0.002, 0.003, 0.011));      // afav_dpipipi_kpi
  CorrData.push_back(dato(-0.045, 0.064, 0.011));      // acp_dkpipi_kk
  CorrData.push_back(dato(-0.054, 0.101, 0.011));      // acp_dkpipi_pipi
  CorrData.push_back(dato(-0.019, 0.011, 0.01));       // acp_dpipipi_kk
  CorrData.push_back(dato(-0.013, 0.016, 0.01));       // acp_dpipipi_pipi
  CorrData.push_back(dato(0.0115, 0.0052, 0.00110));    // rp_dkpipi
  CorrData.push_back(dato(0.00545, 0.0043, 0.0006));   // rm_dkpipi
  CorrData.push_back(dato(0.00432, 0.00053, 0.00021)); // rp_dpipipi
  CorrData.push_back(dato(0.00421, 0.00053, 0.00021)); // rm_dpipipi

  // Correlation Matrix (stat):
  Corr.ResizeTo(11, 11);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.;
  Corr(0, 2) = 0.0;
  Corr(0, 3) = 0.0;
  Corr(0, 4) = 0.0;
  Corr(0, 5) = 0.0;
  Corr(0, 6) = 0.0;
  Corr(0, 7) = 0.;
  Corr(0, 8) = 0.;
  Corr(0, 9) = 0.;
  Corr(0, 10) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.0;
  Corr(1, 3) = 0.0;
  Corr(1, 4) = 0.0;
  Corr(1, 5) = 0.0;
  Corr(1, 6) = 0.0;
  Corr(1, 7) = 0.;
  Corr(1, 8) = 0.;
  Corr(1, 9) = 0.;
  Corr(1, 10) = 0.;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.;
  Corr(2, 4) = 0.;
  Corr(2, 5) = 0.0;
  Corr(2, 6) = 0.0;
  Corr(2, 7) = 0.0;
  Corr(2, 8) = 0.;
  Corr(2, 9) = 0.;
  Corr(2, 10) = 0.;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.2;
  Corr(3, 5) = 0.;
  Corr(3, 6) = 0.;
  Corr(3, 7) = 0.;
  Corr(3, 8) = 0.;
  Corr(3, 9) = 0.;
  Corr(3, 10) = 0.;
  Corr(4, 4) = 1.;
  Corr(4, 5) = 0.0;
  Corr(4, 6) = 0.0;
  Corr(4, 7) = 0.;
  Corr(4, 8) = 0.;
  Corr(4, 9) = 0.;
  Corr(4, 10) = 0.;
  Corr(5, 5) = 1.;
  Corr(5, 6) = 0.08;
  Corr(5, 7) = 0.;
  Corr(5, 8) = 0.;
  Corr(5, 9) = 0.;
  Corr(5, 10) = 0.;
  Corr(6, 6) = 1.;
  Corr(6, 7) = 0.;
  Corr(6, 8) = 0.;
  Corr(6, 9) = 0.;
  Corr(6, 10) = 0.;
  Corr(7, 7) = 1.;
  Corr(7, 8) = 0.0;
  Corr(7, 9) = 0.;
  Corr(7, 10) = 0.0;
  Corr(8, 8) = 1.;
  Corr(8, 9) = 0.;
  Corr(8, 10) = 0.;
  Corr(9, 9) = 1.;
  Corr(9, 10) = 0.0;
  Corr(10, 10) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(11, 11);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.;
  Corr2(0, 2) = 0.;
  Corr2(0, 3) = 0.;
  Corr2(0, 4) = 0.;
  Corr2(0, 5) = 0.;
  Corr2(0, 6) = 0.;
  Corr2(0, 7) = 0.;
  Corr2(0, 8) = 0.;
  Corr2(0, 9) = 0.;
  Corr2(0, 10) = 0.;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.;
  Corr2(1, 3) = 0.;
  Corr2(1, 4) = 0.;
  Corr2(1, 5) = 0.;
  Corr2(1, 6) = 0.;
  Corr2(1, 7) = 0.;
  Corr2(1, 8) = 0.;
  Corr2(1, 9) = 0.;
  Corr2(1, 10) = 0.;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.;
  Corr2(2, 4) = 0.;
  Corr2(2, 5) = 0.;
  Corr2(2, 6) = 0.;
  Corr2(2, 7) = 0.;
  Corr2(2, 8) = 0.;
  Corr2(2, 9) = 0.;
  Corr2(2, 10) = 0.;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = 0.;
  Corr2(3, 5) = 0.;
  Corr2(3, 6) = 0.;
  Corr2(3, 7) = 0.;
  Corr2(3, 8) = 0.;
  Corr2(3, 9) = 0.;
  Corr2(3, 10) = 0.;
  Corr2(4, 4) = 1.;
  Corr2(4, 5) = 0.;
  Corr2(4, 6) = 0.;
  Corr2(4, 7) = 0.;
  Corr2(4, 8) = 0.;
  Corr2(4, 9) = 0.;
  Corr2(4, 10) = 0.;
  Corr2(5, 5) = 1.;
  Corr2(5, 6) = 0.;
  Corr2(5, 7) = 0.;
  Corr2(5, 8) = 0.;
  Corr2(5, 9) = 0.;
  Corr2(5, 10) = 0.;
  Corr2(6, 6) = 1.;
  Corr2(6, 7) = 0.;
  Corr2(6, 8) = 0.;
  Corr2(6, 9) = 0.;
  Corr2(6, 10) = 0.;
  Corr2(7, 7) = 1.;
  Corr2(7, 8) = 0.;
  Corr2(7, 9) = 0.;
  Corr2(7, 10) = 0.;
  Corr2(8, 8) = 1.;
  Corr2(8, 9) = 0.;
  Corr2(8, 10) = 0.;
  Corr2(9, 9) = 1.;
  Corr2(9, 10) = 0.;
  Corr2(10, 10) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID9", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //----------------------------------------------------------------------------------------------------------------------------------

}
// ---------------------------------------------------------

// ---------------------------------------------------------
void MixingModel::Add_NeutralBd_meas()
{

  vector<dato> CorrData;
  TMatrixDSym Corr, Corr2, Corr3;

  //-------------------------------------- B^0 -> DKst0 -------------------------------------------------------------------------

  if (comb == 1)
  { // When treating separately from the Bs counterparts

    // GLW: D -> pipi, D -> KK, D -> 4pi; ADS: D -> Kpi, D -> K3pi
    // https://arxiv.org/pdf/2401.17934
    // Observables 12:
    CorrData.clear();
    CorrData.push_back(dato(0.031, 0.017, 0.015));  // afav_dkstz_kpi
    CorrData.push_back(dato(0.069, 0.013, 0.005));  // rp_dkstz_kpi
    CorrData.push_back(dato(0.093, 0.013, 0.005));  // rm_dkstz_kpi
    CorrData.push_back(dato(-0.012, 0.018, 0.016)); // afav_dkstz_k3pi
    CorrData.push_back(dato(0.060, 0.014, 0.006));  // rp_dkstz_k3pi
    CorrData.push_back(dato(0.038, 0.014, 0.006));  // rm_dkstz_k3pi
    CorrData.push_back(dato(-0.047, 0.063, 0.015)); // acp_dkstz_kk
    CorrData.push_back(dato(0.811, 0.057, 0.017));  // rcp_dkstz_kk
    CorrData.push_back(dato(-0.034, 0.094, 0.016)); // acp_dkstz_pipi
    CorrData.push_back(dato(1.104, 0.111, 0.026));  // rcp_dkstz_pipi
    CorrData.push_back(dato(0.021, 0.087, 0.016));  // acp_dkstz_4pi
    CorrData.push_back(dato(0.882, 0.086, 0.033));  // rcp_dkstz_4pi


    // Correlation Matrix (stat):
    Corr.ResizeTo(12, 12);
    Corr(0, 0) = 1.;
    Corr(0, 1) = 0.09;
    Corr(0, 2) = -0.12;
    Corr(0, 3) = 0.0;
    Corr(0, 4) = 0.0;
    Corr(0, 5) = 0.;
    Corr(0, 6) = 0.;
    Corr(0, 7) = -0.01;
    Corr(0, 8) = 0.;
    Corr(0, 9) = 0.;
    Corr(0, 10) = 0.;
    Corr(0, 11) = 0.;
    Corr(1, 1) = 1.;
    Corr(1, 2) = 0.06;
    Corr(1, 3) = 0.;
    Corr(1, 4) = 0.0;
    Corr(1, 5) = 0.0;
    Corr(1, 6) = 0.0;
    Corr(1, 7) = 0.03;
    Corr(1, 8) = 0.;
    Corr(1, 9) = 0.02;
    Corr(1, 10) = 0.;
    Corr(1, 11) = 0.;
    Corr(2, 2) = 1.;
    Corr(2, 3) = 0.;
    Corr(2, 4) = 0.;
    Corr(2, 5) = 0.0;
    Corr(2, 6) = 0.0;
    Corr(2, 7) = 0.03;
    Corr(2, 8) = 0.0;
    Corr(2, 9) = 0.02;
    Corr(2, 10) = 0.;
    Corr(2, 11) = 0.;
    Corr(3, 3) = 1.;
    Corr(3, 4) = 0.07;
    Corr(3, 5) = -0.05;
    Corr(3, 6) = 0.0;
    Corr(3, 7) = 0.0;
    Corr(3, 8) = 0.0;
    Corr(3, 9) = 0.0;
    Corr(3, 10) = 0.;
    Corr(3, 11) = -0.01;
    Corr(4, 4) = 1.;
    Corr(4, 5) = 0.08;
    Corr(4, 6) = 0.;
    Corr(4, 7) = 0.;
    Corr(4, 8) = 0.;
    Corr(4, 9) = 0.;
    Corr(4, 10) = 0.;
    Corr(4, 11) = 0.02;
    Corr(5, 5) = 1.;
    Corr(5, 6) = 0.0;
    Corr(5, 7) = 0.0;
    Corr(5, 8) = 0.0;
    Corr(5, 9) = 0.0;
    Corr(5, 10) = 0.;
    Corr(5, 11) = 0.01;
    Corr(6, 6) = 1.;
    Corr(6, 7) = 0.04;
    Corr(6, 8) = 0.0;
    Corr(6, 9) = 0.0;
    Corr(6, 10) = 0.0;
    Corr(6, 11) = 0.;
    Corr(7, 7) = 1.;
    Corr(7, 8) = 0.0;
    Corr(7, 9) = 0.05;
    Corr(7, 10) = 0.;
    Corr(7, 11) = 0.;
    Corr(8, 8) = 1.;
    Corr(8, 9) = 0.04;
    Corr(8, 10) = 0.0;
    Corr(8, 11) = 0.0;
    Corr(9, 9) = 1.;
    Corr(9, 10) = 0.0;
    Corr(9, 11) = 0.;
    Corr(10, 10) = 1.;
    Corr(10, 11) = 0.0;
    Corr(11, 11) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr);

    // Correlation Matrix (syst)
    Corr2.ResizeTo(12, 12);
    Corr2(0, 0) = 1.;
    Corr2(0, 1) = 0.;
    Corr2(0, 2) = -0.05;
    Corr2(0, 3) = 0.81;
    Corr2(0, 4) = 0.15;
    Corr2(0, 5) = -0.19;
    Corr2(0, 6) = 0.78;
    Corr2(0, 7) = -0.03;
    Corr2(0, 8) = 0.79;
    Corr2(0, 9) = -0.03;
    Corr2(0, 10) = 0.78;
    Corr2(0, 11) = -0.03;
    Corr2(1, 1) = 1.;
    Corr2(1, 2) = 0.15;
    Corr2(1, 3) = 0.14;
    Corr2(1, 4) = 0.07;
    Corr2(1, 5) = 0.06;
    Corr2(1, 6) = 0.12;
    Corr2(1, 7) = 0.05;
    Corr2(1, 8) = 0.11;
    Corr2(1, 9) = 0.03;
    Corr2(1, 10) = 0.05;
    Corr2(1, 11) = -0.01;
    Corr2(2, 2) = 1.;
    Corr2(2, 3) = -0.19;
    Corr2(2, 4) = 0.07;
    Corr2(2, 5) = 0.04;
    Corr2(2, 6) = -0.26;
    Corr2(2, 7) = 0.14;
    Corr2(2, 8) = -0.25;
    Corr2(2, 9) = 0.04;
    Corr2(2, 10) = -0.26;
    Corr2(2, 11) = -0.01;
    Corr2(3, 3) = 1.;
    Corr2(3, 4) = -0.01;
    Corr2(3, 5) = -0.09;
    Corr2(3, 6) = 0.77;
    Corr2(3, 7) = 0.0;
    Corr2(3, 8) = 0.77;
    Corr2(3, 9) = -0.01;
    Corr2(3, 10) = 0.76;
    Corr2(3, 11) = -0.30;
    Corr2(4, 4) = 1.;
    Corr2(4, 5) = 0.15;
    Corr2(4, 6) = 0.18;
    Corr2(4, 7) = 0.;
    Corr2(4, 8) = 0.17;
    Corr2(4, 9) = 0.01;
    Corr2(4, 10) = 0.19;
    Corr2(4, 11) = 0.04;
    Corr2(5, 5) = 1.;
    Corr2(5, 6) = -0.24;
    Corr2(5, 7) = 0.;
    Corr2(5, 8) = -0.24;
    Corr2(5, 9) = 0.01;
    Corr2(5, 10) = -0.23;
    Corr2(5, 11) = 0.04;
    Corr2(6, 6) = 1.;
    Corr2(6, 7) = -0.03;
    Corr2(6, 8) = 0.94;
    Corr2(6, 9) = -0.02;
    Corr2(6, 10) = 0.93;
    Corr2(6, 11) = 0.;
    Corr2(7, 7) = 1.;
    Corr2(7, 8) = -0.01;
    Corr2(7, 9) = 0.15;
    Corr2(7, 10) = -0.01;
    Corr2(7, 11) = -0.02;
    Corr2(8, 8) = 1.;
    Corr2(8, 9) = -0.03;
    Corr2(8, 10) = 0.95;
    Corr2(8, 11) = 0.;
    Corr2(9, 9) = 1.;
    Corr2(9, 10) = -0.02;
    Corr2(9, 11) = 0.;
    Corr2(10, 10) = 1.;
    Corr2(10, 11) = -0.03;
    Corr2(11, 11) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr2);

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("2401.17934Bd", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));
  }

  if (comb == 1)
  {

    // GGSZ: xpm_Dkstz, ymp_Dkstz
    // https://arxiv.org/abs/2309.05514
    // Observables 4:
    CorrData.clear();
    CorrData.push_back(dato(0.074, 0.086, 0.005, 0.011));  // xp_dkstz
    CorrData.push_back(dato(-0.215, 0.086, 0.004, 0.013)); // xm_dkstz
    CorrData.push_back(dato(-0.336, 0.105, 0.017, 0.009)); // yp_dkstz
    CorrData.push_back(dato(-0.012, 0.128, 0.024, 0.011)); // ym_dkstz

    // Correlation Matrix (stat):
    Corr.ResizeTo(4, 4);
    Corr(0, 0) = 1.;
    Corr(0, 1) = 0.;
    Corr(0, 2) = 0.18;
    Corr(0, 3) = 0.;
    Corr(1, 1) = 1.;
    Corr(1, 2) = 0.;
    Corr(1, 3) = 0.08;
    Corr(2, 2) = 1.;
    Corr(2, 3) = 0.;
    Corr(3, 3) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr);

    // Correlation Matrix (syst)
    Corr2.ResizeTo(4, 4);
    Corr2(0, 0) = 1.;
    Corr2(0, 1) = -0.02;
    Corr2(0, 2) = 0.05;
    Corr2(0, 3) = 0.;
    Corr2(1, 1) = 1.;
    Corr2(1, 2) = 0.05;
    Corr2(1, 3) = 0.04;
    Corr2(2, 2) = 1.;
    Corr2(2, 3) = -0.07;
    Corr2(3, 3) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr2);

    // Correlation Matrix (ext)
    Corr3.ResizeTo(4, 4);
    Corr3(0, 0) = 1.;
    Corr3(0, 1) = -0.14;
    Corr3(0, 2) = 0.34;
    Corr3(0, 3) = -0.09;
    Corr3(1, 1) = 1.;
    Corr3(1, 2) = -0.04;
    Corr3(1, 3) = 0.17;
    Corr3(2, 2) = 1.;
    Corr3(2, 3) = -0.04;
    Corr3(3, 3) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr3);

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("2309.05514", CorrelatedGaussianObservables(CorrData, Corr, Corr2, Corr3)));
  }
  else if (comb == 3)
  {

    // GGSZ LHCb ChargedB + NeutralBd
    // Observables 26

    //-------------------------------------------------  Bpm -> Dhpm -------------------------------------------------------------------------

    // D -> K0spipi, D -> K0sKK
    // https://arxiv.org/pdf/2010.08483
    // Observables 6:
    CorrData.clear();
    CorrData.push_back(dato(0.0568, 0.0096, 0.0020, 0.0023));     // xm_dk
    CorrData.push_back(dato(0.06550, 0.01140, 0.0025, 0.0035)); // ym_dk
    CorrData.push_back(dato(-0.09300, 0.0098, 0.0024, 0.0018));    // xp_dk
    CorrData.push_back(dato(-0.0125, 0.0123, 0.0026, 0.0028));    // yp_dk
    CorrData.push_back(dato(-0.05470, 0.0199, 0.0032, 0.0014));   // xi_x_dpi
    CorrData.push_back(dato(0.00710, 0.02330, 0.0054, 0.0018));   // xi_y_dpi

    //----------------------------------------------------------------------------------------------------------------------------------

    //-------------------------------------------------  Bpm -> D*hpm -------------------------------------------------------------------------


    // GGSZ: D -> K0spipi, D -> K0sKK
    // https://arxiv.org/abs/2310.04277 LHCb
    // Observables 6:
    CorrData.push_back(dato(11.42e-2, 3.16e-2, 1.26e-2, 0.41e-2));  // xp_dstk
    CorrData.push_back(dato(-8.91e-2, 3.55e-2, 2.04e-2, 0.23e-2));  // xm_dstk
    CorrData.push_back(dato(3.60e-2, 4.41e-2, 2.12e-2, 0.30e-2));   // yp_dstk
    CorrData.push_back(dato(-16.75e-2, 3.98e-2, 1.48e-2, 0.64e-2)); // ym_dstk
    CorrData.push_back(dato(0.51e-2, 5.00e-2, 2.66e-2, 0.93e-2));   // Re xi_dstpi
    CorrData.push_back(dato(7.92e-2, 5.04e-2, 3.78e-2, 0.83e-2));   // Im xi_dstpi

    // GGSZ: D -> K0spipi, D -> K0sKK
    // https://arxiv.org/abs/2311.10434 LHCb
    // Observables 6:
    CorrData.push_back(dato(-6.3e-2, 2.9e-2, 1.1e-2, 0.6e-2)); // xm_dstk
    CorrData.push_back(dato(-4.8e-2, 5.7e-2, 1.4e-2, 1.5e-2)); // ym_dstk
    CorrData.push_back(dato(6.0e-2, 2.6e-2, 0.9e-2, 0.2e-2));  // xp_dstk
    CorrData.push_back(dato(5.4e-2, 2.9e-2, 0.9e-2, 0.4e-2));  // yp_dstk
    CorrData.push_back(dato(11.5e-2, 9.4e-2, 3.3e-2, 2.3e-2)); // Re xi_dstpi
    CorrData.push_back(dato(-0.9e-2, 9.7e-2, 2.5e-2, 2.1e-2)); // Im xi_dstpi

    //----------------------------------------------------------------------------------------------------------------------------------

    //-------------------------------------------------  Bpm -> DK*pm -------------------------------------------------------------------------


    // GGSZ: D -> K0spipi, D -> K0sKK
    // LHCB-PAPER-2024-023
    // Observables 4:
    CorrData.push_back(dato(0.135, 0.056, 0.019, 0.009)); // xm_dkst
    CorrData.push_back(dato(-0.170, 0.108, 0.013, 0.039)); // ym_dkst
    CorrData.push_back(dato(0.003, 0.052, 0.018, 0.004)); // xp_dkst
    CorrData.push_back(dato(0.054, 0.061, 0.009, 0.019)); // yp_dkst


    // GGSZ: xpm_Dkstz, ymp_Dkstz
    // https://arxiv.org/abs/2309.05514
    // Observables 4:
    CorrData.push_back(dato(0.074, 0.086, 0.005, 0.011));       // xp_dkstz
    CorrData.push_back(dato(-0.215, 0.086, 0.004, 0.013));      // xm_dkstz
    CorrData.push_back(dato(-0.336, 0.105, 0.017, 0.009));      // yp_dkstz
    CorrData.push_back(dato(-0.012, 0.128, 0.024, 0.011));      // ym_dkstz


    // Correlation Matrix (stat):
    Corr.ResizeTo(26, 26);
    Corr = 0.;

    Corr(0, 0) = 1.;
    Corr(0, 1) = -0.125;
    Corr(0, 2) = -0.013;
    Corr(0, 3) = 0.019;
    Corr(0, 4) = 0.037;
    Corr(0, 5) = -0.161;
    Corr(1, 1) = 1.;
    Corr(1, 2) = -0.011;
    Corr(1, 3) = -0.010;
    Corr(1, 4) = 0.097;
    Corr(1, 5) = 0.041;
    Corr(2, 2) = 1.;
    Corr(2, 3) = 0.105;
    Corr(2, 4) = -0.108;
    Corr(2, 5) = 0.032;
    Corr(3, 3) = 1.;
    Corr(3, 4) = -0.07;
    Corr(3, 5) = -0.147;
    Corr(4, 4) = 1.;
    Corr(4, 5) = 0.15;
    Corr(5, 5) = 1.;

    Corr(6, 6) = 1.;
    Corr(6, 7) = -0.05;
    Corr(6, 8) = -0.08;
    Corr(6, 9) = -0.08;
    Corr(6, 10) = -0.16;
    Corr(6, 11) = -0.14;
    Corr(7, 7) = 1.;
    Corr(7, 8) = 0.03;
    Corr(7, 9) = -0.08;
    Corr(7, 10) = 0.25;
    Corr(7, 11) = 0.06;
    Corr(8, 8) = 1.;
    Corr(8, 9) = -0.09;
    Corr(8, 10) = 0.;
    Corr(8, 11) = -0.19;
    Corr(9, 9) = 1.;
    Corr(9, 10) = 0.;
    Corr(9, 11) = 0.33;
    Corr(10, 10) = 1.;
    Corr(10, 11) = 0.18;
    Corr(11, 11) = 1.;

    Corr(12, 12) = 1.;
    Corr(12, 13) = 0.420;
    Corr(12, 14) = 0.158;
    Corr(12, 15) = 0.105;
    Corr(12, 16) = 0.445;
    Corr(12, 17) = 0.422;
    Corr(13, 13) = 1.;
    Corr(13, 14) = 0.115;
    Corr(13, 15) = 0.232;
    Corr(13, 16) = 0.765;
    Corr(13, 17) = 0.631;
    Corr(14, 14) = 1.;
    Corr(14, 15) = -0.095;
    Corr(14, 16) = 0.012;
    Corr(14, 17) = 0.409;
    Corr(15, 15) = 1.;
    Corr(15, 16) = 0.263;
    Corr(15, 17) = 0.112;
    Corr(16, 16) = 1.;
    Corr(16, 17) = 0.597;
    Corr(17, 17) = 1.;

    Corr(18,18) = 1.;
    Corr(18,19) = 0.448;
    Corr(18,20) = 0.;
    Corr(18,21) = 0.;
    Corr(19,19) = 1.;
    Corr(19,20) = 0.;
    Corr(19,21) = 0.;
    Corr(20,20) = 1.;
    Corr(20,21) = -0.144;
    Corr(21,21) = 1.;

    Corr(22, 22) = 1.;
    Corr(22, 23) = 0.;
    Corr(22, 24) = 0.18;
    Corr(22, 25) = 0.;
    Corr(23, 23) = 1.;
    Corr(23, 24) = 0.;
    Corr(23, 25) = 0.08;
    Corr(24, 24) = 1.;
    Corr(24, 25) = 0.;
    Corr(25, 25) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr);


    // Correlation Matrix (syst)
    Corr2.ResizeTo(26, 26);
    Corr2 = 0.;


    Corr2(0, 0) = 1.;
    Corr2(0, 1) = 0.864;
    Corr2(0, 2) = 0.734;
    Corr2(0, 3) = 0.897;
    Corr2(0, 4) = 0.349;
    Corr2(0, 5) = 0.318;
    Corr2(1, 1) = 1.;
    Corr2(1, 2) = 0.874;
    Corr2(1, 3) = 0.903;
    Corr2(1, 4) = 0.408;
    Corr2(1, 5) = 0.362;
    Corr2(2, 2) = 1.;
    Corr2(2, 3) = 0.771;
    Corr2(2, 4) = 0.563;
    Corr2(2, 5) = 0.447;
    Corr2(3, 3) = 1.;
    Corr2(3, 4) = 0.507;
    Corr2(3, 5) = 0.451;
    Corr2(4, 4) = 1.;
    Corr2(4, 5) = 0.484;
    Corr2(5, 5) = 1.;

    Corr2(6, 6) = 1.;
    Corr2(6, 7) = -0.01;
    Corr2(6, 8) = 0.;
    Corr2(6, 9) = -0.01;
    Corr2(6, 10) = -0.01;
    Corr2(6, 11) = 0.06;
    Corr2(7, 7) = 1.;
    Corr2(7, 8) = 0.02;
    Corr2(7, 9) = 0.01;
    Corr2(7, 10) = -0.02;
    Corr2(7, 11) = -0.05;
    Corr2(8, 8) = 1.;
    Corr2(8, 9) = -0.01;
    Corr2(8, 10) = -0.09;
    Corr2(8, 11) = -0.02;
    Corr2(9, 9) = 1.;
    Corr2(9, 10) = 0.01;
    Corr2(9, 11) = -0.1;
    Corr2(10, 10) = 1.;
    Corr2(10, 11) = 0.02;
    Corr2(11, 11) = 1.;

    Corr2(12, 12) = 1.;
    Corr2(12, 13) = 0.630;
    Corr2(12, 14) = -0.241;
    Corr2(12, 15) = -0.016;
    Corr2(12, 16) = 0.602;
    Corr2(12, 17) = 0.083;
    Corr2(13, 13) = 1.;
    Corr2(13, 14) = 0.008;
    Corr2(13, 15) = 0.154;
    Corr2(13, 16) = 0.735;
    Corr2(13, 17) = 0.230;
    Corr2(14, 14) = 1.;
    Corr2(14, 15) = 0.515;
    Corr2(14, 16) = 0.232;
    Corr2(14, 17) = 0.618;
    Corr2(15, 15) = 1.;
    Corr2(15, 16) = 0.237;
    Corr2(15, 17) = 0.112;
    Corr2(16, 16) = 1.;
    Corr2(16, 17) = -0.201;
    Corr2(17, 17) = 1.;

    Corr2(18,18) = 1.;
    Corr2(18,19) = 0.319;
    Corr2(18,20) = 0.030;
    Corr2(18,21) = -0.158;
    Corr2(19,19) = 1.;
    Corr2(19,20) = 0.022;
    Corr2(19,21) = -0.478;
    Corr2(20,20) = 1.;
    Corr2(20,21) = -0.055;
    Corr2(21,21) = 1.;

    Corr2(22, 22) = 1.;
    Corr2(22, 23) = -0.02;
    Corr2(22, 24) = 0.05;
    Corr2(22, 25) = 0.;
    Corr2(23, 23) = 1.;
    Corr2(23, 24) = 0.05;
    Corr2(23, 25) = 0.04;
    Corr2(24, 24) = 1.;
    Corr2(24, 25) = -0.07;
    Corr2(25, 25) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr2);

    // Correlation Matrix (syst)
    Corr3.ResizeTo(26, 26);
    Corr3 = 0.;


    Corr3(0, 0) = 1.;
    Corr3(0, 1) = -0.047;
    Corr3(0, 2) = -0.490;
    Corr3(0, 3) = 0.322;
    Corr3(0, 4) = 0.189;
    Corr3(0, 5) = 0.144;
    Corr3(0, 6) = 0.27;
    Corr3(0, 7) = 0.07;
    Corr3(0, 8) = 0.;
    Corr3(0, 9) = 0.21;
    Corr3(0, 10) = -0.23;
    Corr3(0, 11) = -0.05;
    Corr3(0, 12) = 0.08;
    Corr3(0, 13) = 0.09;
    Corr3(0, 14) = 0.29;
    Corr3(0, 15) = 0.38;
    Corr3(0, 16) = 0.21;
    Corr3(0, 17) = 0.06;
    Corr3(0, 18) = -0.27;
    Corr3(0, 19) = -0.07;
    Corr3(0, 20) = 0.20;
    Corr3(0, 21) = -0.04;
    Corr3(0, 22) = -0.29;
    Corr3(0, 23) = -0.31;
    Corr3(0, 24) = 0.35;
    Corr3(0, 25) = 0.13;
    Corr3(1, 1) = 1.;
    Corr3(1, 2) = 0.059;
    Corr3(1, 3) = -0.237;
    Corr3(1, 4) = -0.13;
    Corr3(1, 5) = -0.117;
    Corr3(1, 6) = 0.09;
    Corr3(1, 7) = 0.04;
    Corr3(1, 8) = 0.17;
    Corr3(1, 9) = - 0.05;
    Corr3(1, 10) = 0.07;
    Corr3(1, 11) = 0.11;
    Corr3(1, 12) = 0.12;
    Corr3(1, 13) = 0.07;
    Corr3(1, 14) = 0.22;
    Corr3(1, 15) = 0.14;
    Corr3(1, 16) = 0.24;
    Corr3(1, 17) = 0.35;
    Corr3(1, 18) = -0.53;
    Corr3(1, 19) = -0.65;
    Corr3(1, 20) = -0.24;
    Corr3(1, 21) = 0.46;
    Corr3(1, 22) = -0.11;
    Corr3(1, 23) = 0.22;
    Corr3(1, 24) = 0.12;
    Corr3(1, 25) = 0.22;
    Corr3(2, 2) = 1.;
    Corr3(2, 3) = 0.061;
    Corr3(2, 4) = 0.004;
    Corr3(2, 5) = -0.139;
    Corr3(2, 6) = -0.14;
    Corr3(2, 7) = 0.04;
    Corr3(2, 8) = -0.16;
    Corr3(2, 9) = -0.07;
    Corr3(2, 10) = 0.09;
    Corr3(2, 11) = 0.05;
    Corr3(2, 12) = -0.16;
    Corr3(2, 13) = -0.17;
    Corr3(2, 14) = -0.01;
    Corr3(2, 15) = -0.06;
    Corr3(2, 16) = -0.20;
    Corr3(2, 17) = -0.07;
    Corr3(2, 18) = -0.14;
    Corr3(2, 19) = 0.;
    Corr3(2, 20) = -0.16;
    Corr3(2, 21) = 0.01;
    Corr3(2, 22) = -0.06;
    Corr3(2, 23) = 0.48;
    Corr3(2, 24) = 0.03;
    Corr3(2, 25) = -0.15;
    Corr3(3, 3) = 1.;
    Corr3(3, 4) = 0.14;
    Corr3(3, 5) = -0.199;
    Corr3(3, 6) = 0.1;
    Corr3(3, 7) = -0.03;
    Corr3(3, 8) = -0.2;
    Corr3(3, 9) = 0.09;
    Corr3(3, 10) = -0.09;
    Corr3(3, 11) = 0.;
    Corr3(3, 12) = -0.05;
    Corr3(3, 13) = -0.05;
    Corr3(3, 14) = 0.41;
    Corr3(3, 15) = -0.17;
    Corr3(3, 16) = -0.08;
    Corr3(3, 17) = 0.;
    Corr3(3, 18) = -0.42;
    Corr3(3, 19) = -0.25;
    Corr3(3, 20) = 0.27;
    Corr3(3, 21) = -0.08;
    Corr3(3, 22) = -0.06;
    Corr3(3, 23) = -0.49;
    Corr3(3, 24) = 0.27;
    Corr3(3, 25) = -0.01;
    Corr3(4, 4) = 1.;
    Corr3(4, 5) = 0.638;
    Corr3(4, 6) = 0.11;
    Corr3(4, 7) = 0.09;
    Corr3(4, 8) = 0.02;
    Corr3(4, 9) = 0.05;
    Corr3(4, 10) = -0.1;
    Corr3(4, 11) = 0.06;
    Corr3(4, 12) = -0.47;
    Corr3(4, 13) = -0.50;
    Corr3(4, 14) = -0.23;
    Corr3(4, 15) = 0.54;
    Corr3(4, 16) = -0.41;
    Corr3(4, 17) = -0.50;
    Corr3(4, 18) = 0.01;
    Corr3(4, 19) = 0.27;
    Corr3(4, 20) = -0.08;
    Corr3(4, 21) = -0.28;
    Corr3(4, 22) = -0.16;
    Corr3(4, 23) = 0.04;
    Corr3(4, 24) = 0.29;
    Corr3(4, 25) = 0.21;
    Corr3(5, 5) = 1.;
    Corr3(5, 6) = 0.08;
    Corr3(5, 7) = 0.04;
    Corr3(5, 8) = 0.16;
    Corr3(5, 9) = 0.03;
    Corr3(5, 10) = -0.05;
    Corr3(5, 11) = 0.02;
    Corr3(5, 12) = -0.37;
    Corr3(5, 13) = -0.54;
    Corr3(5, 14) = -0.13;
    Corr3(5, 15) = 0.56;
    Corr3(5, 16) = -0.42;
    Corr3(5, 17) = -0.54;
    Corr3(5, 18) = 0.30;
    Corr3(5, 19) = 0.53;
    Corr3(5, 20) = 0.32;
    Corr3(5, 21) = -0.53;
    Corr3(5, 22) = 0.30;
    Corr3(5, 23) = 0.12;
    Corr3(5, 24) = 0.32;
    Corr3(5, 25) = 0.36;

    Corr3(6, 6) = 1.;
    Corr3(6, 7) = 0.26;
    Corr3(6, 8) = 0.44;
    Corr3(6, 9) = 0.07;
    Corr3(6, 10) = -0.59;
    Corr3(6, 11) = 0.24;
    Corr3(6, 12) = 0.03;
    Corr3(6, 13) = 0.01;
    Corr3(6, 14) = 0.11;
    Corr3(6, 15) = 0.16;
    Corr3(6, 16) = 0.07;
    Corr3(6, 17) = 0.04;
    Corr3(6, 18) = -0.03;
    Corr3(6, 19) = -0.05;
    Corr3(6, 20) = 0.01;
    Corr3(6, 21) = 0.04;
    Corr3(6, 22) = -0.05;
    Corr3(6, 23) = -0.14;
    Corr3(6, 24) = 0.15;
    Corr3(6, 25) = 0.04;
    Corr3(7, 7) = 1.;
    Corr3(7, 8) = -0.21;
    Corr3(7, 9) = -0.58;
    Corr3(7, 10) = -0.53;
    Corr3(7, 11) = -0.39;
    Corr3(7, 12) = -0.02;
    Corr3(7, 13) = -0.03;
    Corr3(7, 14) = -0.04;
    Corr3(7, 15) = 0.1;
    Corr3(7, 16) = -0.02;
    Corr3(7, 17) = -0.04;
    Corr3(7, 18) = -0.05;
    Corr3(7, 19) = -0.05;
    Corr3(7, 20) = 0.01;
    Corr3(7, 21) = 0.02;
    Corr3(7, 22) = 0.02;
    Corr3(7, 23) = 0.05;
    Corr3(7, 24) = 0.08;
    Corr3(7, 25) = 0.02;
    Corr3(8, 8) = 1.;
    Corr3(8, 9) = -0.04;
    Corr3(8, 10) = 0.08;
    Corr3(8, 11) = 0.44;
    Corr3(8, 12) = 0.19;
    Corr3(8, 13) = 0.26;
    Corr3(8, 14) = -0.20;
    Corr3(8, 15) = 0.14;
    Corr3(8, 16) = 0.26;
    Corr3(8, 17) = 0.17;
    Corr3(8, 18) = -0.03;
    Corr3(8, 19) = -0.05;
    Corr3(8, 20) = 0.03;
    Corr3(8, 21) = 0.02;
    Corr3(8, 22) = 0.15;
    Corr3(8, 23) = 0.06;
    Corr3(8, 24) = -0.04;
    Corr3(8, 25) = 0.19;
    Corr3(9, 9) = 1.;
    Corr3(9, 10) = -0.24;
    Corr3(9, 11) = -0.03;
    Corr3(9, 12) = 0.04;
    Corr3(9, 13) = 0.05;
    Corr3(9, 14) = 0.05;
    Corr3(9, 15) = 0.1;
    Corr3(9, 16) = 0.06;
    Corr3(9, 17) = 0.03;
    Corr3(9, 18) = 0.02;
    Corr3(9, 19) = 0.02;
    Corr3(9, 20) = -0.02;
    Corr3(9, 21) = -0.01;
    Corr3(9, 22) = -0.06;
    Corr3(9, 23) = -0.04;
    Corr3(9, 24) = 0.11;
    Corr3(9, 25) = 0.01;
    Corr3(10, 10) = 1.;
    Corr3(10, 11) = 0.57;
    Corr3(10, 12) = 0.04;
    Corr3(10, 13) = 0.02;
    Corr3(10, 14) = 0.01;
    Corr3(10, 15) = -0.16;
    Corr3(10, 16) = 0.;
    Corr3(10, 17) = 0.05;
    Corr3(10, 18) = 0.06;
    Corr3(10, 19) = 0.06;
    Corr3(10, 20) = 0.;
    Corr3(10, 21) = -0.02;
    Corr3(10, 22) = 0.04;
    Corr3(10, 23) = 0.09;
    Corr3(10, 24) = -0.11;
    Corr3(10, 25) = 0.05;
    Corr3(11, 11) = 1.;
    Corr3(11, 12) = 0.01;
    Corr3(11, 13) = 0.01;
    Corr3(11, 14) = 0.07;
    Corr3(11, 15) = -0.04;
    Corr3(11, 16) = 0.02;
    Corr3(11, 17) = 0.06;
    Corr3(11, 18) = 0.05;
    Corr3(11, 19) = 0.03;
    Corr3(11, 20) = -0.01;
    Corr3(11, 21) = 0.01;
    Corr3(11, 22) = -0.04;
    Corr3(11, 23) = 0.01;
    Corr3(11, 24) = 0.03;
    Corr3(11, 25) = 0.10;

    Corr3(12, 12) = 1.;
    Corr3(12, 13) = 0.886;
    Corr3(12, 14) = -0.117;
    Corr3(12, 15) = -0.063;
    Corr3(12, 16) = 0.916;
    Corr3(12, 17) = 0.904;
    Corr3(12, 18) = -0.04;
    Corr3(12, 19) = -0.01;
    Corr3(12, 20) = 0.03;
    Corr3(12, 21) = 0.04;
    Corr3(12, 22) = 0.07;
    Corr3(12, 23) = -0.08;
    Corr3(12, 24) = -0.14;
    Corr3(12, 25) = -0.08;
    Corr3(13, 13) = 1.;
    Corr3(13, 14) = -0.124;
    Corr3(13, 15) = -0.115;
    Corr3(13, 16) = 0.935;
    Corr3(13, 17) = 0.889;
    Corr3(13, 18) = -0.03;
    Corr3(13, 19) = 0.;
    Corr3(13, 20) = 0.02;
    Corr3(13, 21) = 0.01;
    Corr3(13, 22) = -0.04;
    Corr3(13, 23) = -0.16;
    Corr3(13, 24) = -0.24;
    Corr3(13, 25) = -0.08;
    Corr3(14, 14) = 1.;
    Corr3(14, 15) = -0.280;
    Corr3(14, 16) = -0.130;
    Corr3(14, 17) = 0.063;
    Corr3(14, 18) = 0.01;
    Corr3(14, 19) = 0.04;
    Corr3(14, 20) = -0.04;
    Corr3(14, 21) = 0.;
    Corr3(14, 22) = -0.16;
    Corr3(14, 23) = -0.18;
    Corr3(14, 24) = 0.05;
    Corr3(14, 25) = 0.25;
    Corr3(15, 15) = 1.;
    Corr3(15, 16) = 0.085;
    Corr3(15, 17) = -0.077;
    Corr3(15, 18) = -0.02;
    Corr3(15, 19) = -0.05;
    Corr3(15, 20) = 0.01;
    Corr3(15, 21) = 0.01;
    Corr3(15, 22) = -0.07;
    Corr3(15, 23) = 0.14;
    Corr3(15, 24) = 0.33;
    Corr3(15, 25) = 0.14;
    Corr3(16, 16) = 1.;
    Corr3(16, 17) = 0.941;
    Corr3(16, 18) = -0.04;
    Corr3(16, 19) = -0.02;
    Corr3(16, 20) = 0.01;
    Corr3(16, 21) = 0.02;
    Corr3(16, 22) = -0.08;
    Corr3(16, 23) = -0.1;
    Corr3(16, 24) = -0.13;
    Corr3(16, 25) = 0.03;
    Corr3(17, 17) = 1.;
    Corr3(17, 18) = -0.04;
    Corr3(17, 19) = -0.01;
    Corr3(17, 20) = 0.;
    Corr3(17, 21) = 0.03;
    Corr3(17, 22) = -0.08;
    Corr3(17, 23) = -0.1;
    Corr3(17, 24) = -0.18;
    Corr3(17, 25) = 0.03;

    Corr3(18,18) = 1.;
    Corr3(18,19) = 0.71;
    Corr3(18,20) = -0.01;
    Corr3(18,21) = -0.02;
    Corr3(18, 22) = 0.05;
    Corr3(18, 23) = -0.03;
    Corr3(18, 24) = 0.03;
    Corr3(18, 25) = -0.03;
    Corr3(19,19) = 1.;
    Corr3(19,20) = -0.05;
    Corr3(19,21) = -0.29;
    Corr3(19, 22) = 0.06;
    Corr3(19, 23) = -0.07;
    Corr3(19, 24) = 0.01;
    Corr3(19, 25) = -0.04;
    Corr3(20,20) = 1.;
    Corr3(20,21) = -0.15;
    Corr3(20, 22) = 0.09;
    Corr3(20, 23) = -0.02;
    Corr3(20, 24) = 0.01;
    Corr3(20, 25) = -0.02;
    Corr3(21,21) = 1.;
    Corr3(21, 22) = -0.03;
    Corr3(21, 23) = 0.08;
    Corr3(21, 24) = -0.02;
    Corr3(21, 25) = 0.01;

    Corr3(22, 22) = 1.;
    Corr3(22, 23) = -0.14;
    Corr3(22, 24) = 0.34;
    Corr3(22, 25) = -0.09;
    Corr3(23, 23) = 1.;
    Corr3(23, 24) = -0.04;
    Corr3(23, 25) = 0.17;
    Corr3(24, 24) = 1.;
    Corr3(24, 25) = -0.04;
    Corr3(25, 25) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr3);

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("DKst0Pcomb", CorrelatedGaussianObservables(CorrData, Corr, Corr2, Corr3)));
  }

  //----------------------------------------------------------------------------------------------------------------------------------

  // B0 -> DK^*0
  // https://arxiv.org/pdf/1509.01098 (Belle)
  // Observables 4:
  CorrData.clear();
  CorrData.push_back(dato(0.2, 0.55, 0.05, 0.1));  // xp_dkstz
  CorrData.push_back(dato(0.2, 0.65, 0.05, 0.1)); // yp_dkstz
  CorrData.push_back(dato(0.55, 0.8, 0.05, 0.)); // xm_dkstz
  CorrData.push_back(dato(-0.65, 0.9, 0.05, 0.1)); // ym_dkstz

  // Correlation Matrix (stat):
  Corr.ResizeTo(4, 4);
  Corr(0, 0) = 1.;
  Corr(0, 1) = -0.075;
  Corr(0, 2) = 0.;
  Corr(0, 3) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.;
  Corr(1, 3) = 0.;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.07;
  Corr(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(4, 4);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.;
  Corr2(0, 2) = 0.;
  Corr2(0, 3) = 0.;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.;
  Corr2(1, 3) = 0.;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.;
  Corr2(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  // Correlation Matrix (ext)
  Corr3.ResizeTo(4, 4);
  Corr3(0, 0) = 1.;
  Corr3(0, 1) = 0.;
  Corr3(0, 2) = 0.;
  Corr3(0, 3) = 0.;
  Corr3(1, 1) = 1.;
  Corr3(1, 2) = 0.;
  Corr3(1, 3) = 0.;
  Corr3(2, 2) = 1.;
  Corr3(2, 3) = 0.;
  Corr3(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr3);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("1509.01098", CorrelatedGaussianObservables(CorrData, Corr, Corr2, Corr3)));



  //-------------------------------------- B0d -> Dpi -------------------------------------------------------------------------

  // Time dependent rates, D -> Kpipi
  // https://arxiv.org/pdf/1805.03448
  // Observables 2:
  CorrData.clear();
  CorrData.push_back(dato(0.058, 0.02, 0.011)); // s_dmpi
  CorrData.push_back(dato(0.038, 0.02, 0.007)); // sb_dmpi

  // Correlation Matrix (stat):
  Corr.ResizeTo(2, 2);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.6;
  Corr(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(2, 2);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = -0.41;
  Corr2(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID12", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //----------------------------------------------------------------------------------------------------------------------------------


  //-------------------------------------- Babar Measurements -------------------------------------------------------------------------

   // https://arxiv.org/pdf/hep-ex/0602049
   // B^0 -> D-pi+, B0 -> D^*(Full rec.)pi, B0 -> Drho
   meas.insert(pair<string, dato>("0602049_aDpi", dato(-0.01, 0.023, 0.007))); // a_Dpi
   meas.insert(pair<string, dato>("0602049_cDpi", dato(-0.033, 0.042, 0.012))); // c_Dpi (lep)
   meas.insert(pair<string, dato>("0602049_aDstarpi", dato(-0.04, 0.023, 0.01))); // a_Dstarpi
   meas.insert(pair<string, dato>("0602049_cDstarpi", dato(0.049, 0.042, 0.015))); // c_Dstarpi (lep)
   meas.insert(pair<string, dato>("0602049_aDrho", dato(-0.024, 0.031, 0.009))); // a_Drho
   meas.insert(pair<string, dato>("0602049_cDrho", dato(-0.098, 0.055, 0.018))); // c_Drho (lep)

   // https://arxiv.org/pdf/hep-ex/0504035
   // B^0 -> Dstar-pi+ (Partial rec.)
   meas.insert(pair<string, dato>("0504035_aDstarpi", dato(-0.034, 0.014, 0.009))); // a_Dstarpi
   meas.insert(pair<string, dato>("0504035_cDstarpi", dato(-0.019, 0.022, 0.013))); // c_Dstarpi (lep)

  //---------------------------------------------------------------------------------------------------------------



  //-------------------------------------- Belle Measurements -------------------------------------------------------------------------

  // https://arxiv.org/pdf/hep-ex/0604013 (HFLAV conversion)
  // B^0 -> D-pi+, B0 -> D^*(Full rec.)pi
  meas.insert(pair<string, dato>("0604013_aDpi", dato(-0.050, 0.021, 0.012))); // a_Dpi
  meas.insert(pair<string, dato>("0604013_cDpi", dato(0.019, 0.021, 0.012))); // c_Dpi (lep)
  meas.insert(pair<string, dato>("0604013_aDstarpi", dato(-0.039, 0.020, 0.013))); // a_Dstarpi
  meas.insert(pair<string, dato>("0604013_cDstarpi", dato(-0.011, 0.02, 0.013))); // c_Dstarpi (lep)

  // https://arxiv.org/pdf/1102.0888
  // B^0 -> Dstar-pi+ (Partial rec.)
  meas.insert(pair<string, dato>("11020888_aDstarpi", dato(-0.046, 0.011, 0.015))); // a_Dstarpi
  meas.insert(pair<string, dato>("11020888_cDstarpi", dato(-0.015, 0.011, 0.015))); // c_Dstarpi (lep)


  //---------------------------------------------------------------------------------------------------------------




  //-------------------------------------  Other useful inputs for the Time Dependent B0 Decay and mixing parameters -------------------------------------------------------------------------

  // kB_Dkstz
  // https://arxiv.org/pdf/1602.03455
  // Observables 1:
  meas.insert(pair<string, dato>("UID25", dato(0.934, 0.0075, 0.024))); // k_dkstz

  // sin(2beta)
  // https://indico.cern.ch/event/1291157/contributions/5903548/attachments/2900988/5087304/bona-utfit.pdf
  // Observables 1:
  meas.insert(pair<string, dato>("UID27", dato(0.699, 0.015, 0.))); // sin(phi_d=2beta)

  //----------------------------------------------------------------------------------------------------------------------------------

}
// ---------------------------------------------------------

// ---------------------------------------------------------
void MixingModel::Add_NeutralBs_meas()
{

  vector<dato> CorrData;
  TMatrixDSym Corr, Corr2;

  //-------------------------------------- Bs^0 -> DKst0bar -------------------------------------------------------------------------

  if (comb == 2)
  { // If treating it separately from the Bd counterpart

    // GLW: D -> pipi, D -> KK, D -> 4pi; ADS: D -> Kpi, D -> K3pi
    // https://arxiv.org/pdf/2401.17934
    // Observables 12:
    CorrData.clear();
    CorrData.push_back(dato(-0.009, 0.011, 0.020)); // afav_dkstz_kpi
    CorrData.push_back(dato(0.004, 0.002, 0.006));  // rp_dkstz_kpi
    CorrData.push_back(dato(0.004, 0.002, 0.006));  // rm_dkstz_kpi
    CorrData.push_back(dato(-0.029, 0.012, 0.021)); // afav_dkstz_k3pi
    CorrData.push_back(dato(0.019, 0.004, 0.007));  // rp_dkstz_k3pi
    CorrData.push_back(dato(0.015, 0.004, 0.007));  // rm_dkstz_k3pi
    CorrData.push_back(dato(0.062, 0.032, 0.021));  // acp_dkstz_kk
    CorrData.push_back(dato(1.00, 0.034, 0.016));   // rcp_dkstz_kk
    CorrData.push_back(dato(-0.001, 0.056, 0.021)); // acp_dkstz_pipi
    CorrData.push_back(dato(0.996, 0.057, 0.023));  // rcp_dkstz_pipi
    CorrData.push_back(dato(0.017, 0.044, 0.022));  // acp_dkstz_4pi
    CorrData.push_back(dato(1.010, 0.048, 0.033));  // rcp_dkstz_4pi

    // Correlation Matrix (stat):
    Corr.ResizeTo(12, 12);
    Corr(0, 0) = 1.;
    Corr(0, 1) = 0.02;
    Corr(0, 2) = -0.02;
    Corr(0, 3) = 0.0;
    Corr(0, 4) = 0.0;
    Corr(0, 5) = 0.;
    Corr(0, 6) = 0.;
    Corr(0, 7) = 0.;
    Corr(0, 8) = 0.;
    Corr(0, 9) = 0.;
    Corr(0, 10) = 0.;
    Corr(0, 11) = 0.;
    Corr(1, 1) = 1.;
    Corr(1, 2) = 0.1;
    Corr(1, 3) = 0.;
    Corr(1, 4) = 0.02;
    Corr(1, 5) = 0.03;
    Corr(1, 6) = 0.0;
    Corr(1, 7) = 0.;
    Corr(1, 8) = 0.;
    Corr(1, 9) = 0.;
    Corr(1, 10) = 0.;
    Corr(1, 11) = 0.;
    Corr(2, 2) = 1.;
    Corr(2, 3) = 0.;
    Corr(2, 4) = 0.02;
    Corr(2, 5) = 0.03;
    Corr(2, 6) = 0.0;
    Corr(2, 7) = 0.0;
    Corr(2, 8) = 0.0;
    Corr(2, 9) = 0.0;
    Corr(2, 10) = 0.0;
    Corr(2, 11) = 0.;
    Corr(3, 3) = 1.;
    Corr(3, 4) = 0.05;
    Corr(3, 5) = -0.04;
    Corr(3, 6) = 0.0;
    Corr(3, 7) = 0.0;
    Corr(3, 8) = 0.0;
    Corr(3, 9) = 0.0;
    Corr(3, 10) = 0.;
    Corr(3, 11) = 0.;
    Corr(4, 4) = 1.;
    Corr(4, 5) = 0.07;
    Corr(4, 6) = 0.;
    Corr(4, 7) = 0.;
    Corr(4, 8) = 0.;
    Corr(4, 9) = 0.;
    Corr(4, 10) = 0.;
    Corr(4, 11) = 0.01;
    Corr(5, 5) = 1.;
    Corr(5, 6) = 0.;
    Corr(5, 7) = 0.0;
    Corr(5, 8) = 0.;
    Corr(5, 9) = 0.;
    Corr(5, 10) = 0.;
    Corr(5, 11) = 0.01;
    Corr(6, 6) = 1.;
    Corr(6, 7) = -0.03;
    Corr(6, 8) = 0.;
    Corr(6, 9) = 0.;
    Corr(6, 10) = 0.;
    Corr(6, 11) = 0.;
    Corr(7, 7) = 1.;
    Corr(7, 8) = 0.;
    Corr(7, 9) = 0.05;
    Corr(7, 10) = 0.;
    Corr(7, 11) = 0.02;
    Corr(8, 8) = 1.;
    Corr(8, 9) = -0.02;
    Corr(8, 10) = 0.;
    Corr(8, 11) = 0.;
    Corr(9, 9) = 1.;
    Corr(9, 10) = 0.0;
    Corr(9, 11) = 0.01;
    Corr(10, 10) = 1.;
    Corr(10, 11) = -0.01;
    Corr(11, 11) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr);

    // Correlation Matrix (syst)
    Corr2.ResizeTo(12, 12);
    Corr2(0, 0) = 1.;
    Corr2(0, 1) = 0.;
    Corr2(0, 2) = -0.05;
    Corr2(0, 3) = 0.81;
    Corr2(0, 4) = 0.15;
    Corr2(0, 5) = -0.19;
    Corr2(0, 6) = 0.78;
    Corr2(0, 7) = -0.03;
    Corr2(0, 8) = 0.79;
    Corr2(0, 9) = -0.03;
    Corr2(0, 10) = 0.78;
    Corr2(0, 11) = -0.03;
    Corr2(1, 1) = 1.;
    Corr2(1, 2) = 0.15;
    Corr2(1, 3) = 0.14;
    Corr2(1, 4) = 0.07;
    Corr2(1, 5) = 0.06;
    Corr2(1, 6) = 0.12;
    Corr2(1, 7) = 0.05;
    Corr2(1, 8) = 0.11;
    Corr2(1, 9) = 0.03;
    Corr2(1, 10) = 0.05;
    Corr2(1, 11) = -0.01;
    Corr2(2, 2) = 1.;
    Corr2(2, 3) = -0.19;
    Corr2(2, 4) = 0.07;
    Corr2(2, 5) = 0.04;
    Corr2(2, 6) = -0.26;
    Corr2(2, 7) = 0.14;
    Corr2(2, 8) = -0.25;
    Corr2(2, 9) = 0.04;
    Corr2(2, 10) = -0.26;
    Corr2(2, 11) = -0.01;
    Corr2(3, 3) = 1.;
    Corr2(3, 4) = -0.01;
    Corr2(3, 5) = -0.09;
    Corr2(3, 6) = 0.77;
    Corr2(3, 7) = 0.0;
    Corr2(3, 8) = 0.77;
    Corr2(3, 9) = -0.01;
    Corr2(3, 10) = 0.76;
    Corr2(3, 11) = -0.30;
    Corr2(4, 4) = 1.;
    Corr2(4, 5) = 0.15;
    Corr2(4, 6) = 0.18;
    Corr2(4, 7) = 0.;
    Corr2(4, 8) = 0.17;
    Corr2(4, 9) = 0.01;
    Corr2(4, 10) = 0.19;
    Corr2(4, 11) = 0.04;
    Corr2(5, 5) = 1.;
    Corr2(5, 6) = -0.24;
    Corr2(5, 7) = 0.;
    Corr2(5, 8) = -0.24;
    Corr2(5, 9) = 0.01;
    Corr2(5, 10) = -0.23;
    Corr2(5, 11) = 0.04;
    Corr2(6, 6) = 1.;
    Corr2(6, 7) = -0.03;
    Corr2(6, 8) = 0.94;
    Corr2(6, 9) = -0.02;
    Corr2(6, 10) = 0.93;
    Corr2(6, 11) = 0.;
    Corr2(7, 7) = 1.;
    Corr2(7, 8) = -0.01;
    Corr2(7, 9) = 0.15;
    Corr2(7, 10) = -0.01;
    Corr2(7, 11) = -0.02;
    Corr2(8, 8) = 1.;
    Corr2(8, 9) = -0.03;
    Corr2(8, 10) = 0.95;
    Corr2(8, 11) = 0.;
    Corr2(9, 9) = 1.;
    Corr2(9, 10) = -0.02;
    Corr2(9, 11) = 0.;
    Corr2(10, 10) = 1.;
    Corr2(10, 11) = -0.03;
    Corr2(11, 11) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr2);

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("2401.17934Bs", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  }
  else if (comb == 3)
  { // this is if you are considering all the observables because Bd and Bs are correlated

    // GLW: D -> pipi, D -> KK, D -> 4pi; ADS: D -> Kpi, D -> K3pi
    // https://arxiv.org/pdf/2401.17934
    // Observables 24:

    // Bd part
    CorrData.clear();
    CorrData.push_back(dato(0.031, 0.017, 0.015));  // afav_dkstz_kpi
    CorrData.push_back(dato(0.069, 0.013, 0.005));  // rp_dkstz_kpi
    CorrData.push_back(dato(0.093, 0.013, 0.005));  // rm_dkstz_kpi
    CorrData.push_back(dato(-0.012, 0.018, 0.016)); // afav_dkstz_k3pi
    CorrData.push_back(dato(0.060, 0.014, 0.006));  // rp_dkstz_k3pi
    CorrData.push_back(dato(0.038, 0.014, 0.006));  // rm_dkstz_k3pi
    CorrData.push_back(dato(-0.047, 0.063, 0.015)); // acp_dkstz_kk
    CorrData.push_back(dato(0.811, 0.057, 0.017));  // rcp_dkstz_kk
    CorrData.push_back(dato(-0.034, 0.094, 0.016)); // acp_dkstz_pipi
    CorrData.push_back(dato(1.104, 0.111, 0.026));  // rcp_dkstz_pipi
    CorrData.push_back(dato(0.021, 0.087, 0.016));  // acp_dkstz_4pi
    CorrData.push_back(dato(0.882, 0.086, 0.033));  // rcp_dkstz_4pi

    // Bs part
    CorrData.push_back(dato(-0.009, 0.011, 0.020)); // afav_dkstz_kpi
    CorrData.push_back(dato(0.004, 0.002, 0.006));  // rp_dkstz_kpi
    CorrData.push_back(dato(0.004, 0.002, 0.006));  // rm_dkstz_kpi
    CorrData.push_back(dato(-0.029, 0.012, 0.021)); // afav_dkstz_k3pi
    CorrData.push_back(dato(0.019, 0.004, 0.007));  // rp_dkstz_k3pi
    CorrData.push_back(dato(0.015, 0.004, 0.007));  // rm_dkstz_k3pi
    CorrData.push_back(dato(0.062, 0.032, 0.021));  // acp_dkstz_kk
    CorrData.push_back(dato(1.00, 0.034, 0.016));   // rcp_dkstz_kk
    CorrData.push_back(dato(-0.001, 0.056, 0.021)); // acp_dkstz_pipi
    CorrData.push_back(dato(0.996, 0.057, 0.023));  // rcp_dkstz_pipi
    CorrData.push_back(dato(0.017, 0.044, 0.022));  // acp_dkstz_4pi
    CorrData.push_back(dato(1.010, 0.048, 0.033));  // rcp_dkstz_4pi

    // Correlation Matrix (stat):
    Corr.ResizeTo(24, 24);

    // B0d part 11X11
    Corr(0, 0) = 1.;
    Corr(0, 1) = 0.09;
    Corr(0, 2) = -0.12;
    Corr(0, 3) = 0.0;
    Corr(0, 4) = 0.0;
    Corr(0, 5) = 0.;
    Corr(0, 6) = 0.;
    Corr(0, 7) = -0.01;
    Corr(0, 8) = 0.;
    Corr(0, 9) = 0.;
    Corr(0, 10) = 0.;
    Corr(0, 11) = 0.;
    Corr(1, 1) = 1.;
    Corr(1, 2) = 0.06;
    Corr(1, 3) = 0.;
    Corr(1, 4) = 0.0;
    Corr(1, 5) = 0.0;
    Corr(1, 6) = 0.0;
    Corr(1, 7) = 0.03;
    Corr(1, 8) = 0.;
    Corr(1, 9) = 0.02;
    Corr(1, 10) = 0.;
    Corr(1, 11) = 0.;
    Corr(2, 2) = 1.;
    Corr(2, 3) = 0.;
    Corr(2, 4) = 0.;
    Corr(2, 5) = 0.0;
    Corr(2, 6) = 0.0;
    Corr(2, 7) = 0.03;
    Corr(2, 8) = 0.0;
    Corr(2, 9) = 0.02;
    Corr(2, 10) = 0.;
    Corr(2, 11) = 0.;
    Corr(3, 3) = 1.;
    Corr(3, 4) = 0.07;
    Corr(3, 5) = -0.05;
    Corr(3, 6) = 0.0;
    Corr(3, 7) = 0.0;
    Corr(3, 8) = 0.0;
    Corr(3, 9) = 0.0;
    Corr(3, 10) = 0.;
    Corr(3, 11) = -0.01;
    Corr(4, 4) = 1.;
    Corr(4, 5) = 0.08;
    Corr(4, 6) = 0.;
    Corr(4, 7) = 0.;
    Corr(4, 8) = 0.;
    Corr(4, 9) = 0.;
    Corr(4, 10) = 0.;
    Corr(4, 11) = 0.02;
    Corr(5, 5) = 1.;
    Corr(5, 6) = 0.0;
    Corr(5, 7) = 0.0;
    Corr(5, 8) = 0.0;
    Corr(5, 9) = 0.0;
    Corr(5, 10) = 0.;
    Corr(5, 11) = 0.01;
    Corr(6, 6) = 1.;
    Corr(6, 7) = 0.04;
    Corr(6, 8) = 0.0;
    Corr(6, 9) = 0.0;
    Corr(6, 10) = 0.0;
    Corr(6, 11) = 0.;
    Corr(7, 7) = 1.;
    Corr(7, 8) = 0.0;
    Corr(7, 9) = 0.05;
    Corr(7, 10) = 0.;
    Corr(7, 11) = 0.;
    Corr(8, 8) = 1.;
    Corr(8, 9) = 0.04;
    Corr(8, 10) = 0.0;
    Corr(8, 11) = 0.0;
    Corr(9, 9) = 1.;
    Corr(9, 10) = 0.0;
    Corr(9, 11) = 0.;
    Corr(10, 10) = 1.;
    Corr(10, 11) = 0.0;
    Corr(11, 11) = 1.;

    // Mixed Matrix B0d and B0s
    Corr(0, 12) = 0.;
    Corr(0, 13) = -0.02;
    Corr(0, 14) = 0.01;
    Corr(0, 15) = 0.;
    Corr(0, 16) = 0.;
    Corr(0, 17) = 0.;
    Corr(0, 18) = 0.;
    Corr(0, 19) = 0.;
    Corr(0, 20) = 0.;
    Corr(0, 21) = 0.;
    Corr(0, 22) = 0.;
    Corr(0, 23) = 0.;

    Corr(1, 12) = -0.02;
    Corr(1, 13) = 0.;
    Corr(1, 14) = 0.;
    Corr(1, 15) = 0.;
    Corr(1, 16) = 0.;
    Corr(1, 17) = 0.;
    Corr(1, 18) = 0.;
    Corr(1, 19) = 0.;
    Corr(1, 20) = 0.;
    Corr(1, 21) = 0.;
    Corr(1, 22) = 0.;
    Corr(1, 23) = -0.01;

    Corr(2, 12) = 0.02;
    Corr(2, 13) = 0.;
    Corr(2, 14) = 0.;
    Corr(2, 15) = 0.;
    Corr(2, 16) = 0.;
    Corr(2, 17) = 0.;
    Corr(2, 18) = 0.;
    Corr(2, 19) = 0.;
    Corr(2, 20) = 0.;
    Corr(2, 21) = 0.;
    Corr(2, 22) = 0.;
    Corr(2, 23) = 0.;

    Corr(3, 12) = 0.;
    Corr(3, 13) = 0.;
    Corr(3, 14) = 0.;
    Corr(3, 15) = 0.;
    Corr(3, 16) = -0.01;
    Corr(3, 17) = 0.02;
    Corr(3, 18) = 0.;
    Corr(3, 19) = 0.;
    Corr(3, 20) = 0.;
    Corr(3, 21) = 0.;
    Corr(3, 22) = 0.;
    Corr(3, 23) = 0.;

    Corr(4, 12) = 0.;
    Corr(4, 13) = 0.;
    Corr(4, 14) = 0.;
    Corr(4, 15) = -0.02;
    Corr(4, 16) = 0.;
    Corr(4, 17) = 0.;
    Corr(4, 18) = 0.;
    Corr(4, 19) = -0.01;
    Corr(4, 20) = 0.;
    Corr(4, 21) = 0.;
    Corr(4, 22) = 0.;
    Corr(4, 23) = 0.;

    Corr(5, 12) = 0.;
    Corr(5, 13) = 0.;
    Corr(5, 14) = 0.;
    Corr(5, 15) = 0.02;
    Corr(5, 16) = 0.;
    Corr(5, 17) = 0.;
    Corr(5, 18) = 0.;
    Corr(5, 19) = -0.01;
    Corr(5, 20) = 0.;
    Corr(5, 21) = 0.;
    Corr(5, 22) = 0.;
    Corr(5, 23) = 0.;

    Corr(6, 12) = 0.;
    Corr(6, 13) = 0.;
    Corr(6, 14) = 0.;
    Corr(6, 15) = 0.;
    Corr(6, 16) = 0.;
    Corr(6, 17) = 0.;
    Corr(6, 18) = 0.02;
    Corr(6, 19) = 0.;
    Corr(6, 20) = 0.;
    Corr(6, 21) = 0.;
    Corr(6, 22) = 0.;
    Corr(6, 23) = 0.;

    Corr(7, 12) = 0.;
    Corr(7, 13) = -0.01;
    Corr(7, 14) = -0.01;
    Corr(7, 15) = 0.;
    Corr(7, 16) = 0.;
    Corr(7, 17) = 0.;
    Corr(7, 18) = 0.;
    Corr(7, 19) = 0.03;
    Corr(7, 20) = 0.;
    Corr(7, 21) = 0.;
    Corr(7, 22) = 0.;
    Corr(7, 23) = 0.;

    Corr(8, 12) = 0.;
    Corr(8, 13) = 0.;
    Corr(8, 14) = 0.;
    Corr(8, 15) = 0.;
    Corr(8, 16) = 0.;
    Corr(8, 17) = 0.;
    Corr(8, 18) = 0.;
    Corr(8, 19) = 0.;
    Corr(8, 20) = 0.02;
    Corr(8, 21) = 0.;
    Corr(8, 22) = 0.;
    Corr(8, 23) = 0.;

    Corr(9, 12) = 0.;
    Corr(9, 13) = -0.01;
    Corr(9, 14) = -0.01;
    Corr(9, 15) = 0.;
    Corr(9, 16) = 0.;
    Corr(9, 17) = 0.;
    Corr(9, 18) = 0.;
    Corr(9, 19) = 0.;
    Corr(9, 20) = 0.;
    Corr(9, 21) = 0.03;
    Corr(9, 22) = 0.;
    Corr(9, 23) = 0.;

    Corr(10, 12) = 0.;
    Corr(10, 13) = 0.;
    Corr(10, 14) = 0.;
    Corr(10, 15) = 0.;
    Corr(10, 16) = 0.;
    Corr(10, 17) = 0.;
    Corr(10, 18) = 0.;
    Corr(10, 19) = 0.;
    Corr(10, 20) = 0.;
    Corr(10, 21) = 0.;
    Corr(10, 22) = 0.02;
    Corr(10, 23) = 0.;

    Corr(11, 12) = 0.;
    Corr(11, 13) = 0.;
    Corr(11, 14) = 0.;
    Corr(11, 15) = 0.;
    Corr(11, 16) = -0.01;
    Corr(11, 17) = -0.01;
    Corr(11, 18) = 0.;
    Corr(11, 19) = 0.;
    Corr(11, 20) = 0.0;
    Corr(11, 21) = 0.;
    Corr(11, 22) = 0.;
    Corr(11, 23) = 0.04;


    // Bs part
    Corr(12, 12) = 1.;
    Corr(12, 13) = 0.02;
    Corr(12, 14) = -0.02;
    Corr(12, 15) = 0.0;
    Corr(12, 16) = 0.0;
    Corr(12, 17) = 0.;
    Corr(12, 18) = 0.;
    Corr(12, 19) = 0.;
    Corr(12, 20) = 0.;
    Corr(12, 21) = 0.;
    Corr(12, 22) = 0.;
    Corr(12, 23) = 0.;
    Corr(13, 13) = 1.;
    Corr(13, 14) = 0.1;
    Corr(13, 15) = 0.;
    Corr(13, 16) = 0.02;
    Corr(13, 17) = 0.03;
    Corr(13, 18) = 0.0;
    Corr(13, 19) = 0.;
    Corr(13, 20) = 0.;
    Corr(13, 21) = 0.;
    Corr(13, 22) = 0.;
    Corr(13, 23) = 0.;
    Corr(14, 14) = 1.;
    Corr(14, 15) = 0.;
    Corr(14, 16) = 0.02;
    Corr(14, 17) = 0.03;
    Corr(14, 18) = 0.0;
    Corr(14, 19) = 0.0;
    Corr(14, 20) = 0.0;
    Corr(14, 21) = 0.0;
    Corr(14, 22) = 0.0;
    Corr(14, 23) = 0.;
    Corr(15, 15) = 1.;
    Corr(15, 16) = 0.05;
    Corr(15, 17) = -0.04;
    Corr(15, 18) = 0.0;
    Corr(15, 19) = 0.0;
    Corr(15, 20) = 0.0;
    Corr(15, 21) = 0.0;
    Corr(15, 22) = 0.;
    Corr(15, 23) = 0.;
    Corr(16, 16) = 1.;
    Corr(16, 17) = 0.07;
    Corr(16, 18) = 0.;
    Corr(16, 19) = 0.;
    Corr(16, 20) = 0.;
    Corr(16, 21) = 0.;
    Corr(16, 22) = 0.;
    Corr(16, 23) = 0.01;
    Corr(17, 17) = 1.;
    Corr(17, 18) = 0.;
    Corr(17, 19) = 0.0;
    Corr(17, 20) = 0.;
    Corr(17, 21) = 0.;
    Corr(17, 22) = 0.;
    Corr(17, 23) = 0.01;
    Corr(18, 18) = 1.;
    Corr(18, 19) = -0.03;
    Corr(18, 20) = 0.;
    Corr(18, 21) = 0.;
    Corr(18, 22) = 0.;
    Corr(18, 23) = 0.;
    Corr(19, 19) = 1.;
    Corr(19, 20) = 0.;
    Corr(19, 21) = 0.05;
    Corr(19, 22) = 0.;
    Corr(19, 23) = 0.02;
    Corr(20, 20) = 1.;
    Corr(20, 21) = -0.02;
    Corr(20, 22) = 0.;
    Corr(20, 23) = 0.;
    Corr(21, 21) = 1.;
    Corr(21, 22) = 0.0;
    Corr(21, 23) = 0.01;
    Corr(22, 22) = 1.;
    Corr(22, 23) = -0.01;
    Corr(23, 23) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr);



    // Correlation Matrix (syst)
    Corr2.ResizeTo(24, 24);

    // Bd part 11X11
    Corr2(0, 0) = 1.;
    Corr2(0, 1) = 0.;
    Corr2(0, 2) = -0.05;
    Corr2(0, 3) = 0.81;
    Corr2(0, 4) = 0.15;
    Corr2(0, 5) = -0.19;
    Corr2(0, 6) = 0.78;
    Corr2(0, 7) = -0.03;
    Corr2(0, 8) = 0.79;
    Corr2(0, 9) = -0.03;
    Corr2(0, 10) = 0.78;
    Corr2(0, 11) = -0.03;
    Corr2(1, 1) = 1.;
    Corr2(1, 2) = 0.15;
    Corr2(1, 3) = 0.14;
    Corr2(1, 4) = 0.07;
    Corr2(1, 5) = 0.06;
    Corr2(1, 6) = 0.12;
    Corr2(1, 7) = 0.05;
    Corr2(1, 8) = 0.11;
    Corr2(1, 9) = 0.03;
    Corr2(1, 10) = 0.05;
    Corr2(1, 11) = -0.01;
    Corr2(2, 2) = 1.;
    Corr2(2, 3) = -0.19;
    Corr2(2, 4) = 0.07;
    Corr2(2, 5) = 0.04;
    Corr2(2, 6) = -0.26;
    Corr2(2, 7) = 0.14;
    Corr2(2, 8) = -0.25;
    Corr2(2, 9) = 0.04;
    Corr2(2, 10) = -0.26;
    Corr2(2, 11) = -0.01;
    Corr2(3, 3) = 1.;
    Corr2(3, 4) = -0.01;
    Corr2(3, 5) = -0.09;
    Corr2(3, 6) = 0.77;
    Corr2(3, 7) = 0.0;
    Corr2(3, 8) = 0.77;
    Corr2(3, 9) = -0.01;
    Corr2(3, 10) = 0.76;
    Corr2(3, 11) = -0.30;
    Corr2(4, 4) = 1.;
    Corr2(4, 5) = 0.15;
    Corr2(4, 6) = 0.18;
    Corr2(4, 7) = 0.;
    Corr2(4, 8) = 0.17;
    Corr2(4, 9) = 0.01;
    Corr2(4, 10) = 0.19;
    Corr2(4, 11) = 0.04;
    Corr2(5, 5) = 1.;
    Corr2(5, 6) = -0.24;
    Corr2(5, 7) = 0.;
    Corr2(5, 8) = -0.24;
    Corr2(5, 9) = 0.01;
    Corr2(5, 10) = -0.23;
    Corr2(5, 11) = 0.04;
    Corr2(6, 6) = 1.;
    Corr2(6, 7) = -0.03;
    Corr2(6, 8) = 0.94;
    Corr2(6, 9) = -0.02;
    Corr2(6, 10) = 0.93;
    Corr2(6, 11) = 0.;
    Corr2(7, 7) = 1.;
    Corr2(7, 8) = -0.01;
    Corr2(7, 9) = 0.15;
    Corr2(7, 10) = -0.01;
    Corr2(7, 11) = -0.02;
    Corr2(8, 8) = 1.;
    Corr2(8, 9) = -0.03;
    Corr2(8, 10) = 0.95;
    Corr2(8, 11) = 0.;
    Corr2(9, 9) = 1.;
    Corr2(9, 10) = -0.02;
    Corr2(9, 11) = 0.;
    Corr2(10, 10) = 1.;
    Corr2(10, 11) = -0.03;
    Corr2(11, 11) = 1.;

    // Mixed part
    Corr2(0, 12) = 1.;
    Corr2(0, 13) = 0.;
    Corr2(0, 14) = -0.05;
    Corr2(0, 15) = 0.81;
    Corr2(0, 16) = 0.15;
    Corr2(0, 17) = -0.19;
    Corr2(0, 18) = 0.78;
    Corr2(0, 19) = -0.03;
    Corr2(0, 20) = 0.79;
    Corr2(0, 21) = -0.03;
    Corr2(0, 22) = 0.78;
    Corr2(0, 23) = -0.03;
    Corr2(1, 12) = 0.;
    Corr2(1, 13) = 1.;
    Corr2(1, 14) = 0.15;
    Corr2(1, 15) = 0.14;
    Corr2(1, 16) = 0.07;
    Corr2(1, 17) = 0.06;
    Corr2(1, 18) = 0.12;
    Corr2(1, 19) = 0.05;
    Corr2(1, 20) = 0.11;
    Corr2(1, 21) = 0.03;
    Corr2(1, 22) = 0.05;
    Corr2(1, 23) = -0.01;
    Corr2(2, 12) = -0.05;
    Corr2(2, 13) = 0.15;
    Corr2(2, 14) = 1.;
    Corr2(2, 15) = -0.19;
    Corr2(2, 16) = 0.07;
    Corr2(2, 17) = 0.04;
    Corr2(2, 18) = -0.26;
    Corr2(2, 19) = 0.14;
    Corr2(2, 20) = -0.25;
    Corr2(2, 21) = 0.04;
    Corr2(2, 22) = -0.26;
    Corr2(2, 23) = -0.01;
    Corr2(3, 12) = 0.81;
    Corr2(3, 13) = 0.14;
    Corr2(3, 14) = -0.19;
    Corr2(3, 15) = 1.;
    Corr2(3, 16) = -0.01;
    Corr2(3, 17) = -0.09;
    Corr2(3, 18) = 0.77;
    Corr2(3, 19) = 0.0;
    Corr2(3, 20) = 0.77;
    Corr2(3, 21) = -0.01;
    Corr2(3, 22) = 0.76;
    Corr2(3, 23) = -0.30;
    Corr2(4, 12) = 0.15;
    Corr2(4, 13) = 0.07;
    Corr2(4, 14) = 0.07;
    Corr2(4, 15) = -0.01;
    Corr2(4, 16) = 1.;
    Corr2(4, 17) = 0.15;
    Corr2(4, 18) = 0.18;
    Corr2(4, 19) = 0.;
    Corr2(4, 20) = 0.17;
    Corr2(4, 21) = 0.01;
    Corr2(4, 22) = 0.19;
    Corr2(4, 23) = 0.04;
    Corr2(5, 12) = -0.19;
    Corr2(5, 13) = 0.06;
    Corr2(5, 14) = 0.04;
    Corr2(5, 15) = -0.09;
    Corr2(5, 16) = 0.15;
    Corr2(5, 17) = 1.;
    Corr2(5, 18) = -0.24;
    Corr2(5, 19) = 0.;
    Corr2(5, 20) = -0.24;
    Corr2(5, 21) = 0.01;
    Corr2(5, 22) = -0.23;
    Corr2(5, 23) = 0.04;
    Corr2(6, 12) = 0.78;
    Corr2(6, 13) = 0.12;
    Corr2(6, 14) = -0.26;
    Corr2(6, 15) = 0.77;
    Corr2(6, 16) = 0.18;
    Corr2(6, 17) = -0.24;
    Corr2(6, 18) = 1.;
    Corr2(6, 19) = -0.03;
    Corr2(6, 20) = 0.94;
    Corr2(6, 21) = -0.02;
    Corr2(6, 22) = 0.93;
    Corr2(6, 23) = 0.;
    Corr2(7, 12) = -0.03;
    Corr2(7, 13) = 0.05;
    Corr2(7, 14) = 0.14;
    Corr2(7, 15) = 0.;
    Corr2(7, 16) = 0.;
    Corr2(7, 17) = 0.;
    Corr2(7, 18) = -0.03;
    Corr2(7, 19) = 1.;
    Corr2(7, 20) = -0.01;
    Corr2(7, 21) = 0.15;
    Corr2(7, 22) = -0.01;
    Corr2(7, 23) = -0.02;
    Corr2(8, 12) = 0.79;
    Corr2(8, 13) = 0.11;
    Corr2(8, 14) = -0.25;
    Corr2(8, 15) = 0.77;
    Corr2(8, 16) = 0.17;
    Corr2(8, 17) = -0.24;
    Corr2(8, 18) = 0.94;
    Corr2(8, 19) = -0.01;
    Corr2(8, 20) = 1.;
    Corr2(8, 21) = -0.03;
    Corr2(8, 22) = 0.95;
    Corr2(8, 23) = 0.;
    Corr2(9, 12) = -0.03;
    Corr2(9, 13) = 0.03;
    Corr2(9, 14) = 0.04;
    Corr2(9, 15) = -0.01;
    Corr2(9, 16) = 0.01;
    Corr2(9, 17) = 0.01;
    Corr2(9, 18) = -0.02;
    Corr2(9, 19) = 0.15;
    Corr2(9, 20) = -0.03;
    Corr2(9, 21) = 1.;
    Corr2(9, 22) = -0.02;
    Corr2(9, 23) = 0.;
    Corr2(10, 12) = 0.78;
    Corr2(10, 13) = 0.05;
    Corr2(10, 14) = -0.26;
    Corr2(10, 15) = 0.76;
    Corr2(10, 16) = 0.19;
    Corr2(10, 17) = -0.23;
    Corr2(10, 18) = 0.93;
    Corr2(10, 19) = -0.01;
    Corr2(10, 20) = 0.95;
    Corr2(10, 21) = -0.02;
    Corr2(10, 22) = 1.;
    Corr2(10, 23) = -0.03;
    Corr2(11, 12) = -0.03;
    Corr2(11, 13) = -0.01;
    Corr2(11, 14) = -0.01;
    Corr2(11, 15) = -0.30;
    Corr2(11, 16) = 0.04;
    Corr2(11, 17) = 0.04;
    Corr2(11, 18) = 0.;
    Corr2(11, 19) = -0.02;
    Corr2(11, 20) = 0.;
    Corr2(11, 21) = 0.;
    Corr2(11, 22) = -0.03;
    Corr2(11, 23) = 1.;


    // Bs Part
    Corr2(12, 12) = 1.;
    Corr2(12, 13) = 0.;
    Corr2(12, 14) = -0.05;
    Corr2(12, 15) = 0.81;
    Corr2(12, 16) = 0.15;
    Corr2(12, 17) = -0.19;
    Corr2(12, 18) = 0.78;
    Corr2(12, 19) = -0.03;
    Corr2(12, 20) = 0.79;
    Corr2(12, 21) = -0.03;
    Corr2(12, 22) = 0.78;
    Corr2(12, 23) = -0.03;
    Corr2(13, 13) = 1.;
    Corr2(13, 14) = 0.15;
    Corr2(13, 15) = 0.14;
    Corr2(13, 16) = 0.07;
    Corr2(13, 17) = 0.06;
    Corr2(13, 18) = 0.12;
    Corr2(13, 19) = 0.05;
    Corr2(13, 20) = 0.11;
    Corr2(13, 21) = 0.03;
    Corr2(13, 22) = 0.05;
    Corr2(13, 23) = -0.01;
    Corr2(14, 14) = 1.;
    Corr2(14, 15) = -0.19;
    Corr2(14, 16) = 0.07;
    Corr2(14, 17) = 0.04;
    Corr2(14, 18) = -0.26;
    Corr2(14, 19) = 0.14;
    Corr2(14, 20) = -0.25;
    Corr2(14, 21) = 0.04;
    Corr2(14, 22) = -0.26;
    Corr2(14, 23) = -0.01;
    Corr2(15, 15) = 1.;
    Corr2(15, 16) = -0.01;
    Corr2(15, 17) = -0.09;
    Corr2(15, 18) = 0.77;
    Corr2(15, 19) = 0.0;
    Corr2(15, 20) = 0.77;
    Corr2(15, 21) = -0.01;
    Corr2(15, 22) = 0.76;
    Corr2(15, 23) = -0.30;
    Corr2(16, 16) = 1.;
    Corr2(16, 17) = 0.15;
    Corr2(16, 18) = 0.18;
    Corr2(16, 19) = 0.;
    Corr2(16, 20) = 0.17;
    Corr2(16, 21) = 0.01;
    Corr2(16, 22) = 0.19;
    Corr2(16, 23) = 0.04;
    Corr2(17, 17) = 1.;
    Corr2(17, 18) = -0.24;
    Corr2(17, 19) = 0.;
    Corr2(17, 20) = -0.24;
    Corr2(17, 21) = 0.01;
    Corr2(17, 22) = -0.23;
    Corr2(17, 23) = 0.04;
    Corr2(18, 18) = 1.;
    Corr2(18, 19) = -0.03;
    Corr2(18, 20) = 0.94;
    Corr2(18, 21) = -0.02;
    Corr2(18, 22) = 0.93;
    Corr2(18, 23) = 0.;
    Corr2(19, 19) = 1.;
    Corr2(19, 20) = -0.01;
    Corr2(19, 21) = 0.15;
    Corr2(19, 22) = -0.01;
    Corr2(19, 23) = -0.02;
    Corr2(20, 20) = 1.;
    Corr2(20, 21) = -0.03;
    Corr2(20, 22) = 0.95;
    Corr2(20, 23) = 0.;
    Corr2(21, 21) = 1.;
    Corr2(21, 22) = -0.02;
    Corr2(21, 23) = 0.;
    Corr2(22, 22) = 1.;
    Corr2(22, 23) = -0.03;
    Corr2(23, 23) = 1.;
    SymmetrizeUpperTriangularMatrix(Corr2);

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("2401.17934", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  }

  //---------------------------------------------------------------------------------------------------------------

  //-------------------------------------- B0s -> DsK -------------------------------------------------------------------------

  // Time dependent B0s, Ds -> KKpi, pipipi, LHCb Run 1 full
  // LHCb-PAPER-2017-04 updated
  // Observables 5:
  CorrData.clear();
  CorrData.push_back(dato(0.75, 0.14, 0.04)); // c_dsk
  CorrData.push_back(dato(0.38, 0.28, 0.15)); // d_dsk
  CorrData.push_back(dato(0.30, 0.28, 0.15)); // db_dsk
  CorrData.push_back(dato(-0.53, 0.21, 0.06)); // s_dsk
  CorrData.push_back(dato(-0.45, 0.20, 0.06)); // sb_dsk

  // Correlation Matrix (stat):
  Corr.ResizeTo(5, 5);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.114;
  Corr(0, 2) = 0.098;
  Corr(0, 3) = 0.018;
  Corr(0, 4) = -0.054;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.546;
  Corr(1, 3) = -0.088;
  Corr(1, 4) = -0.024;
  Corr(2, 2) = 1.;
  Corr(2, 3) = -0.051;
  Corr(2, 4) = -0.024;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.001;
  Corr(4, 4) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(5, 5);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.07;
  Corr2(0, 2) = 0.05;
  Corr2(0, 3) = 0.04;
  Corr2(0, 4) = -0.01;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.53;
  Corr2(1, 3) = 0.02;
  Corr2(1, 4) = 0.02;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.03;
  Corr2(2, 4) = 0.03;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = 0.02;
  Corr2(4, 4) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("BSDSKRun1", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // Time dependent B0s, Ds -> KKpi, pipipi, LHCb Run 1 full
  // LHCb-CONF-2024-020 updated
  // Observables 5:
  CorrData.clear();
  CorrData.push_back(dato(0.791, 0.061, 0.022)); // c_dsk
  CorrData.push_back(dato(-0.051, 0.134, 0.058)); // d_dsk
  CorrData.push_back(dato(-0.303, 0.125, 0.055)); // db_dsk
  CorrData.push_back(dato(-0.571, 0.084, 0.023)); // s_dsk
  CorrData.push_back(dato(-0.503, 0.084, 0.025)); // sb_dsk

  // Correlation Matrix (stat):
  Corr.ResizeTo(5, 5);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.134;
  Corr(0, 2) = 0.130;
  Corr(0, 3) = 0.039;
  Corr(0, 4) = 0.022;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.501;
  Corr(1, 3) = -0.108;
  Corr(1, 4) = -0.036;
  Corr(2, 2) = 1.;
  Corr(2, 3) = -0.056;
  Corr(2, 4) = -0.067;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.006;
  Corr(4, 4) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(5, 5);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.008;
  Corr2(0, 2) = 0.012;
  Corr2(0, 3) = -0.08;
  Corr2(0, 4) = -0.246;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.878;
  Corr2(1, 3) = 0.004;
  Corr2(1, 4) = -0.022;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = -0.002;
  Corr2(2, 4) = -0.022;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = 0.085;
  Corr2(4, 4) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("BSDSKRun2", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //---------------------------------------------------------------------------------------------------------------

  //---------------------------------------------------------------------------------------------------------------

  // 12. PDF: dskpipi (UID11) https://link.springer.com/article/10.1007/JHEP03(2021)137
  // https://arxiv.org/pdf/2011.12041
  // Observables 5:
  CorrData.clear();
  CorrData.push_back(dato(0.631, 0.096, 0.032));  // c_dskpipi
  CorrData.push_back(dato(-0.334, 0.232, 0.097)); // d_dskpipi
  CorrData.push_back(dato(-0.695, 0.215, 0.081)); // db_dskpipi
  CorrData.push_back(dato(-0.424, 0.135, 0.033)); // s_dskpipi
  CorrData.push_back(dato(-0.463, 0.134, 0.031)); // sb_dskpipi

  // Correlation Matrix (stat):
  Corr.ResizeTo(5, 5);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.135;
  Corr(0, 2) = 0.168;
  Corr(0, 3) = 0.04;
  Corr(0, 4) = 0.042;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.504;
  Corr(1, 3) = -0.082;
  Corr(1, 4) = -0.052;
  Corr(2, 2) = 1.;
  Corr(2, 3) = -0.036;
  Corr(2, 4) = -0.113;
  Corr(3, 3) = 1.;
  Corr(3, 4) = 0.014;
  Corr(4, 4) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(5, 5);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.039;
  Corr2(0, 2) = 0.042;
  Corr2(0, 3) = -0.136;
  Corr2(0, 4) = -0.006;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.208;
  Corr2(1, 3) = -0.031;
  Corr2(1, 4) = -0.023;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = -0.021;
  Corr2(2, 4) = -0.035;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = -0.278;
  Corr2(4, 4) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID11", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //---------------------------------------------------------------------------------------------------------------

  //-------------------------------------  Other useful inputs for the Time Dependent B0 Decay and mixing parameters -------------------------------------------------------------------------

  // phis
  // https://hflav-eos.web.cern.ch/hflav-eos/osc/HFLAV_2024/
  // Observables 1:
    meas.insert(pair<string, dato>("UID26", dato(-0.060, 0.014, 0.))); // phis = - 2 beta_s J/psi-phi only
  // phis
  // https://arxiv.org/pdf/2308.01468
  // Observables 1:
  // meas.insert(pair<string, dato>("UID26", dato(-0.031, 0.018, 0.))); // phis = - 2 beta_s

}
// ---------------------------------------------------------

// ---------------------------------------------------------
void MixingModel::Add_time_dependent_Dmeas()
{

  vector<dato> CorrData;
  TMatrixDSym Corr, Corr2, Corr3;

  //-------------------------------------- D -> hh  -------------------------------------------------------------------------

  // Delta ACP; Acp(KK); Run1 semileptonic tagging
  // https://arxiv.org/pdf/1405.2797
  // Observables 4:
  meas.insert(pair<string, dato>("1405.2797_tKKOverTauD_DAcp", dato(1.082, 0.001, 0.004))); // tKK/tau_D DAcp
  meas.insert(pair<string, dato>("1405.2797_tpipiOverTauD_DAcp", dato(1.068, 0.001, 0.004))); // tpipi/tau_D DAcp
  meas.insert(pair<string, dato>("1405.2797_tKKOverTauD_Acp", dato(1.051, 0.001, 0.004))); // tKK/tau_D Acp
  CorrData.clear();
  CorrData.push_back(dato(14e-4, 16e-4, 8e-4)); // Delta Acp
  CorrData.push_back(dato(-6e-4, 15e-4, 10e-4)); // Acp KK
  // ACP(KK); Run1 Hadronic tagging; Correlated due to the removing of the detection asymmetry
  // https://arxiv.org/pdf/1610.09476
  // Observables 1:
  CorrData.push_back(dato(14e-4, 15e-4, 10e-4)); // Acp(KK) pitagged

  Corr.ResizeTo(3, 3);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.23;
  Corr(0, 2) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.36;
  Corr(2, 2) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  Corr2.ResizeTo(3, 3);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.40;
  Corr2(0, 2) = 0.;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 1.;
  Corr2(2, 2) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("1405.2797_1610.09476_Acp", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // tOverTauD; Run1 Hadronic tagging;
  // https://arxiv.org/pdf/1610.09476
  // Observables 1:
  meas.insert(pair<string, dato>("1610.09476_tKKOverTauD", dato(2.2390, 0.0007, 0.0187))); // tKK/tau_D


  // Delta ACP; Run1 Hadronic tagging
  // https://arxiv.org/pdf/1602.03160
  // Observables 3:
  meas.insert(pair<string, dato>("1602.03160_DeltaACPpitagged", dato(-0.0010, 0.0008, 0.0003))); // DeltaACPpitagged
  meas.insert(pair<string, dato>("1602.03160_tKKOverTauD", dato(2.1524, 0.0005, 0.0162))); // tKKpitaggedOverTauD
  meas.insert(pair<string, dato>("1602.03160_tpipiOverTauD",dato(2.0371, 0.0005, 0.0151))); // tpipipitaggedOverTauD


  // Delta ACP; Run2 Hadronic and Semileptonic tagging
  // https://arxiv.org/pdf/1903.08726
  // Observables 6:
  meas.insert(pair<string, dato>("DeltaACPpitagged", dato(-0.00182, 0.00032, 0.00009))); // DeltaACPpitagged
  meas.insert(pair<string, dato>("DeltaACPmutagged", dato(-0.0009, 0.0008, 0.0005))); // DeltaACPmutagged
  meas.insert(pair<string, dato>("tavepitaggedOverTauD", dato(1.74, 0.1, 0.))); // tpitaggedOverTauD
  meas.insert(pair<string, dato>("tavemutaggedOverTauD", dato(1.21, 0.01, 0.))); // tmutaggedOverTauD
  meas.insert(pair<string, dato>("DeltatmutaggedOverTauD",dato(-0.003, 0.001, 0.))); // DeltatmutaggedOverTauD
  // tOverTauD; Run 2 Acp(KK); Correlated through reconstructed mean decay times
  // https://arxiv.org/pdf/2209.03179
  // Observables 2:
  CorrData.clear();
  CorrData.push_back(dato(7.315e-13, 0.02e-13, 0.)); // tKKCDp
  CorrData.push_back(dato(6.868e-13, 0.014e-13, 0.)); // tKKCDs
  CorrData.push_back(dato(0.135, 0.002, 0.));    // DeltatpitaggedOverTauD 1903.08726

  // Correlation Matrix (stat):
  Corr.ResizeTo(3, 3);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.74;
  Corr(0, 2) = 0.23;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.25;
  Corr(2, 2) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(3, 3);
  Corr2.UnitMatrix();

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("tausforDACP", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // Run 2 Acp(KK)
  // https://arxiv.org/pdf/2209.03179
  // Observables 2:
  CorrData.clear();
  CorrData.push_back(dato(13.6e-4, 8.8e-4, 1.6e-4)); // ACPKKDp
  CorrData.push_back(dato(2.8e-4, 6.7e-4, 2.0e-4)); // ACPKKDs

  Corr.ResizeTo(2, 2);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.05;
  Corr(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  Corr2.ResizeTo(2, 2);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.28;
  Corr2(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("ACPKK", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // yCP - yCP(Kpi)
  // HFLAV combo: https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM23/results_mixing.html#kkpipi including latest https://arxiv.org/abs/2202.09106
  meas.insert(pair<string, dato>("UID28", dato(0.00697, 0.00028, 0.))); // ytilde_cp


  // ( DY(KK) + DY(pipi) )/2 LHCb combo Run1 + Run2
  // https://arxiv.org/pdf/2105.09889
  // Observables 1:
  meas.insert(pair<string, dato>("UID29", dato(-0.00019, 0.00013, 0.00004))); // DY // Average of KK and pipi


  // ( DY(KK) - DY(pipi) ) LHCb combo Run1 + Run2
  // https://arxiv.org/pdf/2105.09889
  // Observables 1:
  meas.insert(pair<string, dato>("DYKKmDYpipi", dato(3.3e-4, 2.7e-4, 0.2e-4))); // DYKKmDYpipi


  // Agamma(KK) and Agamma(pipi) Full CDF
  // https://arxiv.org/pdf/1410.5435
  // Observables 2:
  meas.insert(pair<string, dato>("1410.5435_AGammaKK", dato(-0.19e-2, 0.15e-2, 0.04e-2))); // AgammaKK
  meas.insert(pair<string, dato>("1410.5435_AGammapipi", dato(-0.01e-2, 0.18e-2, 0.03e-2))); // Agammapipi


  // tOverTauD CDF
  // https://arxiv.org/pdf/1111.5023
  // Observables 2:
  meas.insert(pair<string, dato>("1111.5023_tKKOverTauD_CDF", dato(2.65, 0.03, 0.))); // tKK/tau_D Acp
  meas.insert(pair<string, dato>("1111.5023_tpipiOverTauD_CDF", dato(2.40, 0.03, 0.))); // tpipi/tau_D Acp


  // Acp(KK), Acp(pipi)
  // https://arxiv.org/pdf/1208.2517
  // Observables 2
  meas.insert(pair<string, dato>("1208.2517_AcpKK_CDF", dato(-0.32e-2, 0.21e-2, 0.))); // Acp(KK)
  meas.insert(pair<string, dato>("1208.2517_Acppipi_CDF", dato(0.31e-2, 0.22e-2, 0.))); // Acp(pipi)


  // Acp(KK), Acp(pipi) Babar
  // https://arxiv.org/pdf/0709.2715
  // Observables 2
  meas.insert(pair<string, dato>("0709.2715_AcpKK_Babar", dato(0., 0.34e-2, 0.13e-2))); // Acp(KK)
  meas.insert(pair<string, dato>("0709.2715_Acppipi_Babar", dato(-0.24e-2, 0.52e-2, 0.22e-2))); // Acp(pipi)


  // Acp(KK), Acp(pipi) Belle
  // https://arxiv.org/pdf/0807.0148
  // Observables 2
  meas.insert(pair<string, dato>("0807.0148_AcpKK_Belle", dato(-0.43e-2, 0.30e-2, 0.11e-2))); // Acp(KK)
  meas.insert(pair<string, dato>("0807.0148_Acppipi_Belle", dato(0.43e-2, 0.52e-2, 0.12e-2))); // Acp(pipi)

  //---------------------------------------------------------------------------------------------------------------


  //-------------------------------------- D -> Kpi -------------------------------------------------------------------------

  // rD_kpi, (x'pm_Kpi)^2, y'pm_Kpi; LHCb Run1 + Run2 duble tag (DT) semileptonic
  // https://arxiv.org/pdf/1611.06143 + https://arxiv.org/pdf/2501.11635
  // Observables 6:
  CorrData.clear();
  CorrData.push_back(dato(347.e-5, 5.1e-5));    // Rd
  CorrData.push_back(dato(5.5e-3, 1.5e-3));     // cKpi
  CorrData.push_back(dato(1.2e-5, 2.4e-5));  //cpKpi
  CorrData.push_back(dato(0.9e-2, 1.5e-2));    // AD
  CorrData.push_back(dato(-1.3e-3, 1.5e-3)); // DcKpi
  CorrData.push_back(dato(1.2e-5, 2.4e-5)); // DcpKpi

  // Correlation Matrix (stat):
  Corr.ResizeTo(6, 6);
  Corr(0, 0) = 1.;
  Corr(0, 1) = -75.2e-2;
  Corr(0, 2) = 60.4e-2;
  Corr(0, 3) = -1.4e-2;
  Corr(0, 4) = 1.3e-2;
  Corr(0, 5) = -0.5e-2;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -92.5e-2;
  Corr(1, 3) = 1.9e-2;
  Corr(1, 4) = -3.5e-2;
  Corr(1, 5) = 2.4e-2;
  Corr(2, 2) = 1.;
  Corr(2, 3) = -1.e-2;
  Corr(2, 4) = 2.4e-2;
  Corr(2, 5) = -2.4e-2;
  Corr(3, 3) = 1;
  Corr(3, 4) = -74.1e-2;
  Corr(3, 5) = 59.7e-2;
  Corr(4, 4) = 1.;
  Corr(4, 5) = -92.5e-2;
  Corr(5, 5) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID30", CorrelatedGaussianObservables(CorrData, Corr)));


  // rD_kpi, c'pm_Kpi, Dc'pm_kpi; LHCb Run1+Run2, pi-tagged
  // https://arxiv.org/abs/2407.18001
  // Observables 9:
  CorrData.clear();

  // Fit in Appendix B
  CorrData.push_back(dato(0.003427, 0.000019, 0.0)); // RKpi
  CorrData.push_back(dato(52.8e-4, 3.3e-4, 0.));    // CKpi
  CorrData.push_back(dato(12.0e-6, 3.5e-6, 0.0));   // CpKpi
  CorrData.push_back(dato(-8.2e-3, 5.9e-3, 0.0));   // AtildeKpi
  CorrData.push_back(dato(3.2e-4, 3.6e-4, 0.0));    // DCtildeKpi
  CorrData.push_back(dato(-2.0e-6, 3.8e-6, 0.0));   // DCtildepKpi
  CorrData.push_back(dato(-0.9e-2, 2.0e-2, 0.0));   // AKpi
  CorrData.push_back(dato(-0.1e-3, 1.0e-3, 0.0));   // DeltaCKpi
  CorrData.push_back(dato(4.6e-6, 9.8e-6, 0.0));    // DeltaCpKpi

  // Correlation Matrix (stat):
  Corr.ResizeTo(9, 9);
  Corr(0, 0) = 1.;
  Corr(0, 1) = -92.7e-2;
  Corr(0, 2) = 80.3e-2;
  Corr(0, 3) = 0.8e-2;
  Corr(0, 4) = -0.7e-2;
  Corr(0, 5) = 0.e-2;
  Corr(0, 6) = 0.3e-2;
  Corr(0, 7) = -0.2e-2;
  Corr(0, 8) = 0.2e-2;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -94.3e-2;
  Corr(1, 3) = -1.4e-2;
  Corr(1, 4) = 1.3e-2;
  Corr(1, 5) = -0.6e-2;
  Corr(1, 6) = -0.5e-2;
  Corr(1, 7) = 0.4e-2;
  Corr(1, 8) = -0.4e-2;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.7e-2;
  Corr(2, 4) = -0.6e-2;
  Corr(2, 5) = 0.e-2;
  Corr(2, 6) = 0.3e-2;
  Corr(2, 7) = -0.3e-2;
  Corr(2, 8) = 0.3e-2;
  Corr(3, 3) = 1.;
  Corr(3, 4) = -93.4e-2;
  Corr(3, 5) = 81.e-2;
  Corr(3, 6) = 0.e-2;
  Corr(3, 7) = 0.e-2;
  Corr(3, 8) = 0.e-2;
  Corr(4, 4) = 1.;
  Corr(4, 5) = -94.3e-2;
  Corr(4, 6) = 0.e-2;
  Corr(4, 7) = 0.e-2;
  Corr(4, 8) = 0.e-2;
  Corr(5, 5) = 1.;
  Corr(5, 6) = 0.e-2;
  Corr(5, 7) = 0.e-2;
  Corr(5, 8) = 0.e-2;
  Corr(6, 6) = 1.;
  Corr(6, 7) = -93.8e-2;
  Corr(6, 8) = 81.1e-2;
  Corr(7, 7) = 1.;
  Corr(7, 8) = -94.3e-2;
  Corr(8, 8) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(9, 9);
  Corr2.UnitMatrix();

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("2407.18001", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // Asymmetries D -> Kpi   BESIII
  // https://arxiv.org/abs/2208.09402
  // First measurements
  // Observables 2:
  CorrData.clear();
  CorrData.push_back(dato(0.132, 0.011, 0.007)); // A_kpi
  CorrData.push_back(dato(0.130, 0.012, 0.008)); // A_kpi ^ pipipi0

  // Correlation Matrix (stat):
  Corr.ResizeTo(2, 2);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.38;
  Corr(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(2, 2);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.16;
  Corr2(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("BESIII_Adk", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  // Second measurements
  // Observables 2:
  CorrData.clear();
  CorrData.push_back(dato(-0.0562, 0.0081, 0.0050, 0.0010)); // rD_kpi cos()
  CorrData.push_back(dato(-0.011, 0.012, 0.007, 0.003));   // rD_kpi sin()

  // Correlation Matrix (stat):
  Corr.ResizeTo(2, 2);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.02;
  Corr(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("BESIII_rDkpi_polar", CorrelatedGaussianObservables(CorrData, Corr)));


  //-------------------------------------- D -> K0spipi -------------------------------------------------------------------------

  // xcp, ycp, dx, dy; LHCb Run1; Prompt + semileptonic
  // https://arxiv.org/pdf/1903.03074
  // Observables 4:
  CorrData.clear();
  CorrData.push_back(dato(0.0027, 0.0016, 0.0004));    // xcp
  CorrData.push_back(dato(0.00740, 0.0036, 0.0011));   // ycp
  CorrData.push_back(dato(-0.00053, 0.0007, 0.00022)); // dx
  CorrData.push_back(dato(0.0006, 0.0016, 0.0003));    // dy

  // Correlation Matrix (stat):
  Corr.ResizeTo(4, 4);
  Corr(0, 0) = 1.;
  Corr(0, 1) = -0.17;
  Corr(0, 2) = 0.04;
  Corr(0, 3) = -0.02;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -0.03;
  Corr(1, 3) = 0.01;
  Corr(2, 2) = 1.;
  Corr(2, 3) = -0.13;
  Corr(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(4, 4);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.15;
  Corr2(0, 2) = 0.01;
  Corr2(0, 3) = -0.02;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = -0.05;
  Corr2(1, 3) = -0.03;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.14;
  Corr2(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID14", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // https://arxiv.org/pdf/2208.06512.pdf
  // Observables 4:
  CorrData.clear();
  CorrData.push_back(dato(0.00401, 0.00045, 0.0002));  // xcp
  CorrData.push_back(dato(0.00551, 0.00116, 0.00059));   // ycp
  CorrData.push_back(dato(-0.00029, 0.00018, 0.00001)); // dx
  CorrData.push_back(dato(0.00031, 0.00035, 0.00013)); // dy

  // Correlation Matrix (stat):
  Corr.ResizeTo(4, 4);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.12;
  Corr(0, 2) = -0.02;
  Corr(0, 3) = -0.02;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -0.01;
  Corr(1, 3) = -0.06;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.07;
  Corr(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(4, 4);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.08;
  Corr2(0, 2) = 0.;
  Corr2(0, 3) = -0.01;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = -0.02;
  Corr2(1, 3) = -0.04;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.33;
  Corr2(3, 3) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("LHCb_kspp_Au2022", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  //---------------------------------------------------------------------------------------------------------------


  }
// ---------------------------------------------------------

// ---------------------------------------------------------
void MixingModel::Add_other_meas()
{

  vector<dato> CorrData;
  TMatrixDSym Corr, Corr2;

  //-------------------------------------- BESIII Branching ratios -------------------------------------------------------------------------
  // https://arxiv.org/pdf/2503.19542
  // Branching ratios DCS/CF
  meas.insert(pair<string, dato>("BESIII_2503.19542_BrDKpi", dato(0.328, 0.027))); // D -> Kpi
  meas.insert(pair<string, dato>("BESIII_2503.19542_BrDK3pi", dato(0.289, 0.028))); // D -> K3pi
  meas.insert(pair<string, dato>("BESIII_2503.19542_BrDKpipi0", dato(0.212, 0.021))); // D -> Kpipi0



  //-------------------------------------- D -> 4pi & D -> hhpi0 & D -> KKpipi -------------------------------------------------------------------------

  // CP-even fractions Cleo-c
  // https://arxiv.org/pdf/1504.05878
  // Observables 2:
  CorrData.clear();
  CorrData.push_back(dato(0.973, 0.017, 0.)); // F_pipipi0
  CorrData.push_back(dato(0.732, 0.055, 0.)); // F_kkpi0

  // Correlation Matrix (stat):
  Corr.ResizeTo(2, 2);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.;
  Corr(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(2, 2);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.;
  Corr2(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID21", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // CP-even fractions BESIII
  // https://arxiv.org/pdf/2409.07197
  // Observables 2:
  CorrData.clear();
  CorrData.push_back(dato(0.9406, 0.0036, 0.0021)); // F_pipipi0
  CorrData.push_back(dato(0.631, 0.014, 0.011)); // F_kkpi0

  // Correlation Matrix (stat):
  Corr.ResizeTo(2, 2);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.;
  Corr(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(2, 2);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.;
  Corr2(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("2409.07197_F_BESIII", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // CP-even fractions, Cleo-c
  // https://arxiv.org/pdf/1504.05878
  // Observables 2:
  // Observables 1:
  meas.insert(pair<string, dato>("UID20", dato(0.737, 0.028, 0.0))); // Fpipipipi


  // CP-even fractioon 4pi, BESIII
  // https://arxiv.org/pdf/2408.16279
  // Observables 1:
  meas.insert(pair<string, dato>("Fpipipipi_BESIII", dato(0.746, 0.010, 0.004))); // F_pipipipi


  // CP-even fraction, BESIII
  // https://arxiv.org/pdf/2502.12873 supersedes https://arxiv.org/pdf/2212.06489
  // Observables 1:
  meas.insert(pair<string, dato>("FKKpipi_BESIII", dato(0.754, 0.010, 0.008))); // F_KKpipi

  //---------------------------------------------------------------------------------------------------------------


  //-------------------------------------- D -> K3pi & D -> Kpipi0 -------------------------------------------------------------------------

  // D -> K3pi decay parameters, LHCb + BESIII + Cleo-C combo
  // https://arxiv.org/pdf/2103.05988
  // Observables 6:
  CorrData.clear();
  CorrData.push_back(dato(0.445, 0.095, 0.));        // kD_k3pi
  CorrData.push_back(dato(166 * d2r, 23 * d2r, 0.)); // dD_k3pi
  CorrData.push_back(dato(0.79, 0.04, 0.));          // kD_kpipi0
  CorrData.push_back(dato(196 * d2r, 11 * d2r, 0.)); // dD_kpipi0
  CorrData.push_back(dato(0.0550, 0.0007, 0.));      // rD_k3pi
  CorrData.push_back(dato(0.0441, 0.0011, 0.));      // rD_kpipi0

  // Correlation Matrix (stat):
  Corr.ResizeTo(6, 6);
  Corr(0, 0) = 1.;
  Corr(0, 1) = -0.75;
  Corr(0, 2) = 0.;
  Corr(0, 3) = -0.07;
  Corr(0, 4) = 0.52;
  Corr(0, 5) = -0.06;
  Corr(1, 1) = 1.;
  Corr(1, 2) = 0.03;
  Corr(1, 3) = 0.17;
  Corr(1, 4) = -0.42;
  Corr(1, 5) = 0.01;
  Corr(2, 2) = 1.;
  Corr(2, 3) = 0.19;
  Corr(2, 4) = -0.01;
  Corr(2, 5) = -0.01;
  Corr(3, 3) = 1.;
  Corr(3, 4) = -0.02;
  Corr(3, 5) = 0.25;
  Corr(4, 4) = 1.;
  Corr(4, 5) = -0.12;
  Corr(5, 5) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(6, 6);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.;
  Corr2(0, 2) = 0.;
  Corr2(0, 3) = 0.;
  Corr2(0, 4) = 0.;
  Corr2(0, 5) = 0.;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.;
  Corr2(1, 3) = 0.;
  Corr2(1, 4) = 0.;
  Corr2(1, 5) = 0.;
  Corr2(2, 2) = 1.;
  Corr2(2, 3) = 0.;
  Corr2(2, 4) = 0.;
  Corr2(2, 5) = 0.;
  Corr2(3, 3) = 1.;
  Corr2(3, 4) = 0.;
  Corr2(3, 5) = 0.;
  Corr2(4, 4) = 1.;
  Corr2(4, 5) = 0.;
  Corr2(5, 5) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID19", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  //---------------------------------------------------------------------------------------------------------------


  //-------------------------------------- D -> K0sKpi -------------------------------------------------------------------------

  // Br(D0 -> K0sKpi) / Br(D0bar -> K0sKpi), restricted to K^*pm region
  // https://arxiv.org/pdf/1509.06628
  // Observables 1:
  meas.insert(pair<string, dato>("UID22", dato(0.37, 0.003, 0.012))); // Br_kskpi


  // Rates and strong phase Cleo-c
  // Observables 3:
  CorrData.clear();
  CorrData.push_back(dato(0.356, 0.034, 0.007)); // RD_kskpi
  CorrData.push_back(dato(-16.6 * d2r, 18.4 * d2r, 0.)); // dD_kskpi
  CorrData.push_back(dato(0.94, 0.12, 0.));              // kD_kskpi

  // Correlation Matrix (stat):
  Corr.ResizeTo(3, 3);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.;
  Corr(0, 2) = 0.;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -0.5;
  Corr(2, 2) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(3, 3);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.;
  Corr2(0, 2) = 0.;
  Corr2(1, 1) = 1.;
  Corr2(1, 2) = 0.;
  Corr2(2, 2) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID23", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));

  //---------------------------------------------------------------------------------------------------------------

}
// ---------------------------------------------------------

// ---------------------------------------------------------
void MixingModel::Add_old_meas()
{

  vector<dato> CorrData;
  TMatrixDSym Corr, Corr2;


  //-------------------------------------- D -> Kpi -------------------------------------------------------------------------

  // rD_kpi, x'p^2, y'p Babar
  // https://arxiv.org/pdf/hep-ex/0703020
  // Observables 3:
  CorrData.clear();

  CorrData.push_back(dato(0.303e-2, 0.016e-2, 0.010e-2));  // RD_exp_babar_kp
  CorrData.push_back(dato(-0.024e-2, 0.043e-2, 0.030e-2)); // xp2plus_exp_babar_kp
  CorrData.push_back(dato(0.98e-2, 0.64e-2, 0.45e-2));     // ypplus_exp_babar_kp

  Corr.ResizeTo(3, 3);

  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.77;
  Corr(0, 2) = -0.87;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -0.94;
  Corr(2, 2) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpi_babar_plus", CorrelatedGaussianObservables(CorrData, Corr)));


  // AD_kpi, x'm^2, y'm Babar
  // https://arxiv.org/pdf/hep-ex/0703020
  // Observables 3:
  CorrData.clear();

  CorrData.push_back(dato(-2.1e-2, 5.2e-2, 1.5e-2));       // AD_exp_babar_kp
  CorrData.push_back(dato(-0.020e-2, 0.041e-2, 0.029e-2)); // xp2minus_exp_babar_kp
  CorrData.push_back(dato(0.96e-2, 0.61e-2, 0.43e-2));     // ypminus_exp_babar_kp

  Corr.ResizeTo(3, 3);

  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.77;
  Corr(0, 2) = -0.87;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -0.94;
  Corr(2, 2) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpi_babar_minus", CorrelatedGaussianObservables(CorrData, Corr)));


  // rD_kpi, x'p^2, y'p Belle
  // https://arxiv.org/pdf/hep-ex/0601029
  // Observables 3:
  CorrData.clear();

  CorrData.push_back(dato(0.364e-2, .018e-2, 0.));  // RD_exp_belle_kp
  CorrData.push_back(dato(0.032e-2, 0.037e-2, 0.)); // xp2plus_exp_belle_kp
  CorrData.push_back(dato(-0.12e-2, 0.58e-2, 0.));  // ypplus_exp_belle_kp

  Corr.ResizeTo(3, 3);

  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.655;
  Corr(0, 2) = -0.834;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -0.909;
  Corr(2, 2) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpi_belle_plus",
                                                              CorrelatedGaussianObservables(CorrData, Corr)));


  // AD_kpi, x'm^2, y'm Belle
  // https://arxiv.org/pdf/hep-ex/0601029
  // Observables 3:
  CorrData.clear();

  CorrData.push_back(dato(2.3e-2, 4.7e-2));     // AD_exp_belle_kp
  CorrData.push_back(dato(0.006e-2, 0.034e-2)); // xp2minus_exp_belle_kp
  CorrData.push_back(dato(0.20e-2, 0.54e-2));   // ypminus_exp_belle_kp

  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.655;
  Corr(0, 2) = -0.834;
  Corr(1, 1) = 1.;
  Corr(1, 2) = -0.909;
  Corr(2, 2) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpi_belle_minus", CorrelatedGaussianObservables(CorrData, Corr)));


  // rD_kpi, x^2, y, cos(d_kpi), sin(d_kpi) NO CPV
  // https://arxiv.org/pdf/1210.0939
  // Observables 5:
  CorrData.clear();

  CorrData.push_back(dato(0.533e-2, 0.107e-2, .045e-2)); // RD_exp_cleoc
  CorrData.push_back(dato(0.06e-2, 0.23e-2, 0.11e-2));   // xsq_exp_cleoc
  CorrData.push_back(dato(4.2e-2, 2.0e-2, 1.0e-2));      // y_exp_cleoc
  CorrData.push_back(dato(0.84, 0.20, 0.06));            // cd_exp_cleoc
  CorrData.push_back(dato(-0.01, 0.41, .04));            // sd_exp_cleoc

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
  SymmetrizeUpperTriangularMatrix(Corr);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("cleoc", CorrelatedGaussianObservables(CorrData, Corr)));


  //---------------------------------------------------------------------------------------------------------------


  //-------------------------------------- D -> K0spipi -------------------------------------------------------------------------


  // x,y, qop, phi Belle
  // https://arxiv.org/pdf/1404.2412
  // Observables 4:
  CorrData.clear();

  CorrData.push_back(dato(0.52e-2, 0.19e-2, 0.06e-2, 0.08e-2));       // x_exp_belle_kpp
  CorrData.push_back(dato(0.28e-2, 0.15e-2, 0.05e-2, 0.05e-2));       // y_exp_belle_kpp
  CorrData.push_back(dato(0.91, 0.16, 0.05, 0.06));                   // qop_exp_belle_kpp
  CorrData.push_back(dato(-6. * d2r, 11. * d2r, 3. * d2r, 4. * d2r)); // phi_exp_belle_kpp

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
  SymmetrizeUpperTriangularMatrix(Corr);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpp_belle", CorrelatedGaussianObservables(CorrData, Corr)));


  // x, y NO CPV
  // https://arxiv.org/pdf/1510.01664
  // 2 Observables:
  CorrData.clear();
  CorrData.push_back(dato(-0.0086, 0.0053, 0.0017)); // x
  CorrData.push_back(dato(0.0003, 0.0046, 0.0013));  // y

  // Correlation Matrix (stat):
  Corr.ResizeTo(2, 2);
  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.37;
  Corr(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  // Correlation Matrix (syst)
  Corr2.ResizeTo(2, 2);
  Corr2(0, 0) = 1.;
  Corr2(0, 1) = 0.;
  Corr2(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr2);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("UID13", CorrelatedGaussianObservables(CorrData, Corr, Corr2)));


  // x, y NO CPV
  // https://arxiv.org/pdf/1004.5053
  // Observables 2:
  CorrData.clear();

  CorrData.push_back(dato(0.16e-2, 0.23e-2, 0.12e-2, 0.08e-2)); // x_exp_babar
  CorrData.push_back(dato(0.57e-2, 0.20e-2, 0.13e-2, 0.07e-2)); // y_exp_babar

  Corr.ResizeTo(2, 2);

  Corr(0, 0) = 1.;
  Corr(0, 1) = 0.0615;
  Corr(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kppkk", CorrelatedGaussianObservables(CorrData, Corr)));


  //---------------------------------------------------------------------------------------------------------------

  //-------------------------------------- D -> pipipi0 -------------------------------------------------------------------------


  // x, y NO CPV Babar
  // https://arxiv.org/pdf/1604.00857
  // Observables 2:
  CorrData.clear();

  CorrData.push_back(dato(1.5e-2, 1.2e-2, 0.6e-2)); // x_exp_babar
  CorrData.push_back(dato(0.2e-2, 0.9e-2, 0.5e-2)); // y_exp_babar

  Corr.ResizeTo(2, 2);

  Corr(0, 0) = 1.;
  Corr(0, 1) = -0.006;
  Corr(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kppp0_Babar", CorrelatedGaussianObservables(CorrData, Corr)));


  //---------------------------------------------------------------------------------------------------------------

  //-------------------------------------- D -> Kpipi0 -------------------------------------------------------------------------

  // xp_kpp, yp_kpp
  // https://arxiv.org/pdf/0807.4544
  // Observables 2:
  CorrData.clear();
  CorrData.push_back(dato(2.48e-2, 0.59e-2, 0.39e-2)); // xpp_kpp_babar_exp
  CorrData.push_back(dato(-0.07e-2, 0.65e-2, 0.5e-2)); // ypp_kpp_babar_exp

  Corr.ResizeTo(2, 2);

  Corr(0, 0) = 1.;
  Corr(0, 1) = -0.69;
  Corr(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpp_babar_plus", CorrelatedGaussianObservables(CorrData, Corr)));

  // xm_kpp, ym_kpp
  // https://arxiv.org/pdf/0807.4544
  // Observables 2:
  CorrData.clear();

  CorrData.push_back(dato(3.5e-2, 0.78e-2, 0.65e-2)); // xpm_kpp_babar_exp
  CorrData.push_back(dato(-0.92e-2, .78e-2, .41e-2)); // ypm_kpp_babar_exp

  Corr.ResizeTo(2, 2);

  Corr(0, 0) = 1.;
  Corr(0, 1) = -0.66;
  Corr(1, 1) = 1.;
  SymmetrizeUpperTriangularMatrix(Corr);

  corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpp_babar_minus", CorrelatedGaussianObservables(CorrData, Corr)));


  //---------------------------------------------------------------------------------------------------------------

  //-------------------------------------- D -> K3pi -------------------------------------------------------------------------

  // x^2 + y^2 NO CPV
  // https://arxiv.org/pdf/1602.07224
  // Observables 1:
  meas.insert(pair<string, dato>("UID18", dato(0.000048, 0.000018, 0.))); // k3pi // informations about x,y 0.25 * (x^2 + y^2)


  //-------------------------------------- D -> Klnu_l -------------------------------------------------------------------------

  // (x^2 + y^2)/2
  // COMBOS average https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM23/results_mixing.html#semileptonic
  meas.insert(pair<string, dato>("RM", dato(0.013e-2, 0.0269e-2, 0.)));

  //---------------------------------------------------------------------------------------------------------------


  }
// ---------------------------------------------------------

// ---------------------------------------------------------
void MixingModel::DefineParameters()
{
  //------------------------------- Adding model parameters for all the combinations ------------------------------------------------------------------------

  if (comb == 0) //-------- Adding model parameters for the charged B modes only combination--------------------------------
  {


    AddParameter("g", 0., M_PI - 1e-12, "#gamma");
    AddParameter("x12", 0.5e-3, 8.e-3, "x_{12}");
    AddParameter("y12", 3e-3, 9e-3, "y_{12}");


    AddParameter("r_dk", 0.07, 0.13, "r_{DK}");
    AddParameter("r_dpi", 0., 0.012, "r_{D#pi}");
    AddParameter("rD_kpi", 0.055, 0.065, "r^{D}_{K#pi}");
    AddParameter("d_dk", 2.3 - M_PI, 2.3 + M_PI - 1e-12, "#delta_{DK}");
    AddParameter("d_dpi", 5.3 - M_PI, 5.3  + M_PI - 1e-12, "#delta_{D#pi}");
    AddParameter("dD_kpi", 3.3 - M_PI, 3.3 +M_PI - 1e-12, "#delta^{D}_{K#pi}");


    AddParameter("rD_k3pi", 0.048, 0.065, "r^{D}_{K3#pi}");
    AddParameter("dD_k3pi", 2.8 - M_PI, 2.8+M_PI - 1e-12, "#delta^{D}_{K3#pi}");
    AddParameter("kD_k3pi", 0., 1., "#kappa^{D}_{K3#pi}");
    AddParameter("F_pipipipi", 0.5, 1., "F_{#pi#pi#pi#pi}");

    AddParameter("rD_kpipi0", 0.03, 0.06, "r^{D}_{K#pi#pi^{0}}");
    AddParameter("dD_kpipi0", 3.5 -M_PI, 3.5 + M_PI - 1e-12, "#delta^{D}_{K#pi#pi^{0}}");
    AddParameter("kD_kpipi0", 0.4, 1., "#kappa^{D}_{K#pi#pi^{0}}");
    AddParameter("F_pipipi0", 0.6, 1., "F_{#pi#pi#pi^{0}}");
    AddParameter("F_kkpi0", 0.2, 1., "F_{KK#pi^0}");

    AddParameter("rD_kskpi", 0.55, 0.7, "r^{D}_{K^{0}_S K #pi}");
    AddParameter("dD_kskpi", 0.3 - M_PI, 0.3 + M_PI - 1e-12, "#delta^{D}_{K^{0}_S K #pi}");
    AddParameter("kD_kskpi", 0.3, 1., "#kappa_{K^{0}_S K #pi}");
    AddParameter("RBRdkdpi", 0.05, 0.1, "#frac{B^{-} ---> D^0 K^{-}}{B^{-} ---> D^0 #pi^{-}}");

    AddParameter("r_dstk", 0., 0.3, "r_{D^{*}K}");
    AddParameter("d_dstk", -0.7 - M_PI, -0.7 + M_PI - 1e-12, "#delta_{D^{*}K}");
    AddParameter("r_dstpi", 0., 0.06, "r_{D^{*} #pi}");
    AddParameter("d_dstpi",  1.5 - M_PI,  M_PI + 1.5 - 1e-12, "#delta_{D^{*} #pi}");

    AddParameter("r_dkst", 0., 0.3, "r_{DK^{*}}");
    AddParameter("d_dkst", 1. - M_PI, 1. + M_PI - 1e-12, "#delta_{DK^{*}}");
    AddParameter("k_dkst", 0.5, 1., "#kappa_{DK^{*}}");

    AddParameter("r_dkpipi", 0., 0.2, "r_{DK #pi #pi}");
    AddParameter("d_dkpipi", -M_PI - 1.5, M_PI - 1.5 - 1e-12, "#delta_{DK #pi #pi}");
    AddParameter("k_dkpipi", 0., 1., "#kappa_{DK #pi #pi}");
    AddParameter("r_dpipipi", 0., 0.2, "r_{D #pi #pi #pi}");
    AddParameter("d_dpipipi", -M_PI - 1.5, M_PI - 1.5 - 1e-12, "#delta_{D #pi #pi #pi}");
    AddParameter("k_dpipipi", 0., 1., "#kappa_{D #pi #pi #pi}");


    AddParameter("PhiM12", 0.034 -2., 0.034 + 2. - 1e-12, "#phi_{M}");
    AddParameter("PhiG12", 0.042 -2., 0.042 + 2. - 1e-12, "#phi_{#Gamma}");


    AddParameter("adKK", -.1, .1, "a_{d}^{KK}");
    AddParameter("adpipi", -.1, .1, "a_{d}^{#pi#pi}");
    AddParameter("DYKKmDYpipi", -0.1, 0.1, "#Delta Y^{KK} - #Delta Y^{#pi#pi}");
    AddParameter("tavepitaggedOverTauD", 1., 3., "t_{ave}^{#pi tagged}/#tau_{D}");
    AddParameter("tavemutaggedOverTauD", 1.1, 1.3, "t_{ave}^{#mu tagged}/#tau_{D}");
    AddParameter("DeltatmutaggedOverTauD", -0.01, 0.01, "#Delta t^{#mu untagged}/#tau_{D}");
    AddParameter("DeltatpitaggedOverTauD", 0.12, 0.15, "#Delta t^{#pi tagged}/#tau_{D}");
    AddParameter("tKKCDp",6.e-13,8e-13,"t_{KK}^{C_{D^{+}}}");
    AddParameter("tKKCDs",6.5e-13,7.3e-13,"t_{KK}^{C_{D_{s}}}");

    AddParameter("F_kkpipi", 0., 1., "F_{K^+K^-#pi^+#pi^-}");

    // tau Run1 semileptonic tagged
    AddParameter("tauKK_DAcp_Run1_sl",1.05,1.11,"#tau^{R1-sl}_{#Delta Acp}(KK)");
    AddParameter("taupipi_DAcp_Run1_sl",1.,1.11,"#tau^{R1-sl}_{#Delta Acp}(#pi#pi)");
    AddParameter("tauKK_Acp_Run1_sl",1.,1.1,"#tau^{R1-sl}_{Acp}(KK)");

    // tau Run2 pi-tagged
    AddParameter("tauKK_DAcp_Run1_pi",2.,2.3,"#tau^{R1-#pi}_{#Delta Acp}(KK)");
    AddParameter("taupipi_DAcp_Run1_pi",1.9,2.2,"#tau^{R1-#pi}_{#Delta Acp}(#pi#pi)");
    AddParameter("tauKK_Acp_Run1_pi",2.1,2.35,"#tau^{R1-#pi}_{Acp}(#pi#pi)");

    // tau CDF
    AddParameter("tauKK_Acp_CDF",2.3,3.,"#tau^{CDF}_{#Delta Acp}(KK)");
    AddParameter("taupipi_Acp_CDF",2.1,2.8,"#tau^{CDF}_{#Delta Acp}(#pi#pi)");

    SetPriorConstantAll();

  }
  else if (comb == 1)
  { //-------- Adding model parameters for the Neutral-Bd modes --------------------------------

    AddParameter("g", -0.8, -0.8+M_PI, "#gamma");
    AddParameter("x12", 0.5e-3, 8.e-3, "x_{12}");
    AddParameter("y12", 3e-3, 9e-3, "y_{12}");

    AddParameter("rD_kpi", 0.055, 0.065, "r^{D}_{K#pi}");
    AddParameter("dD_kpi", 3.3 - M_PI, 3.3 +M_PI - 1e-12, "#delta^{D}_{K#pi}");

    AddParameter("rD_k3pi", 0.048, 0.065, "r^{D}_{K3#pi}");
    AddParameter("dD_k3pi", 2.8 - M_PI, 2.8+M_PI - 1e-12, "#delta^{D}_{K3#pi}");
    AddParameter("kD_k3pi", 0., 1., "#kappa^{D}_{K3#pi}");
    AddParameter("F_pipipipi", 0.5, 1., "F_{#pi#pi#pi#pi}");

    AddParameter("rD_kpipi0", 0.03, 0.06, "r^{D}_{K#pi#pi^{0}}");
    AddParameter("dD_kpipi0", 3.5 -M_PI, 3.5 + M_PI - 1e-12, "#delta^{D}_{K#pi#pi^{0}}");
    AddParameter("kD_kpipi0", 0.4, 1., "#kappa^{D}_{K#pi#pi^{0}}");
    AddParameter("F_pipipi0", 0.6, 1., "F_{#pi#pi#pi^{0}}");
    AddParameter("F_kkpi0", 0.2, 1., "F_{KK#pi^0}");

    AddParameter("rD_kskpi", 0.55, 0.7, "r^{D}_{K^{0}_S K #pi}");
    AddParameter("dD_kskpi", 0.3 - M_PI, 0.3 + M_PI - 1e-12, "#delta^{D}_{K^{0}_S K #pi}");
    AddParameter("kD_kskpi", 0.3, 1., "#kappa_{K^{0}_S K #pi}");

    AddParameter("r_dkstz", 0., 0.5, "r_{D K^{0*}}");
    AddParameter("d_dkstz", 3.2 - M_PI, 3.2 + M_PI - 1e-12, "#delta_{D K^{0*}}");
    AddParameter("k_dkstz", 0.7, 1., "#kappa_{D K^{0*}}");

    AddParameter("l_dmpi", 0., 0.2, "#lambda_{D^{-} #pi}");
    AddParameter("d_dmpi", 0.5 - M_PI, 0.5 + M_PI - 1e-12, "#delta_{D^{-} #pi}");
    AddParameter("phi_d", 0.75-1., 0.75+1., "#phi_d");

    AddParameter("PhiM12", -2., 2., "#phi_{M}");
    AddParameter("PhiG12", -2., 2., "#phi_{#Gamma}");

    AddParameter("adKK", -.1, .1, "a_{d}^{KK}");
    AddParameter("adpipi", -.1, .1, "a_{d}^{#pi#pi}");
    AddParameter("DYKKmDYpipi", -0.1, 0.1, "#Delta Y^{KK} - #Delta Y^{#pi#pi}");
    AddParameter("tavepitaggedOverTauD", 1., 3., "t_{ave}^{#pi tagged}/#tau_{D}");
    AddParameter("tavemutaggedOverTauD", 1.1, 1.3, "t_{ave}^{#mu tagged}/#tau_{D}");
    AddParameter("DeltatmutaggedOverTauD", -0.01, 0.01, "#Delta t^{#mu untagged}/#tau_{D}");
    AddParameter("DeltatpitaggedOverTauD", 0.12, 0.15, "#Delta t^{#pi tagged}/#tau_{D}");
    AddParameter("tKKCDp",6.e-13,8e-13,"t_{KK}^{C_{D^{+}}}");
    AddParameter("tKKCDs",6.5e-13,7.3e-13,"t_{KK}^{C_{D_{s}}}");

    AddParameter("F_kkpipi", 0., 1., "F_{K^+K^-#pi^+#pi^-}");

    // tau Run1 semileptonic tagged
    AddParameter("tauKK_DAcp_Run1_sl",1.05,1.11,"#tau^{R1-sl}_{#Delta Acp}(KK)");
    AddParameter("taupipi_DAcp_Run1_sl",1.,1.11,"#tau^{R1-sl}_{#Delta Acp}(#pi#pi)");
    AddParameter("tauKK_Acp_Run1_sl",1.,1.1,"#tau^{R1-sl}_{Acp}(KK)");

    // tau Run2 pi-tagged
    AddParameter("tauKK_DAcp_Run1_pi",2.,2.3,"#tau^{R1-#pi}_{#Delta Acp}(KK)");
    AddParameter("taupipi_DAcp_Run1_pi",1.9,2.2,"#tau^{R1-#pi}_{#Delta Acp}(#pi#pi)");
    AddParameter("tauKK_Acp_Run1_pi",2.1,2.35,"#tau^{R1-#pi}_{Acp}(#pi#pi)");

    // tau CDF
    AddParameter("tauKK_Acp_CDF",2.3,3.,"#tau^{CDF}_{#Delta Acp}(KK)");
    AddParameter("taupipi_Acp_CDF",2.1,2.8,"#tau^{CDF}_{#Delta Acp}(#pi#pi)");

    AddParameter("l_dstarmpi", 0., 0.2, "#lambda_{D^{*-} #pi}");
    AddParameter("d_dstarmpi", -0.8 - M_PI, -0.8 + M_PI - 1e-12, "#delta_{D^{*-} #pi}");

    AddParameter("l_dmrho", 0., 0.4, "#lambda_{D^{-} #rho}");
    AddParameter("d_dmrho", -1.51 - M_PI, -1.51 + M_PI - 1e-12, "#delta_{D^{-} #rho}");



    SetPriorConstantAll();

  }
  else if (comb == 2)
  { //-------- Adding model parameters for the Neutral_Bs modes --------------------------------

    AddParameter("g", 0., M_PI - 1e-12, "#gamma");
    AddParameter("x12", 0.5e-3, 8.e-3, "x_{12}");
    AddParameter("y12", 3e-3, 9e-3, "y_{12}");

    AddParameter("rD_kpi", 0.055, 0.065, "r^{D}_{K#pi}");
    AddParameter("dD_kpi", 3.3 - M_PI, 3.3 +M_PI - 1e-12, "#delta^{D}_{K#pi}");

    AddParameter("rD_k3pi", 0.048, 0.065, "r^{D}_{K3#pi}");
    AddParameter("dD_k3pi", 2.8 - M_PI, 2.8+M_PI - 1e-12, "#delta^{D}_{K3#pi}");
    AddParameter("kD_k3pi", 0., 1., "#kappa^{D}_{K3#pi}");
    AddParameter("F_pipipipi", 0.5, 1., "F_{#pi#pi#pi#pi}");

    AddParameter("rD_kpipi0", 0.03, 0.06, "r^{D}_{K#pi#pi^{0}}");
    AddParameter("dD_kpipi0", 3.5 -M_PI, 3.5 + M_PI - 1e-12, "#delta^{D}_{K#pi#pi^{0}}");
    AddParameter("kD_kpipi0", 0.4, 1., "#kappa^{D}_{K#pi#pi^{0}}");
    AddParameter("F_pipipi0", 0.6, 1., "F_{#pi#pi#pi^{0}}");
    AddParameter("F_kkpi0", 0.2, 1., "F_{KK#pi^0}");

    AddParameter("rD_kskpi", 0.55, 0.7, "r^{D}_{K^{0}_S K #pi}");
    AddParameter("dD_kskpi", 0.3 - M_PI, 0.3 + M_PI - 1e-12, "#delta^{D}_{K^{0}_S K #pi}");
    AddParameter("kD_kskpi", 0.3, 1., "#kappa_{K^{0}_S K #pi}");

    AddParameter("l_dsk", 0., 1.2, "#lambda_{D^{-}_{s} K}");
    AddParameter("d_dsk", -M_PI, M_PI- 1e-12, "#delta_{D^{-}_{s} K}");
    AddParameter("phis", -M_PI, M_PI -1e-12, "#phi_s");

    AddParameter("l_dskpipi", 0., 1.2, "#lambda_{D^{-}_{s} K #pi #pi}");
    AddParameter("d_dskpipi", -0.3 -M_PI, -0.3 + M_PI - 1e-12, "#delta_{D_s K #pi #pi}");
    AddParameter("k_dskpipi", 0., 1., "#kappa_{D_s K #pi #pi}");

    AddParameter("PhiM12", 0.034 -2., 0.034 + 2. - 1e-12, "#phi_{M}");
    AddParameter("PhiG12", 0.042 -2., 0.042 + 2. - 1e-12, "#phi_{#Gamma}");

    AddParameter("adKK", -.1, .1, "a_{d}^{KK}");
    AddParameter("adpipi", -.1, .1, "a_{d}^{#pi#pi}");
    AddParameter("DYKKmDYpipi", -0.1, 0.1, "#Delta Y^{KK} - #Delta Y^{#pi#pi}");
    AddParameter("tavepitaggedOverTauD", 1., 3., "t_{ave}^{#pi tagged}/#tau_{D}");
    AddParameter("tavemutaggedOverTauD", 1.1, 1.3, "t_{ave}^{#mu tagged}/#tau_{D}");
    AddParameter("DeltatmutaggedOverTauD", -0.01, 0.01, "#Delta t^{#mu untagged}/#tau_{D}");
    AddParameter("DeltatpitaggedOverTauD", 0.12, 0.15, "#Delta t^{#pi tagged}/#tau_{D}");
    AddParameter("tKKCDp",6.e-13,8e-13,"t_{KK}^{C_{D^{+}}}");
    AddParameter("tKKCDs",6.5e-13,7.3e-13,"t_{KK}^{C_{D_{s}}}");

    AddParameter("F_kkpipi", 0., 1., "F_{K^+K^-#pi^+#pi^-}");

    AddParameter("r_dkstzs", 0., 0.5, "r_{B^{0}_s}^{D K^{0*}}");
    AddParameter("d_dkstzs", 4.45 - 2*M_PI + 1e-12, 4.45, "#delta_{B^{0}_s}^{D K^{0*}}");
    AddParameter("k_dkstzs", 0., 1., "#kappa_{B^{0}_s}^{D K^{0*}}");

    // tau Run1 semileptonic tagged
    AddParameter("tauKK_DAcp_Run1_sl",1.05,1.11,"#tau^{R1-sl}_{#Delta Acp}(KK)");
    AddParameter("taupipi_DAcp_Run1_sl",1.,1.11,"#tau^{R1-sl}_{#Delta Acp}(#pi#pi)");
    AddParameter("tauKK_Acp_Run1_sl",1.,1.1,"#tau^{R1-sl}_{Acp}(KK)");

    // tau Run2 pi-tagged
    AddParameter("tauKK_DAcp_Run1_pi",2.,2.3,"#tau^{R1-#pi}_{#Delta Acp}(KK)");
    AddParameter("taupipi_DAcp_Run1_pi",1.9,2.2,"#tau^{R1-#pi}_{#Delta Acp}(#pi#pi)");
    AddParameter("tauKK_Acp_Run1_pi",2.1,2.35,"#tau^{R1-#pi}_{Acp}(#pi#pi)");

    // tau CDF
    AddParameter("tauKK_Acp_CDF",2.3,3.,"#tau^{CDF}_{#Delta Acp}(KK)");
    AddParameter("taupipi_Acp_CDF",2.1,2.8,"#tau^{CDF}_{#Delta Acp}(#pi#pi)");

    SetPriorConstantAll();
  }
  else if (comb == 3)
  { //-------- Adding model parameters for all the modes --------------------------------

    AddParameter("g", 0., M_PI - 1e-12, "#gamma");
    AddParameter("x12", 0.5e-3, 8.e-3, "x_{12}");
    AddParameter("y12", 3e-3, 9e-3, "y_{12}");

    AddParameter("r_dk", 0.07, 0.13, "r_{DK}");
    AddParameter("r_dpi", 0., 0.012, "r_{D#pi}");
    AddParameter("rD_kpi", 0.055, 0.065, "r^{D}_{K#pi}");
    AddParameter("d_dk", 2.3 - M_PI, 2.3 + M_PI - 1e-12, "#delta_{DK}");
    AddParameter("d_dpi", 5.3 - M_PI, 5.3  + M_PI - 1e-12, "#delta_{D#pi}");
    AddParameter("dD_kpi", 3.3 - M_PI, 3.3 +M_PI - 1e-12, "#delta^{D}_{K#pi}");

    AddParameter("rD_k3pi", 0.048, 0.065, "r^{D}_{K3#pi}");
    AddParameter("dD_k3pi", 2.8 - M_PI, 2.8+M_PI - 1e-12, "#delta^{D}_{K3#pi}");
    AddParameter("kD_k3pi", 0., 1., "#kappa^{D}_{K3#pi}");
    AddParameter("F_pipipipi", 0.5, 1., "F_{#pi#pi#pi#pi}");

    AddParameter("rD_kpipi0", 0.03, 0.06, "r^{D}_{K#pi#pi^{0}}");
    AddParameter("dD_kpipi0", 3.5 -M_PI, 3.5 + M_PI - 1e-12, "#delta^{D}_{K#pi#pi^{0}}");
    AddParameter("kD_kpipi0", 0.4, 1., "#kappa^{D}_{K#pi#pi^{0}}");
    AddParameter("F_pipipi0", 0.6, 1., "F_{#pi#pi#pi^{0}}");
    AddParameter("F_kkpi0", 0.2, 1., "F_{KK#pi^0}");

    AddParameter("rD_kskpi", 0.55, 0.7, "r^{D}_{K^{0}_S K #pi}");
    AddParameter("dD_kskpi", 0.3 - M_PI, 0.3 + M_PI - 1e-12, "#delta^{D}_{K^{0}_S K #pi}");
    AddParameter("kD_kskpi", 0.3, 1., "#kappa_{K^{0}_S K #pi}");
    AddParameter("RBRdkdpi", 0.05, 0.1, "#frac{B^{-} ---> D^0 K^{-}}{B^{-} ---> D^0 #pi^{-}}");

    AddParameter("r_dstk", 0., 0.3, "r_{D^{*}K}");
    AddParameter("d_dstk", -0.7 - M_PI, -0.7 + M_PI - 1e-12, "#delta_{D^{*}K}");
    AddParameter("r_dstpi", 0., 0.06, "r_{D^{*} #pi}");
    AddParameter("d_dstpi", 1.5 - M_PI,  M_PI + 1.5 - 1e-12, "#delta_{D^{*} #pi}");

    AddParameter("r_dkst", 0., 0.3, "r_{DK^{*}}");
    AddParameter("d_dkst", 1. - M_PI, 1. + M_PI - 1e-12, "#delta_{DK^{*}}");
    AddParameter("k_dkst", 0.5, 1., "#kappa_{DK^{*}}");

    AddParameter("r_dkstz", 0., 0.5, "r_{D K^{0*}}");
    AddParameter("d_dkstz", 3.2 - M_PI, 3.2 + M_PI - 1e-12, "#delta_{D K^{0*}}");
    AddParameter("k_dkstz", 0.75, 1., "#kappa_{D K^{0*}}");

    AddParameter("r_dkpipi", 0., 0.2, "r_{DK #pi #pi}");
    AddParameter("d_dkpipi", -M_PI - 1.5, M_PI - 1.5 - 1e-12, "#delta_{DK #pi #pi}");
    AddParameter("k_dkpipi", 0., 1., "#kappa_{DK #pi #pi}");
    AddParameter("r_dpipipi", 0., 0.2, "r_{D #pi #pi #pi}");
    AddParameter("d_dpipipi", -M_PI - 1.5, M_PI - 1.5 - 1e-12, "#delta_{D #pi #pi #pi}");
    AddParameter("k_dpipipi", 0., 1., "#kappa_{D #pi #pi #pi}");

    AddParameter("l_dsk", 0., 1.2, "#lambda_{D^{-}_{s} K}");
    AddParameter("d_dsk", -M_PI, M_PI- 1e-12, "#delta_{D^{-}_{s} K}");
    AddParameter("phis", -M_PI, M_PI - 1e-12, "#phi_s");

    AddParameter("l_dskpipi", 0., 1.2, "#lambda_{D^{-}_{s} K #pi #pi}");
    AddParameter("d_dskpipi", -0.3 -M_PI, -0.3 + M_PI - 1e-12, "#delta_{D_s K #pi #pi}");
    AddParameter("k_dskpipi", 0., 1., "#kappa_{D_s K #pi #pi}");

    AddParameter("l_dmpi", 0., 0.4, "#lambda_{D^{-} #pi}");
    AddParameter("d_dmpi", 0.5 - M_PI, 0.5 + M_PI - 1e-12, "#delta_{D^{-} #pi}");
    AddParameter("phi_d", 0.75-1., 0.75+1. , "#phi_d");

    AddParameter("PhiM12", 0.034 - M_PI, 0.034 + M_PI - 1e-12, "#phi_{M}");
    AddParameter("PhiG12", 0.042 - M_PI, 0.042 + M_PI - 1e-12, "#phi_{#Gamma}");

    AddParameter("adKK", -.1, .1, "a_{d}^{KK}");
    AddParameter("adpipi", -.1, .1, "a_{d}^{#pi#pi}");
    AddParameter("DYKKmDYpipi", -0.1, 0.1, "#Delta Y^{KK} - #Delta Y^{#pi#pi}");
    AddParameter("tavepitaggedOverTauD", 1., 3., "t_{ave}^{#pi tagged}/#tau_{D}");
    AddParameter("tavemutaggedOverTauD", 1.1, 1.3, "t_{ave}^{#mu tagged}/#tau_{D}");
    AddParameter("DeltatmutaggedOverTauD", -0.01, 0.01, "#Delta t^{#mu untagged}/#tau_{D}");
    AddParameter("DeltatpitaggedOverTauD", 0.12, 0.15, "#Delta t^{#pi tagged}/#tau_{D}");
    AddParameter("tKKCDp",6.e-13,8e-13,"t_{KK}^{C_{D^{+}}}");
    AddParameter("tKKCDs",6.5e-13,7.3e-13,"t_{KK}^{C_{D_{s}}}");

    AddParameter("F_kkpipi", 0., 1., "F_{K^+K^-#pi^+#pi^-}");

    AddParameter("r_dkstzs", 0., 0.5, "r_{B^{0}_s}^{D K^{0*}}");
    AddParameter("d_dkstzs", 4.45 - 2*M_PI + 1e-12, 4.45, "#delta_{B^{0}_s}^{D K^{0*}}");
    AddParameter("k_dkstzs", 0., 1., "#kappa_{B^{0}_s}^{D K^{0*}}");

    // tau Run1 semileptonic tagged
    AddParameter("tauKK_DAcp_Run1_sl",1.05,1.11,"#tau^{R1-sl}_{#Delta Acp}(KK)");
    AddParameter("taupipi_DAcp_Run1_sl",1.,1.11,"#tau^{R1-sl}_{#Delta Acp}(#pi#pi)");
    AddParameter("tauKK_Acp_Run1_sl",1.,1.1,"#tau^{R1-sl}_{Acp}(KK)");

    // tau Run2 pi-tagged
    AddParameter("tauKK_DAcp_Run1_pi",2.,2.3,"#tau^{R1-#pi}_{#Delta Acp}(KK)");
    AddParameter("taupipi_DAcp_Run1_pi",1.9,2.2,"#tau^{R1-#pi}_{#Delta Acp}(#pi#pi)");
    AddParameter("tauKK_Acp_Run1_pi",2.1,2.35,"#tau^{R1-#pi}_{Acp}(#pi#pi)");

    // tau CDF
    AddParameter("tauKK_Acp_CDF",2.3,3.,"#tau^{CDF}_{#Delta Acp}(KK)");
    AddParameter("taupipi_Acp_CDF",2.1,2.8,"#tau^{CDF}_{#Delta Acp}(#pi#pi)");

    // Babar beauty
    AddParameter("l_dstarmpi", 0., 0.3, "#lambda_{D^{*-} #pi}");
    AddParameter("d_dstarmpi", -0.8-M_PI,  -0.8+ M_PI - 1e-12, "#delta_{D^{*-} #pi}");

    // Babar beauty
    AddParameter("l_dmrho", 0., 0.4, "#lambda_{D^{-} #rho}");
    AddParameter("d_dmrho",   -1.51-M_PI,  -1.51+M_PI, "#delta_{D^{-} #rho}");


    SetPriorConstantAll();

  }
  else if (comb == 4)
  { //-------- Adding model parameters for the Neutral_Bs modes --------------------------------

    AddParameter("x12", 0.5e-3, 8.e-3, "x_{12}");
    AddParameter("y12", 3e-3, 9e-3, "y_{12}");

    AddParameter("rD_kpi", 0.055, 0.065, "r^{D}_{K#pi}");
    AddParameter("dD_kpi", 3.3 - M_PI, 3.3 +M_PI - 1e-12, "#delta^{D}_{K#pi}");

    AddParameter("rD_k3pi", 0.048, 0.065, "r^{D}_{K3#pi}");
    AddParameter("dD_k3pi", 2.8 - M_PI, 2.8+M_PI - 1e-12, "#delta^{D}_{K3#pi}");
    AddParameter("kD_k3pi", 0., 1., "#kappa^{D}_{K3#pi}");
    AddParameter("F_pipipipi", 0.5, 1., "F_{#pi#pi#pi#pi}");

    AddParameter("rD_kpipi0", 0.03, 0.06, "r^{D}_{K#pi#pi^{0}}");
    AddParameter("dD_kpipi0", 3.5 -M_PI, 3.5 + M_PI - 1e-12, "#delta^{D}_{K#pi#pi^{0}}");
    AddParameter("kD_kpipi0", 0.4, 1., "#kappa^{D}_{K#pi#pi^{0}}");
    AddParameter("F_pipipi0", 0.6, 1., "F_{#pi#pi#pi^{0}}");
    AddParameter("F_kkpi0", 0.2, 1., "F_{KK#pi^0}");

    AddParameter("rD_kskpi", 0.55, 0.7, "r^{D}_{K^{0}_S K #pi}");
    AddParameter("dD_kskpi", 0.3 - M_PI, 0.3 + M_PI - 1e-12, "#delta^{D}_{K^{0}_S K #pi}");
    AddParameter("kD_kskpi", 0.3, 1., "#kappa_{K^{0}_S K #pi}");

    AddParameter("PhiM12", 0.034 - M_PI, 0.034 + M_PI - 1e-12, "#phi_{M}");
    AddParameter("PhiG12", 0.042 - M_PI, 0.042 + M_PI - 1e-12, "#phi_{#Gamma}");

    AddParameter("adKK", -.1, .1, "a_{d}^{KK}");
    AddParameter("adpipi", -.1, .1, "a_{d}^{#pi#pi}");
    AddParameter("DYKKmDYpipi", -0.1, 0.1, "#Delta Y^{KK} - #Delta Y^{#pi#pi}");
    AddParameter("tavepitaggedOverTauD", 1., 3., "t_{ave}^{#pi tagged}/#tau_{D}");
    AddParameter("tavemutaggedOverTauD", 1.1, 1.3, "t_{ave}^{#mu tagged}/#tau_{D}");
    AddParameter("DeltatmutaggedOverTauD", -0.01, 0.01, "#Delta t^{#mu untagged}/#tau_{D}");
    AddParameter("DeltatpitaggedOverTauD", 0.12, 0.15, "#Delta t^{#pi tagged}/#tau_{D}");
    AddParameter("tKKCDp",6.e-13,8e-13,"t_{KK}^{C_{D^{+}}}");
    AddParameter("tKKCDs",6.5e-13,7.3e-13,"t_{KK}^{C_{D_{s}}}");

    AddParameter("F_kkpipi", 0., 1., "F_{K^+K^-#pi^+#pi^-}");

    // tau Run1 semileptonic tagged
    AddParameter("tauKK_DAcp_Run1_sl",1.05,1.11,"#tau^{R1-sl}_{#Delta Acp}(KK)");
    AddParameter("taupipi_DAcp_Run1_sl",1.,1.11,"#tau^{R1-sl}_{#Delta Acp}(#pi#pi)");
    AddParameter("tauKK_Acp_Run1_sl",1.,1.1,"#tau^{R1-sl}_{Acp}(KK)");

    // tau Run2 pi-tagged
    AddParameter("tauKK_DAcp_Run1_pi",2.,2.3,"#tau^{R1-#pi}_{#Delta Acp}(KK)");
    AddParameter("taupipi_DAcp_Run1_pi",1.9,2.2,"#tau^{R1-#pi}_{#Delta Acp}(#pi#pi)");
    AddParameter("tauKK_Acp_Run1_pi",2.1,2.35,"#tau^{R1-#pi}_{Acp}(#pi#pi)");

    // tau CDF
    AddParameter("tauKK_Acp_CDF",2.3,3.,"#tau^{CDF}_{#Delta Acp}(KK)");
    AddParameter("taupipi_Acp_CDF",2.1,2.8,"#tau^{CDF}_{#Delta Acp}(#pi#pi)");


    SetPriorConstantAll();
  }

}
// ---------------------------------------------------------

// ---------------------------------------------------------
void MixingModel::DefineHistograms()
{

  //-------------------------------------- Generating 1D and 2D histograms ---------------------------------------
  int j = 0;
  for (int i = 0; i < nVarab.size(); i++)
  { // looping over all the names of the variables
    histos.createH1D(nVarab[i], 200, 1., -1.);
    for (int j = i + 1; j < nVarab.size(); j++)
    {
      histos.createH2D(nVarab[i], nVarab[j], 200, 1., -1., 200, 1., -1.);
    }
  }
}
// ---------------------------------------------------------

double MixingModel::Acp(double rB, double delta_B, double kB, double F_D, double alpha)
{
  return (2 * rB * kB * sin(g) * sin(delta_B) * ((2 * F_D - 1) - alpha * y12)) /
         ((1 + rB * rB) * (1 - alpha * y12 * (2 * F_D - 1)) + 2 * rB * kB * cos(g) * cos(delta_B) * ((2 * F_D - 1) - alpha * y12));
}

// ---------------------------------------------------------

double MixingModel::Afav(double rB, double rD, double delta_B, double delta_D, double kB, double kD, double alpha)
{
  return (2 * rD * rB * kD * kB * sin(g) * sin(delta_B - delta_D) - alpha * rB * kB * sin(g) * (x12 * cos(delta_B) * (1 - rD * rD) + y12 * sin(delta_B) * (1 + rD * rD))) /
         (1 + rD * rD * rB * rB + 2 * rB * rD * kB * kD * cos(g) * cos(delta_B - delta_D) - alpha * y12 * (rD * kD * (1 + rB * rB) * cos(delta_D) + rB * kB * (1 + rD * rD) * cos(g) * cos(delta_B)) + alpha * x12 * (rB * kB * (1 - rD * rD) * cos(g) * sin(delta_B) - rD * kD * (1 - rB * rB) * sin(delta_D)));
}

// ---------------------------------------------------------

double MixingModel::Rcp_h(double rBCP, double delta_BCP, double rBCF, double rD, double delta_BCF, double delta_D, double kB, double kD, double F_D, double alpha)
{
  return ((1 + rBCP * rBCP) * (1 - (2 * F_D - 1) * alpha * y12) + 2 * kB * rBCP * cos(g) * cos(delta_BCP) * ((2 * F_D - 1) - alpha * y12)) /
         (1 + rD * rD * rBCF * rBCF + 2 * rBCF * rD * kB * kD * cos(g) * cos(delta_BCF - delta_D) - alpha * y12 * (rD * kD * (1 + rBCF * rBCF) * cos(delta_D) + kB * rBCF * (1 + rD * rD) * cos(g) * cos(delta_BCF)) + alpha * x12 * (kB * rBCF * (1 - rD * rD) * cos(g) * sin(delta_BCF) - rD * kD * (1 - rBCF * rBCF) * sin(delta_D)));
}

// ---------------------------------------------------------

double MixingModel::Rm(double rB, double rD, double delta_B, double delta_D, double kB, double kD, double alpha)
{
  return (rD * rD + rB * rB + 2 * kB * kD * rB * rD * cos(delta_B - g + delta_D) - alpha * y12 * (rD * kD * (1 + rB * rB) * cos(delta_D) + rB * kB * (1 + rD * rD) * cos(delta_B - g)) - alpha * x12 * (rD * kD * (1 - rB * rB) * sin(-delta_D) + rB * kB * (1 - rD * rD) * sin(delta_B - g))) /
         (1 + rD * rD * rB * rB + 2 * rD * rB * kD * kB * cos(delta_B - g - delta_D) - alpha * y12 * (rD * kD * (1 + rB * rB) * cos(delta_D) + rB * kB * (1 + rD * rD) * cos(delta_B - g)) + alpha * x12 * (rD * kD * (1 - rB * rB) * sin(-delta_D) + rB * kB * (1 - rD * rD) * sin(delta_B - g)));
}

// ---------------------------------------------------------

double MixingModel::Rp(double rB, double rD, double delta_B, double delta_D, double kB, double kD, double alpha)
{
  return (rD * rD + rB * rB + 2 * kB * kD * rB * rD * cos(delta_B + g + delta_D) - alpha * y12 * (rD * kD * (1 + rB * rB) * cos(delta_D) + rB * kB * (1 + rD * rD) * cos(delta_B + g)) - alpha * x12 * (rD * kD * (1 - rB * rB) * sin(-delta_D) + rB * kB * (1 - rD * rD) * sin(delta_B + g))) /
         (1 + rD * rD * rB * rB + 2 * rD * rB * kD * kB * cos(delta_B + g - delta_D) - alpha * y12 * (rD * kD * (1 + rB * rB) * cos(delta_D) + rB * kB * (1 + rD * rD) * cos(delta_B + g)) + alpha * x12 * (rD * kD * (1 - rB * rB) * sin(-delta_D) + rB * kB * (1 - rD * rD) * sin(delta_B + g)));
}

// ---------------------------------------------------------

double MixingModel::Asup(double rB, double rD, double delta_B, double delta_D, double kB, double kD, double alpha)
{
  return (2 * rB * rD * kB * kD * sin(g) * sin(delta_B + delta_D) - alpha * rB * kB * sin(g) * (y12 * (1 + rD * rD) * sin(delta_B) - x12 * (1 - rD * rD) * cos(delta_B))) /
         (rB * rB + rD * rD + 2 * rB * rD * kB * kD * cos(g) * cos(delta_B + delta_D) - alpha * y12 * (rB * kB * (1 + rD * rD) * cos(g) * cos(delta_B) + rD * kD * (1 + rB * rB) * cos(delta_D)) - alpha * x12 * (rB * kB * (1 - rD * rD) * cos(g) * sin(delta_B) + rD * kD * (rB * rB - 1) * sin(delta_D)));
}

// ---------------------------------------------------------

double MixingModel::Rads(double rB, double rD, double delta_B, double delta_D, double kB, double kD, double alpha)
{
  return (rD * rD + rB * rB + 2 * rB * rD * kB * kD * cos(g) * cos(delta_B + delta_D) - alpha * y12 * (rD * kD * (1 + rB * rB) * cos(delta_D) + kB * rB * (1 + rD * rD) * cos(g) * cos(delta_B)) - alpha * x12 * (rD * kD * (1 - rB * rB) * sin(-delta_D) + rB * kB * (1 - rD * rD) * cos(g) * sin(delta_B))) /
         (1 + rD * rD * rB * rB + 2 * rB * rD * kB * kD * cos(g) * cos(delta_B - delta_D) - alpha * y12 * (rD * kD * (1 + rB * rB) * cos(delta_D) + rB * kB * (1 + rD * rD) * cos(g) * cos(delta_B)) + alpha * x12 * (rD * kD * (1 - rB * rB) * sin(-delta_D) + rB * kB * (1 - rD * rD) * cos(g) * sin(delta_B)));
}

// ---------------------------------------------------------

double MixingModel::Rfav(double rB1, double rB2, double rD, double delta_B1, double delta_B2, double delta_D, double BR, double kD, double alpha)
{
  return BR * ((1 + rB1 * rB1 * rD * rD + 2 * kD * rD * rB1 * cos(g) * cos(delta_B1 - delta_D) - alpha * y12 * (rD * kD * (1 + rB1 * rB1) * cos(delta_D) + rB1 * (1 + rD * rD) * cos(g) * cos(delta_B1)) + alpha * x12 * (rD * kD * (1 - rB1 * rB1) * sin(-delta_D) + rB1 * (1 - rD * rD) * cos(g) * sin(delta_B1))) /
               (1 + rB2 * rB2 * rD * rD + 2 * kD * rD * rB2 * cos(g) * cos(delta_B2 - delta_D) - alpha * y12 * (rD * kD * (1 + rB2 * rB2) * cos(delta_D) + rB2 * (1 + rD * rD) * cos(g) * cos(delta_B2)) + alpha * x12 * (rD * kD * (1 - rB2 * rB2) * sin(-delta_D) + rB2 * (1 - rD * rD) * cos(g) * sin(delta_B2))));
}

// ---------------------------------------------------------

double MixingModel::Rsup(double rB1, double rB2, double rD, double delta_B1, double delta_B2, double delta_D, double BR, double kD, double alpha)
{
  return BR * ((rD * rD + rB1 * rB1 + 2 * rD * kD * rB1 * cos(g) * cos(delta_B1 + delta_D) - alpha * y12 * (rD * kD * (1 + rB1 * rB1) * cos(delta_D) + rB1 * (1 + rD * rD) * cos(g) * cos(delta_B1)) - alpha * x12 * (rD * kD * (1 - rB1 * rB1) * sin(-delta_D) + rB1 * (1 - rD * rD) * cos(g) * sin(delta_B1))) /
               (rD * rD + rB2 * rB2 + 2 * rD * kD * rB2 * cos(g) * cos(delta_B2 + delta_D) - alpha * y12 * (rD * kD * (1 + rB2 * rB2) * cos(delta_D) + rB2 * (1 + rD * rD) * cos(g) * cos(delta_B2)) - alpha * x12 * (rD * kD * (1 - rB2 * rB2) * sin(-delta_D) + rB2 * (1 - rD * rD) * cos(g) * sin(delta_B2))));
}

// ---------------------------------------------------------

double MixingModel::y_plus(double delta_D)
{
  return qop * (sin(phi) * (x * cos(delta_D) + y * sin(delta_D)) - cos(phi) * (-x * sin(delta_D) + y * cos(delta_D)));
}

// ---------------------------------------------------------

double MixingModel::y_minus(double delta_D)
{
  return (1. / qop) * (-sin(phi) * (x * cos(delta_D) + y * sin(delta_D)) - cos(phi) * (-x * sin(delta_D) + y * cos(delta_D)));
}

// ---------------------------------------------------------

double MixingModel::x_plus(double delta_D)
{
  return (qop) * (-cos(phi) * (x * cos(delta_D) + y * sin(delta_D)) - sin(phi) * (-x * sin(delta_D) + y * cos(delta_D)));
}

// ---------------------------------------------------------

double MixingModel::x_minus(double delta_D)
{
  return (-1. / qop) * (cos(phi) * (x * cos(delta_D) + y * sin(delta_D)) - sin(phi) * (-x * sin(delta_D) + y * cos(delta_D)));
}

// ---------------------------------------------------------
double MixingModel::Calculate_ChargedB_observables()
{

  double ll1;
  ll1 = 0.;

  //-------------------------------------------------  Bpm -> Dhpm  -------------------------------------------------------------------------

  // GLW: D -> KK, pipi; ADS: D -> Kpi
  // https://arxiv.org/pdf/2012.09903.pdf
  // Observables 8:
  acp_dk_uid0 = Acp(r_dk, d_dk, 1., 1., 2 * 0.523);
  acp_dpi_uid0 = Acp(r_dpi, d_dpi, 1., 1., 2 * 0.523);
  afav_dk_uid0 = Afav(r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 2 * 0.523);
  rcp_uid0 = Rcp_h(r_dk, d_dk, r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 1., 2 * 0.523) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 1., 2 * 0.523);
  rm_dk_uid0 = Rm(r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 2 * 0.523);
  rm_dpi_uid0 = Rm(r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 2 * 0.523);
  rp_dk_uid0 = Rp(r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 2 * 0.523);
  rp_dpi_uid0 = Rp(r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 2 * 0.523);


  // GLW: D -> KKpipi, D -> 4pi
  // https://arxiv.org/abs/2301.10328
  // 6 Observables
  acp_dk_kkpipi_230110328 = Acp(r_dk, d_dk, 1., F_kkpipi, 1.);
  acp_dpi_kkpipi_230110328 = Acp(r_dpi, d_dpi, 1., F_kkpipi, 1.);
  acp_dk_pipipipi_230110328 = Acp(r_dk, d_dk, 1., F_pipipipi, 1.);
  acp_dpi_pipipipi_230110328 = Acp(r_dpi, d_dpi, 1., F_pipipipi, 1.);
  rcp_kpi_kkpipi_230110328 = Rcp_h(r_dk, d_dk, r_dk, rD_k3pi, d_dk, dD_k3pi, 1., kD_k3pi, F_kkpipi, 1.) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_k3pi, d_dpi, dD_k3pi, 1., kD_k3pi, F_kkpipi, 1.);
  rcp_kpi_pipipipi_230110328 = Rcp_h(r_dk, d_dk, r_dk, rD_k3pi, d_dk, dD_k3pi, 1., kD_k3pi, F_pipipipi, 1.) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_k3pi, d_dpi, dD_k3pi, 1., kD_k3pi, F_pipipipi, 1.);


  // GLW: D -> pipipi0, KKpi0; D -> Kpipi0
  //https://arxiv.org/pdf/2112.10617
  // Observables 11:
  rcp_kkpi0_211210617 = Rcp_h(r_dk, d_dk, r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, F_kkpi0, 2 * 0.5) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_kpipi0, d_dpi, dD_kpipi0, 1., kD_kpipi0, F_kkpi0, 2 * 0.5);
  rcp_pipipi0_211210617 = Rcp_h(r_dk, d_dk, r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, F_pipipi0, 2 * 0.5) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_kpipi0, d_dpi, dD_kpipi0, 1., kD_kpipi0, F_pipipi0, 2 * 0.5);
  afav_dk_kpipi0_211210617 = Afav(r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, 2 * 0.5);
  acp_dk_kkpi0_211210617 = Acp(r_dk, d_dk, 1., F_kkpi0, 2 * 0.5);
  acp_dk_pipipi0_211210617 = Acp(r_dk, d_dk, 1., F_pipipi0, 2 * 0.5);
  acp_dpi_kkpi0_211210617 = Acp(r_dpi, d_dpi, 1., F_kkpi0, 2 * 0.5);
  acp_dpi_pipipi0_211210617 = Acp(r_dpi, d_dpi, 1., F_pipipi0, 2 * 0.5);
  rp_dk_211210617 = Rp(r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, 2 * 0.5);
  rm_dk_211210617 = Rm(r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, 2 * 0.5);
  rp_dpi_211210617 = Rp(r_dpi, rD_kpipi0, d_dpi, dD_kpipi0, 1., kD_kpipi0, 2 * 0.5);
  rm_dpi_211210617 = Rm(r_dpi, rD_kpipi0, d_dpi, dD_kpipi0, 1., kD_kpipi0, 2 * 0.5);


  // ADS: D -> K0sKpi
  // https://arxiv.org/pdf/2002.08858
  // Observables 7:
  afav_dpi_kskpi_uid4 = Afav(r_dpi, rD_kskpi, d_dpi, dD_kskpi, 1., kD_kskpi, 1.);
  asup_dpi_kskpi_uid4 = Asup(r_dpi, rD_kskpi, d_dpi, dD_kskpi, 1., kD_kskpi, 1.);
  afav_dk_kskpi_uid4 = Afav(r_dk, rD_kskpi, d_dk, dD_kskpi, 1., kD_kskpi, 1.);
  asup_dk_kskpi_uid4 = Asup(r_dk, rD_kskpi, d_dk, dD_kskpi, 1., kD_kskpi, 1.);
  rfavsup_dpi_kskpi_uid4 = 1. / Rads(r_dpi, rD_kskpi, d_dpi, dD_kskpi, 1., kD_kskpi, 1.);
  rfav_dkdpi_kskpi_uid4 = Rfav(r_dk, r_dpi, rD_kskpi, d_dk, d_dpi, dD_kskpi, RBRdkdpi, kD_kskpi, 1.);
  rsup_dkdpi_kskpi_uid4 = Rsup(r_dk, r_dpi, rD_kskpi, d_dk, d_dpi, dD_kskpi, RBRdkdpi, kD_kskpi, 1.);


  xm_dk_uid3 = r_dk * cos(d_dk - g);
  ym_dk_uid3 = r_dk * sin(d_dk - g);
  xp_dk_uid3 = r_dk * cos(d_dk + g);
  yp_dk_uid3 = r_dk * sin(d_dk + g);
  xi_x_dpi_uid3 = (r_dpi / r_dk) * cos(d_dpi - d_dk);
  xi_y_dpi_uid3 = (r_dpi / r_dk) * sin(d_dpi - d_dk);


  // GLW: D -> KK, D -> pipi; ADS: D -> Kpi
  // https://arxiv.org/abs/2012.09903
  // Observables 18:
  acp_dstk_dg_uid5 = Acp(r_dstk, d_dstk + M_PI, 1., 1., 2 * 0.523);
  acp_dstk_dp_uid5 = Acp(r_dstk, d_dstk, 1., 1., 2 * 0.523);
  afav_dstk_dg_uid5 = Afav(r_dstk, rD_kpi, d_dstk + M_PI, dD_kpi, 1., 1., 2 * 0.523);
  afav_dstk_dp_uid5 = Afav(r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 2 * 0.523);
  rcp_dg_uid5 = Rcp_h(r_dstk, d_dstk + M_PI, r_dstk, rD_kpi, d_dstk + M_PI, dD_kpi, 1., 1., 1., 2 * 0.523) / Rcp_h(r_dstpi, d_dstpi + M_PI, r_dstpi, rD_kpi, d_dstpi + M_PI, dD_kpi, 1., 1., 1., 2 * 0.523);
  rcp_dp_uid5 = Rcp_h(r_dstk, d_dstk, r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 1., 2 * 0.523) / Rcp_h(r_dstpi, d_dstpi, r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 1., 2 * 0.523);
  rm_dstk_dg_uid5 = Rm(r_dstk, rD_kpi, d_dstk + M_PI, dD_kpi, 1., 1., 2 * 0.523);
  rm_dstk_dp_uid5 = Rm(r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 2 * 0.523);
  rp_dstk_dg_uid5 = Rp(r_dstk, rD_kpi, d_dstk + M_PI, dD_kpi, 1., 1., 2 * 0.523);
  rp_dstk_dp_uid5 = Rp(r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 2 * 0.523);
  acp_dstpi_dg_uid5 = Acp(r_dstpi, d_dstpi + M_PI, 1., 1., 2 * 0.523);
  acp_dstpi_dp_uid5 = Acp(r_dstpi, d_dstpi, 1., 1., 2 * 0.523);
  rm_dstpi_dg_uid5 = Rm(r_dstpi, rD_kpi, d_dstpi + M_PI, dD_kpi, 1., 1., 2 * 0.523);
  rm_dstpi_dp_uid5 = Rm(r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 2 * 0.523);
  rp_dstpi_dg_uid5 = Rp(r_dstpi, rD_kpi, d_dstpi + M_PI, dD_kpi, 1., 1., 2 * 0.523);
  rp_dstpi_dp_uid5 = Rp(r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 2 * 0.523);
  afav_dstpi_dg_uid5 = Afav(r_dstpi, rD_kpi, d_dstpi + M_PI, dD_kpi, 1., 1., 2 * 0.523);
  afav_dstpi_dp_uid5 = Afav(r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 2 * 0.523);


  // GLW: D -> KK, D -> pipi, D -> 4pi; ADS: D -> Kpi, D -> K3pi
  // LHCb-PAPER-2024-023
  // Observables 12:
  afav_dkst_kpi = Afav(r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 2 * 0.6);
  acp_dkst_kk = Acp(r_dkst, d_dkst, k_dkst, 1., 2 * 0.6);
  acp_dkst_pipi = Acp(r_dkst, d_dkst, k_dkst, 1., 2 * 0.6);
  asup_dkst_kpi = Asup(r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 2 * 0.6);
  rcp_dkst_kk = Rcp_h(r_dkst, d_dkst, r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 1., 2 * 0.6);
  rcp_dkst_pipi = Rcp_h(r_dkst, d_dkst, r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 1., 2 * 0.6);
  rsup_dkst_kpi = Rads(r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 2 * 0.6);
  afav_dkst_k3pi = Afav(r_dkst, rD_k3pi, d_dkst, dD_k3pi, k_dkst, kD_k3pi, 2 * 0.6);
  acp_dkst_pipipipi = Acp(r_dkst, d_dkst, k_dkst, F_pipipipi, 2 * 0.6);
  asup_dkst_k3pi = Asup(r_dkst, rD_k3pi, d_dkst, dD_k3pi, k_dkst, kD_k3pi, 2 * 0.6);
  rcp_dkst_pipipipi = Rcp_h(r_dkst, d_dkst, r_dkst, rD_k3pi, d_dkst, dD_k3pi, k_dkst, kD_k3pi, F_pipipipi, 2 * 0.6);
  rsup_dkst_k3pi = Rads(r_dkst, rD_k3pi, d_dkst, dD_k3pi, k_dkst, kD_k3pi, 2 * 0.6);


  // GLW: D -> KK, D -> pipi; ADS: D -> Kpi
  // https://arxiv.org/pdf/1505.07044
  // Observables 11:
  rcp_dkpipi_uid9 = Rcp_h(r_dkpipi, d_dkpipi, r_dkpipi, rD_kpi, d_dkpipi, dD_kpi, k_dkpipi, 1., 1., 2 * 0.6) / Rcp_h(r_dpipipi, d_dpipipi, r_dpipipi, rD_kpi, d_dpipipi, dD_kpi, k_dpipipi, 1., 1., 2 * 0.6);
  afav_dkpipi_kpi_uid9 = Afav(r_dkpipi, rD_kpi, d_dkpipi, dD_kpi, k_dkpipi, 1., 2 * 0.6);
  afav_dpipipi_kpi_uid9 = Afav(r_dpipipi, rD_kpi, d_dpipipi, dD_kpi, k_dpipipi, 1., 2 * 0.6);
  acp_dkpipi_kk_uid9 = Acp(r_dkpipi, d_dkpipi, k_dkpipi, 1., 2 * 0.6);
  acp_dkpipi_pipi_uid9 = Acp(r_dkpipi, d_dkpipi, k_dkpipi, 1., 2 * 0.6);
  acp_dpipipi_kk_uid9 = Acp(r_dpipipi, d_dpipipi, k_dpipipi, 1., 2 * 0.6);
  acp_dpipipi_pipi_uid9 = Acp(r_dpipipi, d_dpipipi, k_dpipipi, 1., 2 * 0.6);
  rp_dkpipi_uid9 = Rp(r_dkpipi, rD_kpi, d_dkpipi, dD_kpi, k_dkpipi, 1., 2 * 0.6);
  rm_dkpipi_uid9 = Rm(r_dkpipi, rD_kpi, d_dkpipi, dD_kpi, k_dkpipi, 1., 2 * 0.6);
  rp_dpipipi_uid9 = Rp(r_dpipipi, rD_kpi, d_dpipipi, dD_kpi, k_dpipipi, 1., 2 * 0.6);
  rm_dpipipi_uid9 = Rm(r_dpipipi, rD_kpi, d_dpipipi, dD_kpi, k_dpipipi, 1., 2 * 0.6);


  // Coherence factor kappakstpm
  // https://arxiv.org/pdf/1709.05855
  k_dkst_uid24 = k_dkst;

  //------------------------------------------------------ Calculating the contribution to the LogLikelihood ----------------------------------------------------

  TVectorD corr(8);

  //-------------------------------------------------  Babar measurements  -------------------------------------------------------------------------

  // B -> DK, D -> KK, D -> pipi normalized to D -> Kpi
  // https://journals.aps.org/prd/pdf/10.1103/PhysRevD.82.072004
  // 4 Observables :
  corr.ResizeTo(4);
  corr(0) = Acp(r_dk, d_dk, 1., 1., 2*0.5);
  corr(1) = Acp(r_dk, d_dk, 1., 0., 2*0.5);
  corr(2) =  Rcp_h(r_dk, d_dk, r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 1., 2 * 0.5) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 1., 2 * 0.5);
  corr(3) =  Rcp_h(r_dk, d_dk, r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 0., 2 * 0.5) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 0., 2 * 0.5);
  ll1 += corrmeas.at("Babar_PRD82_072004").logweight(corr);


  //  https://arxiv.org/pdf/0807.2408
  // B -> DstarK
  // D -> KK, pipi + fcp-
  // 4 Observables :
  ll1 += meas.at("Babar_0807.2408_ACPp").logweight(Acp(r_dstk, d_dstk, 1., 1., 2 * 0.5));
  ll1 += meas.at("Babar_0807.2408_ACPm").logweight(Acp(r_dstk, d_dstk, 1., 0., 2 * 0.5));
  ll1 += meas.at("Babar_0807.2408_RCPp").logweight(Rcp_h(r_dstk, d_dstk, r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 1., 2 * 0.5) / Rcp_h(r_dstpi, d_dstpi, r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 1., 2 * 0.5));
  ll1 += meas.at("Babar_0807.2408_RCPm").logweight(Rcp_h(r_dstk, d_dstk, r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 0., 2 * 0.5) / Rcp_h(r_dstpi, d_dstpi, r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 0., 2 * 0.5));


  // https://journals.aps.org/prd/pdf/10.1103/PhysRevD.80.092001
  // B -> DKstar
  // D -> KK, pipi, fcp-
  ll1 += meas.at("Babar_PRD80_092001_ACPp").logweight(Acp(r_dkst, d_dkst, k_dkst, 1., 2 * 0.5));
  ll1 += meas.at("Babar_PRD80_092001_ACPm").logweight(Acp(r_dkst, d_dkst, k_dkst, 0., 2 * 0.5));
  ll1 += meas.at("Babar_PRD80_092001_RCPp").logweight(Rcp_h(r_dkst, d_dkst, r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 1., 2 * 0.5)  );
  ll1 += meas.at("Babar_PRD80_092001_RCPm").logweight(Rcp_h(r_dkst, d_dkst, r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 0., 2 * 0.5)  );
  //https://arxiv.org/pdf/0909.3981
  // B -> DK
  // D -> Kpi
  ll1 += meas.at("Babar_0909.3981_RADS").logweight(Rads(r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 2 * 0.5));
  ll1 += meas.at("Babar_0909.3981_Asup").logweight(Asup(r_dkst, rD_kpi, d_dkst, dD_kpi, k_dkst, 1., 2 * 0.5));


  // https://arxiv.org/pdf/hep-ex/0703037
  // B -> DK
  // GLW D -> pi+pi-pi0
  ll1 += meas.at("Babar_0703037").logweight(Acp(r_dk, d_dk, 1., F_pipipi0, 2 * 0.5));
  corr.ResizeTo(4);
  corr(0) = sqrt( ( r_dk * cos(d_dk + g) - (2 * F_pipipi0 -1 ) ) * ( r_dk * cos(d_dk + g) - (2 * F_pipipi0 -1 ) ) + ( r_dk * sin(d_dk + g) )*( r_dk * sin(d_dk + g) ) );
  double thetap_babar = atan2( r_dk * sin( d_dk + g ), ( r_dk * cos(d_dk + g)  - (2*F_pipipi0 -1) )  );
  if( thetap_babar < 0){
    thetap_babar+= 2*M_PI; // in [0, 2 pi]
  }
  corr(1) = thetap_babar;
  corr(2) = sqrt( ( r_dk * cos(d_dk - g) - (2 * F_pipipi0 -1 ) ) * ( r_dk * cos(d_dk - g) - (2 * F_pipipi0 -1 ) ) + ( r_dk * sin(d_dk - g) )*( r_dk * sin(d_dk - g) ) );
  double thetam_babar = atan2( r_dk * sin( d_dk - g ) , ( r_dk * cos(d_dk - g)  - (2*F_pipipi0 -1)  ) );
  if( thetam_babar < 0){
    thetam_babar+= 2*M_PI;
  }
  corr(3) = thetam_babar;
  ll1 += corrmeas.at("Babar_0703037_rhotheta").logweight(corr);


  // https://arxiv.org/pdf/1006.4241
  // B -> DK
  // D -> Kpi
  ll1 += meas.at("Babar_1006.4241_RDK").logweight( 0.5 * ( Rp(r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 2 * 0.5) + Rm(r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 2 * 0.5) ) );
  ll1 += meas.at("Babar_1006.4241_ADK").logweight( ( - Rp(r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 2 * 0.5) + Rm(r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 2 * 0.5) ) / ( Rp(r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 2 * 0.5) + Rm(r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 2 * 0.5) ) );
  // B -> [Dpi0]_Dstar K
  // D -> Kpi
  ll1 += meas.at("Babar_1006.4241_RDstarKpi0").logweight( 0.5 * (  Rm(r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 2 * 0.5) +  Rp(r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 2 * 0.5)  )  );
  ll1 += meas.at("Babar_1006.4241_ADstarKpi0").logweight( (  Rm(r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 2 * 0.5) -  Rp(r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 2 * 0.5)  ) / (  Rm(r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 2 * 0.5) +  Rp(r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 2 * 0.5)  ) );
  // B -> [Dg]_Dstar K
  // D -> Kpi
  ll1 += meas.at("Babar_1006.4241_RDstarKg").logweight( 0.5 * ( Rm(r_dstk, rD_kpi, d_dstk + M_PI, dD_kpi, 1., 1., 2 * 0.5) + Rp(r_dstk, rD_kpi, d_dstk + M_PI, dD_kpi, 1., 1., 2 * 0.5) ) );
  ll1 += meas.at("Babar_1006.4241_ADstarKg").logweight( ( Rm(r_dstk, rD_kpi, d_dstk + M_PI, dD_kpi, 1., 1., 2 * 0.5) - Rp(r_dstk, rD_kpi, d_dstk + M_PI, dD_kpi, 1., 1., 2 * 0.5) ) / ( Rm(r_dstk, rD_kpi, d_dstk + M_PI, dD_kpi, 1., 1., 2 * 0.5) + Rp(r_dstk, rD_kpi, d_dstk + M_PI, dD_kpi, 1., 1., 2 * 0.5) ) );
  // B -> Dpi
  // D -> Kpi
  ll1 += meas.at("Babar_1006.4241_RDpi").logweight( 0.5 * ( Rm(r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 2 * 0.5) +  Rp(r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 2 * 0.5) ) );
  ll1 += meas.at("Babar_1006.4241_ADpi").logweight( ( Rm(r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 2 * 0.5) -  Rp(r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 2 * 0.5) ) / ( Rm(r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 2 * 0.5) +  Rp(r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 2 * 0.5) ) );
  // B -> [Dpi0]_Dstarpi
  // D -> Kpi
  ll1 += meas.at("Babar_1006.4241_RDstarpipi0").logweight( 0.5 * (  Rm(r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 2 * 0.5) +  Rp(r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 2 * 0.5)  )  );
  ll1 += meas.at("Babar_1006.4241_ADstarpipi0").logweight( (  Rm(r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 2 * 0.5) -  Rp(r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 2 * 0.5)  ) / (  Rm(r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 2 * 0.5) +  Rp(r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 2 * 0.5)  ) );
  // B -> [Dpig]_Dstarpi
  // D -> Kpi
  ll1 += meas.at("Babar_1006.4241_RDstarpig").logweight( 0.5 * ( Rm(r_dstpi, rD_kpi, d_dstpi + M_PI, dD_kpi, 1., 1., 2 * 0.5) + Rp(r_dstpi, rD_kpi, d_dstpi + M_PI, dD_kpi, 1., 1., 2 * 0.5) ) );
  ll1 += meas.at("Babar_1006.4241_ADstarpig").logweight( ( Rm(r_dstpi, rD_kpi, d_dstpi + M_PI, dD_kpi, 1., 1., 2 * 0.5) - Rp(r_dstpi, rD_kpi, d_dstpi + M_PI, dD_kpi, 1., 1., 2 * 0.5) ) / ( Rm(r_dstpi, rD_kpi, d_dstpi + M_PI, dD_kpi, 1., 1., 2 * 0.5) + Rp(r_dstpi, rD_kpi, d_dstpi + M_PI, dD_kpi, 1., 1., 2 * 0.5) ) );


  // https://arxiv.org/pdf/1104.4472
  // B -> DK
  // D -> Kpipi0
  ll1 += meas.at("Babar_1104.4472_Rp_DK_Kpipi0").logweight(Rp(r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, 2 * 0.5));
  ll1 += meas.at("Babar_1104.4472_Rm_DK_Kpipi0").logweight(Rm(r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, 2 * 0.5));


  // https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.105.121801
  // B -> D(star)K(star) BPGGSZ
  // D -> K0spipi + D -> K0sKK
  corr.ResizeTo(12);
  corr(0) = r_dk * cos(d_dk - g);
  corr(1) = r_dk * sin(d_dk - g);
  corr(2) = r_dk * cos(d_dk + g);
  corr(3) = r_dk * sin(d_dk + g);
  corr(4) = r_dstk * cos(d_dstk - g);
  corr(5) = r_dstk * sin(d_dstk - g);
  corr(6) = r_dstk * cos(d_dstk + g);
  corr(7) = r_dstk * sin(d_dstk + g);
  corr(8) = k_dkst * r_dstk * cos(d_dkst - g);
  corr(9) = k_dkst * r_dkst * sin(d_dkst - g);
  corr(10) = k_dkst * r_dkst * cos(d_dkst + g);
  corr(11) = k_dkst * r_dkst * sin(d_dkst + g);
  ll1 += corrmeas.at("Babar_PRL105.121801").logweight(corr);


  //-------------------------------------------------  CDF measurements  -------------------------------------------------------------------------

  // https://journals.aps.org/prd/pdf/10.1103/PhysRevD.81.031105
  // B -> DK
  // D -> KK, D-> pipi
  ll1 += meas.at("CDF_PRD81_031105_ACPp").logweight(Acp(r_dk, d_dk, 1., 1., 2*0.5));
  ll1 += meas.at("CDF_PRD81_031105_RCPp").logweight(Rcp_h(r_dk, d_dk, r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 1., 2 * 0.5) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 1., 2 * 0.5));


  // https://arxiv.org/pdf/1108.5765
  // B -> DK
  // D -> Kpi
  ll1 += meas.at("CDF_1108.5765_RDK").logweight(Rads(r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 1.));
  ll1 += meas.at("CDF_1108.5765_ADK").logweight(Asup(r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 1.));
  // B -> Dpi
  // D -> Kpi
  ll1 += meas.at("CDF_1108.5765_RDpi").logweight(Rads(r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 1.));
  ll1 += meas.at("CDF_1108.5765_ADpi").logweight(Asup(r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 1.));

  //--------------------------------------------------------------------------------------------------------------------------

  //--------------------------------------------------------------------------------------------------------------------------


  // GLW: D -> KK, pipi; ADS: D -> Kpi
  // https://arxiv.org/pdf/2012.09903.pdf
  // Observables 8:
  corr.ResizeTo(8);
  corr(0) = acp_dk_uid0;
  corr(1) = acp_dpi_uid0;
  corr(2) = afav_dk_uid0;
  corr(3) = rcp_uid0;
  corr(4) = rm_dk_uid0;
  corr(5) = rm_dpi_uid0;
  corr(6) = rp_dk_uid0;
  corr(7) = rp_dpi_uid0;
  ll1 += corrmeas.at("UID0").logweight(corr);


  // GLW: D -> KKpipi, D -> 4pi
  // https://arxiv.org/abs/2301.10328
  // 6 Observables
  corr.ResizeTo(6);
  corr(0) = acp_dk_kkpipi_230110328;
  corr(1) = acp_dpi_kkpipi_230110328;
  corr(2) = acp_dk_pipipipi_230110328;
  corr(3) = acp_dpi_pipipipi_230110328;
  corr(4) = rcp_kpi_kkpipi_230110328;
  corr(5) = rcp_kpi_pipipipi_230110328;
  ll1 += corrmeas.at("GLW_2301.10328").logweight(corr);


  // GLW: D -> pipipi0, KKpi0; D -> Kpipi0
  //https://arxiv.org/pdf/2112.10617
  // Observables 11:
  corr.ResizeTo(11);
  corr(0) = rcp_kkpi0_211210617;
  corr(1) = rcp_pipipi0_211210617;
  corr(2) = afav_dk_kpipi0_211210617;
  corr(3) = acp_dk_kkpi0_211210617;
  corr(4) = acp_dk_pipipi0_211210617;
  corr(5) = acp_dpi_kkpi0_211210617;
  corr(6) = acp_dpi_pipipi0_211210617;
  corr(7) = rp_dk_211210617;
  corr(8) = rm_dk_211210617;
  corr(9) = rp_dpi_211210617;
  corr(10) = rm_dpi_211210617;
  ll1 += corrmeas.at("2112.10617").logweight(corr);


  // ADS: D -> K0sKpi
  // https://arxiv.org/pdf/2002.08858
  // Observables 7:
  corr.ResizeTo(7);
  corr(0) = afav_dpi_kskpi_uid4;
  corr(1) = asup_dpi_kskpi_uid4;
  corr(2) = afav_dk_kskpi_uid4;
  corr(3) = asup_dk_kskpi_uid4;
  corr(4) = rfavsup_dpi_kskpi_uid4;
  corr(5) = rfav_dkdpi_kskpi_uid4;
  corr(6) = rsup_dkdpi_kskpi_uid4;
  ll1 += corrmeas.at("UID4").logweight(corr);


  // GLW: D -> KK, D -> K0spi0
  // https://arxiv.org/abs/2308.05048
  // Observables 4:
  corr.ResizeTo(4);
  corr(0) = Acp(r_dk, d_dk, 1., 0., 1.);
  corr(1) = Rcp_h(r_dk, d_dk, r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 0., 1.) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 0., 1.);
  corr(2) = Acp(r_dk, d_dk, 1., 1., 1.);
  corr(3) = Rcp_h(r_dk, d_dk, r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 1., 1.) / Rcp_h(r_dpi, d_dpi, r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 1., 1.);
  ll1 += corrmeas.at("2308.05048").logweight(corr);


  // ADS: D -> Kpipi0
  // https://arxiv.org/pdf/1310.1741
  // Observables 4:
  corr.ResizeTo(4);
  corr(0) = Rads(r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, 1.);
  corr(1) = Asup(r_dk, rD_kpipi0, d_dk, dD_kpipi0, 1., kD_kpipi0, 1.);
  corr(2) = Rads(r_dpi, rD_kpipi0, d_dpi, dD_kpipi0, 1., kD_kpipi0, 1.);
  corr(3) = Asup(r_dpi, rD_kpipi0, d_dpi, dD_kpipi0, 1., kD_kpipi0, 1.);
  ll1 += corrmeas.at("Belle_PRD88_2013").logweight(corr);


  // ADS: D -> Kpi
  // https://arxiv.org/abs/1103.5951
  // Observables 4:
  corr.ResizeTo(4);
  corr(0) = Rads(r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 1.);
  corr(1) = Asup(r_dk, rD_kpi, d_dk, dD_kpi, 1., 1., 1.);
  corr(2) = Rads(r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 1.);
  corr(3) = Asup(r_dpi, rD_kpi, d_dpi, dD_kpi, 1., 1., 1.);
  ll1 += corrmeas.at("Belle_PRL106_2011").logweight(corr);


  // K^*+- region fit: ADS: D -> K0sKpi
  // 2306.02940
  // Observables 7:
  corr.ResizeTo(7);
  corr(0) = afav_dk_kskpi_uid4;
  corr(1) = asup_dk_kskpi_uid4;
  corr(2) = afav_dpi_kskpi_uid4;
  corr(3) = asup_dpi_kskpi_uid4;
  corr(4) = rfav_dkdpi_kskpi_uid4;
  corr(5) = rsup_dkdpi_kskpi_uid4;
  corr(6) = rfavsup_dpi_kskpi_uid4;
  ll1 += corrmeas.at("2306.02940").logweight(corr);


  // GGSZ: D -> K3pi
  // https://arxiv.org/pdf/2209.03692
  // Observables 6:
  corr.ResizeTo(6);
  corr(0) = xp_dk_uid3;
  corr(1) = xm_dk_uid3;
  corr(2) = yp_dk_uid3;
  corr(3) = ym_dk_uid3;
  corr(4) = xi_x_dpi_uid3;
  corr(5) = xi_y_dpi_uid3;
  ll1 += corrmeas.at("2209.03692").logweight(corr);


  // GGSZ D -> K0spipipi0
  // https://arxiv.org/pdf/1908.09499
  // Observables 8:
  corr.ResizeTo(8);
  corr(0) = xm_dk_uid3;
  corr(1) = ym_dk_uid3;
  corr(2) = xp_dk_uid3;
  corr(3) = yp_dk_uid3;
  corr(4) = r_dpi * cos(d_dpi - g);
  corr(5) = r_dpi * sin(d_dpi - g);
  corr(6) = r_dpi * cos(d_dpi + g);
  corr(7) = r_dpi * sin(d_dpi + g);
  ll1 += corrmeas.at("1908.09499").logweight(corr);


  // GGSZ: D -> K0spipi, D -> K0sKK
  // https://arxiv.org/abs/2110.12125
  // Observables 6:
  corr.ResizeTo(6);
  corr(0) = xm_dk_uid3;
  corr(1) = ym_dk_uid3;
  corr(2) = xp_dk_uid3;
  corr(3) = yp_dk_uid3;
  corr(4) = xi_x_dpi_uid3;
  corr(5) = xi_y_dpi_uid3;
  ll1 += corrmeas.at("2110.12125").logweight(corr);


  // GGSZ: D -> K0sKK, D -> K0spipi
  // https://arxiv.org/abs/2301.10328
  // 6 Observables
  corr.ResizeTo(6);
  corr(0) = xm_dk_uid3;
  corr(1) = ym_dk_uid3;
  corr(2) = xp_dk_uid3;
  corr(3) = yp_dk_uid3;
  corr(4) = xi_x_dpi_uid3;
  corr(5) = xi_y_dpi_uid3;
  ll1 += corrmeas.at("GGSZ_2301.10328").logweight(corr);


  // If combiining chargedB modes only
  // If combining all the modes see the neutralBd section since they are all correlated
  // https://arxiv.org/pdf/2010.08483 + https://arxiv.org/abs/2310.04277 +  https://arxiv.org/abs/2311.10434 + LHCB-PAPER-2024-023
  if (comb == 0)
  {

    // 4. GGSZ LHCb ChargedB
    // Observables 22:
    corr.ResizeTo(22);

    // D -> K0spipi, D -> K0sKK
    // https://arxiv.org/pdf/2010.08483
    // Observables 6:
    corr(0) = xm_dk_uid3;
    corr(1) = ym_dk_uid3;
    corr(2) = xp_dk_uid3;
    corr(3) = yp_dk_uid3;
    corr(4) = xi_x_dpi_uid3;
    corr(5) = xi_y_dpi_uid3;

    //----------------------------------------------------------------------------------------------------------------------------------

    //-------------------------------------------------  Bpm -> D*hpm -------------------------------------------------------------------------

    // GGSZ: D -> K0spipi, D -> K0sKK
    // https://arxiv.org/abs/2310.04277 LHCb
    // Observables 6:
    corr(6) = r_dstk * cos(d_dstk + g);
    corr(7) = r_dstk * cos(d_dstk - g);
    corr(8) = r_dstk * sin(d_dstk + g);
    corr(9) = r_dstk * sin(d_dstk - g);
    corr(10) = r_dstpi / r_dstk * cos(d_dstpi - d_dstk);
    corr(11) = r_dstpi / r_dstk * sin(d_dstpi - d_dstk);

    // GGSZ: D -> K0spipi, D -> K0sKK
    // https://arxiv.org/abs/2311.10434 LHCb
    // Observables 6:
    corr(12) = r_dstk * cos(d_dstk - g);
    corr(13) = r_dstk * sin(d_dstk - g);
    corr(14) = r_dstk * cos(d_dstk + g);
    corr(15) = r_dstk * sin(d_dstk + g);
    corr(16) = r_dstpi / r_dstk * cos(d_dstpi - d_dstk);
    corr(17) = r_dstpi / r_dstk * sin(d_dstpi - d_dstk);

    //----------------------------------------------------------------------------------------------------------------------------------

    //-------------------------------------------------  Bpm -> DK*pm -------------------------------------------------------------------------

    // GGSZ: D -> K0spipi, D -> K0sKK
    // LHCB-PAPER-2024-023
    // Observables 4:
    corr(18) = r_dkst * cos(d_dkst - g);
    corr(19) = r_dkst * sin(d_dkst - g);
    corr(20) = r_dkst * cos(d_dkst + g);
    corr(21) = r_dkst * sin(d_dkst + g);
    ll1 += corrmeas.at("GGSZ_LHCb_Cb").logweight(corr);

  }

  //----------------------------------------------------------------------------------------------------------------------------------

  //-------------------------------------------------  Bpm -> D*hpm -------------------------------------------------------------------------

  // GLW: D -> KK, D -> pipi; ADS: D -> Kpi
  // https://arxiv.org/abs/2012.09903
  // Observables 18:
  corr.ResizeTo(18);
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
  corr(16) = afav_dstpi_dg_uid5;
  corr(17) = afav_dstpi_dp_uid5;
  ll1 += corrmeas.at("UID5").logweight(corr);


  // GLW: D -> KK, D -> pipi, D -> K0spi0, ....
  // https://arxiv.org/pdf/hep-ex/0601032
  // Observables 4:
  corr.ResizeTo(4);
  corr(0) = Acp(r_dstk, d_dstk, 1., 0., 1.);
  corr(1) = Rcp_h(r_dstk, d_dstk, r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 0., 1.) / Rcp_h(r_dstpi, d_dstpi, r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 0., 1.);
  corr(2) = Acp(r_dstk, d_dstk, 1., 1., 1.);
  corr(3) = Rcp_h(r_dstk, d_dstk, r_dstk, rD_kpi, d_dstk, dD_kpi, 1., 1., 1., 1.) / Rcp_h(r_dstpi, d_dstpi, r_dstpi, rD_kpi, d_dstpi, dD_kpi, 1., 1., 1., 1.);
  ll1 += corrmeas.at("Belle_PRD73_2006").logweight(corr);


  // GGSZ: D -> K0spipi
  // https://arxiv.org/abs/1003.3360
  // Observables 8:
  corr.ResizeTo(8);
  corr(0) = r_dstk * cos(d_dstk - g);
  corr(1) = r_dstk * sin(d_dstk - g);
  corr(2) = r_dstk * cos(d_dstk + g);
  corr(3) = r_dstk * sin(d_dstk + g);
  corr(4) = -r_dstk * cos(d_dstk - g);
  corr(5) = -r_dstk * sin(d_dstk - g);
  corr(6) = -r_dstk * cos(d_dstk + g);
  corr(7) = -r_dstk * sin(d_dstk + g);
  ll1 += corrmeas.at("Belle_PRD81_2010").logweight(corr);


  //----------------------------------------------------------------------------------------------------------------------------------

  //-------------------------------------------------  Bpm -> DKstpm -------------------------------------------------------------------------

  // Bpm -> DK^*pm
  // https://arxiv.org/pdf/hep-ex/0604054 (Belle)
  corr.ResizeTo(4);
  corr(0) = r_dkst * cos(d_dkst + g);
  corr(1) = r_dkst * sin(d_dkst + g);
  corr(2) = r_dkst * cos(d_dkst - g);
  corr(3) = r_dkst * sin(d_dkst - g);
  ll1 += corrmeas.at("0604054").logweight(corr);

  // GLW: D -> KK, D -> pipi, D -> 4pi; ADS: D -> Kpi, D -> K3pi
  // LHCb-PAPER-2024-023
  // Observables 12:
  corr.ResizeTo(12);
  corr(0) = afav_dkst_kpi;
  corr(1) = acp_dkst_kk;
  corr(2) = acp_dkst_pipi;
  corr(3) = asup_dkst_kpi;
  corr(4) = rcp_dkst_kk;
  corr(5) = rcp_dkst_pipi;
  corr(6) = rsup_dkst_kpi;
  corr(7) = afav_dkst_k3pi;
  corr(8) = acp_dkst_pipipipi;
  corr(9) = asup_dkst_k3pi;
  corr(10) = rcp_dkst_pipipipi;
  corr(11) = rsup_dkst_k3pi;
  ll1 += corrmeas.at("LHCB-PAPER-2024-023-GLWADS").logweight(corr);


  // Coherence factor kappakstpm
  // https://arxiv.org/pdf/1709.05855
  ll1 += meas.at("UID24").logweight(k_dkst_uid24);

  //----------------------------------------------------------------------------------------------------------------------------------

  //-------------------------------------------------  Bpm -> Dhpipi -------------------------------------------------------------------------

  // GLW: D -> KK, D -> pipi; ADS: D -> Kpi
  // https://arxiv.org/pdf/1505.07044
  // Observables 11:
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

  return ll1;
}
// ---------------------------------------------------------

// ---------------------------------------------------------
double MixingModel::Calculate_neutralBdobservables()
{

  double ll2;
  ll2 = 0.;

  // 26. PDF: dkstzcoherence (UID25)
  //  1 Osservabile
  k_dkstz_uid25 = k_dkstz;

  //----------------------------------------------- Calculating time Integrated B0d observables -------------------------------------------------------------------------

  // PDF: glwads-dkst-hh-Kpi-h3pi-dmix (2401.17934)
  // Observables 12:
  acp_dkstz_kk_240117934Bd = Acp(r_dkstz, d_dkstz, k_dkstz, 1., 1.34);
  acp_dkstz_pipi_240117934Bd = Acp(r_dkstz, d_dkstz, k_dkstz, 1., 1.34);
  rcp_dkstz_kk_240117934Bd = Rcp_h(r_dkstz, d_dkstz, r_dkstz, rD_kpi, d_dkstz, dD_kpi, k_dkstz, 1., 1., 1.34);
  rcp_dkstz_pipi_240117934Bd = Rcp_h(r_dkstz, d_dkstz, r_dkstz, rD_kpi, d_dkstz, dD_kpi, k_dkstz, 1., 1., 1.34);
  acp_dkstz_4pi_240117934Bd = Acp(r_dkstz, d_dkstz, k_dkstz, F_pipipipi, 1.34);
  rcp_dkstz_4pi_240117934Bd = Rcp_h(r_dkstz, d_dkstz, r_dkstz, rD_k3pi, d_dkstz, dD_k3pi, k_dkstz, kD_k3pi, F_pipipipi, 1.34);
  rp_dkstz_kpi_240117934Bd = Rp(r_dkstz, rD_kpi, d_dkstz, dD_kpi, k_dkstz, 1., 1.34);
  rm_dkstz_kpi_240117934Bd = Rm(r_dkstz, rD_kpi, d_dkstz, dD_kpi, k_dkstz, 1., 1.34);
  rp_dkstz_k3pi_240117934Bd = Rp(r_dkstz, rD_k3pi, d_dkstz, dD_k3pi, k_dkstz, kD_k3pi, 1.34);
  rm_dkstz_k3pi_240117934Bd = Rm(r_dkstz, rD_k3pi, d_dkstz, dD_k3pi, k_dkstz, kD_k3pi, 1.34);
  afav_dkstz_kpi_240117934Bd = Afav(r_dkstz, rD_kpi, d_dkstz, dD_kpi, k_dkstz, 1., 1.34);
  afav_dkstz_k3pi_240117934Bd = Afav(r_dkstz, rD_k3pi, d_dkstz, dD_k3pi, k_dkstz, kD_k3pi, 1.34);


  xm_dkstz_230905514 = r_dkstz * cos(d_dkstz - g);
  ym_dkstz_230905514 = r_dkstz * sin(d_dkstz - g);
  xp_dkstz_230905514 = r_dkstz * cos(d_dkstz + g);
  yp_dkstz_230905514 = r_dkstz * sin(d_dkstz + g);

  //----------------------------------------------- Calculating time dependent B0d observables -------------------------------------------------------------------------


  s_dmpi_uid12 = -(2 * l_dmpi * sin(d_dmpi - (phi_d + g))) / (1 + l_dmpi * l_dmpi);
  sb_dmpi_uid12 = (2 * l_dmpi * sin(d_dmpi + (phi_d + g))) / (1 + l_dmpi * l_dmpi);

  //-------------------------------------- Babar Measurements -------------------------------------------------------------------------

  // https://arxiv.org/pdf/hep-ex/0602049
  double a_Dpi = -2 * l_dmpi / (1 + l_dmpi * l_dmpi) * sin(phi_d + g) * cos(d_dmpi);
  double c_Dpi = -2 * l_dmpi / (1 + l_dmpi * l_dmpi) * cos(phi_d + g) * sin(d_dmpi);
  double a_Dstarpi = -2 * l_dstarmpi / (1 + l_dstarmpi * l_dstarmpi) * sin(phi_d + g) * cos(d_dstarmpi);
  double c_Dstarpi = -2 * l_dstarmpi / (1 + l_dstarmpi * l_dstarmpi) * cos(phi_d + g) * sin(d_dstarmpi);
  double a_Drho = -2 * l_dmrho / (1 + l_dmrho * l_dmrho) * sin(phi_d + g) * cos(d_dmrho);
  double c_Drho = -2 * l_dmrho / (1 + l_dmrho * l_dmrho) * cos(phi_d + g) * sin(d_dmrho);


  //----------------------------------------------- Contribution to the LogLikelihood -------------------------------------------------------------------------

  TVectorD corr(8);

  if (comb == 1)
  { // If treating it seprately from the Bs counterpart

    corr.ResizeTo(12);
    corr(0) = afav_dkstz_kpi_240117934Bd;
    corr(1) = rp_dkstz_kpi_240117934Bd;
    corr(2) = rm_dkstz_kpi_240117934Bd;
    corr(3) = afav_dkstz_k3pi_240117934Bd;
    corr(4) = rp_dkstz_k3pi_240117934Bd;
    corr(5) = rm_dkstz_k3pi_240117934Bd;
    corr(6) = acp_dkstz_kk_240117934Bd;
    corr(7) = rcp_dkstz_kk_240117934Bd;
    corr(8) = acp_dkstz_pipi_240117934Bd;
    corr(9) = rcp_dkstz_pipi_240117934Bd;
    corr(10) = acp_dkstz_4pi_240117934Bd;
    corr(11) = rcp_dkstz_4pi_240117934Bd;
    ll2 += corrmeas.at("2401.17934Bd").logweight(corr);

    corr.ResizeTo(4);
    corr(0) = xp_dkstz_230905514;
    corr(1) = xm_dkstz_230905514;
    corr(2) = yp_dkstz_230905514;
    corr(3) = ym_dkstz_230905514;
    ll2 += corrmeas.at("2309.05514").logweight(corr);

  }
  else if (comb == 3)
  {

    corr.ResizeTo(26);

    corr(0) = xm_dk_uid3;
    corr(1) = ym_dk_uid3;
    corr(2) = xp_dk_uid3;
    corr(3) = yp_dk_uid3;
    corr(4) = xi_x_dpi_uid3;
    corr(5) = xi_y_dpi_uid3;

    corr(6) = r_dstk * cos(d_dstk + g);
    corr(7) = r_dstk * cos(d_dstk - g);
    corr(8) = r_dstk * sin(d_dstk + g);
    corr(9) = r_dstk * sin(d_dstk - g);
    corr(10) = r_dstpi / r_dstk * cos(d_dstpi - d_dstk);
    corr(11) = r_dstpi / r_dstk * sin(d_dstpi - d_dstk);

    corr(12) = r_dstk * cos(d_dstk - g);
    corr(13) = r_dstk * sin(d_dstk - g);
    corr(14) = r_dstk * cos(d_dstk + g);
    corr(15) = r_dstk * sin(d_dstk + g);
    corr(16) = r_dstpi / r_dstk * cos(d_dstpi - d_dstk);
    corr(17) = r_dstpi / r_dstk * sin(d_dstpi - d_dstk);

    corr(18) = r_dkst * cos(d_dkst - g);
    corr(19) = r_dkst * sin(d_dkst - g);
    corr(20) = r_dkst * cos(d_dkst + g);
    corr(21) = r_dkst * sin(d_dkst + g);

    corr(22) = xp_dkstz_230905514;
    corr(23) = xm_dkstz_230905514;
    corr(24) = yp_dkstz_230905514;
    corr(25) = ym_dkstz_230905514;
    ll2 += corrmeas.at("DKst0Pcomb").logweight(corr);

  }


  //  https://arxiv.org/pdf/1509.01098 (Belle)
  corr.ResizeTo(4);
  corr(0) = r_dkstz * cos( d_dkstz + g );
  corr(1) = r_dkstz * sin( d_dkstz + g );
  corr(2) = r_dkstz * cos( d_dkstz - g );
  corr(3) = r_dkstz * sin( d_dkstz - g );
  ll2 += corrmeas.at("1509.01098").logweight(corr);


  corr.ResizeTo(2);
  corr(0) = s_dmpi_uid12;
  corr(1) = sb_dmpi_uid12;
  ll2 += corrmeas.at("UID12").logweight(corr);


  //-------------------------------------- Babar Measurements -------------------------------------------------------------------------
  // https://arxiv.org/pdf/hep-ex/0602049
  ll2+= meas.at("0602049_aDpi").logweight(a_Dpi);
  ll2+= meas.at("0602049_cDpi").logweight(c_Dpi);
  ll2+= meas.at("0602049_aDstarpi").logweight(a_Dstarpi);
  ll2+=meas.at("0602049_cDstarpi").logweight(c_Dstarpi);
  ll2+=meas.at("0602049_aDrho").logweight(a_Drho);
  ll2+=meas.at("0602049_cDrho").logweight(c_Drho);

  // https://arxiv.org/pdf/hep-ex/0504035
  ll2+=meas.at("0504035_aDstarpi").logweight(a_Dstarpi);
  ll2+=meas.at("0504035_cDstarpi").logweight(c_Dstarpi);


  ll2 += meas.at("UID25").logweight(k_dkstz_uid25);

  ll2 += meas.at("UID27").logweight(sin(phi_d));


  //-------------------------------------- Belle Measurements -------------------------------------------------------------------------

  // https://arxiv.org/pdf/hep-ex/0604013 (HFLAV conversion)
  ll2+= meas.at("0604013_aDpi").logweight(a_Dpi);
  ll2+= meas.at("0604013_cDpi").logweight(c_Dpi);
  ll2+= meas.at("0604013_aDstarpi").logweight(a_Dstarpi);
  ll2+=meas.at("0604013_cDstarpi").logweight(c_Dstarpi);

  // https://arxiv.org/pdf/1102.0888
  ll2+= meas.at("11020888_aDstarpi").logweight(a_Dstarpi);
  ll2+=meas.at("11020888_cDstarpi").logweight(c_Dstarpi);



  return ll2;
}
// ---------------------------------------------------------

// ---------------------------------------------------------
double MixingModel::Calculate_neutralBsobservables()
{

  double ll2;
  ll2 = 0.;

  //----------------------------------------------- Calculating time Integrated B0d observables -------------------------------------------------------------------------

  acp_dkstz_kk_240117934Bs = Acp(r_dkstzs, d_dkstzs, k_dkstzs, 1., 1.34);
  acp_dkstz_pipi_240117934Bs = Acp(r_dkstzs, d_dkstzs, k_dkstzs, 1., 1.34);
  rcp_dkstz_kk_240117934Bs = Rcp_h(r_dkstzs, d_dkstzs, r_dkstzs, rD_kpi, d_dkstzs, dD_kpi, k_dkstzs, 1., 1., 1.34);
  rcp_dkstz_pipi_240117934Bs = Rcp_h(r_dkstzs, d_dkstzs, r_dkstzs, rD_kpi, d_dkstzs, dD_kpi, k_dkstzs, 1., 1., 1.34);
  acp_dkstz_4pi_240117934Bs = Acp(r_dkstzs, d_dkstzs, k_dkstzs, F_pipipipi, 1.34);
  rcp_dkstz_4pi_240117934Bs = Rcp_h(r_dkstzs, d_dkstzs, r_dkstzs, rD_k3pi, d_dkstzs, dD_k3pi, k_dkstzs, kD_k3pi, F_pipipipi, 1.34);
  rp_dkstz_kpi_240117934Bs = Rp(r_dkstzs, rD_kpi, d_dkstzs, dD_kpi, k_dkstzs, 1., 1.34);
  rm_dkstz_kpi_240117934Bs = Rm(r_dkstzs, rD_kpi, d_dkstzs, dD_kpi, k_dkstzs, 1., 1.34);
  rp_dkstz_k3pi_240117934Bs = Rp(r_dkstzs, rD_k3pi, d_dkstzs, dD_k3pi, k_dkstzs, kD_k3pi, 1.34);
  rm_dkstz_k3pi_240117934Bs = Rm(r_dkstzs, rD_k3pi, d_dkstzs, dD_k3pi, k_dkstzs, kD_k3pi, 1.34);
  afav_dkstz_kpi_240117934Bs = Afav(r_dkstzs, rD_kpi, d_dkstzs, dD_kpi, k_dkstzs, 1., 1.34);
  afav_dkstz_k3pi_240117934Bs = Afav(r_dkstzs, rD_k3pi, d_dkstzs, dD_k3pi, k_dkstzs, kD_k3pi, 1.34);

  //----------------------------------------------- Calculating time dependent B0s observables -------------------------------------------------------------------------

  c_dsk_uid10 = (1 - l_dsk * l_dsk) / (1 + l_dsk * l_dsk);
  d_dsk_uid10 = -(2 * l_dsk * cos(d_dsk - (g + phis))) / (1 + l_dsk * l_dsk); // phis = - 2 beta_s
  db_dsk_uid10 = -(2 * l_dsk * cos(d_dsk + (g + phis))) / (1 + l_dsk * l_dsk);
  s_dsk_uid10 = (2 * l_dsk * sin(d_dsk - (g + phis))) / (1 + l_dsk * l_dsk);
  sb_dsk_uid10 = -(2 * l_dsk * sin(d_dsk + (g + phis))) / (1 + l_dsk * l_dsk);


  c_dskpipi_uid11 = (1 - l_dskpipi * l_dskpipi) / (1 + l_dskpipi * l_dskpipi);
  d_dskpipi_uid11 = -(2 * k_dskpipi * l_dskpipi * cos(d_dskpipi - (g + phis))) / (1 + l_dskpipi * l_dskpipi);
  db_dskpipi_uid11 = -(2 * k_dskpipi * l_dskpipi * cos(d_dskpipi + (g + phis))) / (1 + l_dskpipi * l_dskpipi);
  s_dskpipi_uid11 = (2 * k_dskpipi * l_dskpipi * sin(d_dskpipi - (g + phis))) / (1 + l_dskpipi * l_dskpipi);
  sb_dskpipi_uid11 = -(2 * k_dskpipi * l_dskpipi * sin(d_dskpipi + (g + phis))) / (1 + l_dskpipi * l_dskpipi);

  phis_uid26 = phis; // -2 betas

  //----------------------------------------------- Contribution to the LogLikelihood -------------------------------------------------------------------------

  TVectorD corr(8);

  if (comb == 2)
  { // If treating it separately from Bd counterpart

    // PDF: glwads-dkst-hh-Kpi-h3pi-dmix (2401.17934)
    // Observables 12:
    corr.ResizeTo(12);
    corr(0) = afav_dkstz_kpi_240117934Bs;
    corr(1) = rp_dkstz_kpi_240117934Bs;
    corr(2) = rm_dkstz_kpi_240117934Bs;
    corr(3) = afav_dkstz_k3pi_240117934Bs;
    corr(4) = rp_dkstz_k3pi_240117934Bs;
    corr(5) = rm_dkstz_k3pi_240117934Bs;
    corr(6) = acp_dkstz_kk_240117934Bs;
    corr(7) = rcp_dkstz_kk_240117934Bs;
    corr(8) = acp_dkstz_pipi_240117934Bs;
    corr(9) = rcp_dkstz_pipi_240117934Bs;
    corr(10) = acp_dkstz_4pi_240117934Bs;
    corr(11) = rcp_dkstz_4pi_240117934Bs;

    ll2 += corrmeas.at("2401.17934Bs").logweight(corr);
  }
  else if (comb == 3)
  { // When using also the Bd counterpart

    // PDF: glwads-dkst-hh-Kpi-h3pi-dmix (2401.17934)
    // Observables 24:
    corr.ResizeTo(24); // Bd  part
    corr(0) = afav_dkstz_kpi_240117934Bd;
    corr(1) = rp_dkstz_kpi_240117934Bd;
    corr(2) = rm_dkstz_kpi_240117934Bd;
    corr(3) = afav_dkstz_k3pi_240117934Bd;
    corr(4) = rp_dkstz_k3pi_240117934Bd;
    corr(5) = rm_dkstz_k3pi_240117934Bd;
    corr(6) = acp_dkstz_kk_240117934Bd;
    corr(7) = rcp_dkstz_kk_240117934Bd;
    corr(8) = acp_dkstz_pipi_240117934Bd;
    corr(9) = rcp_dkstz_pipi_240117934Bd;
    corr(10) = acp_dkstz_4pi_240117934Bd;
    corr(11) = rcp_dkstz_4pi_240117934Bd;

    // Bs part
    corr(12) = afav_dkstz_kpi_240117934Bs;
    corr(13) = rp_dkstz_kpi_240117934Bs;
    corr(14) = rm_dkstz_kpi_240117934Bs;
    corr(15) = afav_dkstz_k3pi_240117934Bs;
    corr(16) = rp_dkstz_k3pi_240117934Bs;
    corr(17) = rm_dkstz_k3pi_240117934Bs;
    corr(18) = acp_dkstz_kk_240117934Bs;
    corr(19) = rcp_dkstz_kk_240117934Bs;
    corr(20) = acp_dkstz_pipi_240117934Bs;
    corr(21) = rcp_dkstz_pipi_240117934Bs;
    corr(22) = acp_dkstz_4pi_240117934Bs;
    corr(23) = rcp_dkstz_4pi_240117934Bs;

    ll2 += corrmeas.at("2401.17934").logweight(corr);

  }


  corr.ResizeTo(5);
  corr(0) = c_dsk_uid10;
  corr(1) = d_dsk_uid10;
  corr(2) = db_dsk_uid10;
  corr(3) = s_dsk_uid10;
  corr(4) = sb_dsk_uid10;
  ll2 += corrmeas.at("BSDSKRun1").logweight(corr);
  ll2 += corrmeas.at("BSDSKRun2").logweight(corr);


  corr.ResizeTo(5);
  corr(0) = c_dskpipi_uid11;
  corr(1) = d_dskpipi_uid11;
  corr(2) = db_dskpipi_uid11;
  corr(3) = s_dskpipi_uid11;
  corr(4) = sb_dskpipi_uid11;
  ll2 += corrmeas.at("UID11").logweight(corr);


  ll2 += meas.at("UID26").logweight(phis_uid26); // -2betas

  return ll2;
}
// ---------------------------------------------------------

// ---------------------------------------------------------
double MixingModel::Calculate_time_dependent_Dobservables()
{

  double ll3;
  ll3 = 0.;

  xcp_uid14 = xcp;
  ycp_uid14 = ycp;
  dx_uid14 = dx;
  dy_uid14 = dy;

  ycp_uid28 = ycp + 0.5 * rD_kpi * (cos(phi) * (-y * cos(dD_kpi) - x * sin(dD_kpi)) * (qop + 1. / qop) + sin(phi) * (-y * sin(dD_kpi) + x * cos(dD_kpi)) * (qop - 1. / qop));

  DY_uid29 = 0.5 * (-y * cos(phi) * (qop - 1. / qop) + x * sin(phi) * (qop + 1. / qop));

  Rdp_uid30 = rD_kpi * rD_kpi * (1 + AD);
  yp_uid30 = y_plus(dD_kpi);
  xpsq_uid30 = x_plus(dD_kpi) * x_plus(dD_kpi);
  Rdm_uid30 = rD_kpi * rD_kpi * (1 - AD);
  ym_uid30 = y_minus(dD_kpi);
  xmsq_uid30 = x_minus(dD_kpi) * x_minus(dD_kpi);

  Akpi_BESIII = (-2 * rD_kpi * cos(dD_kpi) + y) / (1 + rD_kpi * rD_kpi);
  Akpi_kpipi0_BESIII = (F_pipipi0 * (-2 * rD_kpi * cos(dD_kpi) + y)) / (1 + rD_kpi * rD_kpi + (1 - F_pipipi0) * (2 * rD_kpi * cos(dD_kpi) + y));

  xi_x_BESIII = rD_kpi * cos(dD_kpi);
  xi_y_BESIII = rD_kpi * sin(dD_kpi);

  double tKKpitaggedOverTauD = tavepitaggedOverTauD + 0.5 * DeltatpitaggedOverTauD;
  double tKKmutaggedOverTauD = tavemutaggedOverTauD + 0.5 * DeltatmutaggedOverTauD;
  double tpipipitaggedOverTauD = tavepitaggedOverTauD - 0.5 * DeltatpitaggedOverTauD;
  double tpipimutaggedOverTauD = tavemutaggedOverTauD - 0.5 * DeltatmutaggedOverTauD;
  double DYKK = DY_uid29 + 0.5 * DYKKmDYpipi;
  double DYpipi = DY_uid29 - 0.5 * DYKKmDYpipi;
  double DeltaACP_pitagged = adKK - adpipi + tKKpitaggedOverTauD * DYKK - tpipipitaggedOverTauD * DYpipi;
  double DeltaACP_mutagged = adKK - adpipi + tKKmutaggedOverTauD * DYKK - tpipimutaggedOverTauD * DYpipi;
  double ACPKKDp = adKK + tKKCDp / tauD * DYKK;
  double ACPKKDs = adKK + tKKCDs / tauD * DYKK;

  //----------------------------------------------- Contribution to the LogLikelihood -------------------------------------------------------------------------
  TVectorD corr(8);

  // Delta ACP; Acp(KK); Run1 semileptonic tagging
  // https://arxiv.org/pdf/1405.2797
  // Observables 4:
  ll3 += meas.at("1405.2797_tKKOverTauD_DAcp").logweight(tauKK_DAcp_Run1_sl);
  ll3 += meas.at("1405.2797_tpipiOverTauD_DAcp").logweight(taupipi_DAcp_Run1_sl);
  ll3 += meas.at("1405.2797_tKKOverTauD_Acp").logweight(tauKK_Acp_Run1_sl);
  corr.ResizeTo(3);
  corr(0) = adKK - adpipi + tauKK_DAcp_Run1_sl * DYKK - taupipi_DAcp_Run1_sl * DYpipi; // DAcp Run1 sl
  corr(1) = adKK + tauKK_Acp_Run1_sl * DYKK; // Acp(KK) sl
  // ACP(KK); Run1 Hadronic tagging; Correlated due to the removing of the detection asymmetry
  // https://arxiv.org/pdf/1610.09476
  // Observables 1:
  corr(2) = adKK + tauKK_Acp_Run1_pi * DYKK; // Acp(KK) pi-tagged
  ll3 += corrmeas.at("1405.2797_1610.09476_Acp").logweight(corr);


  // tOverTauD; Run1 Hadronic tagging;
  // https://arxiv.org/pdf/1610.09476
  // Observables 1:
  ll3 += meas.at("1610.09476_tKKOverTauD").logweight(tauKK_Acp_Run1_pi);


  // Delta ACP; Run1 Hadronic tagging
  // https://arxiv.org/pdf/1602.03160
  // Observables 3:
  double DeltaACP_Run1_pitagged = adKK - adpipi + tauKK_DAcp_Run1_pi * DYKK - taupipi_DAcp_Run1_pi * DYpipi;
  ll3 += meas.at("1602.03160_DeltaACPpitagged").logweight(DeltaACP_Run1_pitagged);
  ll3 += meas.at("1602.03160_tKKOverTauD").logweight(tauKK_DAcp_Run1_pi);
  ll3 += meas.at("1602.03160_tpipiOverTauD").logweight(taupipi_DAcp_Run1_pi);


  // Delta ACP; Run2 Hadronic and Semileptonic tagging
  // https://arxiv.org/pdf/1903.08726
  // Observables 6:
  ll3 += meas.at("DeltaACPpitagged").logweight(DeltaACP_pitagged);
  ll3 += meas.at("DeltaACPmutagged").logweight(DeltaACP_mutagged);
  ll3 += meas.at("tavepitaggedOverTauD").logweight(tavepitaggedOverTauD);
  ll3 += meas.at("tavemutaggedOverTauD").logweight(tavemutaggedOverTauD);
  ll3 += meas.at("DeltatmutaggedOverTauD").logweight(DeltatmutaggedOverTauD);
  // tOverTauD; Run 2 Acp(KK); Correlated through reconstructed mean decay times
  // https://arxiv.org/pdf/2209.03179
  // Observables 2:
  corr.ResizeTo(3);
  corr(0) = tKKCDp;
  corr(1) = tKKCDs;
  corr(2) = DeltatpitaggedOverTauD;
  ll3 += corrmeas.at("tausforDACP").logweight(corr);


  // Run 2 Acp(KK)
  // https://arxiv.org/pdf/2209.03179
  // Observables 2:
  corr.ResizeTo(2);
  corr(0) = ACPKKDp;
  corr(1) = ACPKKDs;
  ll3 += corrmeas.at("ACPKK").logweight(corr);


  // yCP - yCP(Kpi)
  // HFLAV combo: https://hflav-eos.web.cern.ch/hflav-eos/charm/CKM23/results_mixing.html#kkpipi including latest https://arxiv.org/abs/2202.09106
  ll3 += meas.at("UID28").logweight(ycp_uid28);


  // ( DY(KK) + DY(pipi) )/2 LHCb combo Run1 + Run2
  // https://arxiv.org/pdf/2105.09889
  // Observables 1:
  ll3 += meas.at("UID29").logweight(DY_uid29);


  // ( DY(KK) - DY(pipi) ) LHCb combo Run1 + Run2
  // https://arxiv.org/pdf/2105.09889
  // Observables 1:
  ll3 += meas.at("DYKKmDYpipi").logweight(DYKKmDYpipi);


  // Agamma(KK) and Agamma(pipi) Full CDF
  // https://arxiv.org/pdf/1410.5435
  // Observables 2:
  ll3 += meas.at("1410.5435_AGammaKK").logweight(-DYKK);
  ll3 += meas.at("1410.5435_AGammapipi").logweight(-DYpipi);


  // tOverTauD CDF
  // https://arxiv.org/pdf/1111.5023
  // Observables 2:
  ll3 += meas.at("1111.5023_tKKOverTauD_CDF").logweight(tauKK_Acp_CDF);
  ll3 += meas.at("1111.5023_tpipiOverTauD_CDF").logweight(taupipi_Acp_CDF);


  // Acp(KK), Acp(pipi)
  // https://arxiv.org/pdf/1208.2517
  // Observables 2
  double ACPKKCDF = adKK + tauKK_Acp_CDF * DYKK;
  double ACPpipiCDF = adpipi + taupipi_Acp_CDF * DYpipi;
  ll3 += meas.at("1208.2517_AcpKK_CDF").logweight(ACPKKCDF);
  ll3 += meas.at("1208.2517_Acppipi_CDF").logweight(ACPpipiCDF);


  // Acp(KK), Acp(pipi) Babar
  // https://arxiv.org/pdf/0709.2715
  // Observables 2
  double ACPKKBfacts = adKK + DYKK; // B factories tau = 1
  double ACPpipiBfacts = adpipi + DYpipi;
  ll3 += meas.at("0709.2715_AcpKK_Babar").logweight(ACPKKBfacts);
  ll3 += meas.at("0709.2715_Acppipi_Babar").logweight(ACPpipiBfacts);


  // Acp(KK), Acp(pipi) Belle
  // https://arxiv.org/pdf/0807.0148
  // Observables 2
  ll3 += meas.at("0807.0148_AcpKK_Belle").logweight(ACPKKBfacts);
  ll3 += meas.at("0807.0148_Acppipi_Belle").logweight(ACPpipiBfacts);


  CKpi = -y12 * cos(PhiG12) * cos(dD_kpi) + x12 * cos(PhiM12) * sin(dD_kpi);
  CpKpi = 1. / 4. * (x12 * x12 + y12 * y12) + 0.25 * Rdp_uid30 * (y12 * y12 - x12 * x12);
  DCKpi = -y12 * sin(PhiG12) * sin(dD_kpi) - x12 * sin(PhiM12) * cos(dD_kpi);
  DCpKpi = 0.5 * x12 * y12 * sin(phi12);
  double AtildeKpi = - 2. * adKK;
  double DCtildeKpi = DCKpi - CKpi * adKK - 2. * rD_kpi * DYKK;
  double DCtildepKpi = DCpKpi - 2. * CpKpi * adKK - 2. * rD_kpi * CKpi * DYKK;

  corr.ResizeTo(9);
  corr(0) = Rdp_uid30;
  corr(1) = CKpi;
  corr(2) = CpKpi;
  corr(3) = AtildeKpi;
  corr(4) = DCtildeKpi;
  corr(5) = DCtildepKpi;
  corr(6) = AD;
  corr(7) = DCKpi;
  corr(8) = DCpKpi;
  ll3 += corrmeas.at("2407.18001").logweight(corr);


  corr.ResizeTo(6);
  corr(0) = Rdp_uid30;
  corr(1) = CKpi;
  corr(2) = CpKpi;
  corr(3) = AD;
  corr(4) = DCKpi;
  corr(5) = DCpKpi;
  ll3 += corrmeas.at("UID30").logweight(corr);


  corr.ResizeTo(2);
  corr(0) = Akpi_BESIII;
  corr(1) = Akpi_kpipi0_BESIII;
  ll3 += corrmeas.at("BESIII_Adk").logweight(corr);


  corr.ResizeTo(2);
  corr(0) = xi_x_BESIII;
  corr(1) = xi_y_BESIII;
  ll3 += corrmeas.at("BESIII_rDkpi_polar").logweight(corr);


  corr.ResizeTo(4);
  corr(0) = xcp_uid14;
  corr(1) = ycp_uid14;
  corr(2) = dx_uid14;
  corr(3) = dy_uid14;

  ll3 += corrmeas.at("UID14").logweight(corr);


  corr.ResizeTo(4);
  corr(0) = xcp;
  corr(1) = ycp;
  corr(2) = dx;
  corr(3) = dy;

  ll3 += corrmeas.at("LHCb_kspp_Au2022").logweight(corr);

  return ll3;
}
// ---------------------------------------------------------

// ---------------------------------------------------------
double MixingModel::Calculate_other_observables()
{

  double ll4;
  ll4 = 0.;

  // 20. PDF: dk3pi_dkpipi0_constraints (UID19)
  //  6 Observables
  kD_k3pi_uid19 = kD_k3pi;
  dD_k3pi_uid19 = dD_k3pi; // because of the phase sign convention
  kD_kpipi0_uid19 = kD_kpipi0;
  dD_kpipi0_uid19 = dD_kpipi0; // because of the phase sign convention
  rD_k3pi_uid19 = rD_k3pi;
  rD_kpipi0_uid19 = rD_kpipi0;

  // 21. PDF: d4pi_dmixing_cleo (UID20)
  //  1 Osservabile
  F_pipipipi_uid20 = F_pipipipi;

  // 22. PDF: CleoDhhpi0Dilution (UID21)
  F_pipipi0_uid21 = F_pipipi0;
  F_kkpi0_uid21 = F_kkpi0;

  //----------------------------------------------- Other Observables-------------------------------------------------------------------------
  // 14. PDF: charm-kspipi-nocpv (UID13)
  // 2 Observables
  x_uid13 = x;
  y_uid13 = y;



  // 23. PDF: dkskpiRWS (UID22)
  //  1 Osservabile
  RD_kskpi_uid22 = (rD_kskpi * rD_kskpi - rD_kskpi * kD_kskpi * (y * cos(dD_kskpi) - x * sin(dD_kskpi))) / (1. - rD_kskpi * kD_kskpi * (y * cos(dD_kskpi) + x * sin(dD_kskpi)));
  // Relative Sign Problem

  // 24. PDF: dkskpi (UID23)
  //  3 Observables
  RD_kskpi_uid23 = (rD_kskpi * rD_kskpi - rD_kskpi * kD_kskpi * (y * cos(dD_kskpi) - x * sin(dD_kskpi))) / (1. - rD_kskpi * kD_kskpi * (y * cos(dD_kskpi) + x * sin(dD_kskpi)));
  dD_kskpi_uid23 = dD_kskpi;
  kD_kskpi_uid23 = kD_kskpi;

  F_pipipipi_BESIII = F_pipipipi;

  //----------------------------------------------- Contribution to the LogLikelihood -------------------------------------------------------------------------
  TVectorD corr(8);

  // https://arxiv.org/pdf/2503.19542
  ll4 += meas.at("BESIII_2503.19542_BrDKpi").logweight( rD_kpi * rD_kpi + rD_kpi * y_plus(dD_kpi) + 0.5 * x_plus(dD_kpi) * x_plus(dD_kpi) * y_plus(dD_kpi) * y_plus(dD_kpi) );
  ll4 += meas.at("BESIII_2503.19542_BrDK3pi").logweight( rD_k3pi * rD_k3pi + kD_k3pi * rD_k3pi * y_plus(dD_k3pi) + 0.5 * x_plus(dD_k3pi) * x_plus(dD_k3pi) * y_plus(dD_k3pi) * y_plus(dD_k3pi) );
  ll4 += meas.at("BESIII_2503.19542_BrDKpipi0").logweight( rD_kpipi0 * rD_kpipi0 + kD_kpipi0 * rD_kpipi0 * y_plus(dD_kpipi0) + 0.5 * x_plus(dD_kpipi0) * x_plus(dD_kpipi0) * y_plus(dD_kpipi0) * y_plus(dD_kpipi0) );



  corr.ResizeTo(2);
  corr(0) = F_pipipi0;
  corr(1) = F_kkpi0;
  ll4 += corrmeas.at("UID21").logweight(corr);


  corr.ResizeTo(2);
  corr(0) = F_pipipi0;
  corr(1) = F_kkpi0;
  ll4 += corrmeas.at("2409.07197_F_BESIII").logweight(corr);


  ll4 += meas.at("UID20").logweight(F_pipipipi_uid20);


  ll4 += meas.at("Fpipipipi_BESIII").logweight(F_pipipipi_BESIII);

  ll4 += meas.at("FKKpipi_BESIII").logweight(F_kkpipi);


  corr.ResizeTo(6);
  corr(0) = kD_k3pi_uid19;
  corr(1) = dD_k3pi_uid19;
  corr(2) = kD_kpipi0_uid19;
  corr(3) = dD_kpipi0_uid19;
  corr(4) = rD_k3pi_uid19;
  corr(5) = rD_kpipi0_uid19;
  ll4 += corrmeas.at("UID19").logweight(corr);


  ll4 += meas.at("UID22").logweight(RD_kskpi_uid22);


  corr.ResizeTo(3);
  corr(0) = RD_kskpi_uid23;
  corr(1) = dD_kskpi_uid23;
  corr(2) = kD_kskpi_uid23;

  ll4 += corrmeas.at("UID23").logweight(corr);


  return ll4;
}
// ---------------------------------------------------------

// ---------------------------------------------------------
double MixingModel::Calculate_old_observables()
{

  double llo;
  llo = 0.;

  Rd = rD_kpi * rD_kpi;

  // 6th Block
  rm = (x * x + y * y) / 2.;

  yp_kpp_plus = y_plus(dD_kpipi0);
  xp_kpp_plus = x_plus(dD_kpipi0);

  yp_kpp_minus = y_minus(dD_kpipi0);
  xp_kpp_minus = x_minus(dD_kpipi0);

  yp_plus = y_plus(dD_kpi);
  xp_plus = x_plus(dD_kpi);
  xp_plus_sq = xp_plus * xp_plus;

  yp_minus = y_minus(dD_kpi);
  xp_minus = x_minus(dD_kpi);
  xp_minus_sq = xp_minus * xp_minus;

  TVectorD corr(4);


  corr.ResizeTo(3);
  corr(0) = Rd;
  corr(1) = xp_plus_sq;
  corr(2) = yp_plus;
  llo += corrmeas.at("kpi_babar_plus").logweight(corr);
  llo += corrmeas.at("kpi_belle_plus").logweight(corr);

  corr(0) = AD;
  corr(1) = xp_minus_sq;
  corr(2) = yp_minus;
  llo += corrmeas.at("kpi_babar_minus").logweight(corr);
  llo += corrmeas.at("kpi_belle_minus").logweight(corr);



  corr.ResizeTo(5);
  corr(0) = Rd;
  corr(1) = x * x;
  corr(2) = y;
  corr(3) = cos(M_PI - dD_kpi);
  corr(4) = sin(M_PI - dD_kpi);
  llo += corrmeas.at("cleoc").logweight(corr);

  // 3rd Block
  double epsI = 2.228 * sin(43.5 * M_PI / 180.) * 1.e-3; // values taken from PDG: https://pdglive.lbl.gov/ParticleGroup.action?init=0&node=MXXX020
  double RCKM = 0.00384 * 0.04120 / 0.2251 / 0.97345 * sin(g); // values taken from https://indico.cern.ch/event/1291157/contributions/5903548/attachments/2900988/5087304/bona-utfit.pdf
  corr.ResizeTo(4);
  corr(0) = x;
  corr(1) = y;
  corr(2) = qop;
  corr(3) = phi - 2*epsI - RCKM; // Because of B-factories
  llo += corrmeas.at("kpp_belle").logweight(corr);


  corr.ResizeTo(2);
  corr(0) = x;
  corr(1) = y;

  llo += corrmeas.at("kppkk").logweight(corr);
  llo += corrmeas.at("kppp0_Babar").logweight(corr);
  llo += corrmeas.at("UID13").logweight(corr);


  corr.ResizeTo(2);

  corr(0) = xp_kpp_plus;
  corr(1) = yp_kpp_plus;
  llo += corrmeas.at("kpp_babar_plus").logweight(corr);

  corr(0) = xp_kpp_minus;
  corr(1) = yp_kpp_minus;
  llo += corrmeas.at("kpp_babar_minus").logweight(corr);


  k3pi_uid18 = 0.25 * (x * x + y * y);
  llo += meas.at("UID18").logweight(k3pi_uid18);


  llo += meas.at("RM").logweight(rm);

  return llo;
}
// ---------------------------------------------------------

double MixingModel::LogLikelihood(const std::vector<double> &parameters)
{

  double ll = 0.;

  if (comb == 0)
  {

    // Variabili globali
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

    // 5. PDF: glwads-dkdpi-kskpi-dmix (UID4)
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

    // 10. PDF: glwads-dhpipi-hh-dmix (UID9)
    r_dkpipi = parameters[29];
    d_dkpipi = parameters[30];
    k_dkpipi = parameters[31];
    r_dpipipi = parameters[32];
    d_dpipipi = parameters[33];
    k_dpipipi = parameters[34];

    // 15. PDF: charm-kspipi (UID14)
    PhiM12 = parameters[35];
    PhiG12 = parameters[36];

    adKK = parameters[37];
    adpipi = parameters[38];
    DYKKmDYpipi = parameters[39];
    tavepitaggedOverTauD = parameters[40];
    tavemutaggedOverTauD = parameters[41];
    DeltatmutaggedOverTauD = parameters[42];
    DeltatpitaggedOverTauD = parameters[43];
    tKKCDp = parameters[44];
    tKKCDs = parameters[45];

    F_kkpipi = parameters[46];

    tauKK_DAcp_Run1_sl = parameters[47];
    taupipi_DAcp_Run1_sl = parameters[48];
    tauKK_Acp_Run1_sl = parameters[49];

    tauKK_DAcp_Run1_pi = parameters[50];
    taupipi_DAcp_Run1_pi = parameters[51];
    tauKK_Acp_Run1_pi = parameters[52];

    tauKK_Acp_CDF = parameters[53];
    taupipi_Acp_CDF = parameters[54];

    // General parameters
    AD = 0.; // NO direct CPV for CF/DCS
    phi12 = remainder(-PhiG12 + PhiM12, 2. * M_PI);
    x = x12;
    y = y12;
    qop = 1. + x12 * y12 * sin(phi12) / (x12 * x12 + y12 * y12);
    phi = atan(-(x12 * x12 * sin(2. * PhiM12) + y12 * y12 * sin(2. * PhiG12)) / (x12 * x12 * cos(2. * PhiM12) + y12 * y12 * cos(2. * PhiG12))) / 2.;
    d = (1. - qop * qop) / (1. + qop * qop);

    xcp = 0.5 * (x * cos(phi) * (qop + 1. / qop) + y * sin(phi) * (qop - 1. / qop));
    ycp = 0.5 * (y * cos(phi) * (qop + 1. / qop) - x * sin(phi) * (qop - 1. / qop));
    dx = 0.5 * (x * cos(phi) * (qop - 1. / qop) + y * sin(phi) * (qop + 1. / qop));
    dy = 0.5 * (y * cos(phi) * (qop - 1. / qop) - x * sin(phi) * (qop + 1. / qop));

    //---------------------------------------------------- Contribution to the LogLikelihood of the charged B measurements  ----------------------------------------------------
    ll += Calculate_ChargedB_observables();
    //------------------------------------------------------  Contribution to the LogLikelihood of the Time Dependent D observables ----------------------------------------------------
    ll += Calculate_time_dependent_Dobservables();
    //----------------------------------------------- Contribution to the LogLikelihood of the Other Observables -------------------------------------------------------------------------
    ll += Calculate_other_observables();
    //-----------------------------------------------  Contribution to the LogLikelihood of the Old Observables -------------------------------------------------------------------------
    ll += Calculate_old_observables();

    //---------------------------------- Fill the Histograms ---------------------

    // Time Integrated B decay chain
    obs["g"] = g * r2d;
    obs["x12"] = x12 * 1000;
    obs["y12"] = y12 * 1000;
    obs["r_dk"] = r_dk * 100;
    obs["r_dpi"] = r_dpi * 1000;
    obs["rD_kpi"] = rD_kpi * 100;
    obs["d_dk"] = d_dk * r2d;
    obs["d_dpi"] = d_dpi * r2d;
    obs["dD_kpi"] = dD_kpi * r2d;
    obs["rD_k3pi"] = rD_k3pi;
    obs["dD_k3pi"] = dD_k3pi * r2d;
    obs["kD_k3pi"] = kD_k3pi;
    obs["F_pipipipi"] = F_pipipipi;
    obs["rD_kpipi0"] = rD_kpipi0;
    obs["dD_kpipi0"] = dD_kpipi0 * r2d;
    obs["kD_kpipi0"] = kD_kpipi0;
    obs["F_pipipi0"] = F_pipipi0;
    obs["F_kkpi0"] = F_kkpi0;
    obs["rD_kskpi"] = rD_kskpi;
    obs["dD_kskpi"] = dD_kskpi * r2d;
    obs["kD_kskpi"] = kD_kskpi;
    obs["RBRdkdpi"] = RBRdkdpi;
    obs["r_dstk"] = r_dstk;
    obs["d_dstk"] = d_dstk * r2d;
    obs["r_dstpi"] = r_dstpi;
    obs["d_dstpi"] = d_dstpi * r2d;
    obs["r_dkst"] = r_dkst;
    obs["d_dkst"] = d_dkst * r2d;
    obs["k_dkst"] = k_dkst;
    obs["r_dkpipi"] = r_dkpipi;
    obs["d_dkpipi"] = d_dkpipi * r2d;
    obs["k_dkpipi"] = k_dkpipi;
    obs["r_dpipipi"] = r_dpipipi;
    obs["d_dpipipi"] = d_dpipipi * r2d;
    obs["k_dpipipi"] = k_dpipipi;
    obs["F_kkpipi"] = F_kkpipi;

    // Time dependent D decay and mixing
    obs["PhiM12"] = PhiM12 * r2d;
    obs["PhiG12"] = PhiG12 * r2d;
    obs["phipphig12"] = remainder(PhiG12 + phi, 2. * M_PI) * r2d;
    obs["phimphig12"] = remainder(-PhiG12 + phi, 2. * M_PI) * r2d;
    obs["AD"] = AD;
    obs["adKK"] = adKK * 1000;
    obs["adpipi"] = adpipi * 1000;
    obs["qopm1"] = (qop - 1) * 100;
    obs["qop"] = qop;
    obs["phi"] = phi * r2d;
    obs["phi12"] = phi12 * r2d;
    obs["delta"] = d;

    obs["x"] = x * 1000;
    obs["y"] = y * 1000;
    obs["M12"] = 0.5 * x12 / tau;
    obs["G12"] = y12 / tau;
    obs["ImM12"] = 0.5 * x12 / tau * sin(PhiM12);
  }
  else if (comb == 1)
  {

    // Variabili globali
    g = parameters[0];
    x12 = parameters[1];
    y12 = parameters[2];

    // 1. PDF: glwads-dh-hh-dmix (UID0)
    rD_kpi = parameters[3];
    dD_kpi = parameters[4];

    //  2. PDF: glwads-dh-h3pi-dmix (UID1)
    rD_k3pi = parameters[5];
    dD_k3pi = parameters[6];
    kD_k3pi = parameters[7];
    F_pipipipi = parameters[8];

    // 3. PDF: glwads-dh-hhpi0-dmix (UID2)
    rD_kpipi0 = parameters[9];
    dD_kpipi0 = parameters[10];
    kD_kpipi0 = parameters[11];
    F_pipipi0 = parameters[12];
    F_kkpi0 = parameters[13];

    // 5. PDF: glwads-dkdpi-kskpi-dmix (UID4)
    rD_kskpi = parameters[14];
    dD_kskpi = parameters[15];
    kD_kskpi = parameters[16];

    //  8. PDF: glwads-dkstz-hh-h3pi-dmix
    r_dkstz = parameters[17];
    d_dkstz = parameters[18];
    k_dkstz = parameters[19];

    // 13. PDF: dmpi (UID12)
    l_dmpi = parameters[20];
    d_dmpi = parameters[21];
    phi_d = parameters[22];

    // 15. PDF: charm-kspipi (UID14)
    PhiM12 = parameters[23];
    PhiG12 = parameters[24];

    adKK = parameters[25];
    adpipi = parameters[26];
    DYKKmDYpipi = parameters[27];
    tavepitaggedOverTauD = parameters[28];
    tavemutaggedOverTauD = parameters[29];
    DeltatmutaggedOverTauD = parameters[30];
    DeltatpitaggedOverTauD = parameters[31];
    tKKCDp = parameters[32];
    tKKCDs = parameters[33];

    F_kkpipi = parameters[34];

    tauKK_DAcp_Run1_sl = parameters[35];
    taupipi_DAcp_Run1_sl = parameters[36];
    tauKK_Acp_Run1_sl = parameters[37];

    tauKK_DAcp_Run1_pi = parameters[38];
    taupipi_DAcp_Run1_pi = parameters[39];
    tauKK_Acp_Run1_pi = parameters[40];

    tauKK_Acp_CDF = parameters[41];
    taupipi_Acp_CDF = parameters[42];

    l_dstarmpi = parameters[43];
    d_dstarmpi = parameters[44];
    l_dmrho = parameters[45];
    d_dmrho = parameters[46];

    // General parameters
    AD = 0.; // NO direct CPV for CF/DCS
    phi12 = remainder(-PhiG12 + PhiM12, 2. * M_PI);
    x = x12;
    y = y12;
    qop = 1. + x12 * y12 * sin(phi12) / (x12 * x12 + y12 * y12);
    phi = atan(-(x12 * x12 * sin(2. * PhiM12) + y12 * y12 * sin(2. * PhiG12)) / (x12 * x12 * cos(2. * PhiM12) + y12 * y12 * cos(2. * PhiG12))) / 2.;
    d = (1. - qop * qop) / (1. + qop * qop);

    xcp = 0.5 * (x * cos(phi) * (qop + 1. / qop) + y * sin(phi) * (qop - 1. / qop));
    ycp = 0.5 * (y * cos(phi) * (qop + 1. / qop) - x * sin(phi) * (qop - 1. / qop));
    dx = 0.5 * (x * cos(phi) * (qop - 1. / qop) + y * sin(phi) * (qop + 1. / qop));
    dy = 0.5 * (y * cos(phi) * (qop - 1. / qop) - x * sin(phi) * (qop + 1. / qop));

    //------------------------------------------------------  Contribution to the LogLikelihood of the Bd observables  ----------------------------------------------------
    ll += Calculate_neutralBdobservables();
    //------------------------------------------------------  Contribution to the LogLikelihood of the Time Dependent D observables ----------------------------------------------------
    ll += Calculate_time_dependent_Dobservables();
    //----------------------------------------------- Contribution to the LogLikelihood of the Other Observables -------------------------------------------------------------------------
    ll += Calculate_other_observables();
    //----------------------------------------------- Contribution to the LogLikelihood of the Old Observables -------------------------------------------------------------------------
    ll += Calculate_old_observables();

    //---------------------------------- Fill the Histograms ---------------------

    // Time Integrated B decay chain
    obs["g"] = g * r2d;
    obs["x12"] = x12 * 1000;
    obs["y12"] = y12 * 1000;
    obs["rD_kpi"] = rD_kpi * 100;
    obs["dD_kpi"] = dD_kpi * r2d;
    obs["rD_k3pi"] = rD_k3pi;
    obs["dD_k3pi"] = dD_k3pi * r2d;
    obs["kD_k3pi"] = kD_k3pi;
    obs["F_pipipipi"] = F_pipipipi;
    obs["rD_kpipi0"] = rD_kpipi0;
    obs["dD_kpipi0"] = dD_kpipi0 * r2d;
    obs["kD_kpipi0"] = kD_kpipi0;
    obs["F_pipipi0"] = F_pipipi0;
    obs["F_kkpi0"] = F_kkpi0;
    obs["rD_kskpi"] = rD_kskpi;
    obs["dD_kskpi"] = dD_kskpi * r2d;
    obs["kD_kskpi"] = kD_kskpi;
    obs["r_dkstz"] = r_dkstz;
    obs["d_dkstz"] = d_dkstz * r2d;
    obs["k_dkstz"] = k_dkstz;

    // Time dependent B decay chain
    obs["l_dmpi"] = l_dmpi;
    obs["d_dmpi"] = d_dmpi * r2d;
    obs["beta"] = phi_d*0.5 * r2d;
    obs["phid"] = phi_d * r2d;

    // Time dependent D decay and mixing
    obs["PhiM12"] = PhiM12 * r2d;
    obs["PhiG12"] = PhiG12 * r2d;
    obs["phipphig12"] = remainder(PhiG12 + phi, 2. * M_PI) * r2d;
    obs["phimphig12"] = remainder(-PhiG12 + phi, 2. * M_PI) * r2d;
    obs["AD"] = AD;
    obs["adKK"] = adKK * 1000;
    obs["adpipi"] = adpipi * 1000;
    obs["qopm1"] = (qop - 1) * 100;
    obs["qop"] = qop;
    obs["phi"] = phi * r2d;
    obs["phi12"] = phi12 * r2d;
    obs["delta"] = d;

    obs["x"] = x * 1000;
    obs["y"] = y * 1000;
    obs["M12"] = 0.5 * x12 / tau ;
    obs["G12"] = y12 / tau ;
    obs["ImM12"] = 0.5 * x12 / tau * sin(PhiM12);

  }
  else if (comb == 2)
  {

    // Variabili globali
    g = parameters[0];
    x12 = parameters[1];
    y12 = parameters[2];

    // 1. PDF: glwads-dh-hh-dmix (UID0)
    rD_kpi = parameters[3];
    dD_kpi = parameters[4];

    //  2. PDF: glwads-dh-h3pi-dmix (UID1)
    rD_k3pi = parameters[5];
    dD_k3pi = parameters[6];
    kD_k3pi = parameters[7];
    F_pipipipi = parameters[8];

    // 3. PDF: glwads-dh-hhpi0-dmix (UID2)
    rD_kpipi0 = parameters[9];
    dD_kpipi0 = parameters[10];
    kD_kpipi0 = parameters[11];
    F_pipipi0 = parameters[12];
    F_kkpi0 = parameters[13];

    // 5. PDF: glwads-dkdpi-kskpi-dmix (UID4)
    rD_kskpi = parameters[14];
    dD_kskpi = parameters[15];
    kD_kskpi = parameters[16];

    // 11. PDF: dsk (UID10)
    l_dsk = parameters[17];
    d_dsk = parameters[18];
    phis = parameters[19];

    // 12. PDF: dskpipi (UID11)
    l_dskpipi = parameters[20];
    d_dskpipi = parameters[21];
    k_dskpipi = parameters[22];

    // 15. PDF: charm-kspipi (UID14)
    PhiM12 = parameters[23];
    PhiG12 = parameters[24];

    adKK = parameters[25];
    adpipi = parameters[26];
    DYKKmDYpipi = parameters[27];
    tavepitaggedOverTauD = parameters[28];
    tavemutaggedOverTauD = parameters[29];
    DeltatmutaggedOverTauD = parameters[30];
    DeltatpitaggedOverTauD = parameters[31];
    tKKCDp = parameters[32];
    tKKCDs = parameters[33];

    F_kkpipi = parameters[34];

    // 2401.17934 Bs part
    r_dkstzs = parameters[35];
    d_dkstzs = parameters[36];
    k_dkstzs = parameters[37];

    tauKK_DAcp_Run1_sl = parameters[38];
    taupipi_DAcp_Run1_sl = parameters[39];
    tauKK_Acp_Run1_sl = parameters[40];

    tauKK_DAcp_Run1_pi = parameters[41];
    taupipi_DAcp_Run1_pi = parameters[42];
    tauKK_Acp_Run1_pi = parameters[43];

    tauKK_Acp_CDF = parameters[44];
    taupipi_Acp_CDF = parameters[45];

    // General parameters
    AD = 0.; // NO direct CPV for CF/DCS
    phi12 = remainder(-PhiG12 + PhiM12, 2. * M_PI);
    x = x12;
    y = y12;
    qop = 1. + x12 * y12 * sin(phi12) / (x12 * x12 + y12 * y12);
    phi = atan(-(x12 * x12 * sin(2. * PhiM12) + y12 * y12 * sin(2. * PhiG12)) / (x12 * x12 * cos(2. * PhiM12) + y12 * y12 * cos(2. * PhiG12))) / 2.;
    d = (1. - qop * qop) / (1. + qop * qop);

    xcp = 0.5 * (x * cos(phi) * (qop + 1. / qop) + y * sin(phi) * (qop - 1. / qop));
    ycp = 0.5 * (y * cos(phi) * (qop + 1. / qop) - x * sin(phi) * (qop - 1. / qop));
    dx = 0.5 * (x * cos(phi) * (qop - 1. / qop) + y * sin(phi) * (qop + 1. / qop));
    dy = 0.5 * (y * cos(phi) * (qop - 1. / qop) - x * sin(phi) * (qop + 1. / qop));

    //------------------------------------------------------  Contribution to the LogLikelihood of the B0s observables  ----------------------------------------------------
    ll += Calculate_neutralBsobservables();
    //------------------------------------------------------  Contribution to the LogLikelihood of the Time Dependent D observables  ----------------------------------------------------
    ll += Calculate_time_dependent_Dobservables();
    //----------------------------------------------- Contribution to the LogLikelihood of the Other Observables -------------------------------------------------------------------------
    ll += Calculate_other_observables();
    //----------------------------------------------- Contribution to the LogLikelihood of the Old Observables -------------------------------------------------------------------------
    ll += Calculate_old_observables();

    //---------------------------------- Fill the Histograms ---------------------

    // Time Integrated B decay chain
    obs["g"] = g * r2d;
    obs["x12"] = x12 * 1000;
    obs["y12"] = y12 * 1000;
    obs["rD_kpi"] = rD_kpi * 100;
    obs["dD_kpi"] = dD_kpi * r2d;
    obs["rD_k3pi"] = rD_k3pi;
    obs["dD_k3pi"] = dD_k3pi * r2d;
    obs["kD_k3pi"] = kD_k3pi;
    obs["F_pipipipi"] = F_pipipipi;
    obs["rD_kpipi0"] = rD_kpipi0;
    obs["dD_kpipi0"] = dD_kpipi0 * r2d;
    obs["kD_kpipi0"] = kD_kpipi0;
    obs["F_pipipi0"] = F_pipipi0;
    obs["F_kkpi0"] = F_kkpi0;
    obs["rD_kskpi"] = rD_kskpi;
    obs["dD_kskpi"] = dD_kskpi * r2d;
    obs["kD_kskpi"] = kD_kskpi;

    // Time integrated Bs chain
    obs["r_dkstzs"] = r_dkstzs;
    obs["d_dkstzs"] = d_dkstzs * r2d;
    obs["k_dkstzs"] = k_dkstzs;

    // Time dependent B decay chain

    obs["l_dsk"] = l_dsk;
    obs["d_dsk"] = d_dsk * r2d;
    obs["phis"] = phis * r2d;
    obs["beta_s"] = -0.5 * phis;
    obs["l_dskpipi"] = l_dskpipi;
    obs["d_dskpipi"] = d_dskpipi * r2d;
    obs["k_dskpipi"] = k_dskpipi;

    // Time dependent D decay and mixing
    obs["PhiM12"] = PhiM12 * r2d;
    obs["PhiG12"] = PhiG12 * r2d;
    obs["phipphig12"] = remainder(PhiG12 + phi, 2. * M_PI) * r2d;
    obs["phimphig12"] = remainder(-PhiG12 + phi, 2. * M_PI) * r2d;
    obs["AD"] = AD;
    obs["adKK"] = adKK * 1000;
    obs["adpipi"] = adpipi * 1000;
    obs["qopm1"] = (qop - 1) * 100;
    obs["qop"] = qop;
    obs["phi"] = phi * r2d;
    obs["phi12"] = phi12 * r2d;
    obs["delta"] = d;

    obs["x"] = x * 1000;
    obs["y"] = y * 1000;
    obs["M12"] = 0.5 * x12 / tau ;
    obs["G12"] = y12 / tau ;
    obs["ImM12"] = 0.5 * x12 / tau * sin(PhiM12);
  }
  else if (comb == 3)
  {

    // Variabili globali
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

    // 5. PDF: glwads-dkdpi-kskpi-dmix (UID4)
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

    //  8. PDF: glwads-dkstz-hh-h3pi-dmix
    r_dkstz = parameters[29];
    d_dkstz = parameters[30];
    k_dkstz = parameters[31];

    // 10. PDF: glwads-dhpipi-hh-dmix (UID9)
    r_dkpipi = parameters[32];
    d_dkpipi = parameters[33];
    k_dkpipi = parameters[34];
    r_dpipipi = parameters[35];
    d_dpipipi = parameters[36];
    k_dpipipi = parameters[37];

    // 11. PDF: dsk (UID10)
    l_dsk = parameters[38];
    d_dsk = parameters[39];
    phis = parameters[40];

    // 12. PDF: dskpipi (UID11)
    l_dskpipi = parameters[41];
    d_dskpipi = parameters[42];
    k_dskpipi = parameters[43];

    // 13. PDF: dmpi (UID12)
    l_dmpi = parameters[44];
    d_dmpi = parameters[45];
    phi_d = parameters[46];

    // 15. PDF: charm-kspipi (UID14)
    PhiM12 = parameters[47];
    PhiG12 = parameters[48];

    adKK = parameters[49];
    adpipi = parameters[50];
    DYKKmDYpipi = parameters[51];
    tavepitaggedOverTauD = parameters[52];
    tavemutaggedOverTauD = parameters[53];
    DeltatmutaggedOverTauD = parameters[54];
    DeltatpitaggedOverTauD = parameters[55];
    tKKCDp = parameters[56];
    tKKCDs = parameters[57];

    // https://arxiv.org/abs/2301.10328
    F_kkpipi = parameters[58];

    // 2401.17934 Bs part
    r_dkstzs = parameters[59];
    d_dkstzs = parameters[60];
    k_dkstzs = parameters[61];

    tauKK_DAcp_Run1_sl = parameters[62];
    taupipi_DAcp_Run1_sl = parameters[63];
    tauKK_Acp_Run1_sl = parameters[64];

    tauKK_DAcp_Run1_pi = parameters[65];
    taupipi_DAcp_Run1_pi = parameters[66];
    tauKK_Acp_Run1_pi = parameters[67];

    tauKK_Acp_CDF = parameters[68];
    taupipi_Acp_CDF = parameters[69];

    l_dstarmpi = parameters[70];
    d_dstarmpi = parameters[71];
    l_dmrho = parameters[72];
    d_dmrho = parameters[73];


    // General parameters
    AD = 0.; // NO direct CPV for CF/DCS
    phi12 = remainder(-PhiG12 + PhiM12, 2. * M_PI);
    x = x12;
    y = y12;
    qop = 1. + x12 * y12 * sin(phi12) / (x12 * x12 + y12 * y12);
    phi = atan(-(x12 * x12 * sin(2. * PhiM12) + y12 * y12 * sin(2. * PhiG12)) / (x12 * x12 * cos(2. * PhiM12) + y12 * y12 * cos(2. * PhiG12))) / 2.;
    d = (1. - qop * qop) / (1. + qop * qop);

    xcp = 0.5 * (x * cos(phi) * (qop + 1. / qop) + y * sin(phi) * (qop - 1. / qop));
    ycp = 0.5 * (y * cos(phi) * (qop + 1. / qop) - x * sin(phi) * (qop - 1. / qop));
    dx = 0.5 * (x * cos(phi) * (qop - 1. / qop) + y * sin(phi) * (qop + 1. / qop));
    dy = 0.5 * (y * cos(phi) * (qop - 1. / qop) - x * sin(phi) * (qop + 1. / qop));

    //---------------------------------------------------- Contribution to the LogLikelihood of the charged B observables ----------------------------------------------------
    ll += Calculate_ChargedB_observables();
    //------------------------------------------------------  Contribution to the LogLikelihood of the neutral Bd observables  ----------------------------------------------------
    ll += Calculate_neutralBdobservables();
    //------------------------------------------------------  Contribution to the LogLikelihood of the neutral Bs observables  ----------------------------------------------------
    ll += Calculate_neutralBsobservables();
    //------------------------------------------------------  Contribution to the LogLikelihood of the Time Dependent D observables  ----------------------------------------------------
    ll += Calculate_time_dependent_Dobservables();
    //----------------------------------------------- Contribution to the LogLikelihood of the Other Observables -------------------------------------------------------------------------
    ll += Calculate_other_observables();
    //----------------------------------------------- Contribution to the LogLikelihood of the Old Observables -------------------------------------------------------------------------
    ll += Calculate_old_observables();

    //---------------------------------- Fill the map "obs", a subset of which will be used to fill the histograms --------------------------------------------------------
    // Time Integrated B decay chain
    obs["g"] = g * r2d;
    obs["x12"] = x12 * 1000;
    obs["y12"] = y12 * 1000;
    obs["r_dk"] = r_dk * 100;
    obs["r_dpi"] = r_dpi * 1000;
    obs["rD_kpi"] = rD_kpi * 100;
    obs["d_dk"] = d_dk * r2d;
    obs["d_dpi"] = d_dpi * r2d;
    obs["dD_kpi"] = dD_kpi * r2d;
    obs["rD_k3pi"] = rD_k3pi;
    obs["dD_k3pi"] = dD_k3pi * r2d;
    obs["kD_k3pi"] = kD_k3pi;
    obs["F_pipipipi"] = F_pipipipi;
    obs["rD_kpipi0"] = rD_kpipi0;
    obs["dD_kpipi0"] = dD_kpipi0 * r2d;
    obs["kD_kpipi0"] = kD_kpipi0;
    obs["F_pipipi0"] = F_pipipi0;
    obs["F_kkpi0"] = F_kkpi0;
    obs["rD_kskpi"] = rD_kskpi;
    obs["dD_kskpi"] = dD_kskpi * r2d;
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
    obs["r_dkstzs"] = r_dkstzs;
    obs["d_dkstzs"] = d_dkstzs * r2d;
    obs["k_dkstzs"] = k_dkstzs;
    obs["r_dkpipi"] = r_dkpipi;
    obs["d_dkpipi"] = d_dkpipi * r2d;
    obs["k_dkpipi"] = k_dkpipi;
    obs["r_dpipipi"] = r_dpipipi;
    obs["d_dpipipi"] = d_dpipipi * r2d;
    obs["k_dpipipi"] = k_dpipipi;
    obs["F_kkpipi"] = F_kkpipi;

    // Time dependent B decay chain

    obs["l_dsk"] = l_dsk;
    obs["d_dsk"] = d_dsk * r2d;
    obs["phis"] = phis * r2d;
    obs["beta_s"] = -0.5 * phis;
    obs["l_dskpipi"] = l_dskpipi;
    obs["d_dskpipi"] = d_dskpipi * r2d;
    obs["k_dskpipi"] = k_dskpipi;
    obs["l_dmpi"] = l_dmpi;
    obs["d_dmpi"] = d_dmpi * r2d;
    obs["beta"] = phi_d*0.5 * r2d;
    obs["phid"] = phi_d * r2d;

    // Time dependent D decay and mixing
    obs["PhiM12"] = PhiM12 * r2d;
    obs["PhiG12"] = PhiG12 * r2d;
    obs["phipphig12"] = remainder(PhiG12 + phi, 2. * M_PI) * r2d;
    obs["phimphig12"] = remainder(-PhiG12 + phi, 2. * M_PI) * r2d;
    obs["AD"] = AD;
    obs["adKK"] = adKK * 1000;
    obs["adpipi"] = adpipi * 1000;
    obs["qopm1"] = (qop - 1) * 100;
    obs["qop"] = qop;
    obs["phi"] = phi * r2d;
    obs["phi12"] = phi12 * r2d;
    obs["delta"] = d;

    obs["x"] = x * 1000;
    obs["y"] = y * 1000;
    obs["M12"] = 0.5 * x12 / tau ; // ps^-1
    obs["G12"] = y12 / tau ;
    obs["ImM12"] = 0.5 * x12 / tau * sin(PhiM12);


    // Other parameters
    obs["DYKKmDYpipi"] = DYKKmDYpipi*1000; // permille
    obs["tavepitaggedOverTauD"] = tavepitaggedOverTauD;
    obs["tavemutaggedOverTauD"] = tavemutaggedOverTauD;
    obs["DeltatmutaggedOverTauD"] = DeltatmutaggedOverTauD;
    obs["DeltatpitaggedOverTauD"] = DeltatpitaggedOverTauD;
    obs["tKKCDp"] = tKKCDp*1e12;
    obs["tKKCDs"] = tKKCDs*1e12;
    obs["tauKK_DAcp_Run1_sl"] = tauKK_DAcp_Run1_sl;
    obs["taupipi_DAcp_Run1_sl"] = taupipi_DAcp_Run1_sl;
    obs["tauKK_Acp_Run1_sl"] = tauKK_Acp_Run1_sl;
    obs["tauKK_DAcp_Run1_pi"] = tauKK_DAcp_Run1_pi;
    obs["taupipi_DAcp_Run1_pi"] = taupipi_DAcp_Run1_pi;
    obs["tauKK_Acp_Run1_pi"] = tauKK_Acp_Run1_pi;
    obs["tauKK_Acp_CDF"] = tauKK_Acp_CDF;
    obs["taupipi_Acp_CDF"] = taupipi_Acp_CDF;
    obs["l_dstarmpi"] = l_dstarmpi;
    obs["d_dstarmpi"] = d_dstarmpi*180./M_PI;
    obs["l_dmrho"] = l_dmrho;
    obs["d_dmrho"] = d_dmrho*180./M_PI;

  }
  else if (comb == 4)
  {

    // Variabili globali
    x12 = parameters[0];
    y12 = parameters[1];

    // 1. PDF: glwads-dh-hh-dmix (UID0)
    rD_kpi = parameters[2];
    dD_kpi = parameters[3];

    //  2. PDF: glwads-dh-h3pi-dmix (UID1)
    rD_k3pi = parameters[4];
    dD_k3pi = parameters[5];
    kD_k3pi = parameters[6];
    F_pipipipi = parameters[7];

    // 3. PDF: glwads-dh-hhpi0-dmix (UID2)
    rD_kpipi0 = parameters[8];
    dD_kpipi0 = parameters[9];
    kD_kpipi0 = parameters[10];
    F_pipipi0 = parameters[11];
    F_kkpi0 = parameters[12];

    // 5. PDF: glwads-dkdpi-kskpi-dmix (UID4)
    rD_kskpi = parameters[13];
    dD_kskpi = parameters[14];
    kD_kskpi = parameters[15];

    // 15. PDF: charm-kspipi (UID14)
    PhiM12 = parameters[16];
    PhiG12 = parameters[17];

    adKK = parameters[18];
    adpipi = parameters[19];
    DYKKmDYpipi = parameters[20];
    tavepitaggedOverTauD = parameters[21];
    tavemutaggedOverTauD = parameters[22];
    DeltatmutaggedOverTauD = parameters[23];
    DeltatpitaggedOverTauD = parameters[24];
    tKKCDp = parameters[25];
    tKKCDs = parameters[26];

    F_kkpipi = parameters[27];

    tauKK_DAcp_Run1_sl = parameters[28];
    taupipi_DAcp_Run1_sl = parameters[29];
    tauKK_Acp_Run1_sl = parameters[30];

    tauKK_DAcp_Run1_pi = parameters[31];
    taupipi_DAcp_Run1_pi = parameters[32];
    tauKK_Acp_Run1_pi = parameters[33];

    tauKK_Acp_CDF = parameters[34];
    taupipi_Acp_CDF = parameters[35];


    // General parameters
    AD = 0.; // NO direct CPV for CF/DCS
    phi12 = remainder(-PhiG12 + PhiM12, 2. * M_PI);
    x = x12;
    y = y12;
    qop = 1. + x12 * y12 * sin(phi12) / (x12 * x12 + y12 * y12);
    phi = atan(-(x12 * x12 * sin(2. * PhiM12) + y12 * y12 * sin(2. * PhiG12)) / (x12 * x12 * cos(2. * PhiM12) + y12 * y12 * cos(2. * PhiG12))) / 2.;
    d = (1. - qop * qop) / (1. + qop * qop);

    xcp = 0.5 * (x * cos(phi) * (qop + 1. / qop) + y * sin(phi) * (qop - 1. / qop));
    ycp = 0.5 * (y * cos(phi) * (qop + 1. / qop) - x * sin(phi) * (qop - 1. / qop));
    dx = 0.5 * (x * cos(phi) * (qop - 1. / qop) + y * sin(phi) * (qop + 1. / qop));
    dy = 0.5 * (y * cos(phi) * (qop - 1. / qop) - x * sin(phi) * (qop + 1. / qop));

    //------------------------------------------------------  Contribution to the LogLikelihood of the Time Dependent D observables  ----------------------------------------------------
    ll += Calculate_time_dependent_Dobservables();
    //----------------------------------------------- Contribution to the LogLikelihood of the Other Observables -------------------------------------------------------------------------
    ll += Calculate_other_observables();
    //----------------------------------------------- Contribution to the LogLikelihood of the Old Observables -------------------------------------------------------------------------
    ll += Calculate_old_observables();

    //---------------------------------- Fill the Histograms ---------------------

    // Time Integrated B decay chain
    obs["x12"] = x12 * 1000;
    obs["y12"] = y12 * 1000;
    obs["rD_kpi"] = rD_kpi * 100;
    obs["dD_kpi"] = dD_kpi * r2d;
    obs["rD_k3pi"] = rD_k3pi;
    obs["dD_k3pi"] = dD_k3pi * r2d;
    obs["kD_k3pi"] = kD_k3pi;
    obs["F_pipipipi"] = F_pipipipi;
    obs["rD_kpipi0"] = rD_kpipi0;
    obs["dD_kpipi0"] = dD_kpipi0 * r2d;
    obs["kD_kpipi0"] = kD_kpipi0;
    obs["F_pipipi0"] = F_pipipi0;
    obs["F_kkpi0"] = F_kkpi0;
    obs["rD_kskpi"] = rD_kskpi;
    obs["dD_kskpi"] = dD_kskpi * r2d;
    obs["kD_kskpi"] = kD_kskpi;
    obs["F_kkpipi"] = F_kkpipi;

    // Time dependent D decay and mixing
    obs["PhiM12"] = PhiM12 * r2d;
    obs["PhiG12"] = PhiG12 * r2d;
    obs["phipphig12"] = remainder(PhiG12 + phi, 2. * M_PI) * r2d;
    obs["phimphig12"] = remainder(-PhiG12 + phi, 2. * M_PI) * r2d;
    obs["AD"] = AD;
    obs["adKK"] = adKK * 1000;
    obs["adpipi"] = adpipi * 1000;
    obs["qopm1"] = (qop - 1) * 100;
    obs["qop"] = qop;
    obs["phi"] = phi * r2d;
    obs["phi12"] = phi12 * r2d;
    obs["delta"] = d;

    obs["x"] = x * 1000;
    obs["y"] = y * 1000;
    obs["M12"] = 0.5 * x12 / tau ;
    obs["G12"] = y12 / tau ;
    obs["ImM12"] = 0.5 * x12 / tau * sin(PhiM12);
  }

  return ll;
}
// ---------------------------------------------------------

void MixingModel::MCMCUserIterationInterface()
{
  std::vector<double> pars;
  for (unsigned int i = 0; i < fMCMCNChains; ++i)
  {
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
