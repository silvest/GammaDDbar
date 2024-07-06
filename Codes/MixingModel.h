#ifndef __MIXINGMODEL__H
#define __MIXINGMODEL__H

#include <BAT/BCModel.h>
#include <BAT/BCH2D.h>
#include <TH2D.h>
#include <map>
#include "histo.h"
#include "dato.h"
#include "CorrelatedGaussianObservables.h"
#include <iostream>
#include <cmath>

// --------------------------------------------------------- This is the class to make the model ----------------------------------------------------

using namespace std;

class MixingModel : public BCModel {
public:

  // Constructors and destructor
  MixingModel(vector<string> nParam, int combination);
  ~MixingModel();

  // Methods to overload, see file MixingModel.cpp
  void DefineParameters(); // Define the parameters
  double LogLikelihood(const std::vector<double> &parameters); // Compute the log likelihood
  void MCMCUserIterationInterface();
  void PrintHistogram(); // This is to print the histograms

  //Methods to insert the measurements
  void Add_ChargedB_meas(); //Charged B decay chain modes
  void Add_NeutralBd_meas(); //Neutral Bd mixing and decay modes
  void Add_NeutralBs_meas(); //B0s mixing and decay modes
  void Add_time_dependent_Dmeas(); // D mixing and decay modes
  void Add_other_meas(); // other inputs
  void Add_old_meas(); // old D analysis measurements

  //Histograms
  void DefineHistograms(); // Function to define the histograms to fill

  //Boolean variables to set the combination
  int comb; // combination variable
  // vector to copy the names of the variables to fill the relative histograms
  std::vector<string> nVarab; // name of the parameters to fill the histograms


  map<string,dato> meas;
  map<string,CorrelatedGaussianObservables> corrmeas;
  map<string,double> obs;
  histo histos;

  //PARAMETERS

  //Combination parameters
  double g, x12, y12, r_dk, r_dpi, rD_kpi, d_dk, d_dpi, dD_kpi, //9
  rD_k3pi, dD_k3pi, kD_k3pi, F_pipipipi, //4
  rD_kpipi0, dD_kpipi0, kD_kpipi0, F_pipipi0, F_kkpi0, //5
  rD_kskpi, dD_kskpi, kD_kskpi, RBRdkdpi, //4
  r_dstk, d_dstk, r_dstpi, d_dstpi, //4
  r_dkst, d_dkst, k_dkst, //3
  r_dkstz, d_dkstz, k_dkstz, //3
  r_dkpipi, d_dkpipi, k_dkpipi, r_dpipipi, d_dpipipi, k_dpipipi, //6
  l_dsk, d_dsk, phis, //3
  l_dskpipi, d_dskpipi, k_dskpipi, //3
  l_dmpi, d_dmpi, beta, // 3
  PhiM12, PhiG12, AD, Delta_totau, DeltaAcp, //5
  //52 total parameters
  F_kkpipi, //airXiv 20301.10328
  // 2401.17934 Bs
  r_dkstzs, d_dkstzs, k_dkstzs; //3

  //old parameter
  double Rd, Rdkp, d; // rdkpi^2

  //auxiliary parameters
  double x,y, phi12, qop, phi; // other parametrization of the mixing parameters

  //OBSERVABLES
  //Combination Observables
  double  acp_dk_uid0, acp_dpi_uid0, afav_dk_uid0, rcp_uid0, rm_dk_uid0, rm_dpi_uid0, rp_dk_uid0, rp_dpi_uid0, //UID0
  aads_dk_k3pi_uid1, aads_dpi_k3pi_uid1, acp_dk_4pi_uid1, acp_dpi_4pi_uid1, afav_dk_k3pi_uid1, rads_dk_k3pi_uid1, rads_dpi_k3pi_uid1, rcp_4pi_uid1, //UID1
  aads_dk_kpipi0_uid2, aads_dpi_kpipi0_uid2, acp_dk_kkpi0_uid2, acp_dk_pipipi0_uid2, acp_dpi_kkpi0_uid2, acp_dpi_pipipi0_uid2, afav_dk_kpipi0_uid2, rads_dk_kpipi0_uid2, rads_dpi_kpipi0_uid2, rcp_kkpi0_uid2, rcp_pipipi0_uid2, //UID2
  xm_dk_uid3, ym_dk_uid3, xp_dk_uid3, yp_dk_uid3, xi_x_dpi_uid3, xi_y_dpi_uid3, //UID3
  afav_dpi_kskpi_uid4, asup_dpi_kskpi_uid4, afav_dk_kskpi_uid4, asup_dk_kskpi_uid4, rfavsup_dpi_kskpi_uid4, rfav_dkdpi_kskpi_uid4, rsup_dkdpi_kskpi_uid4, //UID4
  acp_dstk_dg_uid5, acp_dstk_dp_uid5, afav_dstk_dg_uid5, afav_dstk_dp_uid5, rcp_dg_uid5, rcp_dp_uid5, rm_dstk_dg_uid5, rm_dstk_dp_uid5,  rp_dstk_dg_uid5, rp_dstk_dp_uid5, acp_dstpi_dg_uid5, acp_dstpi_dp_uid5,
  rm_dstpi_dg_uid5, rm_dstpi_dp_uid5, rp_dstpi_dg_uid5, rp_dstpi_dp_uid5, //UID5
  afav_dkst_kpi_uid6, acp_dkst_kk_uid6, acp_dkst_pipi_uid6, rcp_dkst_kk_uid6, rcp_dkst_pipi_uid6, rp_dkst_kpi_uid6, rm_dkst_kpi_uid6, afav_dkst_k3pi_uid6,
  acp_dkst_pipipipi_uid6, rcp_dkst_pipipipi_uid6, rp_dkst_k3pi_uid6, rm_dkst_k3pi_uid6, //UID6
  acp_dkstz_kk_uid7, acp_dkstz_pipi_uid7, rcp_dkstz_kk_uid7, rcp_dkstz_pipi_uid7, acp_dkstz_4pi_uid7, rcp_dkstz_4pi_uid7, rp_dkstz_kpi_uid7,
  rm_dkstz_kpi_uid7, rp_dkstz_k3pi_uid7, rm_dkstz_k3pi_uid7, afav_dkstz_kpi_uid7, afav_dkstz_k3pi_uid7, //UID7
  xm_dkstz_uid8, ym_dkstz_uid8, xp_dkstz_uid8, yp_dkstz_uid8, //UID8
  rcp_dkpipi_uid9, afav_dkpipi_kpi_uid9, afav_dpipipi_kpi_uid9, acp_dkpipi_kk_uid9, acp_dkpipi_pipi_uid9, acp_dpipipi_kk_uid9,
  acp_dpipipi_pipi_uid9, rp_dkpipi_uid9, rm_dkpipi_uid9, rp_dpipipi_uid9, rm_dpipipi_uid9, //UID9
  c_dsk_uid10, d_dsk_uid10, db_dsk_uid10, s_dsk_uid10, sb_dsk_uid10, //UID10
  c_dskpipi_uid11, d_dskpipi_uid11, db_dskpipi_uid11, s_dskpipi_uid11, sb_dskpipi_uid11, //UID11
  s_dmpi_uid12, sb_dmpi_uid12, //UID12
  x_uid13, y_uid13, //UID13
  xcp_uid14, ycp_uid14, dx_uid14, dy_uid14, //UID14
  xcp_uid15, ycp_uid15, dx_uid15, dy_uid15, //UID15
  Rdp_uid16, yp_uid16, xpsq_uid16, Rdm_uid16, ym_uid16, xmsq_uid16, //UID16
  dtotau_uid17, dacp_uid17, //UID17
  k3pi_uid18, //UID18
  kD_k3pi_uid19, dD_k3pi_uid19, kD_kpipi0_uid19, dD_kpipi0_uid19, rD_k3pi_uid19, rD_kpipi0_uid19, //UID19
  F_pipipipi_uid20, //UID20
  F_pipipi0_uid21, F_kkpi0_uid21, //UID21
  RD_kskpi_uid22, //UID22
  RD_kskpi_uid23, dD_kskpi_uid23, kD_kskpi_uid23, //UID23
  k_dkst_uid24,  //UID24
  k_dkstz_uid25, //UID25
  phis_uid26, //UID26
  beta_uid27, //UID27
  ycp_uid28, //UID28
  DY_uid29, //UID29
  Rdp_uid30, yp_uid30, xpsq_uid30, Rdm_uid30, ym_uid30, xmsq_uid30; //UID30

  double xcp, ycp, dx, dy; // D Observables CP final states

  //old_observables
  double rm,
  yp_kpp_plus, xp_kpp_plus,
  yp_kpp_minus, xp_kpp_minus,
  xp_plus, xp_plus_sq, yp_plus,
  xp_minus, xp_minus_sq, yp_minus;

  //https://arxiv.org/abs/2208.09402 observables
  double Akpi_BESIII, Akpi_kpipi0_BESIII,
         xi_x_BESIII, xi_y_BESIII;

  //https://arxiv.org/pdf/2208.10098.pdf
  double F_pipipipi_BESIII;

  //https://arxiv.org/abs/2301.10328
  double acp_dk_kkpipi_230110328, acp_dpi_kkpipi_230110328, acp_dk_pipipipi_230110328, acp_dpi_pipipipi_230110328,
         rcp_kpi_kkpipi_230110328, rcp_kpi_pipipipi_230110328;

  //https://arxiv.org/pdf/2401.17934.pdf
  //Bd observables
  double acp_dkstz_kk_240117934Bd, acp_dkstz_pipi_240117934Bd, rcp_dkstz_kk_240117934Bd, rcp_dkstz_pipi_240117934Bd, acp_dkstz_4pi_240117934Bd, rcp_dkstz_4pi_240117934Bd, rp_dkstz_kpi_240117934Bd,
           rm_dkstz_kpi_240117934Bd, rp_dkstz_k3pi_240117934Bd, rm_dkstz_k3pi_240117934Bd, afav_dkstz_kpi_240117934Bd, afav_dkstz_k3pi_240117934Bd;

  //Bs observables
  double acp_dkstz_kk_240117934Bs, acp_dkstz_pipi_240117934Bs, rcp_dkstz_kk_240117934Bs, rcp_dkstz_pipi_240117934Bs, acp_dkstz_4pi_240117934Bs, rcp_dkstz_4pi_240117934Bs, rp_dkstz_kpi_240117934Bs,
            rm_dkstz_kpi_240117934Bs, rp_dkstz_k3pi_240117934Bs, rm_dkstz_k3pi_240117934Bs, afav_dkstz_kpi_240117934Bs, afav_dkstz_k3pi_240117934Bs;

  // Slides Ribatti
  double CKpi, CpKpi, DCKpi, DCpKpi;

  //Methods to calculate the observables and the Log-Likelihood
  double Calculate_ChargedB_observables();
  double Calculate_neutralBdobservables();
  double Calculate_neutralBsobservables();
  double Calculate_time_dependent_Dobservables();
  double Calculate_other_observables();
  double Calculate_old_observables();

  //General structure of the fit equations
  double Acp(double rB, double delta_B, double kB, double F_D, double alpha );
  double Rcp_h(double rBCP, double delta_BCP, double rBCF, double rD, double delta_BCF, double delta_D, double kB, double kD, double F_D, double alpha);
  double Afav(double rB, double rD, double delta_B, double delta_D, double kB, double kD , double alpha );
  double Asup(double rB, double rD, double delta_B, double delta_D, double kB, double kD , double alpha);
  double Rm(double rB, double rD, double delta_B, double delta_D, double kB, double kD, double alpha);
  double Rp(double rB, double rD, double delta_B, double delta_D, double kB, double kD, double alpha);
  double Rads(double rB, double rD, double delta_B, double delta_D, double kB, double kD, double alpha);
  double Rfav(double rB1, double rB2, double rD, double delta_B1, double delta_B2, double delta_D, double BR, double kD, double alpha);
  double Rsup(double rB1, double rB2, double rD, double delta_B1, double delta_B2, double delta_D, double BR, double kD, double alpha);
  double y_plus(double delta_D);
  double x_plus(double delta_D);
  double y_minus(double delta_D);
  double x_minus(double delta_D);

private:
  double d2r, r2d;
  double tau; // Mean lifetime

};
// ---------------------------------------------------------

#endif
