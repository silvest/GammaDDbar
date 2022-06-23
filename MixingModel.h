#ifndef __MIXINGMODEL__H
#define __MIXINGMODEL__H

#include <BAT/BCModel.h>
#include <BAT/BCH2D.h>
#include <TH2D.h>
#include <map>
#include "histo.h"
#include "dato.h"
#include "CorrelatedGaussianObservables.h"

// ---------------------------------------------------------

using namespace std;

class MixingModel : public BCModel {
public:

    // Constructors and destructor
  MixingModel(bool phig12_i, bool noagamma_i, bool nokspipi_i, bool nok3pi_i, bool nokp_i, bool rad_i, bool epsK_i, bool ckmcorr_i, bool combgamma_all_i, bool combgamma_delta_i);
    ~MixingModel();

    // Methods to overload, see file MixingModel.cpp
    void DefineParameters();
    double LogLikelihood(const std::vector<double> &parameters);
    void MCMCUserIterationInterface(); 
    void PrintHistogram();
        
  bool phig12, noagamma, nokspipi, nok3pi, rad, epsK, ckmcorr, nokp, combgamma_all, combgamma_delta;
    double tauD;
    map<string,dato> meas;
    map<string,CorrelatedGaussianObservables> corrmeas;
    map<string,double> obs;
    histo histos;
    double Rd,Ad,Rdkp,rm,yp_plus,yp_minus,xp_plus,xp_plus_sq,xp_minus,xp_minus_sq,phiM12,
      xp_kpp_plus,xp_kpp_minus,yp_kpp_plus,yp_kpp_minus,x,y,phi,qop,x12,y12,delta_kpi,xp,yp,d,delta_kpipi,
      phiG12,phi12,phi12A,ycp,Ag,xcp,dx,am,ycpksppLHCb,AgksppLHCb;
  double phiKsBfact, phiKsLHCb, RCP, RCP_in, RCP_err;
    private:
    double avexf, errxf, aveyf, erryf, aveqopf, errqopf, avephif, errphif, aveagammaf, erragammaf, aveycpf,
            avexcpf, errxcpf, avedyf, errdyf, avedxf, errdxf, errycpf,
      aveyppf,erryppf,aveypmf,errypmf,avexpp2f,errxpp2f,avexpm2f,errxpm2f,aveRDf,errRDf,r2d,d2r, epsI, avephiKsBfactf, avephiKsLHCbf, aveycpksppLHCbf, aveAgksppLHCbf, avexpf, aveypf, delta_kpi_in, avermf;


};
// ---------------------------------------------------------

#endif

