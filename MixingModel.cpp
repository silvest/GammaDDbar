#include "MixingModel.h"

#include <BAT/BCMath.h>
#include <BAT/BCGaussianPrior.h>

#include <TRandom3.h>

// ---------------------------------------------------------

MixingModel::MixingModel(bool phig12_i, bool noagamma_i, bool nokspipi_i, bool nok3pi_i, bool nokp_i, bool rad_i, bool epsK_i, bool ckmcorr_i, bool combgamma_all_i, bool combgamma_delta_i) : BCModel(), histos(obs)
{
    TH1::SetDefaultBufferSize(1000000); 	
    phig12 = phig12_i;
    noagamma = noagamma_i;
    nokspipi = nokspipi_i;
    nok3pi = nok3pi_i;
    nokp = nokp_i;
    rad = rad_i;
    epsK = epsK_i;
    ckmcorr = ckmcorr_i;
    combgamma_all = combgamma_all_i;
    combgamma_delta = combgamma_delta_i;

    tauD = .410; //ps

    epsI = 2.228 * sin(43.52 / 180. * M_PI);

    vector<dato> CorrData;
    TMatrixDSym Corr, Corr2;

    d2r = M_PI / 180.;
    r2d = (rad ? 1. : 180. / M_PI);

    // UTfit summer 2018 web page A^2 lambda^4 eta
    RCP_in = 0.826 * 0.826 * 0.225 * 0.225 * 0.225 * 0.225 * 0.357;


    // UTfit summer 2018 web page A^2 lambda^4 eta
    RCP_in = 0.826 * 0.826 * 0.225 * 0.225 * 0.225 * 0.225 * 0.357;
    // Yellow report current
    RCP_err = (2. * 1.5e-2 + 4. * 0.12e-2 + 3.e-2) * RCP_in;

    // Experimental data from HFLAV Moriond 2019 + LHCb AGamma 1911.01114 + Belle yCP 1912.10912 (irrilevante)

    if (combgamma_all)
        meas.insert(pair<string, dato>("ycp", dato(0.835e-2, 0.155e-2))); //without LHCb 2019
    else
        meas.insert(pair<string, dato>("ycp", dato(0.719e-2, 0.113e-2))); //HFLAV 2020
    if (!noagamma && !combgamma_all) {
        if (phig12)
            meas.insert(pair<string, dato>("AGamma", dato(1.9e-4, 1.3e-4, 0.4e-4))); //LHCb 2021
        else
            meas.insert(pair<string, dato>("AGamma", dato(1.0e-4, 1.1e-4, 0.3e-4))); //LHCb 2021
    }
    //      meas.insert(pair<string, dato>("AGamma", dato(-0.0309e-2, 0.0204e-2))); //HFLAV 2020

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

    meas.insert(pair<string, dato>("RM", dato(0.013e-2, 0.0269e-2))); //HFLAV

    CorrData.clear();

    CorrData.push_back(dato(0.16e-2, 0.23e-2, 0.12e-2, 0.08e-2)); // x_exp_babar
    CorrData.push_back(dato(0.57e-2, 0.20e-2, 0.13e-2, 0.07e-2)); // y_exp_babar

    Corr.ResizeTo(2, 2);

    Corr(0, 0) = 1.;
    Corr(0, 1) = 0.0615;
    Corr(1, 1) = 1.;

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kppkk",
            CorrelatedGaussianObservables(CorrData, Corr)));

    CorrData.clear();

    CorrData.push_back(dato(2.48e-2, 0.59e-2, 0.39e-2)); // xpp_kpp_babar_exp
    CorrData.push_back(dato(-.07e-2, .65e-2, .5e-2)); // ypp_kpp_babar_exp

    Corr.ResizeTo(2, 2);

    Corr(0, 0) = 1.;
    Corr(0, 1) = -0.69;
    Corr(1, 1) = 1.;

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpp_babar_plus",
            CorrelatedGaussianObservables(CorrData, Corr)));

    CorrData.clear();

    CorrData.push_back(dato(3.5e-2, 0.78e-2, 0.65e-2)); // xpm_kpp_babar_exp
    CorrData.push_back(dato(-0.82e-2, .68e-2, .41e-2)); // ypm_kpp_babar_exp

    Corr.ResizeTo(2, 2);

    Corr(0, 0) = 1.;
    Corr(0, 1) = -0.66;
    Corr(1, 1) = 1.;

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpp_babar_minus",
            CorrelatedGaussianObservables(CorrData, Corr)));

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

    CorrData.clear();

    CorrData.push_back(dato(-2.1e-2, 5.4e-2)); // AD_exp_babar_kp
    CorrData.push_back(dato(-0.020e-2, .050e-2)); // xp2minus_exp_babar_kp
    CorrData.push_back(dato(0.96e-2, .75e-2)); // ypminus_exp_babar_kp

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpi_babar_minus",
            CorrelatedGaussianObservables(CorrData, Corr)));

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

    CorrData.clear();

    CorrData.push_back(dato(2.3e-2, 4.7e-2)); // AD_exp_belle_kp
    CorrData.push_back(dato(0.006e-2, .034e-2)); // xp2minus_exp_belle_kp
    CorrData.push_back(dato(0.20e-2, .54e-2)); // ypminus_exp_belle_kp

    corrmeas.insert(pair<string, CorrelatedGaussianObservables>("kpi_belle_minus",
            CorrelatedGaussianObservables(CorrData, Corr)));

    if (!combgamma_all) {
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


        meas.insert(pair<string, dato>("RMKppp", dato(0.0096e-2, 0.0036e-2))); //LHCb Kppp multiplied by 2

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



    // Observables

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
    if (ckmcorr) {
        histos.createH1D("phiKsLHCb", 200, 1., -1.);
        histos.createH1D("RCP", 200, 1., -1.);
    }
    DefineParameters();

};

// ---------------------------------------------------------

MixingModel::~MixingModel()
{
    // default destructor
};

// ---------------------------------------------------------

void MixingModel::DefineParameters()
{
    // add parameters

    AddParameter("x12", 0., 2.e-2);
    AddParameter("y12", 0., 5.e-2);
    AddParameter("PhiM12", -.2, .2);
    AddParameter("delta_kp", -M_PI, M_PI);
    AddParameter("delta_kpp", -M_PI, M_PI);
    AddParameter("Rd", 0., 1.e-2);
    if (phig12)
        AddParameter("phiG12", -0.2, 0.2);
    else
        AddParameter("phiG12", 0., 0.);

    if (ckmcorr)
        AddParameter("RCP", RCP_in - 5. * RCP_err, RCP_in + 5. * RCP_err);

    SetPriorConstantAll();
    if (ckmcorr)
        GetParameter("RCP").SetPrior(new BCGaussianPrior(RCP_in, RCP_err));


}

// ---------------------------------------------------------

double MixingModel::LogLikelihood(const std::vector<double> &parameters)
{
    double ll = 0.;
    x12 = parameters[0];
    y12 = parameters[1];
    phiM12 = parameters[2];

    delta_kpi = parameters[3];
    delta_kpipi = parameters[4];
    Rd = parameters[5];
    phiG12 = parameters[6];
    if (ckmcorr)
        RCP = parameters[7];

    //compute all observables

    phi12 = remainder(-phiG12 + phiM12, 2. * M_PI);
    x = x12;
    y = y12;
    rm = (x * x + y * y) / 2.;
    qop = 1. + x12 * y12 * sin(phi12) / (x12 * x12 + y12 * y12);
    phi = atan(-(x12 * x12 * sin(2. * phiM12) + y12 * y12 * sin(2. * phiG12)) /
            (x12 * x12 * cos(2. * phiM12) + y12 * y12 * cos(2. * phiG12))) / 2.;
    d = (1. - qop * qop) / (1. + qop * qop);

    phiKsBfact = phi;
    phiKsLHCb = phi;
    if (epsK)
        phiKsBfact += -2 * epsI;
    if (ckmcorr) {
        //    std::cout << RCP << std::endl;
        phiKsBfact -= RCP;
        phiKsLHCb -= RCP;
    }


    xp = x * cos(delta_kpi) + y * sin(delta_kpi);
    yp = -x * sin(delta_kpi) + y * cos(delta_kpi);


    yp_plus = qop * (yp * cos(phi) - xp * sin(phi));
    xp_plus = qop * (xp * cos(phi) + yp * sin(phi));
    xp_plus_sq = xp_plus * xp_plus;

    yp_minus = 1. / qop * (yp * cos(phi) + xp * sin(phi));
    xp_minus = (xp * cos(phi) - yp * sin(phi)) / qop;
    xp_minus_sq = xp_minus * xp_minus;

    double xpkpp = x * cos(delta_kpipi) + y * sin(delta_kpipi);
    double ypkpp = -x * sin(delta_kpipi) + y * cos(delta_kpipi);

    yp_kpp_plus = qop * (ypkpp * cos(phi) - xpkpp * sin(phi));
    xp_kpp_plus = qop * (xpkpp * cos(phi) + ypkpp * sin(phi));

    //no direct CPV in tree decays
    Ad = 0.;

    Rdkp = Rd * (1. + Ad);

    yp_kpp_minus = 1. / qop * (ypkpp * cos(phi) + xpkpp * sin(phi));
    xp_kpp_minus = (xpkpp * cos(phi) - ypkpp * sin(phi)) / qop;

    ycp = (qop + 1. / qop) / 2. * y * cos(phi) - (qop - 1. / qop) / 2. * x * sin(phi);
    Ag = (qop - 1. / qop) / 2. * y * cos(phi)- (qop + 1. / qop) / 2. * x * sin(phi);
    ycpksppLHCb = (qop + 1. / qop) / 2. * y * cos(phiKsLHCb) - (qop - 1. / qop) / 2. * x * sin(phiKsLHCb);
    AgksppLHCb = (qop - 1. / qop) / 2. * y * cos(phiKsLHCb)- (qop + 1. / qop) / 2. * x * sin(phiKsLHCb);
    xcp = (qop + 1. / qop) / 2. * x * cos(phiKsLHCb) + (qop - 1. / qop) / 2. * y * sin(phiKsLHCb);
    dx = (qop - 1. / qop) / 2. * x * cos(phiKsLHCb) + (qop + 1. / qop) / 2. * y * sin(phiKsLHCb);
    am = (qop * qop - 1. / (qop * qop)) / (qop * qop + 1. / (qop * qop));

    obs["y"] = y;
    obs["x"] = x;
    obs["delta_kp"] = delta_kpi * r2d;
    obs["delta_kpp"] = delta_kpipi * r2d;
    obs["phi"] = phi * r2d;
    obs["phiG12"] = phiG12 * r2d;
    obs["phipphig12"] = remainder(phiG12 + phi, 2. * M_PI) * r2d;
    obs["phimphig12"] = remainder(-phiG12 + phi, 2. * M_PI) * r2d;
    obs["phiM12"] = phiM12 * r2d;
    obs["phi12"] = phi12 * r2d;
    obs["x12"] = x12;
    obs["y12"] = y12;
    obs["qopm1"] = qop - 1.;
    obs["delta"] = d;
    obs["RD"] = Rd;
    obs["AD"] = Ad;
    obs["RDkp"] = Rdkp;
    obs["rDkp"] = sqrt(Rdkp);
    obs["Rm"] = rm;
    obs["Am"] = am;
    obs["AGamma"] = Ag;
    obs["ycp"] = ycp;
    obs["dx"] = dx;
    obs["xcp"] = xcp;
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

    TVectorD corr(4);


    corr(0) = x;
    corr(1) = y;
    corr(2) = qop;
    corr(3) = phiKsBfact;

    ll += corrmeas.at("kpp_belle").logweight(corr);

    ll += meas.at("RM").logweight(rm);
    if (!combgamma_all)
        ll += meas.at("RMKppp").logweight(rm);

    corr.ResizeTo(2);

    corr(0) = yp_kpp_plus;
    corr(1) = xp_kpp_plus;

    ll += corrmeas.at("kpp_babar_plus").logweight(corr);

    corr(0) = yp_kpp_minus;
    corr(1) = xp_kpp_minus;

    ll += corrmeas.at("kpp_babar_minus").logweight(corr);

    ll += meas.at("ycp").logweight(ycp);

    //    cout << "ycp " << ycp_ex.weight(ycp) << endl;

    if (!noagamma && !combgamma_all)
        ll += meas.at("AGamma").logweight(Ag);


    //Babar D0 -> Ks pi pi and D0 -> Ks K K 

    corr.ResizeTo(2);
    corr(0) = x;
    corr(1) = y;

    ll += corrmeas.at("kppkk").logweight(corr);


    //CLEOc Psi results

    corr.ResizeTo(5);

    corr(0) = Rd;
    corr(1) = x*x;
    corr(2) = y;
    corr(3) = cos(delta_kpi);
    corr(4) = sin(delta_kpi);

    ll += corrmeas.at("cleoc").logweight(corr);

    // Babar Kpi & Belle Kpi

    corr.ResizeTo(3);

    corr(0) = Rd;
    corr(1) = xp_plus_sq;
    corr(2) = yp_plus;

    ll += corrmeas.at("kpi_babar_plus").logweight(corr);
    ll += corrmeas.at("kpi_belle_plus").logweight(corr);

    corr(0) = Ad;
    corr(1) = xp_minus_sq;
    corr(2) = yp_minus;

    ll += corrmeas.at("kpi_babar_minus").logweight(corr);
    ll += corrmeas.at("kpi_belle_minus").logweight(corr);

    corr(0) = Rd;
    corr(1) = (yp_plus + yp_minus) / 2.;
    corr(2) = (xp_plus_sq + xp_minus_sq + yp_plus * yp_plus + yp_minus * yp_minus) / 2. - corr(1) * corr(1);

    ll += corrmeas.at("kpi_cdf").logweight(corr);

    if (!combgamma_all) {
        corr.ResizeTo(5);

        corr(0) = Rd;
        corr(1) = yp_plus;
        corr(2) = xp_plus_sq;
        corr(3) = yp_minus;
        corr(4) = xp_minus_sq;

        ll += corrmeas.at("kpi_lhcb").logweight(corr);

        corr.ResizeTo(4);

        corr(0) = xcp;
        corr(1) = ycpksppLHCb;
        corr(2) = dx;
        corr(3) = AgksppLHCb;

        ll += corrmeas.at("LHCb_kspp").logweight(corr);
    } else {
        corr.ResizeTo(6);

        corr(0) = x;
        corr(1) = y;
        corr(2) = qop;
        corr(3) = phi / d2r;
        corr(4) = sqrt(Rd);
        corr(5) = delta_kpi / d2r;

        ll += corrmeas.at("LHCb_combo").logweight(corr);
    }
    if (combgamma_delta) {
        corr.ResizeTo(2);

        corr(0) = sqrt(Rd);
        corr(1) = delta_kpi / d2r;

        ll += corrmeas.at("LHCb_combo_delta").logweight(corr);

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
