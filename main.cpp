#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TFile.h>
//#include <mpi.h>
#include "MixingModel.h"

int main(int argc, char ** argv)
{
//    MPI_Init(&argc, &argv);
//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 6) {
        std::cout << argv[0] << " nchains nevents_pre nevents filename future [--phig12 --dream --phigamma12 --noagamma --nokspipi --nok3pi --nokp --rad --epsK --ckmcorr --combgamma --combgamma_delta]" << std::endl;
        exit(0);
    }

    bool phig12 = false;
    bool dream = false, phigamma12 = false, noexpand = true, noagamma = false, nokspipi = false, nok3pi = false, rad = false, epsK = false, ckmcorr = false, nokp = false, combgamma = false, combgamma_delta = false;

    for (int i=6;i<argc;i++){
    if (strncmp(argv[i],"--phig12",8)==0) phig12 = true;
    if (strncmp(argv[i],"--dream",7)==0) dream = true;
    if (strncmp(argv[i],"--phigamma12",12)==0) phigamma12 = true;
    if (strncmp(argv[i],"--noagamma",10)==0) noagamma = true;
    if (strncmp(argv[i],"--nokspipi",10)==0) nokspipi = true;
    if (strncmp(argv[i],"--nok3pi",8)==0) nok3pi = true;
    if (strncmp(argv[i],"--nokp",6)==0) nokp = true;
    if (strncmp(argv[i],"--rad",5)==0) rad = true;
    if (strncmp(argv[i],"--epsK",6)==0) epsK = true;
    if (strncmp(argv[i],"--ckmcorr",9)==0) ckmcorr = true;
    if (strncmp(argv[i],"--combgamma",11)==0) combgamma = true;
    if (strncmp(argv[i],"--combgamma_delta",17)==0) combgamma_delta = true;
 }

    std::cout << "phig12 = " << phig12 << ", future = " << atoi(argv[5]) << ", dream = " << dream << ", phigamma12 = " <<phigamma12 << ", expand = " << !noexpand << ", noagamma = " << noagamma << ", nokspipi = " << nokspipi << ", nok3pi = " << nok3pi << ", nokp = " << nokp << ", rad = " << rad << ", epsK = " << epsK << ", ckmcorr = " << ckmcorr << ", combgamma = " << combgamma <<  ", combgamma_delta = " << combgamma_delta << std::endl;

    // set nicer style for drawing than the ROOT default
  BCAux::SetStyle();

  // open log file
  std::string filename(argv[4]);
  mkdir(filename.c_str(), 0777);
  filename+="/";
  BCLog::OpenLog((filename+"log.txt").c_str(), BCLog::debug, BCLog::debug);
//  BCLog::SetLogLevel(BCLog::debug);

  // create new GaussModel object
  MixingModel m(phig12, noagamma, nokspipi, nok3pi, nokp, rad, epsK, ckmcorr, combgamma, combgamma_delta);

  // set MCMC precision
  m.SetNChains(atoi(argv[1]));
  m.SetNIterationsPreRunMax(atoi(argv[2]));
  m.SetNIterationsRun(atoi(argv[3]));
  //m.SetNIterationsPreRunCheck(5000);
  m.SetProposeMultivariate(true);

  // run MCMC and marginalize posterior wrt. all parameters
  // and all combinations of two parameters
    
  m.MarginalizeAll();

  // if MCMC was run before (MarginalizeAll()) it is
  // possible to use the mode found by MCMC as
  // starting point of Minuit minimization
  //m.FindMode( m.GetBestFitParameters() );

  // draw all marginalized distributions into a PostScript file
  m.PrintAllMarginalized((filename+"parameters.ps").c_str());

  // print results of the analysis into a text file
//  m.PrintResults((filename+"results.txt").c_str());

  TFile out((filename+"results.root").c_str(),"RECREATE");
  //out.WriteMarginalizedDistributions();
  //out.FillAnalysisTree();
  //out.WriteMarkovChain(true);

  m.PrintHistogram();
   
  out.Close();
  
  // close log file
  BCLog::CloseLog();
//  MPI_Finalize();
 
  return 0;

}

