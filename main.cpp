#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <TFile.h>
#include "MixingModel.h"

int main(int argc, char ** argv)
{

  //------------------------------------------------------  How to compile the code ---------------------------------------------------

  if (argc < 7) {
    std::cout << "To compile the code insert the following arguments: " << std::endl ;
    std::cout << argv[0] << " Nchains Nevents_pre Nevents output_filename combination variables_filename" << std::endl << "combination = 0: Charged Beauty" << std::endl << "combination = 1: B0d" << std::endl << "combination = 2: B0s" << std::endl << "combination = 3: All_modes" << std::endl;
    std::cout  << "Add other arguments as follows: " << std::endl << "future [--phig12 --dream --phigamma12 --noagamma --nokspipi --nok3pi --nokp --rad --epsK --ckmcorr --combgamma --combgamma_delta] " << std::endl;
    exit(0);
  }
  // combination = 0 Charged beauty
  // combination = 1 Neutral B0d
  // combination = 2 Neutral B0s
  // combination = 3 All_modes
  // variables_filename = name of the variables for which we want the histograms

  //----------------------------------------------------------------------------------------------------------------------------------

  //-------------------------------  Setting  all the other arguments to zero ---------------------------------------------
  bool phig12 = false;
  bool dream = false, phigamma12 = false, noexpand = true , noagamma = false, nokspipi = false, nok3pi = false, rad = false, epsK = false, ckmcorr = false, nokp = false, combgamma = false, combgamma_delta = false;
  int combination = 0;
  //--------------------------------------------------------------------------------------------------------------------------

  //-------------------------------  Setting all the variables to the values in input ---------------------------------------------

  combination = atoi(argv[5]);

  for (int i=7;i<argc;i++){
    if (strncmp(argv[i],"--phig12",8)==0) phig12 = true; // turn on the phi^Gamma contribution
    if (strncmp(argv[i],"--dream",7)==0) dream = true; // unused in this code
    if (strncmp(argv[i],"--phigamma12",12)==0) phigamma12 = true; // unused in this code
    if (strncmp(argv[i],"--noagamma",10)==0) noagamma = true; // turn off the Agamma contribution
    if (strncmp(argv[i],"--nokspipi",10)==0) nokspipi = true; // unused in this code
    if (strncmp(argv[i],"--nok3pi",8)==0) nok3pi = true; // unused in this code
    if (strncmp(argv[i],"--nokp",6)==0) nokp = true; // unused in this code
    if (strncmp(argv[i],"--rad",5)==0) rad = true; // If rad is true print the root histograms in rad, otherwise they are defined in degrees
    if (strncmp(argv[i],"--epsK",6)==0) epsK = true; // If true takes into account the correction due to epsK
    if (strncmp(argv[i],"--ckmcorr",9)==0) ckmcorr = true; // If true insert the Rcp parameter
    if (strncmp(argv[i],"--combgamma",11)==0) combgamma = true; // use the LHCb combination results
    if (strncmp(argv[i],"--combgamma_delta",17)==0) combgamma_delta = true; // use the LHCb combination results for D -> Kpi channel
    //noexpand unused in this code
  }

  std::cout << "Your choice for the auxiliary variables: " << std::endl;
  std::cout << "phig12 = " << phig12 << ", future = " << atoi(argv[6]) << ", dream = " << dream << ", phigamma12 = " <<phigamma12 << ", expand = " << !noexpand << ", noagamma = " << noagamma << ", nokspipi = " << nokspipi << ", nok3pi = " << nok3pi << ", nokp = " << nokp << ", rad = " << rad << ", epsK = " << epsK << ", ckmcorr = " << ckmcorr << ", combgamma = " << combgamma <<  ", combgamma_delta = " << combgamma_delta << std::endl;

  //----------------------------------------- Setting the variables for the histograms ----------------------------------------

  string Nvarfile;
  std::vector<string> nParameters;
  std::ifstream variables_file;
  Nvarfile = string(argv[6]);
  variables_file.open(Nvarfile);
  string word;
  while ( variables_file >> word) {
    nParameters.push_back(word);
  }
  std::cout << "The histograms will be filled for the following variables: " << std::endl ;
  for(int i = 0; i < nParameters.size(); i++){
    std::cout << nParameters[i] << std::endl;
  }
  variables_file.close();
  //--------------------------------------------------------------------------------------------------------------------------

  //-------------------------------  Performing the fit  ------------------------------------------------------------------

  // open log file
  std::string filename(argv[4]);
  mkdir(filename.c_str(), 0777);
  filename+="/";
  BCLog::OpenLog((filename+"log.txt").c_str(), BCLog::detail, BCLog::detail);

  // create new MixingModel object
  MixingModel m(nParameters, combination, phig12, noagamma, nokspipi, nok3pi, nokp, rad, epsK, ckmcorr, combgamma, combgamma_delta);

  // set MCMC precision
  m.SetNChains(atoi(argv[1]));
  m.SetNIterationsPreRunMax(atoi(argv[2]));
  m.SetNIterationsRun(atoi(argv[3]));
  m.SetProposeMultivariate(true);

  BCLog::OutSummary("Test model created");
  // run MCMC and marginalize posterior wrt. all parameters
  // and all combinations of two parameters

  m.MarginalizeAll();

  // draw all marginalized distributions into a PostScript file
  m.PrintAllMarginalized((filename+"parameters.ps").c_str());


  m.PrintCorrelationMatrix((filename+"coor_matrix.pdf").c_str());
  m.PrintSummary(); //Print all the relevant results in the log


  TFile out((filename+"results.root").c_str(),"RECREATE");

  m.PrintHistogram();

  out.Close();

  // close log file
  BCLog::CloseLog();
  //  MPI_Finalize();

  return 0;

}
