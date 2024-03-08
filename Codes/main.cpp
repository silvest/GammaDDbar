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
    std::cout << argv[0] << " N_chains N_events_pre N_events output_filename combination variables_filename" << std::endl << "combination = 0: Charged Beauty only" << std::endl << "combination = 1: B0d only" << std::endl << "combination = 2: B0s only" << std::endl << "combination = 3: All_modes" << std::endl;
    exit(0);
  }
  // combination = 0 Charged beauty
  // combination = 1 Neutral B0d
  // combination = 2 Neutral B0s
  // combination = 3 All_modes
  // variables_filename = name of the variables for which we want the histograms

  //----------------------------------------------------------------------------------------------------------------------------------


  //----------------------------------------- Setting the variables for the histograms ----------------------------------------
  int Nchains, Npre, Nevents, combination;
  std::string filename, Nvarfile;

  Nchains=atoi(argv[1]); // Number of Markov chains
  Npre=atoi(argv[2]); // Number of iterations to thermalize the MCMC algorithm
  Nevents = atoi(argv[3]); // Number of events
  filename=string(argv[4]); // Name of the output folder
  combination = atoi(argv[5]); // Type of combination
  Nvarfile = string(argv[6]); // Name of the file containing the variable of interest

  std::cout << "Your choice for the parameters: " << std::endl;
  cout << "Number of chains: " << Nchains << endl;
  cout << "Number of events to thermalize: " << Npre << endl;
  cout << "Number of events: " << Nevents << endl;
  cout << "Combination type: " << combination<< endl;
  cout << "Name of the output folder: "<< filename << endl;
  cout << "Path of the file containing the parameters of interest: " << Nvarfile << endl;

  std::vector<string> nParameters; // Names of parameters of interest
  std::ifstream variables_file; // file containing the names of the parameters
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
  mkdir(filename.c_str(), 0777);
  filename+="/";
  BCLog::OpenLog((filename+"log.txt").c_str(), BCLog::detail, BCLog::detail);

  // create new MixingModel object
  MixingModel m(nParameters, combination);

  // set MCMC precision
  m.SetNChains(Nchains);
  m.SetNIterationsPreRunMax(Npre);
  m.SetNIterationsRun(Nevents);
  m.SetProposeMultivariate(true);

  BCLog::OutSummary("Test model created");
  // run MCMC and marginalize posterior wrt. all parameters
  // and all combinations of two parameters

  m.MarginalizeAll();

  // draw all marginalized distributions into a PostScript file
  m.PrintAllMarginalized((filename+"parameters.ps").c_str());


  m.PrintCorrelationMatrix((filename+"coor_matrix.pdf").c_str());
  m.PrintSummary(); //Print all the relevant results in the log


  TFile out((filename+"results.root").c_str(),"RECREATE"); // root file to put our results

  m.PrintHistogram();

  out.Close();

  // close log file
  BCLog::CloseLog();

  return 0;

}
