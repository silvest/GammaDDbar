codes_folder="$PWD/Codes" ## path to the folder containing the script

Nchains=$1 ## Number of chains
Nevents_pre=$2 ## Number of events per chain to run to thermalize the MCMC
Nevents=$3 ## Number of events per chain used to do the statistics
output_filename=$4 ## Name of the output folder
Comb_type=$5 ## 0, 1, 2, 3 Charged, B0d, B0s, All
variables_folder="$PWD/Variables/$6"
Path_to_ROOTSYS="$ROOTSYS/bin/root-config --cflags --libs"
Path_to_BAT_config="bat-config --cflags"
Path_to_BAT_libs="bat-config --libs"

g++ -o main.x "$codes_folder/main.cpp" `$Path_to_ROOTSYS` `$Path_to_BAT_config` `$Path_to_BAT_libs` histo.o CorrelatedGaussianObservables.o MixingModel.o

time ./main.x $Nchains $Nevents_pre $Nevents $output_filename $Comb_type $variables_folder
