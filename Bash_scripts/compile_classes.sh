
wd=$PWD ## Path of the Fit_code folder
codes_folder="$wd/Codes" ## path to the folder containing the classes

g++ -c "$codes_folder/histo.cpp" `$ROOTSYS/bin/root-config --cflags --libs` `bat-config --cflags` `bat-config --libs`
g++ -c "$codes_folder/CorrelatedGaussianObservables.cpp" `$ROOTSYS/bin/root-config --cflags --libs` `bat-config --cflags` `bat-config --libs`
g++ -c "$codes_folder/MixingModel.cpp" `$ROOTSYS/bin/root-config --cflags --libs` `bat-config --cflags` `bat-config --libs`
