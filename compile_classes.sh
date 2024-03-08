
codes_folder="$PWD/Codes" ## path to the folder containing the classes
Path_to_ROOTSYS="$ROOTSYS/bin/root-config --cflags --libs"
Path_to_BAT_config="bat-config --cflags"
Path_to_BAT_libs="bat-config --libs"

g++ -c "$codes_folder/histo.cpp" `$Path_to_ROOTSYS` `$Path_to_BAT_config` `$Path_to_BAT_libs`
g++ -c "$codes_folder/CorrelatedGaussianObservables.cpp" `$Path_to_ROOTSYS` `$Path_to_BAT_config` `$Path_to_BAT_libs`
g++ -c "$codes_folder/MixingModel.cpp" `$Path_to_ROOTSYS` `$Path_to_BAT_config` `$Path_to_BAT_libs`
