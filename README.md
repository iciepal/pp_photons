This is an  example of neutral pion analysis. To run a macro, login to virgo-debian10.hpc.gsi.de and download the project from git:


git clone https://github.com/iciepal/pp_photons.git

Then, in terminal:

 . profile.sh  --> to load the latest HYDRA version

To run the code:

root -l -b loopDST_photons.C++

Event mixing:
Macro invMgg_evMix.C usues rootfile with large statistics: photons_day060_01.root

root -l -b invMgg_evMix.C
output: pp_gg_evMix.root with 2 canvas