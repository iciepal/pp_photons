This is an  example of neutral pion analysis. To run a macro, login to virgo-debian10.hpc.gsi.de and download the project from git:


git clone https://github.com/iciepal/pp_photons.git

Then, in the terminal:

To load the latest HYDRA version:
. /cvmfs/hadessoft.gsi.de/install/debian10/6.24.02/hydra2-6.3/defall.sh

To run the code:
root -l -b loopDST_photons.C+

Event mixing:

Macro invMgg_evMix.C usues rootfile with large statistics: photons_day060_01.root

root -l -b invMgg_evMix.C

output: pp_gg_evMix.root with 2 canvas