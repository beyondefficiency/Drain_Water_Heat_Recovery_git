This repository contains the tools needed to predict the savings from a DWHR device. The main documentation of the methods can be found
in DWHR_SavingsEstimates.py. Important files are:

\Coefficients - This folder contains coefficients for the DWHR performance maps. The values in these performance maps came from previous
DWHR research.

\Profiles - This folder contains the draw profiles that you wish to analyze. When the script finishes analyzing the draw profiles it will
save the results in this folder with '_Analyzed' at the end of the file name.

DWHR_Savings_Estimates.py - This is the analysis script. It references performance map coefficients in \Coefficients and draw profiles
stored in \Profiles to predict the effectiveness and energy savings of a rated DWHR device. See the detailed documentation in the script.

