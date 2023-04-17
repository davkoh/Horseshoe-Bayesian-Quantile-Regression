# Horseshoe-Bayesian-Quantile-Regression

Matlab and R implimentations of the horseshoe prior Bayesian quantile regression model in Kohns & Szendrei (2020). If you use the code in your work please cite as:

Kohns, D., & Szendrei, T. (2020). Horseshoe prior Bayesian quantile regression. arXiv preprint arXiv:2006.07655.

-----

The important codes are in the root directory, and the "functions" directory
*  HSBQR.m (or HSBQR.R) is the main function that you will need if you'd like to use our method.
*  MC_example.m shows how to setup the Monte-Carlo experiments done in the paper, and run the HSBQR on one of experiments (This code uses the files in the MC directory too)
*  GaR_example.m shows how to setup the empirical application as done in the paper, and run the HSBQR with iterative expanding sample sizes (This code uses the files in the GaR Application directory). Updates to the data can be downloaded from https://research.stlouisfed.org/econ/mccracken/fred-databases/. Please note that macroeconomic data can undergo heavy revision over time, so your empirical results might vary. 

-----

This code is free to use for academic purposes only, provided that the paper is cited appropriately. This code comes without technical support of any kind. It is expected to reproduce the results reported in the paper. Under no circumstances will the authors be held responsible for any use (or misuse) of this code in any way.
