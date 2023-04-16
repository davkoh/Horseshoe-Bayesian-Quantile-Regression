# Horseshoe-Bayesian-Quantile-Regression

Matlab and R implimentations of the horseshoe prior Bayesian quantile regression model in Kohns & Szendrei (2020). If you use the code in your work please cite as:

Kohns, D., & Szendrei, T. (2020). Horseshoe prior Bayesian quantile regression. arXiv preprint arXiv:2006.07655.

-----

The important codes are in the root directory, and the "functions" directory
*  HSBQR.m is the main function that you will need if you'd like to use the method.
*  MC_example.m shows how to setup the Monte-Carlo experiments done in the paper, and run the HSBQR on one of experiments (This code uses the files in the MC directory too)

-----
BibTeX citation:

```rb
@article{kohns2020horseshoe,
  title={Horseshoe prior Bayesian quantile regression},
  author={Kohns, David and Szendrei, Tibor},
  journal={arXiv preprint arXiv:2006.07655},
  year={2020}
}
```
