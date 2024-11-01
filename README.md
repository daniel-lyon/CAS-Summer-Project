# CAS Summer Project

This repository contains code and results for the 2023/24 CAS Summer Project "Understanding Ultra Diffuse Galaxies With the Keck Telescopes" at Swinburne University of Technology.

## Introduction
15 Globular Cluster (GC) candidates are analysed to determine if they are true GCs or not. To do this, the python package `pPXF` is used to fit stellar templates to our spectra. First, the recessional velocity is fit to determine if the candidates are associated with their host galaxy. If candidates are not, then they are very likely not GCs. 3/15 candidates were not associated - 5P, 10P, & IKN - being two AGNs and a background galaxy respectively. 12/15 candidates were associated with thier host galaxy and so were each fit again with 256 input parameters combinations (15 kinematic and 15 multiplicated polynomials). The median age and metallicity from the 256 fits are taken to be the results with the 16th and 84th percentiles being the lower and upper error bounds respectively. 

8 GC candidates were confirmed - GC2,4,5,6, B336, H12, PA-41, & Sextans.
2 GC candidates were too young - GC1,3.
2 GC candidates fits could not be adequately determined - DDO190, F8D1.

## Code
All code required to analyse the data is located in `/scripts/`.

`cube_clean_n_stack.ipynb` for cleaning cubes and stacking spectra.

`ppxf_fit.ipynb` for checking recessional velocities.

`Redshift_Spectra_Display.ipynb` for visualising AGN/galaxy spectra.

`mp_stellar_pops.py` for running 256 parameter inputs to get age & metallicity.

`fitting_results_ppxf_explore.py` to visualse 256 parameter results.

Other files include some helper scripts and plotting scripts for paper

`agn_plot_fit.py` for the AGN & galaxy spectra

`spec_plot_fit.py` for the NGC247 & other candidates

## Velocity Fits
Candidates are checked to ensure their velocity aligns with the host galaxy. Results of this are in `/results_GC/`. This directory also contains files that were cleaned and stacked with `cube_clean_n_stack.ipynb` to create single spectra files the fitting was done with. All results of the fitting and input parameters are stored in `results_GC.csv`. Final fits are stored in `/results_GC/fits/`.

## Stellar Pops
Stellar population fits run with `mp_stellar_pops.py` are located in `/results_SP/` and can be visualised with `fitting_results_ppxf_explore.py`. Results and input parameter values are stored in `/results.csv`.

## Images
`/images/` contains images of the candidates. Most are from HST: GC1-6, DDO190, F8D1, PA-41, & B336. H12 image is from the SDSS DR14 Finding Chart. 5P, 10P, IKN, & Sextans are from Legacy Sky Survey DR9.

## Cite
if you used this data in your paper, please be sure to cite us with your favourite reference manager from [NASA ADS](https://ui.adsabs.harvard.edu/abs/2024PASA...41...44F/abstract) or as "Forbes, D. A., Lyon, D., Gannon, J., Romanowsky, A. J., & Brodie, J. P. 2024, Publications of the Astronomical Society of Australia, 41, e044, doi: 10.1017/pasa.2024.41"

## Contact
If would like to contact me, feel free to email anytime: `daniellyon31@gmail.com`
