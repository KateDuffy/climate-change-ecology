# Climate-mediated shifts in temperature fluctuations promote extinction risk 

This repository contains code to perform the analyses and modeling described in the following paper

Duffy, Kate, Tarik C. Gouhier, and Auroop R. Ganguly. "Climate-mediated shifts in temperature fluctuations promote extinction risk." arXiv preprint arXiv:2202.12236 (2022).


## Abstract

Climate-mediated changes in the spatiotemporal distribution of thermal stress can destabilize animal populations and promote extinction risk. Assessment of risks typically uses thermal performance curves, however their reliance on mean conditions limits their sensitivity to temperature fluctuations at finer scales. Using projections from the latest generation of Earth System Models, we show that significant regional differences in the statistical distribution of temperature will emerge over time and give rise to shifts in the mean, variability, and persistence of thermal stress. Integrating these trends into mathematical models that simulate the dynamical and cumulative effects of thermal stress on the performance of 38 ectotherm species distributed around the globe revealed complex regional changes in population stability over the course of the 21st century, with temperate species facing higher risk. However, despite their idiosyncratic effects on stability, projected temperatures universally increase extinction risk. These results show that the effects of climate change may be more extensive than previously predicted based on the statistical relationship between biological performance and average temperature.


## Datasets

The CMIP6 simulation data used in this paper is available via the data portal [https://esgf-node.llnl.gov/search/cmip6/](https://esgf-node.llnl.gov/search/cmip6/). (Sample search constraints: Source ID = BCC-CSM2-MR, Experiment ID = ssp585, Frequency = day, Variable = tas)


The ectotherm thermal performance parameters are available for download as SD1.xls at [https://doi.org/10.1073/pnas.0709472105](https://doi.org/10.1073/pnas.0709472105)


## Getting started


### Climate data analysis

`prepare_data.py` contains functions to interpolate climate data to a 1 degree by 1 degree grid and perform detrending by piecewise linear fit

`quantile_regression.py` contains a function 'analyze_qr' to perform quantile regression with confidence bounds at each geographic location for a given range of tau values

`variance.py` contains a function 'analyze_variance' to calculate variance for the detrended temperature data at each geographic location for 10 year moving windows

`spectral_exponent.py` contains a function 'analyze_SE' to calculate the spectral for the detrended temperature dat at each geographic location for 10 year moving windows

`GLS.py` contains functions to perform generalized least squares (GLS) regression on a parameter (e.g. variance or spectral exponent) at each geographic location 



### Ecology modeling

`extract_climate_ts.py` contains a function 'extract_climatology' to extract and save the temperature time series at the geographic location of each insect species

`simulation.py` contains a function 'run_simulation' that simulates population growth using the r-alpha model and growths rates r determined using each species-specific thermal performance curve

`summary_stats.py` contains code to calculate percent change in performance for each population under a variety of climate scenarios and save the statistics
