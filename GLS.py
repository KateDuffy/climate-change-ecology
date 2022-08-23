import glob
import itertools
import numpy as np
np.warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
import xarray as xr
import os
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
from rpy2.robjects import Formula
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
stats = importr('stats')
base = importr('base')
nlme = importr('nlme')
gls = nlme.gls




def GLS(data_path, parameter, model, ssp, save_dir):


    '''
    Perform generalized least squares (GLS) regression at each geographic location 

    Parameters
    ----------
    data_path: str
        Location of netcdf dataset(s) to analyze
    parameter: str
        Parameter to perform GLS on, such as "variance" or "spectral exponent"
    model: str
        Name of CMIP6 model
    ssp: str
        Name of experiment, e.g. "ssp585"
    save_dir: str
        Location to save dataset

    
    Returns
    ----------
    ds: xr.Dataset
        GLS results - trend and p value at each geographic location

    ''' 

    ds = xr.open_mfdataset(data_path, combine='by_coords')

    years = list(np.arange(1855, 2090, 10))
    years_str = [str(y) for y in years]
    data = np.stack([ds[y].values for y in years_str])
    

    trends, pvals = np.zeros((data.shape[1], data.shape[2])), np.zeros((data.shape[1], data.shape[2]))
    
    for i in range(data.shape[1]):
        if i % 10 == 0:
            print('%s percent done' %(np.round(i/data.shape[1]*100, 3)))
        for j in range(data.shape[2]):
            
            ##### GLS Regression ######
            formula = Formula('exponent ~ time')
            formula2 = Formula('exponent ~ 1')

            env = formula.environment
            env['exponent'] = data[:,i,j]
            env['time'] = np.array(years)

            # specify correlation parameter corAR1(form ~ 1) and check that phi is not zero
            c = nlme.corAR1(form = Formula('~ 1'))
            fit = gls(formula, correlation=c, method="ML")
            fit2 = gls(formula2, correlation=c, method="ML")
                
            trends[i, j] = fit.rx('coefficients')[0][1]
            pvals[i, j] = ro.r.anova(fit, fit2)[-1][1]
            

    sig = np.zeros_like(pvals)
    sig[pvals>=0.05] = np.nan
    sig[pvals<0.05] = 1.
    sigtrends = sig*trends

    out = ds.copy()
    for k in list(out.keys()):
        out = out.drop(k)

    out["trend"] = (("lat", "lon"), trends)
    out["pval"] = (("lat", "lon"), pvals)
    out["sig_trend"] = (("lat", "lon"), sigtrends)
    
    filename = "%s_GLS_%s_%s.nc" %(parameter, model, ssp)
    out.to_netcdf(save_dir + filename)










