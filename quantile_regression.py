import cftime
import dask.array as da
import glob
import itertools
import numpy as np
import os
import platform
import sys
from tqdm import tqdm
import xarray as xr
import time

from rpy2.robjects.packages import importr
base = importr('base')
utils = importr('utils')
#utils.install_packages('quantreg')
quantreg = importr('quantreg')
from rpy2.robjects import Formula
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
from rpy2 import robjects as ro
R = ro.r


def analyze_qr(data_path, model, ssp, save_dir):


    '''
    Perform quantile regression on temperature data

    Parameters
    ----------
    data_path: str
        Location of netcdf dataset(s) to analyze
    model: str
        Name of CMIP6 model
    ssp: str
        Name of experiment, e.g. "ssp585"
    save_dir: str
        Location to save interpolated dataset

    
    Returns
    ----------
    ds: xr.Dataset
        Quantile regression results - slope and confidence bounds for each tau value at each geographic location
    '''        

    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    savefile = os.path.join(save_dir, 'qr_%s_%s.nc' % (model, ssp))

    ds = xr.open_mfdataset(paths=data_path, combine='by_coords')

    data = ds['tas']
    ds.close()
    data = data.chunk({"lat":1, "lon":14, "time":91250})
    
    slope_results, lower_bd, upper_bd = return_results(data)
    print('Saving ...')
    write_results(data, slope_results, lower_bd, upper_bd, savefile)



def return_results(data):
    
    taus = np.arange(0.1, 1, 0.1)
    taus = np.concatenate([np.array([0.025]), taus, np.array([0.975])], axis=0)
    alpha = 0.1

    slope_results = np.empty((data.shape[1], data.shape[2], len(taus)))
    slope_results[:] = np.nan
    lower_bd = np.empty((data.shape[1], data.shape[2], len(taus)))
    lower_bd[:] = np.nan
    upper_bd = np.empty((data.shape[1], data.shape[2], len(taus)))
    upper_bd[:] = np.nan

    
    start = time.time()
    data.load()
    end = time.time()
    print("loaded data ", end-start, " seconds")

    for i in tqdm(range(data.shape[1])):
        for j in tqdm(range(data.shape[2]-1)):
    
            try:
                formula = Formula('tas ~ time')
                env = formula.environment
                env['tas'] = base.as_vector(data[:, i, j].values)
                env['time'] = base.as_vector(np.arange(len(data[:, i, j])))
                
                info = R.plot(R.summary(quantreg.rq(formula, tau = taus, alpha = alpha)))
                info = np.asarray(info)
                    
                slope_results[i, j, :] = info[1,0,:-1]
                lower_bd[i, j, :] = info[1,1,:-1]
                upper_bd[i, j, :] = info[1,2,:-1]

            except:
                print("passed: ", i, j )

    return slope_results, lower_bd, upper_bd


def write_results(data, slope_results, lower_bd, upper_bd, savefile):

    lat = data.lat
    lon = data.lon
    taus = np.arange(0.1, 1, 0.1)
    taus = np.concatenate([np.array([0.025]), taus, np.array([0.975])], axis=0)

    ds = xr.Dataset({"slope": (('lat', 'lon', 'tau'), slope_results),
                    "lower_bd": (('lat', 'lon', 'tau'), lower_bd),
                    "upper_bd": (('lat', 'lon', 'tau'), upper_bd)},
                    coords={'lat': lat, 'lon': lon, 'tau': taus})
    ds.to_netcdf(savefile)

    print("Saved ...", savefile)




