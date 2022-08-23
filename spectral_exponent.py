import cftime
import dask.array as da
import glob
import itertools
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import os
import platform
import sys
from tqdm import tqdm
import xarray as xr

from utils import *


def analyze_SE(data_path, model, ssp, save_dir):


    '''
    Calculate the spectral exponent at each geographic location over a 10 year window 
    
    Parameters
    ----------
    data_path: str
        Location of netcdf dataset(s) to analyze
    model: str
        Name of CMIP6 model
    ssp: str
        Name of experiment, e.g. "ssp585"
    save_dir: str
        Location to save dataset

    
    Returns
    ----------
    ds: xr.Dataset
        Results 

    '''     
    
    
    ds = xr.open_mfdataset(paths=data_path, combine='by_coords',chunks={"time":3650})
    data = ds['tas_detrended']

        for startyear in np.arange(1855, 2090, 10):

        save_file = os.path.join(save_dir, 'spectral_exponent_%s_%s.nc' % (model, ssp, startyear))
        spectral_exponents, end = return_spectral_exponents(data, startyear)
        write_spectral_exponents(ds, spectral_exponents, save_file)

    ds.close()


def return_spectral_exponents(data, startyear):
    '''
    spectral analysis by FFT and calculation of spectral exponent

    '''

    spectral_exponents = {}

    start = cftime.DatetimeNoLeap(startyear, 1, 1)
    end = cftime.DatetimeNoLeap(startyear+9, 12, 31)
    temp = data.sel(time=slice(start, end))
    data.close()
    temp.load()
    print('time period: ', start.year, end.year)

    temp_spectral_exponents = np.empty((temp.shape[1], temp.shape[2]))
    temp_spectral_exponents[:] = np.nan

    for i in tqdm(range(temp.shape[1])):
        for j in range(temp.shape[2]-1):
    
            try:
                timeseries = temp[:, i, j]

                A = np.fft.fft(timeseries)
                power = np.abs(A) ** 2
                n = timeseries.size
                freq = np.fft.fftfreq(n)

                power = power[1:int(n/2)]
                freq = freq[1:int(n/2)]

                fit = np.polyfit(np.log10(freq), np.log10(power), deg=1)

                temp_spectral_exponents[i, j] = fit[0]
            except:
                print("passed ", i, j)


    # write to dictionary
    spectral_exponents[str(start.year)] = temp_spectral_exponents

    return spectral_exponents, end




def write_spectral_exponents(ds, spectral_exponents, save_file):
    
    vname = list(spectral_exponents.keys())[0]
    lat = ds.lat
    lon = ds.lon
    ds_out = xr.Dataset({vname: (('lat' , 'lon'), spectral_exponents[vname])},
                    coords={'lat': lat, 'lon': lon})
    ds_out.to_netcdf(save_file)

    ds_out.close()
    ds_out = None





