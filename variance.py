import cftime
import dask.array as da
import glob
import itertools
import matplotlib.pyplot as plt
import numpy as np
import os
import platform
import sys
from tqdm import tqdm
import xarray as xr
import datetime



def analyze_variance(data_path, model, ssp, save_dir):


    '''
    Calculate variance at each geographic location over a 10 year window 
    
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

    ds = xr.open_mfdataset(paths=data_path, combine='by_coords',
                               chunks={"lat":179, "lon":359, "time":3650})
    data = ds['tas']

    for startyear in np.arange(1855, 2090, 10):

        var_results = return_variance(data, startyear)

        save_file = os.path.join(save_dir, 'variance_%s_%s_%s.nc' % (model, ssp, startyear))
        write_results(ds, var_results, save_file)



def return_variance(data, startyear):
    
    var_results = {}

    try:
        start = cftime.DatetimeNoLeap(startyear, 1  , 1)
        end = cftime.DatetimeNoLeap(startyear+9, 12, 31)
        temp = data.sel(time=slice(start, end))
    except:
        start = np.datetime64(datetime.datetime(startyear, 1, 1))
        end = np.datetime64(datetime.datetime(startyear+9, 12, 31))
        temp = data.sel(time=slice(start, end))

    var_results[str(startyear)] = np.nanvar(temp, axis=0)

    return var_results



def write_results(ds, var_results, save_file):
    
    vname = list(var_results())[0]
    lat = ds.lat
    lon = ds.lon

    ds = xr.Dataset({vname: (('lat', 'lon'), var_results[vname])}, coords={'lat': lat, 'lon': lon})
    ds.to_netcdf(save_file)
        



