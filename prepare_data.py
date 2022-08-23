import cftime
import dask
import dask.array as da
from dask.distributed import Client, progress, as_completed
import datetime as dt
import glob
import itertools
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pwlf
import scipy
from scipy import interpolate
from tqdm import tqdm
import xarray as xr

from utils import *



def interpolate_cmip(file, save_dir='CMIP6_tas_interp', vname='tas'):

    '''
    Interpolate dataset to 1 degree resolution
    
    Parameters
    ----------
    file: str
        Location of netcdf dataset
    save_dir: str
        Location to save interpolated dataset
    vname: str
        Name of variable in dataset
    
    Returns
    ----------
    out_ds: xr.Dataset
        Interoplated dataset
        
    '''   

    save_file = os.path.join(save_dir, 'interp_' + file.split('/')[-1])

    print('Interpolating...')
    dataset = xr.open_mfdataset(file, combine="by_coords")

    # original grid
    x, y, z = dataset.lon, dataset.lat, dataset[vname].data
    X, Y = np.meshgrid(x, y)

    # reference grid
    x2, y2 = np.linspace(1, 359, 359), np.linspace(-89, 89, 179) # 1 degree
    X2, Y2 = np.meshgrid(x2, y2)

    interp_temp = np.zeros((len(dataset.time), len(y2), len(x2)))

    for t in tqdm(range(z.shape[0])):
        z_temp = z[t, :, :]
        out = scipy.interpolate.griddata((X.flatten(), Y.flatten()), z_temp.flatten(),
                                         (X2.flatten(), Y2.flatten()), 'linear').reshape(X2.shape)
        interp_temp[t, :, :] = out

    interp_da = da.from_array(interp_temp)
    out_ds = xr.Dataset(data_vars={vname: (('time', 'lat', 'lon'), interp_da)},
                           coords={'time': dataset.time.data, 'lat': y2.data, 'lon': x2.data})


    print('Saving ', save_file)
    out_ds.to_netcdf(path=save_file)



def linear_detrend_cmip(files, save_dir):

    '''
    Concatenate netcdf datasets along time dimension

    Remove values over 60 deg C (140 deg F)

    Perform piecewise linear detrending

    
    Parameters
    ----------
    files: str
        List of directory locations of netcdf datasets belonging to one climate model/experiment
    save_dir: str
        Location to save detrended dataset

    
    Returns
    ----------
    ds: xr.Dataset
        Detrended dataset
        
    '''   


    ds = xr.open_mfdataset(paths= files, combine='by_coords')

    try:
        x = cftime.JulianDayFromDate(ds.time)
    except: # coerce time to cftime.DatetimeNoLeap
        timestring = [pd.to_datetime(str(v)).strftime('%Y-%m-%d') for v in pd.to_datetime(ds.time.values)]
        dropleap = [x.split('-') for x in timestring if not '-02-29' in x]
        converted = [cftime.DatetimeNoLeap(int(s[0]), int(s[1]), int(s[2]), 12, 0, 0, 0) for s in dropleap]

        ds = ds.sel(time=~((ds.time.dt.month == 2) & (ds.time.dt.day == 29)))
        ds['time'] = converted
        x = cftime.JulianDayFromDate(ds.time)

    # piecewise linear detrending
    slope, intercept = best_fit(x, ds.tas.data)
    x_pred = np.expand_dims(np.expand_dims(x, -1), -1)
    prediction = x_pred * slope + intercept
    ds_detrended = ds.tas - prediction

    ds = ds.drop('tas')
    ds['tas_detrended'] = ds_detrended.astype('float32')


    save_file = os.path.join(save_dir, 'detrended_' + files[0].split('/')[-1])
    print('Saving ', save_file)
    ds.to_netcdf(path=save_file)



