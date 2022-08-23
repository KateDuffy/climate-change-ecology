import cftime
import csv
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scipy.integrate as inte
#from statsmodels.tsa.stattools import adfuller
import xarray as xrs



def run_simulation(s, data_dir, save_dir):

    '''
    Perform temperature-driven dynamical population modeling for a given species under projected temperature changes ("usual"), when controlling for trends in mean temperature ("rm_mean_trend"), when controlling for trends in temperature variability ("rm_sd_trend"), and when controlling for trends in both mean temperature and temperature variability ("rm_mean_sd_trend")

    Parameters
    ----------
    s: str
        Name of insect species
    data_dir: str
        Location of climate time series corresponding to the geographic location of species s
    save_dir: str
        Location to save dataset

    
    Returns
    ----------
    ds: xr.Dataset
        Population time series for the historical ("N_hist"; 1950-1999) and future time periods ("N_fut"; 2050-2099) and under projected and modified climate scenarios

    ''' 


    TPC_data = pd.read_csv("SD1.csv")
    files = glob.glob(data_dir+"*_%s.nc" %s)
    
    r = TPC_data.loc[TPC_data['Species'] == s]
    T_opt, sigma, T_max, lat, lon = r.Topt.item(), r.SIGMA.item(), r.CTmax.item(), r.Lat.item(), r.Long.item()
    
    
    for i, file in enumerate(files):
        ds = xr.open_mfdataset(file, combine="by_coords")
        ds = covert_to_cftime(ds)
        ds_out = ds.copy()
        
        for scenario in ["usual", "rm_mean_trend", "rm_sd_trend", "rm_mean_sd_trend"]:
            print(i, scenario)
            ds = xr.open_mfdataset(file, combine="by_coords")
            ds = covert_to_cftime(ds)
            
            tas_hist, tas_fut = detrend_time_series(ds, scenario)
            
            rs_hist = vfunc(tas_hist.tas.values, T_opt, sigma, T_max)
            rs_fut = vfunc(tas_fut.tas.values, T_opt, sigma, T_max)
            
            N_hist = run_r_alpha(rs_hist, theta=1.)
            N_fut = run_r_alpha(rs_fut, theta=1.)
            
            ds_out["N_hist_"+scenario] = N_hist
            ds_out["N_fut_"+scenario] = N_fut
        
            ds_out["r_hist_"+scenario] = rs_hist
            ds_out["r_fut_"+scenario] = rs_fut

        ds_out.to_netcdf(os.path.join(save_dir, file.split("/")[-1]))



def TPC_growth_rate(t, T_opt, sigma, T_max):
    t -= 273.15
    if t < T_opt:
        r = np.exp( - ( (t - T_opt) / (2 * sigma) ) ** 2)
    elif t >= T_opt:
        r = 1 - ( (t - T_opt)/(T_opt - T_max) ) ** 2
    return r


vfunc = np.vectorize(TPC_growth_rate)


def r_alpha(t, y, pars):
    N = y[0]
    if N < 1e-9:
        dN = 0.
    else:
        loc = np.int(np.floor(t))
        rs = pars['r']
        r_current = rs[loc] + (t-loc) * (rs[loc + 1] - rs[loc])
        dN = N * (r_current - pars['alpha'] * (N ** pars['theta'])) # dN/dt = N * (r - alpha * N^theta)
    return([dN])


def run_r_alpha(r, theta = 1., alpha = 1.):
    t_beg = 0
    t_end = len(r) - 2
    t_span = [t_beg, t_end]
    times = np.linspace(t_beg, t_end, len(r))
    pars = dict(r = r, alpha = alpha, theta = theta)
    N0 = np.max([0.01, 1/alpha * r[0]])
    init_conds = [N0]
    sol = inte.solve_ivp(lambda t,y: r_alpha(t, y, pars),
                         t_span, init_conds, method = 'RK45', t_eval = times)
    return sol.y[0, :]


def rm_mean_trend(ds, m):
    for i in range(5):
        data = ds.tas.values[i*365*10:(i+1)*365*10]
        data = m + (data - data.mean())
        ds.tas[i*365*10:(i+1)*365*10] = data
    return(ds)


def rm_sd_trend(ds, m, sd):
    for i in range(5):
        data = ds.tas.values[i*365*10:(i+1)*365*10]
        data = m + (data - m) * (sd/data.std())
        ds.tas[i*365*10:(i+1)*365*10] = data
    ds.tas[:] = m + (ds.tas.values - m) * (sd/ds.tas.values.std())
    return(ds)


def covert_to_cftime(ds):
    ds.load()
    try:
        x = cftime.JulianDayFromDate(ds.time)
    except: # coerce time to cftime.DatetimeNoLeap
        timestring = [pd.to_datetime(str(v)).strftime('%Y-%m-%d') for v in pd.to_datetime(ds.time.values)]
        dropleap = [x.split('-') for x in timestring if not '-02-29' in x]
        converted = [cftime.DatetimeNoLeap(int(s[0]), int(s[1]), int(s[2]), 12, 0, 0, 0) for s in dropleap]
        ds = ds.sel(time=~((ds.time.dt.month == 2) & (ds.time.dt.day == 29)))
        ds['time'] = converted
    return ds


def detrend_time_series(ds, scenario):
    tas_hist = ds.sel(time=slice(cftime.DatetimeNoLeap(1950, 1, 1, 12),
                                 cftime.DatetimeNoLeap(1999, 12, 31, 12)))
    tas_fut = ds.sel(time=slice(cftime.DatetimeNoLeap(2050, 1, 1, 12),
                                cftime.DatetimeNoLeap(2099, 12, 31, 12)))
    mean_hist, sdev_hist = np.nanmean(tas_hist.tas), np.nanstd(tas_hist.tas)
    mean_fut, sdev_fut = np.nanmean(tas_fut.tas), np.nanstd(tas_fut.tas)
                                 
    if scenario == "usual":
        pass
    elif scenario == "rm_mean_trend":
        tas_fut = rm_mean_trend(tas_fut, mean_hist)
    elif scenario == "rm_sd_trend":
        tas_fut = rm_sd_trend(tas_fut, mean_hist, sdev_hist)
    elif scenario == "rm_mean_sd_trend":
        tas_fut = rm_mean_trend(tas_fut, mean_hist)
        tas_fut = rm_sd_trend(tas_fut, mean_hist, sdev_hist)
                             
    return tas_hist, tas_fut

