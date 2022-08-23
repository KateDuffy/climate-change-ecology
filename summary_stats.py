import cftime
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pwlf
import scipy.integrate as inte
import xarray as xr

def percent_change(n1, n2):
    return (n2 - n1)/n1 * 100
    return (n2 - n1)/n1 * 100


TPC_data = pd.read_csv("SD1.csv")


for scenario in ["usual", "rm_mean_trend", "rm_sd_trend", "rm_mean_sd_trend"]:
    
    df = pd.DataFrame(columns=["species", "model", "lat", "dR", "dN", "dCV", "dSD", "dStability", "r_fut", "r_hist", "N_fut", "N_hist"])
    
    
    for species in np.unique(TPC_data.Species):
        print(species)
        files = glob.glob("./pop_ts/*_%s.nc" % species)
        
        for file in files:
            
            temp = xr.open_dataset(file)
            model = file.split("/")[-1].split("_")[0]
            
            r_hist = temp["r_hist_" + scenario].values
            r_fut = temp["r_fut_" + scenario].values
            
            N_hist = temp["N_hist_" + scenario].values
            N_fut = temp["N_fut_" + scenario].values
            
            if np.sum(N_hist < 1e-5) > 0.:
                dN, dSD, dCV, dS = np.nan, np.nan, np.nan, np.nan
            else:
                dR = percent_change(np.nanmean(r_hist), np.nanmean(r_fut))
                dN = percent_change(np.nanmean(N_hist), np.nanmean(N_fut))
                dSD = percent_change(np.nanstd(N_hist), np.nanstd(N_fut))
                dCV =  percent_change(np.nanstd(N_hist)/np.nanmean(N_hist), np.nanstd(N_fut)/np.nanmean(N_fut))
                dS =  percent_change(np.nanmean(N_hist)/np.nanstd(N_hist), np.nanmean(N_fut)/np.nanstd(N_fut))
            
            row = TPC_data.loc[TPC_data['Species'] == species]
            df = df.append({"species":species, "model":model, "lat":row.Lat.values[0], "dR":dR, "dN":dN,
                           "dSD":dSD, "dCV":dCV, "dStability":dS, "r_hist":r_hist, "r_fut": r_fut,
                           "N_hist":N_hist, "N_fut": N_fut}, ignore_index=True)


    df.to_pickle("./stats/%s.pkl" %(scenario))
