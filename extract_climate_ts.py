import glob
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import os
import tqdm


def extract_climatology(data_path, model, ssp, save_dir):


    '''
    Extract temperature time series at the geographic location of each species in dataset
    
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
        Dataset with temperature time series at lat/lon location of each species 

    '''        
    
    df = pd.read_csv('SD1.csv')

    files = glob.glob('%s*%s*%s*' %(data_dir, model, ssp)) + glob.glob('%s*%s*hist*' %(data_dir, model, 'hist'))
    ds_clim = xr.open_mfdataset(files, combine="by_coords")
    
    for index, row in tqdm.tqdm(df.iterrows()):
        species = row.Species
        lat = row.Lat
        long = row.Long % 360
        
        outfile = os.path.join(save_dir, "%s_%s.nc"%(model, row.Species))
	
        if not os.path.exists(outfile):

            temp = ds_clim.sel(lat=lat, method="nearest").sel(lon=long, method="nearest")
            out = xr.Dataset({'tas':(['time'], temp.tas.values)}, coords={'time': temp.time.values})
            out.to_netcdf(outfile)
            out.close()
            out = None


