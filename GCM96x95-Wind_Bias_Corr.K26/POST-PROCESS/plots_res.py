import os
import argparse
import numpy as np
import xarray as xr
import cmocean

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.util as cutil

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import matplotlib
matplotlib.use('Agg')

def make_plot(data,lon,lat,infos,output):
    # args
    title, cmap, norm, tfs = infos
    data = tfs(data)
    # figure
    plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, zorder=100, edgecolor='k')
    # color map
    pcm = ax.pcolormesh(lon, lat, data, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
    cbar = plt.colorbar(pcm, ax=ax, orientation='vertical', pad=0.05, shrink=0.5)
    plt.title(title)
    # write fig
    plt.savefig(output, bbox_inches='tight')
    plt.close()


def main(filepath, var_name, fig_name, infos, freq):

    # read files
    try:
        ds = xr.open_dataset(filepath)
    except:
        print(f'FAILED TO READ {filepath}')
        return

    print(f'Plotting {var_name}')

    # get fields
    lon = ds.lon.values
    lat = ds.lat.values
    var = getattr(ds,var_name)
    if len(var.shape) == 3:
        var_val = var.values[5,:,:]
    elif len(var.shape) == 4:
        var_val = var.values[5,0,:,:]

    # plot
    plotpath = fig_name + '_' + config +'_' + freq + '.png'
    make_plot(var_val,lon,lat,infos,plotpath)



if __name__=="__main__":

    # Config name
    # -----------
    config = 'R242_PYC'

    # snapshots
    # ---------
    # zonal wind U
    infos = [ 'U (m/s)' , cmocean.cm.balance , colors.Normalize(vmin=-20.0, vmax=20.0), lambda x: x ]
    main( filepath=config+'/ATM/Output/MO/'+config+'_20150101_20151231_1M_histmth.nc' , var_name='u10m' , fig_name='U', infos=infos , freq='y2016m06' )

    # meridional wind V
    infos = [ 'V (m/s)' , cmocean.cm.balance , colors.Normalize(vmin=-20.0, vmax=20.0), lambda x: x ]
    main( filepath=config+'/ATM/Output/MO/'+config+'_20150101_20151231_1M_histmth.nc' , var_name='v10m' , fig_name='V', infos=infos , freq='y2016m06' )
    
    # specific humidity q
    infos = [ 'q (-)' , cmocean.cm.algae , colors.Normalize(vmin=0.0, vmax=0.03), lambda x: x ]
    main( filepath=config+'/ATM/Output/MO/'+config+'_20150101_20151231_1M_histmth.nc' , var_name='q2m' , fig_name='q', infos=infos , freq='y2016m06' )
    
    # potential temperature theta
    infos = [ 'Pot. Temp. (K)' , cmocean.cm.thermal , colors.Normalize(vmin=200.0, vmax=350.0), lambda x: x ]
    main( filepath=config+'/ATM/Output/MO/'+config+'_20150101_20151231_1M_histmth.nc' , var_name='theta' , fig_name='Theta', infos=infos , freq='y2016m06' )

