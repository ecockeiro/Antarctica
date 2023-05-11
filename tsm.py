#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 22:41:04 2023

@author: everson
"""

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cmocean
import metpy.calc as mpcalc
from metpy.units import units
import cartopy.io.shapereader as shpreader # Import shapefiles
import matplotlib.colors as mcolors
import dask.array as da

# sst
file_4 = xr.open_mfdataset('/home/everson/siac_2023/dados/dados_gelo/TSM_ERSST.v5/ANUAL/2010/*.nc4')
file_4 = file_4.assign_coords({'lon' :((file_4.lon.values + 180) % 360) - 180}).sortby('lon')

# Extensão (Extent)
lat_slice = slice(-45.,-90.)
lon_slice = slice(-180.,-10.)

#pega as lat/lon
lats = file_4.lat.sel(lat=lat_slice).values
lons = file_4.lon.sel(lon=lon_slice).values

for i in range(len(file_4['time'])):
    args_1 = dict(
        time = file_4.time[i],
        lat=lat_slice,
        lon=lon_slice
        )

    # flux_hl = file_1.mslhf.metpy.sel(**args_1).metpy.unit_array.squeeze()
    # flux_sl = file_2.msshf.metpy.sel(**args_1).metpy.unit_array.squeeze()
    # sic = file_3.siconc.metpy.sel(**args_1).metpy.unit_array.squeeze()
    sst = file_4.sst.sel(**args_1).squeeze()
    # temp2m =file_5.t2m.metpy.sel(**args_1).metpy.unit_array.squeeze().to('degC')
    
    vtempo = file_4.time.data[i].astype('datetime64[ms]').astype('O')

    ###################### Especificações do Plot #############################
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(15,15))
    
    # usando a projeção da coordenada cilindrica equidistante 
    ax = plt.axes(projection=ccrs.PlateCarree())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      color='gray',
                      alpha=1.0, 
                      linestyle='--', 
                      linewidth=0.5, 
                      xlocs=np.arange(-180, 180, 10), 
                      ylocs=np.arange(-90, 90, 10), 
                      draw_labels=True
                      )
    gl.top_labels = False 
    gl.right_labels = False
    gl.xlabel_style = {'size': 29, 'color': 'black'}
    gl.ylabel_style = {'size': 29, 'color': 'black'}
    
    # Intevalos da adv de temp
    intervalo_min_1 = -1.8
    intervalo_max_1 = 10.8
    interval_1 = 0.8             # de quanto em quanto voce quer que varie
    levels_1 = np.arange(intervalo_min_1, intervalo_max_1, interval_1)
    
    ax.add_feature(cfeature.LAND)
    ax.coastlines(resolution='10m', color='black', linewidth=2)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=2)
    
    # Plota a imagem adv de temp
    sst_sombreado = ax.contourf(lons, 
                            lats, 
                            sst, 
                            cmap='Spectral_r',
                            levels = levels_1,
                            extend = 'neither',
                            )
    
    #adicionando shapefile
    shapefile = list(
        shpreader.Reader(
        '/home/everson/Downloads/GFS-analysis_and_forecast-main/shapefiles/BR_UF_2021/BR_UF_2021.shp'
        ).geometries()
        )
    
    ax.add_geometries(
        shapefile, ccrs.PlateCarree(), 
        edgecolor = 'black', 
        facecolor='none', 
        linewidth=0.5
        )
    
    # adiciona continente e bordas
    ax.add_feature(cfeature.LAND)
    ax.coastlines(resolution='10m', color='black', linewidth=3)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=3)
    
    # adiciona legenda 
    barra_de_cores = plt.colorbar(sst_sombreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
    plt.title('Valid time: {}'.format(vtempo),
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    plt.savefig(f'/home/everson/siac_2023/dados/imagens/sst/sst_{vtempo}.png', bbox_inches='tight')
    plt.close()
