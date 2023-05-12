#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 11 14:31:23 2023

@author: bjerknes
"""

import h5py
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import numpy as np
import cmocean
import xarray as xr
import metpy
from metpy.units import units

# Abrir o arquivo HDF5 em modo de leitura
# arquivo_hdf5 = h5py.File('/home/bjerknes/Downloads/DADOS_GELO/2019_01/636598.anom.oisst-avhrr-v02r01.20190101-ice.oisst-avhrr-v02r01.20190116nc/636598.ice.oisst-avhrr-v02r01.20190114.nc', 'r')
arquivo_hdf5 = xr.open_dataset('/home/bjerknes/Downloads/DADOS_GELO/2019_01/636598.anom.oisst-avhrr-v02r01.20190101-ice.oisst-avhrr-v02r01.20190116nc/636598.ice.oisst-avhrr-v02r01.20190101.nc').metpy.parse_cf()
tsm = xr.open_dataset('/home/bjerknes/Downloads/DADOS_GELO/2019_01/636598.ice.oisst-avhrr-v02r01.20190117-sst.oisst-avhrr-v02r01.20190131.nc/636598.sst.oisst-avhrr-v02r01.20190101.nc').metpy.parse_cf()
tsm = tsm.assign_coords(dict(lon = (((tsm.lon.values + 180) % 360) - 180))).sortby('lon')
ERA5 = xr.open_dataset('/home/bjerknes/Downloads/adaptor.mars.internal-1683921620.7452216-30185-5-0e0f87d6-e17d-4c68-b03e-cef7a8f1f50b.nc').metpy.parse_cf()

lat_slice_1 = slice(-55.,-85.)
lon_slice_1 = slice(-180.,180.)

lat_slice_2 = slice(-85.,-55.)
lon_slice_2 = slice(-180.,180.)

lats = ERA5.latitude.sel(latitude=lat_slice_1).values
lons = ERA5.longitude.sel(longitude=lon_slice_1).values

lats_2 = tsm.lat.sel(lat=lat_slice_2).values
lons_2 = tsm.lon.sel(lon=lon_slice_2).values

# Acessar as variáveis necessárias
temp_2m = ERA5['t2m'][1,:,:]-273.15
sst = tsm['sst'][1,::-1,:]
diff = temp_2m - sst

concentracao_gelo = arquivo_hdf5['ice']
latitude = arquivo_hdf5['lat'][:]
longitude = arquivo_hdf5['lon'][:]
# concentracao_gelo = np.ma.masked_equal(concentracao_gelo, 9.96921e+36)
# Reduzir a dimensão dos dados de concentração de gelo
concentracao_gelo = np.squeeze(concentracao_gelo)

#lat e lon extend da imagem

latmin, latmax = -90, -50
lonmin, lonmax = -95, -15        

# Criar uma figura
fig,ax = plt.subplots(figsize=(15,10),subplot_kw=dict(projection=ccrs.PlateCarree()))
args = dict(color='gray',
            alpha=1.0, 
            linestyle='--', 
            linewidth=0.5,
            xlocs=np.arange(lonmin, lonmax, 10), 
            ylocs=np.arange(latmin, latmax, 10), 
            draw_labels=True)

#ax = plt.axes(projection=ccrs.PlateCarree())
gl = ax.gridlines(crs=ccrs.PlateCarree(), **args)
gl.xlabel_style = {'size': 15, 'color': 'black', 'rotation': 0}
gl.ylabel_style = {'size': 15, 'color': 'Gray', 'rotation': 0}
gl.top_labels = False
gl.right_labels = False
gl.bottom_labels = True
gl.rotate_labels = True
ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.COASTLINE, linewidth=1, edgecolor='black')
ax.add_feature(cfeature.OCEAN)   

# intevalos da concentracao_gelo
intervalo_min2 = 0
intervalo_max2 = 1
interval_2 = 0.05            
levels_2 = np.arange(intervalo_min2, intervalo_max2, interval_2)

# # intevalos da sst
# intervalo_min3 = -2
# intervalo_max3 = 14
# interval_3 = 0.5            
# levels_3 = np.arange(intervalo_min3, intervalo_max3, interval_3)
 
# Criar um mapa de calor com a concentração de gelo marinho
im = ax.contourf(longitude, latitude, concentracao_gelo, cmap='Blues_r', levels = levels_2, extend = 'max', zorder=2)
im2 = ax.contourf(longitude, latitude, diff, cmap='jet', extend = 'neither', zorder=1)

# adiciona legenda 
barra_de_cores = plt.colorbar(im, 
                              orientation = 'horizontal', 
                              pad=0.05, 
                              fraction=0.05
                              )

barra_de_cores_2 = plt.colorbar(im2, 
                              orientation = 'vertical', 
                              pad=0.05, 
                              fraction=0.05
                              )

barra_de_cores.set_label('(%)', fontsize=20,  loc='center')
barra_de_cores_2.set_label('(°C)', fontsize=20,)
barra_de_cores_2.ax.tick_params(labelsize=font_size)

font_size = 20 # Adjust as appropriate.
barra_de_cores.ax.tick_params(labelsize=font_size)

# Adicionar títulos e rótulos dos eixos
ax.set_title('Concentração de Gelo Marinho\n'
             'DEZEMBRO - 2019', fontweight='bold', 
             fontsize=23)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Mostrar o mapa de calor
plt.show()
# plt.savefig(args, kwargs)
