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

# Abrir o arquivo HDF5 em modo de leitura
arquivo_hdf5 = h5py.File('/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Antartica/siac_2023/dados/dados_gelo/SEA_ICE_CONC_MONTHLY/2019/MENSAL/seaice_conc_monthly_sh_f17_201912_v03r01.nc', 'r')

# Acessar as variáveis necessárias
concentracao_gelo = arquivo_hdf5['seaice_conc_monthly_cdr'][:]
latitude = arquivo_hdf5['latitude'][:]
longitude = arquivo_hdf5['longitude'][:]
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
ax.add_feature(cfeature.COASTLINE, linewidth=1, edgecolor='black', zorder=3)
ax.add_feature(cfeature.OCEAN)   

# intevalos da concentracao_gelo
intervalo_min2 = 0
intervalo_max2 = np.amax(np.array(concentracao_gelo))
interval_2 = 5            
levels_2 = np.arange(intervalo_min2, intervalo_max2, interval_2)
 
# Criar um mapa de calor com a concentração de gelo marinho
im = ax.contourf(longitude, latitude, concentracao_gelo, cmap='Blues_r', levels = levels_2, extend = 'max', zorder=2)


# Adicionar uma barra de cores
cbar = fig.colorbar(im)

# Adicionar títulos e rótulos dos eixos
ax.set_title('Concentração de Gelo Marinho')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Mostrar o mapa de calor
plt.show()
