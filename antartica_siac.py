#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 18:52:43 2023

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
from scipy.interpolate import griddata

### DATASETS

# Concetraçao de gelo
CGELO = xr.open_mfdataset('/home/everson/Documentos/ssd_antigo/maq_virtual/PIBIC-Antarctica/Dados_usados/Concentracao_gelo/mensal_2019/*.nc').metpy.parse_cf()

# TSM
TSM = xr.open_mfdataset('/home/everson/Documentos/ssd_antigo/maq_virtual/PIBIC-Antarctica/Dados_usados/SST/SST_2019/*.nc4').metpy.parse_cf()
TSM = TSM.assign_coords(dict(lon = (((TSM.lon.values + 180) % 360) - 180))).sortby('lon')

# ERA5
ERA5 = xr.open_dataset('/home/everson/Documentos/ssd_antigo/maq_virtual/PIBIC-Antarctica/Dados_usados/ERA5/flux_t2m.nc').metpy.parse_cf()


# Extensão (Extent)
lat_slice = slice(-45.,-90.)
lon_slice = slice(-180.,180.)

#pega as lat/lon
lats = ERA5.latitude.sel(latitude=lat_slice).values
lons = ERA5.longitude.sel(longitude=lon_slice).values

for i in range(len(ERA5['time'])):
    args_1 = dict(longitude=lon_slice, latitude=lat_slice, time = ERA5.time[i])


    flux_hl = ERA5.slhf.metpy.sel(**args_1).metpy.unit_array.to('kJ/m**2')
    flux_sl = ERA5.sshf.metpy.sel(**args_1).metpy.unit_array.to('kJ/m**2')
    
    sst = TSM.sst.sel(lon=lon_slice, time= TSM.time[i], lat=lat_slice).metpy.unit_array
    temp2m = ERA5.t2m.metpy.sel(**args_1).metpy.unit_array.to('degC')
    
    NLATS = np.linspace(lats.min(), lats.max(), temp2m.shape[0])
    NLONS = np.linspace(lons.min(), lons.max(), temp2m.shape[1])
    
    grid_lons, grid_lats = np.meshgrid(NLONS, NLATS)
    # Redimensionar sst para o novo shape usando interpolação
    sst_resized = griddata((lons.flatten(), lats.flatten()), sst.flatten(), (grid_lons, grid_lats), method='linear')


    # Redimensionar sst para o shape de temp2m
    sst_resized = mpinterp.interpolate_2d(sst, temp2m.shape)
    diff_temp = sst-temp2m
    
    vtempo = file_4.time.data[i].astype('datetime64[ms]').astype('O')

    ###################### Especificações do Plot #############################
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      color='gray',
                      alpha=1.0, 
                      linestyle='--', 
                      linewidth=0.5,
                      xlocs=np.arange(-180, 180, 10), 
                      ylocs=np.arange(-90, -45, 10), 
                      draw_labels=True)
    gl.xlabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.ylabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    gl.top_labels = True
    gl.right_labels = False
    gl.bottom_labels = False
    gl.rotate_labels = False
    ax.set_extent([-150, -10, -90, -55], ccrs.PlateCarree())
    
    # Intevalos da adv de temp
    intervalo_min_1 = -5
    intervalo_max_1 = 5.5
    interval_1 = 0.5             # de quanto em quanto voce quer que varie
    levels_1 = np.arange(intervalo_min_1, intervalo_max_1, interval_1)
    
    ax.add_feature(cfeature.LAND)
    ax.coastlines(resolution='10m', color='black', linewidth=2)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=2)
    
    # Plota a imagem adv de temp
    sst_sombreado = ax.contourf(lons, 
                            lats, 
                            sst, 
                            cmap='Spectral',
                            levels = levels_1,
                            extend = 'neither',
                            transform=ccrs.PlateCarree()
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
    
    # ###################### Especificações do Plot #############################
    # # escolha o tamanho do plot em polegadas (largura x altura)
    # plt.figure(figsize=(25,25))
    # ax = plt.axes(projection=ccrs.SouthPolarStereo())
    # gl = ax.gridlines(crs=ccrs.PlateCarree(),
    #                   color='gray',
    #                   alpha=1.0, 
    #                   linestyle='--', 
    #                   linewidth=0.5,
    #                   xlocs=np.arange(-180, 180, 10), 
    #                   ylocs=np.arange(-90, -45, 10), 
    #                   draw_labels=True)
    # gl.xlabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    # gl.ylabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
    # gl.top_labels = True
    # gl.right_labels = False
    # gl.bottom_labels = False
    # gl.rotate_labels = False
    # ax.set_extent([-150, -10, -90, -55], ccrs.PlateCarree())
    
    # # Intevalos da adv de temp
    # intervalo_min_1 = -40
    # intervalo_max_1 = 20
    # interval_1 = 2             # de quanto em quanto voce quer que varie
    # levels_1 = np.arange(intervalo_min_1, intervalo_max_1, interval_1)
    
    # ax.add_feature(cfeature.LAND)
    # ax.coastlines(resolution='10m', color='black', linewidth=2)
    # ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=2)
    
    # # Plota a imagem adv de temp
    # temp2m_sombreado = ax.contourf(lons, 
    #                         lats, 
    #                         temp2m, 
    #                         cmap='Spectral_r',
    #                         levels = levels_1,
    #                         extend = 'neither',
    #                         transform=ccrs.PlateCarree()
    #                         )
    
    # #adicionando shapefile
    # shapefile = list(
    #     shpreader.Reader(
    #     '/home/everson/Downloads/GFS-analysis_and_forecast-main/shapefiles/BR_UF_2021/BR_UF_2021.shp'
    #     ).geometries()
    #     )
    
    # ax.add_geometries(
    #     shapefile, ccrs.PlateCarree(), 
    #     edgecolor = 'black', 
    #     facecolor='none', 
    #     linewidth=0.5
    #     )
    
    # # adiciona continente e bordas
    # ax.add_feature(cfeature.LAND)
    # ax.coastlines(resolution='10m', color='black', linewidth=3)
    # ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=3)
    
    # # adiciona legenda 
    # barra_de_cores = plt.colorbar(temp2m_sombreado, 
    #                               orientation = 'horizontal', 
    #                               pad=0.04, 
    #                               fraction=0.04
    #                               )
    # font_size = 20 # Adjust as appropriate.
    # barra_de_cores.ax.tick_params(labelsize=font_size)
    
    # plt.title('Valid time: {}'.format(vtempo),
    #           fontweight='bold', 
    #           fontsize=35, 
    #           loc='left'
    #           )
    
    # plt.savefig(f'/home/everson/Siac/Siac/Dados/Track_validado_2010_2019/dados/imagens/temp_2m/temp_2m_{vtempo}.png', bbox_inches='tight')
    # plt.close()

    '''
        ###################### Especificações do Plot #############################
        # escolha o tamanho do plot em polegadas (largura x altura)
        plt.figure(figsize=(25,25))
        ax = plt.axes(projection=ccrs.SouthPolarStereo())
        gl = ax.gridlines(crs=ccrs.PlateCarree(),
                          color='gray',
                          alpha=1.0, 
                          linestyle='--', 
                          linewidth=0.5,
                          xlocs=np.arange(-180, 180, 10), 
                          ylocs=np.arange(-90, -45, 10), 
                          draw_labels=True)
        gl.xlabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
        gl.ylabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
        gl.top_labels = True
        gl.right_labels = False
        gl.bottom_labels = False
        gl.rotate_labels = False
        ax.set_extent([-150, -10, -90, -55], ccrs.PlateCarree())
        
        # Intevalos da adv de temp
        intervalo_min_1 = -10
        intervalo_max_1 = 10.8
        interval_1 = 0.8             # de quanto em quanto voce quer que varie
        levels_1 = np.arange(intervalo_min_1, intervalo_max_1, interval_1)
        
        ax.add_feature(cfeature.LAND)
        ax.coastlines(resolution='10m', color='black', linewidth=2)
        ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=2)
        
        # Plota a imagem adv de temp
        mslhf_sombreado = ax.contourf(lons, 
                                lats, 
                                flux_hl, 
                                cmap='seismic',
                                levels = levels_1,
                                extend = 'neither',
                                transform=ccrs.PlateCarree()
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
        barra_de_cores = plt.colorbar(mslhf_sombreado, 
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
        
        plt.savefig(f'/home/everson/Siac/Siac/Dados/Track_validado_2010_2019/dados/imagens/calor_latente/calor_latente_{vtempo}.png', bbox_inches='tight')
        plt.close()
        
        ###################### Especificações do Plot #############################
        # escolha o tamanho do plot em polegadas (largura x altura)
        plt.figure(figsize=(25,25))
        ax = plt.axes(projection=ccrs.SouthPolarStereo())
        gl = ax.gridlines(crs=ccrs.PlateCarree(),
                          color='gray',
                          alpha=1.0, 
                          linestyle='--', 
                          linewidth=0.5,
                          xlocs=np.arange(-180, 180, 10), 
                          ylocs=np.arange(-90, -45, 10), 
                          draw_labels=True)
        gl.xlabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
        gl.ylabel_style = {'size': 29, 'color': 'black', 'rotation': 0}
        gl.top_labels = True
        gl.right_labels = False
        gl.bottom_labels = False
        gl.rotate_labels = False
        ax.set_extent([-150, -10, -90, -55], ccrs.PlateCarree())
        
        # Intevalos da adv de temp
        intervalo_min_1 = -10
        intervalo_max_1 = 10.8
        interval_1 = 0.8             # de quanto em quanto voce quer que varie
        levels_1 = np.arange(intervalo_min_1, intervalo_max_1, interval_1)
        
        ax.add_feature(cfeature.LAND)
        ax.coastlines(resolution='10m', color='black', linewidth=2)
        ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=2)
        
        # Plota a imagem adv de temp
        msshf_sombreado = ax.contourf(lons, 
                                lats, 
                                flux_sl, 
                                cmap='BrBG',
                                levels = levels_1,
                                extend = 'neither',
                                transform=ccrs.PlateCarree()
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
        barra_de_cores = plt.colorbar(msshf_sombreado, 
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
        
        plt.savefig(f'/home/everson/Siac/Siac/Dados/Track_validado_2010_2019/dados/imagens/calor_sensivel/calor_sensivel_{vtempo}.png', bbox_inches='tight')
        plt.close()
    '''