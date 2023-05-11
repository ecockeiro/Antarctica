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

### DATASETS

# ERA5
ERA5 = xr.open_dataset('/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Antartica/siac_2023/Dados_usados/ERA5/flux_t2m.nc').metpy.parse_cf()
temp_2m = ERA5.t2m
fluxo_sens = ERA5.sshf
fluxo_late = ERA5.slhf

# SST
SST =xr.open_mfdataset('/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Antartica/siac_2023/Dados_usados/SST/SST_2019/*.nc4').metpy.parse_cf()
SST = SST.assign_coords(dict(lon = (((SST.lon.values + 180) % 360) - 180))).sortby('lon')

# Concentraçao de gelo
CGELO = xr.open_mfdataset('/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Antartica/siac_2023/Dados_usados/Concentracao_gelo/mensal_2019/*.nc').metpy.parse_cf()
concentracao_gelo = CGELO.seaice_conc_monthly_cdr
latitude = CGELO.latitude
longitude = CGELO.longitude

# Extensão (Extent)
lat_slice = slice(-45.,-90.)
lon_slice = slice(-180.,180.)

#pega as lat/lon
lats = ERA5.latitude.sel(latitude=lat_slice).values
lons = ERA5.longitude.sel(longitude=lon_slice).values

for i in range(len(file_4['time'])):
    args_1 = dict(
        time = file_4.time[i],
        lat=lat_slice,
        lon=lon_slice
        )

    # flux_hl = file_1.mslhf.metpy.sel(**args_1).metpy.unit_array.squeeze()
    # flux_sl = file_2.msshf.metpy.sel(**args_1).metpy.unit_array.squeeze()
    # sic = file_3.siconc.metpy.sel(**args_1).metpy.unit_array.squeeze()
    sst = np.array(file_4.sst.sel(**args_1))
    sst= sst[~np.isnan(sst)]
    # temp2m =file_5.t2m.metpy.sel(**args_1).metpy.unit_array.squeeze().to('degC')
    
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
