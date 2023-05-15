#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 10:23:05 2023

@author: everson

Esse script plota todas as trajetória dentro 
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader # Import shapefiles
import imageio.v2 as imageio

# dataset
datafiles = glob.glob("/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Antartica/siac_2023/Dados_usados/Track_mensal_2019/2019-01.csv")
datafiles.sort(reverse=True)

lonmin, lonmax = -130, 0
latmin, latmax = -45, -90

# fig = plt.figure(figsize=(10,6))
# ax = fig.add_subplot(111)
#to read individual data files containing the coordinates of the track for each typhoon
fig = plt.figure(figsize=(17, 15), dpi=100)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

args = dict(draw_labels=True,
            linewidth=1, color='gray', 
            alpha=0.5, linestyle='--', 
            xlocs=np.arange(-180, 180, 10), 
            ylocs=np.arange(-90, -50, 5))


gl = ax.gridlines(crs=ccrs.PlateCarree(), **args)
gl.xlabel_style = {'size': 15, 'color': 'black', 'rotation': 0}
gl.ylabel_style = {'size': 15, 'color': 'black', 'rotation': 0}
gl.top_labels = False
gl.right_labels = False
gl.bottom_labels = True
gl.rotate_labels = True
ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND)
# ax.add_feature(cfeature.OCEAN)
ax.coastlines(resolution='50m', color='black', linewidth=1)
# ax.add_feature(cfeature.OCEAN)

shapefile = list(
    shpreader.Reader(
    '/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Coisas/br_unidades_da_federacao/BR_UF_2019.shp'
    ).geometries()
    )
ax.add_geometries(
    shapefile, ccrs.PlateCarree(), 
    edgecolor = 'black', 
    facecolor='none', 
    linewidth=0.5
    )

# colors_0 = ['black','red','blue','orange','gray','cyan','orangered','green','magenta','midnightblue', 'navy', 'purple', 'brown', 'lime']
# colors_1 = ['red', 'black', 'dimgray', 'dimgrey', 'gray', 'grey', 'darkgray', 'darkgray', 'darkgrey', 'darkgrey', 'silver', 'silver', 'lightgray', 'lightgray']
 #default colors from Python (can be automated if the order is not important)
# mark = ['o', 'v', 's', '^', '*', 'D', '+', 'x', '1', '4','P', '<','>', 'd']

# Quantidade de ciclones para plotar
numero = 12
images = []
for jj in range(len(datafiles)):
    #abre uma planilha por vez
    print(datafiles[jj])
    dff = pd.read_csv(datafiles[jj],sep=',', dtype={'Data': object}) #read data file and time as string
    
    #converte a longitude de 0 - 360 para -180 a 180
    #dff['Longitude'] = (((dff['Longitude'] + 180) % 360) - 180)
    
    # arredonda para resolução de 0.25 graus
    # dff['Latitude'] = round(dff.Latitude / 0.25) * 0.25
    # dff['Longitude'] = round(dff.Longitude / 0.25) * 0.25

    # rodada = datafiles[jj].split('/')[12].split('.')[0] #extract rodada information from the filename
    nome = datafiles[jj].split('/')
    rodada = nome[-1].split('.')[0]
    # track_name = datafiles[0].split('/')[8] #extract track information from the filename
    
    # rodada = datafiles[jj]
    # track_name = datafiles[jj]
    # track_name = track_name.capitalize()
    
    lons = dff['Longitude'].values
    lats = dff['Latitude'].values
   
    # ax.plot(lons, lats,color=colors_1[jj],linestyle=':', marker='.', linewidth=2.2 , ms=8, zorder=7, label=f"TRACK {jj} ({rodada})")
    ax.plot(lons, lats, color='black', linestyle='--', marker='.', linewidth=2.2 , ms=8, zorder=7)
    ax.set_title(f"Tracks: {(rodada)} ", 
                 fontweight='bold', 
                 fontsize=23)
    # ax.legend(loc='upper right', fontsize=18)
    
    plt.savefig(f'/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Antartica/imagens/Trajetorias/trajetória_{(rodada)}.png', bbox_inches='tight')

    # # Salvar a figura atual em um arquivo de imagem
    # filename = f'/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Antartica/imagens/Trajetorias/trajetória_{(rodada)}.png'
    # plt.savefig(filename, dpi=100, bbox_inches='tight')
    
    # # Adicionar o arquivo de imagem à lista de imagens
    # images.append(imageio.imread(filename))
    
    # # Salvar a lista de imagens como um arquivo .gif
    # filename = '/media/bjerknes/HD_todo_pod/Everson/Coqueiro/CENPES/Imagens/trajetória_Yakecan_GFSxERA5.gif'
    # imageio.mimsave(filename, images, format='GIF', duration=0.5, subrectangles=True, subrectangles_palette=256, palettesize=256, quantizer='nq') # duration em ms (milisegundos)
    



    