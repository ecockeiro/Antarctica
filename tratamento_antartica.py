#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 20:20:49 2023

@author: everson
"""

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.feature as cfeature
import numpy as np

#dataset
df = pd.read_csv('/home/everson/Siac/Siac/Dados/Track_validado_2010_2019/dados/2010.csv', sep=',')

#lat e lon extend da imagem
latmin, latmax = -90, -53
lonmin, lonmax = -120, -30

#converte e recorta latitude/longitude
df['Longitude'] = (((df['Longitude']+180) %360) -180)
df = df[(df.Longitude>=-120) & (df.Longitude<=-30) & (df.Latitude<=-55) & (df.Latitude>=-73)].copy()

# Conta a ocorrência de cada item
counts = df['N_do_Ciclone'].value_counts()

# Cria uma lista com os itens a serem removidos
itens_para_remover = counts[counts < 3].index.tolist()

# Remove os itens da lista do DataFrame
df = df[~df['N_do_Ciclone'].isin(itens_para_remover)]

# arredonda para resolução de 0.25 graus
df['Latitude'] = round(df.Latitude / 0.25) * 0.25
df['Longitude'] = round(df.Longitude / 0.25) * 0.25

#salva a planilha recortada
df.to_csv('/home/everson/Siac/Siac/Dados/Track_validado_2010_2019/2010_nova_planilha.csv', index='False')    

#escolhe o tamanho da figura e polegadas (x,y)
plt.figure(figsize=(15,10))

args = dict(color='gray',
            alpha=1.0, 
            linestyle='--', 
            linewidth=0.5,
            xlocs=np.arange(-180, 180, 10), 
            ylocs=np.arange(-90, -45, 10), 
            draw_labels=True)

# Criação da projeção polar centrada na Antártica

ax = plt.axes(projection=ccrs.SouthPolarStereo())

gl = ax.gridlines(crs=ccrs.PlateCarree(), **args)
gl.xlabel_style = {'size': 15, 'color': 'black', 'rotation': 0}
gl.ylabel_style = {'size': 15, 'color': 'Gray', 'rotation': 0}
gl.top_labels = False
gl.right_labels = False
gl.bottom_labels = True
gl.rotate_labels = True
ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND)
ax.coastlines(resolution='50m', color='black', linewidth=1)

#one_ciclone = data[(data.N_do_Ciclone == 200)].copy()
def plot_cyclones(ciclone):

    for i in range(len(ciclone)-1):
        #print(range(len(ciclone)))
        #print(i)
        ciclone1 = ciclone.iloc[i]
        ciclone2 = ciclone.iloc[i+1]
        if i == 0:
            frist_lat, frist_lon = ciclone1.Latitude, ciclone1.Longitude
            second_lat, second_lon = ciclone2.Latitude, ciclone2.Longitude
            plt.plot([frist_lon, second_lon], [frist_lat, second_lat],
                     color='black',linewidth=0.4,
                     transform=ccrs.PlateCarree()
                     )
        else:
            if i == len(ciclone)-2:
                frist_lat, frist_lon = ciclone1.Latitude, ciclone1.Longitude
                second_lat, second_lon = ciclone2.Latitude, ciclone2.Longitude
                plt.plot([frist_lon, second_lon], [frist_lat, second_lat],
                     color='black',linewidth=0.4,
                     transform=ccrs.PlateCarree()
                     )
            else:
                frist_lat, frist_lon = ciclone1.Latitude, ciclone1.Longitude
                second_lat, second_lon = ciclone2.Latitude, ciclone2.Longitude
                plt.plot([frist_lon, second_lon], [frist_lat, second_lat],
                     color='black',linewidth=0.4,
                     transform=ccrs.PlateCarree()
                     )



ciclones = list( dict.fromkeys(df.N_do_Ciclone) )

for Ciclone in ciclones:
    one_ciclone = df[(df.N_do_Ciclone == Ciclone)].copy()
    plot_cyclones(one_ciclone)
    
    
    
plt.show()