#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 17:59:35 2023

@author: everson
"""

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.feature as cfeature
import numpy as np
import glob 
import os

# função que lê as pastas em um diretório
def get_folders(directory):
    return [f for f in os.listdir(directory) if os.path.isdir(os.path.join(directory, f))]
# diretório

directories = '/home/everson/siac_2023/dados_10_anos_mensais/Estações'

# lista de pasta para serem selecionadas
folders = get_folders(directories)
folders.sort()

# range = ao numero de meses, poderia usar o range(len(datafiles))
for i in range(len(folders)):
    
    # Cria uma lista de DataFrames
    dfs = []
    
    # II - Adiciona os rótulos nas colunas, pois quando o track gera as planilhas geram planilhas sem rótulos
    
    #seleciona a planilha
    path_name = os.path.join(directories, folders[i])
    datafiles = glob.glob(f'{path_name}/*.csv')
    datafiles.sort()
    
    for j in  range(len(datafiles)):
    
        df = pd.read_csv(datafiles[j] , sep=',', usecols=[0, 1, 2, 3, 4])
        
        # # Carrega a primeira planilha como referência
        # df_ref = pd.read_csv(datafiles[0], sep=',', usecols=[0, 1, 2, 3, 4])
        
        
        
        estacao = folders[i]
            
        # Adiciona o DataFrame à lista de DataFrames
        dfs.append(df)
    
        # Concatena todos os DataFrames na lista em uma única planilha
        result = pd.concat(dfs, axis=0, ignore_index=True)
        
        # Salva o resultado em um arquivo de saída
        result.to_csv(f'/home/everson/siac_2023/dados_10_anos_mensais/Estações/{(estacao)}/10_{(estacao)}.csv', index=False)
        
        # plota os tracks
        # Caso queira plotar as saídas das planilhas descomente essa parte!
        
        
        
        
        
'''
        #lat e lon extend da imagem
        latmin, latmax = -90, -50
        lonmin, lonmax = -120, -30
        
        fig = plt.figure(figsize=(10,5))
        args = dict(color='gray',
                    alpha=1.0, 
                    linestyle='--', 
                    linewidth=0.5,
                    xlocs=np.arange(lonmin, lonmax, 10), 
                    ylocs=np.arange(latmin, latmax, 10), 
                    draw_labels=True)
        
        
        # Criação da projeção cartesiana

        ax = plt.axes(projection=ccrs.PlateCarree())
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
        
        
        # Criação da projeção polar centrada na Antártica
        
        # ax = plt.axes(projection=ccrs.SouthPolarStereo())

        # gl = ax.gridlines(crs=ccrs.PlateCarree(), **args)
        # gl.xlabel_style = {'size': 15, 'color': 'black', 'rotation': 0}
        # gl.ylabel_style = {'size': 15, 'color': 'Gray', 'rotation': 0}
        # gl.top_labels = False
        # gl.right_labels = False
        # gl.bottom_labels = True
        # gl.rotate_labels = True
        # ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
        # ax.add_feature(cfeature.LAND)
        # ax.coastlines(resolution='50m', color='black', linewidth=1)
        
        
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
                              color='red',linewidth=0.5,
                              transform=ccrs.PlateCarree()
                              )
                else:
                    if i == len(ciclone)-2:
                        frist_lat, frist_lon = ciclone1.Latitude, ciclone1.Longitude
                        second_lat, second_lon = ciclone2.Latitude, ciclone2.Longitude
                        plt.plot([frist_lon, second_lon], [frist_lat, second_lat],
                              color='black',linewidth=0.5,
                              transform=ccrs.PlateCarree()
                              )
                    else:
                        frist_lat, frist_lon = ciclone1.Latitude, ciclone1.Longitude
                        second_lat, second_lon = ciclone2.Latitude, ciclone2.Longitude
                        plt.plot([frist_lon, second_lon], [frist_lat, second_lat],
                              color='blue',linewidth=0.5,
                              transform=ccrs.PlateCarree()
                              )



        ciclones = list( dict.fromkeys(result.N_do_Ciclone) )

        for Ciclone in ciclones:
            one_ciclone = result[(result.N_do_Ciclone == Ciclone)].copy()
            plot_cyclones(one_ciclone)
            
            
            
        plt.show()
        '''