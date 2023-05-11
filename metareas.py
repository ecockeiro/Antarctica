#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 13:16:47 2023

@author: ladsin
"""

import pandas as pd
import numpy as np
import os
import glob 
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import cartopy.feature as cfeature
import xarray as xr
# import regionmask as rmk



    
pasta1 = '/work/archive/Everson/Coqueiro/Antártica/dados/geneses_mensais/Ciclogeneses/'
pasta2 = '/work/archive/Everson/Coqueiro/Antártica/dados/geneses_mensais/Ciclolises/'

# Listando todos os arquivos em cada pasta
lista_arquivos1 = os.listdir(pasta1)
lista_arquivos2 = os.listdir(pasta2)

# Iterando sobre a lista de pares de arquivos
for arquivo1, arquivo2 in zip(lista_arquivos1, lista_arquivos2):
        # print(arquivo1,arquivo2)
    
        path_0 = np.array(pd.read_csv(pasta1 + arquivo1, sep=',', usecols=[0, 1, 2, 3, 4]))
        path_1 = np.array(pd.read_csv(pasta2 + arquivo2, sep=',', usecols=[0, 1, 2, 3, 4]))
        
        # # arredonda para resolução de 0.25 graus
        # path_0['Latitude'] = round(path_0[3] / 0.25) * 0.25
        # path_0['Longitude'] = round(path_0.Longitude / 0.25) * 0.25
        # path_1['Latitude'] = round(path_1.Latitude / 0.25) * 0.25
        # path_1['Longitude'] = round(path_1.Longitude / 0.25) * 0.25
        
        df_ciclogeneses = pd.DataFrame(path_0, columns = ['N_do_Ciclone_Ciclogêneses', 'Data', 'Longitude', 'Latitude', 'Vorticidade'])
    
        df_ciclolises = pd.DataFrame(path_1, columns = ['N_do_Ciclone_Ciclolises', 'Data', 'Longitude', 'Latitude', 'Vorticidade'])
    
        
        
        nome = arquivo1.split('.')[0]
        
        area_1_df_ciclogeneses = df_ciclogeneses[(df_ciclogeneses.Longitude>=-90) & (df_ciclogeneses.Longitude<=-70) & (df_ciclogeneses.Latitude<=-60) & (df_ciclogeneses.Latitude>=-66)]
        print(f'Para {(nome)} a quantidade de ciclogêneses na Area 1 : ' + str(area_1_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count()))
        
        total_1_df_ciclogeneses = str(area_1_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count())
        
        area_2_df_ciclogeneses = df_ciclogeneses[(df_ciclogeneses.Longitude>=-90) & (df_ciclogeneses.Longitude<=-70) & (df_ciclogeneses.Latitude<=-66) & (df_ciclogeneses.Latitude>=-73)]
        print(f'Para {(nome)} a quantidade de ciclogêneses na Area 2: ' + str(area_2_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count()))
        
        total_2_df_ciclogeneses = str(area_2_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count())
        
        area_3_df_ciclogeneses = df_ciclogeneses[(df_ciclogeneses.Longitude>=-50) & (df_ciclogeneses.Longitude<=-20) & (df_ciclogeneses.Latitude<=-60) & (df_ciclogeneses.Latitude>=-68)]
        print(f'Para {(nome)} a quantidade de ciclogêneses na Area 3: ' + str(area_3_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count()))
        
        total_3_df_ciclogeneses = str(area_3_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count())
        
        area_4_df_ciclogeneses = df_ciclogeneses[(df_ciclogeneses.Longitude>=-50) & (df_ciclogeneses.Longitude<=-20) & (df_ciclogeneses.Latitude<=-68) & (df_ciclogeneses.Latitude>=-78)]
        print(f'Para {(nome)} a quantidade de ciclogêneses na Area 4: ' + str( area_4_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count()))
        
        total_4_df_ciclogeneses = str(area_4_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count())
        
        zone_1_df_ciclogeneses = df_ciclogeneses[(df_ciclogeneses.Longitude>=-67.5) & (df_ciclogeneses.Longitude<=-60) & (df_ciclogeneses.Latitude<=-54) & (df_ciclogeneses.Latitude>=-60)]
        print(f'Para {(nome)} a quantidade de ciclogêneses na Zona 1: ' + str(zone_1_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count()))
        
        total_8_df_ciclogeneses = str(area_1_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count())
        
        # zone_2_df_ciclogeneses = df_ciclogeneses[(df_ciclogeneses.Longitude>=-55) & (df_ciclogeneses.Longitude<=-54) & (df_ciclogeneses.Latitude<=-60) & (df_ciclogeneses.Latitude>=-67.50)]
        # print(f'Para {(nome)} a quantidade de ciclogêneses na Zona 2: ' + str(zone_2_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count()))
        
        zone_3_df_ciclogeneses = df_ciclogeneses[(df_ciclogeneses.Longitude>=-70) & (df_ciclogeneses.Longitude<=-60) & (df_ciclogeneses.Latitude<=-60) & (df_ciclogeneses.Latitude>=-66)]
        print(f'Para {(nome)} a quantidade de ciclogêneses na Zona 3: ' + str(zone_3_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count()))
        
        total_5_df_ciclogeneses = str(zone_3_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count())
        
        zone_4_df_ciclogeneses = df_ciclogeneses[(df_ciclogeneses.Longitude>=-70) & (df_ciclogeneses.Longitude<=-50) & (df_ciclogeneses.Latitude<=-90) & (df_ciclogeneses.Latitude>=-66)]
        print(f'Para {(nome)} a quantidade de ciclogêneses na Zona 4: ' + str(zone_4_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count()))
        
        total_6_df_ciclogeneses = str(zone_4_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count())
        
        zone_5_df_ciclogeneses = df_ciclogeneses[(df_ciclogeneses.Longitude>=-60) & (df_ciclogeneses.Longitude<=-50) & (df_ciclogeneses.Latitude<=-60) & (df_ciclogeneses.Latitude>=-66)]
        print(f'Para {(nome)} a quantidade de ciclogêneses na Zona 5: ' + str(zone_5_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count()))
        
        total_7_df_ciclogeneses = str(zone_5_df_ciclogeneses['N_do_Ciclone_Ciclogêneses'].count())    
        
        
        area_1_df_ciclolises = df_ciclolises[(df_ciclolises.Longitude>=-90) & (df_ciclolises.Longitude<=-70) & (df_ciclolises.Latitude<=-60) & (df_ciclolises.Latitude>=-66)]
        print(f'Para {(nome)} a quantidade de ciclolises na Area 1 : ' + str(area_1_df_ciclolises['N_do_Ciclone_Ciclolises'].count()))
        
        total_1_df_ciclolises = str(area_1_df_ciclolises['N_do_Ciclone_Ciclolises'].count())
        
        area_2_df_ciclolises = df_ciclolises[(df_ciclolises.Longitude>=-90) & (df_ciclolises.Longitude<=-70) & (df_ciclolises.Latitude<=-66) & (df_ciclolises.Latitude>=-73)]
        print(f'Para {(nome)} a quantidade de ciclolises na Area 2: ' + str(area_2_df_ciclolises['N_do_Ciclone_Ciclolises'].count()))
        
        total_2_df_ciclolises = str(area_2_df_ciclolises['N_do_Ciclone_Ciclolises'].count())
        
        area_3_df_ciclolises = df_ciclolises[(df_ciclolises.Longitude>=-50) & (df_ciclolises.Longitude<=-20) & (df_ciclolises.Latitude<=-60) & (df_ciclolises.Latitude>=-68)]
        print(f'Para {(nome)} a quantidade de ciclolises na Area 3: ' + str(area_3_df_ciclolises['N_do_Ciclone_Ciclolises'].count()))
        
        total_3_df_ciclolises = str(area_3_df_ciclolises['N_do_Ciclone_Ciclolises'].count())
        
        area_4_df_ciclolises = df_ciclolises[(df_ciclolises.Longitude>=-50) & (df_ciclolises.Longitude<=-20) & (df_ciclolises.Latitude<=-68) & (df_ciclolises.Latitude>=-78)]
        print(f'Para {(nome)} a quantidade de ciclolises na Area 4: ' + str(area_4_df_ciclolises['N_do_Ciclone_Ciclolises'].count()))
        
        total_4_df_ciclolises = str(area_4_df_ciclolises['N_do_Ciclone_Ciclolises'].count())
        
        zone_1_df_ciclolises = df_ciclolises[(df_ciclolises.Longitude>=-67.50) & (df_ciclolises.Longitude<=-60) & (df_ciclolises.Latitude<=-54) & (df_ciclolises.Latitude>=-60)]
        print(f'Para {(nome)} a quantidade de ciclolises na Zona 1: ' + str(zone_1_df_ciclolises['N_do_Ciclone_Ciclolises'].count()))
        
        total_8_df_ciclolises = str(zone_1_df_ciclolises['N_do_Ciclone_Ciclolises'].count())
        
        # zone_2 = df_ciclolises[(df_ciclolises.Longitude>=-67.5) & (df_ciclolises.Longitude<=-60) & (df_ciclolises.Latitude<=-55) & (df_ciclolises.Latitude>=-60)]
        # print(f'Para {(nome)} a quantidade de ciclolises na Zona 2: ' + str(zone_2['N_do_Ciclone_Ciclolises'].count()))
        
        zone_3_df_ciclolises = df_ciclolises[(df_ciclolises.Longitude>=-70) & (df_ciclolises.Longitude<=-60) & (df_ciclolises.Latitude<=-60) & (df_ciclolises.Latitude>=-66)]
        print(f'Para {(nome)} a quantidade de ciclolises na Zona 3: ' + str(zone_3_df_ciclolises['N_do_Ciclone_Ciclolises'].count()))
        
        total_5_df_ciclolises = str(zone_3_df_ciclolises['N_do_Ciclone_Ciclolises'].count())
        
        zone_4_df_ciclolises = df_ciclolises[(df_ciclolises.Longitude>=-70) & (df_ciclolises.Longitude<=-50) & (df_ciclolises.Latitude<=-90) & (df_ciclolises.Latitude>=-66)]
        print(f'Para {(nome)} a quantidade de ciclolises na Zona 4: ' + str(zone_4_df_ciclolises['N_do_Ciclone_Ciclolises'].count()))
        
        total_6_df_ciclolises = str(zone_4_df_ciclolises['N_do_Ciclone_Ciclolises'].count())
        
        zone_5_df_ciclolises = df_ciclolises[(df_ciclolises.Longitude>=-60) & (df_ciclolises.Longitude<=-50) & (df_ciclolises.Latitude<=-60) & (df_ciclolises.Latitude>=-66)]
        print(f'Para {(nome)} a quantidade de ciclolises na Zona 5: ' + str(zone_5_df_ciclolises['N_do_Ciclone_Ciclolises'].count()))
        
        total_7_df_ciclolises = str(zone_5_df_ciclolises['N_do_Ciclone_Ciclolises'].count())
        
        
        # cria uma figura e eixos

        #lat e lon extend da imagem

        latmin, latmax = -90, -50
        lonmin, lonmax = -95, -15

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
        ax.coastlines(resolution='50m', color='black', linewidth=1)

        # define as coordenadas de latitude e longitude do retângulo
        lat_min_area_1, lon_min_area_1 = -66.0, -90.0
        lat_max_area_1, lon_max_area_1 = -60.0, -70.0

        lat_min_area_2, lon_min_area_2 = -73.0, -90.0
        lat_max_area_2, lon_max_area_2 = -66.2, -70.0

        lat_min_area_3, lon_min_area_3 = -68.0, -50.0
        lat_max_area_3, lon_max_area_3 = -60.0, -20.0

        lat_min_area_4, lon_min_area_4 = -78.0, -50.0
        lat_max_area_4, lon_max_area_4 = -68.1, -20.0

        lat_min_zone_1, lon_min_zone_1 = -60.0, -67.5
        lat_max_zone_1, lon_max_zone_1 = -54.0, -60.0

        # lat_min_zone_2, lon_min_zone_2 = -55.0, -67.5
        # lat_max_zone_2, lon_max_zone_2 = -54.0, -60.0

        lat_min_zone_3, lon_min_zone_3 = -66.0, -69.9
        lat_max_zone_3, lon_max_zone_3 = -60.0, -60.0

        lat_min_zone_4, lon_min_zone_4 = -80.0, -69.9
        lat_max_zone_4, lon_max_zone_4 = -66.2, -50.1

        lat_min_zone_5, lon_min_zone_5 = -66.0, -59.9
        lat_max_zone_5, lon_max_zone_5 = -60.0, -50.1

        # converte as coordenadas de latitude e longitude para coordenadas de pixel
        x_min_area_1, y_min_area_1 = ax.projection.transform_point(lon_min_area_1, lat_min_area_1, ccrs.PlateCarree(), ax.transData)
        x_max_area_1, y_max_area_1 = ax.projection.transform_point(lon_max_area_1, lat_max_area_1, ccrs.PlateCarree(), ax.transData)

        x_min_area_2, y_min_area_2 = ax.projection.transform_point(lon_min_area_2, lat_min_area_2, ccrs.PlateCarree(), ax.transData)
        x_max_area_2, y_max_area_2 = ax.projection.transform_point(lon_max_area_2, lat_max_area_2, ccrs.PlateCarree(), ax.transData)

        x_min_area_3, y_min_area_3 = ax.projection.transform_point(lon_min_area_3, lat_min_area_3, ccrs.PlateCarree(), ax.transData)
        x_max_area_3, y_max_area_3 = ax.projection.transform_point(lon_max_area_3, lat_max_area_3, ccrs.PlateCarree(), ax.transData)

        x_min_area_4, y_min_area_4 = ax.projection.transform_point(lon_min_area_4, lat_min_area_4, ccrs.PlateCarree(), ax.transData)
        x_max_area_4, y_max_area_4 = ax.projection.transform_point(lon_max_area_4, lat_max_area_4, ccrs.PlateCarree(), ax.transData)

        x_min_zone_1, y_min_zone_1 = ax.projection.transform_point(lon_min_zone_1, lat_min_zone_1, ccrs.PlateCarree(), ax.transData)
        x_max_zone_1, y_max_zone_1 = ax.projection.transform_point(lon_max_zone_1, lat_max_zone_1, ccrs.PlateCarree(), ax.transData)

        # x_min_zone_2, y_min_zone_2 = ax.projection.transform_point(lon_min_zone_2, lat_min_zone_2, ccrs.PlateCarree(), ax.transData)
        # x_max_zone_2, y_max_zone_2 = ax.projection.transform_point(lon_max_zone_2, lat_max_zone_2, ccrs.PlateCarree(), ax.transData)

        x_min_zone_3, y_min_zone_3 = ax.projection.transform_point(lon_min_zone_3, lat_min_zone_3, ccrs.PlateCarree(), ax.transData)
        x_max_zone_3, y_max_zone_3 = ax.projection.transform_point(lon_max_zone_3, lat_max_zone_3, ccrs.PlateCarree(), ax.transData)

        x_min_zone_4, y_min_zone_4 = ax.projection.transform_point(lon_min_zone_4, lat_min_zone_4, ccrs.PlateCarree(), ax.transData)
        x_max_zone_4, y_max_zone_4 = ax.projection.transform_point(lon_max_zone_4, lat_max_zone_4, ccrs.PlateCarree(), ax.transData)

        x_min_zone_5, y_min_zone_5 = ax.projection.transform_point(lon_min_zone_5, lat_min_zone_5, ccrs.PlateCarree(), ax.transData)
        x_max_zone_5, y_max_zone_5 = ax.projection.transform_point(lon_max_zone_5, lat_max_zone_5, ccrs.PlateCarree(), ax.transData)

        # calcula a largura e altura do retângulo em pixels
        w_area_1, h_area_1 = abs(x_max_area_1 - x_min_area_1), abs(y_max_area_1 - y_min_area_1)
        w_area_2, h_area_2 = abs(x_max_area_2 - x_min_area_2), abs(y_max_area_2 - y_min_area_2)
        w_area_3, h_area_3 = abs(x_max_area_3 - x_min_area_3), abs(y_max_area_3 - y_min_area_3)
        w_area_4, h_area_4 = abs(x_max_area_4 - x_min_area_4), abs(y_max_area_4 - y_min_area_4)
        w_zone_1, h_zone_1 = abs(x_max_zone_1 - x_min_zone_1), abs(y_max_zone_1 - y_min_zone_1)
        # w_zone_2, h_zone_2 = abs(x_max_zone_2 - x_min_zone_2), abs(y_max_zone_2 - y_min_zone_2)
        w_zone_3, h_zone_3 = abs(x_max_zone_3 - x_min_zone_3), abs(y_max_zone_3 - y_min_zone_3)
        w_zone_4, h_zone_4 = abs(x_max_zone_4 - x_min_zone_4), abs(y_max_zone_4 - y_min_zone_4)
        w_zone_5, h_zone_5 = abs(x_max_zone_5 - x_min_zone_5), abs(y_max_zone_5 - y_min_zone_5)

        # cria um retângulo com coordenadas (x, y), largura w e altura h
        area_1 = Rectangle((x_min_area_1, y_min_area_1), w_area_1, h_area_1, linewidth=1.5, edgecolor='black', facecolor='none')
        area_2 = Rectangle((x_min_area_2, y_min_area_2), w_area_2, h_area_2, linewidth=1.5, edgecolor='black', facecolor='none')
        area_3 = Rectangle((x_min_area_3, y_min_area_3), w_area_3, h_area_3, linewidth=1.5, edgecolor='black', facecolor='none')
        area_4 = Rectangle((x_min_area_4, y_min_area_4), w_area_4, h_area_4, linewidth=1.5, edgecolor='black', facecolor='none')
        zone_1 = Rectangle((x_min_zone_1, y_min_zone_1), w_zone_1, h_zone_1, linewidth=1.5, edgecolor='black', facecolor='none')
        # zone_2 = Rectangle((x_min_zone_2, y_min_zone_2), w_zone_2, h_zone_2, linewidth=1.5, edgecolor='black', facecolor='none')
        zone_3 = Rectangle((x_min_zone_3, y_min_zone_3), w_zone_3, h_zone_3, linewidth=1.5, edgecolor='black', facecolor='none')
        zone_4 = Rectangle((x_min_zone_4, y_min_zone_4), w_zone_4, h_zone_4, linewidth=1.5, edgecolor='black', facecolor='none')
        zone_5 = Rectangle((x_min_zone_5, y_min_zone_5), w_zone_5, h_zone_5, linewidth=1.5, edgecolor='black', facecolor='none')

        # adiciona o retângulo na figura
        ax.add_patch(area_1)
        # adiciona o texto dentro do retângulo
        area_1_x = (x_min_area_1 + x_max_area_1) / 2
        area_1_y = (y_min_area_1 + y_max_area_1) / 2
        ax.text(area_1_x, area_1_y, "Área 1", ha='center', va='bottom', fontsize=20)
        ax.text(area_1_x, area_1_y, "Ciclogêneses:" + total_1_df_ciclogeneses + "\n""Ciclolises:" + total_1_df_ciclolises, ha='center', va='top', fontsize=14)

        ax.add_patch(area_2)
        # adiciona o texto dentro do retângulo
        area_2_x = (x_min_area_2 + x_max_area_2) / 2
        area_2_y = (y_min_area_2 + y_max_area_2) / 2
        ax.text(area_2_x, area_2_y, "Área 2", ha='center', va='bottom', fontsize=20)
        ax.text(area_2_x, area_2_y, f"Ciclogêneses:" + total_2_df_ciclogeneses + "\n""Ciclolises:" + total_2_df_ciclolises, ha='center', va='top', fontsize=14)

        ax.add_patch(area_3)
        # adiciona o texto dentro do retângulo
        area_3_x = (x_min_area_3 + x_max_area_3) / 2
        area_3_y = (y_min_area_3 + y_max_area_3) / 2
        ax.text(area_3_x, area_3_y, "Área 3", ha='center', va='bottom', fontsize=20)
        ax.text(area_3_x, area_3_y, "Ciclogêneses:" + total_3_df_ciclogeneses + "\n""Ciclolises:"  + total_3_df_ciclolises, ha='center', va='top', fontsize=14)

        ax.add_patch(area_4)
        # adiciona o texto dentro do retângulo
        area_4_x = (x_min_area_4 + x_max_area_4) / 2
        area_4_y = (y_min_area_4 + y_max_area_4) / 2
        ax.text(area_4_x, area_4_y, "Área 4", ha='center', va='bottom', fontsize=20)
        ax.text(area_4_x, area_4_y, "Ciclogêneses:" + total_4_df_ciclogeneses + "\n""Ciclolises:" + total_4_df_ciclolises, ha='center', va='top', fontsize=14)

        ax.add_patch(zone_1)
        # adiciona o texto dentro do retângulo
        zone_1_x = (x_min_zone_1 + x_max_zone_1) / 2
        zone_1_y = (y_min_zone_1 + y_max_zone_1) / 2
        ax.text(zone_1_x, zone_1_y, "Drake", ha='center', va='bottom', fontsize=17)
        ax.text(zone_1_x, zone_1_y, "Gêneses:" + total_8_df_ciclogeneses + "\n""Lises:"  + total_8_df_ciclolises, ha='center', va='top', fontsize=12)

        # ax.add_patch(zone_2)
        # # adiciona o texto dentro do retângulo
        # zone_2_x = (x_min_zone_2 + x_max_zone_2) / 2
        # zone_2_y = (y_min_zone_2 + y_max_zone_2) / 2
        # ax.text(zone_2_x, zone_2_y, "FDM", ha='center', va='bottom', fontsize=7)
        # ax.text(zone_2_x, zone_2_y, "Gêneses:" + total_5_df_ciclogeneses + "\n""Lises:"  + total_5_df_ciclolises, ha='center', va='top', fontsize=5)


        ax.add_patch(zone_3)
        # adiciona o texto dentro do retângulo
        zone_3_x = (x_min_zone_3 + x_max_zone_3) / 2
        zone_3_y = (y_min_zone_3 + y_max_zone_3) / 2
        ax.text(zone_3_x, zone_3_y, "Zona 3", ha='center', va='bottom', fontsize=17)
        ax.text(zone_3_x, zone_3_y, "Gêneses:" + total_5_df_ciclogeneses + "\n""Lises:"  + total_5_df_ciclolises, ha='center', va='top', fontsize=12)

        ax.add_patch(zone_4)
        # adiciona o texto dentro do retângulo
        zone_4_x = (x_min_zone_4 + x_max_zone_4) / 2
        zone_4_y = (y_min_zone_4 + y_max_zone_4) / 2
        ax.text(zone_4_x, zone_4_y, "Zona 4", ha='center', va='bottom', fontsize=20)
        ax.text(zone_4_x, zone_4_y, "Ciclogêneses:" + total_6_df_ciclogeneses + "\n""Ciclolises:"  + total_6_df_ciclolises, ha='center', va='top', fontsize=15)

        ax.add_patch(zone_5)
        zone_5_x = (x_min_zone_5 + x_max_zone_5) / 2
        zone_5_y = (y_min_zone_5 + y_max_zone_5) / 2
        ax.text(zone_5_x, zone_5_y, "Zona 5", ha='center', va='bottom', fontsize=17)
        ax.text(zone_5_x, zone_5_y, "Gêneses:" + total_7_df_ciclogeneses + "\n""Lises:"  + total_7_df_ciclolises, ha='center', va='top', fontsize=12)
        
        #plt.title('Valid time: {}'.format(vtime), fontsize=35, loc='right')
        #analise
        plt.title('Gênese de ciclones \n'f'{(nome)}', fontsize=35, loc='center')
        # # exibe a figura
        # plt.show()
        
        plt.savefig(f'/work/archive/Everson/Coqueiro/Antártica/imagens/Geneses/{nome}.png', bbox_inches='tight')
        
        
    