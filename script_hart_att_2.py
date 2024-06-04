#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 10:51:12 2024

@author: bjerknes
"""

import pandas as pd
import numpy as np
import xarray as xr
import os
import glob


# # função que lê as pastas em um diretório
def get_folders(directory):
        return [f for f in os.listdir(directory) if os.path.isdir(os.path.join(directory, f))]
    # diretório

directories = '/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Monografia/Dados/Anos'
# lista de pasta para serem selecionadas
folders = get_folders(directories)
folders.sort()

def vento_termico_hart(geo, lat, lon):
    ''' o vento termico em baixos e altos niveis, segundo hart, eh caculado a partir
    da perturbacao de altura geopotencial nos tres niveis da atmosfera: 900hPa, 600hPa e 300hPa
    Essa pertubacao em uma camada pode ser obtida pela diferenca entre a altura geopotencial maxima e minima
    dentro de um raio de 500km a partir do centro do ciclone
    '''
    
    dlon = 5 # primeira aproximacao, 5 graus é cerca de 555 km.
    dlat = 5
    sub = geo.sel(dict(longitude = slice(lon - dlon, lon + dlon), latitude = slice(lat + dlat, lat - dlat))) # Seleciona uma 'quadrado' em torno do ponto
    # Nao me pergunte o porque, mas parece que nesse netCDF, o array de latitude eh decrescente enquanto que o de longitude eh crescente
    
    lat_arr = sub.latitude.values # tamanho 41 se dlat = 5
    lon_arr = sub.longitude.values # tamanho 41  se dlon = 5
    dx = np.abs(lon_arr[1] - lon_arr[0]) # Resolucao em x, graus
    dy = np.abs(lat_arr[1] - lat_arr[0]) # Resolucao em y, graus
    
    xx, yy = np.meshgrid(lon_arr, lat_arr)
    coords = np.column_stack((np.ravel(xx), np.ravel(yy))) # todos os pontos a serem testados. Dimensao (41x41 , 2) = (1681, 2). Se dlat = dlon = 5
    x = 0 # primeira coluna em coords = longitude
    y = 1 # segunda coluna em coords = latitude
    
    # Calculo de distancia para cada ponto, provavelmente a parte mais custosa do script, em termos de eficienica. Imagina testar 1681 pontos em cada passo de tempo?!
    # Se a distancia do centro do ciclone ate cada ponto for menor que 500 km, a gente considera ele na busca do geopotencial maximo e minimo.
    # Estou usando a formula de Haversine para calcular distancia entre dois pontos em uma esfera.
    dlat = np.deg2rad(coords[:, y] - lat) # dlat em radiano, nao faz mal sobrescrever aqueles do inicio
    dlon = np.deg2rad(coords[:, x] - lon) # dlon em radiano
    R = 6371 # Raio da Terra em km
    a = np.power(np.sin(dlat/2), 2) + np.cos(np.deg2rad(lat)) * np.cos(np.deg2rad(coords[:, y])) * np.power(np.sin(dlon/2), 2)
    c = 2 *  np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    d = R * c
    bool_arr = d <= 500 # array de booleans, True se d para cada ponto for menor que 500 km.
    
    # Agora separando o geopotencial para cada uma das tres camadas
    g = 9.80665 # aceleracao da gravidade
    altura_300 = np.ravel(sub.sel(dict(level = 300)).values)/g
    altura_600 = np.ravel(sub.sel(dict(level = 600)).values)/g
    altura_900 = np.ravel(sub.sel(dict(level = 900)).values)/g
    
    # Perturbacao de altura geopotencial
    dz_300 = np.max(altura_300[bool_arr]) - np.min(altura_300[bool_arr])
    dz_600 = np.max(altura_600[bool_arr]) - np.min(altura_600[bool_arr])
    dz_900 = np.max(altura_900[bool_arr]) - np.min(altura_900[bool_arr])
    
    thermal_upper = (dz_300 - dz_600)/(np.log(300) - np.log(600))
    thermal_lower = (dz_600 - dz_900)/(np.log(600) - np.log(900))
    
    return thermal_upper, thermal_lower
    

for ano in folders[8:9]:
    datafiles = glob.glob(f'{directories}/{ano}/*.csv')
    datafiles.sort()
        
        
    # range = ao numero de meses, poderia usar o range(len(datafiles))
    for plan in datafiles[:]:
        # II - Adiciona os rótulos nas colunas, pois quando o track gera as planilhas geram planilhas sem rótulos
        #seleciona a planilha
        track = np.array(pd.read_csv(plan, sep=',', usecols=[0, 1, 2, 3, 4]))
    
        #adiciona os rotulos das colunas
        track = pd.DataFrame(track, columns=['N_do_Ciclone', 'Data', 'Longitude', 'Latitude', 'Vorticidade'])
        
        col_data = track['Data']
        track = track.replace('.', ',', regex=True)
        track['Data'] = col_data
        
        time_list = pd.to_datetime(track['Data'].str.strip(), format='"%Y-%m-%d %H:%M:%S"')
    
        nome = plan.split('/')
        data = nome[-1].split('.')[0]
        ano= nome[9]
    
        # III - Converte as longitudes para -180 a 180 graus
        track['Longitude'] = (((track['Longitude']+180) % 360) - 180)
    
        # Identificar ciclones que começam fora do domínio
        ciclones_fora_do_dominio = track[track.groupby('N_do_Ciclone')['Longitude'].transform('first') < -180]
    
        # # Remover registros desses ciclones do DataFrame original
        track = track[~track['N_do_Ciclone'].isin(ciclones_fora_do_dominio['N_do_Ciclone'])].copy()
    
        # converte e recorta latitude/longitude
        track = track[(track.Longitude >= -120) & (track.Longitude <= -10) & (track.Latitude <= -40) & (track.Latitude >= -78)].copy()
    
        # IV - Remove a ocorrencia de ciclones a partir de um numero pré-definido (exemplo = 4)
        # Conta a ocorrência de cada item
        counts = track['N_do_Ciclone'].value_counts()
        # Cria uma lista com os ciclones a serem removidos
        itens_para_remover = counts[counts < 6].index.tolist()
    
        # Remove os ciclones da lista do DataFrame
        track = track[~track['N_do_Ciclone'].isin(itens_para_remover)]
        
        # arredonda para resolução de 0.25 graus
        latitude = round(track.Latitude / 0.25) * 0.25
        longitude = round(track.Longitude / 0.25) * 0.25
        
        # Dataset netcdf
        netcdf = xr.open_dataset(f'/home/bjerknes/Downloads/geopotencial/{ano}.nc')
        geopotencial = netcdf['z']
        geopotencial = geopotencial.assign_coords(dict(longitude = (((geopotencial.longitude.values + 180) % 360) - 180)))
        track['VTU'] = 0.0  # Adiciona coluna 'VTU' com valor inicial zero
        track['VTL'] = 0.0  # Adiciona coluna 'VTL' com valor inicial zero
        track['fase'] = ''  # Adiciona coluna 'fase' com valor inicial vazio

        for t in track.index[:]: # Intervalo fechado  [ini, fim] -> inclui inicio e fim
                # current_geo = geopotencial.sel(dict(time = time_list[t])) # fixa um tempo
            current_geo = geopotencial.sel(time=time_list[t])
            upper, lower = vento_termico_hart(current_geo, latitude[t], longitude[t]) # calcula o Vt
            track.loc[t,'VTU'] = round(upper,1)
            track.loc[t,'VTL'] = round(lower,1)
            # for i in range(len(track['VTU'])):
            if -180< track['VTU'][t] < 0 and -300 < track['VTL'][t] < 0:
                track.loc[t,'fase'] = 'Extratropical em transicao Tropical'
                
            elif track['VTU'][t] < -180 and -200 < track['VTL'][t] < 0:
                track.loc[t,'fase'] = 'Ciclone Extratropical Ocluso'
                
            elif track['VTU'][t] < -180 and track['VTL'][t] < -200:
                track.loc[t,'fase'] = 'Ciclone Extratropical'
                
            elif -180 < track['VTU'][t] < -10 and track['VTL'][t] > -50:
                track.loc[t,'fase'] = 'Ciclone Subtropical'
                
            elif track['VTU'][t] < -180 and track['VTL'][t] > -50:
                track.loc[t,'fase'] = 'Aprisionamento quente/Hibridos'
                         	
            elif 0 < track['VTU'][t] < 65 and  0 < track['VTL'][t] < 60:
                track.loc[t,'fase'] = 'Tempestade Tropical'
                
            else:
                track.loc[t,'fase'] = 'Tropical'
            track.to_csv(f'/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Monografia/hart/Planilhas_de_fases/{data}.csv', index=False)
            #print("{} || -Vt_u = {:.2f} , -Vt_l = {:.2f}".format(time_list[t], upper, lower))
        #print("Fim!")
        print(track['fase'].value_counts())
        # return None
     
