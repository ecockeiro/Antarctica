#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 11:57:20 2023

@author: everson
"""

from datetime import datetime
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import pandas as pd
import cartopy.feature as cfeature
import glob 
import os
import shutil

# Paths

# # Diretório dos dados .nc de componente u e v
# diretorio_pastas = '/media/bjerknes/HD_todo_pod/Everson/Coqueiro/CENPES/DADOS/Raoni/dados_input_vort'

# Diretorio do indat TRACK
indat = "/home/bjerknes/TRACK-1.5.0/indat" 

# Diretorio do outdat TRACK
outdat = "/home/bjerknes/TRACK-1.5.0/outdat"

# Definindo o caminho para a pasta onde estão os arquivos .gz
caminho_dados_tr_trs_neg = "/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Monografia/Track/master"

# Diretorio ../utils/bin
utils_bin = '/home/bjerknes/TRACK-1.5.0/utils/bin'


# função que lê as pastas em um diretório
# def get_folders(directory):
#     return [f for f in os.listdir(directory) if os.path.isdir(os.path.join(directory, f))]

# # lista de pasta para serem selecionadas
# folders = get_folders(diretorio_pastas)
# folders.sort(reverse=True)

# # Converte as componentes u e v em vorticidade relativa 
# for j in range(len(folders)):
#     # path_name = str(diretorio_pastas) + '/' + str(folders[j])
    
#     path_name = os.path.join(diretorio_pastas, folders[j])
    
#     ## Abre os dados na pasta
#     data = xr.open_mfdataset(f'{path_name}/*grib2.nc')
#     data['u'] = data['U_GRD_L100'].copy()
#     data['v'] = data['V_GRD_L100'].copy()
    
#     Data = str(folders[j])
    
#     # data = data.sel(level=850)
    
#     # Vai preenchendo o np.array...
#     vort_850= mpcalc.vorticity(data.u, data.v)
    
    
#     #sys.exit()    
#     # Cria o datarray dentro da estrutura do data 
#     data['vo'] = data['u'].copy()
#     data['vo'].data = vort_850
#     data['vo'] = xr.DataArray(vort_850, dims = ['time','lat','lon'])   
#     data['vo'].attrs['standard_name'] = 'atmosphere_relative_vorticity'
#     data['vo'].attrs['long_name'] = 'Vorticity (relative)'
#     data['vo'].attrs['units'] = 's**-1'
#     data['lat'].attrs['units'] = 'degrees_north'
#     data['lon'].attrs['units'] = 'degrees_east'
    
    
#     # Deixa somente o vort_850 dentro do data
#     data = data[['vo']]
#     data = data.chunk()

#     # Salva no arquivo netcdf
#     fname = f'/media/bjerknes/HD_todo_pod/Everson/Coqueiro/CENPES/DADOS/Raoni/vort_nc/{(Data)}.nc'
    
#     data.to_netcdf(f'/media/bjerknes/HD_todo_pod/Everson/Coqueiro/CENPES/DADOS/Raoni/vort_nc/{(Data)}.nc')
    
#     os.system("ncap2 -O -s 'time=double(time)' {} {}".format(fname,fname))
#     shutil.copy(fname,indat)
#     print ('finished saving')

# # Diretório dos dados .nc de vorticidade relativa
diretorio_pastas = '/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Monografia/Dados/vorticidade_relativa_850/mensal'

# função que lê as pastas em um diretório
def get_folders(directory):
    return [f for f in os.listdir(directory) if os.path.isdir(os.path.join(directory, f))]

# lista de pasta para serem selecionadas
folders = get_folders(diretorio_pastas)
folders.sort(reverse=True)

# Converte as componentes u e v em vorticidade relativa 
for j in range(len(folders)):
    # path_name = str(diretorio_pastas) + '/' + str(folders[j])
    
    path = os.path.join(diretorio_pastas, folders[j])
    
    # ## Abre os dados na pasta
    # data = xr.open_mfdataset(f'{path_name}/*.nc')

    # path = '/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Monografia/Dados/vorticidade_relativa_850/mensal/2010'
    for fname in os.listdir(path):
        shutil.copy(os.path.join(path, fname), indat)
        print ('finished copy')
        
    
    # Percorra todos os arquivos .nc na pasta e execute o comando bin/track.linux com o filtro especificado em specfiltT42.in
    for arquivo in os.listdir(indat):
        sorted(arquivo)
        if arquivo.endswith(".nc"):
            nome_arquivo_sem_extensao = os.path.splitext(arquivo)[0]
            os.system(f"/home/bjerknes/TRACK-1.5.0/bin/track.linux -i {arquivo} -f filt < /home/bjerknes/TRACK-1.5.0/filtro_era5.in")
            
            caminho_spectral = os.path.join(outdat, "specfil.filt_band000")
            caminho_dat = os.path.join(outdat, nome_arquivo_sem_extensao + ".dat")
            os.system('mv ' + caminho_spectral + ' ' + caminho_dat)
            
            caminho_spectral = caminho_dat
            caminho_dat = os.path.join(indat, nome_arquivo_sem_extensao + ".dat")
            os.system('mv ' + caminho_spectral + ' ' + caminho_dat)
            
    
    # Muda para pasta TRACK
    os.chdir('/home/bjerknes/TRACK-1.5.0/') 
            
    # Percorra todos os arquivos .dat na pasta indat e execute o comando ./master com as opções especificadas
    for arquivo in os.listdir(indat):
        if arquivo.endswith(".dat"):
            nome_arquivo_sem_extensao = os.path.splitext(arquivo)[0]
            os.system(f"./master -c={nome_arquivo_sem_extensao} -d=now -e=track.linux -i={arquivo} -f=vo1991 -n=1,1000,1 -o=/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Monografia/Track/master -r=RUN_AT_ -s=RUNDATIN.VOR -j=RUN_AT.in")
    
    
    
    # Percorrer cada diretório e descompactar o arquivo tr_trs_neg.gz
    for diretorio in os.listdir(caminho_dados_tr_trs_neg):
        if os.path.isdir(os.path.join(caminho_dados_tr_trs_neg, diretorio)):
            caminho_arquivo_gz = os.path.join(caminho_dados_tr_trs_neg, diretorio, "tr_trs_neg.gz")
            if os.path.exists(caminho_arquivo_gz):
                os.system(f"gunzip {caminho_arquivo_gz}")
                
                dia = '01'
                mes = diretorio[5:7]
                ano = diretorio[0:4]
                hora = '00'
                data = ano+mes+dia+hora
                
                caminho_arquivo_descompactado = os.path.join(caminho_dados_tr_trs_neg, diretorio, "tr_trs_neg")
                
                shutil.copy(caminho_arquivo_descompactado, utils_bin)
                
                # Muda para o diretrotio /utils/bin/
                os.chdir('/home/bjerknes/TRACK-1.5.0/utils/bin/')    
                
                os.system(f"./count tr_trs_neg 0 0 5 4 0 {data} 1")
    
                os.system('mv ' + os.path.join(utils_bin,'tr_trs_neg.new') +' '+ os.path.join(utils_bin,diretorio))
    
                os.system(f"tr2csv {diretorio}")
                
    
                os.system('mv ' + os.path.join(utils_bin,'alltr.csv') + ' ' + os.path.join(utils_bin, f'{diretorio}.csv'))
                
                path_ciclone = f'/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Monografia/Track/planilhas_sem_recorte/{ano}'
    
                os.system('cp ' + os.path.join(utils_bin, f'{diretorio}.csv') +' '+ os.path.join(path_ciclone, f'{diretorio}.csv'))
                os.remove(f'{diretorio}')
                
                # usar somente caso os dados nao estejam salvos como 01, 02, 03 e sim como 1, 2, 3..
                """
                if int(dia) >= 10 and int(dia) <=31:
                    os.system(f"./count tr_trs_neg 0 0 5 4 0 {data} 1")
        
                    os.system('mv ' + os.path.join(utils_bin,'tr_trs_neg.new') +' '+ os.path.join(utils_bin,diretorio))
        
                    os.system(f"tr2csv {diretorio}")
                    
        
                    os.system('mv ' + os.path.join(utils_bin,'alltr.csv') + ' ' + os.path.join(utils_bin, f'{diretorio}.csv'))
                    
                    path_ciclone = '/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Monografia/Track/planilhas_sem_recorte'
        
                    os.system('cp ' + os.path.join(utils_bin, f'{diretorio}.csv') +' '+ os.path.join(path_ciclone, f'{diretorio}.csv'))
                    os.remove(f'{diretorio}')
                else:
                    dia = diretorio[:1]
                    mes = diretorio[1:3]
                    ano = diretorio[3:7]
                    hora = '00'
                    data = ano+mes+'0'+dia+hora
                    
                    os.system(f"./count tr_trs_neg 0 0 5 4 0 {data} 3")
        
                    os.system('mv ' + os.path.join(utils_bin,'tr_trs_neg.new') +' '+ os.path.join(utils_bin,diretorio))
        
                    os.system(f"tr2csv {diretorio}")
                    new_name = "0" + f"{diretorio}.csv"
        
                    os.system('mv ' + os.path.join(utils_bin,'alltr.csv') + ' ' + os.path.join(utils_bin, new_name))
                    
                    path_ciclone = '/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Antartica/siac_2023/Dados_usados/trakcs_planilhas'
        
                    os.system('cp ' + os.path.join(utils_bin, new_name) +' '+ os.path.join(path_ciclone, new_name))
                    os.remove(f'{diretorio}')
                    """
    
            
    # I.2 - Carrega a pasta contendo todas as planilhas dos Tracks
    
    datafiles = glob.glob(f'/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Monografia/Dados/vorticidade_relativa_850/todas_planilhas_recortadas/{ano}/*.csv')
    datafiles.sort()
    
    
    # range = ao numero de meses, poderia usar o range(len(datafiles))
    for i in range(len(datafiles)):
        
        
        # II - Adiciona os rótulos nas colunas, pois quando o track gera as planilhas geram planilhas sem rótulos
        
        #seleciona a planilha
        path = np.array(pd.read_csv(datafiles[i], sep=',', usecols=[0, 1, 2, 3, 4]))
        
        #adiciona os rotulos das colunas
        df = pd.DataFrame(path, columns = ['N_do_Ciclone', 'Data', 'Longitude', 'Latitude', 'Vorticidade'])
        
        nome = datafiles[i].split('/')
        data = nome[-1].split('.')[0]
     
        # III - Converte as longitudes para -180 a 180 graus
        
        # converte e recorta latitude/longitude
        df['Longitude'] = (((df['Longitude']+180) %360) - 180)
        
        df = df[(df.Longitude>=-140) & (df.Longitude<=-10) & (df.Latitude<=-40) & (df.Latitude>=-80)].copy()
        
        
        # IV - Remove a ocorrencia de ciclones a partir de um numero pré-definido (exemplo = 4)
        
        # Conta a ocorrência de cada item
        counts = df['N_do_Ciclone'].value_counts()
        
        # Cria uma lista com os ciclones a serem removidos
        itens_para_remover = counts[counts < 6].index.tolist()
        
        # Remove os ciclones da lista do DataFrame
        df = df[~df['N_do_Ciclone'].isin(itens_para_remover)]
        
        
        # V - Salva a planilha dentro da pasta automaticamente a partir do rotulo 'Data'
        
        #salva a planilha recortada
        df.to_csv(f'/media/bjerknes/HD_todo_pod/Everson/Coqueiro/Monografia/Track/planilha_com_recorte/{ano}/{data}.csv', index=False)    
        # df.to_csv(f'/work/archive/Everson/Coqueiro/CENPES/DADOS/track/rodada_1/.csv', index=False) 
        
        
        # plota os tracks
        # Caso queira plotar as saídas das planilhas descomente essa parte!
        
        #lat e lon extend da imagem
        latmin, latmax = -84, -38
        lonmin, lonmax = -140, -10
        
        fig = plt.figure(figsize=(10,10))
        args = dict(color='gray',
                    alpha=1.0, 
                    linestyle='--', 
                    linewidth=0.5,
                    xlocs=np.arange(lonmin, lonmax, 15), 
                    ylocs=np.arange(latmin, latmax, 5), 
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
        ax.set_title(f'{data}')
        
        
        
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
    
    
    
        ciclones = list( dict.fromkeys(df.N_do_Ciclone) )
    
        for Ciclone in ciclones:
            one_ciclone = df[(df.N_do_Ciclone == Ciclone)].copy()
            plot_cyclones(one_ciclone)
            # print(one_ciclone)
            
            
            
            
        plt.show()
    
    
    # # limpa as pastas
    
    # # Remove .nc do indat        
    # for filename_nc in os.listdir(indat):
    #     if filename_nc.endswith(".nc"):
    #         file_path = os.path.join(indat,filename_nc)
    #         os.remove(file_path)
            
    # # Remove os .csv do utils_bin
    # for filename in os.listdir(utils_bin):
    #     if filename.endswith(".csv"):
    #         file_path = os.path.join(utils_bin, filename)
    #         os.remove(file_path)
            
    # # Remove .dat do indat  
    # for filename_dat in os.listdir(indat):
    #     if filename_dat.endswith('.dat'):
    #         file_path = os.path.join(indat,filename_dat)
    #         os.remove(file_path)
