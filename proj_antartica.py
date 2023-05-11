#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 10:45:35 2023

@author: ladsin
"""
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import proplot as pplt
import numpy as np
import cartopy.feature as cfeature

#delimita a lats e lons
latmin, latmax = -90, -55
lonmin, lonmax = -180, 180


#escolhe o tamanho da figura e polegadas (x,y)
plt.figure(figsize=(10,10))

args = dict(color='gray',
            alpha=1.0, 
            linestyle='--', 
            linewidth=0.5,
            xlocs=np.arange(-180, 180, 30), 
            ylocs=np.arange(-90, -45, 10), 
            draw_labels=True)

# Criação da projeção polar centrada na Antártica

ax = plt.axes(projection=ccrs.SouthPolarStereo())

gl = ax.gridlines(crs=ccrs.PlateCarree(), **args)
gl.xlabel_style = {'size': 15, 'color': 'black', 'rotation': 0}
gl.ylabel_style = {'size': 15, 'color': 'Gray', 'rotation': 0}
gl.top_labels = True
gl.right_labels = True
gl.bottom_labels = True
gl.rotate_labels = True
ax.set_extent([lonmin, lonmax, latmin, latmax], ccrs.PlateCarree())
ax.add_feature(cfeature.LAND)
ax.coastlines(resolution='50m', color='black', linewidth=1)




