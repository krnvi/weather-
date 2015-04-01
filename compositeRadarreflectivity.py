#!/usr/bin/env python 

import sys ; import numpy as np ; import scipy as sp ; import datetime
from grads import * ; 
import matplotlib.pyplot as plt; from pylab import * ; from matplotlib import *

from mpl_toolkits.basemap import Basemap, cm, shiftgrid, maskoceans, interp, shapefile
from scipy.io import netcdf 
from netCDF4 import Dataset


nc_file='wrfpost_REFCclm.nc'
fh     =Dataset(nc_file,'r')
lons   =fh.variables['lon'][:]
lats   =fh.variables['lat'][:]
time   =fh.variables['time'][:]
avc    =fh.variables['refcclm'][0:384,:,:]
fh.close()

avc_1=avc[0:95,:,:].mean(axis=0)
avc_2=avc[96:191,:,:].mean(axis=0)
avc_3=avc[192:287,:,:].mean(axis=0)
avc_4=avc[288:383,:,:].mean(axis=0)


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
#ax = fig.add_axes([0.1,0.1,0.8,0.8])

subplot(221)

m=Basemap(projection='merc',lat_ts=10, llcrnrlon=lons.min(),urcrnrlon=lons.max()-3, llcrnrlat=lats.min()+3,urcrnrlat=lats.max(), resolution='c')
#m = Basemap(projection='ortho',lat_0=lats.min(),lon_0=lons.min(), resolution='l')
lon, lat =np.meshgrid(lons, lats)
XX, YY=m(lon,lat)
cmap = cm.StepSeq
#cmap = cm.GMT_drywet
#clevs=[12,16,20,24,28,32]
#cs = m.contourf(XX,YY,avc_1,clevs,cmap=cmap)
#cs1 = m.contour(XX,YY,avc_1,clevs,linewidths=1.0,colors='k',animated=True)
cs=m.pcolormesh(XX,YY,avc_1,shading='flat', cmap=cmap,vmin=0, vmax=32);
#m.drawcoastlines(linewidth=0.25) ; m.drawcountries(linewidth=0.25); m.drawstates(linewidth=0.25) 
m.readshapefile('subdiv/India_subdiv','ind',drawbounds=True, zorder=None, linewidth=1.0, color='k', antialiased=1, ax=None, default_encoding='utf-8')
m.drawmapboundary(fill_color='white') ; m.fillcontinents(color='white',lake_color='white')
m.drawparallels(np.arange(0.,40.,10.), labels=[1,0,0,0])
m.drawmeridians(np.arange(60.,100.,10.), labels=[0,0,0,1])
#cbar=m.colorbar(cs,location='bottom', pad='10%') ; cbar.set_label('dbZ') ; 
#m.bluemarble()
plt.title('Radar Reflectivity (dbz)'); 

########################
subplot(222)
m=Basemap(projection='merc',lat_ts=10, llcrnrlon=lons.min(),urcrnrlon=lons.max()-3, llcrnrlat=lats.min()+3,urcrnrlat=lats.max(), resolution='c')
#m = Basemap(projection='ortho',lat_0=lats.min(),lon_0=lons.min(), resolution='l')
lon, lat =np.meshgrid(lons, lats)
XX, YY=m(lon,lat)
cmap = cm.StepSeq
#cmap = cm.GMT_drywet
#clevs=[12,16,20,24,28,32]
#cs = m.contourf(XX,YY,avc_1,clevs,cmap=cmap)
#cs1 = m.contour(XX,YY,avc_1,clevs,linewidths=1.0,colors='k',animated=True)
cs=m.pcolormesh(XX,YY,avc_2,shading='flat', cmap=cmap,vmin=0, vmax=32);
#m.drawcoastlines(linewidth=0.25) ; m.drawcountries(linewidth=0.25); m.drawstates(linewidth=0.25) 
m.readshapefile('subdiv/India_subdiv','ind',drawbounds=True, zorder=None, linewidth=1.0, color='k', antialiased=1, ax=None, default_encoding='utf-8')
m.drawmapboundary(fill_color='white') ; m.fillcontinents(color='white',lake_color='white')
m.drawparallels(np.arange(0.,40.,10.), labels=[1,0,0,0])
m.drawmeridians(np.arange(60.,100.,10.), labels=[0,0,0,1])
#cbar=m.colorbar(cs,location='bottom', pad='10%') ; cbar.set_label('dbZ') ; 
#m.bluemarble()
plt.title('Radar Reflectivity (dbz)'); 
#################################
subplot(223)
m=Basemap(projection='merc',lat_ts=10, llcrnrlon=lons.min(),urcrnrlon=lons.max()-3, llcrnrlat=lats.min()+3,urcrnrlat=lats.max(), resolution='c')
#m = Basemap(projection='ortho',lat_0=lats.min(),lon_0=lons.min(), resolution='l')
lon, lat =np.meshgrid(lons, lats)
XX, YY=m(lon,lat)
cmap = cm.StepSeq
#cmap = cm.GMT_drywet
#clevs=[12,16,20,24,28,32]
#cs = m.contourf(XX,YY,avc_1,clevs,cmap=cmap)
#cs1 = m.contour(XX,YY,avc_1,clevs,linewidths=1.0,colors='k',animated=True)
cs=m.pcolormesh(XX,YY,avc_3,shading='flat', cmap=cmap,vmin=0, vmax=32);
#m.drawcoastlines(linewidth=0.25) ; m.drawcountries(linewidth=0.25); m.drawstates(linewidth=0.25) 
m.readshapefile('subdiv/India_subdiv','ind',drawbounds=True, zorder=None, linewidth=1.0, color='k', antialiased=1, ax=None, default_encoding='utf-8')
m.drawmapboundary(fill_color='white') ; m.fillcontinents(color='white',lake_color='white')
m.drawparallels(np.arange(0.,40.,10.), labels=[1,0,0,0])
m.drawmeridians(np.arange(60.,100.,10.), labels=[0,0,0,1])
cbar=m.colorbar(cs,location='bottom', pad='10%') ; cbar.set_label('dbZ') ; 
#m.bluemarble()
plt.title('Radar Reflectivity (dbz)'); 
#########################################
subplot(224)
m=Basemap(projection='merc',lat_ts=10, llcrnrlon=lons.min(),urcrnrlon=lons.max()-3, llcrnrlat=lats.min()+3,urcrnrlat=lats.max(), resolution='c')
#m = Basemap(projection='ortho',lat_0=lats.min(),lon_0=lons.min(), resolution='l')
lon, lat =np.meshgrid(lons, lats)
XX, YY=m(lon,lat)
cmap = cm.StepSeq
#cmap = cm.GMT_drywet
#clevs=[12,16,20,24,28,32]
#cs = m.contourf(XX,YY,avc_1,clevs,cmap=cmap)
#cs1 = m.contour(XX,YY,avc_1,clevs,linewidths=1.0,colors='k',animated=True)
cs=m.pcolormesh(XX,YY,avc_4,shading='flat', cmap=cmap,vmin=0, vmax=32);
#m.drawcoastlines(linewidth=0.25) ; m.drawcountries(linewidth=0.25); m.drawstates(linewidth=0.25) 
m.readshapefile('subdiv/India_subdiv','ind',drawbounds=True, zorder=None, linewidth=1.0, color='k', antialiased=1, ax=None, default_encoding='utf-8')
m.drawmapboundary(fill_color='white') ; m.fillcontinents(color='white',lake_color='white')
m.drawparallels(np.arange(0.,40.,10.), labels=[1,0,0,0])
m.drawmeridians(np.arange(60.,100.,10.), labels=[0,0,0,1])
cbar=m.colorbar(cs,location='bottom', pad='10%') ; cbar.set_label('dbZ') ; 
#m.bluemarble()
plt.title('Radar Reflectivity (dbz)'); 

savefig('reflectivity.png', dpi=100); plt.show()