#!/usr/bin/env python 

import sys ; import numpy as np ; import scipy as sp ; import datetime
from grads import * ; 
import matplotlib.pyplot as plt; from pylab import * ; from matplotlib import *

from mpl_toolkits.basemap import Basemap, cm, shiftgrid, maskoceans, interp, shapefile
from scipy.io import netcdf 
from netCDF4 import Dataset


nc_file='wrfpost_UGRD10m.nc'
fh     =Dataset(nc_file,'r')
lons   =fh.variables['lon'][:]
lats   =fh.variables['lat'][:]
time   =fh.variables['time'][:]
u10   =fh.variables['ugrd10m'][24:48,:,:]
fh.close()
nc_file ='wrfpost_VGRD10m.nc'
fh      =Dataset(nc_file,'r')
v10     =fh.variables['vgrd10m'][24:48,:,:]
fh.close()

nc_file ='wrfpost_MSLETmsl.nc'
fh      =Dataset(nc_file,'r')
slp1     =fh.variables['msletmsl'][6:12,:,:]
fh.close

uin=u10.max(axis=0)
vin=v10.max(axis=0)
p=0.01*(slp1.mean(axis=0))

lon, lat = np.meshgrid(lons,lats) 

#m = Basemap(projection='ortho',lat_0=lats.min(),lon_0=lons.min(), resolution='l')
m=Basemap(projection='merc',lat_ts=10, llcrnrlon=lons.min(),urcrnrlon=lons.max()-3, llcrnrlat=lats.min()+3,urcrnrlat=lats.max(), resolution='c')
fig1 = plt.figure(figsize=(8,10)) ; ax = fig1.add_axes([0.1,0.1,0.8,0.8])
x, y = m(lon, lat)

cmap = plt.cm.Reds
cmap = plt.cm.Spectral
cmap = plt.cm.RdBu_r

clevs = np.arange(990,1010,1)
cs1 = m.contour(x,y,p,clevs,linewidths=1.0,colors='k',animated=True)
cs2 = m.contourf(x,y,p,clevs,cmap=cmap,animated=True)

ugrid,newlons = shiftgrid(62.,uin,lons,start=True)
vgrid,newlons = shiftgrid(62.,vin,lons,start=True)

uproj,vproj,xx,yy = m.transform_vector(ugrid,vgrid,newlons,lats,51,51,returnxy=True,masked=True)
Q = m.quiver(xx,yy,uproj,vproj,scale=250,zorder=10,color='r')
qk = plt.quiverkey(Q, 0.1, 0.1, 20, '15 m/s', labelpos='W',coordinates='data',color='r')


m.readshapefile('subdiv/India_subdiv','ind',drawbounds=True, zorder=None, linewidth=1.0, color='k', antialiased=1, ax=None, default_encoding='utf-8')
m.drawmapboundary(fill_color='aqua') ; #m.fillcontinents(color='white',lake_color='white')
m.drawparallels(np.arange(0.,40.,10.), labels=[1,0,0,0])
m.drawmeridians(np.arange(60.,100.,10.), labels=[0,0,0,1])
cbar=m.colorbar(cs2,location='bottom', pad='5%') ; cbar.set_label('hpa') ; 
#plt.clabel(cs2, fontsize=9, inline=1)
m.bluemarble() ;
plt.title('Mean Sea Level Pressure(hpa) & Wind Vectors'); savefig('slp.png', dpi=100);
plt.show()






