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
u10   =fh.variables['ugrd10m'][6:12,:,:]
fh.close()
nc_file ='wrfpost_VGRD10m.nc'
fh      =Dataset(nc_file,'r')
v10     =fh.variables['vgrd10m'][6:12,:,:]
fh.close()

nc_file ='wrfpost_MSLETmsl.nc'
fh      =Dataset(nc_file,'r')
slp1     =fh.variables['msletmsl'][6:12,:,:]
fh.close

u=u10.max(axis=0)
v=v10.max(axis=0)
p=0.01*(slp1.mean(axis=0))

slp = np.zeros((p.shape[0],p.shape[1]+1),np.float)
slp[:,0:-1] = p[::-1]; slp[:,-1] = p[::-1,0]
#lons.append(360.) ; lons = np.array(lons)
########################################################################

fig=plt.figure(figsize=(8,8)); ax = fig.add_axes([0.1,0.1,0.8,0.8])
m=Basemap(projection='merc',lat_ts=10, llcrnrlon=lons.min(),urcrnrlon=lons.max()-3, llcrnrlat=lats.min()+3,urcrnrlat=lats.max(), resolution='c')
#m = Basemap(projection='ortho',lat_0=lats.min(),lon_0=lons.min(), resolution='l')
lon, lat =np.meshgrid(lons, lats)
XX, YY=m(lon,lat)


cmap = plt.cm.Reds
#cmap = plt.cm.RdBu_r
#norm = mpl.colors.Normalize(vmin=5, vmax=10)

clevs = np.arange(990,1010,1)
cs1 = m.contour(XX,YY,p,clevs,linewidths=1.0,colors='k',animated=True)
#cs2 = m.contourf(XX,YY,p,clevs,cmap=cmap,animated=True)


urot,vrot,x,y = m.rotate_vector(u[:,:],v[:,:],lons[:],lats[:],returnxy=True)
Q = m.quiver(lon,lat,u,v,latlon=True) #or specify, e.g., width=0.003, scale=400)
qk = plt.quiverkey(Q, 0.95, 1.05, 25, '25 m/s', labelpos='W')


###########################################################################################################
#ugrid,newlons = shiftgrid(180.,u,lons,start=False)
#vgrid,newlons = shiftgrid(180.,v,lons,start=False)
# transform vectors to projection grid.
#uproj,vproj,xx,yy = m.transform_vector(ugrid,vgrid,newlons,latitudes,31,31,returnxy=True,masked=True)
#Q = m.quiver(XX,YY,u,v,scale=9000)
# make quiver key.
#qk = plt.quiverkey(Q, 0.1, 0.1, 20, '20 m/s', labelpos='W')
#m.barbs(XX,YY,u,v, length=7, color='red')
##############################################################################################################

#m.drawcoastlines(linewidth=0.25) ; m.drawcountries(linewidth=0.25); m.drawstates(linewidth=0.25) 
m.readshapefile('subdiv/India_subdiv','ind',drawbounds=True, zorder=None, linewidth=1.0, color='k', antialiased=1, ax=None, default_encoding='utf-8')
m.drawmapboundary(fill_color='aqua') ; 
m.drawparallels(np.arange(0.,40.,10.), labels=[1,0,0,0])
m.drawmeridians(np.arange(60.,100.,10.), labels=[0,0,0,1])
#cbar=m.colorbar(cs1,location='bottom', pad='5%') ; cbar.set_label('hpa') ; 
plt.clabel(cs1, fontsize=9, inline=1)
#m.bluemarble() ;
m.fillcontinents(color='white',lake_color='white')
plt.title('6 hour Mean Sea Level Pressure(hpa)'); savefig('slp.png', dpi=100);
plt.show()

#############################################################################################################
quit




















