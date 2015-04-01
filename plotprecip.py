#!/usr/bin/env python

import sys ; import numpy as np ; import scipy as sp ; import datetime
from grads import * ; 
import matplotlib.pyplot as plt; from pylab import * ; from matplotlib import *

from mpl_toolkits.basemap import Basemap, cm, shiftgrid, maskoceans, interp, shapefile
from scipy.io import netcdf 
from netCDF4 import Dataset


nc_file='wrfpost_APCPsfc.nc'
fh     =Dataset(nc_file,'r')
lons   =fh.variables['lon'][:]
lats   =fh.variables['lat'][:]
time   =fh.variables['time'][:]
rain   =fh.variables['apcpsfc'][0:95,:,:]
fh.close()
rf_day1=rain.sum(axis=0)


########################################################################

fig=plt.figure(figsize=(8,8)); ax = fig.add_axes([0.1,0.1,0.8,0.8])
m=Basemap(projection='merc',lat_ts=10, llcrnrlon=lons.min(),urcrnrlon=lons.max()-3, llcrnrlat=lats.min()+3,urcrnrlat=lats.max(), resolution='c')
#m = Basemap(projection='ortho',lat_0=lats.min(),lon_0=lons.min(), resolution='l')
lon, lat =np.meshgrid(lons, lats)
XX, YY=m(lon,lat)
my_cmap = matplotlib.cm.get_cmap('rainbow')
my_cmap.set_under('w')

cmap = cm.s3pcpn
cmap = cm.s3pcpn_l
cmap = cm.StepSeq
#cmap = cm.GMT_drywet
#norm = mpl.colors.Normalize(vmin=5, vmax=10)

#m.pcolormesh(XX,YY,rf_day1,shading='flat', cmap=my_cmap);
clevs=[0,1,2.5,5,7.5,10,15,20,30,40,50,70,100,150,200,250,300,400,500,600,750]
cs = m.contourf(XX,YY,rf_day1,clevs,cmap=cmap)

#m.drawcoastlines(linewidth=0.25) ; m.drawcountries(linewidth=0.25); m.drawstates(linewidth=0.25) 
m.readshapefile('subdiv/India_subdiv','ind',drawbounds=True, zorder=None, linewidth=1.0, color='k', antialiased=1, ax=None, default_encoding='utf-8')
m.drawmapboundary(fill_color='aqua') ; 
m.drawparallels(np.arange(0.,40.,10.), labels=[1,0,0,0])
m.drawmeridians(np.arange(60.,100.,10.), labels=[0,0,0,1])
cbar=m.colorbar(cs,location='bottom', pad='5%') ; cbar.set_label('mm') ; m.etopo()

plt.title('24 hour Accumalated Rainfall(mm/day)'); savefig('precipitation.png', dpi=100);
plt.close()

#############################################################################################################
quit




















