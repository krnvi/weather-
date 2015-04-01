#!/usr/bin/env python

import sys ; import numpy as np ; import scipy as sp ; import datetime
from grads import * ; 
import matplotlib.pyplot as plt; from pylab import * ; from matplotlib import *

from mpl_toolkits.basemap import Basemap, cm, shiftgrid, maskoceans, interp, shapefile
from scipy.io import netcdf 
from netCDF4 import Dataset


nc_file='wrfpost_TMPsfc.nc'
fh     =Dataset(nc_file,'r')
lons   =fh.variables['lon'][:]
lats   =fh.variables['lat'][:]
time   =fh.variables['time'][:]
temp   =fh.variables['tmpsfc'][0:95,:,:]
fh.close()
tmp_day1=temp.mean(axis=0)-273.15


########################################################################

fig=plt.figure(figsize=(8,8)); ax = fig.add_axes([0.1,0.1,0.8,0.8])

m=Basemap(projection='merc',lat_ts=10, llcrnrlon=lons.min(),urcrnrlon=lons.max()-3, llcrnrlat=lats.min()+3,urcrnrlat=lats.max(), resolution='c')
#m = Basemap(projection='ortho',lat_0=lats.min(),lon_0=lons.min(), resolution='l')
lon, lat =np.meshgrid(lons, lats)
XX, YY=m(lon,lat)
my_cmap = matplotlib.cm.get_cmap('rainbow')
my_cmap.set_under('w')

#cmap = cm.s3pcpn_l
cmap = cm.StepSeq
cmap = plt.get_cmap('gist_ncar')
#norm = mpl.colors.Normalize(vmin=5, vmax=10)

m.pcolormesh(XX,YY,tmp_day1,shading='flat', cmap=cmap, vmin=16, vmax=35);


#m.drawcoastlines(linewidth=0.25) ; m.drawcountries(linewidth=0.25); m.drawstates(linewidth=0.25) 
m.readshapefile('subdiv/India_subdiv','ind',drawbounds=True, zorder=None, linewidth=1.0, color='k', antialiased=1, ax=None, default_encoding='utf-8')
#m.drawmapboundary(fill_color='aqua') ; 
m.drawparallels(np.arange(0.,40.,10.), labels=[1,0,0,0])
m.drawmeridians(np.arange(60.,100.,10.), labels=[0,0,0,1])
m.drawlsmask(land_color='grey',ocean_color='aqua',lakes=True)
cbar=m.colorbar(location='bottom', pad='5%') ; cbar.set_label('degC') ; #m.etopo()
#m.drawlsmask(land_color='grey',ocean_color='aqua',lakes=True)
m.fillcontinents(color='white',lake_color='white')
plt.title('24 hour Sea Surface Temperature(degC)'); savefig('sst_philin.png', dpi=100);
plt.close()

#############################################################################################################
quit




















