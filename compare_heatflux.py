'''

  Software which uses the MHW definition
  of Hobday et al. (2015) applied to
  select SST time series around the globe

'''

# Load required modules

import numpy as np
from scipy import io
from datetime import date
from Scientific.IO import NetCDF

from matplotlib import pyplot as plt
from matplotlib import dates as mdates
from matplotlib import colors
import mpl_toolkits.basemap as bm

import deseason as ds
import ecoliver as ecj

import marineHeatWaves as mhw

#
# CFSv2 air temperature
#

# load data
yr0 = 2013
file = '/data/MHWs/Tasmania_2015_2016/CFSv2/CFSv2_' + str(yr0) + '_heatflux.nc'
fileobj = NetCDF.NetCDFFile(file, 'r')
lon_CFS = fileobj.variables['x'].getValue()
lat_CFS = fileobj.variables['y'].getValue()
t_CFS = fileobj.variables['t'].getValue()
heatflux = fileobj.variables['heatflux'].getValue()
fileobj.close()

for yr in range(2014, 2016+1):
    file = '/data/MHWs/Tasmania_2015_2016/CFSv2/CFSv2_' + str(yr) + '_heatflux.nc'
    fileobj = NetCDF.NetCDFFile(file, 'r')
    t_CFS = np.append(t_CFS, fileobj.variables['t'].getValue())
    heatflux = np.append(heatflux, fileobj.variables['heatflux'].getValue(), axis=0)
    fileobj.close()

# Build up time and date vectors
t_CFS = np.floor(t_CFS + date(1990,1,1).toordinal()).astype(int)
year_CFS = []
month_CFS = []
dates_CFS = []
for tt in range(len(t_CFS)):
    year_CFS.append(date.fromordinal(t_CFS[tt]).year)
    month_CFS.append(date.fromordinal(t_CFS[tt]).month)
    dates_CFS.append(date.fromordinal(t_CFS[tt]))

year_CFS = np.array(year_CFS)
month_CFS = np.array(month_CFS)
T_CFS = len(t_CFS)

# Create monthly means and anomalies
year_start = 2013
year_end = 2016
Tl = (year_end-year_start+1)*12

heatflux_mth = np.NaN*np.zeros((lat_CFS.shape[0], lat_CFS.shape[1], Tl))
heatflux_mth_abs = np.NaN*np.zeros((lat_CFS.shape[0], lat_CFS.shape[1], Tl))

for ii in range(len(lon_CFS)):
    print ii, len(lon_CFS)
    for jj in range(len(lat_CFS)):
        if np.isnan(heatflux[:,jj,ii]).sum() == T_CFS:
            continue
        heatflux_ds, heatflux_s, heatflux_beta = ds.deseason_harmonic(heatflux[:,jj,ii], 2, 1460)
        tt = 0
        for yr in range(year_start, year_end+1):
            for mth in range(1,12+1):
                heatflux_mth[jj,ii,tt] = np.nanmean(heatflux_ds[(year_CFS==yr) * (month_CFS==mth)])
                heatflux_mth_abs[jj,ii,tt] = np.nanmean(heatflux[(year_CFS==yr) * (month_CFS==mth),jj,ii])
                tt += 1

#
# GODAS
#

yr0 = 2013
file = '/data/MHWs/Tasmania_2015_2016/GODAS/thflx.' + str(yr0) + '.nc'
fileobj = NetCDF.NetCDFFile(file, 'r')
lon_GOD = fileobj.variables['lon'].getValue()
lat_GOD = fileobj.variables['lat'].getValue()
heatflux_GOD_mth_abs = np.swapaxes(np.swapaxes(fileobj.variables['thflx'].getValue(), 0, 1), 1, 2)
fillValue = fileobj.variables['thflx']._FillValue
scale = fileobj.variables['thflx'].scale_factor
offset = fileobj.variables['thflx'].add_offset
time_GOD = np.floor(date(1800,1,1).toordinal() + fileobj.variables['time'].getValue() + 15).astype(int)
fileobj.close()

# Create monthly means and anomalies
year_start = 2013
year_end = 2016
Tl = (year_end-year_start+1)*12

for yr in range(2014, 2016+1):
    file = '/data/MHWs/Tasmania_2015_2016/GODAS/thflx.' + str(yr) + '.nc'
    fileobj = NetCDF.NetCDFFile(file, 'r')
    heatflux_GOD_mth_abs = np.append(heatflux_GOD_mth_abs, np.swapaxes(np.swapaxes(fileobj.variables['thflx'].getValue(), 0, 1), 1, 2), axis=2)
    time_GOD = np.append(time_GOD, np.floor(date(1800,1,1).toordinal() + fileobj.variables['time'].getValue() + 15).astype(int), 0)
    fileobj.close()

# Date vector
dates_GOD = []
T_GOD = len(time_GOD)
for tt in range(len(time_GOD)):
    dates_GOD.append(date.fromordinal(time_GOD[tt]))

# fillvalue, offset, scale
heatflux_GOD_mth_abs = heatflux_GOD_mth_abs.astype(float)
heatflux_GOD_mth_abs[heatflux_GOD_mth_abs==fillValue] = np.nan
heatflux_GOD_mth_abs = heatflux_GOD_mth_abs*scale + offset

heatflux_GOD_mth = np.NaN*np.zeros(heatflux_GOD_mth_abs.shape)

for ii in range(len(lon_GOD)):
    print ii, len(lon_GOD)
    for jj in range(len(lat_GOD)):
        if np.isnan(heatflux_GOD_mth_abs[jj,ii,:]).sum() == heatflux_GOD_mth_abs.shape[2]:
            continue
        heatflux_ds, heatflux_s, heatflux_beta = ds.deseason_harmonic(heatflux_GOD_mth_abs[jj,ii,:], 2, 12)
        heatflux_GOD_mth[jj,ii,:] = np.array(heatflux_ds)[:,0]

#
# OceanMAPS
# 

# load data
file = '/data/MHWs/Tasmania_2015_2016/OMAPS/OMAPS_2012_08_tse.nc'
fileobj = NetCDF.NetCDFFile(file, 'r')
lon = fileobj.variables['x'].getValue()
lat = fileobj.variables['y'].getValue()
z = fileobj.variables['z'].getValue()
fileobj.close()

cnt = 0;
for yr in range(2013, 2016+1):
    print yr
    for mth in range(12):
        file = '/data/MHWs/Tasmania_2015_2016/OMAPS/OMAPS_' + str(yr) + '_' + str(mth+1).zfill(2) + '_tse.nc'
        fileobj = NetCDF.NetCDFFile(file, 'r')
        if cnt == 0:
            t = fileobj.variables['t'].getValue()
            sfc_hflux = fileobj.variables['sfc_hflux'][:,:,:]
        else:
            t = np.append(t, fileobj.variables['t'].getValue())
            sfc_hflux = np.append(sfc_hflux, fileobj.variables['sfc_hflux'][:,:,:], axis=0)
        fileobj.close()
        cnt += 1

# Build up time and date vectors
t = np.floor(t).astype(int) - 366 # Convert from MATLAB datenum format
year = []
month = []
dates = []
for tt in range(len(t)):
    year.append(date.fromordinal(t[tt]).year)
    month.append(date.fromordinal(t[tt]).month)
    dates.append(date.fromordinal(t[tt]))

year = np.array(year)
month = np.array(month)
T = len(t)

# Create monthly means and anomalies
year_start = 2013
year_end = 2016
Tl = (year_end-year_start+1)*12

sfc_hflux_mth = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
sfc_hflux_mth_abs = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))

for ii in range(lon.shape[1]):
    print ii, lon.shape[1]
    for jj in range(lat.shape[0]):
        if ~(np.isnan(sfc_hflux[:,jj,ii]).sum() == T):
            sfc_hflux_ds, sfc_hflux_s, sfc_hflux_beta = ds.deseason_harmonic(sfc_hflux[:,jj,ii], 2, 365.25)
            tt = 0
            for yr in range(year_start, year_end+1):
                for mth in range(1,12+1):
                    sfc_hflux_mth[jj,ii,tt] = np.nanmean(sfc_hflux_ds[(year==yr) * (month==mth)])
                    sfc_hflux_mth_abs[jj,ii,tt] = np.nanmean(sfc_hflux[(year==yr) * (month==mth),jj,ii])
                    tt += 1

#
# Monthly maps
#

domain = [-46, 143, -36, 161]
labels_month = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
labels_year = np.arange(year_start, year_end+1)
wbgyr = colors.ListedColormap(np.genfromtxt('/home/ecoliver/Desktop/python_modules/cmaps/WhiteBlueGreenYellowRed.rgb')/255.)
llon_GOD, llat_GOD = np.meshgrid(lon_GOD, lat_GOD)

# Heat flux
fig = plt.figure(figsize=(11,8))
year = 2016
plt.clf()
cnt = 0
for tt in range((year-2013)*12, (year-2013+1)*12):
    cnt += 1
    # GODAS
    AX = plt.subplot(3,12,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])
    lonproj, latproj = proj(llon_GOD, llat_GOD)
    H = plt.contourf(lonproj, latproj, heatflux_GOD_mth_abs[:,:,tt], levels=np.arange(-1,1+0.1,0.1)*400, cmap=plt.cm.RdBu_r)
    plt.clim(-300,300)
    plt.title(labels_month[np.mod(tt,12)])
    if cnt == 1:
        plt.ylabel('GODAS')
    # CFSv2
    AX = plt.subplot(3,12,12 + cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])
    lonproj, latproj = proj(lon_CFS, lat_CFS)
    H = plt.contourf(lonproj, latproj, heatflux_mth_abs[:,:,tt], levels=np.arange(-1,1+0.1,0.1)*400, cmap=plt.cm.RdBu_r)
    plt.clim(-300,300)
    if cnt == 1:
        plt.ylabel('CFSv2')
    # OceanMAPS
    AX = plt.subplot(3,12,24 + cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])
    lonproj, latproj = proj(lon, lat)
    H = plt.contourf(lonproj, latproj, sfc_hflux_mth_abs[:,:,tt], levels=np.arange(-1,1+0.1,0.1)*400, cmap=plt.cm.RdBu_r)
    plt.clim(-300,300)
    if cnt == 1:
        plt.ylabel('OceanMAPS')

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[W m$^{-2}$]')

plt.savefig('mhw_properties/MHW_SFCHFLUX_Comparisons_Monthly_' + str(year) +'.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

# Heat flux anomalies
fig = plt.figure(figsize=(11,8))
year = 2016
plt.clf()
cnt = 0
for tt in range((year-2013)*12, (year-2013+1)*12):
    cnt += 1
    # GODAS
    AX = plt.subplot(3,12,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])
    lonproj, latproj = proj(llon_GOD, llat_GOD)
    H = plt.contourf(lonproj, latproj, heatflux_GOD_mth[:,:,tt], levels=np.arange(-1,1+0.1,0.1)*150, cmap=plt.cm.RdBu_r)
    plt.clim(-100,100)
    plt.title(labels_month[np.mod(tt,12)])
    if cnt == 1:
        plt.ylabel('GODAS')
    # CFSv2
    AX = plt.subplot(3,12,12 + cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])
    lonproj, latproj = proj(lon_CFS, lat_CFS)
    H = plt.contourf(lonproj, latproj, heatflux_mth[:,:,tt], levels=np.arange(-1,1+0.1,0.1)*150, cmap=plt.cm.RdBu_r)
    plt.clim(-100,100)
    if cnt == 1:
        plt.ylabel('CFSv2')
    # OceanMAPS
    AX = plt.subplot(3,12,24 + cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i') 
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])
    lonproj, latproj = proj(lon, lat)
    H = plt.contourf(lonproj, latproj, sfc_hflux_mth[:,:,tt], levels=np.arange(-1,1+0.1,0.1)*150, cmap=plt.cm.RdBu_r)
    plt.clim(-100,100)
    if cnt == 1:
        plt.ylabel('OceanMAPS')

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[W m$^{-2}$]')

plt.savefig('mhw_properties/MHW_SFCHFLUX_Comparisons_MonthlyAnom_' + str(year) +'.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

#
# Daily, spatially-averaged time series
#

locations = {}
#locations['lon'] = [147, 155]
#locations['lat'] = [-45, -37]
locations['lon'] = [147, 157]
locations['lat'] = [-46, -39]

sfc_hflux_ts = np.zeros((T))
sfc_hflux_ts_ds = np.zeros((T))
heatflux_ts = np.zeros((T_CFS))
heatflux_ts_ds = np.zeros((T_CFS))
heatflux_GOD_ts = np.zeros((T_GOD))
heatflux_GOD_ts_ds = np.zeros((T_GOD))

# OceanMAPS
i1 = np.where(lon[0,:] > locations['lon'][0])[0][0] - 1
i2 = np.where(lon[0,:] > locations['lon'][1])[0][0]
j1 = np.where(lat[:,0] > locations['lat'][0])[0][0] - 1
j2 = np.where(lat[:,0] > locations['lat'][1])[0][0]
#sfc_hflux_ts = np.nanmean(np.nanmean(sfc_hflux[:,j1:j2+1,:][:,:,i1:i2+1], axis=2), axis=1)
#sfc_hflux_ts_ds = np.array(ds.deseason_harmonic(sfc_hflux_ts, 2, 365.25)[0])
sfc_hflux_ts = np.nanmean(np.nanmean(sfc_hflux_mth_abs[j1:j2+1,:,:][:,i1:i2+1,:], axis=1), axis=0)
sfc_hflux_ts_ds = np.array(ds.deseason_harmonic(sfc_hflux_ts, 2, 12)[0])

# CFSv2
i1 = np.where(lon_CFS[0,:] > locations['lon'][0])[0][0] - 1
i2 = np.where(lon_CFS[0,:] > locations['lon'][1])[0][0]
j1 = np.where(lat_CFS[:,0] > locations['lat'][0])[0][0] - 1
j2 = np.where(lat_CFS[:,0] > locations['lat'][1])[0][0]
#heatflux_ts = np.nanmean(np.nanmean(heatflux[:,j1:j2+1,:][:,:,i1:i2+1], axis=2), axis=1)
#heatflux_ts_ds = np.array(ds.deseason_harmonic(heatflux_ts, 2, 1460)[0])
heatflux_ts = np.nanmean(np.nanmean(heatflux_mth_abs[j1:j2+1,:,:][:,i1:i2+1,:], axis=1), axis=0)
heatflux_ts_ds = np.array(ds.deseason_harmonic(heatflux_ts, 2, 12)[0])

# GODAS
i1 = np.where(lon_GOD > locations['lon'][0])[0][0] - 1
i2 = np.where(lon_GOD > locations['lon'][1])[0][0]
j1 = np.where(lat_GOD > locations['lat'][0])[0][0] - 1
j2 = np.where(lat_GOD > locations['lat'][1])[0][0]
heatflux_GOD_ts = np.nanmean(np.nanmean(heatflux_GOD_mth_abs[j1:j2+1,:,:][:,i1:i2+1,:], axis=1), axis=0)
heatflux_GOD_ts_ds = np.array(ds.deseason_harmonic(heatflux_GOD_ts, 2, 12)[0])

#
# Plot some timeseries
#

plt.figure(figsize=(16,8))

ts = date(2013,1,1).toordinal()
te = date(2016,5,1).toordinal()

plt.clf()
plt.subplot(2,1,1)
#plt.plot(dates, sfc_hflux_ts, '-', color='0.5', linewidth=1)
#plt.plot(dates, ecj.runavg(sfc_hflux_ts, 61), 'k-', linewidth=2)
#plt.plot(dates_CFS, heatflux_ts, 'b-')
plt.plot(dates_GOD, sfc_hflux_ts[:T_GOD], 'r-o')
plt.plot(dates_GOD, heatflux_ts[:T_GOD], 'b-o')
plt.plot(dates_GOD, heatflux_GOD_ts, 'k-o')
plt.xlim(ts, te)
plt.ylim(-220, 220)
plt.grid()
plt.ylabel('Surface heat flux [W m$^{-2}$]')
plt.title('Product comparison (2013-2016)')
plt.subplot(2,1,2)
plt.plot(dates_GOD, sfc_hflux_ts_ds[:T_GOD], 'r-o')
plt.plot(dates_GOD, heatflux_ts_ds[:T_GOD], 'b-o')
plt.plot(dates_GOD, heatflux_GOD_ts_ds, 'k-o')
plt.xlim(ts, te)
#plt.ylim(-300, 300)
plt.grid()
plt.ylabel('Surface heat flux anomaly [W m$^{-2}$]')
plt.legend(['OceanMAPS', 'CFSv2', 'GODAS'], loc='lower left')

plt.savefig('mhw_properties/MHW_SFCHFLUX_Comparison_timeSeries.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

