'''

  Software which uses the MHW definition
  of Hobday et al. (2015) applied to
  select SST time series around the globe

'''

# Load required modules

import numpy as np
from scipy import io
from datetime import date
from netCDF4 import Dataset

from matplotlib import pyplot as plt
from matplotlib import dates as mdates
from matplotlib import colors
import mpl_toolkits.basemap as bm

import deseason as ds
import ecoliver as ecj

import marineHeatWaves as mhw

# Some basic parameters

coldSpells = False # Detect coldspells instead of heatwaves
seasonalMeans = False
col_clim = '0.25'
col_thresh = 'g-'
if coldSpells:
    mhwname = 'MCS'
    mhwfullname = 'coldspell'
    col_evMax = (0, 102./255, 204./255)
    col_ev = (153./255, 204./255, 1)
    col_bar = (0.5, 0.5, 1)
    cmap_i = plt.cm.YlGnBu_r
else:
    mhwname = 'MHW'
    mhwfullname = 'heatwave'
    col_evMax = 'r'
    col_ev = (1, 0.6, 0.5)
    col_bar = (1, 0.5, 0.5)
    cmap_i = plt.cm.hot_r

#
# OceanMAPS
# 

# load data
file = '/data/MHWs/Tasmania_2015_2016/OMAPS/OMAPS_2012_08_tse.nc'
fileobj = Dataset(file, 'r')
lon = fileobj.variables['x'][:]
lat = fileobj.variables['y'][:]
z = fileobj.variables['z'][:]
fileobj.close()
file = '/data/MHWs/Tasmania_2015_2016/OMAPS/OMAPS_2012_08_uv.nc'
fileobj = Dataset(file, 'r')
lon_uv = fileobj.variables['x'][:]
lat_uv = fileobj.variables['y'][:]
fileobj.close()

cnt = 0;
for yr in range(2013, 2016+1):
    print yr
    for mth in range(12):
        file = '/data/MHWs/Tasmania_2015_2016/OMAPS/OMAPS_' + str(yr) + '_' + str(mth+1).zfill(2) + '_tse.nc'
        fileobj = Dataset(file, 'r')
        if cnt == 0:
            t = fileobj.variables['t'][:]
            temp_sfc = fileobj.variables['temp'][:,-1,:,:]
            temp_100 = 0.5*fileobj.variables['temp'][:,-14,:,:] + 0.5*fileobj.variables['temp'][:,-15,:,:]
            temp_200 = 0.5*fileobj.variables['temp'][:,-24,:,:] + 0.5*fileobj.variables['temp'][:,-25,:,:]
            temp_350 = 0.5*fileobj.variables['temp'][:,-30,:,:] + 0.5*fileobj.variables['temp'][:,-31,:,:]
            temp_500 = (45./60)*fileobj.variables['temp'][:,-33,:,:] + (15./60)*fileobj.variables['temp'][:,-34,:,:]
            eta = fileobj.variables['eta'][:,:,:]
            sfc_hflux = fileobj.variables['sfc_hflux'][:,:,:]
        else:
            t = np.append(t, fileobj.variables['t'][:])
            temp_sfc = np.append(temp_sfc, fileobj.variables['temp'][:,-1,:,:], axis=0)
            temp_100 = np.append(temp_100, 0.5*fileobj.variables['temp'][:,-14,:,:] + 0.5*fileobj.variables['temp'][:,-15,:,:], axis=0)
            temp_200 = np.append(temp_200, 0.5*fileobj.variables['temp'][:,-24,:,:] + 0.5*fileobj.variables['temp'][:,-25,:,:], axis=0)
            temp_350 = np.append(temp_350, 0.5*fileobj.variables['temp'][:,-30,:,:] + 0.5*fileobj.variables['temp'][:,-31,:,:], axis=0)
            temp_500 = np.append(temp_500, (45./60)*fileobj.variables['temp'][:,-33,:,:] + (15./60)*fileobj.variables['temp'][:,-34,:,:], axis=0)
            eta = np.append(eta, fileobj.variables['eta'][:,:,:], axis=0)
            sfc_hflux = np.append(sfc_hflux, fileobj.variables['sfc_hflux'][:,:,:], axis=0)
        fileobj.close()
        file = '/data/MHWs/Tasmania_2015_2016/OMAPS/OMAPS_' + str(yr) + '_' + str(mth+1).zfill(2) + '_uv.nc'
        fileobj = Dataset(file, 'r')
        if cnt == 0:
            tau_x = fileobj.variables['tau_x'][:,:,:]
            tau_y = fileobj.variables['tau_y'][:,:,:]
        else:
            tau_x = np.append(tau_x, fileobj.variables['tau_x'][:,:,:], axis=0)
            tau_y = np.append(tau_y, fileobj.variables['tau_y'][:,:,:], axis=0)
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

temp_sfc_mth = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
temp_sfc_mth_abs = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
temp_100_mth = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
temp_100_mth_abs = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
temp_200_mth = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
temp_200_mth_abs = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
temp_350_mth = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
temp_350_mth_abs = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
temp_500_mth = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
temp_500_mth_abs = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
eta_mth = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
eta_mth_abs = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
sfc_hflux_mth = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
sfc_hflux_mth_abs = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
tau_x_mth = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
tau_x_mth_abs = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
tau_y_mth = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))
tau_y_mth_abs = np.NaN*np.zeros((lat.shape[0], lat.shape[1], Tl))

for ii in range(lon.shape[1]):
    print ii, lon.shape[1]
    for jj in range(lat.shape[0]):
        if ~(np.isnan(temp_sfc[:,jj,ii]).sum() == T):
            temp_sfc_ds, temp_sfc_s, temp_sfc_beta = ds.deseason_harmonic(temp_sfc[:,jj,ii], 2, 365.25)
            tt = 0
            for yr in range(year_start, year_end+1):
                for mth in range(1,12+1):
                    temp_sfc_mth[jj,ii,tt] = np.nanmean(temp_sfc_ds[(year==yr) * (month==mth)])
                    temp_sfc_mth_abs[jj,ii,tt] = np.nanmean(temp_sfc[(year==yr) * (month==mth),jj,ii])
                    tt += 1
        if ~(np.isnan(temp_100[:,jj,ii]).sum() == T):
            temp_100_ds, temp_100_s, temp_100_beta = ds.deseason_harmonic(temp_100[:,jj,ii], 2, 365.25)
            tt = 0
            for yr in range(year_start, year_end+1):
                for mth in range(1,12+1):
                    temp_100_mth[jj,ii,tt] = np.nanmean(temp_100_ds[(year==yr) * (month==mth)])
                    temp_100_mth_abs[jj,ii,tt] = np.nanmean(temp_100[(year==yr) * (month==mth),jj,ii])
                    tt += 1
        if ~(np.isnan(temp_200[:,jj,ii]).sum() == T):
            temp_200_ds, temp_200_s, temp_200_beta = ds.deseason_harmonic(temp_200[:,jj,ii], 2, 365.25)
            tt = 0
            for yr in range(year_start, year_end+1):
                for mth in range(1,12+1):
                    temp_200_mth[jj,ii,tt] = np.nanmean(temp_200_ds[(year==yr) * (month==mth)])
                    temp_200_mth_abs[jj,ii,tt] = np.nanmean(temp_200[(year==yr) * (month==mth),jj,ii])
                    tt += 1
        if ~(np.isnan(temp_350[:,jj,ii]).sum() == T):
            temp_350_ds, temp_350_s, temp_350_beta = ds.deseason_harmonic(temp_350[:,jj,ii], 2, 365.25)
            tt = 0
            for yr in range(year_start, year_end+1):
                for mth in range(1,12+1):
                    temp_350_mth[jj,ii,tt] = np.nanmean(temp_350_ds[(year==yr) * (month==mth)])
                    temp_350_mth_abs[jj,ii,tt] = np.nanmean(temp_350[(year==yr) * (month==mth),jj,ii])
                    tt += 1
        if ~(np.isnan(temp_500[:,jj,ii]).sum() == T):
            temp_500_ds, temp_500_s, temp_500_beta = ds.deseason_harmonic(temp_500[:,jj,ii], 2, 365.25)
            tt = 0
            for yr in range(year_start, year_end+1):
                for mth in range(1,12+1):
                    temp_500_mth[jj,ii,tt] = np.nanmean(temp_500_ds[(year==yr) * (month==mth)])
                    temp_500_mth_abs[jj,ii,tt] = np.nanmean(temp_500[(year==yr) * (month==mth),jj,ii])
                    tt += 1
        if ~(np.isnan(eta[:,jj,ii]).sum() == T):
            eta_ds, eta_s, eta_beta = ds.deseason_harmonic(eta[:,jj,ii], 2, 365.25)
            tt = 0
            for yr in range(year_start, year_end+1):
                for mth in range(1,12+1):
                    eta_mth[jj,ii,tt] = np.nanmean(eta_ds[(year==yr) * (month==mth)])
                    eta_mth_abs[jj,ii,tt] = np.nanmean(eta[(year==yr) * (month==mth),jj,ii])
                    tt += 1
        if ~(np.isnan(sfc_hflux[:,jj,ii]).sum() == T):
            sfc_hflux_ds, sfc_hflux_s, sfc_hflux_beta = ds.deseason_harmonic(sfc_hflux[:,jj,ii], 2, 365.25)
            tt = 0
            for yr in range(year_start, year_end+1):
                for mth in range(1,12+1):
                    sfc_hflux_mth[jj,ii,tt] = np.nanmean(sfc_hflux_ds[(year==yr) * (month==mth)])
                    sfc_hflux_mth_abs[jj,ii,tt] = np.nanmean(sfc_hflux[(year==yr) * (month==mth),jj,ii])
                    tt += 1
        if ~(np.isnan(tau_x[:,jj,ii]).sum() == T):
            tau_x_ds, tau_x_s, tau_x_beta = ds.deseason_harmonic(tau_x[:,jj,ii], 2, 365.25)
            tt = 0
            for yr in range(year_start, year_end+1):
                for mth in range(1,12+1):
                    tau_x_mth[jj,ii,tt] = np.nanmean(tau_x_ds[(year==yr) * (month==mth)])
                    tau_x_mth_abs[jj,ii,tt] = np.nanmean(tau_x[(year==yr) * (month==mth),jj,ii])
                    tt += 1
        if ~(np.isnan(tau_y[:,jj,ii]).sum() == T):
            tau_y_ds, tau_y_s, tau_y_beta = ds.deseason_harmonic(tau_y[:,jj,ii], 2, 365.25)
            tt = 0
            for yr in range(year_start, year_end+1):
                for mth in range(1,12+1):
                    tau_y_mth[jj,ii,tt] = np.nanmean(tau_y_ds[(year==yr) * (month==mth)])
                    tau_y_mth_abs[jj,ii,tt] = np.nanmean(tau_y[(year==yr) * (month==mth),jj,ii])
                    tt += 1

# Daily, spatially-averaged time series
locations = {}
#locations['lon'] = [147, 155]
#locations['lat'] = [-45, -37]
locations['lon'] = [147, 157]
locations['lat'] = [-46, -39]

temp_sfc_ts = np.zeros((T))
temp_sfc_ts_ds = np.zeros((T))
temp_100_ts = np.zeros((T))
temp_100_ts_ds = np.zeros((T))
temp_200_ts = np.zeros((T))
temp_200_ts_ds = np.zeros((T))
temp_350_ts = np.zeros((T))
temp_350_ts_ds = np.zeros((T))
temp_500_ts = np.zeros((T))
temp_500_ts_ds = np.zeros((T))
sfc_hflux_ts = np.zeros((T))
sfc_hflux_ts_ds = np.zeros((T))
tau_x_ts = np.zeros((T))
tau_x_ts_ds = np.zeros((T))
tau_y_ts = np.zeros((T))
tau_y_ts_ds = np.zeros((T))
i1 = np.where(lon[0,:] > locations['lon'][0])[0][0] - 1
i2 = np.where(lon[0,:] > locations['lon'][1])[0][0]
j1 = np.where(lat[:,0] > locations['lat'][0])[0][0] - 1
j2 = np.where(lat[:,0] > locations['lat'][1])[0][0]

temp_sfc_ts = np.nanmean(np.nanmean(temp_sfc[:,j1:j2+1,:][:,:,i1:i2+1], axis=2), axis=1)
temp_100_ts = np.nanmean(np.nanmean(temp_100[:,j1:j2+1,:][:,:,i1:i2+1], axis=2), axis=1)
temp_200_ts = np.nanmean(np.nanmean(temp_200[:,j1:j2+1,:][:,:,i1:i2+1], axis=2), axis=1)
temp_350_ts = np.nanmean(np.nanmean(temp_350[:,j1:j2+1,:][:,:,i1:i2+1], axis=2), axis=1)
temp_500_ts = np.nanmean(np.nanmean(temp_500[:,j1:j2+1,:][:,:,i1:i2+1], axis=2), axis=1)
sfc_hflux_ts = np.nanmean(np.nanmean(sfc_hflux[:,j1:j2+1,:][:,:,i1:i2+1], axis=2), axis=1)
tau_x_ts = np.nanmean(np.nanmean(tau_x[:,j1:j2+1,:][:,:,i1:i2+1], axis=2), axis=1)
tau_y_ts = np.nanmean(np.nanmean(tau_y[:,j1:j2+1,:][:,:,i1:i2+1], axis=2), axis=1)

temp_sfc_ts_ds = np.array(ds.deseason_harmonic(temp_sfc_ts, 2, 365.25)[0])
temp_100_ts_ds = np.array(ds.deseason_harmonic(temp_100_ts, 2, 365.25)[0])
temp_200_ts_ds = np.array(ds.deseason_harmonic(temp_200_ts, 2, 365.25)[0])
temp_350_ts_ds = np.array(ds.deseason_harmonic(temp_350_ts, 2, 365.25)[0])
temp_500_ts_ds = np.array(ds.deseason_harmonic(temp_500_ts, 2, 365.25)[0])
sfc_hflux_ts_ds = np.array(ds.deseason_harmonic(sfc_hflux_ts, 2, 365.25)[0])
tau_x_ts_ds = np.array(ds.deseason_harmonic(tau_x_ts, 2, 365.25)[0])
tau_y_ts_ds = np.array(ds.deseason_harmonic(tau_y_ts, 2, 365.25)[0])

# Load obs
#outfile = '/data/MHWs/Tasmania_2015_2016/NOAAOISST/mhw_SEAus_NOAAOISST_HadSST'
outfile = '/data/MHWs/Tasmania_2015_2016/NOAAOISST/' + str(locations['lon'][0]) + '_' + str(locations['lon'][1]) + '_' + str(locations['lat'][0]) + '_' + str(locations['lat'][1]) + '/mhw_SEAus_NOAAOISST_HadSST'
data = np.load(outfile + '.npz')
t_sst = data['t']
dates_sst = data['dates']
year_sst = data['year']
month_sst = data['month']
day_sst = data['day']

tt1 = np.where(t_sst == t[0])[0][0]
t_sst, dates_sst, T_sst, year_sst, month_sst, day_sst, doy_sst = ecj.timevector([dates_sst[tt1].year, dates_sst[tt1].month, dates_sst[tt1].day], [dates_sst[-1].year, dates_sst[-1].month, dates_sst[-1].day])

sst = data['sst'][tt1:]
#sst_had = data['sst_had']
#sst_had_clim = data['sst_had_clim']
sst_ds = np.array(ds.deseason_harmonic(sst, 2, 365.25)[0])

#
# Plot some timeseries
#

plt.figure(figsize=(16,8))

ts = date(2013,1,1).toordinal()
te = date(2016,7,1).toordinal()

plt.clf()
plt.subplot(2,1,1)
## Plot SST, seasonal cycle, threshold, etc
plt.plot(dates_sst, sst, color='0.5', linewidth=2)
plt.plot(dates, temp_sfc_ts, 'k-', linewidth=2)
plt.plot(dates, temp_100_ts, 'b-', linewidth=2)
plt.plot(dates, temp_200_ts, 'c-', linewidth=2)
plt.plot(dates, temp_350_ts, 'm-', linewidth=2)
plt.plot(dates, temp_500_ts, 'r-', linewidth=2)
plt.xlim(ts, te)
#plt.grid()
plt.ylabel('Temperature [$^\circ$C]')
plt.title('(A) OceanMAPS (2013-2016)')
plt.legend(['Observed SST', 'Temp. (2.5 m)', 'Temp. (100 m)', 'Temp. (200 m)', 'Temp. (350 m)', 'Temp. (500 m)'], loc='lower left', fontsize=12)
plt.subplot(2,1,2)
## Plot SST anomalies
plt.plot(dates_sst, sst_ds, color='0.5', linewidth=2)
plt.plot(dates, temp_sfc_ts_ds, 'k-', linewidth=2)
plt.plot(dates, temp_100_ts_ds, 'b-', linewidth=2)
plt.plot(dates, temp_200_ts_ds, 'c-', linewidth=2)
plt.plot(dates, temp_350_ts_ds, 'm-', linewidth=2)
plt.plot(dates, temp_500_ts_ds, 'r-', linewidth=2)
plt.xlim(ts, te)
#plt.ylim(-1.5, 2.5)
#plt.grid()
plt.ylabel('Temperature Anomaly [$^\circ$C]')
plt.title('(B)')

plt.savefig('../../documents/14_Tasmania_2015_2016/figures/' + mhwname + '_timeSeries_OMAPS_Temp.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

plt.clf()
plt.subplot(2,1,1)
## Plot SST, seasonal cycle, threshold, etc
plt.plot(dates, sfc_hflux_ts, '-', color='0.5', linewidth=1)
plt.plot(dates, ecj.runavg(sfc_hflux_ts, 61), 'k-', linewidth=2)
plt.xlim(ts, te)
plt.ylim(-250, 250)
plt.grid()
plt.ylabel('Surface heat flux [W m$^{-2}$]')
plt.title('OceanMAPS (2013-2016)')
plt.subplot(2,1,2)
## Plot SST anomalies
plt.plot(dates, sfc_hflux_ts_ds, '-', color='0.5', linewidth=1)
plt.plot(dates, ecj.runavg(sfc_hflux_ts_ds[:,0], 61), 'k-', linewidth=2)
plt.xlim(ts, te)
plt.ylim(-100, 100)
plt.grid()
plt.ylabel('Surface heat flux anomaly [W m$^{-2}$]')
plt.legend(['Surface heat flux', '61-day running mean'], loc='upper left')

plt.savefig('mhw_properties/' + mhwname + '_timeSeries_OMAPS_SFCHFLUX.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

plt.clf()
plt.subplot(2,1,1)
## Plot SST, seasonal cycle, threshold, etc
plt.plot(dates, np.sqrt(tau_x_ts**2+tau_y_ts**2), '-', color='0.5', linewidth=1)
plt.plot(dates, ecj.runavg(np.sqrt(tau_x_ts**2+tau_y_ts**2), 61), 'k-', linewidth=2)
plt.xlim(ts, te)
plt.ylim(0, 0.5)
plt.grid()
plt.ylabel('Surface wind stress [N m$^{-2}$]')
plt.title('OceanMAPS (2013-2016)')
plt.subplot(2,1,2)
## Plot SST anomalies
plt.plot(dates, np.sqrt(tau_x_ts_ds**2+tau_y_ts_ds**2), '-', color='0.5', linewidth=1)
plt.plot(dates, ecj.runavg(np.sqrt(tau_x_ts_ds**2+tau_y_ts_ds**2)[:,0], 61), 'k-', linewidth=2)
plt.xlim(ts, te)
plt.ylim(0, 0.4)
plt.grid()
plt.ylabel('Surface wind stress anomaly [N m$^{-2}$]')
plt.legend(['Surface wind stress', '61-day running mean'], loc='upper left')

plt.savefig('mhw_properties/' + mhwname + '_timeSeries_OMAPS_TAU.png', bbox_inches='tight', pad_inches=0.5, dpi=150)



#
# Monthly maps
#


domain = [-46, 143, -36, 161]
labels_month = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
labels_year = np.arange(year_start, year_end+1)
wbgyr = colors.ListedColormap(np.genfromtxt('/home/ecoliver/Desktop/python_modules/cmaps/WhiteBlueGreenYellowRed.rgb')/255.)

# SSTAs
fig = plt.figure(figsize=(11,8))
plt.clf()
cnt = 0
for tt in range(24+9-1, 36+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(lon, lat)
    H = plt.contourf(lonproj, latproj, temp_sfc_mth[:,:,tt], levels=range(-5,5+1), cmap=plt.cm.RdBu_r)
    plt.clim(-4, 4)
    plt.contour(lonproj, latproj, eta_mth[:,:,tt], levels=np.arange(-1,1+0.1,0.1), colors='k')
    plt.title(labels_month[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_OMAPS_SSTAnom_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

# MHW Intensities
mhw_map = temp_sfc_mth.copy()
mhw_map[mhw_map<=0.] = 0.
plt.clf()
cnt = 0
for tt in range(24+9-1, 36+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(lon, lat)
    H = plt.contourf(lonproj, latproj, mhw_map[:,:,tt], levels=[0.5,1,1.5,2,2.5,3,4,5], cmap=plt.cm.gist_heat_r)
    plt.clim(0.5,4)
    plt.contour(lonproj, latproj, eta_mth[:,:,tt], levels=np.arange(-1,1+0.1,0.1), colors='k')
    plt.title(labels_month[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_OMAPS_MHWInt_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

# Temp' (200m)
plt.clf()
cnt = 0
for tt in range(24+9-1, 36+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(lon, lat)
    H = plt.contourf(lonproj, latproj, temp_200_mth[:,:,tt], levels=range(-5,5+1), cmap=plt.cm.RdBu_r)
    plt.clim(-4, 4)
    plt.contour(lonproj, latproj, eta_mth[:,:,tt], levels=np.arange(-1,1+0.1,0.1), colors='k')
    plt.title(labels_month[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_OMAPS_T200Anom_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

# Temp' (500m)
plt.clf()
cnt = 0
for tt in range(24+9-1, 36+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(lon, lat)
    H = plt.contourf(lonproj, latproj, temp_500_mth[:,:,tt], levels=range(-5,5+1), cmap=plt.cm.RdBu_r)
    plt.clim(-4, 4)
    plt.contour(lonproj, latproj, eta_mth[:,:,tt], levels=np.arange(-1,1+0.1,0.1), colors='k')
    plt.title(labels_month[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_OMAPS_T500Anom_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

# Heat flux
plt.clf()
cnt = 0
for tt in range(24+9-1, 36+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(lon, lat)
    H = plt.contourf(lonproj, latproj, sfc_hflux_mth[:,:,tt], levels=np.arange(-1,1+0.2,0.2)*400, cmap=plt.cm.RdBu_r)
    #H = plt.contourf(lonproj, latproj, sfc_hflux_mth[:,:,tt])#, levels=np.arange(-1,1+0.2,0.2), cmap=plt.cm.RdBu_r)
    plt.title(labels_month[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[W m$^{-2}$]')

plt.savefig('mhw_properties/' + mhwname + '_OMAPS_SFCHFLUXAnom_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

# Wind Stress
plt.clf()
d = 5
cnt = 0
for tt in range(24+9-1, 36+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(lon_uv, lat_uv)
    H = plt.contourf(lonproj, latproj, np.sqrt(tau_x_mth_abs[:,:,tt]**2 + tau_y_mth_abs[:,:,tt]**2), levels=np.arange(0.2,1+0.1,0.1)*0.25, cmap=plt.cm.gnuplot2_r)
    plt.clim(0, 0.4)
    plt.quiver(lonproj[::d,::d], latproj[::d,::d], tau_x_mth_abs[::d,::d,tt], tau_y_mth_abs[::d,::d,tt], scale=2, linewidth=0.4)
    plt.title(labels_month[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[N m$^{-2}$]')

plt.savefig('mhw_properties/' + mhwname + '_OMAPS_TAU_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

# Wind Stress Anomaly
plt.clf()
d = 5
cnt = 0
for tt in range(24+9-1, 36+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(lon_uv, lat_uv)
    H = plt.contourf(lonproj, latproj, np.sqrt(tau_x_mth[:,:,tt]**2 + tau_y_mth[:,:,tt]**2), levels=np.arange(0.2,1+0.1,0.1)*0.25, cmap=plt.cm.gnuplot2_r)
    plt.clim(0, 0.4)
    plt.quiver(lonproj[::d,::d], latproj[::d,::d], tau_x_mth[::d,::d,tt], tau_y_mth[::d,::d,tt], scale=2, linewidth=0.4)
    plt.title(labels_month[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[N m$^{-2}$]')

plt.savefig('mhw_properties/' + mhwname + '_OMAPS_TAUAnom_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=150)


