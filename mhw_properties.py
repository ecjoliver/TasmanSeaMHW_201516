'''

  Software which uses the MHW definition
  of Hobday et al. (2015) applied to
  select SST time series around the globe

'''

# Load required modules

import numpy as np
from scipy import io
from scipy import ndimage
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
# observations
#

pathroot = '/mnt/erebor/'
pathroot = '/mnt/EREBOR/'
#pathroot = '/home/ecoliver/Desktop/'
#pathroot = '/media/ecoliver/DataOne/'
header = pathroot+'data/sst/noaa_oi_v2/avhrr/'
file0 = header + '1982/avhrr-only-v2.19820101.nc'

#
# lat and lons of obs
#

fileobj = Dataset(file0, 'r')
lon = fileobj.variables['lon'][:].astype(float)
lat = fileobj.variables['lat'][:].astype(float)
fill_value = fileobj.variables['sst']._FillValue.astype(float)
scale = fileobj.variables['sst'].scale_factor.astype(float)
offset = fileobj.variables['sst'].add_offset.astype(float)
fileobj.close()

#
# HadISST
#

file = pathroot + '/data/sst/HadSST/HadISST1/HadISST_sst.nc'
#file = pathroot + '/data/sst/HadSST/HadSST3/HadSST.3.1.1.0.median.nc'
fileobj = Dataset(file, 'r')
time = fileobj.variables['time'][:]
lon_had = fileobj.variables['longitude'][:]
lat_had = fileobj.variables['latitude'][:]
sst_monthly = fileobj.variables['sst'][:].data
fillValue = fileobj.variables['sst']._FillValue
sst_monthly[sst_monthly==fillValue] = np.nan
sst_monthly[sst_monthly<=-2] = np.nan
sst_monthly[sst_monthly>=35] = np.nan
# Make time vector
#t_had, dates_had, T_had, year_had, month_had, day_had, doy_had = ecj.timevector([1870, 1, 1], [2016, 1, 30])
#t_had, dates_had, T_had, year_had, month_had, day_had, doy_had = ecj.timevector([1850, 1, 1], [2016, 2, 28])
t_had = np.floor(date(1870,1,1).toordinal() + time).astype(int) # HadISST
#t_had = np.floor(date(1850,1,1).toordinal() + time).astype(int) # HadSST3
fileobj.close()

#t_had = t_had[day_had==15]
#year_had = year_had[day_had==15]
#month_had = month_had[day_had==15]
#T_had = len(t_had)

T_had = len(t_had)
dates_had = []
year_had = np.zeros(T_had)
month_had = np.zeros(T_had)
for mth in range(T_had):
    dates_had.append(date.fromordinal(t_had[mth].astype(int)))
    year_had[mth] = dates_had[mth].year
    month_had[mth] = dates_had[mth].month

#
# Southern hemisphere DJF mean SST
#

file = '/data/MHWs/Tasmania_2015_2016/NOAAOISST/X131.217.255.76.111.22.53.23.nc'
fileobj = Dataset(file, 'r')
sst_DJF = np.mean(fileobj.variables['anom'][:], axis=0)
lon_DJF = fileobj.variables['lon'][:]
lat_DJF = fileobj.variables['lat'][:]
fileobj.close()

#
# OceanCurrent
#

# Load data
#file = '/data/MHWs/Tasmania_2015_2016/OceanCurrent/IMOS-Aggregation-20160507T151727.036+1000.nc'
file = '/data/MHWs/Tasmania_2015_2016/OceanCurrent/IMOS-Aggregation-20160707T154423.429+1000.nc'
fileobj = Dataset(file, 'r')
time_OC = fileobj.variables['TIME'][:].astype(int) + date(1985,1,1).toordinal()
lon_OC = fileobj.variables['LONGITUDE'][:]
lat_OC = fileobj.variables['LATITUDE'][:]
U0 = fileobj.variables['UCUR'][:].data.astype(float)
V0 = fileobj.variables['VCUR'][:].data.astype(float)
U0[U0==fileobj.variables['UCUR']._FillValue] = np.nan
V0[V0==fileobj.variables['VCUR']._FillValue] = np.nan
#U0 = U0 * fileobj.variables['UCUR'].scale_factor + fileobj.variables['UCUR'].add_offset
#V0 = V0 * fileobj.variables['VCUR'].scale_factor + fileobj.variables['VCUR'].add_offset
fileobj.close()

# Some days are missing, so fill out matrix with NaNs for missing days
llon_OC, llat_OC = np.meshgrid(lon_OC, lat_OC)
date_OC_start = date.fromordinal(time_OC[0])
date_OC_end = date.fromordinal(time_OC[-1])
t_OC, dates_OC, T_OC, year_OC, month_OC, day_OC, doy_OC = ecj.timevector([date_OC_start.year, date_OC_start.month, date_OC_start.day], [date_OC_end.year, date_OC_end.month, date_OC_end.day])

U = np.nan*np.zeros((T_OC, len(lat_OC), len(lon_OC)))
V = np.nan*np.zeros((T_OC, len(lat_OC), len(lon_OC)))
for tt in range(T_OC):
    ttt = time_OC == t_OC[tt]
    if ttt.sum() == 0:
        continue
    U[tt,:,:] = U0[ttt,:,:]
    V[tt,:,:] = V0[ttt,:,:]

# Create monthly means and anomalies
year_start = 2012
year_end = 2016
Tl = (year_end-year_start+1)*12

U_mth = np.NaN*np.zeros((len(lat_OC), len(lon_OC), Tl))
V_mth = np.NaN*np.zeros((len(lat_OC), len(lon_OC), Tl))
U_mth_abs = np.NaN*np.zeros((len(lat_OC), len(lon_OC), Tl))
V_mth_abs = np.NaN*np.zeros((len(lat_OC), len(lon_OC), Tl))
EKE_mth = np.NaN*np.zeros((len(lat_OC), len(lon_OC), Tl))

for ii in range(len(lon_OC)):
    print ii, len(lon_OC)
    for jj in range(len(lat_OC)):
        if (np.isnan(U[:,jj,ii]).sum() == T_OC) + (np.isnan(V[:,jj,ii]).sum() == T_OC):
            continue
        U_ds, U_s, U_beta = ds.deseason_harmonic(U[:,jj,ii], 2, 365.25)
        V_ds, V_s, V_beta = ds.deseason_harmonic(V[:,jj,ii], 2, 365.25)
        EKE = 0.5*(np.array(U_ds)**2 + np.array(V_ds)**2)
        tt = 0
        for yr in range(year_start, year_end+1):
            for mth in range(1,12+1):
                U_mth[jj,ii,tt] = np.nanmean(U_ds[(year_OC==yr) * (month_OC==mth)])
                V_mth[jj,ii,tt] = np.nanmean(V_ds[(year_OC==yr) * (month_OC==mth)])
                U_mth_abs[jj,ii,tt] = np.nanmean(U[(year_OC==yr) * (month_OC==mth),jj,ii])
                V_mth_abs[jj,ii,tt] = np.nanmean(V[(year_OC==yr) * (month_OC==mth),jj,ii])
                EKE_mth[jj,ii,tt] = np.nanmean(EKE[(year_OC==yr) * (month_OC==mth)])
                tt += 1

#
# CFSv2 air temperature
#

# load data
yr0 = 2012
file = '/data/MHWs/Tasmania_2015_2016/CFSv2/CFSv2_' + str(yr0) + '.nc'
fileobj = Dataset(file, 'r')
lon_CFS = fileobj.variables['x'][:]
lat_CFS = fileobj.variables['y'][:]
t_CFS = fileobj.variables['t'][:]
air_temp = fileobj.variables['air_temp'][:]
u = fileobj.variables['u'][:]
v = fileobj.variables['v'][:]
fileobj.close()

for yr in range(2013, 2016+1):
    file = '/data/MHWs/Tasmania_2015_2016/CFSv2/CFSv2_' + str(yr) + '.nc'
    fileobj = Dataset(file, 'r')
    t_CFS = np.append(t_CFS, fileobj.variables['t'][:])
    air_temp = np.append(air_temp, fileobj.variables['air_temp'][:], axis=0)
    u = np.append(u, fileobj.variables['u'][:], axis=0)
    v = np.append(v, fileobj.variables['v'][:], axis=0)
    fileobj.close()

# Build up time and date vectors
t_CFS = np.floor(t_CFS + date(1990,1,1).toordinal()).astype(int)
year_CFS = []
month_CFS = []
for tt in range(len(t_CFS)):
    year_CFS.append(date.fromordinal(t_CFS[tt]).year)
    month_CFS.append(date.fromordinal(t_CFS[tt]).month)

year_CFS = np.array(year_CFS)
month_CFS = np.array(month_CFS)
T_CFS = len(t_CFS)

# Create monthly means and anomalies
year_start = 2012
year_end = 2016

air_temp_mth = np.NaN*np.zeros((lat_CFS.shape[0], lat_CFS.shape[1], Tl))
air_temp_mth_abs = np.NaN*np.zeros((lat_CFS.shape[0], lat_CFS.shape[1], Tl))
u_mth = np.NaN*np.zeros((lat_CFS.shape[0], lat_CFS.shape[1], Tl))
u_mth_abs = np.NaN*np.zeros((lat_CFS.shape[0], lat_CFS.shape[1], Tl))
v_mth = np.NaN*np.zeros((lat_CFS.shape[0], lat_CFS.shape[1], Tl))
v_mth_abs = np.NaN*np.zeros((lat_CFS.shape[0], lat_CFS.shape[1], Tl))

for ii in range(len(lon_CFS)):
    print ii, len(lon_CFS)
    for jj in range(len(lat_CFS)):
        if np.isnan(air_temp[:,jj,ii]).sum() == T_CFS:
            continue
        air_temp_ds, air_temp_s, air_temp_beta = ds.deseason_harmonic(air_temp[:,jj,ii], 2, 1460)
        u_ds, u_s, u_beta = ds.deseason_harmonic(u[:,jj,ii], 2, 1460)
        v_ds, v_s, v_beta = ds.deseason_harmonic(v[:,jj,ii], 2, 1460)
        tt = 0
        for yr in range(year_start, year_end+1):
            for mth in range(1,12+1):
                air_temp_mth[jj,ii,tt] = np.nanmean(air_temp_ds[(year_CFS==yr) * (month_CFS==mth)])
                air_temp_mth_abs[jj,ii,tt] = np.nanmean(air_temp[(year_CFS==yr) * (month_CFS==mth),jj,ii])
                u_mth[jj,ii,tt] = np.nanmean(u_ds[(year_CFS==yr) * (month_CFS==mth)])
                u_mth_abs[jj,ii,tt] = np.nanmean(u[(year_CFS==yr) * (month_CFS==mth),jj,ii])
                v_mth[jj,ii,tt] = np.nanmean(v_ds[(year_CFS==yr) * (month_CFS==mth)])
                v_mth_abs[jj,ii,tt] = np.nanmean(v[(year_CFS==yr) * (month_CFS==mth),jj,ii])
                tt += 1

#
# OceanMAPS circulation
#

# load data
file = '/data/MHWs/Tasmania_2015_2016/OMAPS/OMAPS_2012_08_uv.nc'
fileobj = Dataset(file, 'r')
lon_OMAPS = fileobj.variables['x'][:]
lat_OMAPS = fileobj.variables['y'][:]
depth_OMAPS = fileobj.variables['z'][:]
fileobj.close()

# Create monthly means and anomalies
year_start = 2012
year_end = 2016
Tl = (year_end-year_start+1)*12
kH = -14 # -14 -> 95m depth

#u_mth = np.NaN*np.zeros((lat_OMAPS.shape[0], lat_OMAPS.shape[1], Tl))
u_mth_abs = np.NaN*np.zeros((lat_OMAPS.shape[0], lat_OMAPS.shape[1], Tl))
#v_mth = np.NaN*np.zeros((lat_OMAPS.shape[0], lat_OMAPS.shape[1], Tl))
v_mth_abs = np.NaN*np.zeros((lat_OMAPS.shape[0], lat_OMAPS.shape[1], Tl))

yr = 2012
mth = 1
for tt in range(Tl):
    print tt+1, Tl
    if (yr == 2012) * (mth < 8):
        # Up counter
        mth += 1
        if mth == 13:
            yr += 1
            mth = 1
        continue
    # load data, monthly average
    # UV
    file = '/data/MHWs/Tasmania_2015_2016/OMAPS/OMAPS_' + str(yr) + '_' + str(mth).zfill(2) + '_uv.nc'
    fileobj = Dataset(file, 'r')
    u_mth_abs[:,:,tt] = np.nanmean(np.mean(fileobj.variables['u'][:,kH:,:,:], axis=0), axis=0)
    v_mth_abs[:,:,tt] = np.nanmean(np.mean(fileobj.variables['v'][:,kH:,:,:], axis=0), axis=0)
    fileobj.close()
    # Up counter
    mth += 1
    if mth == 13:
        yr += 1
        mth = 1



#
# Load data at locations of interest
#

locations = {}
locations['lon'] = [147, 157] #[147, 155]
locations['lat'] = [-46, -39] #[-45, -37]
locations['lon_map'] = [np.mean(locations['lon'])]
locations['lat_map'] = [np.mean(locations['lat'])]
locations['name'] = ['SE Australia']
print "Area of box:", ecj.latlonArea(locations['lon'][0], locations['lat'][0], locations['lon'][1], locations['lat'][1]), "km^2"

# Time and date vectors
t, dates, T, year, month, day, doy = ecj.timevector([1982,1,1], [2016,7,5])

# Load data
sst = np.zeros((T))
i1 = np.where(lon > locations['lon'][0])[0][0] - 1
i2 = np.where(lon > locations['lon'][1])[0][0]
j1 = np.where(lat > locations['lat'][0])[0][0] - 1
j2 = np.where(lat > locations['lat'][1])[0][0]
for tt in range(T):
    try:
        file = header + str(dates[tt].year) + '/avhrr-only-v2.' + str(dates[tt].year) + str(dates[tt].month).zfill(2) + str(dates[tt].day).zfill(2) + '.nc'
        fileobj = Dataset(file, 'r')
    except:
        file = header + str(dates[tt].year) + '/avhrr-only-v2.' + str(dates[tt].year) + str(dates[tt].month).zfill(2) + str(dates[tt].day).zfill(2) + '_preliminary.nc'
        fileobj = Dataset(file, 'r')
    print str(dates[tt].year) + str(dates[tt].month).zfill(2) + str(dates[tt].day).zfill(2)
    sst0 = fileobj.variables['sst'][0,0,j1:j2+1,i1:i2+1].astype(float).data
    sst0[sst0==fill_value] = np.nan
    sst[tt] = np.nanmean(sst0) #*scale + offset)
    fileobj.close()

# HadISST
sst_had = np.zeros((T_had))
ii = (lon_had >= locations['lon'][0]) * (lon_had <= locations['lon'][1])
jj = (lat_had >= locations['lat'][0]) * (lat_had <= locations['lat'][1])
sst_had = np.nanmean(np.nanmean(sst_monthly[:,:,ii][:,jj,:], axis=2), axis=1)
missing_had = 1. - 1.*np.sum(np.sum(np.isnan(sst_monthly[:,:,ii][:,jj,:]), axis=2), axis=1) / (ii.sum() + jj.sum())
# Climatology
sst_had_clim = np.zeros((T_had))
for tt in range(len(t_had)):
    sst_had_clim[tt] = np.nanmean(sst_had[month_had==month_had[tt]])

# HadSST3 Check some counds of valid months
# 1 - 1.*np.isnan(sst_had[(year_had>=1900)*(year_had<=2016)]).sum()/len(sst_had[(year_had>=1900)*(year_had<=2016)])
# 1 - 1.*np.isnan(sst_had[(year_had>=1900)*(year_had<=1950)]).sum()/len(sst_had[(year_had>=1900)*(year_had<=1950)])
# 1 - 1.*np.isnan(sst_had[(year_had>=1951)*(year_had<=2016)]).sum()/len(sst_had[(year_had>=1951)*(year_had<=2016)])
# 1 - 1.*np.isnan(sst_had[(year_had>=1911)*(year_had<=1940)]).sum()/len(sst_had[(year_had>=1911)*(year_had<=1940)])
# 1 - 1.*np.isnan(sst_had[(year_had>=1881)*(year_had<=1910)]).sum()/len(sst_had[(year_had>=1881)*(year_had<=1910)])
plt.figure()
plt.clf()
plt.plot(np.arange(1855, 2005+1, 10), 100*np.sum(~np.isnan(np.reshape(sst_had[:120*np.floor(len(sst_had)/120.).astype(int)], (16, 120))), axis=1)/120., 'k-o', linewidth=2)
plt.ylim(0, 105)
plt.xlim(1850, 2020)
plt.ylabel('%')
plt.title('Proportion of valid months per decade (HadSST3)')
#plt.savefig('../../documents/14_Tasmania_2015_2016/figures/PropValidMonths_HadSST3.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Mask for low obs count
mask_had = np.ones(sst_had.shape)
mask_had[missing_had<=0.75] = np.nan
sst_had = sst_had*mask_had

#
# Save data
#

outfile = '/data/MHWs/Tasmania_2015_2016/NOAAOISST/' + str(locations['lon'][0]) + '_' + str(locations['lon'][1]) + '_' + str(locations['lat'][0]) + '_' + str(locations['lat'][1]) + '/mhw_SEAus_NOAAOISST_HadSST'

#np.savez(outfile, lon=locations['lon'], lat=locations['lat'], sst=sst, t=t, dates=dates, T=T, year=year, month=month, day=day, doy=doy, sst_had=sst_had, sst_had_clim=sst_had_clim, t_had=t_had, dates_had=dates_had, T_had=T_had, year_had=year_had, month_had=month_had)

# Load data
data = np.load(outfile + '.npz')
sst = data['sst']
dates_had = data['dates_had']
sst_had = data['sst_had']
sst_had_clim = data['sst_had_clim']

#
# Apply Marine Heat Wave definition
#

#mhws, clim = mhw.detect(t, sst, coldSpells=coldSpells)
mhws, clim = mhw.detect(t, sst, coldSpells=coldSpells, climatologyPeriod=[1982,2005])
mhwBlock = mhw.blockAverage(t, mhws, temp=sst)
mean, trend, dtrend = mhw.meanTrend(mhwBlock)

# Print event properties (1982-2005 climatology)
evMax = np.argmax(mhws['duration'])
print 'Maximum intensity:', mhws['intensity_max'][evMax], 'deg. C'
print 'Average intensity:', mhws['intensity_mean'][evMax], 'deg. C'
print 'Cumulative intensity:', mhws['intensity_cumulative'][evMax], 'deg. C-days'
print 'Duration:', mhws['duration'][evMax], 'days'
print 'Start date:', mhws['date_start'][evMax].strftime("%d %B %Y")
print 'End date:', mhws['date_end'][evMax].strftime("%d %B %Y")

# ALTERNATE SHIFT MEAN BASED ON WARMING SINCE EARLY PERIOD
dt = np.nanmean(sst_had[(year_had>=1982)*(year_had<=2005)]) - np.nanmean(sst_had[(year_had>=1911)*(year_had<=1940)])
dt = np.nanmean(sst_had[(year_had>=1982)*(year_had<=2005)]) - np.nanmean(sst_had[(year_had>=1881)*(year_had<=1910)])
mhws_early, clim_early = mhw.detect(t, sst+dt, coldSpells=coldSpells, climatologyPeriod=[1982,2005], alternateClimatology=[t, sst])

# Print event properties (earlier climatology)
evMax = np.argmax(mhws_early['duration'])
print 'Maximum intensity:', mhws_early['intensity_max'][evMax], 'deg. C'
print 'Average intensity:', mhws_early['intensity_mean'][evMax], 'deg. C'
print 'Cumulative intensity:', mhws_early['intensity_cumulative'][evMax], 'deg. C-days'
print 'Duration:', mhws_early['duration'][evMax], 'days'
print 'Start date:', mhws_early['date_start'][evMax].strftime("%d %B %Y")
print 'End date:', mhws_early['date_end'][evMax].strftime("%d %B %Y")

# Plot various summary things

plt.figure(figsize=(15,7))
plt.subplot(2,2,1)
evMax = np.argmax(mhws['duration'])
plt.bar(range(mhws['n_events']), mhws['duration'], width=0.6, color=(0.7,0.7,0.7))
plt.bar(evMax, mhws['duration'][evMax], width=0.6, color=col_bar)
plt.xlim(0, mhws['n_events'])
plt.ylabel('[days]')
plt.title('Duration')
plt.subplot(2,2,2)
evMax = np.argmax(np.abs(mhws['intensity_max']))
plt.bar(range(mhws['n_events']), mhws['intensity_max'], width=0.6, color=(0.7,0.7,0.7))
plt.bar(evMax, mhws['intensity_max'][evMax], width=0.6, color=col_bar)
plt.xlim(0, mhws['n_events'])
plt.ylabel(r'[$^\circ$C]')
plt.title('Maximum Intensity')
plt.subplot(2,2,4)
evMax = np.argmax(np.abs(mhws['intensity_mean']))
plt.bar(range(mhws['n_events']), mhws['intensity_mean'], width=0.6, color=(0.7,0.7,0.7))
plt.bar(evMax, mhws['intensity_mean'][evMax], width=0.6, color=col_bar)
plt.xlim(0, mhws['n_events'])
plt.title('Mean Intensity')
plt.ylabel(r'[$^\circ$C]')
plt.xlabel(mhwname + ' event number')
plt.subplot(2,2,3)
evMax = np.argmax(np.abs(mhws['intensity_cumulative']))
plt.bar(range(mhws['n_events']), mhws['intensity_cumulative'], width=0.6, color=(0.7,0.7,0.7))
plt.bar(evMax, mhws['intensity_cumulative'][evMax], width=0.6, color=col_bar)
plt.xlim(0, mhws['n_events'])
plt.title(r'Cumulative Intensity')
plt.ylabel(r'[$^\circ$C$\times$days]')
plt.xlabel(mhwname + ' event number')
plt.savefig('mhw_properties/' + mhwname + '_list_byNumber.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

ts = date(1982,1,1).toordinal()
te = date(2016,8,1).toordinal()
plt.figure(figsize=(15,7))
plt.subplot(2,2,1)
evMax = np.argmax(mhws['duration'])
plt.bar(mhws['date_peak'], mhws['duration'], width=150, color=(0.7,0.7,0.7))
plt.bar(mhws['date_peak'][evMax], mhws['duration'][evMax], width=150, color=col_bar)
plt.xlim(ts, te)
plt.ylabel('[days]')
plt.title('(A) Duration')
plt.subplot(2,2,2)
evMax = np.argmax(np.abs(mhws['intensity_max']))
plt.bar(mhws['date_peak'], mhws['intensity_max'], width=150, color=(0.7,0.7,0.7))
plt.bar(mhws['date_peak'][evMax], mhws['intensity_max'][evMax], width=150, color=col_bar)
plt.xlim(ts, te)
plt.ylabel(r'[$^\circ$C]')
plt.title('(B) Maximum Intensity')
plt.subplot(2,2,4)
evMax = np.argmax(np.abs(mhws['intensity_mean']))
plt.bar(mhws['date_peak'], mhws['intensity_mean'], width=150, color=(0.7,0.7,0.7))
plt.bar(mhws['date_peak'][evMax], mhws['intensity_mean'][evMax], width=150, color=col_bar)
plt.xlim(ts, te)
plt.title('(C) Mean Intensity')
plt.ylabel(r'[$^\circ$C]')
plt.subplot(2,2,3)
evMax = np.argmax(np.abs(mhws['intensity_cumulative']))
plt.bar(mhws['date_peak'], mhws['intensity_cumulative'], width=150, color=(0.7,0.7,0.7))
plt.bar(mhws['date_peak'][evMax], mhws['intensity_cumulative'][evMax], width=150, color=col_bar)
plt.xlim(ts, te)
plt.title(r'(D) Cumulative Intensity')
plt.ylabel(r'[$^\circ$C$\times$days]')
plt.savefig('../../documents/14_Tasmania_2015_2016/figures/' + mhwname + '_list_byDate.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

plt.clf()
plt.subplot(2,2,1)
evMax = np.argmax(mhws['duration'])
plt.bar(mhws['date_peak'], mhws['duration'], width=150, color=(0.7,0.7,0.7))
plt.bar(mhws['date_peak'][evMax], mhws['duration'][evMax], width=150, color=col_bar)
plt.xlim(ts, te)
plt.ylabel('[days]')
plt.title('(D) Duration')
plt.subplot(2,2,2)
evMax = np.argmax(np.abs(mhws['intensity_max']))
plt.bar(mhws['date_peak'], mhws['intensity_max'], width=150, color=(0.7,0.7,0.7))
plt.bar(mhws['date_peak'][evMax], mhws['intensity_max'][evMax], width=150, color=col_bar)
plt.xlim(ts, te)
plt.ylabel(r'[$^\circ$C]')
plt.title('(E) Maximum Intensity')
plt.savefig('../../documents/14_Tasmania_2015_2016/figures/' + mhwname + '_list_byDate_Dur_iMax.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Plot top 10 events

# Maximum intensity
outfile = open('mhw_properties/' + mhwname + '_topTen_iMax.txt', 'w')
evs = np.argsort(np.abs(mhws['intensity_max']))[-10:]
plt.figure(figsize=(23,16))
for i in range(10):
    ev = evs[-(i+1)]
    plt.subplot(5,2,i+1)
    # Find indices for all ten MHWs before and after event of interest and shade accordingly
    for ev0 in np.arange(max(ev-10,0), min(ev+11,mhws['n_events']-1), 1):
        t1 = np.where(t==mhws['time_start'][ev0])[0][0]
        t2 = np.where(t==mhws['time_end'][ev0])[0][0]
        plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_ev)
    # Find indices for MHW of interest (2011 WA event) and shade accordingly
    t1 = np.where(t==mhws['time_start'][ev])[0][0]
    t2 = np.where(t==mhws['time_end'][ev])[0][0]
    plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_evMax)
    # Plot SST, seasonal cycle, threshold, shade MHWs with main event in red
    plt.plot(dates, sst, 'k-', linewidth=2)
    plt.plot(dates, clim['thresh'], col_thresh, linewidth=2)
    plt.plot(dates, clim['seas'], col_clim, linewidth=2)
    plt.title('Number ' + str(i+1))
    plt.xlim(mhws['time_start'][ev]-150, mhws['time_end'][ev]+150)
    if coldSpells:
        plt.ylim(clim['seas'].min() + mhws['intensity_max'][ev] - 0.5, clim['seas'].max() + 1)
    else:
        plt.ylim(clim['seas'].min() - 1, clim['seas'].max() + mhws['intensity_max'][ev] + 0.5)
    plt.ylabel(r'SST [$^\circ$C]')
    # Save stats
    outfile.write('Number ' + str(i+1) + '\n')
    outfile.write('Maximum intensity: ' + str(mhws['intensity_max'][ev]) + ' deg. C\n')
    outfile.write('Average intensity: '+ str( mhws['intensity_mean'][ev]) + ' deg. C\n')
    outfile.write('Cumulative intensity: ' + str(mhws['intensity_cumulative'][ev]) + ' deg. C-days\n')
    outfile.write('Duration: ' + str(mhws['duration'][ev]) + ' days\n')
    outfile.write('Start date: ' + str(mhws['date_start'][ev].strftime("%d %B %Y")) + '\n')
    outfile.write('End date: ' + str(mhws['date_end'][ev].strftime("%d %B %Y")) + '\n')
    outfile.write('\n')

plt.legend(['SST', 'threshold', 'seasonal climatology'], loc=4)
outfile.close()
plt.savefig('mhw_properties/' + mhwname + '_topTen_iMax.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Mean intensity
outfile = open('mhw_properties/' + mhwname + '_topTen_iMean.txt', 'w')
evs = np.argsort(np.abs(mhws['intensity_mean']))[-10:]
plt.clf()
for i in range(10):
    ev = evs[-(i+1)]
    plt.subplot(5,2,i+1)
    # Find indices for all ten MHWs before and after event of interest and shade accordingly
    for ev0 in np.arange(max(ev-10,0), min(ev+11,mhws['n_events']-1), 1):
        t1 = np.where(t==mhws['time_start'][ev0])[0][0]
        t2 = np.where(t==mhws['time_end'][ev0])[0][0]
        plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_ev)
    # Find indices for MHW of interest (2011 WA event) and shade accordingly
    t1 = np.where(t==mhws['time_start'][ev])[0][0]
    t2 = np.where(t==mhws['time_end'][ev])[0][0]
    plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_evMax)
    # Plot SST, seasonal cycle, threshold, shade MHWs with main event in red
    plt.plot(dates, sst, 'k-', linewidth=2)
    plt.plot(dates, clim['thresh'], col_thresh, linewidth=2)
    plt.plot(dates, clim['seas'], col_clim, linewidth=2)
    plt.title('Number ' + str(i+1))
    plt.xlim(mhws['time_start'][ev]-150, mhws['time_end'][ev]+150)
    if coldSpells:
        plt.ylim(clim['seas'].min() + mhws['intensity_max'][ev] - 0.5, clim['seas'].max() + 1)
    else:
        plt.ylim(clim['seas'].min() - 1, clim['seas'].max() + mhws['intensity_max'][ev] + 0.5)
    plt.ylabel(r'SST [$^\circ$C]')
    # Save stats
    outfile.write('Number ' + str(i+1) + '\n')
    outfile.write('Maximum intensity: ' + str(mhws['intensity_max'][ev]) + ' deg. C\n')
    outfile.write('Average intensity: '+ str( mhws['intensity_mean'][ev]) + ' deg. C\n')
    outfile.write('Cumulative intensity: ' + str(mhws['intensity_cumulative'][ev]) + ' deg. C-days\n')
    outfile.write('Duration: ' + str(mhws['duration'][ev]) + ' days\n')
    outfile.write('Start date: ' + str(mhws['date_start'][ev].strftime("%d %B %Y")) + '\n')
    outfile.write('End date: ' + str(mhws['date_end'][ev].strftime("%d %B %Y")) + '\n')
    outfile.write('\n')

plt.legend(['SST', 'threshold', 'seasonal climatology'], loc=4)
outfile.close()
plt.savefig('mhw_properties/' + mhwname + '_topTen_iMean.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Cumulative intensity
outfile = open('mhw_properties/' + mhwname + '_topTen_iCum.txt', 'w')
evs = np.argsort(np.abs(mhws['intensity_cumulative']))[-10:]
plt.clf()
for i in range(10):
    ev = evs[-(i+1)]
    plt.subplot(5,2,i+1)
    # Find indices for all ten MHWs before and after event of interest and shade accordingly
    for ev0 in np.arange(max(ev-10,0), min(ev+11,mhws['n_events']-1), 1):
        t1 = np.where(t==mhws['time_start'][ev0])[0][0]
        t2 = np.where(t==mhws['time_end'][ev0])[0][0]
        plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_ev)
    # Find indices for MHW of interest (2011 WA event) and shade accordingly
    t1 = np.where(t==mhws['time_start'][ev])[0][0]
    t2 = np.where(t==mhws['time_end'][ev])[0][0]
    plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_evMax)
    # Plot SST, seasonal cycle, threshold, shade MHWs with main event in red
    plt.plot(dates, sst, 'k-', linewidth=2)
    plt.plot(dates, clim['thresh'], col_thresh, linewidth=2)
    plt.plot(dates, clim['seas'], col_clim, linewidth=2)
    plt.title('Number ' + str(i+1))
    plt.xlim(mhws['time_start'][ev]-150, mhws['time_end'][ev]+150)
    if coldSpells:
        plt.ylim(clim['seas'].min() + mhws['intensity_max'][ev] - 0.5, clim['seas'].max() + 1)
    else:
        plt.ylim(clim['seas'].min() - 1, clim['seas'].max() + mhws['intensity_max'][ev] + 0.5)
    plt.ylabel(r'SST [$^\circ$C]')
    # Save stats
    outfile.write('Number ' + str(i+1) + '\n')
    outfile.write('Maximum intensity: ' + str(mhws['intensity_max'][ev]) + ' deg. C\n')
    outfile.write('Average intensity: '+ str( mhws['intensity_mean'][ev]) + ' deg. C\n')
    outfile.write('Cumulative intensity: ' + str(mhws['intensity_cumulative'][ev]) + ' deg. C-days\n')
    outfile.write('Duration: ' + str(mhws['duration'][ev]) + ' days\n')
    outfile.write('Start date: ' + str(mhws['date_start'][ev].strftime("%d %B %Y")) + '\n')
    outfile.write('End date: ' + str(mhws['date_end'][ev].strftime("%d %B %Y")) + '\n')
    outfile.write('\n')

plt.legend(['SST', 'threshold', 'seasonal climatology'], loc=4)
outfile.close()
plt.savefig('mhw_properties/' + mhwname + '_topTen_iCum.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Duration
outfile = open('mhw_properties/' + mhwname + '_topTen_Dur.txt', 'w')
evs = np.argsort(mhws['duration'])[-10:]
plt.clf()
for i in range(10):
    ev = evs[-(i+1)]
    plt.subplot(5,2,i+1)
    # Find indices for all ten MHWs before and after event of interest and shade accordingly
    for ev0 in np.arange(max(ev-10,0), min(ev+11,mhws['n_events']-1), 1):
        t1 = np.where(t==mhws['time_start'][ev0])[0][0]
        t2 = np.where(t==mhws['time_end'][ev0])[0][0]
        plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_ev)
    # Find indices for MHW of interest (2011 WA event) and shade accordingly
    t1 = np.where(t==mhws['time_start'][ev])[0][0]
    t2 = np.where(t==mhws['time_end'][ev])[0][0]
    plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_evMax)
    # Plot SST, seasonal cycle, threshold, shade MHWs with main event in red
    plt.plot(dates, sst, 'k-', linewidth=2)
    plt.plot(dates, clim['thresh'], col_thresh, linewidth=2)
    plt.plot(dates, clim['seas'], col_clim, linewidth=2)
    plt.title('Number ' + str(i+1))
    plt.xlim(mhws['time_start'][ev]-150, mhws['time_end'][ev]+150)
    if coldSpells:
        plt.ylim(clim['seas'].min() + mhws['intensity_max'][ev] - 0.5, clim['seas'].max() + 1)
    else:
        plt.ylim(clim['seas'].min() - 1, clim['seas'].max() + mhws['intensity_max'][ev] + 0.5)
    plt.ylabel(r'SST [$^\circ$C]')
    # Save stats
    outfile.write('Number ' + str(i+1) + '\n')
    outfile.write('Maximum intensity: ' + str(mhws['intensity_max'][ev]) + ' deg. C\n')
    outfile.write('Average intensity: '+ str( mhws['intensity_mean'][ev]) + ' deg. C\n')
    outfile.write('Cumulative intensity: ' + str(mhws['intensity_cumulative'][ev]) + ' deg. C-days\n')
    outfile.write('Duration: ' + str(mhws['duration'][ev]) + ' days\n')
    outfile.write('Start date: ' + str(mhws['date_start'][ev].strftime("%d %B %Y")) + '\n')
    outfile.write('End date: ' + str(mhws['date_end'][ev].strftime("%d %B %Y")) + '\n')
    outfile.write('\n')

plt.legend(['SST', 'threshold', 'seasonal climatology'], loc=4)
outfile.close()
plt.savefig('mhw_properties/' + mhwname + '_topTen_Dur.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Annual averages
years = mhwBlock['years_centre']
plt.figure(figsize=(13,7))
plt.subplot(2,2,2)
plt.plot(years, mhwBlock['count'], 'k-')
plt.plot(years, mhwBlock['count'], 'ko')
if np.abs(trend['count']) - dtrend['count'] > 0:
     plt.title('Frequency (trend = ' + '{:.2}'.format(10*trend['count']) + '* per decade)')
else:
     plt.title('Frequency (trend = ' + '{:.2}'.format(10*trend['count']) + ' per decade)')
plt.ylabel('[count per year]')
plt.grid()
plt.xlim(years.min()-1, years.max()+1)
plt.subplot(2,2,1)
plt.plot(years, mhwBlock['duration'], 'k-')
plt.plot(years, mhwBlock['duration'], 'ko')
if np.abs(trend['duration']) - dtrend['duration'] > 0:
    plt.title('Duration (trend = ' + '{:.2}'.format(10*trend['duration']) + '* per decade)')
else:
    plt.title('Duration (trend = ' + '{:.2}'.format(10*trend['duration']) + ' per decade)')
plt.ylabel('[days]')
plt.grid()
plt.xlim(years.min()-1, years.max()+1)
plt.subplot(2,2,4)
plt.plot(years, mhwBlock['intensity_max'], '-', color=col_evMax)
plt.plot(years, mhwBlock['intensity_mean'], 'k-')
plt.plot(years, mhwBlock['intensity_max'], 'o', color=col_evMax)
plt.plot(years, mhwBlock['intensity_mean'], 'ko')
plt.legend(['Max', 'mean'], loc=2)
if (np.abs(trend['intensity_max']) - dtrend['intensity_max'] > 0) * (np.abs(trend['intensity_mean']) - dtrend['intensity_mean'] > 0):
    plt.title('Intensity (trend = ' + '{:.2}'.format(10*trend['intensity_max']) + '* (max), ' + '{:.2}'.format(10*trend['intensity_mean'])  + '* (mean) per decade)')
elif (np.abs(trend['intensity_max']) - dtrend['intensity_max'] > 0):
    plt.title('Intensity (trend = ' + '{:.2}'.format(10*trend['intensity_max']) + '* (max), ' + '{:.2}'.format(10*trend['intensity_mean'])  + ' (mean) per decade)')
elif (np.abs(trend['intensity_mean']) - dtrend['intensity_mean'] > 0):
    plt.title('Intensity (trend = ' + '{:.2}'.format(10*trend['intensity_max']) + ' (max), ' + '{:.2}'.format(10*trend['intensity_mean'])  + '* (mean) per decade)')
else:
    plt.title('Intensity (trend = ' + '{:.2}'.format(10*trend['intensity_max']) + ' (max), ' + '{:.2}'.format(10*trend['intensity_mean'])  + ' (mean) per decade)')
plt.ylabel(r'[$^\circ$C]')
plt.grid()
plt.xlim(years.min()-1, years.max()+1)
plt.subplot(2,2,3)
plt.plot(years, mhwBlock['intensity_cumulative'], 'k-')
plt.plot(years, mhwBlock['intensity_cumulative'], 'ko')
if np.abs(trend['intensity_cumulative']) - dtrend['intensity_cumulative'] > 0:
    plt.title('Cumulative intensity (trend = ' + '{:.2}'.format(10*trend['intensity_cumulative']) + '* per decade)')
else:
    plt.title('Cumulative intensity (trend = ' + '{:.2}'.format(10*trend['intensity_cumulative']) + ' per decade)')
plt.ylabel(r'[$^\circ$C$\times$days]')
plt.grid()
plt.xlim(years.min()-1, years.max()+1)
plt.savefig('mhw_properties/' + mhwname + '_annualAverages_meanTrend.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Save results as text data
outfile = 'mhw_properties/' + mhwname + '_data'

# Event data
csvfile = open(outfile +'.events.csv', 'w')
#csvfile.write('# Marine heatwave statistics for individual detected events at [' + str(lon[i]) + ' E ' + str(lat[j]) + ' N] from NOAA OI AVHRR V2 SST data (1982-2014)\n')
csvfile.write('# Marine ' + mhwfullname + ' statistics for individual detected events at ' + locations['name'][0] + ' from NOAA OI AVHRR V2 SST data (1982-2014)\n')
csvfile.write('Event number, Start year, Start month, Start day, Peak year, Peak month, Peak day, End year, End month, End day, Duration [days], Maximum intensity [deg C], Mean intensity [deg C], Cumulative intensity [deg C x days], Intensity variability [deg C], Rate of onset [deg C / days], Rate of decline [deg C / days], Maximum intensity (rel. thresh.) [deg C], Mean intensity (rel. thresh.) [deg C], Cumulative intensity (rel. thresh.) [deg C x days], Intensity variability (rel. thresh.) [deg C], Maximum intensity (absolute) [deg C], Mean intensity (absolute) [deg C], Cumulative intensity (absolute) [deg C x days], Intensity variability (absolute) [deg C], Maximum intensity (normalized) [unitless], Mean intensity (normalized) [unitless]\n')
for ev in range(mhws['n_events']):
    csvfile.write(str(ev+1) + ', ' + str(mhws['date_start'][ev].year) + ', ' + str(mhws['date_start'][ev].month) + ', ' + str(mhws['date_start'][ev].day) + ', ' + str(mhws['date_peak'][ev].year) + ', ' + str(mhws['date_peak'][ev].month) + ', ' + str(mhws['date_peak'][ev].day) + ', ' + str(mhws['date_end'][ev].year) + ', ' + str(mhws['date_end'][ev].month) + ', ' + str(mhws['date_end'][ev].day) + ', ' + str(mhws['duration'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_max'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_mean'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_cumulative'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_var'][ev]) + ', ' + '{:.4}'.format(mhws['rate_onset'][ev]) + ', ' + '{:.4}'.format(mhws['rate_decline'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_max_relThresh'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_mean_relThresh'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_cumulative_relThresh'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_var_relThresh'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_max_abs'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_mean_abs'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_cumulative_abs'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_var_abs'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_max_norm'][ev]) + ', ' + '{:.4}'.format(mhws['intensity_mean_norm'][ev]) + '\n')
csvfile.close()

# Annual average data
csvfile = open(outfile +'.annual.csv', 'w')
#csvfile.write('# Annual average marine heatwave statistics at [' + str(lon[i]) + ' E ' + str(lat[j]) + ' N] from NOAA OI AVHRR V2 SST data (1982-2014)\n')
csvfile.write('# Annual average marine ' + mhwfullname + ' statistics at ' + locations['name'][0] + ' from NOAA OI AVHRR V2 SST data (1982-2014)\n')
csvfile.write('# A value of nan indicates missing data. This should correspond to a year with no MHW events (count = 0)\n')
csvfile.write('Year, ' + mhwname + ' event count [number], Duration [days], Maximum intensity [deg C], Mean intensity [deg C], Cumulative intensity [deg C x days], Total ' + mhwname + ' days [days], Total cumulative intensity [deg C x days], Intensity variability [deg C], Rate of onset [deg C / days], Rate of decline [deg C / days], Maximum intensity (rel. thresh.) [deg C], Mean intensity (rel. thresh.) [deg C], Cumulative intensity (rel. thresh.) [deg C x days], Intensity variability (rel. thresh.) [deg C], Maximum intensity (absolute) [deg C], Mean intensity (absolute) [deg C], Cumulative intensity (absolute) [deg C x days], Intensity variability (absolute) [deg C], Maximum intensity (normalized) [unitless], Mean intensity (normalized) [unitless], Mean temperature [deg C], Max temperature [deg C], Min temperature [deg C]\n')
for yr in range(len(mhwBlock['years_centre'])):
    csvfile.write(str(mhwBlock['years_centre'][yr].astype(int)) + ', ' + str(mhwBlock['count'][yr].astype(int)) + ', ' + '{:.4}'.format(mhwBlock['duration'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_max'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_mean'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_cumulative'][yr]) + ', ' + '{:.4}'.format(mhwBlock['total_days'][yr]) + ', ' + '{:.4}'.format(mhwBlock['total_icum'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_var'][yr]) + ', ' + '{:.4}'.format(mhwBlock['rate_onset'][yr]) + ', ' + '{:.4}'.format(mhwBlock['rate_decline'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_max_relThresh'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_mean_relThresh'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_cumulative_relThresh'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_var_relThresh'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_max_abs'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_mean_abs'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_cumulative_abs'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_var_abs'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_max_norm'][yr]) + ', ' + '{:.4}'.format(mhwBlock['intensity_mean_norm'][yr]) + ', ' + '{:.4}'.format(mhwBlock['temp_mean'][yr]) + ', ' + '{:.4}'.format(mhwBlock['temp_max'][yr]) + ', ' + '{:.4}'.format(mhwBlock['temp_min'][yr]) + '\n')
csvfile.close()

# Climatology
csvfile = open(outfile +'.climatology.csv', 'w')
#csvfile.write('# SST climatology at [' + str(lon[i]) + ' E ' + str(lat[j]) + ' N] from NOAA OI AVHRR V2 SST data (1982-2014)\n')
csvfile.write('# SST climatology at ' + locations['name'][0] + ' from NOAA OI AVHRR V2 SST data (1982-2014)\n')
csvfile.write('# Climatology includes the seasonal cycle and the 90th percentile threshold\n')
csvfile.write('Day, Seasonal cycle [deg C], Threshold [deg C]\n')
for tt in range(len(clim['seas'])):
    csvfile.write(str(doy[tt]) + ', ' + '{:.4}'.format(clim['seas'][tt]) + ', ' + '{:.4}'.format(clim['thresh'][tt]) + '\n')
csvfile.close()

#
# Plot some timeseries
#

plt.figure(figsize=(16,8))

ts = date(1880,1,1).toordinal()
te = date(2016,8,1).toordinal()

plt.clf()
ax = plt.subplot(2,1,1)
plt.plot(dates_had, sst_had, 'k-')
plt.plot(dates, sst, 'b-')
ax.set_xticks([date(yr,1,1).toordinal() for yr in range(1880,2050,10)])
plt.xlim(ts, te)
#plt.grid()
plt.ylabel('SST [$^\circ$C]')
plt.title('(A) Full period (1900-2015)')
ax = plt.subplot(2,1,2)
plt.plot(dates, sst - clim['seas'], 'b-')
plt.plot(dates_had, sst_had - sst_had_clim, 'k-')
ax.set_xticks([date(yr,1,1).toordinal() for yr in range(1880,2050,10)])
plt.xlim(ts, te)
#plt.grid()
plt.ylabel('SST Anomaly [$^\circ$C]')
#plt.legend(['HadISST', 'NOAA OI SST'], loc='upper left')
plt.legend(['NOAA OI SST', 'HadISST'], loc='upper left')
plt.title('(B)')

plt.savefig('../../documents/14_Tasmania_2015_2016/figures/' + mhwname + '_timeSeries_Full.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

ev = np.argsort(mhws['duration'])[-1]
ts = date(2012,1,1).toordinal()
te = t[-1] #date(2016,6,1).toordinal()

plt.clf()
plt.subplot(2,1,1)
# Find indices for all ten MHWs before and after event of interest and shade accordingly
for ev0 in range(0, mhws['n_events']):
    t1 = np.where(t==mhws['time_start'][ev0])[0][0]
    t2 = np.where(t==mhws['time_end'][ev0])[0][0]
    plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_ev)
# Find indices for MHW of interest (2011 WA event) and shade accordingly
t1 = np.where(t==mhws['time_start'][ev])[0][0]
t2 = np.where(t==mhws['time_end'][ev])[0][0]
plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_evMax)
## Plot SST, seasonal cycle, threshold, etc
plt.plot(dates, sst, 'k-', linewidth=2)
plt.plot(dates, clim['thresh'], '-', linewidth=2, color=(0.4,0.4,1))
plt.plot(dates, clim['seas'], '0.5', linewidth=2)
plt.plot(dates_had, sst_had, 'ko', markerfacecolor='w', markeredgewidth=2)
plt.xlim(ts, te)
plt.ylim(12, 20)
#plt.grid()
plt.ylabel('SST [$^\circ$C]')
#plt.title('(B) Recent years (2012-2016)')
plt.subplot(2,1,2)
# Find indices for all ten MHWs before and after event of interest and shade accordingly
for ev0 in range(0, mhws['n_events']):
    t1 = np.where(t==mhws['time_start'][ev0])[0][0]
    t2 = np.where(t==mhws['time_end'][ev0])[0][0]
    plt.fill_between(dates[t1:t2+1], (sst - clim['seas'])[t1:t2+1], 0., color=col_ev)
# Find indices for MHW of interest (2011 WA event) and shade accordingly
t1 = np.where(t==mhws['time_start'][ev])[0][0]
t2 = np.where(t==mhws['time_end'][ev])[0][0]
plt.fill_between(dates[t1:t2+1], ((sst - clim['seas']))[t1:t2+1], 0., color=col_evMax)
## Plot SST anomalies
plt.plot(dates, sst - clim['seas'], 'k-', linewidth=2)
plt.plot(dates, np.nan*clim['thresh'], '-', linewidth=2, color=(0.4,0.4,1))
plt.plot(dates, np.nan*clim['seas'], '0.5', linewidth=2)
plt.plot(dates_had, sst_had - sst_had_clim, 'ko', markerfacecolor='w', markeredgewidth=2)
plt.xlim(ts, te)
plt.ylim(-1, 3)
#plt.grid()
plt.ylabel('SST Anomaly [$^\circ$C]')
#plt.legend(['NOAA OI SST', 'Threshold', 'Climatology', 'HadISST'], loc='upper left')
#plt.title('(C)')

plt.savefig('../../documents/14_Tasmania_2015_2016/figures/' + mhwname + '_timeSeries_Recent.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Running standard deviation
years_had = np.arange(1870, 2016+1)
sst_std = np.nan*np.zeros(len(years_had))
sst_std_daily = np.sqrt(ecj.runavg((sst_had-sst_had_clim-ecj.runavg(sst_had-sst_had_clim,31))**2, 31))
cnt = 0
for yr in years_had:
    sst_std[cnt] = np.nanmean(sst_std_daily[year_had==yr])
    cnt += 1

plt.figure()
plt.clf()
#plt.plot(years_had, sst_std, 'k-')
#plt.plot(years_had, sst_std, 'k.')
plt.plot(years_had, sst_std, 'k-o', linewidth=1)
plt.ylim(0, 0.55)
plt.xlim(1880, 2020)
#plt.grid()
plt.ylabel('SST 30-year running std. dev. [$^\circ$C]')

plt.savefig('../../documents/14_Tasmania_2015_2016/figures/RunningStdDev.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

#
# Figure in style of Conversation piece
#

plt.figure(figsize=(12,7))
plt.clf()

ts = date(1982,1,1).toordinal()
te = date(2016,8,1).toordinal()
# Bar plot: Duration
plt.subplot(2,2,3)
evMax = np.argmax(mhws['duration'])
plt.bar(mhws['date_peak'], mhws['duration'], width=150, color=(0.7,0.7,0.7))
plt.bar(mhws['date_peak'][evMax], mhws['duration'][evMax], width=150, color='r')
plt.xlim(ts, te)
plt.ylabel('[days]')
plt.title('Duration of all events')
# Bar plot: Max intensity
plt.subplot(2,2,4)
evMax = np.argmax(np.abs(mhws['intensity_max']))
plt.bar(mhws['date_peak'], mhws['intensity_max'], width=150, color=(0.7,0.7,0.7))
plt.bar(mhws['date_peak'][evMax], mhws['intensity_max'][evMax], width=150, color='r')
plt.xlim(ts, te)
plt.ylabel(r'[$^\circ$C]')
plt.title('Intensity of all events')
# Map
plt.subplot(2,2,1)
domain = [-46.5, 142, -35, 165]
llon_DJF, llat_DJF = np.meshgrid(lon_DJF, lat_DJF)
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
proj.drawcoastlines()
proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[True,False,False,False])
proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,True])
lonproj, latproj = proj(llon_DJF, llat_DJF)
#H = plt.contourf(lonproj, latproj, sst_DJF, levels=[-4,-3,-2.5,-2,-1.5,-1,1,1.5,2,2.5,3,4], cmap=plt.cm.RdBu_r)
H = plt.contourf(lonproj, latproj, sst_DJF, levels=np.arange(0,4.0 + 0.5,0.5), cmap=plt.cm.gist_heat_r)
lonproj, latproj = proj([locations['lon'][0], locations['lon'][0], locations['lon'][1], locations['lon'][1], locations['lon'][0]], [locations['lat'][0], locations['lat'][1], locations['lat'][1], locations['lat'][0], locations['lat'][0]])
plt.plot(lonproj, latproj, 'k-', linewidth=2)
h = plt.colorbar()
h.set_label('$^\circ$C')
plt.title('Average SSTa over DJF')
# Time series
ax = plt.subplot(2,2,2)
evs = np.argsort(mhws['duration'])[-10:]
ev = evs[-1]
for ev0 in np.arange(max(ev-10,0), min(ev+11,mhws['n_events']-1), 1):
    t1 = np.where(t==mhws['time_start'][ev0])[0][0]
    t2 = np.where(t==mhws['time_end'][ev0])[0][0]
    plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_ev)
t1 = np.where(t==mhws['time_start'][ev])[0][0]
t2 = np.where(t==mhws['time_end'][ev])[0][0]
plt.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], color=col_evMax)
plt.plot(dates, sst, 'k-', linewidth=2)
plt.plot(dates, clim['thresh'], col_thresh, linewidth=2)
plt.plot(dates, clim['seas'], col_clim, linewidth=2)
plt.title('2015-2016 Marine Heatwaves')
plt.xlim(mhws['time_start'][ev]-150, mhws['time_end'][ev]+30)
if coldSpells:
    plt.ylim(clim['seas'].min() + mhws['intensity_max'][ev] - 0.5, clim['seas'].max() + 1)
else:
    plt.ylim(clim['seas'].min() - 1, clim['seas'].max() + mhws['intensity_max'][ev] + 0.5)
plt.ylabel(r'SST [$^\circ$C]')
plt.legend(['SST', 'Threshold', 'Climatology'], loc='upper left', fontsize=12)
ax.set_xticks(ax.get_xticks()[1::2])

plt.savefig('mhw_properties/ETas_MHW_2015_2016_Properties.png', bbox_inches='tight', pad_inches=0.5, dpi=300)


#
# Monthly SST stats...
#

# HadSST months that were part of the event
ev = np.argsort(mhws['duration'])[-1] # longest event
tt_ev = np.where((t_had >= mhws['date_start'][ev].toordinal()) * (t_had <= mhws['date_end'][ev].toordinal()))[0]
for ttt in range(len(tt_ev)):
    print (year_had[tt_ev[ttt]], month_had[tt_ev[ttt]])

# Generate SSTa and SSTa over event
ssta_had = sst_had - sst_had_clim
ssta_had_ev = ssta_had[tt_ev]

# Number of previous months which exceed max monthly anomaly
tt_max = np.where(ssta_had > np.max(ssta_had_ev))[0]
print len(tt_max)
for ttt in range(len(tt_max)):
    print (year_had[tt_max[ttt]], month_had[tt_max[ttt]])

# Number of previous N-month runs which exceed this run, in terms of avg. anomaly
ssta_had_run = ecj.runavg(ssta_had, len(ssta_had_ev))
tt_max = np.where(ssta_had_run > np.mean(ssta_had_ev))[0]
print len(tt_max)
for ttt in range(len(tt_max)):
    print (year_had[tt_max[ttt]], month_had[tt_max[ttt]])

# Rank of months in this event
for ttt in tt_ev:
    print len(ecj.nonans(ssta_had)) - np.where(np.argsort(ssta_had)==ttt)[0]

#
# MHW stats
#

ev = np.argmax(mhws['duration']) # Find longest event
print 'Maximum intensity:', mhws['intensity_max'][ev], 'deg. C'
print 'Average intensity:', mhws['intensity_mean'][ev], 'deg. C'
print 'Cumulative intensity:', mhws['intensity_cumulative'][ev], 'deg. C-days'
print 'Duration:', mhws['duration'][ev], 'days'
print 'Start date:', mhws['date_start'][ev].strftime("%d %B %Y")
print 'End date:', mhws['date_end'][ev].strftime("%d %B %Y")

#
# Monthy SST maps (NOAA OI SST)
#

# Generate monthly maps for 2012-2014
year_start = 2012
year_end = 2016

Nly = 35*4 /2
Nlx = 60*4 /2
Tl = (year_end-year_start+1)*12

sst_map_daily = np.zeros((2*Nly,2*Nlx,T))
mhw_map_daily = np.zeros((2*Nly,2*Nlx,T))
sst_map = np.NaN*np.zeros((2*Nly,2*Nlx,Tl))
sst_map_abs = np.NaN*np.zeros((2*Nly,2*Nlx,Tl))
mhw_map = np.zeros((2*Nly,2*Nlx,Tl))
i = np.where(lon > locations['lon_map'])[0][0]
j = np.where(lat > locations['lat_map'])[0][0]
lon_map = lon[(i-Nlx):(i+Nlx)]
lat_map = lat[(j-Nly):(j+Nly)]

# Load data
for tt in range(T):
    try:
        file = header + str(dates[tt].year) + '/avhrr-only-v2.' + str(dates[tt].year) + str(dates[tt].month).zfill(2) + str(dates[tt].day).zfill(2) + '.nc'
        fileobj = Dataset(file, 'r')
    except:
        file = header + str(dates[tt].year) + '/avhrr-only-v2.' + str(dates[tt].year) + str(dates[tt].month).zfill(2) + str(dates[tt].day).zfill(2) + '_preliminary.nc'
        fileobj = Dataset(file, 'r')
    print str(dates[tt].year) + str(dates[tt].month).zfill(2) + str(dates[tt].day).zfill(2)
    sst0 = fileobj.variables['sst'][0,0,(j-Nly):(j+Nly),(i-Nlx):(i+Nlx)].astype(float).data
    sst0[sst0==fill_value] = np.nan
    sst_map_daily[:,:,tt] = sst0 #*scale + offset
    fileobj.close()

# Generate monthly anomalies
for ii in range(2*Nlx):
    print ii, 2*Nlx
    for jj in range(2*Nly):
        if ~np.isnan(np.sum(sst_map_daily[jj,ii,:])):
            mhws_map, clim_map = mhw.detect(t, sst_map_daily[jj,ii,:], coldSpells=coldSpells, climatologyPeriod=[1982,2005])
            #sst_ds, sst_s, sst_beta = ds.deseason_harmonic(sst_map_daily[jj,ii,:], 2, 365.25)
            for ev0 in range(mhws_map['n_events']):
                mhw_map_daily[jj,ii,mhws_map['index_start'][ev0]:mhws_map['index_end'][ev0]+1] = 1.
            sst_ds = sst_map_daily[jj,ii,:] - clim_map['seas']
            tt = 0
            for yr in range(year_start, year_end+1):
                for mth in range(1,12+1):
                    sst_map[jj,ii,tt] = np.nanmean(sst_ds[(year==yr) * (month==mth)])
                    sst_map_abs[jj,ii,tt] = np.nanmean(sst_map_daily[jj,ii,(year==yr) * (month==mth)])
                    mhw_map[jj,ii,tt] = mhw_map_daily[jj,ii,(year==yr) * (month==mth)].sum() / len(mhw_map_daily[jj,ii,(year==yr) * (month==mth)])
                    tt += 1

# outfile = '/data/MHWs/Tasmania_2015_2016/NOAAOISST/sst_map'
# np.savez(outfile, lon_map=lon_map, lat_map=lat_map, sst_map_daily=sst_map_daily, mhw_map_daily=mhw_map_daily, sst_map=sst_map, sst_map_abs=sst_map_abs, mhw_map=mhw_map)

# Load data
data = np.load(outfile + '.npz')
sst_map_daily = data['sst_map_daily']
mhw_map_daily = data['mhw_map_daily']
sst_map = data['sst_map']
sst_map_abs = data['sst_map_abs']
mhw_map = data['mhw_map']
lon_map = data['lon_map']
lat_map = data['lat_map']


#domain = [np.nanmin(lat_map), np.nanmin(lon_map), np.nanmax(lat_map), np.nanmax(lon_map)]
domain = [-46, 143, -36, 161]
labels_month = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
labels_month_ABC = ['(E) Jan', '(F) Feb', '(G) Mar', '(H) Apr', '(I) May', 'Jun', 'Jul', 'Aug', '(A) Sep', '(B) Oct', '(C) Nov', '(D) Dec']
labels_year = np.arange(year_start, year_end+1)
llon_map, llat_map = np.meshgrid(lon_map, lat_map)
wbgyr = colors.ListedColormap(np.genfromtxt('/home/ecoliver/Desktop/python_modules/cmaps/WhiteBlueGreenYellowRed.rgb')/255.)

# SSTAs
fig = plt.figure(figsize=(27,8))
plt.clf()
for tt in range(sst_map.shape[2]):
    AX = plt.subplot(5,12,tt+1)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(llon_map, llat_map)
    H = plt.contourf(lonproj, latproj, sst_map[:,:,tt], levels=range(-5,5+1), cmap=plt.cm.RdBu_r)
    plt.clim(-4, 4)
    if tt+1 <= 12:
        plt.title(labels_month[tt])
    if np.mod(tt+1,12) == 1:
        plt.ylabel(labels_year[(tt+1)/12])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_SSTAnom_Monthly.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# MHW intensities 
sstapos_map = sst_map.copy()
sstapos_map[sstapos_map<=0.] = 0.
plt.clf()
for tt in range(sst_map.shape[2]):
    AX = plt.subplot(5,12,tt+1)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(llon_map, llat_map)
    H = plt.contourf(lonproj, latproj, sstapos_map[:,:,tt], levels=[0.5,1,1.5,2,2.5,3,4,5], cmap=plt.cm.gist_heat_r)
    plt.clim(0.5,4)
    if tt+1 <= 12:
        plt.title(labels_month[tt])
    if np.mod(tt+1,12) == 1:
        plt.ylabel(labels_year[(tt+1)/12])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_MHWInt_Monthly.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# SSTs
plt.clf()
d = 3
for tt in range(sst_map.shape[2]):
    AX = plt.subplot(5,12,tt+1)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(llon_map, llat_map)
    H = plt.contourf(lonproj, latproj, sst_map_abs[:,:,tt], levels=range(9,25+1), cmap=wbgyr)
    #lonproj, latproj = proj(llon_OC, llat_OC)
    #plt.quiver(lonproj[::d,::d], latproj[::d,::d], U_mth_abs[::d,::d,tt], V_mth_abs[::d,::d,tt], scale=5, linewidth=0.2)
    #plt.quiver(lonproj[::d,::d], latproj[::d,::d], U_mth[::d,::d,tt], V_mth[::d,::d,tt], scale=5, linewidth=0.2)
    if tt+1 <= 12:
        plt.title(labels_month[tt])
    if np.mod(tt+1,12) == 1:
        plt.ylabel(labels_year[(tt+1)/12])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_SST_Monthly.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Maps over Sep-Mar only

# SSTAs
fig = plt.figure(figsize=(25,10))
plt.clf()
cnt = 0
for tt in range(4,sst_map.shape[2]-6):
    if (np.mod(tt+1,12) > 3) * (np.mod(tt+1,12) < 9):
        continue
    cnt += 1
    AX = plt.subplot(4,7,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(llon_map, llat_map)
    H = plt.contourf(lonproj, latproj, sst_map[:,:,tt], levels=range(-5,5+1), cmap=plt.cm.RdBu_r)
    plt.clim(-4, 4)
    if tt+1 <= 15:
        plt.title(labels_month[np.mod(tt,12)])
    if np.mod(tt+1,12) == 9:
        plt.ylabel(labels_year[(tt+1)/12])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_SSTAnom_Monthly_SONDJFM.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# MHW Intensities
sstapos_map = sst_map.copy()
sstapos_map[sstapos_map<=0.] = 0.
plt.clf()
cnt = 0
for tt in range(4,sst_map.shape[2]-6):
    if (np.mod(tt+1,12) > 3) * (np.mod(tt+1,12) < 9):
        continue
    cnt += 1
    AX = plt.subplot(4,7,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(llon_map, llat_map)
    H = plt.contourf(lonproj, latproj, sstapos_map[:,:,tt], levels=[0.5,1,1.5,2,2.5,3,4,5], cmap=plt.cm.gist_heat_r)
    plt.clim(0.5,4)
    if tt+1 <= 15:
        plt.title(labels_month[np.mod(tt,12)])
    if np.mod(tt+1,12) == 9:
        plt.ylabel(labels_year[(tt+1)/12])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_MHWInt_Monthly_SONDJFM.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# SSTs
plt.clf()
cnt = 0
d = 3
for tt in range(4,sst_map.shape[2]-6):
    if (np.mod(tt+1,12) > 3) * (np.mod(tt+1,12) < 9):
        continue
    cnt += 1
    AX = plt.subplot(4,7,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(llon_map, llat_map)
    H = plt.contourf(lonproj, latproj, sst_map_abs[:,:,tt], levels=range(9,25+1), cmap=wbgyr)
    #lonproj, latproj = proj(llon_OC, llat_OC)
    #plt.quiver(lonproj[::d,::d], latproj[::d,::d], U_mth_abs[::d,::d,tt], V_mth_abs[::d,::d,tt], scale=10, linewidth=0.4)
    #plt.quiver(lonproj[::d,::d], latproj[::d,::d], U_mth[::d,::d,tt], V_mth[::d,::d,tt], scale=10, linewidth=0.4)
    if tt+1 <= 15:
        plt.title(labels_month[np.mod(tt,12)])
    if np.mod(tt+1,12) == 9:
        plt.ylabel(labels_year[(tt+1)/12])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_SST_Monthly_SONDJFM.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Currents
plt.clf()
cnt = 0
d = 3
for tt in range(4,sst_map.shape[2]-6):
    if (np.mod(tt+1,12) > 3) * (np.mod(tt+1,12) < 9):
        continue
    cnt += 1
    AX = plt.subplot(4,7,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(llon_OC, llat_OC)
    H = plt.contourf(lonproj, latproj, np.sqrt(U_mth_abs[:,:,tt]**2 + V_mth_abs[:,:,tt]**2), levels=np.arange(0.2,1.6+0.2,0.2), cmap=plt.cm.gnuplot2_r)
    plt.clim(0, 2)
    plt.quiver(lonproj[::d,::d], latproj[::d,::d], U_mth_abs[::d,::d,tt], V_mth_abs[::d,::d,tt], scale=10, linewidth=0.4)
    if tt+1 <= 15:
        plt.title(labels_month[np.mod(tt,12)])
    if np.mod(tt+1,12) == 9:
        plt.ylabel(labels_year[(tt+1)/12])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[m s$^{-1}$]')

plt.savefig('mhw_properties/' + mhwname + '_Currents_Monthly_SONDJFM.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# EKE
plt.clf()
cnt = 0
for tt in range(4,EKE_mth.shape[2]-6):
    if (np.mod(tt+1,12) > 3) * (np.mod(tt+1,12) < 9):
        continue
    cnt += 1
    AX = plt.subplot(4,7,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(llon_OC, llat_OC)
    H = plt.contourf(lonproj, latproj, EKE_mth[:,:,tt], levels=1e-8*np.arange(0.1,1.0+0.1,0.1), cmap=plt.cm.gnuplot2_r)
    #plt.clim(0, 2)
    if tt+1 <= 15:
        plt.title(labels_month[np.mod(tt,12)])
    if np.mod(tt+1,12) == 9:
        plt.ylabel(labels_year[(tt+1)/12])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[m$^2$ s$^k{-2}$]')

plt.savefig('mhw_properties/' + mhwname + '_EKE_Monthly_SONDJFM.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Air temperature
plt.clf()
cnt = 0
for tt in range(4,sst_map.shape[2]-6):
    if (np.mod(tt+1,12) > 3) * (np.mod(tt+1,12) < 9):
        continue
    cnt += 1
    AX = plt.subplot(4,7,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(lon_CFS, lat_CFS)
    H = plt.contourf(lonproj, latproj, air_temp_mth_abs[:,:,tt], levels=np.append(np.append(0, np.arange(4,26+1)), 32), cmap=wbgyr)
    plt.clim(4,26)
    if tt+1 <= 15:
        plt.title(labels_month[np.mod(tt,12)])
    if np.mod(tt+1,12) == 9:
        plt.ylabel(labels_year[(tt+1)/12])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_SAT_Monthly_SONDJFM.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Air temperature anomalies
plt.clf()
d = 4
cnt = 0
for tt in range(4,sst_map.shape[2]-6):
    if (np.mod(tt+1,12) > 3) * (np.mod(tt+1,12) < 9):
        continue
    cnt += 1
    AX = plt.subplot(4,7,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(lon_CFS, lat_CFS)
    H = plt.contourf(lonproj, latproj, air_temp_mth[:,:,tt], levels=[-5,-4,-3,-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5,3,4,5], cmap=plt.cm.RdBu_r)
    plt.clim(-4, 4)
    plt.streamplot(lonproj[::d,::d], latproj[::d,::d], u_mth_abs[::d,::d,tt], v_mth_abs[::d,::d,tt], linewidth=np.sqrt(u_mth_abs[::d,::d,tt]**2+v_mth_abs[::d,::d,tt]**2)*0.2, density=1, color='k')
    if tt+1 <= 15:
        plt.title(labels_month[np.mod(tt,12)])
    if np.mod(tt+1,12) == 9:
        plt.ylabel(labels_year[(tt+1)/12])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_SATAnom_Monthly_SONDJFM.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Monthly, 2016

# SSTAs
fig = plt.figure(figsize=(11,8))
plt.clf()
cnt = 0
for tt in range(36+9-1, 48+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    if np.mod(cnt,3) == 1:
        proj.drawparallels(np.arange(np.floor(domain[0]/3)*3,np.ceil(domain[2]/5)*3+1,3), labels=[True,False,False,False])
    else:
        proj.drawparallels(np.arange(np.floor(domain[0]/3)*3,np.ceil(domain[2]/3)*3+1,3), labels=[False,False,False,False])
    if cnt <= 6:
        proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])
    else:
        proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,True])
    lonproj, latproj = proj(llon_map, llat_map)
    H = plt.contourf(lonproj, latproj, sst_map[:,:,tt], levels=range(-5,5+1), cmap=plt.cm.RdBu_r)
    plt.clim(-4, 4)
    #plt.contour(lonproj, latproj, mhw_map[:,:,tt], levels=[0.5], colors='k', linewidth=2)
    plt.text(proj(145.1,-37)[0], proj(145.1,-37)[1], labels_month_ABC[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_SSTAnom_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# MHW Intensities
sstapos_map = sst_map.copy()
sstapos_map[sstapos_map<=0.] = 0.
plt.clf()
cnt = 0
for tt in range(36+9-1, 48+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    if np.mod(cnt,3) == 1:
        proj.drawparallels(np.arange(np.floor(domain[0]/3)*3,np.ceil(domain[2]/5)*3+1,3), labels=[True,False,False,False], dashes=[6,900])
    else:
        proj.drawparallels(np.arange(np.floor(domain[0]/3)*3,np.ceil(domain[2]/3)*3+1,3), labels=[False,False,False,False], linewidth=0.)
    if cnt <= 6:
        proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False], linewidth=0.)
    else:
        proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,True], dashes=[6,900])
    lonproj, latproj = proj(llon_map, llat_map)
    H = plt.contourf(lonproj, latproj, sstapos_map[:,:,tt], levels=[0.5,1,1.5,2,2.5,3,4,5], cmap=plt.cm.gist_heat_r)
    plt.clim(0.5,4)
    plt.contour(lonproj, latproj, mhw_map[:,:,tt], levels=[0.90], colors='k', linewidth=2)
    plt.text(proj(145.1,-37)[0], proj(145.1,-37)[1], labels_month_ABC[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_MHWInt_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.savefig('../../documents/14_Tasmania_2015_2016/figures/' + mhwname + '_MHWInt_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# SSTs
plt.clf()
cnt = 0
for tt in range(36+9-1, 48+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(llon_map, llat_map)
    H = plt.contourf(lonproj, latproj, sst_map_abs[:,:,tt], levels=range(9,25+1), cmap=wbgyr)
    plt.title(labels_month[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_SST_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Currents
plt.clf()
d = 3
cnt = 0
for tt in range(36+9-1, 48+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    if np.mod(cnt,3) == 1:
        proj.drawparallels(np.arange(np.floor(domain[0]/3)*3,np.ceil(domain[2]/5)*3+1,3), labels=[True,False,False,False], dashes=[6,900])
    else:
        proj.drawparallels(np.arange(np.floor(domain[0]/3)*3,np.ceil(domain[2]/3)*3+1,3), labels=[False,False,False,False], linewidth=0.)
    if cnt <= 6:
        proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False], linewidth=0.)
    else:
        proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,True], dashes=[6,900])
    lonproj, latproj = proj(llon_OC, llat_OC)
    H = plt.contourf(lonproj, latproj, np.nanmean(np.sqrt(U_mth_abs[:,:,tt:tt+1+1]**2 + V_mth_abs[:,:,tt:tt+1+1]**2), axis=2), levels=np.arange(0.2,1.2+0.2,0.2), cmap=plt.cm.gnuplot2_r)
    plt.clim(0, 2)
    plt.quiver(lonproj[::d,::d], latproj[::d,::d], np.nanmean(U_mth_abs[::d,::d,tt:tt+1+1], axis=2), np.nanmean(V_mth_abs[::d,::d,tt:tt+1+1], axis=2), scale=10, linewidth=0.4, pivot='middle')
    if cnt == 1:
        plt.quiver(proj(144.2,-37.2)[0], proj(144.2,-37.2)[1], 0.5, 0., scale=10, linewidth=0.1, pivot='middle')
    plt.text(proj(145.1,-37)[0], proj(145.1,-37)[1], labels_month_ABC[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[m s$^{-1}$]')

plt.savefig('mhw_properties/' + mhwname + '_Currents_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.savefig('../../documents/14_Tasmania_2015_2016/figures/' + mhwname + '_Currents_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# EKE
plt.clf()
cnt = 0
for tt in range(36+9-1, 48+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    if np.mod(cnt,3) == 1:
        proj.drawparallels(np.arange(np.floor(domain[0]/3)*3,np.ceil(domain[2]/5)*3+1,3), labels=[True,False,False,False])
    else:
        proj.drawparallels(np.arange(np.floor(domain[0]/3)*3,np.ceil(domain[2]/3)*3+1,3), labels=[False,False,False,False])
    if cnt <= 6:
        proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])
    else:
        proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,True])
    lonproj, latproj = proj(llon_OC, llat_OC)
    H = plt.contourf(lonproj, latproj, np.nanmean(EKE_mth[:,:,tt:tt+1+1], axis=2), levels=1e-8*np.arange(0.1,1.0+0.1,0.1), cmap=plt.cm.gnuplot2_r)
    #plt.clim(0, 2)
    plt.text(proj(145.1,-37)[0], proj(145.1,-37)[1], labels_month_ABC[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[m$^2$ s$^{-2}$]')

plt.savefig('mhw_properties/' + mhwname + '_EKE_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Air temperature
plt.clf()
cnt = 0
for tt in range(36+9-1, 48+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,5), labels=[False,False,False,False])#labels=[True,False,False,False])
    proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])#labels=[False,False,False,True])
    lonproj, latproj = proj(lon_CFS, lat_CFS)
    H = plt.contourf(lonproj, latproj, air_temp_mth_abs[:,:,tt], levels=np.append(np.append(0, np.arange(4,26+1)), 32), cmap=wbgyr)
    plt.clim(4,26)
    plt.title(labels_month[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_SAT_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Air temperature anomalies
plt.clf()
cnt = 0
d = 5
for tt in range(36+9-1, 48+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    if np.mod(cnt,3) == 1:
        proj.drawparallels(np.arange(np.floor(domain[0]/3)*3,np.ceil(domain[2]/5)*3+1,3), labels=[True,False,False,False], dashes=[6,900])
    else:
        proj.drawparallels(np.arange(np.floor(domain[0]/3)*3,np.ceil(domain[2]/3)*3+1,3), labels=[False,False,False,False], linewidth=0.)
    if cnt <= 6:
        proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False], linewidth=0.)
    else:
        proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,True], dashes=[6,900])
    lonproj, latproj = proj(lon_CFS, lat_CFS)
    H = plt.contourf(lonproj, latproj, air_temp_mth[:,:,tt], levels=[-5,-3.5,-2,-1,-0.5,0.5,1,2,3.5,5], cmap=plt.cm.RdBu_r)
    plt.clim(-4.5, 4.5)
    plt.quiver(lonproj[::d,::d], latproj[::d,::d], u_mth_abs[::d,::d,tt], v_mth_abs[::d,::d,tt], scale=75, linewidth=0.4, pivot='middle')
    if cnt == 1:
        plt.quiver(proj(144.2,-37.2)[0], proj(144.2,-37.2)[1], 5., 0., scale=75, linewidth=0.1, pivot='middle')
    #plt.streamplot(lonproj[::d,::d], latproj[::d,::d], u_mth_abs[::d,::d,tt], v_mth_abs[::d,::d,tt], linewidth=np.sqrt(u_mth_abs[::d,::d,tt]**2+v_mth_abs[::d,::d,tt]**2)*0.2, density=1, color='k')
    plt.text(proj(145.1,-37)[0], proj(145.1,-37)[1], labels_month_ABC[np.mod(tt,12)])

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[$^\circ$C]')

plt.savefig('mhw_properties/' + mhwname + '_SATAnom_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=300)
plt.savefig('../../documents/14_Tasmania_2015_2016/figures/' + mhwname + '_SATAnom_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# OceanMAPS Currents
domain = [-46.5, 145, -37, 160]
fig = plt.figure()
d = 4
cnt = 0
for tt in range(36+9-1, 48+4+1):
    cnt += 1
    AX = plt.subplot(3,3,cnt)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
    proj.drawcoastlines()
    if np.mod(cnt,3) == 1:
        proj.drawparallels(np.arange(np.floor(domain[0]/3)*3,np.ceil(domain[2]/5)*3+1,3), labels=[True,False,False,False])
    else:
        proj.drawparallels(np.arange(np.floor(domain[0]/3)*3,np.ceil(domain[2]/3)*3+1,3), labels=[False,False,False,False])
    if cnt <= 6:
        proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,False])
    else:
        proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,5), labels=[False,False,False,True])
    lonproj, latproj = proj(lon_OMAPS, lat_OMAPS)
    H = plt.contourf(lonproj, latproj, np.nanmean(np.sqrt(u_mth_abs[:,:,tt:tt+1+1]**2 + v_mth_abs[:,:,tt:tt+1+1]**2), axis=2), levels=np.arange(0.2,1.6+0.2,0.2), cmap=plt.cm.gnuplot2_r)
    plt.clim(0, 2)
    plt.quiver(lonproj[::d,::d], latproj[::d,::d], np.nanmean(u_mth_abs[::d,::d,tt:tt+1+1], axis=2), np.nanmean(v_mth_abs[::d,::d,tt:tt+1+1], axis=2), scale=5, linewidth=0.4, pivot='middle')
    if cnt == 1:
        plt.quiver(proj(144.2,-37.2)[0], proj(144.2,-37.2)[1], 0.5, 0., scale=10, linewidth=0.1, pivot='middle')
    plt.text(proj(145.1,-37)[0], proj(145.1,-37)[1], labels_month_ABC[np.mod(tt,12)])
    lonproj, latproj = proj([locations['lon'][0], locations['lon'][0], locations['lon'][1], locations['lon'][1], locations['lon'][0]], [locations['lat'][0], locations['lat'][1], locations['lat'][1], locations['lat'][0], locations['lat'][0]])
    plt.plot(lonproj, latproj, 'k-', linewidth=2)

AXPOS = AX.get_position()
CAX = fig.add_axes([AXPOS.x1+0.015, AXPOS.y0, 0.01, AXPOS.y1-AXPOS.y0])
HB = plt.colorbar(H, CAX, orientation='vertical')
HB.set_label(r'[m s$^{-1}$]')

plt.savefig('mhw_properties/' + mhwname + '_OMAPS_Currents_Monthly2016.png', bbox_inches='tight', pad_inches=0.5, dpi=300)



#
# Australia map of MHW SSTa
#

domain = [-50, 100, -5, 180]
domain_Tas = [-44, 144, -39, 150]
llon_DJF, llat_DJF = np.meshgrid(lon_DJF, lat_DJF)

#mhw_DJF = sst_DJF.copy()
#mhw_DJF[mhw_DJF<=0.] = 0.

fig = plt.figure(figsize=(27,8))
plt.clf()
proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='i')
#proj = bm.Basemap(projection='nsper',lon_0=145,lat_0=-35, satellite_height=3000*1000., resolution='i')
#proj = bm.Basemap(width=8000000,height=7000000, resolution='l',projection='aea', lat_1=-50.,lat_2=-20,lon_0=145,lat_0=-35)
proj.drawcoastlines(color='0.1')
proj.fillcontinents(color='0.1')
proj.drawparallels(np.arange(np.floor(domain[0]/5)*5,np.ceil(domain[2]/5)*5+1,20), labels=[True,False,False,False], dashes=[6,900])
proj.drawmeridians(np.arange(np.floor(domain[1]/5)*5,np.ceil(domain[3]/5)*5+1,20), labels=[False,False,False,True], dashes=[6,900])
lonproj, latproj = proj(llon_DJF, llat_DJF)
#H = plt.contourf(lonproj, latproj, mhw_DJF, levels=[1,1.5,2,2.5,3,4], cmap=plt.cm.gist_heat_r)
H = plt.contourf(lonproj, latproj, sst_DJF, levels=[-3.5,-2.5,-2,-1.5,-1,1,1.5,2,2.5,3.5], cmap=plt.cm.RdBu_r)
lonproj, latproj = proj([locations['lon'][0], locations['lon'][0], locations['lon'][1], locations['lon'][1], locations['lon'][0]], [locations['lat'][0], locations['lat'][1], locations['lat'][1], locations['lat'][0], locations['lat'][0]])
plt.plot(lonproj, latproj, 'k-', linewidth=2)
lonproj, latproj = proj(148.233505711, -42.5918436917)
plt.plot(lonproj, latproj, 'o', markerfacecolor='w', markersize=8, markeredgewidth=2)

HB = plt.colorbar()
HB.set_label(r'[$^\circ$C]')

plt.title('(A) Mean DJF Sea Surface Temperature Anomaly')

#plt.axes([0.38, 0.4, 0.2, 0.2])
#proj = bm.Basemap(projection='merc', llcrnrlat=domain_Tas[0], llcrnrlon=domain_Tas[1], urcrnrlat=domain_Tas[2], urcrnrlon=domain_Tas[3], resolution='i')
#proj.drawcoastlines(color='0.1')
#proj.fillcontinents(color='0.1')
#lonproj, latproj = proj(llon_DJF, llat_DJF)
##H = plt.contourf(lonproj, latproj, mhw_DJF, levels=[1,1.5,2,2.5,3,4], cmap=plt.cm.gist_heat_r)
#plt.contourf(lonproj, latproj, sst_DJF, levels=[-4,-3,-2.5,-2,-1.5,-1,1,1.5,2,2.5,3,4], cmap=plt.cm.RdBu_r)
#lonproj, latproj = proj(148.233505711, -42.5918436917)
#plt.plot(lonproj, latproj, 'o', markerfacecolor='w', markersize=8, markeredgewidth=2)

plt.savefig('../../documents/14_Tasmania_2015_2016/figures/' + mhwname + '_Map_AUSNZ_SSTAnom_DJF.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# Area of anomaly regions
dA = np.zeros(sst_DJF.shape) # Area of each grid cell, in km^2
for j in range(sst_DJF.shape[0]-1):
    for i in range(sst_DJF.shape[1]-1):
        dlon = ecj.latlon2km(lon_DJF[i], lat_DJF[j], lon_DJF[i+1], lat_DJF[j])
        dlat = ecj.latlon2km(lon_DJF[i], lat_DJF[j], lon_DJF[i], lat_DJF[j+1])
        dA[j,i] = dlon*dlat

# Region of interest, contains this cell
i0 = 595
j0 = 48
A_Tas = 68401.
A_UK = 243610.

# Are of region that includes that cell
regions, nregions = ndimage.label( (sst_DJF>1.).astype(int) )
(regions == regions[j0,i0])
A1 = dA[(regions == regions[j0,i0])].sum()
regions, nregions = ndimage.label( (sst_DJF>2.).astype(int) )
(regions == regions[j0,i0])
A2 = dA[(regions == regions[j0,i0])].sum()
regions, nregions = ndimage.label( (sst_DJF>3.).astype(int) )
(regions == regions[j0,i0])
A3 = dA[(regions == regions[j0,i0])].sum()

A1/A_Tas
A1/A_UK

A2/A_Tas
A2/A_UK

A3/A_Tas
A3/A_UK


#
# Avg EKE over SEAus
#

EKE_avg = EKE_mth.copy()
EKE_avg[EKE_avg==0] = np.nan
EKE_avg[EKE_avg>1] = np.nan
EKE_avg[U_mth_abs==-32766.0] = np.nan

i_SEAus = (lon_OC >= locations['lon'][0]) * (lon_OC <= locations['lon'][1])
j_SEAus = (lat_OC >= locations['lat'][0]) * (lat_OC <= locations['lat'][1])

EKE_SEAus = np.nanmean(np.nanmean(EKE_avg[j_SEAus,:,:][:,i_SEAus,:], axis=0), axis=0)

plt.clf()
cnt = 0
for year in range(2012,2015+1,1):
    #plt.plot(EKE_SEAus[(cnt*12):(cnt+1)*12])
    mths = range(cnt*12, (cnt+1)*12+1)
    plt.plot(EKE_SEAus[mths[8]-1:mths[9]+7+1], '-o')
    cnt += 1

plt.xlim(0.5, 7.5)
#plt.grid()

plt.legend(['2012/13', '2013/14', '2014/15', '2015/16'], loc='upper left')
plt.ylabel(r'Eddy Kinetic Energy [m$^2$ s$^{-2}$]')
plt.xticks(range(1,7+1), np.append(labels_month[8:], labels_month[:3]))

plt.savefig('../../documents/14_Tasmania_2015_2016/figures/EKE_SEAus.png', bbox_inches='tight', pad_inches=0.5, dpi=300)




