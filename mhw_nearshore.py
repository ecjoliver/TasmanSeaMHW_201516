'''

  Maria Island time series

'''

# Load required modules

import numpy as np
from scipy import io
from datetime import date

from matplotlib import pyplot as plt
from matplotlib import dates as mdates
from matplotlib import colors
import mpl_toolkits.basemap as bm

import deseason as ds
import ecoliver as ecj

import marineHeatWaves as mhw

#
# Globals
#

header = '/data/MHWs/Tasmania_2015_2016/'

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
# Load Craig's temperature logger data
#

outfile = header + 'CraigMundy_temperature_loggers/CraigMundyTempLoggers.npz'
data = np.load(outfile)
temp_CM = data['sst'].item()
t_CM = data['t'].item()
dates_CM = data['dates'].item()
lon_CM = data['lon'].item()
lat_CM = data['lat'].item()
depth_CM = data['depth'].item()

# Sites with long enough records

sites2016 = ['Cape_Peron', 'Schouten_Island', 'Coles_Bay', 'Swansea', 'Bicheno', 'Magistrates_Point', 'Wineglass_Bay', 'George_III_Reef', 'Iron_Pot', 'Mouldies_Hole', 'One_Tree_Point', 'Wedge_Island']
mhws_CM = {}
clim_CM = {}
for site in sites2016:
    temp_CM[site] = ecj.runavg(temp_CM[site], 3)
    missing = np.isnan(temp_CM[site])
    mhws_CM[site], clim_CM[site] = mhw.detect(t_CM[site], temp_CM[site], climatologyPeriod=[dates_CM[site][0].year+1,dates_CM[site][-1].year-1])
    temp_CM[site][missing] = np.nan

#
# Load MITS
#

# BGC data
file = header + 'IMOS_NRSMAI/IMOS_NRSMAI_BGC_daily.mat'
matobj = io.loadmat(file)
t_temp = matobj['t_d'].flatten()
dates_temp = matobj['dates_d']
depth_temp = np.nanmean(matobj['depth_d'], axis=0)
temp = matobj['temp_d']
salt = matobj['salt_d']
dox1 = matobj['dox1_d']
dox2 = matobj['dox2_d']
cphl = matobj['cphl_d']
turb = matobj['turb_d']
t_temp, dates_temp, T_temp, year_temp, month_temp, day_temp, doy_temp = ecj.timevector([dates_temp[0,0], dates_temp[0,1], dates_temp[0,2]], [dates_temp[-1,0], dates_temp[-1,1], dates_temp[-1,2]])

# VEL data
file = header + 'IMOS_NRSMAI/IMOS_NRSMAI_VEL_daily_deTided.mat'
matobj = io.loadmat(file)
t_u = matobj['t_d'].flatten()
dates_u = matobj['dates_d']
depth_u = np.nanmean(matobj['depth_d'], axis=0)
u = matobj['u_d']
v = matobj['v_d']
t_u, dates_u, T_u, year_u, month_u, day_u, doy_u = ecj.timevector([dates_u[0,0], dates_u[0,1], dates_u[0,2]], [dates_u[-1,0], dates_u[-1,1], dates_u[-1,2]])

#
# Calculate deseasonalized data
#

temp_ds = np.nan*np.ones(temp.shape)
salt_ds = np.nan*np.ones(salt.shape)
dox1_ds = np.nan*np.ones(dox1.shape)
dox2_ds = np.nan*np.ones(dox2.shape)
turb_ds = np.nan*np.ones(turb.shape)
cphl_ds = np.nan*np.ones(cphl.shape)
u_ds = np.nan*np.ones(u.shape)
v_ds = np.nan*np.ones(v.shape)

for k in range(temp.shape[1]):
    temp_ds[:,k] = np.array(ds.deseason_harmonic(temp[:,k], 3, 365.25)[0]).flatten()
    salt_ds[:,k] = np.array(ds.deseason_harmonic(salt[:,k], 3, 365.25)[0]).flatten()
    dox1_ds[:,k] = np.array(ds.deseason_harmonic(dox1[:,k], 3, 365.25)[0]).flatten()
    dox2_ds[:,k] = np.array(ds.deseason_harmonic(dox2[:,k], 3, 365.25)[0]).flatten()
    turb_ds[:,k] = np.array(ds.deseason_harmonic(turb[:,k], 3, 365.25)[0]).flatten()
    cphl_ds[:,k] = np.array(ds.deseason_harmonic(cphl[:,k], 3, 365.25)[0]).flatten()

for k in range(u.shape[1]):
    u_ds[:,k] = np.array(ds.deseason_harmonic(u[:,k], 3, 365.25)[0]).flatten()
    v_ds[:,k] = np.array(ds.deseason_harmonic(v[:,k], 3, 365.25)[0]).flatten()

# Some stats
years = np.unique(year_temp)[0:-1]
for year in years:
    DJ = range(np.where((year_temp == year) * (month_temp == 12) * (day_temp == 1))[0][0], np.where((year_temp == year+1) * (month_temp == 1) * (day_temp == 31))[0][0]+1)
    print np.nanmean(temp[DJ,0]), np.nanmean(temp_ds[DJ,0])

years = np.unique(year_u)[0:-1]
for year in years:
    DJ = ((year_u == year) * (month_u == 12)) + ((year_u == year+1) * (month_u == 1))
    print np.nanmean(np.nanmean(v[:,::-3], axis=1)[DJ]), np.nanmean(np.nanmean(v_ds[:,::-3], axis=1)[DJ])

mhws, clim = mhw.detect(t_temp, temp[:,0], climatologyPeriod=[2009,2015], maxPadLength=5)
ev = np.argsort(mhws['duration'])[-1]
mhws['duration'][ev]
mhws['date_start'][ev]
mhws['date_end'][ev]

mhws['intensity_cumulative'][ev]
mhws['intensity_max'][ev]
mhws['intensity_mean'][ev]

#
# Monthly averages
#

# Monthly time/date vectors
date_min = date.fromordinal(min(min(t_u), min(t_temp)))
date_max = date.fromordinal(max(max(t_u), max(t_temp)))
t_mth, dates_mth, T_mth, year_mth, month_mth, day_mth, doy_mth = ecj.timevector([date_min.year, date_min.month, 1], [date_max.year, date_max.month, 28])
t_mth = t_mth[day_mth == 15]
T_mth = len(t_mth)
dates_mth = []
for tt in range(len(t_mth)):
    dates_mth.append(date.fromordinal(t_mth[tt]))

# Calculate monthly means
temp_mth = np.nan*np.zeros((T_mth, len(depth_temp)))
salt_mth = np.nan*np.zeros((T_mth, len(depth_temp)))
dox1_mth = np.nan*np.zeros((T_mth, len(depth_temp)))
dox2_mth = np.nan*np.zeros((T_mth, len(depth_temp)))
turb_mth = np.nan*np.zeros((T_mth, len(depth_temp)))
cphl_mth = np.nan*np.zeros((T_mth, len(depth_temp)))
temp_mth_ds = np.nan*np.zeros((T_mth, len(depth_temp)))
salt_mth_ds = np.nan*np.zeros((T_mth, len(depth_temp)))
dox1_mth_ds = np.nan*np.zeros((T_mth, len(depth_temp)))
dox2_mth_ds = np.nan*np.zeros((T_mth, len(depth_temp)))
turb_mth_ds = np.nan*np.zeros((T_mth, len(depth_temp)))
cphl_mth_ds = np.nan*np.zeros((T_mth, len(depth_temp)))
u_mth = np.nan*np.zeros((T_mth, len(depth_u)))
v_mth = np.nan*np.zeros((T_mth, len(depth_u)))
u_mth_ds = np.nan*np.zeros((T_mth, len(depth_u)))
v_mth_ds = np.nan*np.zeros((T_mth, len(depth_u)))

for k in range(len(depth_temp)):
    for tt in range(T_mth):
        temp_mth[tt,k] = np.nanmean(temp[(year_temp == dates_mth[tt].year) * (month_temp == dates_mth[tt].month), k])
        salt_mth[tt,k] = np.nanmean(salt[(year_temp == dates_mth[tt].year) * (month_temp == dates_mth[tt].month), k])
        dox1_mth[tt,k] = np.nanmean(dox1[(year_temp == dates_mth[tt].year) * (month_temp == dates_mth[tt].month), k])
        dox2_mth[tt,k] = np.nanmean(dox2[(year_temp == dates_mth[tt].year) * (month_temp == dates_mth[tt].month), k])
        turb_mth[tt,k] = np.nanmean(turb[(year_temp == dates_mth[tt].year) * (month_temp == dates_mth[tt].month), k])
        cphl_mth[tt,k] = np.nanmean(cphl[(year_temp == dates_mth[tt].year) * (month_temp == dates_mth[tt].month), k])
        temp_mth_ds[tt,k] = np.nanmean(temp_ds[(year_temp == dates_mth[tt].year) * (month_temp == dates_mth[tt].month), k])
        salt_mth_ds[tt,k] = np.nanmean(salt_ds[(year_temp == dates_mth[tt].year) * (month_temp == dates_mth[tt].month), k])
        dox1_mth_ds[tt,k] = np.nanmean(dox1_ds[(year_temp == dates_mth[tt].year) * (month_temp == dates_mth[tt].month), k])
        dox2_mth_ds[tt,k] = np.nanmean(dox2_ds[(year_temp == dates_mth[tt].year) * (month_temp == dates_mth[tt].month), k])
        turb_mth_ds[tt,k] = np.nanmean(turb_ds[(year_temp == dates_mth[tt].year) * (month_temp == dates_mth[tt].month), k])
        cphl_mth_ds[tt,k] = np.nanmean(cphl_ds[(year_temp == dates_mth[tt].year) * (month_temp == dates_mth[tt].month), k])

for k in range(len(depth_u)):
    for tt in range(T_mth):
        u_mth[tt,k] = np.nanmean(u[(year_u == dates_mth[tt].year) * (month_u == dates_mth[tt].month), k])
        v_mth[tt,k] = np.nanmean(v[(year_u == dates_mth[tt].year) * (month_u == dates_mth[tt].month), k])
        u_mth_ds[tt,k] = np.nanmean(u_ds[(year_u == dates_mth[tt].year) * (month_u == dates_mth[tt].month), k])
        v_mth_ds[tt,k] = np.nanmean(v_ds[(year_u == dates_mth[tt].year) * (month_u == dates_mth[tt].month), k])

# Along-shore velocity
th = np.arctan2(29.63, 103.8) # Angle of coastline
vel = (u*np.sin(th) + v*np.cos(th))
vel_ds = (u_ds*np.sin(th) + v_ds*np.cos(th))
vel_mth = (u_mth*np.sin(th) + v_mth*np.cos(th))
vel_mth_ds = (u_mth_ds*np.sin(th) + v_mth_ds*np.cos(th))

#
# Plots
#

ev = np.argsort(mhws['duration'])[-1]
ts = date(2008,7,31).toordinal()
te = dates_temp[-1].toordinal() #date(2016,4,1).toordinal()

plt.figure(figsize=(16,10))
plt.clf()
plt.subplot(3,1,1)
# Find indices for all ten MHWs before and after event of interest and shade accordingly
for ev0 in range(0, mhws['n_events']):
    t1 = np.where(t_temp==mhws['time_start'][ev0])[0][0]
    t2 = np.where(t_temp==mhws['time_end'][ev0])[0][0]
    plt.fill_between(dates_temp[t1:t2+1], temp[t1:t2+1,0], clim['thresh'][t1:t2+1], color=col_ev)
# Find indices for MHW of interest (2011 WA event) and shade accordingly
t1 = np.where(t_temp==mhws['time_start'][ev])[0][0]
t2 = np.where(t_temp==mhws['time_end'][ev])[0][0]
plt.fill_between(dates_temp[t1:t2+1], temp[t1:t2+1,0], clim['thresh'][t1:t2+1], color=col_evMax)
## Plot SST, seasonal cycle, threshold, etc
plt.plot(dates_temp, temp[:,0], 'k-', linewidth=1)
plt.plot(dates_temp, clim['thresh'], 'g-', linewidth=1)
plt.plot(dates_temp, clim['seas'], '0.5', linewidth=1)
plt.plot(dates_mth, temp_mth[:,0], 'o', markerfacecolor='w', markeredgewidth=1)
plt.xlim(ts, te)
plt.ylim(11.5, 19.5)
plt.grid()
plt.ylabel('T (20 m) [$^\circ$C]')
plt.title('(A) Maria Island Time Series')
plt.subplot(3,1,2)
# Find indices for all ten MHWs before and after event of interest and shade accordingly
for ev0 in range(0, mhws['n_events']):
    t1 = np.where(t_temp==mhws['time_start'][ev0])[0][0]
    t2 = np.where(t_temp==mhws['time_end'][ev0])[0][0]
    plt.fill_between(dates_temp[t1:t2+1], (temp[:,0] - clim['seas'])[t1:t2+1], 0., color=col_ev)
# Find indices for MHW of interest (2011 WA event) and shade accordingly
t1 = np.where(t_temp==mhws['time_start'][ev])[0][0]
t2 = np.where(t_temp==mhws['time_end'][ev])[0][0]
plt.fill_between(dates_temp[t1:t2+1], ((temp[:,0] - clim['seas']))[t1:t2+1], 0., color=col_evMax)
## Plot SST anomalies
plt.plot(dates_temp, temp[:,0] - clim['seas'], 'k-', linewidth=1)
plt.plot(dates_mth, temp_mth_ds[:,0], 'o', markerfacecolor='w', markeredgewidth=1)
#plt.plot(dates_temp, np.nan*clim['thresh'], 'g-', linewidth=1)
#plt.plot(dates_temp, np.nan*clim['seas'], '0.5', linewidth=1)
plt.xlim(ts, te)
plt.ylim(-1.75, 3.25)
plt.grid()
plt.ylabel('T (20 m) Anomaly [$^\circ$C]')
#plt.legend(['Maria Island Temperature (20 m)', 'Monthly means', 'Threshold', 'Climatology'], loc='upper left', fontsize=12)
plt.legend(['Maria Island Temperature (20 m)', 'Monthly means'], loc='upper left', fontsize=12)
plt.title('(B)')
plt.subplot(3,1,3)
#plt.plot(dates_u, np.nanmean(v_ds[:,::-3], axis=1), '-', color='0.6', linewidth=1)
#plt.plot(dates_mth, np.nanmean(v_mth_ds[:,::-3], axis=1), 'k-o', color='k', linewidth=2, markerfacecolor='w', markeredgewidth=2)
plt.plot(dates_u, np.nanmean(v[:,::-3], axis=1), '-', color='0.6', linewidth=1)
plt.plot(dates_mth, np.nanmean(v_mth[:,::-3], axis=1), 'k-o', color='k', linewidth=2, markerfacecolor='w', markeredgewidth=2)
plt.xlim(ts, te)
plt.ylim(-0.2,0.3)
plt.grid()
plt.ylabel('v (depth-averaged) [m s$^{-1}$]')
plt.legend(['Maria Island Velocity (depth-averaged)', 'Monthly means'], loc='upper left', fontsize=12)
plt.title('(C)')

# plt.savefig('../../documents/14_Tasmania_2015_2016/figures/' + mhwname + '_MITS.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

plt.figure(figsize=(16,10))
plt.clf()
plt.subplot(4,2,1)
plt.plot(dates_temp, salt[:,0], 'k-', linewidth=1)
plt.plot(dates_mth, salt_mth[:,0], 'o', markerfacecolor='w', markeredgewidth=1)
plt.grid()
plt.ylabel('Salinity [PSU]')
plt.subplot(4,2,2)
plt.plot(dates_temp, salt_ds[:,0], 'k-', linewidth=1)
plt.plot(dates_mth, salt_mth_ds[:,0], 'o', markerfacecolor='w', markeredgewidth=1)
plt.grid()
plt.ylabel('Salinity anomaly [PSU]')
plt.subplot(4,2,3)
plt.plot(dates_temp, dox1[:,0], 'k-', linewidth=1)
plt.plot(dates_mth, dox1_mth[:,0], 'o', markerfacecolor='w', markeredgewidth=1)
plt.ylim(200, 300)
plt.grid()
plt.ylabel('Dissolved O$_2$ [umol L$^{-1}$]')
plt.subplot(4,2,4)
plt.plot(dates_temp, dox1_ds[:,0], 'k-', linewidth=1)
plt.plot(dates_mth, dox1_mth_ds[:,0], 'o', markerfacecolor='w', markeredgewidth=1)
plt.ylim(-30, 30)
plt.grid()
plt.ylabel('Dissolved O$_2$ anomaly [umol L$^{-1}$]')
plt.subplot(4,2,5)
plt.plot(dates_temp, turb[:,0], 'k-', linewidth=1)
plt.plot(dates_mth, turb_mth[:,0], 'o', markerfacecolor='w', markeredgewidth=1)
plt.grid()
plt.ylabel('Turbidity [NTU]')
plt.subplot(4,2,6)
plt.plot(dates_temp, turb_ds[:,0], 'k-', linewidth=1)
plt.plot(dates_mth, turb_mth_ds[:,0], 'o', markerfacecolor='w', markeredgewidth=1)
plt.grid()
plt.ylabel('Turbidity anomaly [NTU]')
plt.subplot(4,2,7)
plt.semilogy(dates_temp, cphl[:,0], 'k-', linewidth=1)
plt.semilogy(dates_mth, cphl_mth[:,0], 'o', markerfacecolor='w', markeredgewidth=1)
plt.grid()
plt.ylabel('Chlorophyll [mg m$^{-3}$]')
plt.subplot(4,2,8)
plt.plot(dates_temp, cphl_ds[:,0], 'k-', linewidth=1)
plt.plot(dates_mth, cphl_mth_ds[:,0], 'o', markerfacecolor='w', markeredgewidth=1)
plt.ylim(-1.5, 3.5)
plt.grid()
plt.ylabel('Chlorophyll anomaly [mg m$^{-3}$]')

# plt.savefig('../../documents/14_Tasmania_2015_2016/figures/' + mhwname + '_MITS_SaltO2TurbChl.png', bbox_inches='tight', pad_inches=0.5, dpi=150)


#
# Proper figure
#

t_CM['Maria_Island'] = t_temp
dates_CM['Maria_Island'] = dates_temp
temp_CM['Maria_Island'] = temp[:,0]
mhws_CM['Maria_Island'] = mhws
clim_CM['Maria_Island'] = clim
lon_CM['Maria_Island'] = 148.233505711
lat_CM['Maria_Island'] = 42.5918436917
depth_CM['Maria_Island'] = 20.

# Order sites by increasing latitude
sites0 = mhws_CM.keys()
lats0 = []
for site in sites0:
    lats0.append(lat_CM[site])

sites = list(np.array(sites0)[np.argsort(lats0)])

# DOES NOT WORK QUITE PROPERLY
#sites0 = mhws_CM.keys()
#sites = []
#for site in sites0:
#    if len(sites) == 0:
#        sites.append(site)
#    elif len(sites) == 1:
#        if lat_CM[site] > lat_CM[sites[-1]]:
#            sites.append(site)
#        else:
#            sites.insert(0, site)
#    else:
#        for i in range(len(sites)-1):
#            if (lat_CM[site] > lat_CM[sites[i]]) + (lat_CM[site] < lat_CM[sites[i+1]]):
#                sites.insert(i, site)
#                break
#        #else:
#        #    sites.insert(0, site)

# Some stats
for site in mhws_CM.keys():
    tt1 = np.where(~np.isnan(temp_CM[site]))[0][0]
    print site, dates_CM[site][tt1], dates_CM[site][-1], lon_CM[site], lat_CM[site], depth_CM[site], 100. - 100.*np.isnan(temp_CM[site][tt1:]).sum() / len(temp_CM[site][tt1:])

# MHW stats
for site in mhws_CM.keys():
    evDur = np.argsort(mhws_CM[site]['duration'])[-1]
    print site, mhws_CM[site]['duration'][evDur], mhws_CM[site]['intensity_max'][evDur],  mhws_CM[site]['intensity_mean'][evDur], mhws_CM[site]['date_start'][evDur], mhws_CM[site]['date_end'][evDur]


domain = [-44, 144, -40, 149]
domain_draw = [-44, 144, -40, 149]
panel = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)', '(j)', '(k)', '(l)', '(m)', '(n)', '(o)']

plt.figure(figsize=(16,11))
plt.clf()
year0 = 2008
cnt = 0
for site in sites:
    cnt += 1
    evDur = np.argsort(mhws_CM[site]['duration'])[-1]
    evMax = np.argsort(mhws_CM[site]['intensity_max'])[-1]
    AX = plt.subplot(7,2,cnt)
    AXPOS = AX.get_position()
    # Find indices for all ten MHWs before and after event of interest and shade accordingly
    for ev0 in range(0, mhws_CM[site]['n_events']):
        t1 = np.where(t_CM[site]==mhws_CM[site]['time_start'][ev0])[0][0]
        t2 = np.where(t_CM[site]==mhws_CM[site]['time_end'][ev0])[0][0]
        plt.fill_between(dates_CM[site][t1:t2+1], (temp_CM[site] - clim_CM[site]['seas'])[t1:t2+1], 0., color=col_ev)
    # Find indices for most intense MHW and shade accordingly
    t1 = np.where(t_CM[site]==mhws_CM[site]['time_start'][evMax])[0][0]
    t2 = np.where(t_CM[site]==mhws_CM[site]['time_end'][evMax])[0][0]
    plt.fill_between(dates_CM[site][t1:t2+1], ((temp_CM[site] - clim_CM[site]['seas']))[t1:t2+1], 0., color=col_evMax)
    ## Plot SST anomalies
    plt.plot(dates_CM[site], temp_CM[site] - clim_CM[site]['seas'], 'k-', linewidth=1)
    # Find indices for longest MHW and plot accordingly
    t1 = np.where(t_CM[site]==mhws_CM[site]['time_start'][evDur])[0][0]
    t2 = np.where(t_CM[site]==mhws_CM[site]['time_end'][evDur])[0][0]
    plt.plot(dates_CM[site][t1:t2+1], ((temp_CM[site] - clim_CM[site]['seas']))[t1:t2+1], color=(0.05, 0.05, 204./255), linewidth=1.5)
    plt.xlim(date.toordinal(date(year0,1,1)), date.toordinal(date(2016,6,1)))
    plt.ylim(-2, 3.5)
    #plt.grid()
    plt.text(date.toordinal(date(year0,2,1)), 2.4+0.5, panel[cnt-1] + ' ' + site.replace('_', ' '), fontsize=13)
    if cnt < 13:
        AX.set_xticklabels([])
    if np.mod(cnt,2) == 0:
        AX.set_yticklabels([])
    else:
        plt.ylabel('SST Anomaly [$^\circ$C]')
    # Map
    AXPOS_map = AXPOS
    AXPOS_map.x0 = AXPOS_map.x0 + 0.005
    AXPOS_map.x1 = AXPOS_map.x1 - 0.27
    AXPOS_map.y0 = AXPOS_map.y0 + 0.01
    AXPOS_map.y1 = AXPOS_map.y1 - 0.01
    AX = plt.axes(AXPOS_map)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='l')
    proj.drawcoastlines()
    #proj.drawparallels(np.arange(domain_draw[0],domain_draw[2]+1,5), labels=[True,False,False,False])
    #proj.drawmeridians(np.arange(domain_draw[1],domain_draw[3]+1,5), labels=[False,False,False,True])
    lonproj, latproj = proj(lon_CM[site], -1.*lat_CM[site])
    plt.plot(lonproj, latproj, 'o', markeredgecolor='k', markeredgewidth=1, markerfacecolor='r', markersize=7)

    AX = plt.subplot(7,2,14)
    AXPOS = AX.get_position()
    AX.yaxis.tick_right()
    AX.yaxis.set_label_position('right')
    plt.plot(dates_u, np.nanmean(v[:,::-3], axis=1), '-', color='0.6', linewidth=1)
    plt.plot(dates_mth, np.nanmean(v_mth[:,::-3], axis=1), 'k-o', color='k', linewidth=2, markerfacecolor='w', markeredgewidth=2)
    plt.xlim(date.toordinal(date(year0,1,1)), date.toordinal(date(2016,6,1)))
    plt.ylim(-0.2,0.3)
    #plt.grid()
    plt.ylabel('v [m s$^{-1}$]')
    plt.text(date.toordinal(date(year0,2,1)), 0.24, '(n) Maria Island (ADCP)', fontsize=13)
    # Map
    AXPOS_map = AXPOS
    AXPOS_map.x0 = AXPOS_map.x0 + 0.005
    AXPOS_map.x1 = AXPOS_map.x1 - 0.27
    AXPOS_map.y0 = AXPOS_map.y0 + 0.01
    AXPOS_map.y1 = AXPOS_map.y1 - 0.01
    AX = plt.axes(AXPOS_map)
    proj = bm.Basemap(projection='merc', llcrnrlat=domain[0], llcrnrlon=domain[1], urcrnrlat=domain[2], urcrnrlon=domain[3], resolution='l')
    proj.drawcoastlines()
    #proj.drawparallels(np.arange(domain_draw[0],domain_draw[2]+1,5), labels=[True,False,False,False])
    #proj.drawmeridians(np.arange(domain_draw[1],domain_draw[3]+1,5), labels=[False,False,False,True])
    lonproj, latproj = proj(lon_CM['Maria_Island'], -1.*lat_CM['Maria_Island'])
    plt.plot(lonproj, latproj, 'o', markeredgecolor='k', markeredgewidth=1, markerfacecolor='r', markersize=7)

# plt.savefig('../../documents/14_Tasmania_2015_2016/figures/' + mhwname + '_nearshore.png', bbox_inches='tight', pad_inches=0.5, dpi=300)



