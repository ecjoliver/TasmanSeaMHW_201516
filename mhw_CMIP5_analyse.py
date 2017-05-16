'''

  Generate data regarding extreme SSTs from ACCESS 1.3 GCM

'''

import numpy as np
import scipy as sp
from scipy import stats
from scipy import signal
from scipy import io
import scipy.optimize as opt
from datetime import date
import os
import sys
import ecoliver as ecj
import deseason as ds

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

import marineHeatWaves as mhw

#
# Load data
#

# Which climate models
periods = ['hist', 'rcp']
models = ['ACCESS1-3', 'CSIRO-Mk3-6-0', 'CNRM-CM5', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'CanESM2']
#models.remove(sys.argv[-1])
#models.remove('CNRM-CM5')
#models = [sys.argv[-1]]

#whichBox = '147_155_-45_-37' # Old
whichBox = '147_157_-46_-39' # New

header = '/data/MHWs/Tasmania_2015_2016/CMIP5/' + whichBox + '/'
#header = '/home/ecoliver/Desktop/data/extreme_SSTs/Tasmania/'

# Load
T_ts = {}
t = {}
dates = {}
T = {}
year = {}
month = {}
day = {}
doy = {}
Ly = {}
NRUNS = {}
NENS = {}
nyears = {}
for period in periods:
    T_ts[period] = {}
    t[period] = {}
    dates[period] = {}
    T[period] = {}
    year[period] = {}
    month[period] = {}
    day[period] = {}
    doy[period] = {}
    Ly[period] = {}
    NRUNS[period] = {}
    NENS[period] = {}
    nyears[period] = {}
    for model in models:
        data = np.load(header + 'mhw_' + model + '_' + period + '.npz')
        lon = data['lon']
        lat = data['lat']
        T_ts[period][model] = data['T_ts_' + period]
        t[period][model] = data['t']
        dates[period][model] = [date.fromordinal(tt.astype(int)) for tt in t[period][model]]
        T[period][model] = data['T']
        year[period][model] = data['year']
        month[period][model] = data['month']
        day[period][model] = data['day']
        doy[period][model] = data['doy']
        Ly[period][model] = data['Ly']
        NRUNS[period][model] = int(data['NRUNS']) #T_ts[period][model].shape[1]
        NENS[period][model] = data['NENS'].tolist() #T_ts[period][model].shape[2]
        nyears[period][model] = len(np.unique(year[period]))

NMODS = len(models)
hist = 0
histNat = 1
rcp45 = 0
rcp85 = 1

# Convert Kelvin to Celcius where necessary
for per in periods:
    for model in models:
        for run in range(T_ts[per][model].shape[1]):
            for iens in range(NENS[per][model][run]): #T_ts[per][model].shape[2]):
                print per, model, run, iens
                if np.nanmedian(T_ts[per][model][:,run,iens]) > 150: # Check if units are obviously K
                    T_ts[per][model][:,run,iens] -= 273.15
                    print np.nanmedian(T_ts[per][model][:,run,iens])

# Observations
data = np.load('/data/MHWs/Tasmania_2015_2016/NOAAOISST/' + whichBox + '/mhw_SEAus_NOAAOISST_HadSST.npz')
#data = np.load('/home/ecoliver/Desktop/data/extreme_SSTs/Tasmania/mhw_SEAus_NOAAOISST_HadSST.npz')
t_obs = data['t']
T_obs = data['T']
year_obs = data['year']
month_obs = data['month']
day_obs = data['day']
doy_obs = data['doy']
sst_obs = data['sst']
year_had = data['year_had']
sst_had = data['sst_had']

# Model-specific fixes
# CSIRO-Mk3-6-0
# One errant time series, pad the fucker
if np.in1d(models, 'CSIRO-Mk3-6-0').sum():
    model = 'CSIRO-Mk3-6-0'
    T_ts['hist'][model][:,histNat,8] = ecj.pad(T_ts['hist'][model][:,histNat,8])
# HadGEM2-ES
# historical ensembles 0 and 3 need ~30 days to finish of 2005.  Take from other ensembles
if np.in1d(models, 'HadGEM2-ES').sum():
    model = 'HadGEM2-ES'
    nov30 = np.where((year['hist'][model]==2005) * (month['hist'][model]==11) * (day['hist'][model]==30))[0][0]
    T_ts['hist'][model][nov30+1:,hist,0] = T_ts['hist'][model][nov30+1:,hist,1] + (T_ts['hist'][model][nov30,hist,0] - T_ts['hist'][model][nov30,hist,1])
    T_ts['hist'][model][nov30+1:,hist,3] = T_ts['hist'][model][nov30+1:,hist,2] + (T_ts['hist'][model][nov30,hist,3] - T_ts['hist'][model][nov30,hist,2])
# 'CNRM-CM5'
# Remove ensembles that seem to have CRAZY values (1e35!)
if np.in1d(models, 'CNRM-CM5').sum():
    model = 'CNRM-CM5'
    tmp = np.nan*np.zeros((T['hist'][model], NRUNS['hist'][model], 5))
    tmp[:,hist,0] = T_ts['hist']['CNRM-CM5'][:,hist,7]
    tmp[:,histNat,:] = T_ts['hist']['CNRM-CM5'][:,histNat,1:5+1]
    T_ts['hist']['CNRM-CM5'] = tmp.copy()
    NENS['hist'][model][hist] = 1
    NENS['hist'][model][histNat] = 5
    del(tmp)

#
# Bias correct?
#

# Bias correct non-seasonal, non-linear-trend signal spectrally and by variance
sst_obs_dtr = sst_obs - (ecj.trend(t_obs, sst_obs)[0] + ecj.trend(t_obs, sst_obs)[1]*(t_obs-t_obs.mean()))
sst_obs_ds, s, beta = ds.deseason_harmonic(sst_obs_dtr, 2, 365.25)
sd_obs = np.std(sst_obs_ds[year_obs<=2005])
##fft_obs = np.fft.fft(sst_obs_ds[year_obs<=2005])
#fft_obs = np.fft.fft(sst_obs_ds[year_obs<=2005] - np.nanmean(sst_obs_ds[year_obs<=2005]))
#freqs_obs = np.fft.fftfreq(len(sst_obs_ds[year_obs<=2005]))
#L = 2*((np.floor(0.1*len(freqs_obs)).astype(int)+1)/2) + 1 # bandwidth smoothing window ~ 10% series length
#L = 2*((np.floor(0.3*len(freqs_obs)).astype(int)+1)/2) + 1 # bandwidth smoothing window ~ 30% series length
#L = 2*((np.floor(0.5*len(freqs_obs)).astype(int)+1)/2) + 1

# Compare spectral densities of historical runs against obs, over 1982-2005, and get mean ratio for each model
#bias_gain = {}
#for model in models:
#    bias_gain0 = np.zeros((len(freqs_obs), NENS['hist'][model][hist]))
#    for ens in range(NENS['hist'][model][hist]):
#        T_ts_dtr = T_ts['hist'][model][:,hist,ens] - (ecj.trend(t['hist'][model], T_ts['hist'][model][:,hist,ens])[0] + ecj.trend(t['hist'][model], T_ts['hist'][model][:,hist,ens])[1]*(t['hist'][model]-t['hist'][model].mean()))
#        T_ts_ds, s, beta = ds.deseason_harmonic(T_ts_dtr, 2, Ly[per][model])
#        # spectral gain
#        #fft_mod = np.fft.fft(T_ts_ds[(year['hist'][model]>=1982)*(year['hist'][model]<=2005)])
#        fft_mod = np.fft.fft(T_ts_ds[(year['hist'][model]>=1982)*(year['hist'][model]<=2005)] - np.nanmean(T_ts_ds[(year['hist'][model]>=1982)*(year['hist'][model]<=2005)]))
#        freqs_mod = np.fft.fftfreq(len(T_ts_ds[(year['hist'][model]>=1982)*(year['hist'][model]<=2005)]))
#        bias_gain0[:,ens] = np.sqrt(ecj.runavg_periodic(np.abs(fft_obs)[:,0]**2, L)/sp.interpolate.interp1d(freqs_mod, ecj.runavg_periodic(np.abs(fft_mod)[:,0]**2, L), fill_value="extrapolate")(freqs_obs)) # actually sqrt(gain)
#    # Ensemble means
#    bias_gain[model] = np.mean(bias_gain0, axis=1)

## FFT TEST
#model = 'CNRM-CM5'
#ens = 0
#T_ts_dtr = T_ts['hist'][model][:,hist,ens] - (ecj.trend(t['hist'][model], T_ts['hist'][model][:,hist,ens])[0] + ecj.trend(t['hist'][model], T_ts['hist'][model][:,hist,ens])[1]*(t['hist'][model]-t['hist'][model].mean()))
#T_ts_ds, s, beta = ds.deseason_harmonic(T_ts_dtr, 2, Ly[per][model])
#ts_obs = sst_obs_ds[year_obs<=2005]
#ts_mod = T_ts_ds[(year['hist'][model]>=1982)*(year['hist'][model]<=2005)]
#fft_obs = np.fft.fft(ts_obs)
#fft_mod = np.fft.fft(ts_mod)
#freqs = np.fft.fftfreq(len(ts_mod))
#L = 801
#gain = np.zeros(fft_mod.shape)
#gain[:,0] = ecj.runavg_periodic(np.abs(fft_obs)[:,0], L)/ecj.runavg_periodic(np.abs(fft_mod)[:,0], L) # actually sqrt(gain)
#fft_modBC = gain*fft_mod
#ts_modBC = np.real(np.fft.ifft(fft_modBC)) # BC = bias-corrected
#plt.clf()
#plt.plot(t_obs[year_obs<=2005], ts_obs, 'k-')
#plt.plot(t['hist'][model][(year['hist'][model]>=1982)*(year['hist'][model]<=2005)], ts_mod, 'b-')
#plt.plot(t['hist'][model][(year['hist'][model]>=1982)*(year['hist'][model]<=2005)], ts_modBC, 'r-')
##
#per = 'hist'
#run = 0
#T_ts_dtr = T_ts[per][model][:,run,ens] - (ecj.trend(t[per][model], T_ts[per][model][:,run,ens])[0] + ecj.trend(t[per][model], T_ts[per][model][:,run,ens])[1]*(t[per][model]-t[per][model].mean()))
#T_ts_ds, s, beta = ds.deseason_harmonic(T_ts_dtr, 2, Ly[per][model])
#ts_mod_full = T_ts_ds.copy()
#fft_mod_full = np.fft.fft(ts_mod_full)
#freqs_full = np.fft.fftfreq(len(ts_mod_full))
#gain_full = np.zeros(fft_mod_full.shape)
#gain_full[:,0] = sp.interpolate.interp1d(freqs, gain[:,0], fill_value="extrapolate")(freqs_full)
#fft_modBC_full = gain_full*fft_mod_full
#ts_modBC_full = np.real(np.fft.ifft(fft_modBC_full)) # BC = bias-corrected
#plt.clf()
#plt.plot(t_obs[year_obs<=2005], ts_obs, 'k-')
#plt.plot(t['hist'][model][(year['hist'][model]>=1982)*(year['hist'][model]<=2005)], ts_mod, 'c.')
#plt.plot(t['hist'][model][(year['hist'][model]>=1982)*(year['hist'][model]<=2005)], ts_modBC, 'r.')
#plt.plot(t['hist'][model], ts_mod_full, 'b-')
#plt.plot(t['hist'][model], ts_modBC_full, 'r-')
#plt.plot(t['hist'][model], ts_modBC_full*np.std(ts_obs)/np.std(ts_modBC_full[(year['hist'][model]>=1982)*(year['hist'][model]<=2005)]), 'm-')
##

# Correct all model runs (spectral)
#for per in periods:
#    for model in models:
#        for run in range(NRUNS[per][model]):
#            for ens in range(NENS[per][model][run]):
#                # First calculate linear trend and seasonal cycle
#                T_ts_dtr = T_ts[per][model][:,run,ens] - (ecj.trend(t[per][model], T_ts[per][model][:,run,ens])[0] + ecj.trend(t[per][model], T_ts[per][model][:,run,ens])[1]*(t[per][model]-t[per][model].mean()))
#                T_ts_ds, s, beta = ds.deseason_harmonic(T_ts_dtr, 2, Ly[per][model])
#                # Bias correct spectral density
#                fft_mod_full = np.fft.fft(T_ts_ds)
#                freqs_full = np.fft.fftfreq(len(T_ts_ds))
#                gain_full = np.zeros(fft_mod_full.shape)
#                gain_full[:,0] = sp.interpolate.interp1d(freqs_obs, bias_gain[model], fill_value="extrapolate")(freqs_full)
#                #gain_full[np.abs(freqs_full) < np.abs(freqs_obs[freqs_obs!=0.]).min()] = 1.
#                #gain_full[np.abs(freqs_full) < 0.28] = 1.
#                #gain_full[:,0] = ecj.runavg_periodic(gain_full, 1001)
#                fft_modBC_full = gain_full*fft_mod_full
#                T_ts_ds_BC = np.real(np.fft.ifft(fft_modBC_full))[:,0] # BC = bias-corrected
#                T_ts_ds_BC = T_ts_ds_BC - np.nanmean(T_ts_ds_BC) + np.nanmean(np.array(T_ts_ds)[:,0])
#                # Next redefine time series as sum of trend, seasonal cycle, and bias-corrected de-seasoned remainder
#                T_ts[per][model][:,run,ens] = ecj.trend(t[per][model], T_ts[per][model][:,run,ens])[0] + ecj.trend(t[per][model], T_ts[per][model][:,run,ens])[1]*(t[per][model]-t[per][model].mean()) + np.array(s)[:,0] + T_ts_ds_BC

# Compare std dev of historical runs against obs, over 1982-2005, and get mean ratio for each model
bias_std = {}
for model in models:
    bias_std0 = np.zeros(NENS['hist'][model][hist])
    for ens in range(NENS['hist'][model][hist]):
        T_ts_dtr = T_ts['hist'][model][:,hist,ens] - (ecj.trend(t['hist'][model], T_ts['hist'][model][:,hist,ens])[0] + ecj.trend(t['hist'][model], T_ts['hist'][model][:,hist,ens])[1]*(t['hist'][model]-t['hist'][model].mean()))
        T_ts_ds, s, beta = ds.deseason_harmonic(T_ts_dtr, 2, Ly[per][model])
        # std. dev.
        bias_std0[ens] = np.std(sst_obs_ds[year_obs<=2005]) / np.std(T_ts_ds[(year['hist'][model]>=1982)*(year['hist'][model]<=2005)])
    # Ensemble means
    bias_std[model] = bias_std0.mean()

# Correct all model runs (variance)
for per in periods:
    for model in models:
        for run in range(NRUNS[per][model]):
            for ens in range(NENS[per][model][run]):
                # First calculate linear trend and seasonal cycle
                T_ts_dtr = T_ts[per][model][:,run,ens] - (ecj.trend(t[per][model], T_ts[per][model][:,run,ens])[0] + ecj.trend(t[per][model], T_ts[per][model][:,run,ens])[1]*(t[per][model]-t[per][model].mean()))
                T_ts_ds, s, beta = ds.deseason_harmonic(T_ts_dtr, 2, Ly[per][model])
                # Bias correct variance
                T_ts_ds_BC = bias_std[model]*np.array(T_ts_ds)[:,0]
                # Next redefine time series as sum of trend, seasonal cycle, and bias-corrected de-seasoned remainder
                T_ts[per][model][:,run,ens] = ecj.trend(t[per][model], T_ts[per][model][:,run,ens])[0] + ecj.trend(t[per][model], T_ts[per][model][:,run,ens])[1]*(t[per][model]-t[per][model].mean()) + np.array(s)[:,0] + T_ts_ds_BC

# Attempt at skewness and kurtosis
#sig_obs = np.std(sst_obs_ds[year_obs<=2005])
#sig_mod = np.std(T_ts_ds[(year['hist'][model]>=1982)*(year['hist'][model]<=2005)])
#gam_obs = stats.skew(sst_obs_ds[year_obs<=2005])
#gam_mod = stats.skew(T_ts_ds[(year['hist'][model]>=1982)*(year['hist'][model]<=2005)])
#kap_obs = stats.kurtosis(sst_obs_ds[year_obs<=2005])
#kap_mod = stats.kurtosis(T_ts_ds[(year['hist'][model]>=1982)*(year['hist'][model]<=2005)])
#
#obs = sst_obs_ds[year_obs<=2005] - np.mean(sst_obs_ds[year_obs<=2005])
#x = np.array(T_ts_ds[(year['hist'][model]>=1982)*(year['hist'][model]<=2005)]) - np.mean(np.array(T_ts_ds[(year['hist'][model]>=1982)*(year['hist'][model]<=2005)]))
#z = x*sig_obs/sig_mod # scale variance
#z = np.sign(x)*np.abs(x**3 * (sig_mod**3/sig_obs**3) * (gam_mod/gam_obs))**(1./3) # scale skewness (changes variance)
#z = np.sign(x)*np.abs(x**4 * (sig_mod**4/sig_obs**4) * (kap_mod/kap_obs))**(1./4) # scale kurtosis (changes variance)

#
# DIAGNOSE
#
plt.figure()
plt.clf()
per = 'hist'
run = 0 # historical
for m in range(NMODS):
    plt.subplot(4,2,m+1)
    plt.plot(t[per][models[m]], T_ts[per][models[m]][:,run,:])
    plt.legend(np.arange(NENS[per][models[m]][run]).astype(str))
    plt.title(models[m] + ', run = ' + per + '-' + str(run))
plt.clf()
per = 'hist'
run = 1 # historicalNat
for m in range(NMODS):
    plt.subplot(4,2,m+1)
    plt.plot(t[per][models[m]], T_ts[per][models[m]][:,run,:])
    plt.legend(np.arange(NENS[per][models[m]][run]).astype(str))
    plt.title(models[m] + ', run = ' + per + '-' + str(run))
plt.clf()
per = 'rcp'
run = 1 # RCP8.5
for m in range(NMODS):
    plt.subplot(4,2,m+1)
    plt.plot(t[per][models[m]], T_ts[per][models[m]][:,run,:])
    plt.legend(np.arange(NENS[per][models[m]][run]).astype(str))
    plt.title(models[m] + ', run = ' + per + '-' + str(run))
#
model = 'CNRM-CM5'
plt.clf()
per = 'hist'
run = 0 # historical
for ens in range(NENS[per][model][run]):
    plt.subplot(5,2,ens+1)
    plt.plot(t[per][model], T_ts[per][model][:,run,ens])
    plt.title(str(ens))
# plt.savefig('mhw_CMIP5/ts_' + model + '_hist' + str(run) + '.png', bbox_inches='tight', pad_inches=0.5, dpi=150)
#
# END
#

#
# Analysis
#

# MHW approach

#mhws_obs, clim = mhw.detect(t_obs[year_obs<=2005], sst_obs[year_obs<=2005], climatologyPeriod=[1982,2005])
#mhws_obs_2016, clim = mhw.detect(t_obs, sst_obs, climatologyPeriod=[1982,2005])

climPeriod = [1911, 1940]
climPeriod = [1881, 1910]
#climPeriod = [1961, 1990]
#climPeriod = [1982, 2005]

# MHW detection, using (1911-1940)->(1982-2005) warming signal frmo HadISST (dt) to offset time series, and (1982-2005) as bsae period
dt = np.nanmean(sst_had[(year_had>=1982)*(year_had<=2005)]) - np.nanmean(sst_had[(year_had>=climPeriod[0])*(year_had<=climPeriod[1])])
mhws_obs, clim = mhw.detect(t_obs[year_obs<=2005], sst_obs[year_obs<=2005]+dt, climatologyPeriod=[1982,2005], alternateClimatology=[t_obs, sst_obs])
mhws_obs_2016, clim = mhw.detect(t_obs, sst_obs+dt, climatologyPeriod=[1982,2005], alternateClimatology=[t_obs, sst_obs])

## King et al. 2015 correction to histNat
#for model in models:
#    tt_climPeriod = (year['hist'][model]>=climPeriod[0]) * (year['hist'][model]<=climPeriod[1])
#    tt_histPeriod = (year['hist'][model]>=1911) * (year['hist'][model]<=1940)
#    dT = (np.nanmean(T_ts['hist'][model][tt_histPeriod, hist, :]) - np.nanmean(T_ts['hist'][model][tt_climPeriod, hist, :])) - (np.nanmean(T_ts['hist'][model][tt_histPeriod, histNat, :]) - np.nanmean(T_ts['hist'][model][tt_climPeriod, histNat, :]))
#    T_ts['hist'][model][:,histNat,:] = T_ts['hist'][model][:,histNat,:] + dT
#    print model, dT

per = 'hist'
mhws_hist_recent = {}
mhws_hist = {}
mhws_histNat = {}
clim_hist_recent = {}
clim_hist = {}
clim_histNat = {}
which_pre2005 = {}
for model in models:
    mhws_hist_recent[model] = {}
    mhws_hist[model] = {}
    mhws_histNat[model] = {}
    clim_hist_recent[model] = {}
    clim_hist[model] = {}
    clim_histNat[model] = {}
    which_pre2005[model] = (year[per][model] <= 2005)
    for iens in range(NENS[per][model][hist]):
        mhws_hist_recent[model][iens], clim_hist_recent[model][iens] = mhw.detect(t[per][model][(year[per][model]>=1982) * which_pre2005[model]], T_ts[per][model][(year[per][model]>=1982) * which_pre2005[model],hist,iens], climatologyPeriod=[1982,2005])
        mhws_hist[model][iens], clim_hist[model][iens] = mhw.detect(t[per][model][which_pre2005[model]], T_ts[per][model][which_pre2005[model],hist,iens], climatologyPeriod=climPeriod)
    for iens in range(NENS[per][model][histNat]):
        mhws_histNat[model][iens], clim_histNat[model][iens] = mhw.detect(t[per][model][which_pre2005[model]], T_ts[per][model][which_pre2005[model],histNat,iens], climatologyPeriod=climPeriod)

per = 'rcp'
mhws_rcp45 = {}
mhws_rcp85 = {}
clim_rcp45 = {}
clim_rcp85 = {}
which_pre2020 = {}
for model in models:
    mhws_rcp45[model] = {}
    mhws_rcp85[model] = {}
    clim_rcp45[model] = {}
    clim_rcp85[model] = {}
    which_pre2020[model] = (year[per][model] <= 2020)
    for iens in range(NENS[per][model][rcp45]):
        if iens > NENS['hist'][model][hist]:
            iens0 = NENS['hist'][model][hist]
        else:
            iens0 = iens*1
        mhws_rcp45[model][iens], clim_rcp45[model][iens] = mhw.detect(t[per][model][which_pre2020[model]], T_ts[per][model][which_pre2020[model],rcp45,iens], climatologyPeriod=climPeriod, alternateClimatology=[t['hist'][model], T_ts['hist'][model][:,hist,iens0]])
    for iens in range(NENS[per][model][rcp85]):
        if iens > NENS['hist'][model][hist]-1:
            iens0 = NENS['hist'][model][hist]-1
        else:
            iens0 = iens*1
        mhws_rcp85[model][iens], clim_rcp85[model][iens] = mhw.detect(t[per][model][which_pre2020[model]], T_ts[per][model][which_pre2020[model],rcp85,iens], climatologyPeriod=climPeriod, alternateClimatology=[t['hist'][model], T_ts['hist'][model][:,hist,iens0]])

# Model selection
# If KS_p value is large, then we can't reject that the model ensemble and obs may come from the same distribution
# I.e. a p-value of 0.02 means that at the 2% significance level the obs and model come from different distributions!
KS_stat = {}
KS_p = {}
#for key in ['intensity_max', 'duration']:
for key in ['intensity_mean', 'duration']:
    KS_stat[key] = {}
    KS_p[key] = {}
    # Historical run, climatology over the full period but pdf of events over obs period
    for model in models:
        KS_stat[key][model] = {}
        KS_p[key][model] = {}
        for iens in range(NENS['hist'][model][hist]):
            KS_stat[key][model][iens] = {}
            KS_p[key][model][iens] = {}
            which = np.array(mhws_hist[model][iens]['time_end']) >= date(1982,1,1).toordinal()
            KS_stat[key][model][iens], KS_p[key][model][iens] = sp.stats.ks_2samp(np.array(mhws_hist[model][iens][key])[which], mhws_obs[key])

# Generate single list of all model-ensembles, to be used for resampling
MODELS_ENS = {}
NENS_tot = {}
for per in periods:
    MODELS_ENS[per] = {}
    NENS_tot[per] = {}
    for run in range(NRUNS[per][model]):
        MODELS_ENS[per][run] = []
        for model in models:
            for ens in range(NENS[per][model][run]):
                MODELS_ENS[per][run].append([model, ens])
        # Total number of model-ensembles
        NENS_tot[per][run] = len(MODELS_ENS[per][run])

# Parameters for bootstrap resampling
NBOOT = 10000 # Total number of bootstrap resamples
kBOOT = NENS_tot['hist'][0]/2 # Subset of ensembles to sample out
# n = 13 # total number of ensembles
# k  = 7
# n = NENS_tot['hist'][0]
# k = NENS_tot['hist'][0]/2 
# 1.*np.math.factorial(n) / (np.math.factorial(k) * np.math.factorial(n-k)) # number of recombinations

pENS = {}
for per in periods:
    pENS[per] = {}
    for run in range(NRUNS[per][model]):
        # Weight ensemble members according to inv. of # ens per model
        pENS[per][run] = np.zeros(NENS_tot[per][run])
        for ens in range(NENS_tot[per][run]):
            model = MODELS_ENS[per][run][ens][0]
            pENS[per][run][ens] = 1./NENS[per][model][run]
        pENS[per][run] /= pENS[per][run].sum()

#np.random.choice(NENS_tot['rcp'][rcp85], size=kBOOT)
# equiv to: np.floor(NENS_tot['rcp'][rcp85]*np.random.rand(kBOOT)).astype(int)

# Generate pdfs of MHW properties
pdf = {}
pdf_boot = {}
x = {}
#x['intensity_mean'] = np.arange(0,3.,0.01)
x['intensity_max'] = np.arange(0,5.,0.01)
x['duration'] = np.arange(0,500,0.25)

for key in x.keys():
    pdf[key] = {}
    pdf_boot[key] = {}
    # Observations
    pdf[key]['obs'] = sp.stats.gaussian_kde(mhws_obs[key]).evaluate(x[key])
    ## Historical run, over the observations period
    #tmp = []
    #for model in models:
    #    for iens in range(NENS['hist'][model][hist]):
    #        tmp.extend(mhws_hist_recent[model][iens][key])
    #pdf[key]['hist_obs'] = sp.stats.gaussian_kde(tmp).evaluate(x[key])
    # Historical run, over the full period
    tmp = []
    for model in models:
        for iens in range(NENS['hist'][model][hist]):
            tmp.extend(mhws_hist[model][iens][key])
    pdf[key]['hist'] = sp.stats.gaussian_kde(tmp).evaluate(x[key])
    # Historical-Nat run, over the full period
    tmp = []
    for model in models:
        for iens in range(NENS['hist'][model][histNat]):
            tmp.extend(mhws_histNat[model][iens][key])
    pdf[key]['histNat'] = sp.stats.gaussian_kde(tmp).evaluate(x[key])
    # Bootstrap
    pdf_boot[key]['histNat'] = np.zeros((len(x[key]), NBOOT))
    for i in range(NBOOT):
        print key, 'histNat', i, NBOOT
        which_modEns = np.random.choice(NENS_tot['hist'][histNat], size=kBOOT).tolist()
        #np.floor(NENS_tot['hist'][histNat]*np.random.rand(kBOOT)).astype(int).tolist()
        tmp = []
        for ii in range(kBOOT):
            model = MODELS_ENS['hist'][histNat][which_modEns[ii]][0]
            iens = MODELS_ENS['hist'][histNat][which_modEns[ii]][1]
            tmp.extend(mhws_histNat[model][iens][key])
        pdf_boot[key]['histNat'][:,i] = sp.stats.gaussian_kde(tmp).evaluate(x[key])
    # Historical run, climatology over the full period but pdf of events over obs period
    tmp = []
    for model in models:
        for iens in range(NENS['hist'][model][hist]):
            which = np.array(mhws_hist[model][iens]['time_end']) >= date(1982,1,1).toordinal()
            tmp.extend(np.array(mhws_hist[model][iens][key])[which].tolist())
    pdf[key]['hist_recent'] = sp.stats.gaussian_kde(tmp).evaluate(x[key])
    # Bootstrap
    pdf_boot[key]['hist_recent'] = np.zeros((len(x[key]), NBOOT))
    for i in range(NBOOT):
        print key, 'hist_recent', i, NBOOT
        #which_modEns = np.floor(NENS_tot['hist'][hist]*np.random.rand(kBOOT)).astype(int).tolist()
        which_modEns = np.random.choice(NENS_tot['hist'][hist], size=kBOOT).tolist()
        tmp = []
        for ii in range(kBOOT):
            model = MODELS_ENS['hist'][hist][which_modEns[ii]][0]
            iens = MODELS_ENS['hist'][hist][which_modEns[ii]][1]
            which = np.array(mhws_hist[model][iens]['time_end']) >= date(1982,1,1).toordinal()
            tmp.extend(np.array(mhws_hist[model][iens][key])[which].tolist())
        pdf_boot[key]['hist_recent'][:,i] = sp.stats.gaussian_kde(tmp).evaluate(x[key])
    # Historical run, climatology over the full period but pdf of events over 1950-1980
    tmp = []
    for model in models:
        for iens in range(NENS['hist'][model][hist]):
            which = (np.array(mhws_hist[model][iens]['time_end']) >= date(1950,1,1).toordinal()) * (np.array(mhws_hist[model][iens]['time_end']) <= date(1979,12,31).toordinal())
            tmp.extend(np.array(mhws_hist[model][iens][key])[which].tolist())
    pdf[key]['hist_old'] = sp.stats.gaussian_kde(tmp).evaluate(x[key])
    # Bootstrap
    pdf_boot[key]['hist_old'] = np.zeros((len(x[key]), NBOOT))
    for i in range(NBOOT):
        print key, 'hist_old', i, NBOOT
        #which_modEns = np.floor(NENS_tot['hist'][hist]*np.random.rand(kBOOT)).astype(int).tolist()
        which_modEns = np.random.choice(NENS_tot['hist'][hist], size=kBOOT).tolist()
        tmp = []
        for ii in range(kBOOT):
            model = MODELS_ENS['hist'][hist][which_modEns[ii]][0]
            iens = MODELS_ENS['hist'][hist][which_modEns[ii]][1]
            which = (np.array(mhws_hist[model][iens]['time_end']) >= date(1950,1,1).toordinal()) * (np.array(mhws_hist[model][iens]['time_end']) <= date(1979,12,31).toordinal())
            tmp.extend(np.array(mhws_hist[model][iens][key])[which].tolist())
        pdf_boot[key]['hist_old'][:,i] = sp.stats.gaussian_kde(tmp).evaluate(x[key])
    # RCP4.5 run, climatology over the full historical period but pdf of events over recent period
    #tmp = []
    #for model in models:
    #    for iens in range(NENS['rcp'][model]):
    #        which = np.array(mhws_rcp45[model][iens]['time_start']) < date(2020,1,1).toordinal()
    #        tmp.extend(np.array(mhws_rcp45[model][iens][key])[which].tolist())
    #pdf[key]['rcp45'] = sp.stats.gaussian_kde(tmp).evaluate(x[key])
    # Bootstrap
    #pdf_boot[key]['rcp45'] = np.zeros((len(x[key]), NBOOT))
    #for i in range(NBOOT):
    #    print key, 'rcp45', i, NBOOT
    #    which_modEns = np.floor(NENS_tot['rcp']*np.random.rand(kBOOT)).astype(int).tolist()
    #    tmp = []
    #    for ii in range(kBOOT):
    #        model = MODELS_ENS['rcp'][which_modEns[ii]][0]
    #        iens = MODELS_ENS['rcp'][which_modEns[ii]][1]
    #        which = np.array(mhws_rcp45[model][iens]['time_start']) < date(2020,1,1).toordinal()
    #        tmp.extend(np.array(mhws_rcp45[model][iens][key])[which].tolist())
    #    pdf_boot[key]['rcp45'][:,i] = sp.stats.gaussian_kde(tmp).evaluate(x[key])
    # RCP8.5 run, climatology over the full historical period but pdf of events over recent period
    tmp = []
    for model in models:
        for iens in range(NENS['rcp'][model][rcp85]):
            which = np.array(mhws_rcp85[model][iens]['time_start']) < date(2020,1,1).toordinal()
            tmp.extend(np.array(mhws_rcp85[model][iens][key])[which].tolist())
    pdf[key]['rcp85'] = sp.stats.gaussian_kde(tmp).evaluate(x[key])
    # Bootstrap
    pdf_boot[key]['rcp85'] = np.zeros((len(x[key]), NBOOT))
    for i in range(NBOOT):
        print key, 'rcp85', i, NBOOT
        #which_modEns = np.floor(NENS_tot['rcp'][rcp85]*np.random.rand(kBOOT)).astype(int).tolist()
        which_modEns = np.random.choice(NENS_tot['rcp'][rcp85], size=kBOOT).tolist()
        tmp = []
        for ii in range(kBOOT):
            model = MODELS_ENS['rcp'][rcp85][which_modEns[ii]][0]
            iens = MODELS_ENS['rcp'][rcp85][which_modEns[ii]][1]
            which = np.array(mhws_rcp85[model][iens]['time_start']) < date(2020,1,1).toordinal()
            tmp.extend(np.array(mhws_rcp85[model][iens][key])[which].tolist())
        pdf_boot[key]['rcp85'][:,i] = sp.stats.gaussian_kde(tmp).evaluate(x[key])
    # RCP8.5 run, climatology over the full historical period but pdf of events over future period
    tmp = []
    for model in models:
        for iens in range(NENS['rcp'][model][rcp85]):
            which = (np.array(mhws_rcp85[model][iens]['time_start']) >= date(2020,1,1).toordinal()) * (np.array(mhws_rcp85[model][iens]['time_start']) <= date(2039,12,31).toordinal())
            tmp.extend(np.array(mhws_rcp85[model][iens][key])[which].tolist())
    pdf[key]['rcp85_future'] = sp.stats.gaussian_kde(tmp).evaluate(x[key])
    # Bootstrap
    pdf_boot[key]['rcp85_future'] = np.zeros((len(x[key]), NBOOT))
    for i in range(NBOOT):
        print key, 'rcp85_future', i, NBOOT
        #which_modEns = np.floor(NENS_tot['rcp'][rcp85]*np.random.rand(kBOOT)).astype(int).tolist()
        which_modEns = np.random.choice(NENS_tot['rcp'][rcp85], size=kBOOT).tolist()
        tmp = []
        for ii in range(kBOOT):
            model = MODELS_ENS['rcp'][rcp85][which_modEns[ii]][0]
            iens = MODELS_ENS['rcp'][rcp85][which_modEns[ii]][1]
            which = (np.array(mhws_rcp85[model][iens]['time_start']) >= date(2020,1,1).toordinal()) * (np.array(mhws_rcp85[model][iens]['time_start']) <= date(2039,12,31).toordinal())
            tmp.extend(np.array(mhws_rcp85[model][iens][key])[which].tolist())
        pdf_boot[key]['rcp85_future'][:,i] = sp.stats.gaussian_kde(tmp).evaluate(x[key])

outfile = '/data/MHWs/Tasmania_2015_2016/CMIP5/' + whichBox + '/mhw_CMIP5_analyse_pdfs_' + str(climPeriod[0]) + \
          '_' + str(climPeriod[1])
# np.savez(outfile, pdf=pdf, pdf_boot=pdf_boot)
data = np.load(outfile + '.npz')
pdf = data['pdf'].item()
pdf_boot = data['pdf_boot'].item()


#
# Fraction of attributable risk (FAR)
#

# Calculate FAR from pdfs
FAR = {}
FAR_boot = {}
for run in ['hist_old', 'hist_recent', 'rcp85', 'rcp85_future']:
    FAR[run] = {}
    FAR_boot[run] = {}
    for key in pdf.keys():
        FAR[run][key] = 1 - (np.cumsum(pdf[key]['histNat'][::-1])[::-1] / np.sum(pdf[key]['histNat']))/(np.cumsum(pdf[key][run][::-1])[::-1] / np.sum(pdf[key][run]))
        FAR_boot[run][key] = 1 - (np.cumsum(pdf_boot[key]['histNat'][::-1,:], axis=0)[::-1,:] / np.sum(pdf_boot[key]['histNat'], axis=0))/(np.cumsum(pdf_boot[key][run][::-1,:], axis=0)[::-1,:] / np.sum(pdf_boot[key][run], axis=0))

# Get FAR values for particular critical event thresholds
#duration0 = 50
#intensity0 = 1.30

#duration0 = 90
#intensity0 = 1.8

#duration0 = 66
#intensity0 = 1.63

#duration0 = 329
#intensity0 = 2.90

#duration0 = np.sort(mhws_obs['duration'])[-2]
#intensity0 = np.sort(mhws_obs['intensity_max'])[-2]

duration0 = np.sort(mhws_obs_2016['duration'])[-2]
intensity0 = np.sort(mhws_obs_2016['intensity_max'])[-2]

# When did they occur?
tt0 = np.where(np.array(mhws_obs_2016['duration'])==duration0)[0][0]
mhws_obs_2016['date_start'][tt0]
mhws_obs_2016['date_end'][tt0]
tt0 = np.where(np.array(mhws_obs_2016['intensity_max'])==intensity0)[0][0]
mhws_obs_2016['date_peak'][tt0]

#FAR_ev = {}
FAR_ev_boot = {}
#FAR_ev['intensity_mean'] = np.interp(intensity0, x['intensity_mean'], FAR['intensity_mean'])
#FAR_ev['duration'] = np.interp(duration0, x['duration'], FAR['duration'])
for run in FAR_boot.keys():
    FAR_ev_boot[run] = {}
    FAR_ev_boot[run]['intensity_max'] = np.zeros(NBOOT)
    FAR_ev_boot[run]['duration'] = np.zeros(NBOOT)
    for i in range(NBOOT):
        FAR_ev_boot[run]['intensity_max'][i] = np.interp(intensity0, x['intensity_max'], FAR_boot[run]['intensity_max'][:,i])
        FAR_ev_boot[run]['duration'][i] = np.interp(duration0, x['duration'], FAR_boot[run]['duration'][:,i])

#for key in FAR_ev.keys():
#    print key, FAR_ev[key], 1./(1 - FAR_ev[key])

for run in FAR_ev_boot.keys():
    for key in FAR_ev_boot[run].keys():
        # All FAR values:
        #print run, key, np.percentile(FAR_ev_boot[run][key], [5, 50, 95]), 1./(1 - np.percentile(FAR_ev_boot[run][key], [5, 50, 95])), 1 - 1.*(FAR_ev_boot[run][key] > 0).sum() / len(FAR_ev_boot[run][key])
        # Only FAR values < 1:
        #tmp = FAR_ev_boot[run][key][FAR_ev_boot[run][key]<=0.99]
        tmp = FAR_ev_boot[run][key] #[FAR_ev_boot[run][key]<0.99]
        #print run, key, np.percentile(tmp, [5, 50, 95]), 1./(1 - np.percentile(tmp, [5, 50, 95])), 1 - 1.*(tmp > 0).sum() / len(tmp)
        print run, key, np.percentile(FAR_ev_boot[run][key], [50, 10]), 1./(1 - np.percentile(FAR_ev_boot[run][key], [50, 10])), 1 - 1.*(FAR_ev_boot[run][key] > 0).sum() / len(FAR_ev_boot[run][key])

# Plot it up
#from statsmodels.nonparametric.kde import KDEUnivariate
#
#def kde_statsmodels_u(x, x_grid, bandwidth=0.2, **kwargs):
#    """Univariate Kernel Density Estimation with Statsmodels"""
#    kde = KDEUnivariate(x)
#    kde.fit(bw=bandwidth, **kwargs)
#    return kde.evaluate(x_grid)
#
#x_FAR = np.arange(-10, 2, 0.01)
#hist_FAR = {}
#pdf_FAR = {}
#for run in FAR_ev_boot.keys():
#    hist_FAR[run] = {}
#    pdf_FAR[run] = {}
#    for key in FAR_ev_boot[run].keys():
#        # All FAR values:
#        #hist_FAR[run][key], bins = np.histogram(FAR_ev_boot[run][key], bins=40, range=(-0.5,1.5), density=True)
#        # Only FAR values < 1:
#        #hist_FAR[run][key], bins = np.histogram(FAR_ev_boot[run][key][FAR_ev_boot[run][key]<1.], bins=40, range=(-0.5,1.5), density=True)
#        #pdf_FAR[run][key] = sp.stats.gaussian_kde(FAR_ev_boot[run][key][FAR_ev_boot[run][key]<=0.95], bw_method=1e-15).evaluate(x_FAR)
#        pdf_FAR[run][key] = kde_statsmodels_u(FAR_ev_boot[run][key][FAR_ev_boot[run][key]<0.99], x_FAR, bandwidth=0.05)
#        #pdf_FAR[run][key] = kde_statsmodels_u(FAR_ev_boot[run][key][FAR_ev_boot[run][key]<1.], x_FAR, bandwidth=0.05)
#
##hist_FAR['intensity_mean'], bins = np.histogram(FAR_ev_boot['intensity_mean'], bins=40, range=(-0.5,1.5))
##bins = bins[0:-1] + np.diff(bins)

# Figure for paper

fig = plt.figure(figsize=(15,10))
plt.clf()
#bins = bins[:-1] #+ np.diff(bins)

key = 'intensity_max'
plt.subplot(3,2,1)
plt.plot(x[key], pdf[key]['obs'], 'k-', linewidth=2)
plt.plot(x[key], pdf[key]['histNat'], 'b-')
plt.plot(x[key], pdf[key]['hist_recent'], 'k-')
plt.plot(x[key], pdf[key]['hist_old'], 'k--')
plt.plot(x[key], pdf[key]['rcp85'], 'r-')
plt.plot(x[key], pdf[key]['rcp85_future'], 'r--')
plt.plot(intensity0, -0.01, '^', markerfacecolor='k', markeredgecolor='k', markersize=10, clip_on=False, zorder=100)
plt.plot(np.sort(mhws_obs_2016[key])[-1], -0.01, '^', markerfacecolor=(1,0.3,0.3), markeredgecolor='k', markersize=10, clip_on=False, zorder=100)
plt.xlim(0, 4)
plt.ylim(0, 1.2)
#plt.grid()
plt.title('(a)')
plt.ylabel('Probability density')
plt.xlabel('MHW max intensity [deg. C]')

key = 'duration'
plt.subplot(3,2,2)
plt.semilogx(x[key], pdf[key]['obs'], 'k-', linewidth=2)
plt.semilogx(x[key], pdf[key]['histNat'], 'b-')
plt.semilogx(x[key], pdf[key]['hist_recent'], 'k-')
plt.semilogx(x[key], pdf[key]['hist_old'], 'k--')
plt.semilogx(x[key], pdf[key]['rcp85'], 'r-')
plt.semilogx(x[key], pdf[key]['rcp85_future'], 'r--')
plt.semilogx(duration0, -0.0005, '^', markerfacecolor='k', markeredgecolor='k', markersize=10, clip_on=False, zorder=100)
plt.semilogx(np.sort(mhws_obs_2016[key])[-1], -0.0005, '^', markerfacecolor=(1,0.3,0.3), markeredgecolor='k', markersize=10, clip_on=False, zorder=100)
plt.xlim(10, 1000)
plt.ylim(0, 0.035)
#plt.grid()
plt.legend(['Observations (1982-2005)', 'histNat (1850-2005)', 'hist (1982-2005)', 'hist (1950-1980)', 'rcp85 (2006-2020)', 'rcp85 (2020-2040)'], fontsize=10)
plt.title('(b)')
plt.xlabel('MHW duration [days]')

nBins = 75
key = 'intensity_max'
ax1a = plt.subplot(3,2,5)
ax1b = ax1a.twinx()
ax1c = ax1a.twinx()
ax1d = ax1a.twinx()
ax2 = ax1a.twiny()
na, bins, patches = ax1a.hist(FAR_ev_boot['hist_old'][key], bins=nBins, range=(-0.5,0.99), normed=True, histtype='step', color='0.5', linewidth=1)
    #ax1a.plot(np.percentile(FAR_ev_boot['hist_old'][key],10), np.interp(np.percentile(FAR_ev_boot['hist_old'][key],10), bins[:-1]+0.5*np.diff(bins)[0], na), 'o', markerfacecolor='w', markeredgewidth=2, markersize=7, markeredgecolor='0.5')
    #ax1a.plot(np.percentile(FAR_ev_boot['hist_old'][key],50), np.interp(np.percentile(FAR_ev_boot['hist_old'][key],50), bins[:-1]+0.5*np.diff(bins)[0], na), 'o', markerfacecolor='0.5', markeredgewidth=2, markersize=7, markeredgecolor='0.5')
nb, bins, patches = ax1b.hist(FAR_ev_boot['hist_recent'][key], bins=nBins, range=(-0.5,0.99), normed=True, histtype='step', color='k', linewidth=2)
    ax1b.plot(np.percentile(FAR_ev_boot['hist_recent'][key],10), np.interp(np.percentile(FAR_ev_boot['hist_recent'][key],10), bins[:-1]+0.5*np.diff(bins)[0], nb), 'o', markerfacecolor='w', markeredgewidth=2, markersize=7, markeredgecolor='k')
    ax1b.plot(np.percentile(FAR_ev_boot['hist_recent'][key],50), np.interp(np.percentile(FAR_ev_boot['hist_recent'][key],50), bins[:-1]+0.5*np.diff(bins)[0], nb), 'o', markerfacecolor='k', markeredgewidth=2, markersize=7, markeredgecolor='k')
nc, bins, patches = ax1c.hist(FAR_ev_boot['rcp85'][key], bins=nBins, range=(-0.5,0.99), normed=True, histtype='step', color='r', linewidth=2)
    ax1c.plot(np.percentile(FAR_ev_boot['rcp85'][key],10), np.interp(np.percentile(FAR_ev_boot['rcp85'][key],10), bins[:-1]+0.5*np.diff(bins)[0], nc), 'o', markerfacecolor='w', markeredgewidth=2, markersize=7, markeredgecolor='r')
    ax1c.plot(np.percentile(FAR_ev_boot['rcp85'][key],50), np.interp(np.percentile(FAR_ev_boot['rcp85'][key],50), bins[:-1]+0.5*np.diff(bins)[0], nc), 'o', markerfacecolor='r', markeredgewidth=2, markersize=7, markeredgecolor='r')
nd, bins, patches = ax1d.hist(FAR_ev_boot['rcp85_future'][key], bins=nBins, range=(-0.5,0.99), normed=True, histtype='step', color=(1,0.3,0.3), linewidth=1)
    #ax1d.plot(np.percentile(FAR_ev_boot['rcp85_future'][key],10), np.interp(np.percentile(FAR_ev_boot['rcp85_future'][key],10), bins[:-1]+0.5*np.diff(bins)[0], nd), 'o', markerfacecolor='w', markeredgewidth=2, markersize=7, markeredgecolor=(1,0.3,0.3))
    #ax1d.plot(np.percentile(FAR_ev_boot['rcp85_future'][key],50), np.interp(np.percentile(FAR_ev_boot['rcp85_future'][key],50), bins[:-1]+0.5*np.diff(bins)[0], nd), 'o', markerfacecolor=(1,0.3,0.3), markeredgewidth=2, markersize=7, markeredgecolor=(1,0.3,0.3))
ax1a.set_xlim(0, 1)
ax1b.set_xlim(0, 1)
ax1c.set_xlim(0, 1)
ax1d.set_xlim(0, 1)
ax2.set_xlim(0, 1)
ax1a.set_ylim(0, 1.1*na.max())
ax1b.set_ylim(0, 1.1*nb.max())
ax1c.set_ylim(0, 1.1*nc.max())
ax1d.set_ylim(0, 1.1*nd.max())
#ax2.set_ylim(0, 1)
ax2.set_xticks(np.array([0, 0.5, 0.75, 0.90]))
ax2.set_xticklabels(['(c) 1x (None)', '2x', '4x', '10x'])
ax1b.set_yticks([])
ax1c.set_yticks([])
ax1d.set_yticks([])
#ax1a.grid()
#ax2.grid()
ax1a.set_ylabel('Normalised frequency')
ax1a.set_xlabel('Fraction of Attributable Risk (FAR)')

nBins = 500
key = 'duration'
ax1a = plt.subplot(3,2,6)
ax1b = ax1a.twinx()
ax1c = ax1a.twinx()
ax1d = ax1a.twinx()
ax2 = ax1a.twiny()
na, bins, patches = ax1a.hist(FAR_ev_boot['hist_old'][key], bins=nBins, range=(-0.5,0.999), normed=True, histtype='step', color='0.5', linewidth=1)
    #ax1a.plot(np.percentile(FAR_ev_boot['hist_old'][key],10), np.interp(np.percentile(FAR_ev_boot['hist_old'][key],10), bins[:-1]+0.5*np.diff(bins)[0], na), 'o', markerfacecolor='w', markeredgewidth=2, markersize=7, markeredgecolor='0.5')
    #ax1a.plot(np.percentile(FAR_ev_boot['hist_old'][key],50), np.interp(np.percentile(FAR_ev_boot['hist_old'][key],50), bins[:-1]+0.5*np.diff(bins)[0], na), 'o', markerfacecolor='0.5', markeredgewidth=2, markersize=7, markeredgecolor='0.5')
nb, bins, patches = ax1b.hist(FAR_ev_boot['hist_recent'][key], bins=nBins, range=(-0.5,0.999), normed=True, histtype='step', color='k', linewidth=2)
    ax1b.plot(np.percentile(FAR_ev_boot['hist_recent'][key],10), np.interp(np.percentile(FAR_ev_boot['hist_recent'][key],10), bins[:-1]+0.5*np.diff(bins)[0], nb), 'o', markerfacecolor='w', markeredgewidth=2, markersize=7, markeredgecolor='k')
    ax1b.plot(np.percentile(FAR_ev_boot['hist_recent'][key],50), np.interp(np.percentile(FAR_ev_boot['hist_recent'][key],50), bins[:-1]+0.5*np.diff(bins)[0], nb), 'o', markerfacecolor='k', markeredgewidth=2, markersize=7, markeredgecolor='k')
nc, bins, patches = ax1c.hist(FAR_ev_boot['rcp85'][key], bins=nBins, range=(-0.5,0.999), normed=True, histtype='step', color='r', linewidth=2)
    ax1c.plot(np.percentile(FAR_ev_boot['rcp85'][key],10), np.interp(np.percentile(FAR_ev_boot['rcp85'][key],10), bins[:-1]+0.5*np.diff(bins)[0], nc), 'o', markerfacecolor='w', markeredgewidth=2, markersize=7, markeredgecolor='r')
    ax1c.plot(np.percentile(FAR_ev_boot['rcp85'][key],50), np.interp(np.percentile(FAR_ev_boot['rcp85'][key],50), bins[:-1]+0.5*np.diff(bins)[0], nc), 'o', markerfacecolor='r', markeredgewidth=2, markersize=7, markeredgecolor='r')
nd, bins, patches = ax1d.hist(FAR_ev_boot['rcp85_future'][key], bins=nBins, range=(-0.5,0.999), normed=True, histtype='step', color=(1,0.3,0.3), linewidth=1)
    #ax1d.plot(np.percentile(FAR_ev_boot['rcp85_future'][key],10), np.interp(np.percentile(FAR_ev_boot['rcp85_future'][key],10), bins[:-1]+0.5*np.diff(bins)[0], nd), 'o', markerfacecolor='w', markeredgewidth=2, markersize=7, markeredgecolor=(1,0.3,0.3))
    #ax1d.plot(np.percentile(FAR_ev_boot['rcp85_future'][key],50), np.interp(np.percentile(FAR_ev_boot['rcp85_future'][key],50), bins[:-1]+0.5*np.diff(bins)[0], nd), 'o', markerfacecolor=(1,0.3,0.3), markeredgewidth=2, markersize=7, markeredgecolor=(1,0.3,0.3))
ax1a.set_xlim(0.9, 1)
ax1b.set_xlim(0.9, 1)
ax1c.set_xlim(0.9, 1)
ax1d.set_xlim(0.9, 1)
ax2.set_xlim(0.9, 1)
ax1a.set_ylim(0, 1.1*na.max())
ax1b.set_ylim(0, 1.1*nb.max())
ax1c.set_ylim(0, 1.1*nc.max())
ax1d.set_ylim(0, 1.1*nd.max())
ax2.set_xticks(np.array([0.9, 0.95, 0.98, 0.99]))
ax2.set_xticklabels(['(d) 10x', '20x', '50x', '100x'])
ax1b.set_yticks([])
ax1c.set_yticks([])
ax1d.set_yticks([])
#ax1a.grid()
#ax2.grid()
ax1a.set_xlabel('Fraction of Attributable Risk (FAR)')
#ax1a.legend(['hist (1950-1980)'], loc='upper left', fontsize=10, framealpha=0.8)
#ax1b.legend(['hist (1982-2005)'], loc='upper left', fontsize=10, framealpha=0.8)
#ax1c.legend(['rcp85 (2006-2020)'], loc='upper left', fontsize=10, framealpha=0.8)
#ax1d.legend(['rcp85 (2020-2040)'], loc='upper left', fontsize=10, framealpha=0.8)

# plt.savefig('../../documents/14_Tasmania_2015_2016/figures/DA_SEAus_MHWs_AllModels.png', bbox_inches='tight', pad_inches=0.5, dpi=300)



# plt.savefig('../../documents/14_Tasmania_2015_2016/figures/DA_SEAus_MHWs_AllModels_woKing_wo' + sys.argv[-1] + '.png', bbox_inches='tight', pad_inches=0.5, dpi=150)
# plt.savefig('../../documents/14_Tasmania_2015_2016/figures/DA_SEAus_MHWs_AllModels_woKing_only' + sys.argv[-1] + '.png', bbox_inches='tight', pad_inches=0.5, dpi=150)

