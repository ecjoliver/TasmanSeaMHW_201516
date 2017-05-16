import numpy as np

def deseason_harmonic(dat, K, L):
    '''
    deseasoned_data, season, beta = deseason_harmonic(dat, K, L)

    Subtracts the seasonal cycle (season) from the data (data). Season
    is calculated by fitting K harmonics of the annual cycle (as well as the
    mean) to the data. Assumes on year has L elements (i.e., 365 for daily data,
    73 for pentad data, 52 for weekly data, etc.).
    Outputs the deseasonalized data, the season, and the fitting parameters (beta)

    Handles missing values as np.nan's

    Written by Eric Oliver, Dalhousie University, 2007-2011
    Adapted from original MATLAB script on 28 November 2012
    '''

#   Valid (non-NaN) elements
    valid = ~np.isnan(dat)

#   ensure that mat is a matrix and a column vector
    dat = np.mat(dat)
    if dat.shape[1]!=0:
        dat = dat.T

#   length of time series and generate time vector
    n = len(dat)
    time = np.mat(np.arange(1,n+1)/(1.*L))

#   set up mean and harmonics to fit data
    P = np.mat(np.ones((n,1)))
    for k in range(1,K+1):
        P = np.concatenate((P, np.cos(k*2*np.pi*time.T)), 1)
        P = np.concatenate((P, np.sin(k*2*np.pi*time.T)), 1)

#   Remove seasonal cycle by harmonic regression
    beta = (np.linalg.inv(P[valid,:].T*P[valid,:])*P[valid,:].T)*dat[valid]
    season = P*beta
    dat_ds = dat - season

    return dat_ds, season, beta

def deseason_harmonic_2D(dat, K, L, detrend=False):
    '''
    deseasoned_data, season, beta = deseason_harmonic_2D(dat, K, L)

    Subtracts the seasonal cycle (season) from the data (data). Season
    is calculated by fitting K harmonics of the annual cycle (as well as the
    mean) to the data. Assumes on year has L elements (i.e., 365 for daily data,
    73 for pentad data, 52 for weekly data, etc.).
    Outputs the deseasonalized data, the season, the trend, and the fitting parameters (beta)

    First dimension of dat must be time dimension.
    Does not handle missing values.

    Written by Eric Oliver, Dalhousie University, 2007-2011
    Adapted from original MATLAB script on 28 November 2012
    '''

#   Valid (non-NaN) elements
    valid = ~np.isnan(dat)

#   ensure that mat is a matrix and a column vector
    dat = np.mat(dat)
    #if dat.shape[1]!=0:
    #    dat = dat.T

#   length of time series and generate time vector
    n = dat.shape[0]
    time = np.mat(np.arange(1,n+1)/(1.*L))

#   set up mean and harmonics to fit data
    P = np.mat(np.ones((n,1)))
    for k in range(1,K+1):
        P = np.concatenate((P, np.cos(k*2*np.pi*time.T)), 1)
        P = np.concatenate((P, np.sin(k*2*np.pi*time.T)), 1)
    if detrend:
        P = np.concatenate((P, time.T-time.mean()), 1)

#   Remove seasonal cycle by harmonic regression
    beta = (np.linalg.inv(P.T*P)*P.T)*dat
    season = P[:,0:K*2+1]*beta[0:K*2+1,:]
    if detrend:
        trend = P[:,-1]*beta[-1,:]
    else:
        trend = 0.*season
    dat_ds = dat - season - trend

    return dat_ds, season, trend, beta



