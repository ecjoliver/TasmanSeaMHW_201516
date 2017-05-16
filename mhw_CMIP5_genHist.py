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
from netCDF4 import Dataset
import os
import ecoliver as ecj
import deseason as ds

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

import marineHeatWaves as mhw

#
# Basic variables
#

# Location of interest
lon = [147, 157]
lat = [-46, -39]
# Old box
#lon = [147, 155]
#lat = [-45, -37]

# Info on runs and ensemble members
pathroot = '/mnt/insect/'

# Model
model = 'ACCESS1-3'
header = [[pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/historical/day/ocean/day/r1i1p1/files/tos_20121112/',
           pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/historical/day/ocean/day/r2i1p1/files/tos_20121112/',
           pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/historical/day/ocean/day/r3i1p1/files/tos_20121112/'],
          [pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/historicalNat/day/ocean/day/r1i1p1/files/tos_20130906/',
           pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/historicalNat/day/ocean/day/r2i1p1/files/tos_20130912/',
           pathroot+'data/CMIP5/CSIRO-BOM/ACCESS1-3/historicalNat/day/ocean/day/r3i1p1/files/tos_20140501/']]
# Time and date vectors
Ly = 365.25
t, tmp, T, year, month, day, doy = ecj.timevector([1850,1,1], [2005,12,31])

# Model
model = 'CSIRO-Mk3-6-0'
header = [[pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r1i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r2i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r3i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r4i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r5i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r6i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r7i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r8i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r9i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historical/day/ocean/day/r10i1p1/v20111222/tos/'],
          [pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r1i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r2i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r3i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r4i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r5i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r6i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r7i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r8i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r9i1p1/v20111222/tos/',
           pathroot+'data/CMIP5/CSIRO-QCCCE/CSIRO-Mk3-6-0/historicalNat/day/ocean/day/r10i1p1/v20111222/tos/']]
# Time and date vectors
Ly = 365
t, tmp, T, year, month, day, doy = ecj.timevector([1850,1,1], [2012,12,31])
feb29s = ~((month==2) * (day==29))
t = t[feb29s]
year = year[feb29s]
month = month[feb29s]
day = day[feb29s]
doy = doy[feb29s]
T = len(t)

# Model
model = 'CNRM-CM5'
header = [[pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r1i1p1/v20111213/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r2i1p1/v20111114/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r3i1p1/v20111114/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r4i1p1/v20111115/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r5i1p1/v20111115/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r6i1p1/v20111115/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r7i1p1/v20111125/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r8i1p1/v20111117/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r9i1p1/v20111115/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historical/day/ocean/day/r10i1p1/v20111122/tos/'],
          [pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historicalNat/day/ocean/day/r1i1p1/v20120702/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historicalNat/day/ocean/day/r2i1p1/v20111213/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historicalNat/day/ocean/day/r3i1p1/v20111213/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historicalNat/day/ocean/day/r4i1p1/v20111213/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historicalNat/day/ocean/day/r5i1p1/v20111213/tos/',
           pathroot+'data/CMIP5/CNRM-CERFACS/CNRM-CM5/historicalNat/day/ocean/day/r8i1p1/v20111213/tos/']]
# Time and date vectors
Ly = 365.25
t, tmp, T, year, month, day, doy = ecj.timevector([1850,1,1], [2012,12,31])

# Model
model = 'HadGEM2-ES'
header = [[pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historical/day/ocean/day/r1i1p1/v20110131/tos/',
           pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historical/day/ocean/day/r2i1p1/v20110418/tos/',
           pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historical/day/ocean/day/r3i1p1/v20110418/tos/',
           pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historical/day/ocean/day/r4i1p1/v20110418/tos/',],
          [pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historicalNat/day/ocean/day/r1i1p1/v20110728/tos/',
           pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historicalNat/day/ocean/day/r2i1p1/v20110609/tos/',
           pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historicalNat/day/ocean/day/r3i1p1/v20110609/tos/',
           pathroot+'data/CMIP5/MOHC/HadGEM2-ES/historicalNat/day/ocean/day/r4i1p1/v20110609/tos/']]
# Time and date vectors
Ly = 360.
t, tmp, T, year, month, day, doy = ecj.timevector([1859,12,1], [2019,12,30])
feb29s = ~((month==2) * (day==29))
t = t[feb29s]
year = year[feb29s]
month = month[feb29s]
day = day[feb29s]
doy = doy[feb29s]
y360 = ~((day > 30) * (month > 3))
t = t[y360]
year = year[y360]
month = month[y360]
day = day[y360]
doy = doy[y360]
T = len(t)

# Model
model = 'IPSL-CM5A-LR'
header = [[pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historical/day/ocean/day/r1i1p1/v20111010/tos/',
           pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historical/day/ocean/day/r2i1p1/v20111010/tos/',
           pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historical/day/ocean/day/r3i1p1/v20111010/tos/',
           pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historical/day/ocean/day/r4i1p1/v20111010/tos/',
           pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historical/day/ocean/day/r5i1p1/v20111119/tos/',
           pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historical/day/ocean/day/r6i1p1/v20120526/tos/'],
          [pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historicalNat/day/ocean/day/r1i1p1/v20130506/tos/',
           pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historicalNat/day/ocean/day/r2i1p1/v20130506/tos/',
           pathroot+'data/CMIP5/IPSL/IPSL-CM5A-LR/historicalNat/day/ocean/day/r3i1p1/v20130506/tos/']]
# Time and date vectors
Ly = 365
t, tmp, T, year, month, day, doy = ecj.timevector([1850,1,1], [2012,12,31])
feb29s = ~((month==2) * (day==29))
t = t[feb29s]
year = year[feb29s]
month = month[feb29s]
day = day[feb29s]
doy = doy[feb29s]
T = len(t)

# Model
model = 'IPSL-CM5A-MR'
header = [[pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/historical/day/ocean/day/r1i1p1/v20111119/tos/',
           pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/historical/day/ocean/day/r2i1p1/v20120430/tos/',
           pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/historical/day/ocean/day/r3i1p1/v20120804/tos/'],
          [pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/historicalNat/day/ocean/day/r1i1p1/v20120804/tos/',
           pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/historicalNat/day/ocean/day/r2i1p1/v20120804/tos/',
           pathroot+'data/CMIP5/IPSL/IPSL-CM5A-MR/historicalNat/day/ocean/day/r3i1p1/v20120804/tos/']]
# Time and date vectors
Ly = 365
t, tmp, T, year, month, day, doy = ecj.timevector([1850,1,1], [2012,12,31])
feb29s = ~((month==2) * (day==29))
t = t[feb29s]
year = year[feb29s]
month = month[feb29s]
day = day[feb29s]
doy = doy[feb29s]
T = len(t)

# Model
model = 'CanESM2'
header = [[pathroot+'data/CMIP5/CCCma/CanESM2/historical/day/ocean/day/r1i1p1/tos/1/'],
          [pathroot+'data/CMIP5/CCCma/CanESM2/historicalNat/day/ocean/day/r1i1p1/tos/1/',
           pathroot+'data/CMIP5/CCCma/CanESM2/historicalNat/day/ocean/day/r3i1p1/tos/1/',
           pathroot+'data/CMIP5/CCCma/CanESM2/historicalNat/day/ocean/day/r5i1p1/tos/1/']]
# Time and date vectors
Ly = 365
t, tmp, T, year, month, day, doy = ecj.timevector([1850,1,1], [2012,12,31])
feb29s = ~((month==2) * (day==29))
t = t[feb29s]
year = year[feb29s]
month = month[feb29s]
day = day[feb29s]
doy = doy[feb29s]
T = len(t)


NRUNS = len(header)
NENS = []
for run in range(NRUNS):
    NENS.append(len(header[run]))

#NENS = len(header[0])
hist = 0
histNat = 1
nyears = len(np.unique(year))

# Temperature array
sst_offset_hist = 0
T_ts_hist = np.NaN*np.zeros([T,NRUNS,np.max(NENS)])

#
# Load data from each run/ensemble member
#

for run in range(NRUNS):
    for ens in range(NENS[run]):
        print 'Run', run, 'Ensemble member', ens
        # Files
        files = []
        os.chdir(header[run][ens])
        for file in os.listdir('.'):
            if file.endswith('.nc'):
                files.append(header[run][ens] + file)
        N_files = len(files)
        file0 = files[0]
        # lat and lons of obs
        fileobj = Dataset(file0, mode='r')
        lon_full = fileobj.variables['lon'][:].astype(float)
        lat_full = fileobj.variables['lat'][:].astype(float)
        fill_value = fileobj.variables['tos']._FillValue.astype(float)
        fileobj.close()
        #   load SST
        t_ts = np.nan*np.ones((T,))
        t0 = 0
        if (model == 'ACCESS1-3') + (model == 'IPSL-CM5A-LR') + (model == 'IPSL-CM5A-MR'):
            i_ll, j_ll = ecj.findxy(lon_full, lat_full, (lon[0], lat[0]))
            i_ur, j_ur = ecj.findxy(lon_full, lat_full, (lon[1], lat[1]))
        elif (model == 'CSIRO-Mk3-6-0') + (model == 'HadGEM2-ES') + (model == 'CanESM2'):
            i_ll = ecj.find_nearest(lon_full, lon[0])[1]
            i_ur = ecj.find_nearest(lon_full, lon[1])[1]
            j_ll = ecj.find_nearest(lat_full, lat[0])[1]
            j_ur = ecj.find_nearest(lat_full, lat[1])[1]
        elif model == 'CNRM-CM5':
            #mask = (lon_full >= lon[0]) * (lon_full <= lon[1]) * (lat_full >= lat[0]) * (lat_full <= lat[1])
            j_ll = np.where((lat_full >= lat[0]) * (lat_full <= lat[1]))[0][0]
            j_ur = np.where((lat_full >= lat[0]) * (lat_full <= lat[1]))[0][-1]
            i_ll = np.where((lon_full[j_ll,:] >= lon[0]) * (lon_full[j_ll,:] <= lon[1]))[0][0]
            i_ur = np.where((lon_full[j_ll,:] >= lon[0]) * (lon_full[j_ll,:] <= lon[1]))[0][-1]
        for file in files:
            fileobj = Dataset(file, mode='r')
            t_file = fileobj.variables['time'][:]
            T_file = len(t_file)
            sst_file_map = fileobj.variables['tos'][:,j_ll:j_ur+1,:][:,:,i_ll:i_ur+1]
            sst_file_map = sst_file_map.data
            sst_file_map[sst_file_map==fill_value] = np.nan
            sst_file = np.nanmean(np.nanmean(sst_file_map, axis=2), axis=1)
            T_ts_hist[t0:(t0+T_file),run,ens] = sst_file
            t_ts[range(t0, t0 + T_file)] = t_file
            t0 = t0 + T_file
            fileobj.close()
        #   Sort based on time-vector and add offset (if required)
        tt = np.argsort(t_ts)
        T_ts_hist[:,run,ens] = T_ts_hist[tt,run,ens] + sst_offset_hist

#
# Save data
#

# Restrict data to end at 2005
#tt_2005 = np.where(t == date(2005, 12, 31).toordinal())[0][0]

outfile = '/data/MHWs/Tasmania_2015_2016/CMIP5/mhw_' + model + '_hist'
np.savez(outfile, lon=lon, lat=lat, T_ts_hist=T_ts_hist, t=t, T=T, year=year, month=month, day=day, doy=doy, Ly=Ly, header=header, NRUNS=NRUNS, NENS=NENS)

