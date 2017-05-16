
clear all

addpath(genpath('../include'));
rehash;

% set-up vars
%header_out = '/mnt/compass/data/MHWs/Tasmania_2015_2016/OMAPS/';
header_out = '/data/MHWs/Tasmania_2015_2016/OMAPS/';
header = 'http://opendap.bom.gov.au:8080/thredds/dodsC/oceanmaps_an_datasets/version_';

% Read in information about versions, and temporal range
% Note data/file_info.txt.{version} must be updated by hand based on the dates of available data on the opendap server
versions = {'2.0', '2.1', '2.2', '2.2.1'};
N_versions = length(versions);
N_ens = 4;
for vers = 1:N_versions
  dat = dlmread(['data/file_info.txt.' versions{vers}]);
  for ens = 1:N_ens,
    time_ens{vers}{ens} = datenum(num2str(dat(ens,2)), 'yyyymmdd'):datenum(num2str(dat(ens,3)), 'yyyymmdd');
    dates_ens{vers}{ens} = datevec(time_ens{vers}{ens});
  end
  % Find latest start time and earliest end time for each version
  t1 = Inf; t2 = -Inf;
  for ens = 1:N_ens,
    t1 = min(t1, time_ens{vers}{ens}(1));
    t2 = max(t2, time_ens{vers}{ens}(end));
  end
  time_vers{vers} = t1:t2;
  dates_vers{vers} = datevec(time_vers{vers});
end

% Full time range for generation of OMAPS forcing data
t1 = datenum([1993 1 1]) + 7151 + 1; t2 = datenum([2016 6 30]); time = t1:t2;
dates = datevec(time); years = unique(dates(:,1))';

% Generate weights for averaging across difference versions and across ensemble members
w = zeros(length(time), N_versions);
for vers = 1:N_versions
%  if vers == 2 % special case for version 2.1, want to skip first 10 overlaps
%    tt = find(time >= time_vers{vers}(1)+10 & time <= time_vers{vers}(end));
%  else
    tt = find(time >= time_vers{vers}(1) & time <= time_vers{vers}(end));
%  end
  w(tt,vers) = 1;
  % find overlap with previous version, and linearly scale weights across overlap zone
  if vers>1
    tt = find(w(:,vers) == 1 & w(:,vers-1) == 1);
    w_linear = linspace(0, 1, length(tt)+2)'; % need +2 for extra 0 and beginning and 1 at end
    w(tt,vers) = w_linear(2:end-1);
    w(tt,vers-1) = 1-w_linear(2:end-1);
  end
  % For each version, find overlap across ensembles and create weights for geenration of ensemble mean
  w_ens{vers} = zeros(length(time), N_ens);
  tt = find(time >= time_vers{vers}(1) & time <= time_vers{vers}(end));
  for ttt = tt;
    for ens = 1:N_ens
      ens_exist = length(find(time_ens{vers}{ens} == time(ttt)));
      if ens_exist
        w_ens{vers}(ttt,ens) = 1;
      end
    end
  end
  w_ens{vers} = w_ens{vers} ./ repmat(sum(w_ens{vers},2), [1 N_ens]);
  w_ens{vers}(isnan(w_ens{vers})) = 0;
end
% Special case for Apr/2016
tt = find((dates(:,1)==2016) & (dates(:,2)==4) & (dates(:,3)==5)); % 5/4/2014
w_ens{4}(tt,:) = [1 0 1 1]/3.; % 5/4/2014
w_ens{4}(tt+1,:) = [1 0 0 1]/2.; % 6/4/2014
w_ens{4}(tt+2,:) = [1 0 0 0]/1.; % 7/4/2014
w_ens{4}(tt+3,:) = [1 0 0 0]/1.; % 8/4/2014
w_ens{4}(tt+4,:) = [1 1 0 0]/2.; % 9/4/2014
w_ens{4}(tt+5,:) = [1 1 1 0]/3.; % 10/4/2014

% Load some basic variables
vers = 1; ens = 1;
file0 = [versions{vers} '/ocean_an0' num2str(ens-1) '_20110427_u.nc'];
lon = ncread([header file0], 'xu_ocean');
lat = ncread([header file0], 'yu_ocean');
depth = ncread([header file0], 'st_ocean');
%load data/OMAPS_BRAN_Offset_eta.mat eta_offset % BRAN3 / OMaps sea level offset

% Domain
i1 = findnearest(143, lon, -1);
i2 = findnearest(161, lon, +1);
j1 = findnearest(-46, lat, -1);
j2 = findnearest(-30, lat, +1);
k1 = 1; % 0 m
k2 = 40; % >1000 m

% Coordinates
x = lon(i1:i2); X = length(x);
y = lat(j1:j2); Y = length(y);
[x, y] = meshgrid(x, y);
x = x'; y = y';
z = depth(k1:k2); Z = length(z);
z = -z; z = flipud(z);
t_full = time; % datenum([1990 1 1]); % SHOC time format

for year = years, for month = [1:12],
  % Initialize variables
  tt = find((dates(:,1)==year) & (dates(:,2)==month));
  if length(tt)==0, continue, end;
  t = t_full(tt); T = length(t);
  u = NaN*zeros(X,Y,Z,T);
  v = NaN*zeros(X,Y,Z,T);
  temp = NaN*zeros(X,Y,Z,T);
  salt = NaN*zeros(X,Y,Z,T);
  eta = NaN*zeros(X,Y,T);
  tau_x = NaN*zeros(X,Y,T);
  tau_y = NaN*zeros(X,Y,T);
  sfc_hflux = NaN*zeros(X,Y,T);
  pme_sbc = NaN*zeros(X,Y,T);
  evap = NaN*zeros(X,Y,T);
  lprec = NaN*zeros(X,Y,T);
  swflx = NaN*zeros(X,Y,T);
  evap_heat = NaN*zeros(X,Y,T);
  sens_heat = NaN*zeros(X,Y,T);
  lw_heat = NaN*zeros(X,Y,T);
  % Loop over days
  tcnt = 0;
  for ttt = tt', tcnt = tcnt + 1; [ttt dates(ttt,1:3)]
    day = dates(ttt,3);
    % Load data, using a try-catch system to keep trying if there are errors in loading data
    read_fail = 1; % assume a priori that the data is not readable
    while read_fail
      try % try to read each of the data files in turn
        U = zeros(X,Y,Z);
        V = zeros(X,Y,Z);
        TEMP = zeros(X,Y,Z);
        SALT = zeros(X,Y,Z);
        ETA = zeros(X,Y);
        TAU_X = zeros(X,Y);
        TAU_Y = zeros(X,Y);
        SFC_HFLUX = zeros(X,Y);
        PME_SBC = zeros(X,Y);
        EVAP = zeros(X,Y);
        LPREC = zeros(X,Y);
        SWFLX = zeros(X,Y);
        EVAP_HEAT = zeros(X,Y);
        SENS_HEAT = zeros(X,Y);
        LW_HEAT = zeros(X,Y);
        for vers = 1:N_versions, if w(ttt,vers) > 0
          for ens = 1:N_ens,
            file_u = [versions{vers} '/ocean_an0' num2str(ens-1) '_' num2str(year) num2str(month, '%.2d') num2str(day, '%.2d') '_u.nc'];
            file_v = [versions{vers} '/ocean_an0' num2str(ens-1) '_' num2str(year) num2str(month, '%.2d') num2str(day, '%.2d') '_v.nc'];
            file_temp = [versions{vers} '/ocean_an0' num2str(ens-1) '_' num2str(year) num2str(month, '%.2d') num2str(day, '%.2d') '_temp.nc'];
            file_salt = [versions{vers} '/ocean_an0' num2str(ens-1) '_' num2str(year) num2str(month, '%.2d') num2str(day, '%.2d') '_salt.nc'];
            file_eta = [versions{vers} '/ocean_an0' num2str(ens-1) '_' num2str(year) num2str(month, '%.2d') num2str(day, '%.2d') '_eta.nc'];
            file_force = [versions{vers} '/ocean_an0' num2str(ens-1) '_' num2str(year) num2str(month, '%.2d') num2str(day, '%.2d') '_force.nc'];
            if w_ens{vers}(ttt,ens) > 0
              U = U + w(ttt,vers)*w_ens{vers}(ttt,ens)*ncread([header file_u], 'u', [i1 j1 k1 1], [i2-i1+1 j2-j1+1 k2 1]);
              V = V + w(ttt,vers)*w_ens{vers}(ttt,ens)*ncread([header file_v], 'v', [i1 j1 k1 1], [i2-i1+1 j2-j1+1 k2 1]);
              TEMP = TEMP + w(ttt,vers)*w_ens{vers}(ttt,ens)*ncread([header file_temp], 'temp', [i1 j1 k1 1], [i2-i1+1 j2-j1+1 k2 1]);
              SALT = SALT + w(ttt,vers)*w_ens{vers}(ttt,ens)*ncread([header file_salt], 'salt', [i1 j1 k1 1], [i2-i1+1 j2-j1+1 k2 1]);
              ETA = ETA + w(ttt,vers)*w_ens{vers}(ttt,ens)*ncread([header file_eta], 'eta_t', [i1 j1 1], [i2-i1+1 j2-j1+1 1]);
              TAU_X = TAU_X + w(ttt,vers)*w_ens{vers}(ttt,ens)*mean(ncread([header file_force], 'tau_x', [i1 j1 1], [i2-i1+1 j2-j1+1 8]),3);
              TAU_Y = TAU_Y + w(ttt,vers)*w_ens{vers}(ttt,ens)*mean(ncread([header file_force], 'tau_y', [i1 j1 1], [i2-i1+1 j2-j1+1 8]),3);
              SFC_HFLUX = SFC_HFLUX + w(ttt,vers)*w_ens{vers}(ttt,ens)*mean(ncread([header file_force], 'sfc_hflux', [i1 j1 1], [i2-i1+1 j2-j1+1 8]),3);
              PME_SBC = PME_SBC + w(ttt,vers)*w_ens{vers}(ttt,ens)*mean(ncread([header file_force], 'pme_sbc', [i1 j1 1], [i2-i1+1 j2-j1+1 8]),3);
              EVAP = EVAP + w(ttt,vers)*w_ens{vers}(ttt,ens)*mean(ncread([header file_force], 'evap', [i1 j1 1], [i2-i1+1 j2-j1+1 8]),3);
              LPREC = LPREC + w(ttt,vers)*w_ens{vers}(ttt,ens)*mean(ncread([header file_force], 'lprec', [i1 j1 1], [i2-i1+1 j2-j1+1 8]),3);
              SWFLX = SWFLX + w(ttt,vers)*w_ens{vers}(ttt,ens)*mean(ncread([header file_force], 'swflx', [i1 j1 1], [i2-i1+1 j2-j1+1 8]),3);
              EVAP_HEAT = EVAP_HEAT + w(ttt,vers)*w_ens{vers}(ttt,ens)*mean(ncread([header file_force], 'evap_heat', [i1 j1 1], [i2-i1+1 j2-j1+1 8]),3);
              SENS_HEAT = SENS_HEAT + w(ttt,vers)*w_ens{vers}(ttt,ens)*mean(ncread([header file_force], 'sens_heat', [i1 j1 1], [i2-i1+1 j2-j1+1 8]),3);
              LW_HEAT = LW_HEAT + w(ttt,vers)*w_ens{vers}(ttt,ens)*mean(ncread([header file_force], 'lw_heat', [i1 j1 1], [i2-i1+1 j2-j1+1 8]),3);
            end
          end
        end, end
        read_fail = 0; % If we get this far than all the data is readable and we should exit the while loop
      catch % if one of the files failed to load, wait 2 seconds and try again
        disp('Error loading data, trying again in 2s...')
        pause(2)
      end
    end
    % Add BRAN3 - OMAPS sea level offset
    %ETA = ETA + eta_offset;
    % Fill in land/dry cells
    %ETA = inpaint_nans(ETA, 4);
    %for k = 1:Z
    %  U(:,:,k) = inpaint_nans(U(:,:,k), 4);
    %  V(:,:,k) = inpaint_nans(V(:,:,k), 4);
    %  TEMP(:,:,k) = inpaint_nans(TEMP(:,:,k), 4);
    %  SALT(:,:,k) = inpaint_nans(SALT(:,:,k), 4);
    %end
    % Insert monthly data into variable variable
    u(:,:,:,tcnt) = flipdim(U, 3);
    v(:,:,:,tcnt) = flipdim(V, 3);
    temp(:,:,:,tcnt) = flipdim(TEMP, 3);
    salt(:,:,:,tcnt) = flipdim(SALT, 3);
    eta(:,:,tcnt) = ETA;
    tau_x(:,:,tcnt) = TAU_X;
    tau_y(:,:,tcnt) = TAU_Y;
    sfc_hflux(:,:,tcnt) = SFC_HFLUX;
    pme_sbc(:,:,tcnt) = PME_SBC;
    evap(:,:,tcnt) = EVAP;
    lprec(:,:,tcnt) = LPREC;
    swflx(:,:,tcnt) = SWFLX;
    evap_heat(:,:,tcnt) = EVAP_HEAT;
    sens_heat(:,:,tcnt) = SENS_HEAT;
    lw_heat(:,:,tcnt) = LW_HEAT;
  end
  % Save NetCDF file - u, v
  file = [header_out 'OMAPS_' num2str(year) '_' num2str(month,'%.2d') '_uv.nc'];
  % x
  nccreate(file, 'x', 'Dimensions', {'ni', X, 'nj', Y});
  ncwriteatt(file, 'x', 'long_name', 'X coordinate (longitude)');
  ncwriteatt(file, 'x', 'coordinate_type', 'longitude');
  ncwriteatt(file, 'x', 'units', 'degrees_east');
  ncwrite(file, 'x', x)
  % y
  nccreate(file, 'y', 'Dimensions', {'ni', X, 'nj', Y});
  ncwriteatt(file, 'y', 'long_name', 'Y coordinate (latitude)');
  ncwriteatt(file, 'y', 'coordinate_type', 'latitude');
  ncwriteatt(file, 'y', 'units', 'degrees_north');
  ncwrite(file, 'y', y)
  % z
  nccreate(file, 'z', 'Dimensions', {'nk', Z});
  ncwriteatt(file, 'z', 'long_name', 'Z coordinate');
  ncwriteatt(file, 'z', 'coordinate_type', 'Z');
  ncwriteatt(file, 'z', 'units', 'metres');
  ncwrite(file, 'z', z)
  % t
  t_12h = t + 0.5; % add 12 hours since this is the analysis time stamp
  nccreate(file, 't', 'Dimensions', {'record', T});
  ncwriteatt(file, 't', 'long_name', 'Time');
  ncwriteatt(file, 't', 'coordinate_type', 'time');
  ncwriteatt(file, 't', 'units', 'days since 1990-01-01 00:00:00 10');
  ncwrite(file, 't', t_12h);
  % u-vel
  nccreate(file, 'u', 'Dimensions', {'ni', X, 'nj', Y, 'nk', Z, 'record', T});
  ncwriteatt(file, 'u', 'long_name', 'East component of current');
  ncwriteatt(file, 'u', 'coordinates', 'x, y, z, t');
  ncwriteatt(file, 'u', 'units', 'm s-1');
  ncwrite(file, 'u', u);
  % v-vel
  nccreate(file, 'v', 'Dimensions', {'ni', X, 'nj', Y, 'nk', Z, 'record', T});
  ncwriteatt(file, 'v', 'long_name', 'North component of current');
  ncwriteatt(file, 'v', 'coordinates', 'x, y, z, t');
  ncwriteatt(file, 'v', 'units', 'm s-1');
  ncwrite(file, 'v', v);
  % tau_x
  nccreate(file, 'tau_x', 'Dimensions', {'ni', X, 'nj', Y, 'record', T});
  ncwriteatt(file, 'tau_x', 'long_name', 'i-directed wind stress');
  ncwriteatt(file, 'tau_x', 'coordinates', 'x, y, t');
  ncwriteatt(file, 'tau_x', 'units', 'N/m^2');
  ncwrite(file, 'tau_x', tau_x);
  % tau_y
  nccreate(file, 'tau_y', 'Dimensions', {'ni', X, 'nj', Y, 'record', T});
  ncwriteatt(file, 'tau_y', 'long_name', 'j-directed wind stress');
  ncwriteatt(file, 'tau_y', 'coordinates', 'x, y, t');
  ncwriteatt(file, 'tau_y', 'units', 'N/m^2');
  ncwrite(file, 'tau_y', tau_y);
  % Save NetCDF file - temp, salt, eta
  file = [header_out 'OMAPS_' num2str(year) '_' num2str(month,'%.2d') '_tse.nc'];
  % x
  nccreate(file, 'x', 'Dimensions', {'ni', X, 'nj', Y});
  ncwriteatt(file, 'x', 'long_name', 'X coordinate (longitude)');
  ncwriteatt(file, 'x', 'coordinate_type', 'longitude');
  ncwriteatt(file, 'x', 'units', 'degrees_east');
  ncwrite(file, 'x', x-0.05)
  % y
  nccreate(file, 'y', 'Dimensions', {'ni', X, 'nj', Y});
  ncwriteatt(file, 'y', 'long_name', 'Y coordinate (latitude)');
  ncwriteatt(file, 'y', 'coordinate_type', 'latitude');
  ncwriteatt(file, 'y', 'units', 'degrees_north');
  ncwrite(file, 'y', y-0.05)
  % z
  nccreate(file, 'z', 'Dimensions', {'nk', Z});
  ncwriteatt(file, 'z', 'long_name', 'Z coordinate');
  ncwriteatt(file, 'z', 'coordinate_type', 'Z');
  ncwriteatt(file, 'z', 'units', 'metres');
  ncwrite(file, 'z', z)
  % t
  nccreate(file, 't', 'Dimensions', {'record', T});
  ncwriteatt(file, 't', 'long_name', 'Time');
  ncwriteatt(file, 't', 'coordinate_type', 'time');
  ncwriteatt(file, 't', 'units', 'days since 1990-01-01 00:00:00 10');
  ncwrite(file, 't', t_12h);
  % eta
  nccreate(file, 'eta', 'Dimensions', {'ni', X, 'nj', Y, 'record', T});
  ncwriteatt(file, 'eta', 'long_name', 'Surface Elevation');
  ncwriteatt(file, 'eta', 'coordinates', 'x, y, t');
  ncwriteatt(file, 'eta', 'units', 'metres');
  ncwrite(file, 'eta', eta);
  % temp
  nccreate(file, 'temp', 'Dimensions', {'ni', X, 'nj', Y, 'nk', Z, 'record', T});
  ncwriteatt(file, 'temp', 'long_name', 'Temperature');
  ncwriteatt(file, 'temp', 'coordinates', 'x, y, z, t');
  ncwriteatt(file, 'temp', 'units', 'degrees C');
  ncwrite(file, 'temp', temp);
  % salt
  nccreate(file, 'salt', 'Dimensions', {'ni', X, 'nj', Y, 'nk', Z, 'record', T});
  ncwriteatt(file, 'salt', 'long_name', 'Salinity');
  ncwriteatt(file, 'salt', 'coordinates', 'x, y, z, t');
  ncwriteatt(file, 'salt', 'units', 'PSU');
  ncwrite(file, 'salt', salt);
  % sfc_hflux
  nccreate(file, 'sfc_hflux', 'Dimensions', {'ni', X, 'nj', Y, 'record', T});
  ncwriteatt(file, 'sfc_hflux', 'long_name', 'total surface heat flux');
  ncwriteatt(file, 'sfc_hflux', 'coordinates', 'x, y, t');
  ncwriteatt(file, 'sfc_hflux', 'units', 'Watts/m^2');
  ncwrite(file, 'sfc_hflux', sfc_hflux);
  % pme_sbc
  nccreate(file, 'pme_sbc', 'Dimensions', {'ni', X, 'nj', Y, 'record', T});
  ncwriteatt(file, 'pme_sbc', 'long_name', 'precip-evap via sbc (liquid, frozen, evaporation)');
  ncwriteatt(file, 'pme_sbc', 'coordinates', 'x, y, t');
  ncwriteatt(file, 'pme_sbc', 'units', '(kg/m^3)*(m/sec)');
  ncwrite(file, 'pme_sbc', pme_sbc);
  % evap
  nccreate(file, 'evap', 'Dimensions', {'ni', X, 'nj', Y, 'record', T});
  ncwriteatt(file, 'evap', 'long_name', 'evaporative mass flux (>0 leaves ocean)');
  ncwriteatt(file, 'evap', 'coordinates', 'x, y, t');
  ncwriteatt(file, 'evap', 'units', '(kg/m^3)*(m/sec)');
  ncwrite(file, 'evap', evap);
  % lprec
  nccreate(file, 'lprec', 'Dimensions', {'ni', X, 'nj', Y, 'record', T});
  ncwriteatt(file, 'lprec', 'long_name', 'liquid precip into ocean (>0 enters ocean)');
  ncwriteatt(file, 'lprec', 'coordinates', 'x, y, t');
  ncwriteatt(file, 'lprec', 'units', '(kg/m^3)*(m/sec)');
  ncwrite(file, 'lprec', lprec);
  % swflx
  nccreate(file, 'swflx', 'Dimensions', {'ni', X, 'nj', Y, 'record', T});
  ncwriteatt(file, 'swflx', 'long_name', 'shortwave flux into ocean (>0 heats ocean)');
  ncwriteatt(file, 'swflx', 'coordinates', 'x, y, t');
  ncwriteatt(file, 'swflx', 'units', 'W/m^2');
  ncwrite(file, 'swflx', swflx);
  % evap_heat
  nccreate(file, 'evap_heat', 'Dimensions', {'ni', X, 'nj', Y, 'record', T});
  ncwriteatt(file, 'evap_heat', 'long_name', 'latent heat flux into ocean (<0 cools ocean)');
  ncwriteatt(file, 'evap_heat', 'coordinates', 'x, y, t');
  ncwriteatt(file, 'evap_heat', 'units', 'W/m^2');
  ncwrite(file, 'evap_heat', evap_heat);
  % sens_heat
  nccreate(file, 'sens_heat', 'Dimensions', {'ni', X, 'nj', Y, 'record', T});
  ncwriteatt(file, 'sens_heat', 'long_name', 'sensible heat into ocean (<0 cools ocean)');
  ncwriteatt(file, 'sens_heat', 'coordinates', 'x, y, t');
  ncwriteatt(file, 'sens_heat', 'units', 'W/m^2');
  ncwrite(file, 'sens_heat', sens_heat);
  % lw_heat
  nccreate(file, 'lw_heat', 'Dimensions', {'ni', X, 'nj', Y, 'record', T});
  ncwriteatt(file, 'lw_heat', 'long_name', 'longwave flux into ocean (<0 cools ocean)');
  ncwriteatt(file, 'lw_heat', 'coordinates', 'x, y, t');
  ncwriteatt(file, 'lw_heat', 'units', 'W/m^2');
  ncwrite(file, 'lw_heat', lw_heat);
end, end



%[is, js] = meshgrid([1:X], [1:Y]);
%is = is'; js = js';

%ETA_vec = ETA(:,:,1);
%i_land = find(isnan(ETA_vec(:)));
%for i = 1:length(i_land)
%  x_vec = x(:); y_vec = y(:); %is_vec = is(:);
%  x_wet = x_vec; y_wet = y_vec;
%  x_wet(i_land) = []; y_wet(i_land) = [];
%  x_land = x_vec(i_land); y_land = y_vec(i_land);
%  dist = (x_wet - x_land(i)).^2 + (y_wet - y_land(i)).^2;
%  [dist_min, i_closest] = min(dist);
%end



