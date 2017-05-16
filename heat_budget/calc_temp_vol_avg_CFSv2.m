%Jessica Benthuysen
%26 April 2016
%calc_temp_vol_avg calculates the flux terms at the edges
%of a specified area, the volume-averaged temperature,
%and the surface flux averaged over the surface area.
%
%Input: area bounds defined by xw, xe, yn, and ys.
%For simplicitly of calculation, use values on the xu,yu-grid 
%i.e. on the face of the temperature-cell, and away from land points
%Currently,ensure that temperature cells are located one grid cell away 
%from land points, down to 60 m depth.
%To 60 m depth and ys = -45 N, yn = -40 N, then xe >= 149.00 deg E
%
%Note that: xu(1,1:3)  %  143.0000  143.1000  143.2000
%           yu(1:3,1)' %  -46.0000  -45.9000  -45.8000
%
%Input: year, month: of oceanmaps file
%
%The flux terms are then depth-averaged over 60 m (can be changed later).
%
%outputs the volume-averaged temperature tendency term from horiz. advection,
%calculated from applying Gauss' theorem as follows
% = -1/h int^0_-h (1/A int^dA div_h dot (uT) dA) dz
% = -1/h int^0_-h (1/A {(uT)(eastern face) dy(eastern) - 
%                      (uT)(western face) dy(western) + 
%                      (vT)(northern face) dx(northern) -
%                      (vT)(southern face) dx(southern)}) dz
%
%Area defined by bounds is A = R^2 (xw - xe)*(sin(yn) - sin(ys))
%
%Output is units of deg.C s^-1 for the volume averaged advection term;
%deg. C for the volume-averaged temperature term itself,
%deg. C s^-1 for the area averaged surface flux term
%Output is for all days in the specified month.

function f = calc_temp_vol_avg(xw, xe, yn, ys, h, year, month, header_OMAPS, data_shf);

%%%%%%%%%%%%%%%%%%%%%%%%
%%Example.
%%Consider bounds, slightly modified (see below).
%xw = 149; xe = 160; ys = -45; yn = -40; 
%year = 2015; month = 11;
%%%%%%%%%%%%%%%%%%%%%%%%%
filename_ts = [header_OMAPS 'OMAPS_' num2str(year) '_' num2str(month, '%.2d') '_tse.nc']; 
filename_uv = [header_OMAPS 'OMAPS_' num2str(year) '_' num2str(month, '%.2d') '_uv.nc'];
if strcmp(data_shf, 'OceanMAPS')
  filename_hf = [header_OMAPS 'OMAPS_' num2str(year) '_' num2str(month, '%.2d') '_tse.nc'];
elseif strcmp(data_shf, 'CFSv2')
  filename_hf = [header_OMAPS '../' data_shf '/' data_shf '_' num2str(year) '_heatflux.nc'];
elseif strcmp(data_shf, 'ERAInterim')
  filename_hf = [header_OMAPS '../' data_shf '/' data_shf '_' num2str(year) '_heatflux.nc'];
elseif strcmp(data_shf, 'GODAS')
  filename_hf = [header_OMAPS '../' data_shf '/thflx.' num2str(year) '.nc'];
end

% depth of mixed layer
if h == 60
  izh = 10;
elseif h == 100
  izh = 14;
elseif h == 150;
  izh = 19;
elseif h == 200;
  izh = 24;
else
  error('Unsupported depth');
end

R = 6371*1000; %m, earth radius specified in MOM4p1 (noted in mom4p0_manual)
A = R*R*(xe - xw)*(2*pi/360)*(sind(yn) - sind(ys)); %m^2, area of domain

%load the depth-bounds on the cell (these are the same for u,v,T).
load OFAM_grid_data_SEAustralia.mat
nz = length(OFAM_GRID.zt); 
z_bnds = OFAM_GRID.z_bnds(nz-izh:nz,1); z_bnds(1) = h; %changed from 61.0355 m
DZ = -diff(z_bnds); %DZ is the depth of each cell

%Determine indices corresponding to xw,xe,yn,ys;
%note that bounds must exclude land points for flux calculations

xt = OFAM_GRID.xt; yt = OFAM_GRID.yt; 
xu = OFAM_GRID.xu; yu = OFAM_GRID.yu; 

iyn = find(yu(:,1) == yn); iys = find(yu(:,1) == ys); 
ixe = find(xu(1,:) == xe); ixw = find(xu(1,:) == xw); 
if isnan(iyn) == 1; display('yn not in u-grid'); end; 
if isnan(iys) == 1; display('ys not in u-grid'); end; 
if isnan(ixw) == 1; display('xw not in u-grid'); end; 
if isnan(ixe) == 1; display('xe not in u-grid'); end; 

%Eric's original bounds include land points at depth.
%xw = 145; xe = 160; ys = -45; yn = -40; 
%figure; pcolor(squeeze(temp2(nz-10,:,:)));
%From -45 to -40 N, and down to 60 m depth, xt >= 149.00 (ixe >= 61 on u-grid)
%find(isnan(OFAM_GRID.mask_ts(nz-10,10:(10 + 5*10 + 1),59)))
%find(isnan(OFAM_GRID.mask_uv(nz-10,10:(10 + 5*10 + 1),59)))

%figure;pcolor(xt,yt,temp2); shading flat; hold on;
%plot(xt(iyn,ixw:ixe),yt(iyn,ixw:ixe), 'r')

cdf_ts = netcdf(filename_ts, 'r'); 
cdf_uv = netcdf(filename_uv, 'r'); 

time = cdf_ts{'t'}(:); nt = length(time); 

%Define temperature grid for domain:
%IXt = ixw:ixe-1;  IYt = iys:iyn-1; %OLD 26 April
IXt = ixw+1:ixe;  IYt = iys+1:iyn; %new grid
XT = xt(IYt,IXt); YT = yt(IYt,IXt); 
[ny,nx] = size(XT); 
%Calculate the area of each T grid cell, [dA] is m^2.
for iy = 1:ny; for ix = 1:nx;
  dA(iy,ix) = R*R*.1*(2*pi/360)*(sind(YT(iy,ix) + .05) - sind(YT(iy,ix) - .05));  
end; end;

%Calculate width of T grid cell face along edges of domain
% (spacing is 0.1 deg N and 0.1 deg E):
dy_w = R*(.1*2*pi/360); %m, meridional length
dy_e = R*(.1*2*pi/360); %m, meridional length
dx_n = R*cosd(yn)*(.1*2*pi/360); %m, zonal width 
dx_s = R*cosd(ys)*(.1*2*pi/360); %m, zonal width 
int_dy_w = sum(ny*dy_w); int_dy_e = sum(ny*dy_e); 
int_dx_n = sum(nx*dx_n); int_dx_s = sum(nx*dx_s); 

%Load temperature (deg C) and surface heat flux (W/m^2, > 0 heats ocean)
temp = cdf_ts{'temp'}(:,nz-(izh-1):nz,IYt,IXt); %nt x (depth of top 60 m) x ny x nx
%
% EDITS FOR CFSv2 HEAT FLUX
%
if strcmp(data_shf, 'OceanMAPS')
  sfc_hflux = cdf_ts{'sfc_hflux'}(:,IYt,IXt); %nt x ny x nx
  mask_hf = ones(length(IYt), length(IXt)); mask_hf(isnan(squeeze(mean(sfc_hflux,1)))) = 0;
elseif strcmp(data_shf, 'CFSv2')
  cdf_hf = netcdf(filename_hf, 'r');
  time_hf = cdf_hf{'t'}(:); dates_hf = datevec(time_hf + datenum([1990 1 1]));
  t1 = find(dates_hf(:,2) == month, 1, 'first'); t2 = find(dates_hf(:,2) == month, 1, 'last');
  x_hf = cdf_hf{'x'}(1,:); i1 = find(x_hf >= xw, 1, 'first'); i2 = find(x_hf <= xe, 1, 'last');
  y_hf = cdf_hf{'y'}(:,1)'; j1 = find(y_hf >= ys, 1, 'first'); j2 = find(y_hf <= yn, 1, 'last');
  sfc_hflux = cdf_hf{'heatflux'}(t1:t2,j1:j2,i1:i2); % Pull out surface heat flux for month of interest, and restricted to domain
  dates_hf = dates_hf(t1:t2,:); x_hf = x_hf(i1:i2); y_hf = y_hf(j1:j2);
  for iy = 1:length(y_hf); for ix = 1:length(x_hf)
    dA_hf(iy,ix) = R*R*.1*(2*pi/360)*(sind(y_hf(iy) + .05) - sind(y_hf(iy) - .5));
  end, end
  mask_hf = ones(length(y_hf), length(x_hf)); mask_hf(isnan(squeeze(mean(sfc_hflux,1)))) = 0;
elseif strcmp(data_shf, 'ERAInterim')
  cdf_hf = netcdf(filename_hf, 'r');
  time_hf = cdf_hf{'t'}(:); dates_hf = datevec(time_hf + datenum([1990 1 1]));
  t1 = find(dates_hf(:,2) == month, 1, 'first'); t2 = find(dates_hf(:,2) == month, 1, 'last');
  x_hf = cdf_hf{'x'}(1,:); i1 = find(x_hf >= xw, 1, 'first'); i2 = find(x_hf <= xe, 1, 'last');
  y_hf = cdf_hf{'y'}(:,1)'; j1 = find(y_hf >= ys, 1, 'first'); j2 = find(y_hf <= yn, 1, 'last');
  sfc_hflux = cdf_hf{'heatflux'}(t1:t2,j1:j2,i1:i2); % Pull out surface heat flux for month of interest, and restricted to domain
  dates_hf = dates_hf(t1:t2,:); x_hf = x_hf(i1:i2); y_hf = y_hf(j1:j2);
  for iy = 1:length(y_hf); for ix = 1:length(x_hf)
    dA_hf(iy,ix) = R*R*.1*(2*pi/360)*(sind(y_hf(iy) + .05) - sind(y_hf(iy) - .5));
  end, end
  mask_hf = ones(length(y_hf), length(x_hf)); mask_hf(isnan(squeeze(mean(sfc_hflux,1)))) = 0;
elseif strcmp(data_shf, 'GODAS')
  cdf_hf = netcdf(filename_hf, 'r');
  x_hf = cdf_hf{'lon'}; i1 = find(x_hf >= xw, 1, 'first'); i2 = find(x_hf <= xe, 1, 'last');
  y_hf = cdf_hf{'lat'}; j1 = find(y_hf >= ys, 1, 'first'); j2 = find(y_hf <= yn, 1, 'last');
  x_hf = x_hf(i1:i2); y_hf = y_hf(j1:j2);
  sfc_hflux = cdf_hf{'thflx'}(month,j1:j2,i1:i2); % Pull out surface heat flux for month of interest, and restricted to domain
  sfc_hflux(sfc_hflux==32767) = NaN;
  sfc_hflux = sfc_hflux*0.0855;
  for iy = 1:length(y_hf); for ix = 1:length(x_hf)
    dA_hf(iy,ix) = R*R*.1*(2*pi/360)*(sind(y_hf(iy) + .05) - sind(y_hf(iy) - .5));
  end, end
  mask_hf = ones(length(y_hf), length(x_hf)); mask_hf(isnan(squeeze(mean(sfc_hflux,1)))) = 0;
end
%
% END EDITS
%

%Determine fluxes on the faces of the domain, where the faces are given by the
%u-grid and temperature points are in the interior.
%if h == 60; 
  ue1 = cdf_uv{'u'}(:,nz-(izh-1):nz,iys:iyn,ixe); %on east face of domain 
  uw1 = cdf_uv{'u'}(:,nz-(izh-1):nz,iys:iyn,ixw); %on west face of domain
  vn1 = squeeze(cdf_uv{'v'}(:,nz-(izh-1):nz,iyn,ixw:ixe)); %on north face of domain
  vs1 = squeeze(cdf_uv{'v'}(:,nz-(izh-1):nz,iys,ixw:ixe)); %on south face of domain
  close(cdf_uv)
  %average velocity onto the center of each T grid cell face
  for iy = 1:length(iys:iyn - 1) %number of T grid cells in the meridional direction
    ue(:,:,iy) = mean(ue1(:,:,iy:iy+1),3); 
    uw(:,:,iy) = mean(uw1(:,:,iy:iy+1),3); 
  end; 
  for ix = 1:length(ixw:ixe - 1) %number of T grid cells in the zonal direction
    vn(:,:,ix) = mean(vn1(:,:,ix:ix+1),3); 
    vs(:,:,ix) = mean(vs1(:,:,ix:ix+1),3); 
  end; 
  %
  %average temperature onto the faces of T grid cell:
  %mean(xt(iys,ixw-1:ixw)) %xw %mean(xt(iys,ixe-1:ixe)) %xe %OLD 26 April
  %mean(xt(iys,ixw:ixw+1)) %xw %new grid
  %mean(xt(iys,ixe:ixe+1)) %xe
  %tw = squeeze(mean(cdf_ts{'temp'}(:,nz-9:nz,iys:iyn-1,ixw-1:ixw),4)); %OLD 
  %te = squeeze(mean(cdf_ts{'temp'}(:,nz-9:nz,iys:iyn-1,ixe-1:ixe),4)); %OLD
  tw = squeeze(mean(cdf_ts{'temp'}(:,nz-(izh-1):nz,iys+1:iyn,ixw:ixw+1),4)); 
  te = squeeze(mean(cdf_ts{'temp'}(:,nz-(izh-1):nz,iys+1:iyn,ixe:ixe+1),4)); 
  %mean(yt(iys-1:iys,ixw)) %ys %mean(yt(iyn-1:iyn,ixw)) %yn %OLD grid 26 April
  %mean(yt(iys:iys+1,ixw)) %ys %new grid
  %mean(yt(iyn:iyn+1,ixw)) %yn
  %ts = squeeze(mean(cdf_ts{'temp'}(:,nz-9:nz,iys-1:iys,ixw:ixe-1),3)); %OLD 
  %tn = squeeze(mean(cdf_ts{'temp'}(:,nz-9:nz,iyn-1:iyn,ixw:ixe-1),3)); 
  ts = squeeze(mean(cdf_ts{'temp'}(:,nz-(izh-1):nz,iys:iys+1,ixw+1:ixe),3)); 
  tn = squeeze(mean(cdf_ts{'temp'}(:,nz-(izh-1):nz,iyn:iyn+1,ixw+1:ixe),3)); 
%end;
%figure; pcolor(time*ones(1,10),ones(nt,1)*OFAM_GRID.zu(nz-9:nz)',squeeze(vn(:,:,1)))
%colormap('redblue'); caxis([-1 1]); datetick('x', 'ddmmmyy'); 
close(cdf_ts)

%
% NEW CODE - INTEGRATES OVER Z, DEPENDING ON LOCAL DEPTH, THEN INTEGRATES OVER X/Y
%
% Depth-average terms - note that depth varies with (x,y) so denominator must change
int_utw_d = zeros(nt,size(uw,3));
for i = 1:size(uw,3)
  utw = uw(:,:,i).*tw(:,:,i);
  wet = ~isnan(utw(1,:));
  int_utw_d(:,i) = nansum(utw.*repmat(DZ', [nt 1]), 2)/sum(DZ(wet)); % deg C x m/s; nt x ny
end
%
int_ute_d = zeros(nt,size(ue,3));
for i = 1:size(ue,3)
  ute = ue(:,:,i).*te(:,:,i);
  wet = ~isnan(ute(1,:));
  int_ute_d(:,i) = nansum(ute.*repmat(DZ', [nt 1]), 2)/sum(DZ(wet)); % deg C x m/s; nt x ny
end
%
int_vtn_d = zeros(nt,size(vn,3));
for i = 1:size(vn,3)
  vtn = vn(:,:,i).*tn(:,:,i);
  wet = ~isnan(vtn(1,:));
  int_vtn_d(:,i) = nansum(vtn.*repmat(DZ', [nt 1]), 2)/sum(DZ(wet)); % deg C x m/s; nt x nx
end
%
int_vts_d = zeros(nt,size(vs,3));
for i = 1:size(vs,3)
  vts = vs(:,:,i).*ts(:,:,i);
  wet = ~isnan(vts(1,:));
  int_vts_d(:,i) = nansum(vts.*repmat(DZ', [nt 1]), 2)/sum(DZ(wet)); % deg C x m/s; nt x nx
end

% Integrate zonally and meridionally along domain edges:
int_utw = nansum(int_utw_d*dy_w,2); % deg C x m^2/s; nt (number of days)
int_ute = nansum(int_ute_d*dy_e,2); % deg C x m^2/s; nt (number of days)
int_vtn = nansum(int_vtn_d*dx_n,2); % deg C x m^2/s; nt (number of days)
int_vts = nansum(int_vts_d*dx_s,2); % deg C x m^2/s; nt (number of days)

%int_adv_d is the contribution to the volume-averaged temperature tendency equation,
%on the right hand side, i.e. note the minus sign.
int_adv_d = -(1/A)*(int_ute - int_utw + int_vtn - int_vts); %deg C x s^-1
%
% OLD CODE - INTEGRATES OVER X/Y AND LATER OVER Z
%
%Integrate zonally and meridionally along domain edges:
%for it = 1:nt
%  utw = squeeze(uw(it,:,:)).*squeeze(tw(it,:,:)); 
%  int_utw(it,:) = sum(utw*dy_w,2)'; %deg C x m^2/s; nt x (depth over top 60 m) 
%  %
%  ute = squeeze(ue(it,:,:)).*squeeze(te(it,:,:)); 
%  int_ute(it,:) = sum(ute*dy_e,2)'; %nt x (depth over top 60 m) 
%  %
%  vtn = squeeze(vn(it,:,:)).*squeeze(tn(it,:,:));  
%  int_vtn(it,:) = sum(vtn*dx_n,2)'; %nt x (depth over top 60 m) 
%  %
%  vts = squeeze(vs(it,:,:)).*squeeze(ts(it,:,:));  
%  int_vts(it,:) = sum(vts*dx_s,2)'; %nt x (depth over top 60 m) 
%end; 

%Depth average terms, sum for the total contribution from horizontal advection,
%and divide by total area of domain.
%int_utw_d = sum(int_utw.*(ones(nt,1)*DZ'),2)/sum(DZ); %deg C x m^2/s; nt (number of days)
%int_ute_d = sum(int_ute.*(ones(nt,1)*DZ'),2)/sum(DZ); %nt (number of days)
%int_vtn_d = sum(int_vtn.*(ones(nt,1)*DZ'),2)/sum(DZ); %nt (number of days)
%int_vts_d = sum(int_vts.*(ones(nt,1)*DZ'),2)/sum(DZ); %nt (number of days)

%int_adv_d is the contribution to the volume-averaged temperature tendency equation,
%on the right hand side, i.e. note the minus sign.
%int_adv_d = -(1/A)*(int_ute_d - int_utw_d + int_vtn_d - int_vts_d); %deg C x s^-1
%
% OLD CODE END
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Volume-average temperature within domain and area-average surface heat flux.
%(note that sum(sum(dA)) is very close to total area 'A' itself)

DDZ = repmat(DZ, [1 ny nx]);
for it = 1:nt; 
  %
  % EDITS FOR CFSv2 HEAT FLUX
  %
  if strcmp(data_shf, 'OceanMAPS')
    sfc_hflux_avg(it,1) = nansum(nansum(squeeze(sfc_hflux(it,:,:)).*dA))/sum(sum(dA.*mask_hf)); %W/m^2 
  elseif strcmp(data_shf, 'CFSv2')
    [year month day] = datevec(time(it)); % get year month day of OceanMAPS, and find which indices from 6-hourly CFSv2 are on that day
    which = find((dates_hf(:,1) == year) .* (dates_hf(:,2) == month) .* (dates_hf(:,3) == day));
    tmp = squeeze(mean(sfc_hflux(which,:,:),1)); sfc_hflux_avg(it,1) = nansum(nansum(tmp.*dA_hf))/sum(sum(dA_hf.*mask_hf)); % Daily and spatial average
  elseif strcmp(data_shf, 'ERAInterim')
    [year month day] = datevec(time(it)); % get year month day of OceanMAPS, and find which indices from 6-hourly CFSv2 are on that day
    which = find((dates_hf(:,1) == year) .* (dates_hf(:,2) == month) .* (dates_hf(:,3) == day));
    tmp = squeeze(mean(sfc_hflux(which,:,:),1)); sfc_hflux_avg(it,1) = nansum(nansum(tmp.*dA_hf))/sum(sum(dA_hf.*mask_hf)); % Daily and spatial average
  elseif strcmp(data_shf, 'GODAS')
    sfc_hflux_avg(it,1) = nansum(nansum(sfc_hflux.*dA_hf))/sum(sum(dA_hf.*mask_hf)); %W/m^2 
  end
  %
  % END EDITS
  %
  %
  % NEW DEPTH-AREA AVG TEMP
  %
  % Depth-average temp - note that depth varies with (x,y) so denominator must change
  wet = squeeze(~isnan(temp(it,:,:,:)));
  temp_d = squeeze(nansum(squeeze(temp(it,:,:,:)).*DDZ,1)) ./ squeeze(sum(DDZ.*wet, 1));
  % Area average temp - note that dry cells must be ignored
  mask_temp = ones(size(temp_d)); mask_temp(isnan(temp_d)) = 0;
  temp_avg(it,:) = nansum(nansum(temp_d.*dA))/sum(sum(dA.*mask_temp));
  %
  % OLD CODE
  %
  %  for iz = 1:length(DZ)
  %    temp_avg_tmp(1,iz) = sum(sum(squeeze(temp(it,iz,:,:)).*dA))/sum(sum(dA)); %deg C, area avg
  %  end; %iz 
  %  temp_avg(it,:) = sum(temp_avg_tmp.*DZ')/sum(DZ); %deg C, volume avg
end; %it

%Convert surface heat flux units:
cp = 3990; %J x deg C^-1 x kg^-1, specific heat at constant pressure
rho_ref = 1035; %kg/m^3, reference density applied in MOM4p1

Qsnet = sfc_hflux_avg/(h*cp*rho_ref); %deg C x s^-1 
%Note units:   {W x m^2}/ {(m) x (J x deg C^-1 x kg^-1) x (kg x m^-3)}
%            = {kg x s^-3} / {(m) x (kg m^2 s^-2 degC^-1 kg^-1) x (kg m^-3)}
%            = {s^-3} x {s^2 deg C}
%            = deg C x s^-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Save variables

f.time      = time; %days of the month, year
f.int_adv_d = int_adv_d; %deg C x s^-1, volume-averaged horizontal advection term on RHS of T-eqn
f.temp_avg  = temp_avg; %deg C, volume averaged temperature in domain; need to take the time 
%difference of this term for left hand side of temperature tendency equation (for units of deg C x s^-1)
f.Qsnet     = Qsnet; %deg C x s^-1, area averaged surface heat flux ( > 0 for heat into the ocean)

% Individual componentes of advection contribution
f.int_adv_de = -(1/A)*int_ute;
f.int_adv_dw =  (1/A)*int_utw;
f.int_adv_dn = -(1/A)*int_vtn;
f.int_adv_ds =  (1/A)*int_vts;

%Volume averaged temperature tendency, dT/dt = int_adv_d (horizontal advection) + Qsnet (surface heat flux)

%return;

%Plot example
%figure(20); 
%plot(f.time, f.int_adv_d, 'c'); hold on;
%plot(f.time, f.Qsnet, 'r'); 
%dTdt = (f.temp_avg(2:end) - f.temp_avg(1:end-1))/(24*60*60); 
%time2 = .5*(f.time(2:end) + f.time(1:end-1)); 
%plot(time2, dTdt, 'k'); datetick('x', 'ddmmm', 'keepticks'); grid on; 
%ylabel('\circC s^{-1}'); title('dT/dt (black), horiz adv (blue), surf heat flux (red)'); 

%For Nov 2015, it appears that advection and surface heat flux follow 
%similar temporal pattern;
%net surface heat flux into ocean and the warming is being offset by advection
%as well as some other process that is part of the 'residual' that is not incorporated 
%into the calculation.
%
%If we time-integrate the above contributions to the temperature tendency equation
%as in Figure 12 of the CSR paper, we will find a less noisy/net effect of the processes
%on the total volume averaged temperature. 












