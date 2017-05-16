%Jessica Benthuysen
%21 April 2016
%grid_setup sets up the grid for OFAM/SE Australia.

clear all

addpath(genpath('/home/ecoliver/matlab/netcdf_toolbox/'));
addpath(genpath('/home/ecoliver/matlab/mexcdf/'));
addpath(genpath('/home/ecoliver/Desktop/include/'));
rehash

%path(path, '/mnt/share/jbenthuy/Desktop');

cdfname1 = 'OMAPS_2015_11_tse.nc';  
cdfname2 = 'OMAPS_2015_11_uv.nc';  

cdf = netcdf(cdfname1, 'r'); 

%  'x'    'y'    'z'    't'    'eta'    'temp'    'salt'    'sfc_hflux'    'pme_sbc'    'evap'
%  'lprec'    'swflx'    'evap_heat'    'sens_heat'    'lw_heat'

%MOM4p1, p 241: Heat flux [W/m2]= Qlongwave radiation leaving the ocean - Qshortwave radiation
%entering the ocean + latent heat of vaporization + Q sensible heat transfer proportional
%to the difference between atmosphere and ocean temperatures.
%
%nc{'sfc_hflux'}.long_name = ncchar(''total surface heat flux'');
%nc{'swflx'}.long_name = ncchar(''shortwave flux into ocean (>0 heats ocean)'');
%nc{'evap_heat'}.long_name = ncchar(''latent heat flux into ocean (<0 cools ocean)'');
%nc{'sens_heat'}.long_name = ncchar(''sensible heat into ocean (<0 cools ocean)'');
%nc{'lw_heat'}.long_name = ncchar(''longwave flux into ocean (<0 cools ocean)'');

sfc_hflux = cdf{'sfc_hflux'}(1,:,:);
xt = cdf{'x'}(:); 
yt = cdf{'y'}(:); 
zt = cdf{'z'}(:); 
t = cdf{'t'}(:); %actual time, does not reference to 1 Jan 1990.  
 
%swflx = cdf{'swflx'}(1,:,:);
%evap_heat = cdf{'evap_heat'}(1,:,:);
%sens_heat = cdf{'sens_heat'}(1,:,:);
%lw_heat = cdf{'lw_heat'}(1,:,:);

temp = cdf{'temp'}(:); %is a NaN over land
[nt, nzt, nyt, nxt] = size(temp); %30, 40, 161, 181

%yt(60,1) is -40.05;

figure; 
pcolor(xt,yt,squeeze(temp(1,end,:,:)))

figure;
pcolor(ones(nzt,1)*xt(60,:), zt*ones(1,nxt), squeeze(temp(1,:,60,:))) 
shading flat; 

%xt(60,70): -40.05 and 149.95
%vertical structure of velocity
figure; plot(squeeze(temp(1,:,60,:)), zt)

%xt(1,1:3) %    143.0500  143.1500  143.2500
%yt(1:3,1)' %    -45.9500  -45.8500  -45.7500

%I think we chose 60 m from the ARGO decay scale of the warming. 

%%%%%%%%%%%%%%%%%%
close(cdf); 
cdf = netcdf(cdfname2, 'r'); 

xu = cdf{'x'}(:); 
yu = cdf{'y'}(:); 
zu = cdf{'z'}(:); 

u = cdf{'u'}(:); 
v = cdf{'v'}(:);

tau_x = cdf{'tau_x'}(:);  
tau_y = cdf{'tau_y'}(:);  

[ntu, nzu, nyu, nxu] = size(u); %30, 40, 161, 181 

%xu(1,1:3) %  143.0000  143.1000  143.2000
%yu(1:3,1)' %  -46.0000  -45.9000  -45.8000

figure;
pcolor(ones(nzu,1)*xu(60,:), zu*ones(1,nxu), squeeze(v(1,:,60,:))) 
shading flat; colormap('redblue'); caxis([-1 1]) 

%Vertical grid bounds are from BRAN- note that they are not flipped 
%at present: st_bnds_u, st_bnds_v, st_bnds_ts.
%save('BRAN_grid_bounds.mat', 'st_bnds_u', 'st_bnds_v', 'st_bnds_ts')
%They are the same vertical bounds for the different variables.
load('BRAN_grid_bounds.mat'); 
st_bnds_u2 = flipud(st_bnds_u); 
%[st_bnds_u2(1:10,:) mean(st_bnds_u2(1:10,:),2) -zu(1:10)] 

%DZ is the height of the vertical cell centered about the point;
%same for temperature and zonal velocity.
for iz = 1:nzu; DZ(iz) = diff(st_bnds_u2(iz,:)); end; 

%%Create land mask
temp2   = squeeze(temp(1,:,:,:)); 
u2      = squeeze(u(1,:,:,:)); 
mask_ts = ones(nzt,nyt,nxt); 
mask_uv = ones(nzu,nyu,nxu); 
%Note that u,v,temp are NaN over land.
for iz = 1:nzu; for iy = 1:nyu; for ix = 1:nxu;
if isnan(temp2(iz,iy,ix)) == 1; mask_ts(iz,iy,ix) = NaN; end; 
if isnan(u2(iz,iy,ix)) == 1; mask_u(iz,iy,ix) = NaN; end; 
end; end; end;

%save grid here....
OFAM_GRID.z_bnds = st_bnds_u2; 
OFAM_GRID.xt = xt;
OFAM_GRID.yt = yt; 
OFAM_GRID.zt = zt; 
OFAM_GRID.xu = xu;
OFAM_GRID.yu = yu;
OFAM_GRID.zu = zu; 
OFAM_GRID.mask_ts = mask_ts; %1 over water, NaN over land
OFAM_GRID.mask_uv = mask_uv; 

save('OFAM_grid_data_SEAustralia.mat', 'OFAM_GRID'); 



return; 


