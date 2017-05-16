%Jessica Benthuysen
%April 27, 2016
%calc_heatbudget calculates the terms in the heat budget from OceanCurrents
%using: dT/dt = - div_H dot (uT) + d/dz(K_V * dT/dz) + R,
%where R is the residual term. The temperature tendency equation is volume
%averaged (see Peter et al 2006, JGR for a mixed layer example).
%
%The code outputs the terms as follows:
%Volume averaged temperature tendency, dT/dt = int_adv_d (horizontal advection) + Qsnet (surface heat flux)

clear all

addpath(genpath('/home/ecoliver/matlab/netcdf_toolbox/'));
addpath(genpath('/home/ecoliver/matlab/mexcdf/'));
addpath(genpath('/home/ecoliver/Desktop/include/'));
addpath(genpath('heat_budget/'));
rehash

header_OMAPS = '/data/MHWs/Tasmania_2015_2016/OMAPS/';
data_shf = 'CFSv2'; %'ERAInterim'; %'CFSv2';
%data_shf = 'GODAS'; %'OceanMAPS';
%data_shf = 'OceanMAPS';

%xw = 147; xe = 155; ys = -45; yn = -37; h = 100; % Old box
xw = 147; xe = 157; ys = -46; yn = -39; h = 100; % New box
years = 2012:2015; months = [9 10 11 12 1+12 2+12 3+12];
years = 2012:2015; months = [9 10 11 12 1+12 2+12 3+12 4+12 5+12];

time = []; T_adv = []; T_Q = []; T_tot = [];
T_adve = []; T_advw = []; T_advn = []; T_advs = [];

for year0 = years, year0

  t = []; dTdt_adv = []; dTdt_Q = []; T = [];
  dTdt_adve = []; dTdt_advw = []; dTdt_advn = []; dTdt_advs = [];

  for month0 = months
    if month0 > 12
      month = month0 - 12; year = year0 + 1;
    else
      month = month0; year = year0;
    end

    F = calc_temp_vol_avg_CFSv2(xw, xe, yn, ys, h, year, month, header_OMAPS, data_shf);
    t = [t; F.time];
    dTdt_adv = [dTdt_adv; F.int_adv_d];
    dTdt_adve = [dTdt_adve; F.int_adv_de];
    dTdt_advw = [dTdt_advw; F.int_adv_dw];
    dTdt_advn = [dTdt_advn; F.int_adv_dn];
    dTdt_advs = [dTdt_advs; F.int_adv_ds];
    dTdt_Q = [dTdt_Q; F.Qsnet];
    T = [T; F.temp_avg];
  end

  %Time-integrate the terms in order to see their net effect on the change in volume integrated temperature.
  time = [time; t];
  T_adv = [T_adv; cumtrapz(t, dTdt_adv)*(24*60*60)];
  T_adve = [T_adve; cumtrapz(t, dTdt_adve)*(24*60*60)];
  T_advw = [T_advw; cumtrapz(t, dTdt_advw)*(24*60*60)];
  T_advn = [T_advn; cumtrapz(t, dTdt_advn)*(24*60*60)];
  T_advs = [T_advs; cumtrapz(t, dTdt_advs)*(24*60*60)];
  T_Q = [T_Q; cumtrapz(t, dTdt_Q)*(24*60*60)];
  T_tot = [T_tot; T - T(1)];
end

figure; set(gcf, 'position', [2099          66         593         829]); clf;
ymin = floor(min([T_tot; T_adv; T_Q; T_tot - T_adv - T_Q])*2)/2;
ymax = ceil(max([T_tot; T_adv; T_Q; T_tot - T_adv - T_Q])*2)/2;
for year0 = years
  subaxis(length(years), 1, year0-years(1)+1); title([num2str(year0) ' / ' num2str(year0+1)]); hold on;
  plot(time, T_tot, 'k', 'linewidth', 2); %change in temperature from first day
  plot(time, T_adv, 'b', 'linewidth', 2); %change in temp. due to horiz. adv.
  plot(time, T_Q,   'r', 'linewidth', 2); %change in temp. due to sfc heat flux
  %plot(time, T_adv+T_Q, 'k--', 'linewidth', 1); %change in temp. due to horiz. adv. + sfc heat flux
  plot(time, T_tot - T_adv - T_Q, 'k--', 'linewidth', 1); % residual
  ylabel('deg C'); 
  grid on;
  xlim([datenum([year0 9 1]) datenum([year0+1 3+1 1])])
  ylim([ymin ymax]);
  datetick('x', 'mmm', 'keeplimits');
  %if year0 == years(1), legend('Total', 'Advection', 'Surface het flux', 'Adv. + Sfc. heat flux', 'location', 'NorthWest'); end
  if year0 == years(1), legend('Total', 'Advection', 'Surface heat flux', 'Residual', 'location', 'NorthWest'); end
  if year0 < years(end), set(gca, 'xticklabel', []); end
end
% set(gcf, 'color', 'w'); export_fig -png ../../documents/14_Tasmania_2015_2016/figures/heat_budget_CFSv2_100m.png -r150 -a1
% set(gcf, 'color', 'w'); export_fig -pdf ../../documents/14_Tasmania_2015_2016/figures/heat_budget_CFSv2_100m.pdf


T_tot_ds = deseason_harmonic_irreg(T_tot, time, 3, 365.25);
T_adv_ds = deseason_harmonic_irreg(T_adv, time, 3, 365.25);
T_Q_ds = deseason_harmonic_irreg(T_Q, time, 3, 365.25);
figure; set(gcf, 'position', [2099          66         593         829]); clf;
for year0 = years
  subaxis(length(years), 1, year0-years(1)+1); title([num2str(year0) ' / ' num2str(year0+1)]); hold on;
  plot(time, T_tot_ds, 'k', 'linewidth', 2); %change in temperature from first day
  plot(time, T_adv_ds, 'b', 'linewidth', 2); %change in temp. due to horiz. adv.
  plot(time, T_Q_ds,   'r', 'linewidth', 2); %change in temp. due to sfc heat flux
  plot(time, T_tot_ds - T_adv_ds - T_Q_ds,   'k--', 'linewidth', 1); % residual
  ylabel('deg C');
  grid on;
  xlim([datenum([year0 9 1]) datenum([year0+1 3+1 1])])
  ylim([-1.75 1.75]);
  datetick('x', 'mmm', 'keeplimits');
  if year0 == years(1), legend('Total', 'Advection', 'Surface heat flux', 'Residual', 'location', 'NorthWest'); end
  if year0 < years(end), set(gca, 'xticklabel', []); end
end

% set(gcf, 'color', 'w'); export_fig -png ../../documents/14_Tasmania_2015_2016/figures/heat_budget_CFSv2_100m_ds.png -r150 -a1
% set(gcf, 'color', 'w'); export_fig -pdf ../../documents/14_Tasmania_2015_2016/figures/heat_budget_CFSv2_100m_ds.pdf

% Climatology of heat budget
w = 5;
dates = datevec(time);
T_tot_clim = NaN*T_tot;
T_adv_clim = NaN*T_adv;
T_adve_clim = NaN*T_adve;
T_advw_clim = NaN*T_advw;
T_advn_clim = NaN*T_advn;
T_advs_clim = NaN*T_advs;
T_Q_clim = NaN*T_Q;
for tt = 1:length(time)
  tt_clim0 = find((dates(:,2)==dates(tt,2)) & (dates(:,3)==dates(tt,3))); % & (time < datenum([2015 9 1])));
  tt_clim = [];
  for iw = -w:w
    tt_clim = [tt_clim; tt_clim0 + iw];
  end
  tt_clim(tt_clim < 1) = [];
  tt_clim(tt_clim > length(time)) = [];
  T_tot_clim(tt) = mean(T_tot(tt_clim));
  T_adv_clim(tt) = mean(T_adv(tt_clim));
  T_adve_clim(tt) = mean(T_adve(tt_clim));
  T_advw_clim(tt) = mean(T_advw(tt_clim));
  T_advn_clim(tt) = mean(T_advn(tt_clim));
  T_advs_clim(tt) = mean(T_advs(tt_clim));
  T_Q_clim(tt) = mean(T_Q(tt_clim));
end
% Fix feb29 issue
tt_feb29 = find((dates(:,2) == 2) & (dates(:,3) == 29));
T_tot_clim(tt_feb29) = 0.5*T_tot_clim(tt_feb29-1) + 0.5*T_tot_clim(tt_feb29+1);
T_adv_clim(tt_feb29) = 0.5*T_adv_clim(tt_feb29-1) + 0.5*T_adv_clim(tt_feb29+1);
T_adve_clim(tt_feb29) = 0.5*T_adve_clim(tt_feb29-1) + 0.5*T_adve_clim(tt_feb29+1);
T_advw_clim(tt_feb29) = 0.5*T_advw_clim(tt_feb29-1) + 0.5*T_advw_clim(tt_feb29+1);
T_advn_clim(tt_feb29) = 0.5*T_advn_clim(tt_feb29-1) + 0.5*T_advn_clim(tt_feb29+1);
T_advs_clim(tt_feb29) = 0.5*T_advs_clim(tt_feb29-1) + 0.5*T_advs_clim(tt_feb29+1);
T_Q_clim(tt_feb29) = 0.5*T_Q_clim(tt_feb29-1) + 0.5*T_Q_clim(tt_feb29+1);
% Kill first and last parts of climatology
tt = find((dates(:,2) == months(1)) & (dates(:,3) == 1));
for ttt = tt'
  T_tot_clim(ttt:ttt+w) = NaN;
  T_adv_clim(ttt:ttt+w) = NaN;
  T_adve_clim(ttt:ttt+w) = NaN;
  T_advw_clim(ttt:ttt+w) = NaN;
  T_advn_clim(ttt:ttt+w) = NaN;
  T_advs_clim(ttt:ttt+w) = NaN;
  T_Q_clim(ttt:ttt+w) = NaN;
end
tt = find((dates(:,2) == mod(months(end),12)) & (dates(:,3) == 31));
for ttt = tt'
  T_tot_clim(ttt-w:ttt) = NaN;
  T_adv_clim(ttt-w:ttt) = NaN;
  T_adve_clim(ttt-w:ttt) = NaN;
  T_advw_clim(ttt-w:ttt) = NaN;
  T_advn_clim(ttt-w:ttt) = NaN;
  T_advs_clim(ttt-w:ttt) = NaN;
  T_Q_clim(ttt-w:ttt) = NaN;
end

% Calculate anomalies
T_tot_ds = T_tot - T_tot_clim;
T_adv_ds = T_adv - T_adv_clim;
T_adve_ds = T_adve - T_adve_clim;
T_advw_ds = T_advw - T_advw_clim;
T_advn_ds = T_advn - T_advn_clim;
T_advs_ds = T_advs - T_advs_clim;
T_Q_ds = T_Q - T_Q_clim;

figure; set(gcf, 'position', [2099          66         593         829]); clf;
for year0 = years
  subaxis(length(years), 1, year0-years(1)+1); title([num2str(year0) ' / ' num2str(year0+1)]); hold on;
  plot(time, T_tot_ds, 'k', 'linewidth', 2); %change in temperature from first day
  plot(time, T_adv_ds, 'b', 'linewidth', 2); %change in temp. due to horiz. adv.
  plot(time, T_Q_ds,   'r', 'linewidth', 2); %change in temp. due to sfc heat flux
  plot(time, T_tot_ds - T_adv_ds - T_Q_ds,   'k--', 'linewidth', 1); % residual
  ylabel('deg C');
  grid on;
  xlim([datenum([year0 9 1]) datenum([year0+1 3+1 1])])
  ylim([-1.75 1.75]);
  datetick('x', 'mmm', 'keeplimits');
  if year0 == years(1), legend('Total', 'Advection', 'Surface heat flux', 'Residual', 'location', 'NorthWest'); end
  if year0 < years(end), set(gca, 'xticklabel', []); end
end
% set(gcf, 'color', 'w'); export_fig -png ../../documents/14_Tasmania_2015_2016/figures/heat_budget_CFSv2_100m_ds.png -r150 -a1
% set(gcf, 'color', 'w'); export_fig -pdf ../../documents/14_Tasmania_2015_2016/figures/heat_budget_CFSv2_100m_ds.pdf

figure; set(gcf, 'position', [2098         254         593         627]); clf;
% Climatology
subaxis(3,1,1); hold on; title('Climatology (Sep 2012 - March 2015)');
plot(time, T_tot_clim, 'k', 'linewidth', 2); %change in temperature from first day
plot(time, T_adv_clim, 'b', 'linewidth', 2); %change in temp. due to horiz. adv.
plot(time, T_Q_clim,   'r', 'linewidth', 2); %change in temp. due to sfc heat flux
plot(time, T_tot_clim - T_adv_clim - T_Q_clim,   'k--', 'linewidth', 1); % residual
ylabel('deg C');
grid on;
xlim([datenum([2012 9 1]) datenum([2012+1 3+1 1])])
ylim([-1.5 5.5]);
datetick('x', 'mmm', 'keeplimits');
set(gca, 'xticklabel', []);
set(gca, 'ytick', [-1:5]);
legend('Total', 'Advection', 'Surface heat flux', 'Residual', 'location', 'NorthWest');
% 2015-2016
subaxis(3,1,2); hold on; title('2015-2016');
plot(time, T_tot, 'k', 'linewidth', 2); %change in temperature from first day
plot(time, T_adv, 'b', 'linewidth', 2); %change in temp. due to horiz. adv.
plot(time, T_Q,   'r', 'linewidth', 2); %change in temp. due to sfc heat flux
plot(time, T_tot - T_adv - T_Q,   'k--', 'linewidth', 1); % residual
ylabel('deg C');
grid on;
xlim([datenum([2015 9 1]) datenum([2015+1 3+1 1])])
ylim([-1.5 5.5]);
datetick('x', 'mmm', 'keeplimits');
set(gca, 'xticklabel', []);
set(gca, 'ytick', [-1:5]);
% 2015-2016 Anomaly
subaxis(3,1,3); hold on; title('2015-2016 Anomaly');
plot(time, T_tot_ds, 'k', 'linewidth', 2); %change in temperature from first day
plot(time, T_adv_ds, 'b', 'linewidth', 2); %change in temp. due to horiz. adv.
plot(time, T_Q_ds,   'r', 'linewidth', 2); %change in temp. due to sfc heat flux
plot(time, T_tot_ds - T_adv_ds - T_Q_ds,   'k--', 'linewidth', 1); % residual
ylabel('deg C');
grid on;
xlim([datenum([2015 9 1]) datenum([2015+1 3+1 1])])
ylim([-1.75 1.75]);
datetick('x', 'mmm', 'keeplimits');
% set(gcf, 'color', 'w'); export_fig -png ../../documents/14_Tasmania_2015_2016/figures/heat_budget_CFSv2_100m_clim_20152016_anom.png -r150 -a1
% set(gcf, 'color', 'w'); export_fig -pdf ../../documents/14_Tasmania_2015_2016/figures/heat_budget_CFSv2_100m_clim_20152016_anom.pdf
% set(gcf, 'color', 'w'); export_fig -pdf ../../documents/14_Tasmania_2015_2016/figures/heat_budget_GODAS_100m_clim_20152016_anom.pdf

figure; set(gcf, 'position', [2098         254         593         627]); clf;
% Ratio (Climatology)
subaxis(3,1,1); hold on; title('Climatology (Sep 2012 - March 2015)');
plot(time, T_tot_clim, 'k', 'linewidth', 2); %change in temperature from first day
plot(time, 100.*T_adv_clim ./ T_tot_clim, 'b', 'linewidth', 2); %change in temp. due to horiz. adv.
plot(time, 100.*T_Q_clim ./ T_tot_clim,   'r', 'linewidth', 2); %change in temp. due to sfc heat flux
ylabel('deg C');
grid on;
xlim([datenum([2012 9 1]) datenum([2012+1 3+1 1])])
ylim([-100 100])
datetick('x', 'mmm', 'keeplimits');
set(gca, 'xticklabel', []);
%set(gca, 'ytick', [-1:5]);
%legend('Total', 'Advection', 'Surface heat flux', 'Residual', 'location', 'NorthWest');
% Ratio (2015-2016)
subaxis(3,1,2); hold on; title('2015-2016');
plot(time, T_tot, 'k', 'linewidth', 2); %change in temperature from first day
plot(time, 100.*T_adv ./ T_tot, 'b', 'linewidth', 2); %change in temp. due to horiz. adv.
plot(time, 100.*T_Q ./ T_tot,   'r', 'linewidth', 2); %change in temp. due to sfc heat flux
ylabel('deg C');
grid on;
xlim([datenum([2015 9 1]) datenum([2015+1 3+1 1])])
ylim([-100 100])
datetick('x', 'mmm', 'keeplimits');
%set(gca, 'xticklabel', []);
%set(gca, 'ytick', [-1:5]);

% Horizontal advection components
figure; set(gcf, 'position', [2098         254         593         627]); clf;
% Climatology
subaxis(3,1,1); hold on; title('Climatology (Sep 2012 - March 2015)');
plot(time, T_advn_clim, 'b-', 'linewidth', 2);
plot(time, T_advs_clim, 'b--', 'linewidth', 2);
plot(time, T_advw_clim, 'r-', 'linewidth', 2);
plot(time, T_adve_clim, 'r--', 'linewidth', 2);
ylabel('deg C');
grid on;
xlim([datenum([2012 9 1]) datenum([2012+1 3+1 1])])
%ylim([-1.5 5.5]);
datetick('x', 'mmm', 'keeplimits');
set(gca, 'xticklabel', []);
%set(gca, 'ytick', [-1:5]);
% 2015-2016
subaxis(3,1,2); hold on; title('2015-2016');
plot(time, T_advn, 'b-', 'linewidth', 2);
plot(time, T_advs, 'b--', 'linewidth', 2);
plot(time, T_advw, 'r-', 'linewidth', 2);
plot(time, T_adve, 'r--', 'linewidth', 2);
ylabel('deg C');
grid on;
xlim([datenum([2015 9 1]) datenum([2015+1 3+1 1])])
%ylim([-1.5 5.5]);
datetick('x', 'mmm', 'keeplimits');
set(gca, 'xticklabel', []);
%set(gca, 'ytick', [-1:5]);
% 2015-2016 Anomaly
subaxis(3,1,3); hold on; title('2015-2016 Anomaly');
plot(time, T_advn_ds, 'b-', 'linewidth', 2);
plot(time, T_advs_ds, 'b--', 'linewidth', 2);
plot(time, T_advw_ds, 'r-', 'linewidth', 2);
plot(time, T_adve_ds, 'r--', 'linewidth', 2);
ylabel('deg C');
grid on;
xlim([datenum([2015 9 1]) datenum([2015+1 3+1 1])])
ylim([-6 6]);
datetick('x', 'mmm', 'keeplimits');
legend('Adv. North', 'Adv. South', 'Adv. West', 'Adv. East', 'location', 'NorthWest');
% set(gcf, 'color', 'w'); export_fig -pdf ../../documents/14_Tasmania_2015_2016/figures/heat_budget_CFSv2_100m_clim_20152016_anom_advDecomp.pdf

