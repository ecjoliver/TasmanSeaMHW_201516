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

addpath(genpath('/home/ecoliver/matlab/netcdf_toolbox/'))
addpath(genpath('/home/ecoliver/matlab/mexcdf/'))
rehash

%Example.
%Consider bounds to 60 m depth within the mixed layer.
xw = 149; xe = 160; ys = -45; yn = -40; 
year = 2015; month = 11;

F = calc_temp_vol_avg(xw,xe,yn,ys,year,month);

%Time-integrate the terms in order to see their net effect on the
%change in volume integrated temperature.

figure(21); 
plot(F.time, F.temp_avg - F.temp_avg(1), 'k'); %change in temperature from first day
hold on; grid on; datetick('x', 'ddmmm', 'keepticks'); 
plot(F.time, cumsum(F.int_adv_d)*(24*60*60), 'c'); %change in temp. due to horiz. adv.
plot(F.time, cumsum(F.Qsnet)*(24*60*60), 'r'); %change in temp. due to surface heat flux
Residual = F.temp_avg - F.temp_avg(1) - (cumsum(F.int_adv_d)*(24*60*60)) - ...
            cumsum(F.Qsnet)*(24*60*60); %residual contribution 
%to the change in volume averaged temperature not explained by the above two terms
plot(F.time, Residual, 'k', 'color', [.5 .5 .5]); 
ylabel('deg C'); 
title('Change in volume averaged temperature (time-integral of dT/dt equation)'); 
legend('dT', 'horiz adv', 'Qsnet', 'residual', 'location', 'northwest');  




