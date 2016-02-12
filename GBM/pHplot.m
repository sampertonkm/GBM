% pHplotter.m
% Script by K. Samperton, Dept. of Geosciences, Princeton University
% Created on 2014-02-20, Thursday, 12:00, Princeton NJ USA
% Updated on 2014-02-20, Thursday, 

% GOALS: 1. Create a script to plot carbonate ion speciation as a function
%           pH at constant alkalinity.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 1. Reset the workspace:
close all; clear; clc;

% STEP 2. Define the system alkalinity (ALK) and dissolved inorganic 
% carbon (DIC):
ALK = 0.0022; % units of mol/kg
DIC = 0.0021; % units of mol/kg

% STEP 3. Set temperature and salinity and solve for K values using my
% equilk.m function (see equilk.m for details). Default values are 25 deg
% Centigrade and 35, respectively:
Tc = 25; % temperature
S  = 35; % salinty

[K1,K2,Kw,Ksp,KH] = equilK(Tc,S);

% STEP 3. Set desired pH range for plot:
pH = 1:0.01:13;

% STEP 4. At a given pH, solve for H, OH, CO3, HCO3, CO2aq and CO2at_ppm:
for i = 1:numel(pH)
    H(i)     = 10^-pH(i);
    OH(i)    = Kw/H(i);
    CO3(i)   = DIC/(1+(H(i)/K2)+((H(i)^2)/K1*K2));
    HCO3(i)  = DIC/(1+(H(i)/K1)+(K2/H(i)));
    CO2aq(i) = DIC/(1+(K1/H(i))+(K1*K2/(H(i)^2)));
    CO2at(i) = CO2aq(i)/KH;
end

pK1    = -log10(K1);
pK2    = -log10(K2);
prange = [10^-6 10^1];
SWpH   = 8.09;

% STEP 5. Plot up the different pH arrays:
plot(pH,H, 'm-'); hold on
plot(pH,OH, 'm-.'); hold on
plot(pH,CO3, 'k'); hold on
plot(pH,HCO3, 'r'); hold on
plot(pH,CO2aq, 'g'); hold on
plot(pH,CO2at, 'y'); hold on
plot([pK1 pK1], prange, 'b'); hold on
plot([pK2 pK2], prange, 'b'); hold on
plot([SWpH SWpH], prange, 'c'); hold on

% STEP 6. Set figure properties:
title('Fig 1. Bjerrum plot of carbonate speciation in seawater');
xlabel('pH');
ylabel('Concentration (mol kg^{-1})');
set(gcf, 'color', 'w');
set(gca, 'YScale', 'log');
xlim([0 14])
ylim([10^-5 10^0])
legend('H^+', 'OH^-', 'CO_3', 'HCO_3', 'CO_{2}aq', 'CO_{2,atm}', ...
    'pK_1 = 5.86', 'pK_2 = 8.92', 'pK_{SW} = 8.09', 'Location', ...
    'EastOutside');
grid on