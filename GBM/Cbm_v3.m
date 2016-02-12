% Cbm_v3.m
% Script by K. Samperton, Dept. of Geosciences, Princeton University
% Created on 2014-02-27, THU, 10:30, Princeton NJ USA
% Updated on 2014-04-03, THU, 10:00, Princeton NJ USA

% GOALS: 1. Create a box model of the carbon cycle through time that tracks
%           changes in MT_CO2, C speciation, temperature, etc. of the 
%           system by varying carbon inputs and outputs/sources and sinks.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 1. Reset the workspace and close all existing figures:
close all; clear; clc;

% STEP 2. Define system starting conditions (ALK, DIC, MT_CO2, Voc, Vat, 
% and MT_ALK), starting fluxes (Vvolc, SilW, and Corg), the model run time 
% in years (time), and the starting temperature and salinity (Tc and S):

Voc    = 1.3*10^21;   % units of kg
Vat    = 1.8*10^20;   % units of mol  
ALK    = 0.0022;      % units of mol/kg
MT_CO2 = 0.002*Voc;   % units of mol C
MT_ALK = ALK*Voc;     % units of mol, mass of total alkalinity

Vvolc = 2.0*10^12;    % units of mol C/yr, mass of total volcanic C input
SilW  = 2.0*10^12;    % units of mol C/yr, mass of silicate weathering
Corg  = 0;            % units of mol C/yr, mass of Corg burial

time  = 2500000;      % units of yr, total model run time

Tc = 25;              % units of degrees Centigrade, temperature
S  = 35;              % salinty

% STEP 3. Solve for K values using my equilk.m function (see equilk.m for 
% details and corresponding documentation). Default starting values of 
% temperature and salinity are 25 degrees Centigrade and 35, respectively:
[K1,K2,Kw,Ksp,KH] = equilK(Tc,S);

% STEP 4. Define constants to streamline arithmetic for Hplus solution (A 
% and B):
A = (KH*MT_CO2)/Vat;
B = (KH*Voc)/Vat;

% STEP 5. Define intial conditions for the model (J0):
J0 = [MT_CO2 MT_ALK Vvolc SilW Corg Tc]';

% STEP 6. Solve the corresponding ODE, as defined by the dxdt function, 
% using the initial conditions in J0:
[t,J] = ode15s('dxdt', time, J0);

% STEP 7. Extract parameters of interest to solve for carbon speciation
% at different time steps:

tsMT_CO2 = J(:,1); % generate a vector of the evolving total mass of C

Cspeciation = zeros(9,size(tsMT_CO2,1));
temp = zeros(size(tsMT_CO2,1),1);
for i = 1:numel(tsMT_CO2)
   A = (KH*(tsMT_CO2(i)))/Vat;
   B = ((KH*Voc)/Vat);
   term1 = -(B+1);
   term2 = (-B*K1) - (ALK*(B+1));
   term3 = -(B*K1*K2) - (B*ALK*K1) + (A*K1) + (Kw*(B+1));
   term4 = (B*K1*Kw) - (B*ALK*K1*K2) + (2*A*K1*K2);
   term5 = B*K1*K2*Kw;
   Hpolynomial = [term1 term2 term3 term4 term5];
   Hroots = roots(Hpolynomial);
   Cspeciation(1,i) = Hroots(find(Hroots > 0));  %#ok<FNDSB> % H
   Cspeciation(2,i) = -log10(Cspeciation(1,i));  % pH  
   Cspeciation(3,i) = Kw/Cspeciation(1,i);       % OH
   Cspeciation(4,i) = (ALK - Cspeciation(3,i) + Cspeciation(1,i))/...
       ((Cspeciation(1,i)/K2)+2);                % CO3
   Cspeciation(5,i) = Cspeciation(4,i)*Cspeciation(1,i)/K2;  % HCO3
   Cspeciation(6,i) = Cspeciation(5,i)*Cspeciation(1,i)/K1;  % CO2aq
   Cspeciation(7,i) = Cspeciation(6,i)/KH*(10^6);            % CO2at_ppm
   Cspeciation(8,i) = Cspeciation(5,i) + Cspeciation(6,i) + ...
       Cspeciation(4,i); % DIC
   Cspeciation(9,i) = Cspeciation(5,i) + 2*Cspeciation(4,i) + ...
       Cspeciation(3,i) - Cspeciation(1,i); % ALK
   temp(i) = 5*log(Cspeciation(7,i)/345)+12;   % temperature, units of deg C
end

% STEP 8. Plot the model outputs:
figure;
subplot(2,3,1); plot(t',Cspeciation(2,:)) % pH
    ylabel('pH'); xlabel('Time'); grid on
subplot(2,3,2); plot(t',Cspeciation(7,:)) % pCO2_atm
    ylabel('pCO_2atm (ppm)'); xlabel('Time'); grid on
subplot(2,3,3); plot(t',Cspeciation(8,:)) % HCO3
    ylabel('DIC (mol/kg)'); xlabel('Time'); grid on
subplot(2,3,4); plot(t',Cspeciation(9,:)) % CO2aq
    ylabel('ALK (mol/kg)'); xlabel('Time'); grid on
subplot(2,3,5); plot(t',temp) % temp
    ylabel('Temperature (degrees C)'); xlabel('Time'); grid on
set(gcf, 'color', 'w');