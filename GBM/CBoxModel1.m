% CBoxModel1.m
% Script by K. Samperton, Dept. of Geosciences, Princeton University
% Created on 2014-02-27, THU, 10:30, Princeton NJ USA
% Updated on 2014-03-13, THU, 09:00, Princeton NJ USA

% GOALS: 1. Create a box model of carbon speciation through time and track
%           changes in MT_CO2 of the system by varying inputs/outputs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 1. Reset the workspace:
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

time  = 10000;     % units of yr, total model run time

Tc = 25;              % units of deg C, temperature
S  = 35;              % salinty

% STEP 3. Solve for K values using my equilk.m function (see equilk.m for 
% details). Default values are 25 deg Centigrade and 35, respectively:
[K1,K2,Kw,Ksp,KH] = equilK(Tc,S);

% STEP 4. Define new constants to assist in algebra to solve for Hplus:
A   = (KH*MT_CO2)/Vat;
B   = ((KH*Voc)/Vat);

% STEP 5. Define the coefficients for the polynomial equation terms 1-5:
term1 = -(B+1);
term2 = (-B*K1) - (ALK*(B+1));
term3 = -(B*K1*K2) - (B*ALK*K1) + (A*K1) + (Kw*(B+1));
term4 = (B*K1*Kw) - (B*ALK*K1*K2) + (2*A*K1*K2);
term5 = B*K1*K2*Kw;

% STEP 6. Assemble the CO2 speciation polynomial, solve for HC species:
Hpolynomial = [term1 term2 term3 term4 term5];
Hroots = roots(Hpolynomial);
H         = Hroots(find(Hroots > 0)); %#ok<FNDSB>
pH        = -log10(H);
OH        = Kw/H;
CO3       = (ALK - OH + H)/((H/K2)+2);
HCO3      = CO3*H/K2;
CO2aq     = HCO3*H/K1;
CO2at_ppm = CO2aq/KH*(10^6);

% STEP 7. Set intial conditions for the model (J0):
J0 = [MT_CO2 MT_ALK Vvolc SilW Corg]';

% STEP 8. Solve the corresponding ODE as defined by the dxdt function:
[t,J] = ode15s('dxdt', time, J0);

% STEP 9. Extract parameters of interest to solve for carbon speciation
% at different time steps:
tsMT_CO2 = J(:,1);

Cspeciation = zeros(9,11);
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
end

% STEP 10. Plot the model outputs:
figure;
subplot(2,2,1); plot(t',Cspeciation(2,:)) % pH
    ylabel('pH'); xlabel('time'); grid on
subplot(2,2,2); plot(t',Cspeciation(7,:)) % pCO2_atm
    ylabel('pCO2atm'); xlabel('time'); grid on
subplot(2,2,3); plot(t',Cspeciation(8,:)) % HCO3
    ylabel('DIC'); xlabel('time'); grid on
subplot(2,2,4); plot(t',Cspeciation(9,:)) % CO2aq
    ylabel('ALK'); xlabel('time'), ylim([0 10^-2]); grid on
set(gcf, 'color', 'w');