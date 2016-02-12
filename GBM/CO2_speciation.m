% CO2_speciation.m
% Script by K. Samperton, Dept. of Geosciences, Princeton University
% Created on 2014-02-13, Thursday, 10:30, Princeton NJ USA
% Updated on 2014-02-20, Thursday, 11:50, Princeton NJ USA

% GOALS: 1. Create a program with inputs of alkalinity (ALK), total mass 
%        of CO2 in the ocean+atmosphere (MT_CO2), volume of the ocean 
%        (Voc), and volume of atmosphere (Vat), and determine the 
%        corresponding H+ concentration (H) for the system. 
%        2. Use the calculated H+ to determine C ion speciation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 1. Reset the workspace:
close all; clear; clc;

% STEP 2. Define five system constants (ALK, DIC, MT_CO2, Voc, and Vat):
Voc    = 1.3*10^21;   % units of kg
Vat    = 1.8*10^20;   % units of mol  
ALK    = 0.0022;      % units of mol/kg
MT_CO2 = 0.002*Voc;   % units of mol C

% STEP 3. Set temperature and salinity and solve for K values using my
% equilk.m function (see equilk.m for details). Default values are 25 deg
% Centigrade and 35, respectively:
Tc = 25; % temperature
S  = 35; % salinty

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
H         = Hroots(find(Hroots > 0))
pH        = -log10(H)
OH        = Kw/H
CO3       = (ALK - OH + H)/((H/K2)+2)
HCO3      = CO3*H/K2
CO2aq     = HCO3*H/K1
CO2at_ppm = CO2aq/KH*(10^6)