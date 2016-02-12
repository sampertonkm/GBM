function [K1,K2,Kw,Ksp,KH] = equilK(Tc,S)

% equilk.m
% Function by K. Samperton, Dept. of Geosciences, Princeton University
% Created on 2014-02-11, Tuesday, 10:45, Princeton NJ USA
% Updated on 2014-02-11, Tuesday, 13:45, Princeton NJ USA

% GOAL: 1. Create a script where you input a T,S (temperature and 
%          salinity...ignore P variabliity for now), and the output is C
%          speciation equilibrium coefficients (K1, K2, KH2O, and Ksp).

% STEP 1. Convert input temperature from Centigrade to K. Note that all
% equations are calibrated for temperature in K:
T = Tc + 273.15;

% STEP 2. Define the empirical K1, K2, KH2O, and Ksp relationships as 
% defined in Zeebe and Wolf-Gladrow (2001), where lnK = f(T,S). Note that 
% the equations for K1, K2, and KH2O are based on Roy et al. (1993) and 
% DOE (1994), while the equation for Ksp is based on Mucci (1983):

ln_K1 = 2.83655 - (2307.1266/T) - (1.5529413*log(T)) - ...
    ((0.207608410 + (4.0484/T))*(S^0.5)) + ...
    (0.0846834*S) - (0.00654208*(S^1.5)) + log(1-(0.001005*S));

ln_K2 = -9.226508 - (3351.6106/T) - (0.2005743*log(T)) - ...
    ((0.106901773 + (23.9722/T))*(S^0.5)) + ...
    (0.1130822*S) - (0.00846934*(S^1.5)) + log(1-(0.001005*S));

ln_Kw = 148.96502 - (13847.26/T) - (23.6521*log(T)) + ...
    (((118.67/T) - 5.977 + (1.0495*log(T)))*(S^0.5)) - (0.01615*S);

log_Ksp = -171.9065 - (0.077993*T) + (2839.319/T) + ...
    (71.595*log10(T)) + ...
    ((-0.77712 + (0.0028426*T) + (178.34/T))*(S^0.5)) - ...
    (0.07711*S) + (0.0041249*(S^(1.5))); 

ln_KH = (9345.17/T) - 60.2409 + (23.3585*(log(T/100))) + ...
    S*(0.023517 - (0.00023656*T) + (0.0047036*((T/100)^2)));

K1  = exp(ln_K1);
K2  = exp(ln_K2);
Kw  = exp(ln_Kw);
Ksp = 10^(log_Ksp);
KH  = exp(ln_KH);

end