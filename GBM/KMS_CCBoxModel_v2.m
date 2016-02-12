% KMS_CCBoxModel_v2.m
% Script by K.M. Samperton, Department of Geosciences, Princeton University
% Created on 2014-05-12, MON, 10:30, Princeton NJ USA
% Updated on 2014-05-12, MON, 10:30, Princeton NJ USA

% GOALS: 1) Quantify the compositional/isotopic evolution of the Earth 
%           system from 4.5 to 0 Ga using a box-modeling approach.
%        2) Track five global reservoirs through time, including:
%           a. Mantle
%           b. Continental crust
%           c. Oceanic crust
%           d. Seawater-atmosphere
%           e. Sediments
%        3) The isotopic systems that will be investigated are:
%           a. Rb-Sr
%           b. Sm-Nd
%           c. Lu-Hf
%           d. U-Th-Pb

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 1. Time model, clear workspace, close windows, clear command line:
tic; clear; close all; clc;

% STEP 2. Set system constants and boundary conditions:

  % Model run time (units: years) and initial reservoir masses (units: 
  % kilograms except where noted):
    time  = 4.5*10^9;   % age of the Earth (i.e., model run time)
    iMman = 4.28*10^24; % initial mass of mantle
    iMcc  = 0;          % initial mass of continental crust
    iMoc  = 0;          % initial mass of oceanic crust
    iMsw  = 1.3*10^21;  % initial mass of the ocean
    iMat  = 1.8*10^20;  % initial mass of the atmosphere (units: mol)
    iMsed = 0;          % initial mass of sediments
    
    ALK = 0.0022;        % alkalinity of seawater (units: mol/kg)
    MT_CO2 = 0.002*iMsw; % mass of total carbon (units: mol C)
    MT_ALK = ALK*iMsw;   % mass of total alkalinity (units: mol)

  % Define isotopic decay constants (units: years^-1):
    lamba_Rb_87  = 1.42*10^-11;    % Rb-87
    lamba_Sm_147 = 6.54*10^-12;    % Sm-147
    lamba_Lu_176 = 1.96*10^-11;    % Lu-176
    lamba_U_235  = 9.8485*10^-10;  % U-235
    lamba_U_238  = 1.55125*10^-10; % U-238
    lamba_Th_232 = 4.9475*10^-11;  % Th-232

  % Define starting fluxes:
    Vvolc = 2.0*10^12;    % units of mol C/yr...flux of total volcanic C input
    SilW  = 2.0*10^12;    % units of mol C/yr...flux of silicate weathering
    Corg  = 0;            % units of mol C/yr...flux of carbon burial

% STEP 3. Define an intial conditions array for the model (J0):
J0 = [iMcc iMsed Vvolc SilW Corg]';

% STEP 4. Solve the corresponding ODE, as defined by the dxdt function, 
% using the initial conditions in J0:
[t,J] = ode15s('timestep', time, J0);
trev = fliplr(t');

% STEP 5. Plot the model outputs:

  % Plot #1, mass of continental crust:
    subplot(1,3,1)
    plot(trev, J(:,1)/max(J(:,1)))
    ylabel('Mass of continental crust (normalized)', 'fontsize', 14)
    xlabel('Time (years)', 'fontsize', 14)
    set(gca, 'XDir', 'reverse')
    set(gca, 'fontsize', 14)
    xlim([0 4.5*10^9])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    grid on
    
  % Plot #2, mass of sediments:
    subplot(1,3,2)
    plot(trev, J(:,2)/max(J(:,2)))
    ylabel('Mass of sediments (normalized)', 'fontsize', 14)
    xlabel('Time (years)', 'fontsize', 14)
    set(gca, 'XDir', 'reverse')
    set(gca, 'fontsize', 14)
    xlim([0 4.5*10^9])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    grid on
  
  % Set figure properties:
    set(gcf, 'color', 'w');
    
    toc