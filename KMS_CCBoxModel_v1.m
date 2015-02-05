% KMS_CCBoxModel_v1.m
% Script by K.M. Samperton, Department of Geosciences, Princeton University
% Created on 2014-04-22, TUE, 11:00, Princeton NJ USA
% Updated on 2014-05-12, MON, 09:30, Princeton NJ USA

% GOALS: 1) Model the growth and evolution of the continental crust from
%           4.5 to 0 Ga using a box modeling-approach.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 1. Reset the workspace and close all figures:
close all; clear; clc;

% STEP 2. Define the model run time, system starting conditions, and 
% starting fluxes:

time  = 4.5*10^9;   % units of yr...total model run time

iMcc  = 0;          % units of kg...initial mass of continental crust
iMsed = 0;          % units of kg...initial mass of sediments

ALK    = 0.0022;      % units of mol/kg.....alkalinity
Voc    = 1.3*10^21;   % units of kg.........mass of the ocean
Vat    = 1.8*10^20;   % units of mol........mass of the atmosphere   
MT_CO2 = 0.002*Voc;   % units of mol C......mass of total carbon
MT_ALK = ALK*Voc;     % units of mol........mass of total alkalinity

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
    plot(trev, J(:,1))
    ylabel('Total mass of continental crust (kg)', 'fontsize', 14)
    xlabel('Time (years)', 'fontsize', 14)
    set(gca, 'XDir', 'reverse')
    set(gca, 'fontsize', 14)
    xlim([0 4.5*10^9])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    grid on
    
  % Plot #2, mass of sediments:
    subplot(1,3,2)
    plot(trev, J(:,2))
    ylabel('Total mass of sediments (kg)', 'fontsize', 14)
    xlabel('Time (years)', 'fontsize', 14)
    set(gca, 'XDir', 'reverse')
    set(gca, 'fontsize', 14)
    xlim([0 4.5*10^9])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    grid on
  
  % Set figure properties:
    set(gcf, 'color', 'w');