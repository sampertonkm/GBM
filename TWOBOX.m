clear; close all; clc

% Test a two-box model of crust-mantle evolution of the Sm-Nd system

time = 4.5*10^9;           % yr, model run time (age of Earth)

Mmantlei = 4.28*10^24;  % kg, initial mass of mantle (= whole system)
Mcrusti  = 0;           % kg, initial mass of crust (doesn't exist)
P = 1;                  % kg/yr, rate of crust production

Sm_mani = 10^10;                % ppm, Sm in mantle at t=0
Nd_mani = 0;                % ppm, Nd in mantle at t=0
Sm_crui = 0;                % ppm, Sm in crust at t=0
Nd_crui = 0;                % ppm, Nd in crust at t=0

lambda147 = 6.54*10^-12; % yr^-1, Sm-147 decay constant
Enrich   = 2;           % enrichment factor

% Set initial conditions array for the model
J0 = [Mmantlei Mcrusti Sm_mani Nd_mani Sm_crui Nd_crui];
tspan = 0:10^6:time;

[t,J] = ode15s('EARTHHISTORY',tspan,J0);
trev = fliplr(t');
tGa  = trev.*10^-9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MODEL VISUALIZATION

% Figure 1, Panel 1. Mass of continental crust vs. time:
subplot(2,3,1)
plot(tGa, J(:,2))
ylabel('Mass of continental crust (kg)', 'fontsize', 14)
xlabel('Time (Ga)', 'fontsize', 14)
xlim([0 max(tGa)])
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
set(gca, 'XDir', 'reverse')
set(gca, 'fontsize', 14)
set(gcf, 'color', 'w')
grid on

% Figure 1, Panel 2. Mass of mantle vs. time:
subplot(2,3,2)
plot(tGa, J(:,1))
ylabel('Mass of mantle (kg)', 'fontsize', 14)
xlabel('Time (Ga)', 'fontsize', 14)
xlim([0 max(tGa)])
ylim([min(J(:,1)) max(J(:,1))])
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
set(gca, 'XDir', 'reverse')
set(gca, 'fontsize', 14)
grid on

% Figure 1, Panel 3. Sm-147_mantle decay vs. time:
subplot(2,3,3)
plot(tGa, J(:,3))
ylabel('^{147}Sm_{mantle} (moles)', 'fontsize', 14)
xlabel('Time (Ga)', 'fontsize', 14)
xlim([0 max(tGa)])
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
set(gca, 'XDir', 'reverse')
set(gca, 'fontsize', 14)
grid on

% Figure 1, Panel 4. Nd-142_mantle growth vs. time:
subplot(2,3,4)
plot(tGa, J(:,4))
ylabel('^{143}Nd_{mantle} (moles)', 'fontsize', 14)
xlabel('Time (Ga)', 'fontsize', 14)
xlim([0 max(tGa)])
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
set(gca, 'XDir', 'reverse')
set(gca, 'fontsize', 14)
grid on

% Figure 1, Panel 5. Sm-147_crust decay vs. time:
subplot(2,3,5)
plot(tGa, J(:,5))
ylabel('^{147}Sm_{crust} (moles)', 'fontsize', 14)
xlabel('Time (Ga)', 'fontsize', 14)
xlim([0 max(tGa)])
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
set(gca, 'XDir', 'reverse')
set(gca, 'fontsize', 14)
grid on

% Figure 1, Panel 6. Sm-147_crust decay vs. time:
subplot(2,3,6)
plot(tGa, J(:,6))
ylabel('^{147}Sm_{crust} (moles)', 'fontsize', 14)
xlabel('Time (Ga)', 'fontsize', 14)
xlim([0 max(tGa)])
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
set(gca, 'XDir', 'reverse')
set(gca, 'fontsize', 14)
grid on
