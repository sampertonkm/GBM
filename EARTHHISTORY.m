function dJdt = EARTHHISTORY(t,J)

dJdt = zeros(6,1);

lambda147 = 6.54*10^-12; % yr^-1, Sm-147 decay constant
tdel = 10^6; % timestep in years

P = 1;    % kg/yr, rate of crust production/mantle loss (arc basalts)
D = 0.01; % kg/yr, rate of crustal delamination/mantle gain

%Sm_mol = J(3)*(1/10^6)*(10^24)*(1/147);
%Sm_atm = Sm_mol*(6.022*10^23);
%Nd_atm = (Sm_atm)*(exp(lambda147*tdel)-1);
%Nd_mol = Nd_atm*(1/(6.022*10^23));
%Nd_ppm = Nd_mol*(142)*(10^6)*(1/(10^24));

%J(3) = (Sm_atm-Nd_atm)*(1/(6.022*10^23))*(147)*(10^6)*(1/(10^24)); % ppm Sm, mantle
%J(4) = Nd_ppm; % ppm Nd, mantle

% STEP 4. Solve the corresponding ODE:
dJdt(1) = -P+D;   % change in mantle mass
dJdt(2) = P-D;    % change in crust mass
dJdt(3) = -lambda147*J(3);
dJdt(4) = lambda147*J(3);
dJdt(5) = -lambda147*J(3);
dJdt(6) = lambda147*J(3);

%J(3) = (Sm_atm-Nd_atm)*(1/(6.022*10^23))*(147)*(10^6)*(1/(10^24)); % ppm Sm, mantle
%J(4) = Nd_ppm; % ppm Nd, mantle

end