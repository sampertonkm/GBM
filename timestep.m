function dJdt = timestep(time,J)

% STEP 1. Set system constants:

% STEP 2. Preallocate a dJdt matrix:
dJdt = zeros(5,1);

% STEP 3. Set fluxes for the constant input model:
Cbur = J(4);   % set carbonate burial = silicate weathering

% STEP 4. Solve the corresponding ODE:
dJdt(1) = 1;    % change in Mcc balance
dJdt(2) = 0.01*dJdt(1);    % change in Msed balance
dJdt(3) = 0;       % change in volcanic outgassing flux
dJdt(4) = 0;       % change in silicate weathering flux
dJdt(5) = 0;       % change in carbonate burial flux

end