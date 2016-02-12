function dJdt = dxdt(time,J)

% STEP 1. Set system constants:
Voc    = 1.3*10^21; % mass of the ocean, units of kg
Vat    = 1.8*10^20; % mass of the atmosphere, units of mol
ALK    = 0.0022;    % units of mol/kg

% STEP 2. Determine K values (calculate for constant Tc for now):
[K1,K2,Kw,Ksp,KH] = equilK(25,35);

% STEP 4. Define new constants to assist in algebra to solve for Hplus:
A = (KH*J(1))/Vat;
B = (KH*Voc)/Vat;

% STEP 5. Define the coefficients for the polynomial equation terms 1-5:
term1 = -(B+1);
term2 = (-B*K1) - (ALK*(B+1));
term3 = -(B*K1*K2) - (B*ALK*K1) + (A*K1) + (Kw*(B+1));
term4 = (B*K1*Kw) - (B*ALK*K1*K2) + (2*A*K1*K2);
term5 = B*K1*K2*Kw;

% STEP 6. Assemble the CO2 speciation polynomial, solve for HC species:
Hpolynomial = [term1 term2 term3 term4 term5];
Hroots      = roots(Hpolynomial);
H           = Hroots(find(Hroots > 0));
pH          = -log10(H);
OH          = Kw/H;
CO3         = (ALK - OH + H)/((H/K2)+2);
HCO3        = CO3*H/K2;
CO2aq       = HCO3*H/K1;
CO2at_ppm   = CO2aq/KH*(10^6);

temp = 5*log(CO2at_ppm/345)+12;   % temperature, units of deg C
J(4) = 29e11*exp((temp-13)/12.5); % dependence of SilW on Tc
J(6) = temp;

% STEP 2. Preallocate a dJdt matrix:
dJdt = zeros(6,1);

% STEP 3. Set fluxes for the constant input model:
Cbur = J(4);   % set carbonate burial = silicate weathering

% STEP 4. Solve the corresponding ODE:
dJdt(1) = J(3) - J(5) - Cbur;    % MT_CO2 balance
dJdt(2) = (2*J(4)) - (2*Cbur);   % MT_ALK balance
dJdt(3) = 0;                     % Constant volcanic outgassing flux
dJdt(4) = 0;                     % Constant silicate weathering flux
dJdt(5) = 0;                     % Constant carbonate burial rate
dJdt(6) = 0;

%J(6) = 5*log(CO2at_ppm/280e-6)+12; % temperature, units of deg C

end