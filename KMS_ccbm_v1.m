% KMS_ccbm_v1.m
% Script by K.M. Samperton, Dept. of Geosciences, Princeton Univ.
% Created on 2014-05-20, TUE, 08:30, Princeton NJ USA
% Updated on 2014-05-26, MON, 12:45, Princeton NJ USA

% GOALS: 1) Quantify the mass/compositional/isotopic evolution of the 
%           Earth system from 4.5 to 0 Ga using a box-modeling approach.
%        2) Constrain five global reservoirs through time, including:
%           a. Mantle................Reservoir #1
%           b. Continental crust.....Reservoir #2
%           c. Oceanic crust.........Reservoir #3
%           d. Seawater-atmosphere...Reservoir #4
%           e. Sediments.............Reservoir #5
%        3) The radiogenic isotope systems that will be investigated 
%           initially include:
%           a. Sm147-Nd143   (t1/2 = 106  Gyr)
%           b. Rb087-Sr087   (t1/2 = 48.8 Gyr)
%           c. Lu176-Hf176   (t1/2 = 35.4 Gyr)
%           d. U238-Pb206    (t1/2 = 4.47 Gyr)
%           e. U235-Pb207    (t1/2 = 704  Myr)
%           f. Th232-Pb208   (t1/2 = 14.0 Gyr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART I. ESTABLISH BOUNDARY CONDITIONS AND RESERVOIR CHARACTERISTICS

    % Step 1. Time model, clear workspace, close windows, and clc:
        tic; clear; close all; clc;

    % Step 2. Set global variables to be used by differential equations:
        global EarthAge AvNum lambda_Sm147 lambda_Rb087 lambda_Lu176...
               lambda_U238 lambda_U235 lambda_Th232 Prodc Prodo Delam ...
               Subd SilW Serp Sed Anatx CONDITIONS

        global e12Sm e12Nd e21Sm e21Nd e13Sm e13Nd e31Sm e31Nd e24Sm ...
               e24Nd e43Sm e43Nd e45Sm e45Nd e52Sm e52Nd ...
               e12Rb e12Sr e21Rb e21Sr e13Rb e13Sr e31Rb e31Sr e24Rb ...
               e24Sr e43Rb e43Sr e45Rb e45Sr e52Rb e52Sr ...
               e12Lu e12Hf e21Lu e21Hf e13Lu e13Hf e31Lu e31Hf e24Lu ...
               e24Hf e43Lu e43Hf e45Lu e45Hf e52Lu e52Hf ...
               e12U e12Th e12Pb e21U e21Th e21Pb e13U e13Th e13Pb ...
               e31U e31Th e31Pb e24U e24Th e24Pb e43U e43Th e43Pb ...
               e45U e45Th e45Pb e52U e52Th e52Pb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 3. Set system constants and boundary conditions, including
    % decay constants, initial reservoir masses and mixing terms:
        EarthAge     = 4.5*10^9;       % [yr]: model run time/age of Earth
        AvNum = 6.0221413*10^23;       % [atoms/mol]: Avogadro constant

        lambda_Sm147 = 6.54*10^-12;    % [yr^-1]: Sm-147 decay constant
        lambda_Rb087 = 1.42*10^-11;    % [yr^-1]: Rb-087 decay constant
        lambda_Lu176 = 1.96*10^-11;    % [yr^-1]: Lu-176 decay constant
        lambda_U238  = 1.55125*10^-10; % [yr^-1]: U-238 decay constant
        lambda_U235  = 9.8485*10^-10;  % [yr^-1]: U-235 decay constant    
        lambda_Th232 = 4.9475*10^-11;  % [yr^-1]: Th-232 decay constant

        Mmanti  = 4.28*10^24; % [kg]: initial mass of mantle
        Mcontci = 1;          % [kg]: initial mass of continental crust
        Moceaci = 1;          % [kg]: initial mass of oceanic crust
        Mseawi  = 1.30*10^21; % [kg]: initial mass of the ocean
        Matmoi  = 1.80*10^20; % [mol]: initial mass of the atmosphere
        Msedi   = 1;          % [kg]: initial mass of sediments

        Prodc  = 1;     % [kg/yr]: rate of continental crust production
        Prodo  = 1;     % [kg/yr]: rate of oceanic production
        Delam  = 0.1;   % [kg/yr]: rate of crustal delamination
        Subd   = 0.1;   % [kg/yr]: rate of subduction
        SilW   = 0.01;  % [kg/yr]: rate of silicate weathering
        Serp   = 0.001; % [kg/yr]: rate of seawater-crust alteration
        Sed    = 0.005; % [kg/yr]: rate of crustal sedimentation
        Anatx  = 0.001; % [kg/yr]: rate of sediment melting (anatexis)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % Step 4. Set initial concentrations for each element in the mantle:

        % Sm147-Nd143
        Sm147_manti  = 1; % [ppm Sm-147]: in mantle at t=0
        Nd143_manti  = 0; % [ppm Nd-143]: in mantle at t=0
        Sm147_contci = 0; % [ppm Sm-147]: in continental crust at t=0
        Nd143_contci = 0; % [ppm Nd-143]: in continental crust at t=0
        Sm147_oceaci = 0; % [ppm Sm-147]: in oceanic crust at t=0
        Nd143_oceaci = 0; % [ppm Nd-143]: in oceanic crust at t=0
        Sm147_seawi  = 0; % [ppm Sm-147]: in the ocean at t=0
        Nd143_seawi  = 0; % [ppm Nd-143]: in the ocean at t=0
        Sm147_sedi   = 0; % [ppm Sm-147]: in sediments at t=0
        Nd143_sedi   = 0; % [ppm Nd-143]: in sediments at t=0
        % Rb087-Sr087
        Rb087_manti  = 1; % [ppm Rb-087]: in mantle at t=0
        Sr087_manti  = 0; % [ppm Sr-087]: in mantle at t=0
        Rb087_contci = 0; % [ppm Rb-087]: in continental crust at t=0
        Sr087_contci = 0; % [ppm Sr-087]: in continental crust at t=0
        Rb087_oceaci = 0; % [ppm Rb-087]: in oceanic crust at t=0
        Sr087_oceaci = 0; % [ppm Sr-087]: in oceanic crust at t=0
        Rb087_seawi  = 0; % [ppm Rb-087]: in the ocean at t=0
        Sr087_seawi  = 0; % [ppm Sr-087]: in the ocean at t=0
        Rb087_sedi   = 0; % [ppm Rb-087]: in sediments at t=0
        Sr087_sedi   = 0; % [ppm Sr-087]: in sediments at t=0
        % Lu176-Hf176
        Lu176_manti  = 1; % [ppm Lu-176]: in mantle at t=0
        Hf176_manti  = 0; % [ppm Hf-176]: in mantle at t=0
        Lu176_contci = 0; % [ppm Lu-176]: in continental crust at t=0
        Hf176_contci = 0; % [ppm Hf-176]: in continental crust at t=0
        Lu176_oceaci = 0; % [ppm Lu-176]: in oceanic crust at t=0
        Hf176_oceaci = 0; % [ppm Hf-176]: in oceanic crust at t=0
        Lu176_seawi  = 0; % [ppm Lu-176]: in the ocean at t=0
        Hf176_seawi  = 0; % [ppm Hf-176]: in the ocean at t=0
        Lu176_sedi   = 0; % [ppm Lu-176]: in sediments at t=0
        Hf176_sedi   = 0; % [ppm Hf-176]: in sediments at t=0
        % U238-Pb206
        U238_manti   = 1; % [ppm U-238]:  in mantle at t=0
        Pb206_manti  = 0; % [ppm Pb-206]: in mantle at t=0
        U238_contci  = 0; % [ppm U-238]:  in continental crust at t=0
        Pb206_contci = 0; % [ppm Pb-206]: in continental crust at t=0
        U238_oceaci  = 0; % [ppm U-238]:  in oceanic crust at t=0
        Pb206_oceaci = 0; % [ppm Pb-206]: in oceanic crust at t=0
        U238_seawi   = 0; % [ppm U-238]:  in the ocean at t=0
        Pb206_seawi  = 0; % [ppm Pb-206]: in the ocean at t=0
        U238_sedi    = 0; % [ppm U-238]:  in sediments at t=0
        Pb206_sedi   = 0; % [ppm Pb-206]: in sediments at t=0
        % U235-Pb207
        U235_manti   = 1; % [ppm U-235]:  in mantle at t=0
        Pb207_manti  = 0; % [ppm Pb-207]: in mantle at t=0
        U235_contci  = 0; % [ppm U-235]:  in continental crust at t=0
        Pb207_contci = 0; % [ppm Pb-207]: in continental crust at t=0
        U235_oceaci  = 0; % [ppm U-235]:  in oceanic crust at t=0
        Pb207_oceaci = 0; % [ppm Pb-207]: in oceanic crust at t=0
        U235_seawi   = 0; % [ppm U-235]:  in the ocean at t=0
        Pb207_seawi  = 0; % [ppm Pb-207]: in the ocean at t=0
        U235_sedi    = 0; % [ppm U-235]:  in sediments at t=0
        Pb207_sedi   = 0; % [ppm Pb-207]: in sediments at t=0
        % Th232-Pb208
        Th232_manti  = 1; % [ppm Th-232]: in mantle at t=0
        Pb208_manti  = 0; % [ppm Pb-208]: in mantle at t=0
        Th232_contci = 0; % [ppm Th-232]: in continental crust at t=0
        Pb208_contci = 0; % [ppm Pb-208]: in continental crust at t=0
        Th232_oceaci = 0; % [ppm Th-232]: in oceanic crust at t=0
        Pb208_oceaci = 0; % [ppm Pb-208]: in oceanic crust at t=0
        Th232_seawi  = 0; % [ppm Th-232]: in the ocean at t=0
        Pb208_seawi  = 0; % [ppm Pb-208]: in the ocean at t=0
        Th232_sedi   = 0; % [ppm Th-232]: in sediments at t=0
        Pb208_sedi   = 0; % [ppm Pb-208]: in sediments at t=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 5. Set enrichment factors for each element for each mixing term:
        
        % Sm-Nd
        e12Sm = 1; % Sm, mantle-to-continental crust (Prodc)
        e12Nd = 2; % Nd, mantle-to-continental crust (Prodc)
        e21Sm = 1; % Sm, continental crust-to-mantle (Delam)
        e21Nd = 0.5; % Nd, continental crust-to-mantle (Delam)
        e13Sm = 1; % Sm, mantle-to-oceanic crust (Prodo)
        e13Nd = 1.5; % Nd, mantle-to-oceanic crust (Prodo)
        e31Sm = 1; % Sm, oceanic crust-to-mantle (Subd)
        e31Nd = 1; % Nd, oceanic crust-to-mantle (Subd)
        e24Sm = 1; % Sm, continental crust-to-the ocean (SilW)
        e24Nd = 1; % Nd, continental crust-to-the ocean (SilW)
        e43Sm = 1; % Sm, the ocean-to-oceanic crust (Serp)
        e43Nd = 1; % Nd, the ocean-to-oceanic crust (Serp)
        e45Sm = 1; % Sm, the ocean-to-sediments (Sed)
        e45Nd = 1; % Nd, the ocean-to-sediments (Sed)
        e52Sm = 2; % Sm, sediments-to-continental crust (Anatx)
        e52Nd = 1; % Nd, sediments-to-continental crust (Anatx)
        % Rb-Sr
        e12Rb = 1; % Rb, mantle-to-continental crust (Prodc)
        e12Sr = 1; % Sr, mantle-to-continental crust (Prodc)
        e21Rb = 1; % Rb, continental crust-to-mantle (Delam)
        e21Sr = 1; % Sr, continental crust-to-mantle (Delam)
        e13Rb = 1; % Rb, mantle-to-oceanic crust (Prodo)
        e13Sr = 1; % Sr, mantle-to-oceanic crust (Prodo)
        e31Rb = 1; % Rb, oceanic crust-to-mantle (Subd)
        e31Sr = 1; % Sr, oceanic crust-to-mantle (Subd)
        e24Rb = 1; % Rb, continental crust-to-the ocean (SilW)
        e24Sr = 1; % Sr, continental crust-to-the ocean (SilW)
        e43Rb = 1; % Rb, the ocean-to-oceanic crust (Serp)
        e43Sr = 1; % Sr, the ocean-to-oceanic crust (Serp)
        e45Rb = 1; % Rb, the ocean-to-sediments (Sed)
        e45Sr = 1; % Sr, the ocean-to-sediments (Sed)
        e52Rb = 1; % Rb, sediments-to-continental crust (Anatx)
        e52Sr = 1; % Sr, sediments-to-continental crust (Anatx)
        % Lu-Hf
        e12Lu = 1; % Lu, mantle-to-continental crust (Prodc)
        e12Hf = 1; % Hf, mantle-to-continental crust (Prodc)
        e21Lu = 1; % Lu, continental crust-to-mantle (Delam)
        e21Hf = 1; % Hf, continental crust-to-mantle (Delam)
        e13Lu = 1; % Lu, mantle-to-oceanic crust (Prodo)
        e13Hf = 1; % Hf, mantle-to-oceanic crust (Prodo)
        e31Lu = 1; % Lu, oceanic crust-to-mantle (Subd)
        e31Hf = 1; % Hf, oceanic crust-to-mantle (Subd)
        e24Lu = 1; % Lu, continental crust-to-the ocean (SilW)
        e24Hf = 1; % Hf, continental crust-to-the ocean (SilW)
        e43Lu = 1; % Lu, the ocean-to-oceanic crust (Serp)
        e43Hf = 1; % Hf, the ocean-to-oceanic crust (Serp)
        e45Lu = 1; % Lu, the ocean-to-sediments (Sed)
        e45Hf = 1; % Hf, the ocean-to-sediments (Sed)
        e52Lu = 1; % Lu, sediments-to-continental crust (Anatx)
        e52Hf = 1; % Hf, sediments-to-continental crust (Anatx)
        % U-Th-Pb
        e12U  = 1; % U,  mantle-to-continental crust (Prodc)
        e12Th = 1; % Th, mantle-to-continental crust (Prodc)
        e12Pb = 1; % Pb, mantle-to-continental crust (Prodc)
        e21U  = 1; % U,  continental crust-to-mantle (Delam)
        e21Th = 1; % Th, continental crust-to-mantle (Delam)
        e21Pb = 1; % Pb, continental crust-to-mantle (Delam)
        e13U  = 1; % U,  mantle-to-oceanic crust (Prodo)
        e13Th = 1; % Th, mantle-to-oceanic crust (Prodo)
        e13Pb = 1; % Pb, mantle-to-oceanic crust (Prodo)
        e31U  = 1; % U,  oceanic crust-to-mantle (Subd)
        e31Th = 1; % Th, oceanic crust-to-mantle (Subd)
        e31Pb = 1; % Pb, oceanic crust-to-mantle (Subd)
        e24U  = 1; % U,  continental crust-to-the ocean (SilW)
        e24Th = 1; % Tg, continental crust-to-the ocean (SilW)
        e24Pb = 1; % Pb, continental crust-to-the ocean (SilW)
        e43U  = 1; % U,  the ocean-to-oceanic crust (Serp)
        e43Th = 1; % Th, the ocean-to-oceanic crust (Serp)
        e43Pb = 1; % Pb, the ocean-to-oceanic crust (Serp)
        e45U  = 1; % U,  the ocean-to-sediments (Sed)
        e45Th = 1; % Th, the ocean-to-sediments (Sed)
        e45Pb = 1; % Pb, the ocean-to-sediments (Sed)
        e52U  = 1; % U,  sediments-to-continental crust (Anatx)
        e52Th = 1; % Th, sediments-to-continental crust (Anatx)
        e52Pb = 1; % Pb, sediments-to-continental crust (Anatx)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART II. PERFORM VARIOUS CALCULATIONS/CONVERSIONS PRIOR TO MODEL RUN

    % Step 1. Convert initial concentrations in mantle from ppm to atoms:
        % Sm147-Nd143
        Sm147_manti_atoms  = Sm147_manti*(1/10^6)*(Mmanti)*(1/147)*AvNum;
        Nd143_manti_atoms  = Nd143_manti*(1/10^6)*(Mmanti)*(1/143)*AvNum;
        Sm147_contci_atoms = Sm147_contci*(1/10^6)*(Mmanti)*(1/147)*AvNum;
        Nd143_contci_atoms = Nd143_contci*(1/10^6)*(Mmanti)*(1/143)*AvNum;
        Sm147_oceaci_atoms = Sm147_oceaci*(1/10^6)*(Mmanti)*(1/147)*AvNum;
        Nd143_oceaci_atoms = Nd143_oceaci*(1/10^6)*(Mmanti)*(1/143)*AvNum;
        Sm147_seawi_atoms  = Sm147_seawi*(1/10^6)*(Mmanti)*(1/147)*AvNum;
        Nd143_seawi_atoms  = Nd143_seawi*(1/10^6)*(Mmanti)*(1/143)*AvNum;
        Sm147_sedi_atoms   = Sm147_sedi*(1/10^6)*(Mmanti)*(1/147)*AvNum;
        Nd143_sedi_atoms   = Nd143_sedi*(1/10^6)*(Mmanti)*(1/143)*AvNum;
        % Rb087-Sr087
        Rb087_manti_atoms  = Rb087_manti*(1/10^6)*(Mmanti)*(1/87)*AvNum;
        Sr087_manti_atoms  = Sr087_manti*(1/10^6)*(Mmanti)*(1/87)*AvNum;
        Rb087_contci_atoms = Rb087_contci*(1/10^6)*(Mmanti)*(1/87)*AvNum;
        Sr087_contci_atoms = Sr087_contci*(1/10^6)*(Mmanti)*(1/87)*AvNum;
        Rb087_oceaci_atoms = Rb087_oceaci*(1/10^6)*(Mmanti)*(1/87)*AvNum;
        Sr087_oceaci_atoms = Sr087_oceaci*(1/10^6)*(Mmanti)*(1/87)*AvNum;
        Rb087_seawi_atoms  = Rb087_seawi*(1/10^6)*(Mmanti)*(1/87)*AvNum;
        Sr087_seawi_atoms  = Sr087_seawi*(1/10^6)*(Mmanti)*(1/87)*AvNum;
        Rb087_sedi_atoms   = Rb087_sedi*(1/10^6)*(Mmanti)*(1/87)*AvNum;
        Sr087_sedi_atoms   = Sr087_sedi*(1/10^6)*(Mmanti)*(1/87)*AvNum;
        % Lu176-Hf176
        Lu176_manti_atoms  = Lu176_manti*(1/10^6)*(Mmanti)*(1/176)*AvNum;
        Hf176_manti_atoms  = Hf176_manti*(1/10^6)*(Mmanti)*(1/176)*AvNum;
        Lu176_contci_atoms = Lu176_contci*(1/10^6)*(Mmanti)*(1/176)*AvNum;
        Hf176_contci_atoms = Hf176_contci*(1/10^6)*(Mmanti)*(1/176)*AvNum;
        Lu176_oceaci_atoms = Lu176_oceaci*(1/10^6)*(Mmanti)*(1/176)*AvNum;
        Hf176_oceaci_atoms = Hf176_oceaci*(1/10^6)*(Mmanti)*(1/176)*AvNum;
        Lu176_seawi_atoms  = Lu176_seawi*(1/10^6)*(Mmanti)*(1/176)*AvNum;
        Hf176_seawi_atoms  = Hf176_seawi*(1/10^6)*(Mmanti)*(1/176)*AvNum;
        Lu176_sedi_atoms   = Lu176_sedi*(1/10^6)*(Mmanti)*(1/176)*AvNum;
        Hf176_sedi_atoms   = Hf176_sedi*(1/10^6)*(Mmanti)*(1/176)*AvNum;
        % U238-Pb206
        U238_manti_atoms   = U238_manti*(1/10^6)*(Mmanti)*(1/238)*AvNum;
        Pb206_manti_atoms  = Pb206_manti*(1/10^6)*(Mmanti)*(1/206)*AvNum;
        U238_contci_atoms  = U238_contci*(1/10^6)*(Mmanti)*(1/238)*AvNum;
        Pb206_contci_atoms = Pb206_contci*(1/10^6)*(Mmanti)*(1/206)*AvNum;
        U238_oceaci_atoms  = U238_oceaci*(1/10^6)*(Mmanti)*(1/238)*AvNum;
        Pb206_oceaci_atoms = Pb206_oceaci*(1/10^6)*(Mmanti)*(1/206)*AvNum;
        U238_seawi_atoms   = U238_seawi*(1/10^6)*(Mmanti)*(1/238)*AvNum;
        Pb206_seawi_atoms  = Pb206_seawi*(1/10^6)*(Mmanti)*(1/206)*AvNum;
        U238_sedi_atoms    = U238_sedi*(1/10^6)*(Mmanti)*(1/238)*AvNum;
        Pb206_sedi_atoms   = Pb206_sedi*(1/10^6)*(Mmanti)*(1/206)*AvNum;
        % U235-Pb207
        U235_manti_atoms   = U235_manti*(1/10^6)*(Mmanti)*(1/235)*AvNum;
        Pb207_manti_atoms  = Pb207_manti*(1/10^6)*(Mmanti)*(1/207)*AvNum;
        U235_contci_atoms  = U235_contci*(1/10^6)*(Mmanti)*(1/235)*AvNum;
        Pb207_contci_atoms = Pb207_contci*(1/10^6)*(Mmanti)*(1/207)*AvNum;
        U235_oceaci_atoms  = U235_oceaci*(1/10^6)*(Mmanti)*(1/235)*AvNum;
        Pb207_oceaci_atoms = Pb207_oceaci*(1/10^6)*(Mmanti)*(1/207)*AvNum;
        U235_seawi_atoms   = U235_seawi*(1/10^6)*(Mmanti)*(1/235)*AvNum;
        Pb207_seawi_atoms  = Pb207_seawi*(1/10^6)*(Mmanti)*(1/207)*AvNum;
        U235_sedi_atoms    = U235_sedi*(1/10^6)*(Mmanti)*(1/235)*AvNum;
        Pb207_sedi_atoms   = Pb207_sedi*(1/10^6)*(Mmanti)*(1/207)*AvNum;
        % Th232-Pb208
        Th232_manti_atoms  = Th232_manti*(1/10^6)*(Mmanti)*(1/232)*AvNum;
        Pb208_manti_atoms  = Pb208_manti*(1/10^6)*(Mmanti)*(1/208)*AvNum;
        Th232_contci_atoms = Th232_contci*(1/10^6)*(Mmanti)*(1/232)*AvNum;
        Pb208_contci_atoms = Pb208_contci*(1/10^6)*(Mmanti)*(1/208)*AvNum;
        Th232_oceaci_atoms = Th232_oceaci*(1/10^6)*(Mmanti)*(1/232)*AvNum;
        Pb208_oceaci_atoms = Pb208_oceaci*(1/10^6)*(Mmanti)*(1/208)*AvNum;
        Th232_seawi_atoms  = Th232_seawi*(1/10^6)*(Mmanti)*(1/232)*AvNum;
        Pb208_seawi_atoms  = Pb208_seawi*(1/10^6)*(Mmanti)*(1/208)*AvNum;
        Th232_sedi_atoms   = Th232_sedi*(1/10^6)*(Mmanti)*(1/232)*AvNum;
        Pb208_sedi_atoms   = Pb208_sedi*(1/10^6)*(Mmanti)*(1/208)*AvNum;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART III. RUN THE MODEL

    % Step 1. Create array of initial conditions and set time interval:
    CONDITIONS = [Mmanti Mcontci Moceaci Mseawi Msedi...
      Sm147_manti_atoms Nd143_manti_atoms Sm147_contci_atoms ...
      Nd143_contci_atoms Sm147_oceaci_atoms Nd143_oceaci_atoms ...
      Sm147_seawi_atoms Nd143_seawi_atoms Sm147_sedi_atoms ...
      Nd143_sedi_atoms ...
      Rb087_manti_atoms Sr087_manti_atoms Rb087_contci_atoms ...
      Sr087_contci_atoms Rb087_oceaci_atoms Sr087_oceaci_atoms ...
      Rb087_seawi_atoms Sr087_seawi_atoms Rb087_sedi_atoms ...
      Sr087_sedi_atoms ...
      Lu176_manti_atoms Hf176_manti_atoms Lu176_contci_atoms ...
      Hf176_contci_atoms Lu176_oceaci_atoms Hf176_oceaci_atoms ...
      Lu176_seawi_atoms Hf176_seawi_atoms Lu176_sedi_atoms ...
      Hf176_sedi_atoms ...
      U238_manti_atoms Pb206_manti_atoms U238_contci_atoms ...
      Pb206_contci_atoms U238_oceaci_atoms Pb206_oceaci_atoms ...
      U238_seawi_atoms Pb206_seawi_atoms U238_sedi_atoms ...
      Pb206_sedi_atoms ...
      U235_manti_atoms Pb207_manti_atoms U235_contci_atoms ...
      Pb207_contci_atoms U235_oceaci_atoms Pb207_oceaci_atoms ...
      U235_seawi_atoms Pb207_seawi_atoms U235_sedi_atoms ...
      Pb207_sedi_atoms ...
      Th232_manti_atoms Pb208_manti_atoms Th232_contci_atoms ...
      Pb208_contci_atoms Th232_oceaci_atoms Pb208_oceaci_atoms ...
      Th232_seawi_atoms Pb208_seawi_atoms Th232_sedi_atoms ...
      Pb208_sedi_atoms];

    tspan = 0:10^6:EarthAge;

    % Step 2. Run the model and convert time units from years to Gyr:
    [t,EH] = ode15s('EarthHistory1', tspan, CONDITIONS);
    trev = fliplr(t');
    tGa  = trev.*10^-9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PART IV. DATA VISUALIZATION OF MODEL OUTPUTS

    % STEP 1. Plot the model outputs:
    
    % Figure 1, Panel 1. Mass of mantle vs. time:
    figure; set(gcf, 'Position', [1 1 1500 700])
    subplot(2,5,1); plot(tGa, EH(:,1));
    ylabel('Mantle (kg)', 'fontsize', 14)
    xlabel('Time (Ga)', 'fontsize', 14); xlim([0 max(tGa)])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    set(gca, 'XDir', 'reverse'); set(gca, 'fontsize', 14); grid on
    set(gcf, 'color', 'w')
    
    % Figure 1, Panel 2. Mass of continental crust vs. time:
    subplot(2,5,2); plot(tGa, EH(:,2))
    ylabel('Continental crust (kg)', 'fontsize', 14)
    xlabel('Time (Ga)', 'fontsize', 14); xlim([0 max(tGa)])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    set(gca, 'XDir', 'reverse'); set(gca, 'fontsize', 14); grid on

    % Figure 1, Panel 3. Mass of oceanic crust vs. time:
    subplot(2,5,3); plot(tGa, EH(:,3))
    ylabel('Oceanic crust (kg)', 'fontsize', 14)
    xlabel('Time (Ga)', 'fontsize', 14); xlim([0 max(tGa)])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    set(gca, 'XDir', 'reverse'); set(gca, 'fontsize', 14); grid on

    % Figure 1, Panel 4. Mass of seawater vs. time:
    subplot(2,5,4); plot(tGa, EH(:,4))
    ylabel('Seawater (kg)', 'fontsize', 14)
    xlabel('Time (Ga)', 'fontsize', 14); xlim([0 max(tGa)])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    set(gca, 'XDir', 'reverse'); set(gca, 'fontsize', 14); grid on
    
    % Figure 1, Panel 5. Mass of sediments vs. time:
    subplot(2,5,5); plot(tGa, EH(:,5))
    ylabel('Sediments (kg)', 'fontsize', 14)
    xlabel('Time (Ga)', 'fontsize', 14); xlim([0 max(tGa)])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    set(gca, 'XDir', 'reverse'); set(gca, 'fontsize', 14); grid on
    
    % Figure 1.6) Mantle Sm-Nd, Rb-Sr, Lu-Hf, U-Th-Pb vs. time:
    subplot(2,5,6);
    plot(tGa, EH(:,6), 'b'); hold on
    plot(tGa, EH(:,7), 'b'); hold on
    plot(tGa, EH(:,16),'r'); hold on
    plot(tGa, EH(:,17),'r'); hold on
    plot(tGa, EH(:,26),'g'); hold on
    plot(tGa, EH(:,27),'g'); hold on
    plot(tGa, EH(:,36),'m'); hold on
    plot(tGa, EH(:,37),'m'); hold on
    plot(tGa, EH(:,46),'k'); hold on
    plot(tGa, EH(:,47),'k'); hold on
    plot(tGa, EH(:,56),'c'); hold on
    plot(tGa, EH(:,57),'c');
    ylabel('Mantle Sm-Nd, Rb-Sr, Lu-Hf, U-Th-Pb (atoms)', 'fontsize', 14)
    xlabel('Time (Ga)', 'fontsize', 14); xlim([0 max(tGa)])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    set(gca, 'XDir', 'reverse'); set(gca, 'fontsize', 14); grid on

    % Figure 1.7) Continental crust Sm-Nd, Rb-Sr, Lu-Hf, U-Th-Pb vs. time:
    subplot(2,5,7);
    plot(tGa, EH(:,8), 'b'); hold on
    plot(tGa, EH(:,9), 'b'); hold on
    plot(tGa, EH(:,18),'r'); hold on
    plot(tGa, EH(:,19),'r'); hold on
    plot(tGa, EH(:,28),'g'); hold on
    plot(tGa, EH(:,29),'g'); hold on
    plot(tGa, EH(:,38),'m'); hold on
    plot(tGa, EH(:,39),'m'); hold on
    plot(tGa, EH(:,48),'k'); hold on
    plot(tGa, EH(:,49),'k'); hold on
    plot(tGa, EH(:,58),'c'); hold on
    plot(tGa, EH(:,59),'c');
    ylabel('Continental crust Sm-Nd, Rb-Sr, Lu-Hf, U-Th-Pb (atoms)', ...
      'fontsize', 14)
    xlabel('Time (Ga)', 'fontsize', 14); xlim([0 max(tGa)])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    set(gca, 'XDir', 'reverse'); set(gca, 'fontsize', 14); grid on
    
    % Figure 1.8) Oceanic crust Sm-Nd, Rb-Sr, Lu-Hf, U-Th-Pb vs. time:
    subplot(2,5,8);
    plot(tGa, EH(:,10), 'b'); hold on
    plot(tGa, EH(:,11), 'b'); hold on
    plot(tGa, EH(:,20), 'r'); hold on
    plot(tGa, EH(:,21), 'r'); hold on
    plot(tGa, EH(:,30), 'g'); hold on
    plot(tGa, EH(:,31), 'g'); hold on
    plot(tGa, EH(:,40), 'm'); hold on
    plot(tGa, EH(:,41), 'm'); hold on
    plot(tGa, EH(:,50), 'k'); hold on
    plot(tGa, EH(:,51), 'k'); hold on
    plot(tGa, EH(:,60), 'c'); hold on
    plot(tGa, EH(:,61), 'c');
    ylabel('Oceanic crust Sm-Nd, Rb-Sr, Lu-Hf, U-Th-Pb (atoms)', ...
      'fontsize', 14)
    xlabel('Time (Ga)', 'fontsize', 14); xlim([0 max(tGa)])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    set(gca, 'XDir', 'reverse'); set(gca, 'fontsize', 14); grid on
    
    % Figure 1.9) Seawater Sm-Nd, Rb-Sr, Lu-Hf, U-Th-Pb vs. time:
    subplot(2,5,9);
    plot(tGa, EH(:,12), 'b'); hold on
    plot(tGa, EH(:,13), 'b'); hold on
    plot(tGa, EH(:,22), 'r'); hold on
    plot(tGa, EH(:,23), 'r'); hold on
    plot(tGa, EH(:,32), 'g'); hold on
    plot(tGa, EH(:,33), 'g'); hold on
    plot(tGa, EH(:,42), 'm'); hold on
    plot(tGa, EH(:,43), 'm'); hold on
    plot(tGa, EH(:,52), 'k'); hold on
    plot(tGa, EH(:,53), 'k'); hold on
    plot(tGa, EH(:,62), 'c'); hold on
    plot(tGa, EH(:,63), 'c');
    ylabel('Seawater Sm-Nd, Rb-Sr, Lu-Hf, U-Th-Pb (atoms)', ...
      'fontsize', 14)
    xlabel('Time (Ga)', 'fontsize', 14); xlim([0 max(tGa)])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    set(gca, 'XDir', 'reverse'); set(gca, 'fontsize', 14); grid on

    % Figure 1.10) Sediments Sm-Nd, Rb-Sr, Lu-Hf, U-Th-Pb vs. time:
    subplot(2,5,10);
    plot(tGa, EH(:,14), 'b'); hold on
    plot(tGa, EH(:,15), 'b'); hold on
    plot(tGa, EH(:,24), 'r'); hold on
    plot(tGa, EH(:,25), 'r'); hold on
    plot(tGa, EH(:,34), 'g'); hold on
    plot(tGa, EH(:,35), 'g'); hold on
    plot(tGa, EH(:,44), 'm'); hold on
    plot(tGa, EH(:,45), 'm'); hold on
    plot(tGa, EH(:,54), 'k'); hold on
    plot(tGa, EH(:,55), 'k'); hold on
    plot(tGa, EH(:,64), 'c'); hold on
    plot(tGa, EH(:,65), 'c');
    ylabel('Sediments Sm-Nd, Rb-Sr, Lu-Hf, U-Th-Pb (atoms)', ...
      'fontsize', 14)
    xlabel('Time (Ga)', 'fontsize', 14); xlim([0 max(tGa)])
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    set(gca, 'XDir', 'reverse'); set(gca, 'fontsize', 14); grid on

    % Plot Nd143/Sm147 and see how these are sensitive to enrichment
    % factors:
    figure; 
    plot(tGa, EH(:,7)./EH(:,6),   'k'); hold on
    plot(tGa, EH(:,9)./EH(:,8),   'c'); hold on
    plot(tGa, EH(:,11)./EH(:,10), 'b'); hold on
    plot(tGa, EH(:,13)./EH(:,12), 'r'); hold on
    plot(tGa, EH(:,15)./EH(:,14), 'g');
    set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on')
    set(gca, 'XDir', 'reverse'); set(gca, 'fontsize', 14); grid on
    set(gcf, 'color', 'w')
    legend('Mant', 'CC', 'OC', 'SW', 'Seds', 'Location', 'NorthWest')
    xlabel('Time (Ga)', 'fontsize', 14); xlim([0 max(tGa)])
    ylabel('^{143}Nd/^{147}Sm', 'fontsize', 14)

    toc