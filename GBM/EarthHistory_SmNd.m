function dJdt = EarthHistory_SmNd(~,EH_SmNd_SmNd)

% Step 1. Match global variables set by primary box model script:
global EarthAge AvNum lambda_Sm147 lambda_Rb087 lambda_Lu176 ...
       lambda_U238 lambda_U235 lambda_Th232 Prodc Prodo Delam Subd ...
       SilW Serp Sed Anatx CONDITIONS

global e12Sm e12Nd e21Sm e21Nd e13Sm e13Nd e31Sm e31Nd e24Sm e24Nd...
           e43Sm e43Nd e45Sm e45Nd e52Sm e52Nd
   
% STEP 2. Preallocate and solve the corresponding ODE:
dJdt = zeros(numel(CONDITIONS),1);

dJdt(1) = -Prodc-Prodo+Delam+Subd; % change in mass of mantle
dJdt(2) = -Delam-SilW+Prodc+Anatx; % change in mass of continental crust
dJdt(3) = -Subd+Prodo+Serp;        % change in mass of oceanic crust
dJdt(4) = -Sed-Serp+SilW;          % change in mass of the ocean
dJdt(5) = -Anatx+Sed;              % change in mass of sediments

dJdt(6)  = -lambda_Sm147*EH_SmNd(6) -((Prodc/EH_SmNd(1))*EH_SmNd(6)*e12Sm)...
            -((Prodo/EH_SmNd(1))*EH_SmNd(6)*e13Sm) +((Delam/EH_SmNd(2))*EH_SmNd(8)*e21Sm)...
            +((Subd/EH_SmNd(3))*EH_SmNd(10)*e31Sm); 
            % change in Sm-147, mantle
            
dJdt(7)  =  lambda_Sm147*EH_SmNd(6) -((Prodc/EH_SmNd(1))*EH_SmNd(7)*e12Nd)...
            -((Prodo/EH_SmNd(1))*EH_SmNd(7)*e13Nd) +((Delam/EH_SmNd(2))*EH_SmNd(9)*e21Nd)...
            +((Subd/EH_SmNd(3))*EH_SmNd(11)*e31Nd);
            % change in Nd-143, mantle
            
dJdt(8)  = -lambda_Sm147*EH_SmNd(8) -((SilW/EH_SmNd(2))*EH_SmNd(8)*e24Sm)...
            -((Delam/EH_SmNd(2))*EH_SmNd(8)*e21Sm) +((Prodc/EH_SmNd(1))*EH_SmNd(6)*e12Sm)...
            +((Anatx/EH_SmNd(5))*EH_SmNd(14)*e52Sm); 
            % change in Sm-147, continental crust
            
dJdt(9)  =  lambda_Sm147*EH_SmNd(8) -((SilW/EH_SmNd(2))*EH_SmNd(9)*e24Nd)...
            -((Delam/EH_SmNd(2))*EH_SmNd(9)*e21Nd) +((Prodc/EH_SmNd(1))*EH_SmNd(7)*e12Nd)...
            +((Anatx/EH_SmNd(5))*EH_SmNd(15)*e52Nd);
            % change in Nd-143, continental crust
            
dJdt(10) = -lambda_Sm147*EH_SmNd(10) -((Subd/EH_SmNd(3))*EH_SmNd(10)*e31Sm)...
            +((Prodo/EH_SmNd(1))*EH_SmNd(6)*e13Sm) +((Serp/EH_SmNd(4))*EH_SmNd(12)*e43Sm); 
            % change in Sm-147, oceanic crust
            
dJdt(11) =  lambda_Sm147*EH_SmNd(10) -((Subd/EH_SmNd(3))*EH_SmNd(11)*e31Nd)...
            +((Prodo/EH_SmNd(1))*EH_SmNd(7)*e13Nd) +((Serp/EH_SmNd(4))*EH_SmNd(13)*e43Nd); 
            % change in Nd-143, oceanic crust
            
dJdt(12) = -lambda_Sm147*EH_SmNd(12) -((Serp/EH_SmNd(4))*EH_SmNd(12)*e43Sm)...
            -((Sed/EH_SmNd(4))*EH_SmNd(12)*e45Sm) +((SilW/EH_SmNd(2))*EH_SmNd(8)*e24Sm);
            % change in Sm-147, the ocean
            
dJdt(13) =  lambda_Sm147*EH_SmNd(12) -((Serp/EH_SmNd(4))*EH_SmNd(13)*e43Nd)...
           -((Sed/EH_SmNd(4))*EH_SmNd(13)*e45Nd) +((SilW/EH_SmNd(2))*EH_SmNd(9)*e24Nd);
           % change in Nd-143, the ocean
           
dJdt(14) = -lambda_Sm147*EH_SmNd(14) -((Anatx/EH_SmNd(5))*EH_SmNd(14)*e52Sm)...
            +((Sed/EH_SmNd(4))*EH_SmNd(12)*e45Sm); 
            % change in Sm-147, sediments
            
dJdt(15) =  lambda_Sm147*EH_SmNd(14) -((Anatx/EH_SmNd(5))*EH_SmNd(15)*e52Nd)...
            +((Sed/EH_SmNd(4))*EH_SmNd(13)*e45Nd);
            % change in Nd-143, sediments

%dJdt(8)  = -lambda_Rb087*EH_SmNd(8);    % change in Rb-087, mantle
%dJdt(9)  =  lambda_Rb087*EH_SmNd(8);    % change in Sr-087, mantle
%dJdt(10) = -lambda_Lu176*EH_SmNd(10);   % change in Lu-176, mantle
%dJdt(11) =  lambda_Lu176*EH_SmNd(10);   % change in Hf-176, mantle
%dJdt(12) = -lambda_U238*EH_SmNd(12);    % change in U-238,  mantle
%dJdt(13) =  lambda_U238*EH_SmNd(12);    % change in Pb-206, mantle
%dJdt(14) = -lambda_U235*EH_SmNd(14);    % change in U-235,  mantle
%dJdt(15) =  lambda_U235*EH_SmNd(14);    % change in Pb-207, mantle
%dJdt(16) = -lambda_Th232*EH_SmNd(16);   % change in Th-232, mantle
%dJdt(17) =  lambda_Th232*EH_SmNd(16);   % change in Pb-208, mantle

% Miscellaneous lines/calculations/notes below (NOT USED):
%Sm_mol = J(3)*(1/10^6)*(10^24)*(1/147);
%Sm_atm = Sm_mol*(6.022*10^23);
%Nd_atm = (Sm_atm)*(exp(lambda147*tdel)-1);
%Nd_mol = Nd_atm*(1/(6.022*10^23));
%Nd_ppm = Nd_mol*(142)*(10^6)*(1/(10^24));
%J(3) = (Sm_atm-Nd_atm)*(1/(6.022*10^23))*(147)*...
%       (10^6)*(1/(10^24)); % ppm Sm, mantle
%J(4) = Nd_ppm; % ppm Nd, mantle

end