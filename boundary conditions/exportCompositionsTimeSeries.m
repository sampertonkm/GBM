%% Load files
if ~exist('mcigncn1','var'); load mcigncn1; end
if ~exist('igncn1','var'); load igncn1; end
if ~exist('msed','var'); load mcsed; end
if ~exist('sed','var'); load sed; end

mcigncn1=findOIBs(mcigncn1);

%% Elements to simulate
ChemicalElements={'SiO2';'TiO2';'Al2O3';'Fe2O3';'Fe2O3T';'FeO';'FeOT';'MgO';'CaO';'Na2O';'K2O';'P2O5';'MnO';'H2O_Total';'La';'Ce';'Pr';'Nd';'Sm';'Eu';'Gd';'Tb';'Dy';'Ho';'Er';'Tm';'Yb';'Lu';'Li';'Be';'B';'C';'CO2';'F';'Cl';'Sc';'Ti';'V';'Cr';'Co';'Ni';'Cu';'Zn';'Ga';'Zr';'Os';'Rb';'Bi';'Hg';'Ba';'Y';'Pb';'Te';'Nb';'Sr87_Sr86';'Tl';'Pt';'Sn';'Cd';'As';'Pd';'Sr';'Se';'S';'Au';'Ta';'Mo';'U';'Cs';'Sb';'Ag';'W';'Th';'Re';'Hf';'Ir';};

% Select elements that are present in both Igneous and Sed datasets.
elements={};
for e=ChemicalElements'
    if sum(~isnan(igncn1.(e{:})))>5000
        if isfield(sed,e{:}) && sum(~isnan(sed.(e{:})))>1000
            elements=[elements; e];
        end
    end
end


%% Timeseries properties
agemin=0;
agemax=3850;
timestep=1;
nbins=(agemax-agemin)/timestep;


%% Cumulative averages for exposed continental crust

Xcc=struct;
SiRange=[43,80];
for elem = elements';
    test=~isnan(mcigncn1.(elem{:})) & mcigncn1.SiO2>SiRange(1) & mcigncn1.SiO2<SiRange(2) & mcigncn1.Elevation>-100 & ~mcigncn1.oibs;
    [c,m,e]=bincumulative(mcigncn1.Age(test),mcigncn1.(elem{:})(test),agemin,agemax,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins,agemax);
    Xcc.Age=c;
    Xcc.(elem{:})=m;
    Xcc.err.(elem{:})=e;
end
save Xcc Xcc

%% Binned average for basaltic input
Xbslt=struct;
SiRange=[43,56];
for elem = elements';
    test=~isnan(mcigncn1.(elem{:})) & mcigncn1.SiO2>SiRange(1) & mcigncn1.SiO2<SiRange(2) & mcigncn1.Elevation>-100 & ~mcigncn1.oibs;
    [c,m,e]=bin(mcigncn1.Age(test),mcigncn1.(elem{:})(test),agemin,agemax,length(mcigncn1.SiO2)./length(igncn1.SiO2),nbins,100/timestep);
    Xbslt.Age=c;
    Xbslt.(elem{:})=m;
    Xbslt.err.(elem{:})=e;
end
save Xbslt Xbslt

%% Binned average for sediments
Xsed=struct;
for elem = elements';
    test=~isnan(mcsed.(elem{:}));
    [c,m,e]=bin(mcsed.Age(test),mcsed.(elem{:})(test),agemin,agemax,length(mcsed.SiO2)./length(sed.SiO2),nbins,100/timestep);
    Xsed.Age=c;
    Xsed.(elem{:})=m;
    Xsed.err.(elem{:})=e;
end
save Xsed Xsed

