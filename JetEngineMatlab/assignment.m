clear all;close all;clc;
warning off;
%% To make sure that matlab will find the functions. You must change it to your situation 
relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Nasa polynomials are loaded and globals are set. 
%% values should not be changed. These are used by all Nasa Functions. 
global Runiv Pref
Runiv=8.314472;
Pref=1e5;       % Reference pressure = standard-state pressure of the NASA polynomials (1 bar)
                % (template had 1.01235e5 "1 atm"; the s0 of H2,O2,CO2,H2O in the database match 1 bar data)
Tref=298.15;    % Reference Temperature
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% ===================== EXPORT SECTION =====================
%  All case-specific values live here and ONLY here. They are read from
%  the group's case file (values already in SI units) and exported as
%  named variables used throughout the rest of the script. To run a
%  different case, point caseFile at a different Groep0XX.txt - nothing
%  else in this script needs to change.
caseFile = 'Groep052.txt';
Case = ReadCaseData(caseFile);

cFuel    = Case.Fuel;                                                       % Fuel name (must match an entry in Sp.Name)
Tamb     = Case.Tamb;                                                       % Ambient temperature [K]
P3overP2 = Case.P3overP2;                                                   % Compressor pressure ratio [-]
Pamb     = Case.Pamb;                                                       % Ambient pressure [Pa]
mfurate  = Case.mfurate;                                                    % Fuel mass flow rate [kg/s]
AF       = Case.AF;                                                         % Air-to-fuel mass ratio [-]
v1       = Case.v1;                                                         % Inlet (flight) velocity [m/s]

% Fuel supply temperature: not in the case file, the fuel is assumed to be
% delivered to the combustor at ambient temperature
Tfuel    = Tamb;                                                            % Fuel temperature at combustor inlet [K]

% Fuel elemental composition CxHy, needed by Combustor.m for the
% combustion stoichiometry. Read from the fuel's own database entry
% (Sp.Elcomp, element order given by El) so it always matches the molar
% mass and polynomials used, e.g. 'Gasoline' is C7.76H13.1, not C8H18.
iFuel = find(strcmp({Sp.Name},cFuel),1);
if isempty(iFuel)
    error('Assignment:UnknownFuel','Fuel "%s" not found in the NASA database (Sp.Name)',cFuel);
end
cEl = {El.Name};
if any(Sp(iFuel).Elcomp(~ismember(cEl,{'C','H'})))
    error('Assignment:FuelNotCxHy','Combustor.m assumes a CxHy fuel, but "%s" contains other elements',cFuel);
end
nCfuel = Sp(iFuel).Elcomp(strcmp(cEl,'C'));                                 % C atoms per fuel molecule
nHfuel = Sp(iFuel).Elcomp(strcmp(cEl,'H'));                                 % H atoms per fuel molecule

% Print every exported value so it's clear at a glance what case is loaded.
fprintf('\nExported case data\n');
fprintf('-------------------------------------------------\n');
PrintVar('caseFile',caseFile);
PrintVar('cFuel',cFuel);
PrintVar('Tamb',Tamb,'K');
PrintVar('P3overP2',P3overP2,'-');
PrintVar('Pamb',Pamb/kPa,'kPa');
PrintVar('mfurate',mfurate,'kg/s');
PrintVar('AF',AF,'kg air/kg fuel');
PrintVar('v1',v1,'m/s');
PrintVar('Tfuel',Tfuel,'K');
PrintVar('nCfuel',nCfuel,'-');
PrintVar('nHfuel',nHfuel,'-');
fprintf('Constants\n');
PrintVar('Runiv',Runiv,'J/(mol K)');
PrintVar('Pref',Pref/kPa,'kPa');
PrintVar('Tref',Tref,'K');
fprintf('-------------------------------------------------\n\n');
%% ============================================================
%% Select species for the case at hand
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air composition
Xair = [0 0.21 0 0 0.79];                                                   % Order is important. Note that these are molefractions
MAir = Xair*Mi';                                                            % Row times Column = inner product 
Yair = Xair.*Mi/MAir;                                                       % Vector. times vector is Matlab's way of making an elementwise multiplication
sMixAir = MixingEntropy(Yair,Mi,Runiv);                                     % Entropy of mixing of air, -sum(Y_i*R_i*ln(X_i)) [J/kg/K]
%% Fuel composition
Yfuel = [1 0 0 0 0];                                                        % Only fuel
%% Range of enthalpies/thermal part of entropy of species
TR = [200:1:3000];NTR=length(TR);
for i=1:NSp                                                                 % Compute properties for all species for temperature range TR 
    hia(:,i) = HNasa(TR,SpS(i));                                            % hia is a NTR by 5 matrix
    sia(:,i) = SNasa(TR,SpS(i));                                            % sia is a NTR by 5 matrix
end
hair_a= Yair*hia';                                                          % Matlab 'inner product': 1x5 times 5xNTR matrix muliplication, 1xNTR resulT -> enthalpy of air for range of T 
sair_a= Yair*sia';                                                          % same but this thermal part of entropy of air for range of T
% whos hia sia hair_a sair_a                                                  % Shows dimensions of arrays on commandline
%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the interpolation method
% Bisection is in the next 'cell'
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using INTERPOLATION
cMethod = 'Interpolation Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
T2 = interp1(hair_a,TR,h2);                                                 % Interpolate h2 on h2air_a to approximate T2. Pretty accurate
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
h2check = Yair*hi2';                                                        % Single value (1x5 times 5x1). Why do I do compute this h2check value? Any ideas?
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral part of th eentropy)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref) + sMixAir;                                % Total specific entropy (incl. entropy of mixing)
S2  = s2thermal - Rg*log(P2/Pref) + sMixAir;
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
T2int = T2;

%% Two methods are presented to 'solve' the conservation equations for the Diffusor
%-------------------------------------------------------------------------
% ----> This part shows the Bisection method
%-------------------------------------------------------------------------
% [1-2] Diffusor :: Example approach using bisection (https://en.wikipedia.org/wiki/Bisection_method)
cMethod = 'Bisection Method';
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2 = 0;
h2 = h1+0.5*v1^2-0.5*v2^2;                                                  % Enhalpy at stage: h2 > h1 due to kinetic energy
TL = T1;
TH = 1000;                                                                  % A guess for the TH (must be too high)
iter = 0;
while abs(TH-TL) > 0.01
    iter = iter+1;
    Ti = (TL+TH)/2;
    for i=1:NSp
        hi2(i)    = HNasa(Ti,SpS(i));
    end
    h2i = Yair*hi2';                                                        % Single value (1x5 times 5x1). Intermediate value
    if h2i > h2
        TH = Ti; % new right boundary
    else
        TL = Ti; % new left boundary
    end
end
T2 = (TH+TL)/2;
T2bis = T2;
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;                                            % ln(P2/P1) = (s2-s1)/Rg , see lecture (s2 are only the temperature integral)
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref) + sMixAir;                                % Total entropy stage 1 (incl. entropy of mixing)
S2  = s2thermal - Rg*log(P2/Pref) + sMixAir;                                % Total entropy stage 2 (incl. entropy of mixing)
% Print to screen
fprintf('\n%14s\n',cMethod);
fprintf('Stage  ||%14s        [unit]\n      NR|%9i %9i\n',sPart,1,2);
fprintf('-------------------------------------\n');
fprintf('%8s| %9.2f %9.2f  [K]\n','Temp',T1,T2);
fprintf('%8s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%8s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
fprintf('---  H/S    -------------------------\n');
fprintf('%8s| %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ);
fprintf('%8s| %9.2f %9.2f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ);
%% Difference between two approaches: so close but not identical
fprintf('----------------------------------------------\n%8s| %9.4f %9.4f  [K]\n----------------------------------------------\n','T2-int vs T2-bis',T2int,T2bis);
%% Here starts your part (compressor,combustor,turbine and nozzle). ...
% Make a choice for which type of solution method you want to use.
method = 'bisection';                                                      % or 'interp'. Each stage below is its own file so you can work on them independently.

%% [2-3] Compressor
v3 = v2;                                                                    % Velocity inside the engine is negligible (as in Turns)
[T3,P3,h2,h3,S2,S3,Wc] = Compressor(T2,P2,v2,v3,SpS,Yair,P3overP2,method,Runiv,Pref);

%% [3-4] Combustor
dPloss = 0;                                                                 % TODO: set a pressure loss fraction if the case requires it
v4 = v3;                                                                    % No velocity change across the combustor
[T4,P4,Yprod,h3mix,h4] = Combustor(T3,P3,Tfuel,SpS,Yair,Yfuel,AF,dPloss,method,nCfuel,nHfuel);   % h3mix: air (T3) + fuel (Tfuel) reactants, not the compressor's h3

%% [4-5] Turbine
v5 = v4;                                                                    % Turbine.m's energy balance assumes no kinetic energy change (v5 = v4)
mfratio = (mfurate*AF+mfurate)/(mfurate*AF);                               % (air+fuel)/air mass flow ratio
[T5,P5,h4check,h5,S4,S5] = Turbine(T4,P4,SpS,Yprod,Wc,mfratio,method,Runiv,Pref);

%% [5-6] Nozzle
[T6,P6,v6,h5check,h6,S5check,S6] = Nozzle(T5,P5,v5,SpS,Yprod,Pamb,method,Runiv,Pref);

%% Engine performance
mair = mfurate*AF;                                                          % Air mass flow rate [kg/s]
mgas = mair+mfurate;                                                        % Exhaust gas mass flow rate [kg/s]
F    = mgas*v6 - mair*v1;                                                   % Thrust from momentum balance; pressure thrust (P6-Pamb)*A6 = 0 since P6 = Pamb
Fspec = F/mair;                                                             % Specific thrust [N/(kg/s)]
TSFC  = mfurate/F;                                                          % Thrust-specific fuel consumption [kg/(N s)]
% Heat of combustion straight from the NASA database: reactants minus
% products, both at Tref (no tabulated LHV needed)
for i=1:NSp
    hiref(i) = HNasa(Tref,SpS(i));
end
dHcomb = (Yfuel + AF*Yair)*hiref' - (1+AF)*Yprod*hiref';                   % [J/kg fuel]
Qin    = mfurate*dHcomb;                                                    % Chemical power released [W]
dEkin  = 0.5*mgas*v6^2 - 0.5*mair*v1^2;                                     % Kinetic energy gain of the flow [W]
Pprop  = F*v1;                                                              % Propulsive power [W]
eta_th   = dEkin/Qin;                                                       % Thermal efficiency
eta_prop = Pprop/dEkin;                                                     % Propulsive efficiency
eta_tot  = Pprop/Qin;                                                       % Overall efficiency = eta_th*eta_prop
% Entropy generated in the combustor, the only irreversible stage. Entropy
% balance (steady, adiabatic): Sgen = mgas*S4 - mair*S3 - mfurate*Sfuel.
% S3 (per kg air) and S4 (per kg gas) include the entropy of mixing; the
% fuel is pure and assumed to be injected at Tfuel and combustor pressure P3.
for i=1:NSp
    siF(i) = SNasa(Tfuel,SpS(i));
end
RgFuel = Runiv*sum(Yfuel./Mi);
Sfuel  = Yfuel*siF' - RgFuel*log(P3/Pref) + MixingEntropy(Yfuel,Mi,Runiv);  % [J/kg fuel/K] (mixing term is 0 for a pure fuel)
Sgen   = mgas*S4 - mair*S3 - mfurate*Sfuel;                                 % Entropy generation rate [W/K]

%% Print overview of all stations (one row per station)
cStation = {'1 inlet (ambient)','2 diffusor exit','3 compressor exit','4 combustor exit','5 turbine exit','6 nozzle exit'};
Tst = [T1 T2 T3 T4 T5 T6];   Pst = [P1 P2 P3 P4 P5 P6];   vst = [v1 v2 v3 v4 v5 v6];
hst = [h1 h2 h3 h4 h5 h6];   Sst = [S1 S2 S3 S4 S5 S6];
fprintf('\nFull cycle\n');
fprintf('%-18s %9s %9s %9s %10s %11s\n','Station','T [K]','P [kPa]','v [m/s]','h [kJ/kg]','S [kJ/kg/K]');
fprintf('---------------------------------------------------------------------\n');
for k=1:6
    fprintf('%-18s %9.2f %9.2f %9.2f %10.2f %11.4f\n',cStation{k},Tst(k),Pst(k)/kPa,vst(k),hst(k)/kJ,Sst(k)/kJ);
end
fprintf('---------------------------------------------------------------------\n');
fprintf('  (h,S at stations 1-3 are per kg air, 4-6 per kg combustion gas,\n');
fprintf('   so S3->S4 is NOT the combustor entropy generation: see Sgen below)\n');

%% Print engine performance
fprintf('\nPerformance\n');
fprintf('-------------------------------------------------\n');
PrintVar('mair',mair,'kg/s');
PrintVar('mgas',mgas,'kg/s');
PrintVar('Wc',Wc/kJ,'kJ/kg air');
PrintVar('dHcomb (LHV)',dHcomb/1e6,'MJ/kg fuel');
PrintVar('Thrust F',F/kN,'kN');
PrintVar('F/mair',Fspec,'N/(kg/s)');
PrintVar('TSFC',TSFC*1e6,'g/(kN s)');
PrintVar('eta_th',eta_th,'-');
PrintVar('eta_prop',eta_prop,'-');
PrintVar('eta_tot',eta_tot,'-');
PrintVar('Sgen combustor',Sgen/kJ,'kW/K');
PrintVar('Sgen/mair',Sgen/mair/kJ,'kJ/(kg air K)');
fprintf('-------------------------------------------------\n');
