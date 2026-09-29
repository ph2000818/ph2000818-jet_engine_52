clear all;close all;clc;
warning off
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
Pref=1.01235e5; % Reference pressure, 1 atm!
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

% Fuel elemental composition CxHy, needed by Combustor.m for the
% combustion stoichiometry. Add a case here if you ever change fuel.
switch cFuel
    case 'Gasoline'
        nCfuel = 8;  nHfuel = 18;                                           % isooctane C8H18 surrogate
    case 'H2'
        nCfuel = 0;  nHfuel = 2;                                            % hydrogen
    otherwise
        error('Assignment:UnknownFuel','No combustion stoichiometry (nC,nH) defined for fuel "%s"',cFuel);
end

% Print every exported value so it's clear at a glance what case is loaded.
fprintf('\n%14s\n','Exported case data');
fprintf('-------------------------------------------------\n');
fprintf('%12s: %s\n','caseFile',caseFile);
fprintf('%12s: %s\n','cFuel',cFuel);
fprintf('%12s: %10.4f  [K]\n','Tamb',Tamb);
fprintf('%12s: %10.4f  [-]\n','P3overP2',P3overP2);
fprintf('%12s: %10.4f  [Pa]\n','Pamb',Pamb);
fprintf('%12s: %10.4f  [kg/s]\n','mfurate',mfurate);
fprintf('%12s: %10.4f  [-]\n','AF',AF);
fprintf('%12s: %10.4f  [m/s]\n','v1',v1);
fprintf('%12s: %10d  [-]\n','nCfuel',nCfuel);
fprintf('%12s: %10d  [-]\n','nHfuel',nHfuel);
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
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total specific entropy
S2  = s2thermal - Rg*log(P2/Pref);
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
S1  = s1thermal - Rg*log(P1/Pref);                                          % Total entropy stage 1
S2  = s2thermal - Rg*log(P2/Pref);                                          % Total entropy stage 2
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
[T4,P4,Yprod,h3mix,h4] = Combustor(T3,P3,SpS,Yair,Yfuel,AF,dPloss,method,nCfuel,nHfuel);   % h3mix: air+fuel mixture, not the compressor's h3

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

%% Print overview of all stations
fprintf('\n%14s\n','Full cycle');
fprintf('Stage  ||%9s %9s %9s %9s %9s %9s\n','1','2','3','4','5','6');
fprintf('-------------------------------------------------------------------\n');
fprintf('%8s| %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f  [K]\n','Temp',T1,T2,T3,T4,T5,T6);
fprintf('%8s| %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa,P3/kPa,P4/kPa,P5/kPa,P6/kPa);
fprintf('%8s| %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f  [m/s]\n','v',v1,v2,v3,v4,v5,v6);
fprintf('---  H/S    -------------------------------------------------------\n');
fprintf('%8s| %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f  [kJ/kg]\n','h',h1/kJ,h2/kJ,h3/kJ,h4/kJ,h5/kJ,h6/kJ);
fprintf('%8s| %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f  [kJ/kg/K]\n','Total S',S1/kJ,S2/kJ,S3/kJ,S4/kJ,S5/kJ,S6/kJ);
fprintf('  (h,S at stations 1-3 are per kg air, 4-6 per kg combustion gas)\n');

%% Print engine performance
fprintf('\n%14s\n','Performance');
fprintf('-------------------------------------------------\n');
fprintf('%12s: %10.4f  [kg/s]\n','mair',mair);
fprintf('%12s: %10.4f  [kg/s]\n','mgas',mgas);
fprintf('%12s: %10.4f  [kJ/kg air]\n','Wc',Wc/kJ);
fprintf('%12s: %10.4f  [MJ/kg fuel]\n','dHcomb',dHcomb/1e6);
fprintf('%12s: %10.4f  [kN]\n','Thrust F',F/kN);
fprintf('%12s: %10.4f  [N/(kg/s)]\n','F/mair',Fspec);
fprintf('%12s: %10.4f  [g/(kN s)]\n','TSFC',TSFC*1e6);
fprintf('%12s: %10.4f  [-]\n','eta_th',eta_th);
fprintf('%12s: %10.4f  [-]\n','eta_prop',eta_prop);
fprintf('%12s: %10.4f  [-]\n','eta_tot',eta_tot);
fprintf('-------------------------------------------------\n');
