function [T4,P4,Yprod,h3,h4] = Combustor(T3,P3,SpS,Yair,Yfuel,AF,dPloss,method)
%COMBUSTOR  Station [3-4]: constant-(approx.)pressure combustion of fuel in air
%   Unlike Diffusor/Compressor/Turbine/Nozzle, this stage is NOT isentropic
%   and the species composition changes (air+fuel -> combustion products).
%   T4 follows from an adiabatic energy balance between reactants (air at
%   T3 + fuel) and products (at T4), not from an entropy relation.
%
%   Input:
%     T3,P3   - inlet (compressor exit) temperature [K] / pressure [Pa]
%     SpS     - NASA species struct array, order {Fuel,O2,CO2,H2O,N2}
%               (matches iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'}) in Assignment.m)
%     Yair    - air mass fractions, order matching SpS
%     Yfuel   - fuel mass fractions, order matching SpS (e.g. [1 0 0 0 0])
%     AF      - air-to-fuel mass ratio
%     dPloss  - fractional pressure loss across combustor (0 if ignored)
%     method  - 'interp' or 'bisection', for solving T4 from the energy balance
%
%   Output:
%     T4,P4   - exit temperature/pressure
%     Yprod   - product mass fraction vector, order matching SpS
%     h3,h4   - inlet (reactants at T3) / exit (products at T4) specific
%               enthalpy [J/kg], on a per-kg-of-mixture basis
%
%   Self-contained: no globals, all dependencies are function arguments.

NSp = length(SpS);
Mi  = [SpS.Mass];

P4 = P3*(1-dPloss);

% Species order is {Fuel,O2,CO2,H2O,N2} (matches iSp in Assignment.m).
iF=1; iO2=2; iCO2=3; iH2O=4; iN2=5;

% Stoichiometry: fuel modeled as isooctane C8H18 (standard "Gasoline"
% surrogate). Adjust nC,nH below if your case uses a different fuel.
nC = 8; nH = 18;
nFuel      = 1/Mi(iF);                                                     % kmol fuel per kg fuel
nO2stoich  = nFuel*(nC+nH/4);                                              % kmol O2 needed for complete combustion
nCO2       = nFuel*nC;
nH2O       = nFuel*nH/2;

nO2in = AF*Yair(iO2)/Mi(iO2);                                              % kmol O2 supplied (per kg fuel)
nN2in = AF*Yair(iN2)/Mi(iN2);                                              % kmol N2 supplied (passes through)
nO2ex = nO2in-nO2stoich;                                                   % kmol excess O2 (unreacted)
if nO2ex < 0
    error('Combustor:RichMixture','AF is below the stoichiometric ratio; this model assumes complete (lean) combustion');
end

mtotal = 1+AF;                                                             % kg mixture per kg fuel
Yprod = zeros(1,NSp);
Yprod(iO2)  = nO2ex*Mi(iO2)/mtotal;
Yprod(iCO2) = nCO2*Mi(iCO2)/mtotal;
Yprod(iH2O) = nH2O*Mi(iH2O)/mtotal;
Yprod(iN2)  = nN2in*Mi(iN2)/mtotal;

% Reactants enthalpy (per kg of mixture = 1 kg fuel + AF kg air):
for i=1:NSp
    hi3(i) = HNasa(T3,SpS(i));
end
h3 = (Yfuel*1 + Yair*AF)*hi3'/mtotal;

% Adiabatic energy balance: h_products(T4) = h3
switch method
    case 'interp'
        TRl = 200:1:3000;
        for i=1:NSp
            hia(:,i) = HNasa(TRl,SpS(i));
        end
        hthermal_a = Yprod*hia';
        T4 = interp1(hthermal_a,TRl,h3);
    case 'bisection'
        TL = T3; TH = 3000;                                                % combustion only raises the temperature
        while abs(TH-TL) > 0.01
            Ti = (TL+TH)/2;
            for i=1:NSp
                hii(i) = HNasa(Ti,SpS(i));
            end
            hi = Yprod*hii';
            if hi > h3
                TH = Ti;
            else
                TL = Ti;
            end
        end
        T4 = (TH+TL)/2;
    otherwise
        error('Combustor:UnknownMethod','method must be ''interp'' or ''bisection''');
end

for i=1:NSp
    hi4(i) = HNasa(T4,SpS(i));
end
h4 = Yprod*hi4';
end
