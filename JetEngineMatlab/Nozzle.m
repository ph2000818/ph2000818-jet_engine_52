function [T6,P6,v6,h5,h6,S5,S6] = Nozzle(T5,P5,v5,SpS,Yprod,Pamb,method,Runiv,Pref)
%NOZZLE  Station [5-6]: isentropic expansion of combustion gas to ambient pressure
%   Structurally the same as the Diffusor example in Assignment.m: P6 is
%   known (=Pamb), T6 is solved from the isentropic relation, and the exit
%   velocity follows from the energy balance (the reverse role of
%   T2/v2 in the Diffusor, where velocity was known and T2 was solved for).
%   Self-contained: no globals, all dependencies are function arguments.
%
%   Input:
%     T5,P5   - inlet (turbine exit) temperature [K] / pressure [Pa]
%     v5      - inlet velocity [m/s]
%     SpS     - NASA species struct array for the combustion products
%     Yprod   - product mass fractions, order matching SpS
%     Pamb    - ambient (exit) pressure [Pa]
%     method  - 'interp' or 'bisection' (see Assignment.m Diffusor example)
%     Runiv,Pref - universal gas constant [kJ/kmol/K], reference pressure [Pa]
%
%   Output:
%     T6,P6   - exit temperature/pressure (P6 = Pamb)
%     v6      - exit velocity [m/s]
%     h5,h6   - inlet/exit specific enthalpy [J/kg]
%     S5,S6   - inlet/exit total specific entropy [J/kg/K]

NSp = length(SpS);
Mi  = [SpS.Mass];
Rg  = Runiv/(Yprod*Mi');

P6 = Pamb;

for i=1:NSp
    hi5(i) = HNasa(T5,SpS(i));
    si5(i) = SNasa(T5,SpS(i));
end
h5 = Yprod*hi5';
s5thermal = Yprod*si5';
S5 = s5thermal - Rg*log(P5/Pref);

% Isentropic assumption: s6thermal(T6) - s5thermal(T5) = Rg*log(P6/P5)
starget = s5thermal + Rg*log(P6/P5);
switch method
    case 'interp'
        TRl = 200:1:3000;
        for i=1:NSp
            sia(:,i) = SNasa(TRl,SpS(i));
        end
        sthermal_a = Yprod*sia';
        T6 = interp1(sthermal_a,TRl,starget);
    case 'bisection'
        TL = 200; TH = T5;                                                 % expansion to lower P cools the gas, so T6 < T5
        while abs(TH-TL) > 0.01
            Ti = (TL+TH)/2;
            for i=1:NSp
                sii(i) = SNasa(Ti,SpS(i));
            end
            si = Yprod*sii';
            if si > starget
                TH = Ti;
            else
                TL = Ti;
            end
        end
        T6 = (TH+TL)/2;
    otherwise
        error('Nozzle:UnknownMethod','method must be ''interp'' or ''bisection''');
end

for i=1:NSp
    hi6(i) = HNasa(T6,SpS(i));
    si6(i) = SNasa(T6,SpS(i));
end
h6 = Yprod*hi6';
s6thermal = Yprod*si6';
S6 = s6thermal - Rg*log(P6/Pref);

% Energy balance: h5 + 0.5*v5^2 = h6 + 0.5*v6^2
v6 = sqrt(v5^2 + 2*(h5-h6));
end
