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
%     Runiv,Pref - universal gas constant [J/(mol*K)], reference pressure [Pa]
%
%   Output:
%     T6,P6   - exit temperature/pressure (P6 = Pamb)
%     v6      - exit velocity [m/s]
%     h5,h6   - inlet/exit specific enthalpy [J/kg]
%     S5,S6   - inlet/exit total specific entropy [J/kg/K]

NSp = length(SpS);
Mi  = [SpS.Mass];
Rg  = Runiv*sum(Yprod./Mi);                                                % Mixture gas constant: 1/Mmix = sum(Y_i/M_i)

P6 = Pamb;

%% Debug: everything going INTO this stage
fprintf('\n[Nozzle 5-6] ---- inputs ----\n');
fprintf('  T5 = %9.4f K     P5 = %11.4f Pa    v5 = %8.4f m/s    Pamb(=P6) = %11.4f Pa    method = %s\n',T5,P5,v5,Pamb,method);
fprintf('  Runiv = %.6f J/mol/K   Pref = %.4f Pa   Rg = %.6f J/kg/K\n',Runiv,Pref,Rg);
fprintf('  Yprod : '); for i=1:NSp, fprintf('%s=%.4f  ',SpS(i).Name,Yprod(i)); end; fprintf('\n');

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

%% Debug: everything coming OUT of this stage
fprintf('[Nozzle 5-6] ---- outputs ----\n');
fprintf('  T6 = %9.4f K     P6 = %11.4f Pa    v6 = %9.4f m/s\n',T6,P6,v6);
fprintf('  h5 = %9.4f kJ/kg  h6 = %9.4f kJ/kg\n',h5/1e3,h6/1e3);
fprintf('  S5 = %9.4f kJ/kg/K  S6 = %9.4f kJ/kg/K\n',S5/1e3,S6/1e3);
fprintf('---------------------------------------------\n');
end
