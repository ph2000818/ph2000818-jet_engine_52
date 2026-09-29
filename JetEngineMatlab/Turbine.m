function [T5,P5,h4,h5,S4,S5] = Turbine(T4,P4,SpS,Yprod,Wc,mfratio,method,Runiv,Pref)
%TURBINE  Station [4-5]: adiabatic expansion of combustion gas
%   The turbine's only job is to drive the compressor, so the specific
%   turbine work is fixed by a power balance with the compressor. T5
%   follows from an energy balance, then P5 from the isentropic relation
%   (same structure as the Diffusor example, but here T5 is found first
%   and used to get P5, instead of the other way around).
%   Self-contained: no globals, all dependencies are function arguments.
%
%   Input:
%     T4,P4    - inlet (combustor exit) temperature [K] / pressure [Pa]
%     SpS      - NASA species struct array for the combustion products
%     Yprod    - product mass fractions, order matching SpS
%     Wc       - specific compressor work to be delivered [J/kg air],
%                from Compressor.m
%     mfratio  - ratio of gas mass flow to air mass flow, (mfa+mff)/mfa,
%                needed because turbine work is per kg of gas while Wc is
%                per kg of air (1 if this correction is ignored)
%     method   - 'interp' or 'bisection' (see Assignment.m Diffusor example)
%     Runiv,Pref - universal gas constant [J/(mol*K)], reference pressure [Pa]
%
%   Output:
%     T5,P5    - exit temperature/pressure
%     h4,h5    - inlet/exit specific enthalpy [J/kg]
%     S4,S5    - inlet/exit total specific entropy [J/kg/K]

NSp = length(SpS);
Mi  = [SpS.Mass];
Rg  = Runiv*sum(Yprod./Mi);                                                % Mixture gas constant: 1/Mmix = sum(Y_i/M_i)

%% Debug: everything going INTO this stage
fprintf('\n[Turbine 4-5] ---- inputs ----\n');
fprintf('  T4 = %9.4f K     P4 = %11.4f Pa    method = %s\n',T4,P4,method);
fprintf('  Wc = %9.4f kJ/kg air   mfratio(gas/air) = %.4f\n',Wc/1e3,mfratio);
fprintf('  Runiv = %.6f J/mol/K   Pref = %.4f Pa   Rg = %.6f J/kg/K\n',Runiv,Pref,Rg);
fprintf('  Yprod : '); for i=1:NSp, fprintf('%s=%.4f  ',SpS(i).Name,Yprod(i)); end; fprintf('\n');

for i=1:NSp
    hi4(i) = HNasa(T4,SpS(i));
    si4(i) = SNasa(T4,SpS(i));
end
h4 = Yprod*hi4';
s4thermal = Yprod*si4';
S4 = s4thermal - Rg*log(P4/Pref);

% Power balance: turbine work (per kg gas) = compressor work (per kg air)
Wt = Wc/mfratio;
h5 = h4-Wt;
fprintf('  Wt = %9.4f kJ/kg gas   (Wc/mfratio)     h4 = %9.4f kJ/kg   h5(target) = %9.4f kJ/kg\n',Wt/1e3,h4/1e3,h5/1e3);

% Solve T5 from the energy balance: Yprod*HNasa(T5,SpS)' = h5
switch method
    case 'interp'
        TRl = 200:1:3000;
        for i=1:NSp
            hia(:,i) = HNasa(TRl,SpS(i));
        end
        hthermal_a = Yprod*hia';
        T5 = interp1(hthermal_a,TRl,h5);
    case 'bisection'
        TL = 200; TH = T4;                                                 % expansion cools the gas, so T5 < T4
        while abs(TH-TL) > 0.01
            Ti = (TL+TH)/2;
            for i=1:NSp
                hii(i) = HNasa(Ti,SpS(i));
            end
            hi = Yprod*hii';
            if hi > h5
                TH = Ti;
            else
                TL = Ti;
            end
        end
        T5 = (TH+TL)/2;
    otherwise
        error('Turbine:UnknownMethod','method must be ''interp'' or ''bisection''');
end

for i=1:NSp
    si5(i) = SNasa(T5,SpS(i));
end
s5thermal = Yprod*si5';
% Isentropic assumption: s5thermal - s4thermal = Rg*log(P5/P4)
lnPr = (s5thermal-s4thermal)/Rg;
P5 = P4*exp(lnPr);
S5 = s5thermal - Rg*log(P5/Pref);

%% Debug: everything coming OUT of this stage
fprintf('[Turbine 4-5] ---- outputs ----\n');
fprintf('  T5 = %9.4f K     P5 = %11.4f Pa\n',T5,P5);
fprintf('  h4 = %9.4f kJ/kg  h5 = %9.4f kJ/kg\n',h4/1e3,h5/1e3);
fprintf('  S4 = %9.4f kJ/kg/K  S5 = %9.4f kJ/kg/K\n',S4/1e3,S5/1e3);
fprintf('---------------------------------------------\n');
end
