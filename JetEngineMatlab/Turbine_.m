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
%     S4,S5    - inlet/exit total specific entropy, incl. entropy of
%                mixing [J/kg/K]

NSp = length(SpS);
Mi  = [SpS.Mass];
Rg  = Runiv*sum(Yprod./Mi);                                                % Mixture gas constant: 1/Mmix = sum(Y_i/M_i)
smix = MixingEntropy(Yprod,Mi,Runiv);                                      % Entropy of mixing (constant: composition does not change)

%% Debug: everything going INTO this stage
fprintf('\n[Turbine 4-5] inputs\n');
PrintVar('T4',T4,'K');
PrintVar('P4',P4/1e3,'kPa');
PrintVar('Wc',Wc/1e3,'kJ/kg air');
PrintVar('mfratio (gas/air)',mfratio,'-');
PrintVar('method',method);
PrintVar('Rg (gas)',Rg,'J/(kg K)');

for i=1:NSp
    hi4(i) = HNasa(T4,SpS(i));
    si4(i) = SNasa(T4,SpS(i));
end
h4 = Yprod*hi4';
s4thermal = Yprod*si4';
S4 = s4thermal - Rg*log(P4/Pref) + smix;

% Power balance: turbine work (per kg gas) = compressor work (per kg air)
Wt = Wc/mfratio;
h5 = h4-Wt;
PrintVar('Wt = Wc/mfratio',Wt/1e3,'kJ/kg gas');

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
% Isentropic assumption (S5 = S4, smix cancels): s5thermal - s4thermal = Rg*log(P5/P4)
lnPr = (s5thermal-s4thermal)/Rg;
P5 = P4*exp(lnPr);
S5 = s5thermal - Rg*log(P5/Pref) + smix;

%% Debug: everything coming OUT of this stage
fprintf('[Turbine 4-5] outputs\n');
PrintVar('T5',T5,'K');
PrintVar('P5',P5/1e3,'kPa');
PrintVar('h4',h4/1e3,'kJ/kg');
PrintVar('h5 = h4 - Wt',h5/1e3,'kJ/kg');
PrintVar('S4',S4/1e3,'kJ/(kg K)');
PrintVar('S5',S5/1e3,'kJ/(kg K)');
fprintf('---------------------------------------------\n');
end
