function [T3,P3,h2,h3,S2,S3,Wc] = Compressor(T2,P2,v2,v3,SpS,Yair,PRc,method,Runiv,Pref)
%COMPRESSOR  Station [2-3]: adiabatic compression of air
%   Mirrors the Diffusor example in Assignment.m (interpolation/bisection
%   on the entropy relation), but here P3 is known (P3 = P2*PRc) and T3 is
%   the unknown solved from the isentropic assumption.
%   Self-contained: no globals, all dependencies are function arguments.
%
%   Input:
%     T2,P2   - inlet temperature [K] and pressure [Pa]
%     v2,v3   - inlet/exit velocity [m/s] (0 if kinetic energy is ignored)
%     SpS     - NASA species struct array for air, order matching Yair
%     Yair    - mass fraction vector for air (matches SpS order)
%     PRc     - compressor pressure ratio, P3/P2
%     method  - 'interp' or 'bisection' (see Assignment.m Diffusor example)
%     Runiv,Pref - universal gas constant [kJ/kmol/K], reference pressure [Pa]
%
%   Output:
%     T3,P3   - exit temperature/pressure
%     h2,h3   - inlet/exit specific enthalpy [J/kg]
%     S2,S3   - inlet/exit total specific entropy [J/kg/K]
%     Wc      - specific compressor work, Wc = h3-h2 [J/kg]

NSp = length(SpS);
Mi  = [SpS.Mass];
Rg  = Runiv/(Yair*Mi');

P3 = P2*PRc;

for i=1:NSp
    hi2(i) = HNasa(T2,SpS(i));
    si2(i) = SNasa(T2,SpS(i));
end
h2 = Yair*hi2';
s2thermal = Yair*si2';
S2 = s2thermal - Rg*log(P2/Pref);

% Energy balance across the compressor (fill in if v2/v3 are not negligible):
% h3 = h2 + 0.5*v2^2 - 0.5*v3^2

% Isentropic assumption: s3thermal(T3) - s2thermal(T2) = Rg*log(P3/P2)
starget = s2thermal + Rg*log(P3/P2);
switch method
    case 'interp'
        TRl = 200:1:3000;
        for i=1:NSp
            sia(:,i) = SNasa(TRl,SpS(i));
        end
        sthermal_a = Yair*sia';
        T3 = interp1(sthermal_a,TRl,starget);
    case 'bisection'
        TL = 200; TH = 3000;                                                % T3 must lie between these bounds
        while abs(TH-TL) > 0.01
            Ti = (TL+TH)/2;
            for i=1:NSp
                sii(i) = SNasa(Ti,SpS(i));
            end
            si = Yair*sii';
            if si > starget
                TH = Ti;
            else
                TL = Ti;
            end
        end
        T3 = (TH+TL)/2;
    otherwise
        error('Compressor:UnknownMethod','method must be ''interp'' or ''bisection''');
end

for i=1:NSp
    hi3(i) = HNasa(T3,SpS(i));
    si3(i) = SNasa(T3,SpS(i));
end
h3 = Yair*hi3';
s3thermal = Yair*si3';
S3 = s3thermal - Rg*log(P3/Pref);
Wc = h3-h2;
end
