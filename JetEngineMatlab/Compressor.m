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
%     Runiv,Pref - universal gas constant [J/(mol*K)], reference pressure [Pa]
%
%   Output:
%     T3,P3   - exit temperature/pressure
%     h2,h3   - inlet/exit specific enthalpy [J/kg]
%     S2,S3   - inlet/exit total specific entropy [J/kg/K]
%     Wc      - specific compressor work [J/kg air],
%               Wc = (h3-h2) + 0.5*(v3^2-v2^2)

NSp = length(SpS);
Mi  = [SpS.Mass];
Rg  = Runiv*sum(Yair./Mi);                                                 % Mixture gas constant: 1/Mmix = sum(Y_i/M_i)

%% Debug: everything going INTO this stage
fprintf('\n[Compressor 2-3] ---- inputs ----\n');
fprintf('  T2 = %9.4f K     P2 = %11.4f Pa    v2 = %8.4f m/s   v3 = %8.4f m/s\n',T2,P2,v2,v3);
fprintf('  PRc(P3/P2) = %8.4f     method = %s\n',PRc,method);
fprintf('  Runiv = %.6f J/mol/K   Pref = %.4f Pa   Rg = %.6f J/kg/K\n',Runiv,Pref,Rg);
fprintf('  Yair : ');
for i=1:NSp, fprintf('%s=%.4f  ',SpS(i).Name,Yair(i)); end
fprintf('\n');

P3 = P2*PRc;

for i=1:NSp
    hi2(i) = HNasa(T2,SpS(i));
    si2(i) = SNasa(T2,SpS(i));
end
h2 = Yair*hi2';
s2thermal = Yair*si2';
S2 = s2thermal - Rg*log(P2/Pref);

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
% Energy balance (adiabatic, work IN): h2 + 0.5*v2^2 + Wc = h3 + 0.5*v3^2
Wc = (h3-h2) + 0.5*(v3^2-v2^2);

%% Debug: everything coming OUT of this stage
fprintf('[Compressor 2-3] ---- outputs ----\n');
fprintf('  starget = %.6f kJ/kg/K   (s2thermal + Rg*ln(P3/P2))\n',starget/1e3);
fprintf('  T3 = %9.4f K     P3 = %11.4f Pa\n',T3,P3);
fprintf('  h2 = %9.4f kJ/kg  h3 = %9.4f kJ/kg  Wc = %9.4f kJ/kg\n',h2/1e3,h3/1e3,Wc/1e3);
fprintf('  S2 = %9.4f kJ/kg/K  S3 = %9.4f kJ/kg/K\n',S2/1e3,S3/1e3);
fprintf('---------------------------------------------\n');
end
