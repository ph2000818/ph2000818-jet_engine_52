function smix = MixingEntropy(Y,Mi,Runiv)
%MIXINGENTROPY  Specific entropy of mixing of an ideal-gas mixture
%   In a mixture each species i sits at its partial pressure X_i*P, so its
%   entropy is s0_i(T) - R_i*ln(X_i*P/Pref), with R_i = Runiv/M_i. Summing
%   Y_i times this over all species gives the total mixture entropy
%       S = sum(Y_i*s0_i(T)) - Rg*ln(P/Pref) + smix
%   where smix = -sum(Y_i*R_i*ln(X_i)) >= 0 is returned by this function.
%   smix depends on the composition only, so it cancels in the isentropic
%   relations of the stages where the composition is fixed. It is needed
%   for the absolute entropies and the entropy generation in the combustor.
%   Self-contained: no globals, all dependencies are function arguments.
%
%   Input:
%     Y      - mass fraction vector (order matching Mi)
%     Mi     - molar masses [kg/mol]
%     Runiv  - universal gas constant [J/(mol*K)]
%
%   Output:
%     smix   - specific entropy of mixing [J/kg/K]

X  = (Y./Mi)/sum(Y./Mi);                                                   % Mole fractions: X_i = (Y_i/M_i)/sum(Y_j/M_j)
ip = Y > 0;                                                                % Absent species contribute nothing (Y*ln(X) -> 0)
smix = -sum(Y(ip).*(Runiv./Mi(ip)).*log(X(ip)));
end
