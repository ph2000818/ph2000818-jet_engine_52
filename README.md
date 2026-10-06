# ph2000818-jet_engine_52

## Structure

`JetEngineMatlab/Assignment.m` runs the full turbojet cycle station by
station: [1-2] Diffusor, [2-3] Compressor, [3-4] Combustor, [4-5] Turbine,
[5-6] Nozzle. The Diffusor is a fully worked example (interpolation and
bisection methods, side by side). The other four stages are implemented
as separate, self-contained function files so they can be read, tested,
or modified independently:

- `Compressor.m` — isentropic compression, `P3` known, solves for `T3`
- `Combustor.m` — adiabatic combustion (fuel + air -> products), not
  isentropic; also derives the product composition (`Yprod`) from the
  fuel's elemental formula (CxHy) and the air-fuel ratio
- `Turbine.m` — expansion sized by a power balance against the compressor
  (`Wc`), solves for `T5` then `P5`
- `Nozzle.m` — isentropic expansion to ambient pressure, solves for `T6`
  and the exit velocity `v6`
- `MixingEntropy.m` — small helper returning the entropy of mixing of a
  gas mixture, used for the total entropies `S`
- `PrintVar.m` — small helper printing one quantity per line (name,
  value, unit), so all console output is a vertical, aligned list

Each of these takes every value it needs as a function argument (no
`global`s), so one person can work on `Turbine.m` without needing to know
what `Compressor.m` or `Assignment.m` set up. `General/` is untouched
course material (NASA polynomial helpers `HNasa`, `SNasa`, `CpNasa`,
`CvNasa`, `UNasa`, `myfind`) shared by every stage.

## Governing equations

`sthermal(T) = Y·SNasa(T,SpS)'` and `h(T) = Y·HNasa(T,SpS)'` are the
mixture entropy (temperature-dependent part only) and enthalpy, computed
from the NASA polynomials for whichever composition `Y` (`Yair`,
`Yprod`, ...) applies at that stage. `Rg = Runiv·sum(Y./Mi)` (i.e. `Runiv/Mmix` with `1/Mmix = sum(Y_i/M_i)`) is the mixture's
specific gas constant. Every stage's total specific entropy is
`S = sthermal(T) - Rg*ln(P/Pref) + smix`, where
`smix = -sum(Y_i*(Runiv/M_i)*ln(X_i))` is the entropy of mixing
(`MixingEntropy.m`; each species sits at its partial pressure `X_i*P`).
`smix` only depends on the composition, so it cancels in the isentropic
relations below. `Pref = 1 bar`, the standard-state pressure of the NASA
polynomials in the database.

**Compressor [2-3]** — isentropic, `P3` known, solve for `T3`:
```
P3 = P2 * (P3/P2)
s3thermal(T3) - s2thermal(T2) = Rg * ln(P3/P2)      -> solve for T3
Wc = (h3 - h2) + 0.5*(v3^2 - v2^2)                   (specific compressor work, per kg air)
```

**Combustor [3-4]** — adiabatic, not isentropic; composition changes:
```
Stoichiometry (fuel CxHy):  CxHy + (x + y/4) O2  ->  x CO2 + (y/2) H2O
  x, y      = C and H atoms of the fuel, read from its database entry (Sp.Elcomp)
  nFuel     = 1 / Mfuel                             [mol fuel / kg fuel] (Mfuel in kg/mol)
  nO2stoich = nFuel * (x + y/4)
  nCO2      = nFuel * x
  nH2O      = nFuel * y/2
  nO2ex     = AF*Yair(O2)/M(O2) - nO2stoich          (must be >= 0: lean mixture)
  nN2       = AF*Yair(N2)/M(N2)                      (passes through unreacted)
  Yprod     = [0, nO2ex*M(O2), nCO2*M(CO2), nH2O*M(H2O), nN2*M(N2)] / (1+AF)

Energy balance (adiabatic flame temperature):
  h3mix = (1·Yfuel·h(Tfuel) + AF·Yair·h(T3))/(1+AF)  (air at T3 + fuel at Tfuel, NOT the compressor's h3)
  Yprod · h(T4) = h3mix                              -> solve for T4
```

**Turbine [4-5]** — sized to exactly drive the compressor, then isentropic:
```
Wt = Wc / mfratio                (mfratio = (air+fuel)/air mass flow)
h5 = h4 - Wt                                         -> solve T5 from Yprod·h(T5) = h5
s5thermal(T5) - s4thermal(T4) = Rg * ln(P5/P4)       -> solve for P5
```

**Nozzle [5-6]** — isentropic expansion to ambient pressure:
```
P6 = Pamb
s6thermal(T6) - s5thermal(T5) = Rg * ln(P6/P5)       -> solve for T6
v6 = sqrt(v5^2 + 2*(h5-h6))                          (energy balance)
```
Velocities inside the engine are negligible, as in Turns: `v2 = v3 = v4
= v5 = 0`. `Turbine.m` assumes no kinetic energy change, so `v5 = v4`.

**Performance** — computed at the end of `Assignment.m`:
```
mair = mfurate*AF,   mgas = mair + mfurate
F    = mgas*v6 - mair*v1                             (thrust; (P6-Pamb)*A6 = 0)
TSFC = mfurate / F
dHcomb = [(Yfuel + AF*Yair) - (1+AF)*Yprod] · h(Tref)   (per kg fuel, from NASA data)
eta_th   = (0.5*mgas*v6^2 - 0.5*mair*v1^2) / (mfurate*dHcomb)
eta_prop = F*v1 / (0.5*mgas*v6^2 - 0.5*mair*v1^2)
eta_tot  = eta_th * eta_prop
Sgen  = mgas*S4 - mair*S3 - mfurate*Sfuel             (combustor entropy generation)
Sfuel = sthermal_fuel(Tfuel) - Rfuel*ln(P3/Pref)     (pure fuel injected at P3)
```
The combustor is the only irreversible stage. The jump in the `S` column
from station 3 to 4 is not its entropy generation, because `S3` is per kg
air and `S4` per kg gas and the fuel's entropy is not in `S3`; use `Sgen`.
The final table prints one row per station (T, P, v, h and S), followed
by the performance block. h and S at stations 1-3 are per kg air, at 4-6
per kg combustion gas. Each stage also prints its inputs and outputs as
a vertical list; the combustor adds the stoichiometry, the equivalence
ratio and a mass-fraction table (air, fuel, reactants, products).

Each "solve for T" step is done with either `interp1` on a precomputed
h(T) or s(T) table, or bisection — selected by the `method` variable in
`Assignment.m`, exactly mirroring the Diffusor's worked example.

## Case data

The group's case values (fuel, `Tamb`, `P3overP2`, `Pamb`, `mfurate`,
`AF`, `v1`) live in `JetEngineMatlab/Groep052.txt`, not hardcoded in the
script. `ReadCaseData.m` parses that file into a struct, and
`Assignment.m` has a single **export section** near the top that reads
it and exports every case variable used by the rest of the script. It
also sets the fuel supply temperature `Tfuel = Tamb` (not in the case
file) and reads the fuel's `nC`/`nH` composition, needed by
`Combustor.m`, from the fuel's NASA database entry (`Sp.Elcomp`), so it
always matches the molar mass and polynomials used.
Running `Assignment.m` prints all of these exported values to the
console first, so it's clear at a glance which case is loaded. To run a
different group's case, change `caseFile` in that section — nothing
else needs to change.
