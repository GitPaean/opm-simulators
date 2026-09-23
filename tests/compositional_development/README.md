# Compositional development cases

These 25 small decks exercise compositional well controls, multiple completions,
CO2 streams, water flow, and EOS and equilibration regions. They are kept with
`opm-simulators` under `tests/compositional_development` on the
`comp-development-cases-support` branch. The original copies in
`compositional_cases/cases_0923` remain available. Each case uses `FULLIMP` in
`RUNSPEC` and writes unified summary and restart output. The `LICENSES` block
is present but commented out in `RUNSPEC`, as requested; Flow ignores the
active keyword with a warning.

| Folder | Cases | Design check |
| --- | --- | --- |
| `co2_injection` | Pure CO2 and 80/20 CO2/methane | `WELLSTRE` and `WINJGAS` set the injector stream; injector-cell CO2/methane fractions follow the requested mix. |
| `water_phase_wells` | Mobile water and water-rate limit | Water injection raises injector-side saturation; producer WRAT control caps water production. |
| `multiple_connection_wells` | Two-completion waterflood and ORAT, with first/second-completion controls | Both completions contribute; the two-completion state differs from each first-only control, and water shifts toward the second injector cell. |
| `rate_control_wells` | Producer ORAT, WRAT, GRAT, LRAT, LRAT limit, RESV; gas and water injectors | Producer and injector rates attain targets or documented BHP limits. `WVPR` checks the RESV target. |
| `multiple_eos_regions` | PR/SRK region map, PR-only control, SRK-only control | Region 1 gas density follows the PR control, region 2 follows the SRK control; pressure and rates respond. |
| `surface_eos_regions` | Two surface EOS regions and region-one control | Surface gas rate changes for the producer in region 2 while reservoir fields stay the same. |
| `multiple_equilibration_regions` | Two `EQLNUM` regions and region-one control | Region 2 retains about 25% CO2 at the bottom; region 1 remains nearly CO2-free. |

All injectors request `RATE` in `WCONINJE`. Gas injectors request 700 sm3/day,
water injectors normally request 0.4 sm3/day, and the dedicated water-injector
case requests 0.2 sm3/day. Injector BHP values are safety limits. The runner
checks final well rates against either the rate target or the active BHP limit.
The RESV case can have recovered nonlinear timestep cutbacks; the runner
requires the simulation to finish and all final physical/rate checks to pass.

The EOS region cases use PR with zero BIC in region 1 and SRK with BIC
`0.20, 0.15, 0.10` in region 2. Their maps put cells 16–30 in region 2;
paired controls use region 1 everywhere. The surface region pair uses the
same contrast for `EOSS`/`BICS`. Both inject an 80/20 CO2/methane stream so
binary interactions matter. The equilibration pair has separate `EQUIL` and
`ZMFVD` definitions for cells 1–10 and 11–20. These differences are checked
against independent expectations as well as the paired controls.

## Run

From this directory:

```bash
./run.sh /path/to/flow_comp /tmp/compositional-development-runs --require-observable-effects
```

The runner needs Python 3 and NumPy. It writes each simulator log and Eclipse
output in the output directory. `verify_rate_controls.py` checks producer and
injector controls. `verify_outputs.py` checks summary/restart content,
physical bounds, composition, and paired case behavior. The optional
`--require-observable-effects` flag makes missing EOSNUM, SURFNUM, or
multi-completion effects fail the run.

The compositional path accepts `SURFNUM` and `TABDIMS(NMEOSS) > 1` with normal
parser strictness. Other Flow models retain their existing keyword validation.

These checks establish sensible behavior for the designed cases. They do not
replace comparison with a trusted simulator or a broader conservation and
parallel test campaign.
