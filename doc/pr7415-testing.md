# PR 7415 phase-rate test branches

Use branch `codex/pr7415-phase-rate-consistency` in **both** `opm-simulators`
and `opm-common`. The simulator branch includes PR 7415 at `a5bf860`; the
common branch includes both commits from PR 5376 at `75c78f603`.
Both branches start from the locally checked-out master revisions.

## Behavior to test

- Standard and multisegment wells compute connection free/solution splits
  through `RatioCalculator`. Producing connections retain their free phase
  fluxes before adding dissolved/vaporized components, including injector
  backflow. Injecting connections select exactly one mixture model, with
  unmixed defaults for single-phase and oil-water configurations.
- Tracers use these unallocated connection splits. Well efficiency is applied
  at the tracer call site, as before.
- Well-level sums are reduced once per accepted timestep. Connection gas
  solution output includes gas dissolved in water as well as gas in oil.
- For producers without crossflow, `WellState::report()` allocates the raw
  free/solution fractions to the accepted wellhead component total. This
  changes only the output object, never well primary variables, connection
  rates, or tracer transport. A zero accepted total gives a zero reported
  split, including for stopped wells with internal crossflow.
- Crossflow is detected from the gathered connection component rates, so an
  injecting connection on another MPI rank is visible to the reporting rank.
  A nonzero total with crossflow, mixed-sign contributions, invalid data, or
  an undefined fraction retains the raw split. Consequently exact well-level
  additivity is deliberately **not** guaranteed for those cases. This is a
  reporting allocation experiment, not a new wellbore transport model.
- The MSW flux model retains its existing lack of Rsw/Rvw mixing.

## Build and focused tests

For the existing parent directory containing the three OPM repositories and
`opm-tests`, with its build directory at `../build`:

```sh
cmake -S . -B ../build
cmake --build ../build --target flow_blackoil test_RatioCalculator test_WellRateAllocation test_wellstate test_Summary -j 4
ctest --test-dir ../build --output-on-failure -R '^(RatioCalculator|WellRateAllocation|wellstate|Summary)$'
```

The new tests exercise the three-phase crossflow overwrite, oil-water and
single-phase crossflow, injector backflow, gas-water mixing, invalid-ratio
fallback, and preservation of automatic derivatives. The reporting tests
cover total/fraction preservation, zero free gas, zero wellhead flow,
unsupported allocations, and connection-state preservation across STOP and
reopen.

`flow_blackoil` exercises the standard-well, MSW, and black-oil tracer paths
changed here. Build the `flow` target as well to test the general executable
with all configured physics variants.

## Suggested simulation comparisons

Run these from the same parent directory. Use separate output directories so
results from each setting remain available for comparison.

```sh
../build/opm-simulators/bin/flow_blackoil opm-tests/model1/MSW_MODEL_1.DATA \
  --newton-min-iterations=1 --solver-max-time-step-in-days=5 \
  --threads-per-process=1 --output-dir=/tmp/pr7415-msw-min1

../build/opm-simulators/bin/flow_blackoil opm-tests/model1/MSW_MODEL_1.DATA \
  --newton-min-iterations=2 --solver-max-time-step-in-days=5 \
  --threads-per-process=1 --output-dir=/tmp/pr7415-msw-min2

mpirun -np 2 ../build/opm-simulators/bin/flow_blackoil opm-tests/model1/MSW_MODEL_1.DATA \
  --newton-min-iterations=1 --solver-max-time-step-in-days=5 \
  --allow-distributed-wells=true --threads-per-process=1 \
  --output-dir=/tmp/pr7415-msw-mpi2

../build/opm-simulators/bin/flow_blackoil opm-tests/tracer/GAS_TRACER-01.DATA \
  --newton-min-iterations=1 --threads-per-process=1 \
  --output-dir=/tmp/pr7415-gas-tracer
```

Check `WGPRF`, `WGPRS`, `WOPRF`, `WOPRS`, their group/field aggregates and
cumulatives, and free/solution tracer results. For ordinary production check
that well-level free plus solution equals the accepted component total.
For connection output check the actual signed component balance, including
crossflow; do not impose nonnegative production on injecting connections.
Different Newton settings need not produce bitwise-identical simulations.

To force perforations onto different ranks in this 6-by-8-by-7 fully active
model, supply an explicit alternating-layer partition:

```sh
python3 - <<'PYTHON'
from pathlib import Path
Path('/tmp/pr7415-msw-layers.partition').write_text(
    ''.join(f'{(i // 48) % 2} {i} 0\n' for i in range(336)))
PYTHON

mpirun -np 2 ../build/opm-simulators/bin/flow_blackoil opm-tests/model1/MSW_MODEL_1.DATA \
  --newton-min-iterations=1 --solver-max-time-step-in-days=5 \
  --allow-distributed-wells=true --external-partition=/tmp/pr7415-msw-layers.partition \
  --threads-per-process=1 --output-dir=/tmp/pr7415-msw-split-wells
```

## Validation performed

- Built `flow_blackoil` and all four test targets listed above.
- Passed `RatioCalculator` (7 cases), `WellRateAllocation` (4 cases),
  `wellstate`, and `Summary`.
- Completed `MSW_MODEL_1` with minimum Newton iterations 1 and 2, the
  default two-rank partition, a simple two-rank partition, and the explicit
  partition above. The explicit partition put five wells, including PROD-3,
  on both ranks.
- Checked 948 well gas/oil summary samples in each MSW run: no nonfinite
  values and no negative free-production rates. The maximum relative
  free-plus-solution mismatch was about 0.24% for the one-iteration serial
  case. The largest mismatches occur in PROD-3's raw reporting path; restart
  connection rates confirm this well has genuine crossflow. This is the
  documented fallback policy, not a guarantee of additivity for crossflow.
- Completed `GAS_TRACER-01` with minimum Newton iterations 1. This is a
  smoke test; no claim of full tracer accuracy or regression-reference
  equivalence is made.

The general `flow` executable and the full regression suite were not rebuilt
or run as part of this focused validation.
