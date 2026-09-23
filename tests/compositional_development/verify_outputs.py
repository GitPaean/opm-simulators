#!/usr/bin/env python3
"""Smoke-check output files and report observable development-case differences."""

from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "tools"))
import eclio


DECK_DIR = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise SystemExit(message)


def check_case(output, deck):
    name = deck.stem
    log = output / (name + ".log")
    require(log.is_file(), f"{name}: missing simulator log")
    log_text = log.read_text()
    require("End of simulation" in log_text, f"{name}: simulation did not complete")
    issues = [line for line in log_text.splitlines()
              if ((line.startswith(("Problem:", "Error:"))
                   and "Solver convergence failure - Iteration limit reached" not in line)
                  or "unhandled for output to restart file" in line)]
    require(not issues, f"{name}: simulator reported: {issues[:3]}")
    recoveries = log_text.count("Solver convergence failure - Iteration limit reached")
    if recoveries:
        print(f"{name}: {recoveries} recovered timestep cutbacks")
    base = output / name / name
    for suffix in ("EGRID", "INIT", "SMSPEC", "UNSMRY", "UNRST"):
        path = Path(str(base) + "." + suffix)
        require(path.is_file() and path.stat().st_size > 0, f"{name}: missing {suffix}")

    labels, values = eclio.summary(str(base) + ".SMSPEC", str(base) + ".UNSMRY")
    require(len(values) > 1, f"{name}: summary has no report steps")
    vertical = name in ("EQLNUM_TWO_REGIONS", "EQLNUM_REGION_ONE_CONTROL")
    region_pair = name.startswith(("EOSNUM_", "SURFNUM_"))
    cells = (1, 10, 11, 20) if vertical else ((1, 15, 16, 30) if region_pair
                                           else (1, 15, 30))
    water = "\nWATER\n" in deck.read_text().split("\nGRID\n", 1)[0]
    summary_keys = {"FPR"}
    for prefix in ("BPR", "BGSAT") + (("BWSAT",) if water else ()):
        summary_keys.update(f"{prefix}:{cell}" for cell in cells)
    if not vertical:
        summary_keys.update(("WBHP:INJ", "WBHP:PROD", "WOPR:PROD", "WGPR:PROD"))
        if water:
            summary_keys.update(("FWIR", "WWIR:INJ"))
        else:
            summary_keys.update(("FGIR", "FGIT", "WGIR:INJ"))
    missing = summary_keys.difference(labels)
    require(not missing, f"{name}: missing summary vectors {sorted(missing)}")
    for key in summary_keys:
        require(np.isfinite(values[:, labels.index(key)]).all(),
                f"{name}: non-finite summary vector {key}")

    restart = eclio.steps(str(base) + ".UNRST")
    require(len(restart) > 1, f"{name}: restart has no later report steps")
    final = restart[max(restart)]
    restart_keys = {"PRESSURE", "SGAS", "SOIL", "DENG", "DENO",
                    "ZMF1", "ZMF2", "ZMF3",
                    "XMF1", "XMF2", "XMF3", "YMF1", "YMF2", "YMF3"}
    if water:
        restart_keys.update(("SWAT", "DENW"))
    else:
        restart_keys.add("VMF")
    missing = restart_keys.difference(final)
    require(not missing, f"{name}: missing restart arrays {sorted(missing)}")
    n_cells = 20 if vertical else 30
    for key in restart_keys:
        require(len(final[key]) == n_cells and np.isfinite(final[key]).all(),
                f"{name}: invalid restart array {key}")
    pressure = final["PRESSURE"]
    require(np.all(pressure > 0), f"{name}: nonpositive reservoir pressure")
    saturation = final["SGAS"] + final["SOIL"]
    if water:
        saturation = saturation + final["SWAT"]
    require(np.all(np.abs(saturation - 1.) < 2e-4),
            f"{name}: phase saturations do not sum to one")
    for key in ("SGAS", "SOIL") + (("SWAT",) if water else ()):
        require(np.all(final[key] >= -1e-6) and np.all(final[key] <= 1 + 1e-6),
                f"{name}: saturation {key} outside [0,1]")
    total_z = sum(final[f"ZMF{i}"] for i in (1, 2, 3))
    require(np.all(np.abs(total_z - 1.) < 2e-4),
            f"{name}: overall mole fractions do not sum to one")
    for i in (1, 2, 3):
        require(np.all(final[f"ZMF{i}"] >= -1e-6)
                and np.all(final[f"ZMF{i}"] <= 1 + 1e-6),
                f"{name}: ZMF{i} outside [0,1]")
    print(f"{name}: {len(summary_keys)} key summary vectors, "
          f"{len(restart_keys)} key restart arrays, {len(restart)} restart steps")
    return final


def check_stream_compositions(final_states):
    for name, co2, methane in (("CO2_PURE_RATE", 1.0, 0.0),
                               ("CO2_METHANE_80_20_RATE", 0.8, 0.2)):
        state = final_states[name]
        for keyword, target in (("ZMF1", co2), ("ZMF2", methane)):
            actual = float(state[keyword][0])
            require(abs(actual - target) < 0.01,
                    f"{name}: injector-cell {keyword} is {actual}, expected {target}")
        print(f"{name}: injector-cell CO2/methane mole fractions "
              f"{state['ZMF1'][0]:.3f}/{state['ZMF2'][0]:.3f}")



def max_summary_difference(output, two_regions, control):
    def read_summary(name):
        base = output / name / name
        return eclio.summary(str(base) + ".SMSPEC", str(base) + ".UNSMRY")

    labels_two, values_two = read_summary(two_regions)
    labels_one, values_one = read_summary(control)
    require(labels_two == labels_one and values_two.shape == values_one.shape,
            f"{two_regions}: region pair has mismatched summary vectors")
    keys = ("FPR", "FGIR", "FGIT", "FGPR", "FGPT", "WGIR:INJ",
            "WGPR:PROD", "WOPR:PROD")
    return max(float(np.max(np.abs(values_two[:, labels_two.index(key)] -
                                   values_one[:, labels_one.index(key)])))
               for key in keys if key in labels_two)


def report_region_effect(states, output):
    missing_effects = []
    pr = states["EOSNUM_REGION_ONE_CONTROL"]["PRESSURE"]
    srk = states["EOSNUM_SRK_PRIMARY_CONTROL"]["PRESSURE"]
    eos_parameter_effect = float(np.max(np.abs(pr - srk)))
    require(eos_parameter_effect > 0.5,
            "EOS positive control: tuned SRK properties do not change pressure")
    print(f"EOS positive control: tuned SRK properties change final pressure by "
          f"up to {eos_parameter_effect:.3f} bar when placed in the primary slot")
    two = states["EOSNUM_TWO_REGIONS"]
    base = states["EOSNUM_REGION_ONE_CONTROL"]
    srk = states["EOSNUM_SRK_PRIMARY_CONTROL"]
    # Region 1 should retain PR-like gas density, while region 2 should
    # acquire SRK-like density. Density is more local than pressure here.
    for cell, expected in ((5, base), (20, srk)):
        index = cell - 1
        density = float(two["DENG"][index])
        expected_density = float(expected["DENG"][index])
        opposite = srk if cell == 5 else base
        opposite_density = float(opposite["DENG"][index])
        require(abs(density - expected_density) < abs(density - opposite_density),
                f"EOSNUM: cell {cell} gas density does not follow assigned EOS region")
    print("EOSNUM: region 1/2 gas densities follow the PR/SRK controls")

    for label, two_regions, control in (
        ("EOSNUM", "EOSNUM_TWO_REGIONS", "EOSNUM_REGION_ONE_CONTROL"),
        ("SURFNUM", "SURFNUM_TWO_REGIONS", "SURFNUM_REGION_ONE_CONTROL"),
        ("EQLNUM", "EQLNUM_TWO_REGIONS", "EQLNUM_REGION_ONE_CONTROL"),
    ):
        fields = ("PRESSURE", "SGAS", "ZMF1", "ZMF2", "ZMF3")
        difference = max(float(np.max(np.abs(states[two_regions][key] -
                                             states[control][key]))) for key in fields)
        summary_difference = max_summary_difference(output, two_regions, control)
        if label == "EQLNUM":
            co2_difference = abs(float(states[two_regions]["ZMF1"][10] -
                                       states[control]["ZMF1"][10]))
            require(co2_difference > 0.1,
                    f"EQLNUM: region 2 has no observable CO2 effect at cell 11 "
                    f"(difference {co2_difference:.6g})")
            require(float(states[two_regions]["ZMF1"][0]) < 0.01
                    and abs(float(states[two_regions]["ZMF1"][19]) - 0.25) < 0.02,
                    "EQLNUM: final CO2 profile does not follow region targets")
            print(f"EQLNUM: observed CO2 mole-fraction difference "
                  f"{co2_difference:.3f} at cell 11 versus region-one control")
        elif label == "SURFNUM":
            require(difference < 1e-5,
                    "SURFNUM: changing surface EOS unexpectedly changed reservoir fields")
            base = output / two_regions / two_regions
            ctrl = output / control / control
            labels, values = eclio.summary(str(base) + ".SMSPEC",
                                           str(base) + ".UNSMRY")
            _, ctrl_values = eclio.summary(str(ctrl) + ".SMSPEC",
                                           str(ctrl) + ".UNSMRY")
            rate = float(values[-1, labels.index("WGPR:PROD")])
            ctrl_rate = float(ctrl_values[-1, labels.index("WGPR:PROD")])
            require(rate > ctrl_rate + 0.1,
                    "SURFNUM: region 2 should increase producer surface gas rate")
            print(f"SURFNUM: producer surface gas rate rises {ctrl_rate:.3f} to "
                  f"{rate:.3f} sm3/day at unchanged reservoir state")
        elif difference < 1e-6 and summary_difference < 1e-6:
            print(f"{label}: NO OBSERVABLE REGION EFFECT in restart or summary fields")
            missing_effects.append(label)
        else:
            print(f"{label}: region pair differs (restart {difference:.6g}, "
                  f"summary {summary_difference:.6g}); "
                  "with bounded physical fields")
    return missing_effects


def report_connection_effect(states, output):
    missing_effects = []
    def difference(multiple, single):
        fields = ("PRESSURE", "SWAT", "SGAS", "ZMF1", "ZMF2", "ZMF3")
        restart = max(float(np.max(np.abs(states[multiple][key] -
                                          states[single][key]))) for key in fields)
        first = output / multiple / multiple
        second = output / single / single
        first_labels, first_values = eclio.summary(str(first) + ".SMSPEC",
                                                    str(first) + ".UNSMRY")
        second_labels, second_values = eclio.summary(str(second) + ".SMSPEC",
                                                      str(second) + ".UNSMRY")
        require(first_labels == second_labels and first_values.shape == second_values.shape,
                f"{multiple}: connection controls have mismatched summary vectors")
        keys = ("FPR", "FWIR", "FWPR", "FOPR", "WWIR:INJ", "WWPR:PROD",
                "WOPR:PROD", "WBHP:INJ", "WBHP:PROD", "BPR:1", "BPR:30")
        summary = max(float(np.max(np.abs(first_values[:, first_labels.index(key)] -
                                          second_values[:, second_labels.index(key)])))
                      for key in keys if key in first_labels)
        return restart, summary

    first = states["FIRST_CONN_WATERFLOOD_CONTROL"]["SWAT"]
    both = states["MULTICONN_WATERFLOOD"]["SWAT"]
    require(float(both[0]) < float(first[0]) - 0.05
            and float(both[1]) > float(both[0]) + 0.05,
            "MULTICONN_WATERFLOOD: water placement does not shift to second injector completion")
    require(float(both[29]) < float(first[29]) - 0.05,
            "MULTICONN_WATERFLOOD: second producer completion has no local water effect")
    print("MULTICONN_WATERFLOOD: water shifts toward injector cell 2 and "
          "production changes cell 30")

    for multiple, first_only in (
        ("MULTICONN_WATERFLOOD", "FIRST_CONN_WATERFLOOD_CONTROL"),
        ("MULTICONN_ORAT", "FIRST_CONN_ORAT_CONTROL"),
    ):
        restart, summary = difference(multiple, first_only)
        if restart < 1e-6 and summary < 1e-6:
            print(f"{multiple}: NO OBSERVABLE SECOND-CONNECTION EFFECT "
                  "versus first-connection-only control")
            missing_effects.append(multiple)
        else:
            print(f"{multiple}: differs from first-connection-only control "
                  f"(restart {restart:.6g}, summary {summary:.6g}); "
                  "consistent with active additional completions")

    for label, other in (("injector", "INJ_SECOND_CONN_CONTROL"),
                         ("producer", "PROD_SECOND_CONN_CONTROL")):
        restart, summary = difference(other, "FIRST_CONN_WATERFLOOD_CONTROL")
        require(restart > 1e-6 or summary > 1e-6,
                f"Connection positive control: selecting the other {label} "
                "completion alone does not change the solution")
        print(f"Connection positive control: selecting the other {label} completion "
              f"changes restart fields by {restart:.6g} and summary vectors by "
              f"{summary:.6g}")
    return missing_effects


def main():
    arguments = [arg for arg in sys.argv[1:] if arg != "--require-observable-effects"]
    require(len(arguments) <= 1, "Usage: verify_outputs.py [output-dir] [--require-observable-effects]")
    output = Path(arguments[0]) if arguments else Path("/tmp/compositional-cases-0923-runs")
    strict = "--require-observable-effects" in sys.argv[1:]
    states = {}
    for deck in sorted(DECK_DIR.glob("*/*.DATA")):
        states[deck.stem] = check_case(output, deck)
    check_stream_compositions(states)
    missing_effects = report_region_effect(states, output)
    missing_effects += report_connection_effect(states, output)
    print(f"Smoke-checked summary and restart output for {len(states)} decks")
    if strict:
        require(not missing_effects,
                "No observable effect for development targets: " + ", ".join(missing_effects))


if __name__ == "__main__":
    main()
