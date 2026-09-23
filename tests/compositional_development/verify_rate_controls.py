#!/usr/bin/env python3
"""Check well-rate controls in the compositional development decks and runs."""

from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "tools"))
import eclio


DECK_DIR = Path(__file__).resolve().parent


def final_rates(output, name):
    base = output / name / name
    labels, values = eclio.summary(str(base) + ".SMSPEC", str(base) + ".UNSMRY")
    return {label: float(values[-1, i]) for i, label in enumerate(labels)}


def check(name, actual, target):
    tolerance = max(1e-4, 0.005 * target)
    if abs(actual - target) > tolerance:
        raise SystemExit(f"{name}: expected {target}, got {actual:.6g}")
    print(f"{name}: {actual:.6g} (target {target})")


def injector_controls():
    controls = []
    for deck in sorted(DECK_DIR.rglob("*.DATA")):
        content = deck.read_text()
        blocks = re.findall(r"^WCONINJE[ \t]*\n(.*?)^[ \t]*/[ \t]*$", content,
                            flags=re.MULTILINE | re.DOTALL)
        if "WCONINJE" in content and not blocks:
            raise SystemExit(f"{deck}: cannot parse WCONINJE")
        for block in blocks:
            for line in block.splitlines():
                record = line.split("--", 1)[0].strip()
                if not record:
                    continue
                fields = record.removesuffix("/").split()
                if len(fields) < 7 or fields[3] != "RATE":
                    raise SystemExit(f"{deck}: injector must declare RATE: {record}")
                if fields[1] not in ("GAS", "WATER"):
                    raise SystemExit(f"{deck}: unexpected injector phase: {record}")
                controls.append((deck, fields[0], fields[1], float(fields[4]),
                                 float(fields[6])))
    return controls


def check_injectors(output):
    controls = injector_controls()
    for deck, well, phase, target, bhp_limit in controls:
        rates = final_rates(output, deck.stem)
        keyword = ("WGIR" if phase == "GAS" else "WWIR") + ":" + well
        actual = rates[keyword]
        if actual <= 0:
            raise SystemExit(f"{deck.stem} {keyword}: expected positive injection, got {actual}")
        tolerance = max(1e-4, 0.005 * target)
        if abs(actual - target) <= tolerance:
            print(f"{deck.stem} {keyword}: {actual:.6g} (RATE target {target})")
        elif abs(rates["WBHP:" + well] - bhp_limit) <= 0.05:
            print(f"{deck.stem} {keyword}: {actual:.6g} (BHP limit {bhp_limit} active)")
        else:
            raise SystemExit(f"{deck.stem} {keyword}: rate {actual:.6g} missed target "
                             f"{target} without reaching BHP limit {bhp_limit}")
    print(f"Checked {len(controls)} injector RATE declarations")


def main():
    output = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/tmp/compositional-cases-0923-runs")
    checks = {
        "PROD_ORAT": ("WOPR:PROD", 0.08),
        "PROD_WRAT": ("WWPR:PROD", 0.10),
        "PROD_GRAT": ("WGPR:PROD", 100.0),
        "INJ_GAS_RATE": ("WGIR:INJ", 700.0),
        "INJ_WATER_RATE": ("WWIR:INJ", 0.20),
        "MULTICONN_ORAT": ("WOPR:PROD", 0.04),
    }
    for name, (keyword, target) in checks.items():
        check(name + " " + keyword, final_rates(output, name)[keyword], target)
    for name in ("PROD_LRAT", "PROD_LRAT_LIMIT"):
        liquid = final_rates(output, name)
        check(name + " WOPR+WWPR",
              liquid["WOPR:PROD"] + liquid["WWPR:PROD"], 0.20)
    check("PROD_RESV WVPR:PROD", final_rates(output, "PROD_RESV")["WVPR:PROD"], 0.20)
    check_injectors(output)


if __name__ == "__main__":
    main()
