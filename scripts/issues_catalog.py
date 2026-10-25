"""File the catalog entries a human should write as GitHub issues (label `catalog`). Dry run by default.

    uv run python scripts/issues_catalog.py            # print
    uv run python scripts/issues_catalog.py --apply    # needs `gh auth login`
"""
import json
import subprocess
import sys

ENTRIES = [
    ("oxide-traps-fluctuators", "Oxide traps and two-level fluctuators", "readout", "Charge noise spectrum (amplitude at 1 Hz, 1/f slope), telegraph switchers; distribution of switching rates; MOS vs SiGe."),
    ("gate-misalignment", "Gate misalignment and width error", "electrostatics", "Lithographic offset and CD error as transforms of the lever-arm basis; the dot is not under the gate the layout says."),
    ("gate-overlap-leakage", "Gate-to-gate overlap and dielectric leakage", "layout", "Leaky barriers and current paths between gates; what a leak looks like in an IV or a stability diagram."),
    ("gate-open-short", "Open or shorted gate", "layout", "A gate that does nothing, or two gates that move together; trivial physics, common in practice."),
    ("ge-fluctuation-valley", "Ge-concentration fluctuation and valley splitting", "valley", "Distribution of valley splitting vs Ge profile and dot size; the fraction of dots that cannot do PSB."),
    ("interface-steps", "Interface steps and roughness", "valley", "Atomic steps under the dot: valley-phase changes with dot position; when to model explicitly."),
    ("strain-variation", "Strain variation", "spin", "g-factor and level-spacing shifts across a wafer and under gates."),
    ("residual-dopants", "Residual dopants and unintentional doping", "disorder", "Isolated charged centres: unintentional dots, hopping sites, leakage."),
    ("dielectric-thickness", "Dielectric thickness and 2DEG depth variation", "electrostatics", "Lever-arm scale variation; closed-form validity vs depth."),
    ("twodeg-density-mobility", "2DEG density and mobility by material class", "electrostatics", "Screening and disorder scale; what MOS and SiGe defaults should be."),
    ("sensor-coupling-asymmetry", "Sensor-dot coupling asymmetry", "readout", "Why one dot is visible and the other barely; capacitive coupling distributions."),
    ("slow-charge-drift", "Slow charge drift", "readout", "Drift on minutes-to-hours timescales; what causes it; how it looks in repeated scans."),
    ("material-class-summary", "Material classes: SiGe vs MOS as summary traits", "layout", "The defaults each class sets across layers, in one place, with sources."),
    ("ge-holes-later", "Germanium holes (later)", "spin", "Whether and when hole devices join; what changes (spin–orbit, g-tensor)."),
]


def main():
    issues = [{"title": f"Catalog: {name}", "body": f"Write `docs/catalog/entries/{eid}.yaml` per `docs/catalog/SCHEMA.md`.\n\nLayer: `{layer}`.\n\n{hint}\n\nStatus `reviewed` is set by a human; cite sources; no real device values.", "labels": ["catalog", f"layer:{layer}"]} for eid, name, layer, hint in ENTRIES]
    if "--apply" not in sys.argv:
        print(json.dumps([i["title"] for i in issues], indent=1)); print(f"\n{len(issues)} issues (dry run)", file=sys.stderr); return
    subprocess.run(["gh", "label", "create", "catalog", "--color", "0e8a16", "--force"], capture_output=True, check=False)
    for i in issues:
        for lbl in i["labels"]:
            subprocess.run(["gh", "label", "create", lbl, "--color", "ededed", "--force"], capture_output=True, check=False)
        subprocess.run(["gh", "issue", "create", "--title", i["title"], "--body", i["body"], *[a for x in i["labels"] for a in ("--label", x)]], check=True)
        print("created:", i["title"])


if __name__ == "__main__":
    main()
