"""Render the gallery (ADR 0008) and check sign-offs and budgets (ADR 0009).

    uv run python scripts/gallery.py                      # all protocols -> build/gallery/
    uv run python scripts/gallery.py --only noise-psd     # one protocol
    uv run python scripts/gallery.py --strict             # wall-time factor 1 (local benchmarking)
    uv run python scripts/gallery.py --no-render          # arrays, hashes and ledger only (no matplotlib)

Exit status is non-zero when a signed protocol's array hash has changed, or a budget is exceeded.
"""

from __future__ import annotations

import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from gallery.runner import CI_TIME_FACTOR, ROOT, load_protocols, run_protocol, write_outputs


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=pathlib.Path, default=ROOT / "build/gallery")
    ap.add_argument("--only", nargs="*", default=None, help="protocol ids")
    ap.add_argument("--strict", action="store_true", help=f"wall-time factor 1 instead of {CI_TIME_FACTOR:g}")
    ap.add_argument("--no-render", action="store_true")
    ap.add_argument("--allow-mismatch", action="store_true", help="report but do not fail on a sign-off mismatch")
    args = ap.parse_args()
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")

    protocols = load_protocols(args.only)
    if not protocols:
        print("gallery: no protocols found")
        return 1
    factor = 1.0 if args.strict else CI_TIME_FACTOR
    results = []
    failed = False
    print(f"{'protocol':34} {'points':>7} {'compile s':>10} {'run ms/pt':>10} {'evals/pt':>9} {'budget':>7} {'sign-off':>9}")
    for p in protocols:
        r = run_protocol(p, factor)
        results.append(r)
        per_pt = 1e3 * r.run_s / max(r.points, 1)
        evals = r.counts["model_evals"] / max(r.points, 1)
        print(f"{p.id:34} {r.points:>7} {r.compile_s:>10.3f} {per_pt:>10.4f} {evals:>9g} {'ok' if r.budget_ok else 'OVER':>7} {r.signoff_status:>9}")
        for note in r.budget_notes:
            print(f"    budget: {note}")
        if not r.budget_ok:
            failed = True
        if r.signoff_status == "MISMATCH" and not args.allow_mismatch:
            failed = True
            print(f"    sign-off: array hash {r.sha256[:16]}… differs from signed {str(r.signoff.get('array_sha256'))[:16]}…")
    write_outputs(results, args.out, render=not args.no_render)
    print(f"gallery: wrote {args.out} ({len(results)} protocols){'; FAIL' if failed else ''}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
