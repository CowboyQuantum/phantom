"""ADR 0003/0004/0024: no qdt or instrument SDK imports anywhere; layers import only downward; the
consumer surface (`phantom.api`, `phantom.contract`) sits above every layer and nothing below imports it."""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src/phantom"
ORDER = ["layout", "electrostatics", "disorder", "valley", "charge", "readout", "spin", "device", "optimum", "compile"]
SURFACE = ["api", "contract"]
"""Modules above `compile`: they may import every layer; no layer and no root module imports them,
except that `contract` imports `api`."""
FORBIDDEN = re.compile(r"^\s*(from|import)\s+(qdt_\w+|qdt\b|qcodes|pyvisa|qm\b|qm\.|qdevil)")
IMPORT = re.compile(r"^\s*(?:from|import)\s+phantom\.(\w+)")
FROM_ROOT = re.compile(r"^\s*from\s+phantom\s+import\s+\(?\s*([\w\s,]+)")
"""`from phantom import api, compile` names modules too."""


def owner(p: pathlib.Path) -> str:
    """The layer a file belongs to, or its module name for a file at the package root."""
    rel = p.relative_to(SRC)
    return rel.parts[0] if len(rel.parts) > 1 else rel.stem


def violations(src: str, dst: str, where: str) -> list[str]:
    out = []
    if dst in SURFACE and src != dst and not (src == "contract" and dst == "api"):
        out.append(f"{where}: {src} imports phantom.{dst} (the consumer surface; nothing below imports it)")
    if src in ORDER and dst in ORDER and src != dst:
        if src != "compile" and ORDER.index(dst) >= ORDER.index(src):
            out.append(f"{where}: {src} imports {dst} (not downward)")
        if dst == "compile":
            out.append(f"{where}: {src} imports compile (only the consumer surface does)")
    if src not in ORDER and src not in SURFACE and dst in ORDER:
        out.append(f"{where}: root module {src} imports layer {dst}")
    return out


bad = []
for p in sorted(SRC.rglob("*.py")):
    src = owner(p)
    for i, line in enumerate(p.read_text(errors="ignore").splitlines(), 1):
        where = f"{p.relative_to(ROOT)}:{i}"
        if FORBIDDEN.search(line):
            bad.append(f"{where}: forbidden import: {line.strip()}")
        m = IMPORT.match(line)
        r = FROM_ROOT.match(line)
        targets = [m.group(1)] if m else [n.strip().split()[0] for n in r.group(1).split(",") if n.strip()] if r else []
        for dst in targets:
            bad.extend(violations(src, dst, where))
if bad:
    print("check_imports: FAIL\n  " + "\n  ".join(bad))
    sys.exit(1)
print("check_imports: OK")
