"""ADR 0003/0004: no qdt or instrument SDK imports anywhere; layers import only downward."""
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
ORDER = ["layout", "electrostatics", "disorder", "valley", "charge", "readout", "spin", "optimum", "compile"]
FORBIDDEN = re.compile(r"^\s*(from|import)\s+(qdt_\w+|qdt\b|qcodes|pyvisa|qm\b|qm\.|qdevil)")
bad = []
for layer in ORDER:
    for p in (ROOT / "src/phantom" / layer).rglob("*.py"):
        for i, line in enumerate(p.read_text(errors="ignore").splitlines(), 1):
            if FORBIDDEN.search(line):
                bad.append(f"{p.relative_to(ROOT)}:{i}: forbidden import: {line.strip()}")
            m = re.match(r"^\s*(?:from|import)\s+phantom\.(\w+)", line)
            if m and m.group(1) in ORDER:
                src, dst = ORDER.index(layer), ORDER.index(m.group(1))
                if layer != "compile" and dst >= src and m.group(1) != layer:
                    bad.append(f"{p.relative_to(ROOT)}:{i}: {layer} imports {m.group(1)} (not downward)")
                if m.group(1) == "compile" and layer != "compile":
                    bad.append(f"{p.relative_to(ROOT)}:{i}: nothing imports compile")
if bad:
    print("check_imports: FAIL\n  " + "\n  ".join(bad)); sys.exit(1)
print("check_imports: OK")
