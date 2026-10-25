"""ADR 0007: the reference library is the only standard; every citation names an entry that exists.

Checks `docs/references/entries/*.yaml` against `docs/references/SCHEMA.md`, every
`@pytest.mark.literature("<id>", item="<label>")` in `tests/`, and the `sources` and `rate.source` fields
of `docs/catalog/entries/*.yaml`; then regenerates the index table in `docs/references/README.md`.

    uv run python scripts/check_references.py            # check and rewrite the index
    uv run python scripts/check_references.py --check    # CI: fail if the index is stale
"""

from __future__ import annotations

import pathlib
import re
import sys

import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
ENTRIES = ROOT / "docs/references/entries"
DATA = ROOT / "docs/references/data"
README = ROOT / "docs/references/README.md"
CATALOG = ROOT / "docs/catalog/entries"
TESTS = ROOT / "tests"

ID_RE = re.compile(r"^[a-z][a-z0-9]*[0-9]{4}(-[a-z0-9]+)+$")
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")
LAYERS = {"layout", "electrostatics", "disorder", "valley", "charge", "readout", "spin"}
KINDS = {"figure", "table", "value", "range", "formula"}
MARK_RE = re.compile(r"pytest\.mark\.literature\(\s*[\"']([^\"']+)[\"']\s*(?:,\s*item\s*=\s*[\"']([^\"']+)[\"'])?")
START, END = "<!-- index:start -->", "<!-- index:end -->"


def load_entries(problems: list[str]) -> dict[str, dict]:
    entries: dict[str, dict] = {}
    for path in sorted(ENTRIES.glob("*.yaml")):
        try:
            e = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        except yaml.YAMLError as exc:
            problems.append(f"{path.name}: not valid YAML ({exc})")
            continue
        eid = e.get("id")
        if eid != path.stem:
            problems.append(f"{path.name}: id {eid!r} must equal the file name")
        if not isinstance(eid, str) or not ID_RE.match(eid):
            problems.append(f"{path.name}: id must look like <firstauthor><year>-<short>")
        for key in ("title", "authors", "year", "layers", "provides", "status"):
            if key not in e:
                problems.append(f"{path.name}: missing field {key!r}")
        doi = e.get("doi")
        if not (isinstance(doi, str) and DOI_RE.match(doi)) and not e.get("arxiv"):
            problems.append(f"{path.name}: doi {doi!r} is not a DOI and no arxiv given")
        if e.get("status") not in {"proposed", "confirmed"}:
            problems.append(f"{path.name}: status must be proposed or confirmed")
        if e.get("status") == "confirmed" and not (e.get("confirmed_by") and e.get("confirmed_on")):
            problems.append(f"{path.name}: confirmed entries need confirmed_by and confirmed_on")
        if not set(e.get("layers") or []) <= LAYERS:
            problems.append(f"{path.name}: layers must be a subset of {sorted(LAYERS)}")
        labels: set[str] = set()
        for item in e.get("provides") or []:
            label, kind = item.get("label"), item.get("kind")
            if not label or label in labels:
                problems.append(f"{path.name}: provides items need unique labels (got {label!r})")
            labels.add(label or "")
            if kind not in KINDS:
                problems.append(f"{path.name}: item {label!r} kind {kind!r} not in {sorted(KINDS)}")
            if "quantity" not in item:
                problems.append(f"{path.name}: item {label!r} needs a quantity")
            if kind in {"figure", "table"}:
                data = item.get("data")
                if not data or not (DATA / data).is_file():
                    problems.append(f"{path.name}: item {label!r} data {data!r} not found under data/")
        if isinstance(eid, str):
            entries[eid] = e
    return entries


def check_tests(entries: dict[str, dict], problems: list[str]) -> None:
    for path in sorted(TESTS.rglob("*.py")):
        if path.name == "conftest.py":
            continue
        text = path.read_text(encoding="utf-8")
        for m in MARK_RE.finditer(text):
            eid, label = m.group(1), m.group(2)
            if eid not in entries:
                problems.append(f"{path.relative_to(ROOT)}: literature id {eid!r} has no entry")
            elif label and label not in {i.get("label") for i in entries[eid].get("provides") or []}:
                problems.append(f"{path.relative_to(ROOT)}: entry {eid!r} has no item {label!r}")
            elif not label:
                problems.append(f"{path.relative_to(ROOT)}: literature marker for {eid!r} needs item=")
        if "literature" in text and "pytest.mark.literature" not in text and "conftest" not in path.name:
            pass  # a word, not a marker


def check_catalog(entries: dict[str, dict], problems: list[str], warnings: list[str]) -> None:
    for path in sorted(CATALOG.glob("*.yaml")):
        e = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        sink = warnings if e.get("status") == "draft" else problems
        cited = list(e.get("sources") or [])
        rate = e.get("rate")
        if isinstance(rate, dict) and rate.get("source") not in (None, "unknown"):
            cited.append(rate["source"])
        for s in cited:
            if s not in entries:
                sink.append(f"catalog {path.name}: source {s!r} is not a reference id")
        if "calibration" not in e:
            sink.append(f"catalog {path.name}: missing calibration field (reference | none)")


def render_index(entries: dict[str, dict]) -> str:
    rows = ["| id | year | layers | status | items |", "|---|---|---|---|---|"]
    for eid, e in sorted(entries.items()):
        items = ", ".join(f"`{i.get('label')}`" for i in e.get("provides") or [])
        rows.append(
            f"| [`{eid}`](entries/{eid}.yaml) | {e.get('year')} | {', '.join(e.get('layers') or [])} "
            f"| {e.get('status')} | {items} |"
        )
    return "\n".join(rows)


def update_readme(entries: dict[str, dict], check_only: bool) -> bool:
    text = README.read_text(encoding="utf-8")
    if START not in text or END not in text:
        print(f"check_references: {README.relative_to(ROOT)} lacks {START} / {END} markers")
        return False
    head, rest = text.split(START, 1)
    _, tail = rest.split(END, 1)
    new = f"{head}{START}\n{render_index(entries)}\n{END}{tail}"
    if new == text:
        return True
    if check_only:
        print("check_references: index table is stale; run scripts/check_references.py")
        return False
    README.write_text(new, encoding="utf-8", newline="\n")
    return True


def main(argv: list[str]) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    check_only = "--check" in argv
    problems: list[str] = []
    warnings: list[str] = []
    entries = load_entries(problems)
    check_tests(entries, problems)
    check_catalog(entries, problems, warnings)
    ok = update_readme(entries, check_only)
    for w in warnings:
        print(f"  warning (draft): {w}")
    if problems or not ok:
        print("check_references: FAIL\n  " + "\n  ".join(problems))
        return 1
    confirmed = sum(e.get("status") == "confirmed" for e in entries.values())
    print(f"check_references: OK ({len(entries)} entries, {confirmed} confirmed)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
