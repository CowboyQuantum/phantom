"""ADR 0007: a `literature` test runs only against a `confirmed` reference entry.

`@pytest.mark.literature("<id>", item="<label>")` names an entry in `docs/references/entries/`. A
missing entry or item fails collection; a `proposed` entry skips the test with a reason, so a P1 row
never passes on an unconfirmed standard.
"""

from __future__ import annotations

import pathlib

import pytest
import yaml

ENTRIES = pathlib.Path(__file__).resolve().parents[1] / "docs/references/entries"


def _entry(eid: str) -> dict | None:
    path = ENTRIES / f"{eid}.yaml"
    return yaml.safe_load(path.read_text(encoding="utf-8")) if path.is_file() else None


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    for item in items:
        mark = item.get_closest_marker("literature")
        if mark is None:
            continue
        if not mark.args:
            raise pytest.UsageError(f"{item.nodeid}: literature marker needs a reference id")
        eid, label = mark.args[0], mark.kwargs.get("item")
        entry = _entry(eid)
        if entry is None:
            raise pytest.UsageError(f"{item.nodeid}: no reference entry {eid!r} (ADR 0007)")
        labels = {i.get("label") for i in entry.get("provides") or []}
        if label not in labels:
            raise pytest.UsageError(f"{item.nodeid}: entry {eid!r} has no item {label!r}")
        if entry.get("status") != "confirmed":
            item.add_marker(
                pytest.mark.skip(reason=f"reference {eid} is {entry.get('status')}, not confirmed (ADR 0007)")
            )
