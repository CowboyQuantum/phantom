"""The runtime gallery protocols' stored devices (ADR 0033): every stored device belongs to a protocol,
was compiled by the code that runs now, and equals a fresh compile on this platform within 1e-9 in every
number (against the largest of its column, as the gallery compares arrays), so storing it changes
nothing but the platform's last digits."""

import json
import math
import pathlib
import re
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from gallery.protocols.step7 import DEVICES_DIR, _compile_fresh, device_key, stale
from gallery.runner import load_protocols
from phantom.device import Device

STORED = sorted(DEVICES_DIR.glob("*.json"))


def _protocol_for(key: str):
    return next((p for p in load_protocols() if p.runner.split(":")[0] in ("step7", "step7b", "step7c", "step8") and device_key(p) == key), None)


def _numbers(x, path=""):
    """(path, value) for every number in a JSON tree."""
    if isinstance(x, dict):
        for k, v in x.items():
            yield from _numbers(v, f"{path}.{k}")
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from _numbers(v, f"{path}[{i}]")
    elif isinstance(x, (int, float)) and not isinstance(x, bool):
        yield path, float(x)


def test_there_is_a_stored_device():
    assert STORED, "run `uv run python scripts/gallery.py --record-devices`"


@pytest.mark.parametrize("path", STORED, ids=[p.stem for p in STORED])
def test_a_stored_device_is_a_fresh_compile_on_this_platform(path):
    p = _protocol_for(path.stem)
    assert p is not None, f"{path.name} belongs to no runtime protocol: remove it or re-record"
    stored = Device.from_json(path)
    assert not stale(stored), f"{path.name} is stale: {stale(stored)}"
    fresh = json.loads(_compile_fresh(p).dumps())
    old = dict(_numbers(json.loads(stored.dumps())))
    new = dict(_numbers(fresh))
    assert set(old) == set(new)
    # as the gallery compares arrays (ADR 0014): each number against the largest of its column, since
    # the shape grid stores changes, differences of nearly equal numbers whose own size is no scale
    scale: dict[str, float] = {}
    for k, a in old.items():
        col = _column(k)
        scale[col] = max(scale.get(col, 0.0), abs(a) if not math.isnan(a) else 0.0)
    worst = max((abs(a - new[k]) / max(scale[_column(k)], 1e-300), k) for k, a in old.items() if not (math.isnan(a) and math.isnan(new[k])))
    assert worst[0] <= 1e-9, worst


def _column(path: str) -> str:
    """The path with every index but the last blanked: one quantity across the nodes and gates."""
    parts = re.split(r"(\[\d+\])", path)
    indices = [i for i, x in enumerate(parts) if re.fullmatch(r"\[\d+\]", x)]
    return "".join("[]" if i in indices[:-1] else x for i, x in enumerate(parts))
