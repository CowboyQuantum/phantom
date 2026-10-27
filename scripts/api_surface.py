"""ADR 0024: snapshot the consumer surface (`phantom.api.__all__`) into docs/api/surface.json.

    uv run python scripts/api_surface.py            # rewrite the snapshot
    uv run python scripts/api_surface.py --check    # CI: fail if the surface differs from the snapshot

Each public name records its kind, and: a function's signature; a class's dataclass fields with their
annotations and whether they have a default, its constructor signature, and its public methods and
properties with their signatures; a value's type and value. A change to the surface shows up here as
a diff, and the pull request that makes it lists it in CHANGELOG.md (a removal or an incompatible
change is a breaking release).
"""

from __future__ import annotations

import dataclasses
import inspect
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/api/surface.json"


def _signature(obj) -> str:
    try:
        return str(inspect.signature(obj))
    except (TypeError, ValueError):
        return "(?)"


def _class(cls: type) -> dict:
    entry: dict = {"kind": "dataclass" if dataclasses.is_dataclass(cls) else "class"}
    if dataclasses.is_dataclass(cls):
        entry["fields"] = [
            {
                "name": f.name,
                "type": f.type if isinstance(f.type, str) else getattr(f.type, "__name__", str(f.type)),
                "default": f.default is not dataclasses.MISSING or f.default_factory is not dataclasses.MISSING,
            }
            for f in dataclasses.fields(cls)
            if f.init
        ]
        entry["frozen"] = cls.__dataclass_params__.frozen
    else:
        entry["init"] = _signature(cls.__init__)
    members: dict = {}
    for name, value in sorted(vars(cls).items()):
        if name.startswith("_"):
            continue
        if isinstance(value, property):
            members[name] = {"kind": "property"}
        elif isinstance(value, classmethod | staticmethod):
            members[name] = {"kind": type(value).__name__, "signature": _signature(value.__func__)}
        elif inspect.isfunction(value):
            members[name] = {"kind": "method", "signature": _signature(value)}
    entry["members"] = members
    return entry


def surface() -> dict:
    from phantom import api

    names: dict = {}
    for name in sorted(api.__all__):
        obj = getattr(api, name)
        if inspect.isclass(obj):
            names[name] = _class(obj)
        elif inspect.isfunction(obj):
            names[name] = {"kind": "function", "signature": _signature(obj)}
        else:
            names[name] = {"kind": "value", "type": type(obj).__name__, "value": obj}
    return {"module": "phantom.api", "names": names}


def render() -> str:
    return json.dumps(surface(), indent=1, sort_keys=True, ensure_ascii=False) + "\n"


def main() -> int:
    text = render()
    if "--check" in sys.argv[1:]:
        current = OUT.read_text(encoding="utf-8") if OUT.exists() else ""
        if current != text:
            print("api_surface: FAIL — phantom.api changed; run `uv run python scripts/api_surface.py`, review the diff and list it in CHANGELOG.md")
            return 1
        print("api_surface: OK")
        return 0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(text, encoding="utf-8")
    print(f"api_surface: wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
