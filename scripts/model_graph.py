"""ADR 0011: regenerate the module dependency graph in docs/MODELS.md from the imports in src/phantom.

    uv run python scripts/model_graph.py            # rewrite the generated section
    uv run python scripts/model_graph.py --check    # CI: fail if the section is stale

The generated section is the import graph between phantom modules (one node per module that defines a
model, one edge per `from phantom.x.y import ...`). The hand-drawn measurement-flow diagram above it in
the same file is maintained by whoever adds or rewires a model.
"""

from __future__ import annotations

import ast
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = ROOT / "src/phantom"
DOC = ROOT / "docs/MODELS.md"
START, END = "<!-- module-graph:start -->", "<!-- module-graph:end -->"
LAYERS = ["layout", "electrostatics", "disorder", "valley", "charge", "readout", "spin", "optimum", "compile"]


def module_name(path: pathlib.Path) -> str:
    rel = path.relative_to(SRC).with_suffix("")
    parts = [p for p in rel.parts if p != "__init__"]
    return ".".join(parts)


def scan() -> tuple[dict[str, str], set[tuple[str, str]]]:
    nodes: dict[str, str] = {}
    edges: set[tuple[str, str]] = set()
    for path in sorted(SRC.rglob("*.py")):
        if path.name == "__init__.py":
            continue
        name = module_name(path)
        tree = ast.parse(path.read_text(encoding="utf-8"))
        doc = ast.get_docstring(tree) or ""
        nodes[name] = doc.strip().splitlines()[0] if doc else ""
        for node in ast.walk(tree):
            if isinstance(node, ast.ImportFrom) and node.module and node.module.startswith("phantom."):
                target = node.module.removeprefix("phantom.")
                if target.count(".") == 0:
                    target_layer = target
                    for alias in node.names:
                        edges.add((name, f"{target_layer}.{alias.name}") if False else (name, target))
                else:
                    edges.add((name, target))
    return nodes, edges


def render(nodes: dict[str, str], edges: set[tuple[str, str]]) -> str:
    lines = ["```mermaid", "flowchart LR"]
    by_layer: dict[str, list[str]] = {}
    for n in nodes:
        by_layer.setdefault(n.split(".")[0], []).append(n)
    for layer in LAYERS:
        if layer not in by_layer:
            continue
        lines.append(f"  subgraph {layer}")
        for n in sorted(by_layer[layer]):
            short = n.split(".", 1)[1] if "." in n else n
            lines.append(f'    {n.replace(".", "_")}["{short}"]')
        lines.append("  end")
    for src, dst in sorted(edges):
        if dst in nodes:
            lines.append(f"  {src.replace('.', '_')} --> {dst.replace('.', '_')}")
    lines.append("```")
    lines.append("")
    lines.append("| module | what it owns |")
    lines.append("|---|---|")
    for n in sorted(nodes):
        lines.append(f"| `{n}` | {nodes[n]} |")
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    nodes, edges = scan()
    text = DOC.read_text(encoding="utf-8")
    if START not in text or END not in text:
        print(f"model_graph: {DOC.relative_to(ROOT)} lacks the {START} / {END} markers")
        return 1
    head, rest = text.split(START, 1)
    _, tail = rest.split(END, 1)
    new = f"{head}{START}\n{render(nodes, edges)}\n{END}{tail}"
    if new == text:
        print(f"model_graph: OK ({len(nodes)} modules, {len(edges)} edges)")
        return 0
    if "--check" in argv:
        print("model_graph: docs/MODELS.md module graph is stale; run scripts/model_graph.py")
        return 1
    DOC.write_text(new, encoding="utf-8", newline="\n")
    print(f"model_graph: rewrote ({len(nodes)} modules, {len(edges)} edges)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
