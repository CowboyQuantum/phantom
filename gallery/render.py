"""PNG per protocol and a static index page. Cosmetic only: never compared, never committed."""

from __future__ import annotations

import html
import pathlib
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from gallery.runner import Result


def render_protocol(r: Result, path: pathlib.Path) -> None:
    spec = r.render_spec if hasattr(r, "render_spec") else r.protocol.render
    panels = spec.get("panels") or [spec]
    fig, axes = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 4.2), squeeze=False)
    for ax, panel in zip(axes[0], panels):
        kind = panel.get("kind", "line")
        x = r.arrays[panel["x"]]
        if kind == "line":
            ys = panel["y"] if isinstance(panel["y"], list) else [panel["y"]]
            style = {"marker": "o", "ms": 4, "ls": "none"} if panel.get("markers") else {"lw": 1.0}
            for key in ys:
                ax.plot(x, r.arrays[key], label=key if len(ys) > 1 else None, **style)
            if len(ys) > 1:
                ax.legend(fontsize=8)
            if panel.get("log_x"):
                ax.set_xscale("log")
            if panel.get("log_y"):
                ax.set_yscale("log")
        elif kind == "map":
            y = r.arrays[panel["y"]]
            z = r.arrays[panel["z"]]
            limits = {}
            if panel.get("robust"):
                # clip the colour scale to the 1st–99th percentile so a few spikes do not hide the rest;
                # symmetric about zero for a diverging map
                lo, hi = np.percentile(z, [1, 99])
                if "RdBu" in panel.get("cmap", "") or "coolwarm" in panel.get("cmap", ""):
                    m_abs = max(abs(lo), abs(hi))
                    lo, hi = -m_abs, m_abs
                limits = {"vmin": float(lo), "vmax": float(hi)}
            m = ax.pcolormesh(x, y, z.T, shading="auto", cmap=panel.get("cmap", "viridis"), rasterized=True, **limits)
            fig.colorbar(m, ax=ax, label=panel.get("zlabel", panel["z"]))
            if panel.get("overlay"):
                poly = r.arrays[panel["overlay"]]  # NaN-separated polyline (m, 2)
                ax.plot(poly[:, 0], poly[:, 1], color="white", lw=0.8, alpha=0.9)
                ax.set_aspect("equal")
        else:
            raise ValueError(f"unknown render kind {kind!r}")
        ax.set_xlabel(panel.get("xlabel", panel["x"]))
        ax.set_ylabel(panel.get("ylabel", panel.get("y", "")))
        if panel.get("title"):
            ax.set_title(panel["title"], fontsize=10)
    fig.suptitle(f"{r.protocol.id} — {r.protocol.name}  [{r.signoff_status}]", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def render_index(results: list[Result], rows: list[dict[str, Any]], meta: dict[str, Any], path: pathlib.Path) -> None:
    def esc(x: Any) -> str:
        return html.escape(str(x))

    style = (
        "<style>body{font-family:system-ui,sans-serif;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#222}"
        "table{border-collapse:collapse;font-size:0.85rem}td,th{border:1px solid #ccc;padding:0.25rem 0.5rem;text-align:left}"
        ".ok{color:#1a7f37}.bad{color:#b60205;font-weight:600}.unsigned{color:#9a6700}img{max-width:100%;border:1px solid #ddd}"
        "section{margin:2rem 0}code{font-size:0.85em}</style></head><body>"
    )
    header = (
        f"<h1>phantom gallery</h1><p>commit <code>{esc(meta['commit'])}</code> · {esc(meta['platform'])} · "
        f"python {esc(meta['python'])} · numpy {esc(meta['numpy'])}. Rendered by CI; nothing here is stored in git (ADR 0008).</p>"
    )
    ledger_head = (
        "<h2>Budget ledger (ADR 0009)</h2><table><tr><th>protocol</th><th>points</th><th>compile s / budget</th>"
        "<th>run ms/pt / budget</th><th>evals/pt / declared</th><th>budget</th><th>sign-off</th></tr>"
    )
    parts = [
        "<!doctype html><html><head><meta charset='utf-8'><title>phantom gallery</title>",
        style,
        header,
        ledger_head,
    ]
    for row in rows:
        b_cls = "ok" if row["budget_ok"] else "bad"
        s_cls = {"signed": "ok", "unsigned": "unsigned"}.get(row["signoff"], "bad")
        parts.append(
            f"<tr><td><a href='#{esc(row['id'])}'><code>{esc(row['id'])}</code></a></td><td>{row['points']}</td>"
            f"<td>{row['compile_s']} / {row['compile_budget_s']}</td>"
            f"<td>{row['run_ms_per_point']} / {row['run_budget_ms_per_point']}</td>"
            f"<td>{row['model_evals_per_point']:g} / {row['declared_model_evals_per_point']}</td>"
            f"<td class='{b_cls}'>{'ok' if row['budget_ok'] else esc('; '.join(row['budget_notes']))}</td>"
            f"<td class='{s_cls}'>{esc(row['signoff'])}</td></tr>"
        )
    parts.append("</table>")
    for r in results:
        p = r.protocol
        must = "".join(f"<li>{esc(m)}</li>" for m in p.must_show) or "<li><em>none recorded yet</em></li>"
        parts.append(
            f"<section id='{esc(p.id)}'><h2><code>{esc(p.id)}</code> — {esc(p.name)}</h2>"
            f"<p>{esc(p.description)}</p><p>layers: {esc(', '.join(p.layers))} · step {esc(p.step)} · "
            f"sign-off: <span class='{ {'signed': 'ok', 'unsigned': 'unsigned'}.get(r.signoff_status, 'bad') }'>{esc(r.signoff_status)}</span> · "
            f"sha256 <code>{esc(r.sha256[:16])}…</code></p>"
            f"<img src='{esc(p.id)}.png' alt='{esc(p.name)}'><p>Must show:</p><ul>{must}</ul></section>"
        )
    parts.append("</body></html>")
    path.write_text("\n".join(parts), encoding="utf-8")
