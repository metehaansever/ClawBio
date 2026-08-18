#!/usr/bin/env python3
"""
Struct Predictor Foldseek — Structural homology search with Foldseek.

Queries a protein structure (CIF or PDB) against major structural databases
and ranks hits by TM-score.

Usage:
    # Search against PDB (default)
    python skills/struct-predictor-foldseek/struct_predictor_foldseek.py \
        --input structure.cif --output /tmp/foldseek_out

    # Search multiple databases
    python skills/struct-predictor-foldseek/struct_predictor_foldseek.py \
        --input structure.cif --databases pdb,afdb --output /tmp/foldseek_out

    # Demo (Trp-cage miniprotein, no input needed)
    python skills/struct-predictor-foldseek/struct_predictor_foldseek.py \
        --demo --output /tmp/foldseek_demo

    # Chain after struct-predictor
    python skills/struct-predictor-foldseek/struct_predictor_foldseek.py \
        --input /tmp/boltz_out/predictions/Trpcage/Trpcage_model_0.cif \
        --output /tmp/foldseek_out
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)

DEMO_STRUCTURE = _SKILL_DIR / "demo_data" / "trpcage.cif"
DEMO_NAME = "Trpcage"

# Default local database cache directory
DEFAULT_DB_CACHE = Path.home() / ".foldseek_dbs"

# Foldseek column order for easy-search TSV output
# Reference: https://github.com/steineggerlab/foldseek/blob/master/README.md
_FOLDSEEK_COLS = [
    "query", "target", "pident", "alnlen", "mismatch",
    "gapopen", "qstart", "qend", "tstart", "tend",
    "evalue", "bits", "alntmscore", "lddt",
]

# Internal alias → canonical column name for TM-score
_TMSCORE_COL = "alntmscore"

# Supported database aliases → foldseek database name
_DB_ALIASES: dict[str, str] = {
    "pdb":    "pdb",
    "afdb":   "afdb50",
    "esm":    "esmatlas",
}

# Default minimum TM-score for the filtered table and report
DEFAULT_MIN_TMSCORE = 0.3

# Max hits shown in report and result.json
DEFAULT_MAX_HITS = 50


# ---------------------------------------------------------------------------
# Foldseek availability check
# ---------------------------------------------------------------------------


def _check_foldseek() -> str:
    """Return the foldseek executable path or raise RuntimeError."""
    exe = shutil.which("foldseek")
    if exe is None:
        raise RuntimeError(
            "foldseek not found on PATH.\n"
            "Install via conda:  conda install -c bioconda -c conda-forge foldseek\n"
            "Then ensure the conda env is active or the binary is on PATH."
        )
    return exe


def _foldseek_version(exe: str) -> str:
    """Return foldseek version string (best-effort)."""
    try:
        out = subprocess.run(
            [exe, "version"],
            capture_output=True, text=True, timeout=10,
        )
        return out.stdout.strip() or out.stderr.strip() or "unknown"
    except Exception:
        return "unknown"


# ---------------------------------------------------------------------------
# Database helpers
# ---------------------------------------------------------------------------


def _resolve_databases(db_arg: str) -> list[tuple[str, str]]:
    """Resolve comma-separated db aliases to list of (alias, foldseek_name) tuples."""
    requested = [d.strip().lower() for d in db_arg.split(",") if d.strip()]
    resolved: list[tuple[str, str]] = []
    unknown = []
    for alias in requested:
        if alias in _DB_ALIASES:
            resolved.append((alias, _DB_ALIASES[alias]))
        else:
            unknown.append(alias)
    if unknown:
        raise ValueError(
            f"Unknown database alias(es): {unknown}. "
            f"Supported: {list(_DB_ALIASES.keys())}"
        )
    return resolved


def _get_db_path(alias: str, foldseek_name: str, db_cache: Path) -> Path:
    """Return the local path where this database should exist."""
    return db_cache / alias / foldseek_name


# ---------------------------------------------------------------------------
# Foldseek runner
# ---------------------------------------------------------------------------


def _run_foldseek_search(
    exe: str,
    query_path: Path,
    db_path: Path,
    results_tsv: Path,
    tmp_dir: Path,
    threads: int = 1,
) -> subprocess.CompletedProcess:
    """Run foldseek easy-search and return the completed process."""
    cmd = [
        exe, "easy-search",
        str(query_path),
        str(db_path),
        str(results_tsv),
        str(tmp_dir),
        "--format-output", ",".join(_FOLDSEEK_COLS),
        "--alignment-type", "1",   # TM-align mode
        "--threads", str(threads),
        "-v", "1",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result


# ---------------------------------------------------------------------------
# TSV parsing
# ---------------------------------------------------------------------------


def _parse_hits(tsv_path: Path) -> list[dict]:
    """Parse Foldseek TSV output into a list of hit dicts."""
    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        return []

    hits: list[dict] = []
    with open(tsv_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < len(_FOLDSEEK_COLS):
                # Pad missing columns with empty string
                parts += [""] * (len(_FOLDSEEK_COLS) - len(parts))
            hit: dict = {}
            for i, col in enumerate(_FOLDSEEK_COLS):
                raw = parts[i]
                if col in ("pident", "tmscore", "rmsd", "evalue", "bits"):
                    try:
                        hit[col] = float(raw)
                    except ValueError:
                        hit[col] = 0.0
                elif col in ("alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend"):
                    try:
                        hit[col] = int(raw)
                    except ValueError:
                        hit[col] = 0
                else:
                    hit[col] = raw
            hits.append(hit)
    return hits


def _filter_and_rank_hits(
    hits: list[dict],
    min_tmscore: float,
    max_hits: int,
) -> list[dict]:
    """Filter by TM-score threshold and return top N, ranked by TM-score descending."""
    filtered = [h for h in hits if h.get(_TMSCORE_COL, 0.0) >= min_tmscore]
    filtered.sort(key=lambda h: h.get(_TMSCORE_COL, 0.0), reverse=True)
    return filtered[:max_hits]


# ---------------------------------------------------------------------------
# Report + figures
# ---------------------------------------------------------------------------


def _write_tmscore_figure(hits: list[dict], figures_dir: Path) -> Path:
    """Write a horizontal bar chart of top hits by TM-score."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        return figures_dir / "tmscore_hits.png"  # skip silently if no matplotlib

    top = hits[:20]  # cap chart to 20 bars
    if not top:
        return figures_dir / "tmscore_hits.png"

    # Color bars by TM-score band
    colors = []
    for s in scores:
        if s >= 0.9:
            colors.append("#2ecc71")
        elif s >= 0.7:
            colors.append("#27ae60")
        elif s >= 0.5:
            colors.append("#f39c12")
        elif s >= 0.3:
            colors.append("#e67e22")
        else:
            colors.append("#e74c3c")

    labels = [h["target"] for h in top]
    scores = [h[_TMSCORE_COL] for h in top]

    fig, ax = plt.subplots(figsize=(9, max(3, 0.35 * len(top) + 1.5)))
    y_pos = np.arange(len(top))
    bars = ax.barh(y_pos, scores, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("TM-score", fontsize=10)
    ax.set_title("Top Foldseek Hits by TM-score", fontsize=12, fontweight="bold")
    ax.set_xlim(0, 1.0)
    ax.axvline(x=0.5, color="gray", linestyle="--", linewidth=0.8, alpha=0.6)

    # Annotate values
    for bar, score in zip(bars, scores):
        ax.text(
            min(score + 0.01, 0.97), bar.get_y() + bar.get_height() / 2,
            f"{score:.3f}", va="center", ha="left", fontsize=7, color="#333333",
        )

    plt.tight_layout()
    out_path = figures_dir / "tmscore_hits.png"
    figures_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _write_hits_tsv(hits: list[dict], path: Path) -> None:
    """Write hits list to a TSV file with header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(_FOLDSEEK_COLS) + "\n")
        for h in hits:
            fh.write("\t".join(str(h.get(c, "")) for c in _FOLDSEEK_COLS) + "\n")


def _tmscore_band(score: float) -> str:
    if score >= 0.9:
        return "Near-identical fold"
    elif score >= 0.7:
        return "Very similar fold"
    elif score >= 0.5:
        return "Likely same fold family"
    elif score >= 0.3:
        return "Possibly related"
    else:
        return "Likely unrelated"


def _generate_report(
    output_dir: Path,
    query_path: Path,
    databases_used: list[str],
    hits_raw: list[dict],
    hits_filtered: list[dict],
    min_tmscore: float,
    max_hits: int,
    cmd: str,
    foldseek_version: str,
    demo: bool,
) -> dict:
    """Write report.md, result.json, figures, tables, and reproducibility files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"
    repro_dir = output_dir / "reproducibility"
    figures_dir.mkdir(exist_ok=True)
    tables_dir.mkdir(exist_ok=True)
    repro_dir.mkdir(exist_ok=True)

    # Write raw TSV
    _write_hits_tsv(hits_raw, tables_dir / "hits.tsv")
    # Write filtered TSV
    _write_hits_tsv(hits_filtered, tables_dir / "hits_filtered.tsv")

    # Figure
    fig_path = _write_tmscore_figure(hits_filtered, figures_dir)

    # --- Build report ---
    top_hit = hits_filtered[0] if hits_filtered else None
    query_name = query_path.name

    lines: list[str] = []
    lines.append("# Struct Predictor Foldseek Report\n")
    if demo:
        lines.append("> **Demo mode** — Trp-cage miniprotein (PDB 1L2Y), synthetic CIF.\n")
    lines.append(f"**Query:** `{query_name}`  \n")
    lines.append(f"**Databases:** {', '.join(databases_used)}  \n")
    lines.append(f"**Foldseek version:** {foldseek_version}  \n")
    lines.append(f"**TM-score threshold:** {min_tmscore}  \n")
    lines.append(f"**Hits (raw / filtered):** {len(hits_raw)} / {len(hits_filtered)}\n")
    lines.append("")

    if top_hit:
        lines.append("## Top Hit\n")
        lines.append(f"| Field | Value |")
        lines.append(f"|-------|-------|")
        lines.append(f"| Target | `{top_hit['target']}` |")
        lines.append(f"| TM-score | {top_hit[_TMSCORE_COL]:.4f} ({_tmscore_band(top_hit[_TMSCORE_COL])}) |")
        lines.append(f"| LDDT | {top_hit['lddt']:.3f} |")
        lines.append(f"| Seq. identity | {top_hit['pident'] * 100:.1f}% |")
        lines.append(f"| E-value | {top_hit['evalue']:.2e} |")
        lines.append(f"| Alignment length | {top_hit['alnlen']} residues |")
        lines.append("")
    else:
        lines.append("## Top Hit\n")
        lines.append(
            f"> No hits passed the TM-score threshold of {min_tmscore}. "
            "The query structure may be novel or unlike any structure in the searched databases."
        )
        lines.append("")

    # Hits table
    if hits_filtered:
        lines.append("## Top Structural Hits\n")
        lines.append("| Rank | Target | TM-score | LDDT | Seq. ID | E-value | Interpretation |")
        lines.append("|------|--------|----------|------|---------|---------|----------------|")
        for i, h in enumerate(hits_filtered[:20], 1):
            lines.append(
                f"| {i} | `{h['target']}` | {h[_TMSCORE_COL]:.4f} | {h['lddt']:.3f} "
                f"| {h['pident']*100:.1f}% | {h['evalue']:.2e} | {_tmscore_band(h[_TMSCORE_COL])} |"
            )
        lines.append("")

    # Figure
    if fig_path.exists():
        lines.append("## TM-score Chart\n")
        lines.append(f"![TM-score hits](figures/tmscore_hits.png)\n")
        lines.append("")

    # Reproducibility
    lines.append("## Reproducibility\n")
    lines.append("```bash")
    lines.append(cmd)
    lines.append("```\n")

    # Disclaimer
    lines.append("---\n")
    lines.append(f"*{DISCLAIMER}*\n")

    report_text = "\n".join(lines)
    (output_dir / "report.md").write_text(report_text, encoding="utf-8")

    # Reproducibility files
    (repro_dir / "commands.sh").write_text(
        f"#!/bin/bash\n# Foldseek search command\n{cmd}\n",
        encoding="utf-8",
    )
    (repro_dir / "environment.txt").write_text(
        f"foldseek {foldseek_version}\n",
        encoding="utf-8",
    )

    # result.json
    result = {
        "query": query_name,
        "databases": databases_used,
        "n_hits_raw": len(hits_raw),
        "n_hits_filtered": len(hits_filtered),
        "min_tmscore_threshold": min_tmscore,
        "top_hit": {
            "target": top_hit["target"],
            "tmscore": top_hit[_TMSCORE_COL],
            "lddt": top_hit["lddt"],
            "seqid": top_hit["pident"],
            "evalue": top_hit["evalue"],
        } if top_hit else None,
        "hits": [
            {
                "rank": i + 1,
                "target": h["target"],
                "tmscore": h[_TMSCORE_COL],
                "lddt": h["lddt"],
                "seqid": h["pident"],
                "evalue": h["evalue"],
                "alnlen": h["alnlen"],
            }
            for i, h in enumerate(hits_filtered[:max_hits])
        ],
        "demo": demo,
        "chat_summary_lines": _build_chat_summary(
            query_name, databases_used, hits_filtered, top_hit, demo
        ),
        "suggested_actions": _build_suggested_actions(top_hit),
        "disclaimer": DISCLAIMER,
    }
    (output_dir / "result.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )

    return result


def _build_chat_summary(
    query_name: str,
    databases: list[str],
    hits: list[dict],
    top_hit: dict | None,
    demo: bool,
) -> list[str]:
    lines = [f"Foldseek search of `{query_name}` against {', '.join(databases)}:"]
    if top_hit:
        lines.append(
            f"  Top hit: **{top_hit['target']}** "
            f"(TM-score {top_hit[_TMSCORE_COL]:.3f}, "
            f"seqid {top_hit['pident']*100:.1f}%)"
        )
        lines.append(
            f"  Interpretation: {_tmscore_band(top_hit[_TMSCORE_COL])}"
        )
        lines.append(f"  {len(hits)} hits passed the TM-score threshold.")
    else:
        lines.append("  No hits passed the TM-score threshold — structure may be novel.")
    if demo:
        lines.append("  (Demo mode: Trp-cage miniprotein CIF)")
    return lines


def _build_suggested_actions(top_hit: dict | None) -> list[str]:
    if top_hit is None:
        return ["Run with --min-tmscore 0.1 to see low-confidence hits"]
    actions = []
    if top_hit[_TMSCORE_COL] >= 0.5:
        pdb_id = top_hit["target"].split("_")[0].upper()
        actions.append(
            f"Fetch PDB entry {pdb_id} for detailed comparison: "
            f"https://www.rcsb.org/structure/{pdb_id}"
        )
    if top_hit[_TMSCORE_COL] < 0.7:
        actions.append("Consider searching additional databases with --databases pdb,afdb")
    return actions


# ---------------------------------------------------------------------------
# Top-level pipeline
# ---------------------------------------------------------------------------


def run_foldseek_search(
    input_path: Path | None,
    output_dir: Path,
    databases: str = "pdb",
    min_tmscore: float = DEFAULT_MIN_TMSCORE,
    max_hits: int = DEFAULT_MAX_HITS,
    db_cache: Path | None = None,
    threads: int = 1,
    demo: bool = False,
) -> dict:
    """Run the full Foldseek structural search pipeline.

    Args:
        input_path: Path to query CIF or PDB. Required unless demo=True.
        output_dir: Where to write the report and artefacts.
        databases: Comma-separated database aliases (pdb, afdb, esm).
        min_tmscore: Minimum TM-score for filtered hit table.
        max_hits: Maximum number of hits reported.
        db_cache: Local directory holding Foldseek databases. Defaults to ~/.foldseek_dbs/.
        threads: Number of CPU threads for Foldseek.
        demo: Use the bundled Trp-cage CIF for a quick offline test.

    Returns:
        result dict (same content as result.json).

    Raises:
        ValueError: If neither input_path nor demo is supplied.
        RuntimeError: If foldseek is not on PATH.
    """
    if not demo and input_path is None:
        raise ValueError("Provide --input or --demo.")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if db_cache is None:
        db_cache = DEFAULT_DB_CACHE

    # Check foldseek availability
    exe = _check_foldseek()
    version = _foldseek_version(exe)
    print(f"  Foldseek: {exe} ({version})")

    # Resolve query path
    if demo:
        query_path = DEMO_STRUCTURE
        print(f"  Demo mode: {DEMO_NAME} (Trp-cage miniprotein, PDB 1L2Y)")
    else:
        query_path = Path(input_path)  # type: ignore[arg-type]
        print(f"  Query: {query_path}")

    if not query_path.exists():
        raise FileNotFoundError(f"Query structure not found: {query_path}")

    # Resolve databases
    db_list = _resolve_databases(databases)
    db_names = [alias for alias, _ in db_list]
    print(f"  Databases: {', '.join(db_names)}")

    # Run search against each database; merge results
    all_hits_raw: list[dict] = []

    with tempfile.TemporaryDirectory(prefix="foldseek_tmp_") as tmpdir:
        tmp_path = Path(tmpdir)

        for db_idx, (alias, foldseek_name) in enumerate(db_list):
            db_path = _get_db_path(alias, foldseek_name, db_cache)

            if not db_path.exists():
                print(
                    f"\n  WARNING: Database '{alias}' not found at {db_path}.\n"
                    f"  To download it, run:\n"
                    f"    foldseek databases {foldseek_name} {db_path} /tmp/foldseek_dl_tmp\n"
                    f"  Skipping '{alias}' database.\n"
                )
                continue

            print(f"  Searching against '{alias}' ({db_path.name})...")
            results_tsv = tmp_path / f"results_{alias}.tsv"

            proc = _run_foldseek_search(
                exe=exe,
                query_path=query_path,
                db_path=db_path,
                results_tsv=results_tsv,
                tmp_dir=tmp_path / f"fs_tmp_{db_idx}",
                threads=threads,
            )

            if proc.returncode != 0:
                err = (proc.stderr or proc.stdout or "").strip()
                raise RuntimeError(
                    f"foldseek exited with code {proc.returncode} for database '{alias}'.\n{err}"
                )

            hits = _parse_hits(results_tsv)
            # Tag each hit with its source database
            for h in hits:
                h["_database"] = alias
            all_hits_raw.extend(hits)
            print(f"    Found {len(hits)} raw hits in '{alias}'.")

    # Filter and rank
    hits_filtered = _filter_and_rank_hits(all_hits_raw, min_tmscore, max_hits)
    print(
        f"  Total: {len(all_hits_raw)} raw hits → "
        f"{len(hits_filtered)} after TM-score ≥ {min_tmscore} filter"
    )

    # Build the command string for reproducibility
    cmd = _build_cmd(query_path, output_dir, databases, min_tmscore, max_hits, demo)

    print(f"  Writing report to {output_dir}/")
    result = _generate_report(
        output_dir=output_dir,
        query_path=query_path,
        databases_used=db_names,
        hits_raw=all_hits_raw,
        hits_filtered=hits_filtered,
        min_tmscore=min_tmscore,
        max_hits=max_hits,
        cmd=cmd,
        foldseek_version=version,
        demo=demo,
    )

    print(f"\n  Report: {output_dir / 'report.md'}")
    print(f"  Full output: {output_dir}/")
    print(f"\n  {DISCLAIMER}")
    return result


def _build_cmd(
    query_path: Path,
    output_dir: Path,
    databases: str,
    min_tmscore: float,
    max_hits: int,
    demo: bool,
) -> str:
    parts = ["python skills/struct-predictor-foldseek/struct_predictor_foldseek.py"]
    if demo:
        parts.append("--demo")
    else:
        parts.append(f"--input {query_path}")
    parts.append(f"--databases {databases}")
    parts.append(f"--min-tmscore {min_tmscore}")
    parts.append(f"--max-hits {max_hits}")
    parts.append(f"--output {output_dir}")
    return " \\\n  ".join(parts)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Struct Predictor Foldseek — structural homology search"
    )
    parser.add_argument(
        "--input", "-i",
        help="Input structure file (CIF or PDB)",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--databases", "-d",
        default="pdb",
        help=(
            "Comma-separated database aliases to search. "
            "Available: pdb, afdb, esm. Default: pdb"
        ),
    )
    parser.add_argument(
        "--min-tmscore",
        type=float,
        default=DEFAULT_MIN_TMSCORE,
        help=f"Minimum TM-score for filtered hit table (default: {DEFAULT_MIN_TMSCORE})",
    )
    parser.add_argument(
        "--max-hits",
        type=int,
        default=DEFAULT_MAX_HITS,
        help=f"Maximum number of hits to report (default: {DEFAULT_MAX_HITS})",
    )
    parser.add_argument(
        "--db-cache",
        default=str(DEFAULT_DB_CACHE),
        help=f"Directory containing local Foldseek databases (default: {DEFAULT_DB_CACHE})",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of CPU threads for Foldseek (default: 1)",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run demo with Trp-cage miniprotein CIF (no input needed)",
    )
    return parser


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    if not args.input and not args.demo:
        parser.print_help()
        print("\nError: provide --input or --demo")
        sys.exit(1)

    print("Struct Predictor Foldseek")
    print("=" * 60)
    print()

    run_foldseek_search(
        input_path=Path(args.input) if args.input else None,
        output_dir=Path(args.output),
        databases=args.databases,
        min_tmscore=args.min_tmscore,
        max_hits=args.max_hits,
        db_cache=Path(args.db_cache),
        threads=args.threads,
        demo=args.demo,
    )


if __name__ == "__main__":
    main()
