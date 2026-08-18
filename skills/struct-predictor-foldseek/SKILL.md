---
name: struct-predictor-foldseek
description: Structural homology search with Foldseek. Queries a predicted or experimental protein structure (CIF/PDB)
  against major structural databases (PDB, AlphaFold/Swiss-Prot, ESMAtlas), ranks hits by TM-score, and writes a
  markdown report with a similarity heatmap and hit table.
license: MIT
metadata:
  version: 0.1.0
  openclaw:
    requires:
      bins:
      - python3
      anyBins:
      - foldseek
    always: false
    emoji: 🔍
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    install:
    - kind: conda
      channel: bioconda
      package: foldseek
      bins:
      - foldseek
      comment: 'conda install -c bioconda -c conda-forge foldseek'
    - kind: pip
      package: matplotlib==3.9.4
    - kind: pip
      package: numpy==2.2.6
    - kind: pip
      package: requests==2.32.3
---

# Struct Predictor — Foldseek

You are the **Struct Predictor Foldseek**, a specialised agent for protein structural homology search using Foldseek.
Given a predicted or experimental structure file (CIF or PDB), you query it against major structural databases and
return ranked hits with TM-scores, RMSD, and sequence identity.

## Core Capabilities

1. **Structural Search**: Run `foldseek easy-search` against one or more target databases (PDB100, AF/Swiss-Prot, ESMAtlas)
2. **Hit Ranking**: Parse TSV output, rank by TM-score (descending), filter by configurable threshold
3. **Report Generation**: Markdown report with summary table, TM-score bar chart, and reproducibility bundle
4. **Demo Mode**: Uses the Trp-cage miniprotein CIF produced by the struct-predictor demo — no input file required
5. **Chaining**: Designed to receive the CIF output from `struct-predictor` and pass it directly into this skill

## CLI Reference

```bash
# Search a structure against PDB100 (default)
python skills/struct-predictor-foldseek/struct_predictor_foldseek.py \
  --input /path/to/structure.cif --output /tmp/foldseek_out

# Search against multiple databases
python skills/struct-predictor-foldseek/struct_predictor_foldseek.py \
  --input structure.pdb --databases pdb,afdb --output /tmp/foldseek_out

# Demo (Trp-cage CIF from demo_data/, no input needed)
python skills/struct-predictor-foldseek/struct_predictor_foldseek.py \
  --demo --output /tmp/foldseek_demo

# Filter by minimum TM-score
python skills/struct-predictor-foldseek/struct_predictor_foldseek.py \
  --input structure.cif --min-tmscore 0.5 --output /tmp/foldseek_out

# Limit number of reported hits
python skills/struct-predictor-foldseek/struct_predictor_foldseek.py \
  --input structure.cif --max-hits 20 --output /tmp/foldseek_out
```

### Plain Text Examples

Search a Boltz-2 prediction against PDB:

    python skills/struct-predictor-foldseek/struct_predictor_foldseek.py \
      --input output/predictions/MyProtein/MyProtein_model_0.cif \
      --output /tmp/foldseek_out

Run the built-in Trp-cage demo:

    python skills/struct-predictor-foldseek/struct_predictor_foldseek.py \
      --demo --output /tmp/foldseek_demo

Chain with struct-predictor (predict then search):

    python skills/struct-predictor/struct_predictor.py --demo --output /tmp/boltz_out
    python skills/struct-predictor-foldseek/struct_predictor_foldseek.py \
      --input /tmp/boltz_out/predictions/Trpcage/Trpcage_model_0.cif \
      --output /tmp/foldseek_out

## Workflow

```
Input CIF/PDB
     │
     ▼
foldseek easy-search <query> <db> results.tsv tmp/
     │
     ▼
Parse TSV → rank by TM-score → filter by --min-tmscore
     │
     ▼
Generate report.md + figures/tmscore_hits.png + result.json
```

## Supported Databases

| Alias     | Foldseek target name   | Description                                      |
|-----------|------------------------|--------------------------------------------------|
| `pdb`     | `pdb`                  | RCSB PDB (experimental structures)               |
| `afdb`    | `afdb50`               | AlphaFold DB / Swiss-Prot representative set     |
| `esm`     | `esmatlas`             | ESM Metagenomic Atlas (large, requires download) |

The `pdb` database is searched by default. Use `--databases pdb,afdb` to search multiple.

## Output Structure

```
output_dir/
  report.md                        # primary markdown report
  result.json                      # machine-readable summary
  figures/
    tmscore_hits.png               # horizontal bar chart of top hits by TM-score
  tables/
    hits.tsv                       # raw Foldseek output (all columns)
    hits_filtered.tsv              # filtered + ranked hits
  reproducibility/
    commands.sh                    # exact foldseek command used
    environment.txt                # foldseek version snapshot
```

## result.json Fields

```json
{
  "query": "MyProtein_model_0.cif",
  "databases": ["pdb"],
  "n_hits_raw": 250,
  "n_hits_filtered": 42,
  "min_tmscore_threshold": 0.3,
  "top_hit": {
    "target": "1L2Y_A",
    "tmscore": 0.97,
    "rmsd": 0.42,
    "seqid": 1.0,
    "evalue": 1.2e-15
  },
  "hits": [...],
  "demo": false,
  "disclaimer": "ClawBio is a research and educational tool..."
}
```

## TSV Column Reference

Foldseek `easy-search` outputs these columns (BLAST-like + structural):

| Column     | Description                              |
|------------|------------------------------------------|
| `query`    | Query structure ID                       |
| `target`   | Database hit identifier                  |
| `pident`   | Sequence identity (0–1)                  |
| `alnlen`   | Alignment length                         |
| `mismatch` | Mismatches                               |
| `gapopen`  | Gap openings                             |
| `qstart`   | Query alignment start                    |
| `qend`     | Query alignment end                      |
| `tstart`   | Target alignment start                   |
| `tend`     | Target alignment end                     |
| `evalue`   | E-value                                  |
| `bits`     | Bit score                                |
| `tmscore`  | TM-score (structural similarity, 0–1)    |
| `rmsd`     | RMSD of alignment (Å)                    |

## TM-Score Interpretation

| TM-score | Interpretation                              |
|----------|---------------------------------------------|
| > 0.9    | Near-identical fold                         |
| 0.7–0.9  | Very similar fold, same topology            |
| 0.5–0.7  | Likely same fold family                     |
| 0.3–0.5  | Possibly related, weak structural homology  |
| < 0.3    | Likely unrelated (random structural noise)  |

Default `--min-tmscore` threshold is **0.3**.

## Demo Data

| Item        | Value                                                             |
|-------------|-------------------------------------------------------------------|
| File        | `skills/struct-predictor-foldseek/demo_data/trpcage.cif`         |
| Structure   | Trp-cage miniprotein (synthetic CIF for offline demo)             |
| Length      | 20 residues                                                       |
| PDB ref     | 1L2Y                                                              |

The demo CIF is a minimal synthetic structure. In a real run, pass the CIF output from `struct-predictor`.

## Dependencies

```bash
# Foldseek — via conda (recommended)
conda install -c bioconda -c conda-forge foldseek

# Python packages (already in ClawBio requirements)
pip install matplotlib numpy requests
```

> **Note:** Foldseek is not available on PyPI. Conda (bioconda) is the supported install method.
> If you are using only a venv, install Foldseek system-wide or via conda into a separate env and
> ensure the `foldseek` binary is on `PATH`. The skill will detect and report if `foldseek` is missing.

## Chaining with struct-predictor

This skill is designed to run immediately after `struct-predictor`:

```
struct-predictor  →  CIF output  →  struct-predictor-foldseek  →  hit report
```

The Bio Orchestrator will chain these automatically when the user asks:
*"Predict the structure of this protein and find structural homologs."*

## Safety

- Foldseek databases are downloaded on first use to a local cache directory (`~/.foldseek_dbs/`).
- No patient/user sequence data is sent to any external server. `easy-search` runs fully locally.
- The skill checks for `foldseek` on PATH before running and exits with a clear error if missing.

## Citations

- van Kempen M et al. (2024) *Fast and accurate protein structure search with Foldseek*. Nature Biotechnology. doi:10.1038/s41587-023-01773-0
- Steinegger M, Söding J (2017) *MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets*. Nature Biotechnology. doi:10.1038/nbt.3988
