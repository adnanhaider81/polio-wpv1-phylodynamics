# polio-wpv1-phylodynamics

[![DOI](https://zenodo.org/badge/1226527658.svg)](https://zenodo.org/badge/latestdoi/1226527658)
[![R demo](https://github.com/adnanhaider81/polio-wpv1-phylodynamics/actions/workflows/r-demo.yml/badge.svg)](https://github.com/adnanhaider81/polio-wpv1-phylodynamics/actions/workflows/r-demo.yml)

Reproducible R workflow for exploratory WPV1 phylogeographic analysis. The repository includes code, documentation, and a small synthetic dataset so the pipeline can be tested without restricted surveillance data.

No raw surveillance sequences, sample identifiers, or line-list metadata are included.

## Portfolio quick view

This repository shows how VP1 sequence data and metadata can be turned into reproducible exploratory phylogeographic summaries: tree building, time scaling, stochastic character mapping, movement tables, lineage persistence summaries, and optional maps. The included synthetic data make the workflow reviewable without exposing restricted surveillance records.

```mermaid
flowchart LR
  A["VP1 alignment"] --> B["Metadata QC"]
  B --> C["Regional grouping"]
  C --> D["Tree building"]
  D --> E["Time scaling"]
  E --> F["Ancestral-state mapping"]
  F --> G["Movement summaries"]
  G --> H["Tables and figures"]
```

## Public repository checklist

| Item | Status |
| --- | --- |
| README, license, citation metadata | Present |
| Reproducible environment | R install script and `DESCRIPTION` |
| Tests or smoke checks | Synthetic demo run |
| Example or synthetic data | `data/synthetic/` |
| Documentation | `docs/` plus README |
| Data privacy note | Present; restricted sequences and line lists are excluded |
| GitHub Actions badge | Present |
| Container recipe | `Dockerfile` |
| Zenodo DOI | [10.5281/zenodo.20257428](https://doi.org/10.5281/zenodo.20257428) |

## What the pipeline does

- Reads a VP1 FASTA alignment and sample metadata.
- Assigns samples to a coarse Pakistan/Afghanistan regional scheme.
- Builds a tree using a pure-R distance mode, FastTree, or IQ-TREE.
- Time-scales the tree with `treedater`.
- Runs stochastic character mapping with `phytools`.
- Summarizes inferred movement corridors, import/export balance, local transmission lineages, cluster persistence, and root-to-tip diagnostics.
- Optionally builds GIS maps when boundary files are available.

## Quick start

Install R packages:

```bash
Rscript scripts/install_packages.R
```

Run the synthetic demo:

```bash
Rscript run_pipeline.R --root . --mode distance --nsim 5 --seed 1
```

The main outputs are written to:

- `results/tables/`
- `results/figures/`
- `results/trees/`

You can also run:

```bash
make install
make demo
```

## Input files

The default synthetic inputs are:

- `data/synthetic/synthetic_alignment.fasta`
- `data/synthetic/synthetic_metadata.csv`

Metadata must contain these columns:

| Column | Description |
| --- | --- |
| `FILENAME` | FASTA header / strain name |
| `Country` | Country code such as `PAK` or `AFG` |
| `StateProv` | Province or state |
| `locality` | District, city, or local label |
| `onsetdate` | Sample date in `YYYY-MM-DD` format |
| `IsENV` | Environmental sample indicator |
| `Cluster` | Program or genetic cluster label |

CSV, TSV, and Excel metadata files are supported. Excel input requires the optional `readxl` package.

## Tree modes

Distance mode is easiest for testing:

```bash
Rscript run_pipeline.R --root . --mode distance --nsim 5
```

FastTree mode requires `FastTree` on `PATH`:

```bash
Rscript run_pipeline.R --root . --mode fasttree --nsim 25
```

IQ-TREE mode requires `iqtree2`, `iqtree3`, or `iqtree` on `PATH`:

```bash
Rscript run_pipeline.R --root . --mode iqtree --ufboot 1000 --nsim 25
```

## Optional maps

Maps require additional geospatial packages and downloaded boundaries:

```bash
Rscript scripts/install_packages.R --optional
Rscript run_pipeline.R --root . --mode distance --make_maps --download_boundaries
```

Boundary downloads are cached under `data/gis/`, which is ignored by git.

## Data policy

This repository is intended as a public methods template. Restricted sequence data, operational metadata, and generated private analysis products should remain outside the repository. Use `data/raw/` locally for private inputs; it is ignored by git.

## Methods status

This workflow is suitable for exploratory analysis and reproducible demonstration. Manuscript-grade phylodynamic inference should add temporal-signal testing, sampling-bias sensitivity analyses, model comparison, uncertainty checks, and an explicit data-use statement.

## Citation

Please cite the archived Zenodo release when using this workflow:

Haider, S. A. (2026). polio-wpv1-phylodynamics (v0.1.1). Zenodo. https://doi.org/10.5281/zenodo.20257428

The all-version Zenodo concept DOI is https://doi.org/10.5281/zenodo.20257427.
