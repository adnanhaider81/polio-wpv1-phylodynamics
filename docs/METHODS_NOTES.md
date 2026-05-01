# Methods Notes

This repository provides a lightweight WPV1 phylogeographic workflow that can be run on a laptop with synthetic data.

## Core workflow

1. Validate FASTA headers against metadata.
2. Assign each sample to a coarse geographic region.
3. Build a phylogeny with distance, FastTree, or IQ-TREE mode.
4. Time-scale the tree using sampling dates.
5. Infer discrete geographic histories using stochastic character mapping.
6. Summarize movement corridors, import/export balance, lineage persistence, and cluster status.

## Interpretation

Transition counts are model-based summaries from stochastic maps. They are not direct observations of infections or travel. Results should be interpreted alongside sampling intensity, surveillance design, and epidemiologic context.

## Recommended additions for manuscripts

- Root-to-tip regression and date-randomization tests.
- Sensitivity to region definitions.
- Sensitivity to environmental versus clinical sampling composition.
- IQ-TREE model-selection and bootstrap convergence checks.
- Comparison against posterior tree sets when available.
- Versioned analysis environment and successful clean-run log.
