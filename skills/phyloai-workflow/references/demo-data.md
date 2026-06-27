# Demo Data

Demo paths are resolved from the installed package, not hardcoded.

## End-to-End Dataset

- `phyloai/demo_data/end_to_end/raw/`
- Small amino-acid FASTA files for walkthroughs.

## Per-Step Data

- `phyloai/demo_data/per_step/aligned/`
- `phyloai/demo_data/per_step/trimmed/`
- `phyloai/demo_data/per_step/concat/`
- `phyloai/demo_data/per_step/gen_trees/`

Demo runs must write outputs to a user run directory, never into the demo data directory.
