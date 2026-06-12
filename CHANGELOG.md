# nexvirome: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - 2026-06-12

Initial public release of nexvirome — a Nextflow/nf-core pipeline for viral
metagenomics from host-dominated samples. Built using the
[nf-core](https://nf-co.re/) template.

### `Added`

- End-to-end pipeline: FastQC → fastp/cutadapt → Bowtie2 host removal (CHM13 T2T)
  → MMseqs2 viral alignment → masking/breadth-aware classification → OTU merge → MultiQC.
- `virome_classifier` package and local modules (classify / coverage / otu_merge).
- Masking- and breadth-aware false-positive control for human-rich matrices
  (human/vector contamination; retroviral *gag/pol/env* and oncogene homology;
  flavivirus NS; herpes core; rRNA), shipped as an encrypted masked-region file.
- Locked "Method B" classification defaults (best-hit, unmasked-breadth ≥ 0.01,
  per-taxon read floor n ≥ 3, rel-abundance gate off; fident 0.85 / alnlen 60 /
  qcov 0.5 / e 1e-3), wired 1:1 from Nextflow params to the classifier CLI.
- Phage → host-genus rollup as an additional OTU table (`merge_otu --phage-host`).
- Reference databases distributed via Zenodo ([10.5281/zenodo.20652876](https://doi.org/10.5281/zenodo.20652876)).

### `Fixed`

- conda profile: expose `virome_classifier` via `PYTHONPATH` instead of a relative
  `pip -e`, which broke under Nextflow's task working directory.

### `Dependencies`

- Python ≥ 3.9 (pandas, numpy, psutil) for the classifier modules.
