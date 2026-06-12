<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nexvirome_logo_dark.png">
    <img alt="nexvirome" src="docs/images/nexvirome_logo_light.png">
  </picture>
</h1>

[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.20652876-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.20652876)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**nexvirome** is a Nextflow/nf-core pipeline for **viral metagenomics from
host-dominated samples** (e.g. sputum, respiratory swabs). It takes paired-end
short reads, removes human host reads, aligns the remainder against a curated
RefSeq viral database with [MMseqs2](https://github.com/soedinglab/MMseqs2), and
classifies hits with a masking- and breadth-aware filter that suppresses the
false positives typical of human-rich matrices (human/vector contamination and
conserved cross-mapping regions — retroviral *gag/pol/env* and oncogene homology,
flavivirus NS, herpes core genes, rRNA). It outputs per-sample virus calls and
merged sample × taxon OTU tables at genus / species / phage-host-genus level.

The pipeline steps are:

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter / quality trimming ([`fastp`](https://github.com/OpenGene/fastp) / [`cutadapt`](https://cutadapt.readthedocs.io/))
3. Human host-read removal ([`Bowtie2`](https://bowtie-bio.sourceforge.net/bowtie2/) against CHM13 T2T)
4. Viral alignment ([`MMseqs2`](https://github.com/soedinglab/MMseqs2) `easy-search` vs RefSeq viral)
5. **Virome classification** — mask conserved/contaminant regions, apply an unmasked-breadth gate and a per-taxon read floor, assign taxonomy (LCA / best-hit)
6. **OTU merge** — sample × taxon tables at genus / species / phage→host-genus
7. Aggregate QC ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

### 1. Prepare a samplesheet

`samplesheet.csv` — one row per sample (paired-end FASTQ):

```csv
sample,fastq_1,fastq_2
SAMPLE_01,/path/SAMPLE_01_R1.fastq.gz,/path/SAMPLE_01_R2.fastq.gz
SAMPLE_02,/path/SAMPLE_02_R1.fastq.gz,/path/SAMPLE_02_R2.fastq.gz
```

### 2. Download the reference databases

The viral MMseqs2 index, taxonomy SQLite, and masked-region file are distributed
via Zenodo (the human host genome is fetched separately from NCBI):

```bash
# Reference DB bundle (≈1.9 GB unpacked)
wget https://zenodo.org/records/20652876/files/nexvirome_db_v20260603.tar.gz
tar -xzf nexvirome_db_v20260603.tar.gz
DB=$PWD/release_db_v20260603

# Human host genome (CHM13 T2T) from NCBI
wget -O CHM13.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
gunzip CHM13.fna.gz
```

> The mask file `nexvirome_db.idx` is a self-describing container the pipeline
> reads directly via `--mask_bed`; no extra step is needed.

### 3. Run the pipeline

```bash
nextflow run nexvirome \
   -profile <docker/singularity/conda> \
   --input            samplesheet.csv \
   --outdir           results \
   --host_fasta       CHM13.fna \
   --mmseqs_database  $DB/mmseqs_db/viral_20260525 \
   --taxonomy_db      $DB/tax_seq_v20260526_MSL41.db \
   --mask_bed         $DB/nexvirome_db.idx
```

To start from already host-removed / trimmed reads, add `--skip_host_removal --skip_cutadapt`.

**Default classification parameters** are the locked "Method B" settings used in
the paper (best-hit assignment, unmasked-breadth ≥ 0.01, per-taxon read floor
n ≥ 3, relative-abundance gate off; HitQuality fident ≥ 0.85 / alnlen ≥ 60 /
qcov ≥ 0.5 / e ≤ 1e-3). Override any of them via params (e.g. `--min_unmasked_coverage`,
`--min_read_count`).

> [!WARNING]
> Provide pipeline parameters via the CLI or a Nextflow `-params-file`. Custom `-c`
> config files may set any configuration _**except parameters**_.

## Pipeline output

Results are written to `--outdir`:

| Path | Contents |
|---|---|
| `virome_classification/<sample>/` | per-sample virus calls: `<sample>_lca_classification.csv`, `.kreport`, `_abundance.tsv` |
| `virome_classification/otu_tables/` | merged sample × taxon tables: `otu_table_{genus,species,family}.csv` and `otu_table_phage_host.csv` (phage rolled up to host genus) |
| `multiqc/multiqc_report.html` | aggregate QC (FastQC, trimming, host-removal stats) |
| `pipeline_info/` | execution logs, software versions, resource usage |

See [`docs/output.md`](docs/output.md) for a full description of each file.

## Credits

nexvirome was developed by Sangchul Park (Korea University). It was built using
the [nf-core](https://nf-co.re/) pipeline template; nexvirome is an independent
project and is not an officially nf-core-curated pipeline.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use nexvirome, please cite the pipeline and its reference databases via the
Zenodo DOI:

> NexVirome reference databases (v20260603). Zenodo. doi: [10.5281/zenodo.20652876](https://doi.org/10.5281/zenodo.20652876)

An extensive list of references for the tools used by the pipeline can be found in
the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
