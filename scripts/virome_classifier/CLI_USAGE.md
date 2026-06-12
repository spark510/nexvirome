# Virome Classifier CLI Usage

Complete command-line interface for viral metagenomics classification.

## Installation

```bash
cd <repo>/scripts
export PYTHONPATH="${PYTHONPATH}:$(pwd)"
```

## Quick Start

### Paired-End Mode

```bash
python -m virome_classifier \
    --r1 sample_R1.result \
    --r2 sample_R2.result \
    --taxonomy /path/to/taxonomy.db \
    --mask /path/to/masked_regions.bed \
    --output results/ \
    --sample sample1 \
    --verbose
```

### Single-End Mode

```bash
python -m virome_classifier \
    --input sample.result \
    --taxonomy /path/to/taxonomy.db \
    --mask /path/to/masked_regions.bed \
    --output results/ \
    --sample sample1
```

## Command-Line Options

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--r1 FILE` | R1 alignment file (paired-end mode) |
| `--r2 FILE` | R2 alignment file (requires --r1) |
| `--input FILE` | Single alignment file (single-end mode) |
| `--taxonomy FILE` | Taxonomy database (SQLite) |
| `--mask FILE` | Masked regions BED file |
| `--output DIR` | Output directory |
| `--sample NAME` | Sample name for output files |

### Filtering Parameters

| Argument | Default | Description |
|----------|---------|-------------|
| `--min-identity FLOAT` | 0.8 | Minimum alignment identity (0-1) |
| `--min-length INT` | 50 | Minimum alignment length (bp) |
| `--max-evalue FLOAT` | 1e-3 | Maximum E-value |
| `--min-query-coverage FLOAT` | 0.5 | Minimum query coverage (0-1) |
| `--min-unmasked-coverage FLOAT` | 0.1 | Minimum unmasked genome coverage |

### Taxonomy Options

| Argument | Default | Description |
|----------|---------|-------------|
| `--root-taxid INT` | 10239 | Root taxonomy ID (10239=Viruses) |
| `--virus-root` | True | Use Viruses (10239) as root in reports |

### Other Options

| Argument | Description |
|----------|-------------|
| `--verbose, -v` | Enable verbose logging |
| `--no-kraken` | Skip Kraken report generation |
| `--help, -h` | Show help message |

## Example Workflows

### 1. Basic Paired-End Analysis

```bash
python -m virome_classifier \
    --r1 MagNA_1_R1.result \
    --r2 MagNA_1_R2.result \
    --taxonomy <repo>/resources/db/custom/tax_seq.db \
    --mask <repo>/notebooks/refseq_masked_250901.bed \
    --output ./results \
    --sample MagNA_1
```

### 2. Strict Filtering

```bash
python -m virome_classifier \
    --r1 sample_R1.result \
    --r2 sample_R2.result \
    --taxonomy taxonomy.db \
    --mask masked.bed \
    --output results/ \
    --sample sample1 \
    --min-identity 0.95 \
    --min-length 100 \
    --min-unmasked-coverage 0.25 \
    --verbose
```

### 3. Single-End with No Kraken Reports

```bash
python -m virome_classifier \
    --input sample.result \
    --taxonomy taxonomy.db \
    --mask masked.bed \
    --output results/ \
    --sample sample1 \
    --no-kraken
```

## Output Files

The pipeline generates the following files in the output directory:

### Classification Results

| File | Format | Description |
|------|--------|-------------|
| `{sample}_lca_classification.csv` | CSV | LCA classification for each query |
| `{sample}_final_hits.csv` | CSV | All filtered alignment hits |
| `{sample}_coverage_stats.csv` | CSV | Coverage statistics per target |

### Kraken Reports (unless --no-kraken)

| File | Format | Description |
|------|--------|-------------|
| `{sample}.kraken` | TSV | Per-read Kraken classification |
| `{sample}.kreport` | TSV | Hierarchical Kraken report |
| `{sample}.abundance.tsv` | TSV | Flat abundance table |

## File Format Examples

### LCA Classification CSV

```csv
query,lca_taxid,lca_name,lca_rank,qlen,read_count,n_hits,n_unique_taxids,all_taxids
read_001,10359,Human betaherpesvirus 5,no rank,151,1,5,2,"10359,10358"
read_002,10310,Human alphaherpesvirus 2,no rank,151,1,3,1,"10310"
```

### Kraken Output

```
C  read_001  10359  151  10239|2731341|...|10359
C  read_002  10310  151  10239|2731341|...|10310
U  read_003  0      151  0
```

### Kraken Report

```
100.00  141246  114     -   10239   Viruses
 84.01  118579  118465  -   -         unclassified
 15.91   22451     0    D   2731341     Duplodnaviria
  7.92   11180   11180   G   10358       Cytomegalovirus
```

## Alternative Execution Methods

### Method 1: Python Module (Recommended)

```bash
python -m virome_classifier [options]
```

### Method 2: Direct Script Execution

```bash
python <repo>/scripts/virome_classifier/cli/classify.py [options]
```

### Method 3: Standalone Script

```bash
# Make symlink (one time only)
ln -s <repo>/scripts/virome_classifier/cli/virome-classify ~/bin/

# Run
virome-classify [options]
```

## Pipeline Steps

The CLI performs the following steps automatically:

1. **Load Taxonomy**: Load NCBI taxonomy database
2. **Parse Alignments**: Parse MMseqs2/BLAST alignment files
3. **Quality Filter**: Apply identity, length, E-value, coverage filters
4. **Masking Filter**: Remove hits to repetitive/low-complexity regions
5. **LCA Classification**: Compute Lowest Common Ancestor for each query
6. **Export Results**: Save classification and coverage results
7. **Generate Reports**: Create Kraken-format reports (optional)

## Troubleshooting

### ImportError: No module named 'virome_classifier'

Make sure PYTHONPATH is set:
```bash
export PYTHONPATH="<repo>/scripts:${PYTHONPATH}"
```

### FileNotFoundError: Taxonomy database not found

Check the path to your taxonomy database:
```bash
ls -lh /path/to/taxonomy.db
```

### Empty output files

Try enabling verbose mode to see what's happening:
```bash
python -m virome_classifier --verbose [other options]
```

### Out of memory

For large datasets, consider:
- Processing samples separately
- Increasing system memory
- Using stricter filtering thresholds

## Performance Tips

1. **Use SSD storage** for database files
2. **Increase filtering stringency** for large datasets
3. **Process samples in parallel** (separate runs)
4. **Monitor memory usage** with `htop` or `top`

## Citation

If you use this pipeline, please cite:
- virome_classifier
- MMseqs2 (if using MMseqs2 alignments)
- Kraken2 (if using Kraken report format)
- NCBI Taxonomy

## Support

For issues or questions:
- Check the main README
- Review example notebooks in `examples/`
- Check GitHub issues
