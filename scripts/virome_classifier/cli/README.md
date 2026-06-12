# virome_classifier CLI Tools

Command-line interface tools for the virome_classifier package.

## Available Commands

### 1. Taxonomy Database Setup

**Auto-download and build NCBI taxonomy database:**

```bash
python -m virome_classifier.cli.setup_taxonomy_db \
    --output taxonomy.db \
    --verbose
```

See: [QUICK_START_TAXONOMY_DB.md](QUICK_START_TAXONOMY_DB.md)

---

### 2. LCA-Based Classification

**Classify reads using Lowest Common Ancestor algorithm:**

```bash
python -m virome_classifier.cli.classify \
    --taxonomy taxonomy.db \
    --r1 sample_R1.tsv \
    --r2 sample_R2.tsv \
    --mask masked.bed \
    --output results/ \
    --sample MySample
```

**Features:**
- Fast LCA-based classification
- Multi-mapping read handling
- Kraken-compatible output
- Abundance tables

---

### 3. Coverage-Based Classification

**Classify using genome coverage thresholds:**

```bash
python -m virome_classifier.cli.classify_coverage \
    --taxonomy taxonomy.db \
    --r1 sample_R1.tsv \
    --r2 sample_R2.tsv \
    --mask masked.bed \
    --output results/ \
    --sample MySample \
    --min-coverage 0.05
```

**Features:**
- Breadth-of-coverage filtering
- Reduces false positives
- Segmented virus support
- Genome-size-corrected abundance

See: [classify_coverage.py](classify_coverage.py)

---

### 4. Build Taxonomy Database (Low-level)

**Build database from existing NCBI dump files:**

```bash
python -m virome_classifier.cli.build_taxonomy_db \
    --nodes nodes.dmp \
    --names names.dmp \
    --output taxonomy.db
```

**Note:** Most users should use `setup_taxonomy_db` instead (auto-downloads).

---

### 5. OTU Table Merging

**Merge multiple sample OTU tables:**

```bash
python -m virome_classifier.cli.merge_otu \
    --input sample1.tsv sample2.tsv sample3.tsv \
    --output merged_otu.tsv
```

---

### 6. Pipeline Profiling

**Analyze pipeline performance:**

```bash
python -m virome_classifier.cli.profile_pipeline \
    --input pipeline_metrics.json \
    --output profile_report.html
```

---

## Quick Start Workflow

### Complete Classification Pipeline

```bash
# 1. Setup taxonomy database (once)
python -m virome_classifier.cli.setup_taxonomy_db \
    --output taxonomy.db \
    --verbose

# 2. Classify your samples (LCA method)
python -m virome_classifier.cli.classify \
    --taxonomy taxonomy.db \
    --r1 sample_R1_mmseqs2.tsv \
    --r2 sample_R2_mmseqs2.tsv \
    --mask masked.bed \
    --output results/ \
    --sample MySample

# OR: Use coverage-based method for higher confidence
python -m virome_classifier.cli.classify_coverage \
    --taxonomy taxonomy.db \
    --r1 sample_R1_mmseqs2.tsv \
    --r2 sample_R2_mmseqs2.tsv \
    --mask masked.bed \
    --segment-info segments.csv \
    --output results/ \
    --sample MySample \
    --min-coverage 0.05 \
    --min-identity 0.85

# 3. Merge multiple samples (if needed)
python -m virome_classifier.cli.merge_otu \
    --input results/sample1_abundance.tsv \
            results/sample2_abundance.tsv \
            results/sample3_abundance.tsv \
    --output results/merged_abundance.tsv
```

## Output Files

### Classification Outputs

Both `classify` and `classify_coverage` produce:

```
results/
├── MySample.kraken              # Kraken-format classifications
├── MySample.kreport             # Kraken report
├── MySample_abundance.tsv       # Species abundance table
└── MySample_lca_result.csv      # Detailed read assignments
```

### Database Outputs

`setup_taxonomy_db` produces:

```
taxonomy.db                      # SQLite database
taxdump_temp/                    # Downloaded files (optional)
├── taxdump.tar.gz
├── nodes.dmp
└── names.dmp
```

## Documentation

- **Quick Start**: [QUICK_START_TAXONOMY_DB.md](QUICK_START_TAXONOMY_DB.md)
- **Complete Guide**: [TAXONOMY_DB_SETUP.md](TAXONOMY_DB_SETUP.md)
- **Coverage Method**: [../../COVERAGE_BASED_CLASSIFICATION.md](../../COVERAGE_BASED_CLASSIFICATION.md)

## Testing

```bash
# Test all CLI tools
cd <repo>/scripts
python3 test_taxonomy_db_setup.py

# Test imports
python3 test_virome_classifier_imports.py
```

## Getting Help

Each command has detailed help:

```bash
python -m virome_classifier.cli.setup_taxonomy_db --help
python -m virome_classifier.cli.classify --help
python -m virome_classifier.cli.classify_coverage --help
python -m virome_classifier.cli.build_taxonomy_db --help
python -m virome_classifier.cli.merge_otu --help
```

## Common Options

Most classification commands support:

- `--taxonomy`: Path to taxonomy database (required)
- `--output`: Output directory (required)
- `--sample`: Sample ID/name (required)
- `--verbose`: Enable detailed logging
- `--min-identity`: Minimum alignment identity (default: 0.8)
- `--min-length`: Minimum alignment length (default: 30)

## Integration

### Python API

```python
from virome_classifier.cli.setup_taxonomy_db import setup_taxonomy_database
from virome_classifier import TaxonomyDB

# Setup database
setup_taxonomy_database("taxonomy.db")

# Load and use
tax = TaxonomyDB.from_sqlite("taxonomy.db")
```

### Nextflow

```groovy
process VIROME_CLASSIFY {
    script:
    """
    python -m virome_classifier.cli.classify \\
        --taxonomy ${taxonomy_db} \\
        --r1 ${r1} \\
        --r2 ${r2} \\
        --output . \\
        --sample ${meta.id}
    """
}
```

See: [../../subworkflows/local/virome_classification.nf](../../subworkflows/local/virome_classification.nf)

## File Formats

### Input: MMseqs2 Results

Tab-separated format with columns:
```
query_id    target_id    identity    aln_length    qstart    qend    tstart    tend
```

### Input: Masked Regions (BED)

Standard BED format:
```
target_id    start    end
```

### Input: Segment Info (CSV)

For segmented viruses:
```
species_taxid,genome_size,segment_name,segment_size
123456,10000,segment_1,5000
123456,10000,segment_2,5000
```

### Output: Abundance Table

Tab-separated:
```
taxid    species_name    rank    read_count    percentage    lineage
```

### Output: Kraken Format

Standard Kraken output:
```
C    read_id    taxid    length    taxid_list
```

## Best Practices

1. **Setup taxonomy once**: Run `setup_taxonomy_db` once, reuse the database
2. **Update regularly**: Refresh taxonomy database monthly/quarterly
3. **Use coverage method**: For high-confidence results, especially with segmented viruses
4. **Keep temp files**: Use `--taxdump-dir` to preserve downloaded files for inspection
5. **Enable verbose**: Use `--verbose` when troubleshooting
6. **Batch processing**: Process multiple samples in parallel using Nextflow

## Performance Tips

- **Taxonomy DB**: Build once, reuse across all samples (~5 min build time)
- **Classification**: ~1-5 minutes per sample (depends on alignment size)
- **Memory**: ~500 MB for taxonomy loading + alignment size
- **Disk**: ~500 MB for taxonomy database + outputs

## Troubleshooting

### "ModuleNotFoundError: virome_classifier"

Make sure you're in the scripts directory:
```bash
cd <repo>/scripts
```

### "Database not found"

Build the database first:
```bash
python -m virome_classifier.cli.setup_taxonomy_db --output taxonomy.db
```

### "Permission denied"

Check output directory permissions or use a different path:
```bash
python -m virome_classifier.cli.classify --output ~/results/ ...
```

### "Download failed"

Use manual download or existing files:
```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
python -m virome_classifier.cli.build_taxonomy_db \
    --nodes nodes.dmp --names names.dmp --output taxonomy.db
```

## Support

For more help, see:
- Source code with detailed comments
- Test scripts for examples
- Documentation in each module
- Jupyter notebook tutorials in [../examples/](../examples/)
