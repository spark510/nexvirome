# Quick Start: NCBI Taxonomy Database Setup

## TL;DR - One Command Setup

```bash
cd <repo>/scripts

# Auto-download and build NCBI taxonomy database
python -m virome_classifier.cli.setup_taxonomy_db \
    --output ../resources/db/custom/taxonomy.db \
    --verbose
```

That's it! The script will:
1. Download taxdump.tar.gz from NCBI (~60 MB)
2. Extract nodes.dmp and names.dmp
3. Build SQLite database
4. Verify it works

Time: ~5-10 minutes total

## What You Get

After running the command, you'll have:

```
resources/db/custom/
├── taxonomy.db              # ← Your SQLite database (~200-500 MB)
└── taxdump_temp/            # Downloaded files (can be deleted)
    ├── taxdump.tar.gz
    ├── nodes.dmp
    └── names.dmp
```

## Using the Database

### In Python

```python
from virome_classifier import TaxonomyDB

# Load database
tax = TaxonomyDB.from_sqlite("resources/db/custom/taxonomy.db")

# Get information
name = tax.get_name(10239)        # "Viruses"
rank = tax.get_rank(10239)        # "acellular root"
lineage = tax.get_lineage(10239)  # [10239]

print(f"Loaded {len(tax):,} taxa")
```

### In CLI

```bash
# LCA-based classification
python -m virome_classifier.cli.classify \
    --taxonomy resources/db/custom/taxonomy.db \
    --r1 sample_R1.tsv \
    --r2 sample_R2.tsv \
    --output results/

# Coverage-based classification
python -m virome_classifier.cli.classify_coverage \
    --taxonomy resources/db/custom/taxonomy.db \
    --r1 sample_R1.tsv \
    --r2 sample_R2.tsv \
    --output results/
```

### In Nextflow

```groovy
process VIROME_CLASSIFY {
    input:
    tuple val(meta), path(r1), path(r2)
    path taxonomy_db  // Pass taxonomy.db here

    script:
    """
    python -m virome_classifier.cli.classify \\
        --taxonomy ${taxonomy_db} \\
        --r1 ${r1} \\
        --r2 ${r2} \\
        --output .
    """
}
```

## Common Use Cases

### 1. First Time Setup

```bash
# Create database directory
mkdir -p <repo>/resources/db/custom

# Download and build
python -m virome_classifier.cli.setup_taxonomy_db \
    --output <repo>/resources/db/custom/taxonomy.db \
    --verbose
```

### 2. Update Existing Database

NCBI updates taxonomy regularly. To get the latest:

```bash
# Force rebuild
python -m virome_classifier.cli.setup_taxonomy_db \
    --output taxonomy.db \
    --force \
    --verbose
```

### 3. Download Only (Inspect First)

```bash
# Just download, don't build yet
python -m virome_classifier.cli.setup_taxonomy_db \
    --download-only \
    --taxdump-dir ./ncbi_taxdump

# Inspect the files
head -20 ./ncbi_taxdump/nodes.dmp
head -20 ./ncbi_taxdump/names.dmp

# Then build when ready
python -m virome_classifier.cli.build_taxonomy_db \
    --nodes ./ncbi_taxdump/nodes.dmp \
    --names ./ncbi_taxdump/names.dmp \
    --output taxonomy.db
```

### 4. Save Disk Space

```bash
# Don't keep downloaded files after building
python -m virome_classifier.cli.setup_taxonomy_db \
    --output taxonomy.db \
    --no-keep-files
```

## Verify It Works

```bash
# Run test script
python3 test_taxonomy_db_setup.py
```

Expected output:
```
✅ PASS: Module Imports
✅ PASS: CLI Help
✅ PASS: Existing Database
✅ PASS: File Locations

Total: 4/4 tests passed

🎉 All tests passed!
```

## Troubleshooting

### "Download failed"

**Problem**: Network issues or firewall blocking NCBI FTP

**Solution**:
```bash
# Download manually
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

# Extract
tar -xzf taxdump.tar.gz

# Build from local files
python -m virome_classifier.cli.build_taxonomy_db \
    --nodes nodes.dmp \
    --names names.dmp \
    --output taxonomy.db
```

### "Permission denied"

**Problem**: No write access to output directory

**Solution**: Use a different output path or fix permissions

```bash
# Use home directory
python -m virome_classifier.cli.setup_taxonomy_db \
    --output ~/taxonomy.db
```

### "Disk space full"

**Problem**: Not enough space (~500 MB needed)

**Solution**: Clean up space or use different location

```bash
# Check available space
df -h .

# Use different location with more space
python -m virome_classifier.cli.setup_taxonomy_db \
    --output /path/with/more/space/taxonomy.db
```

## Full Documentation

For complete documentation, see:
- [TAXONOMY_DB_SETUP.md](TAXONOMY_DB_SETUP.md) - Complete guide
- [build_taxonomy_db.py](build_taxonomy_db.py) - Low-level building
- [setup_taxonomy_db.py](setup_taxonomy_db.py) - High-level automation

## Example Session

```bash
$ cd <repo>/scripts

$ python -m virome_classifier.cli.setup_taxonomy_db --output test.db --verbose

======================================================================
NCBI TAXONOMY DATABASE SETUP
======================================================================
======================================================================
STEP 1: DOWNLOAD NCBI TAXONOMY
======================================================================
📥 Downloading from: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
   Saving to: taxdump_temp/taxdump.tar.gz
   [==================================================] 100.0% (60.5/60.5 MB)
✅ Download complete: taxdump.tar.gz

📦 Extracting: taxdump.tar.gz
✅ Extraction complete
   Extracted 2 files:
     - nodes.dmp (234.5 MB)
     - names.dmp (189.2 MB)

======================================================================
STEP 2: BUILD SQLITE DATABASE
======================================================================
📊 Building taxonomy database...
   Processing names.dmp...
   ✓ Loaded 2,749,286 taxonomy names

   Processing nodes.dmp...
   ✓ Loaded 2,749,286 taxonomy nodes

   Creating SQLite database: test.db
   ✓ Inserted 2,749,286 entries in 1.2 minutes
   ✓ Database optimized

✅ Database built successfully

======================================================================
✅ SETUP COMPLETE
======================================================================
Database: test.db
Size: 423.45 MB

Taxonomy files kept in: taxdump_temp/

Usage:
  from virome_classifier import TaxonomyDB
  tax = TaxonomyDB.from_sqlite('test.db')


$ python3 -c "from virome_classifier import TaxonomyDB; tax = TaxonomyDB.from_sqlite('test.db'); print(f'Taxa: {len(tax):,}')"

[INFO] 📚 Loading taxonomy from test.db...
[INFO] ✅ Loaded 2,749,286 taxa
Taxa: 2,749,286
```

## Next Steps

After setting up the taxonomy database:

1. **Test it**: Run [test_taxonomy_db_setup.py](../../test_taxonomy_db_setup.py)

2. **Use in classification**:
   - LCA method: [classify.py](classify.py)
   - Coverage-based: [classify_coverage.py](classify_coverage.py)

3. **Integrate into Nextflow**: See [virome_classification.nf](../../../subworkflows/local/virome_classification.nf)

4. **Try the tutorial**: [coverage_based_tutorial_kr.ipynb](../examples/coverage_based_tutorial_kr.ipynb)

Happy classifying! 🧬🦠
