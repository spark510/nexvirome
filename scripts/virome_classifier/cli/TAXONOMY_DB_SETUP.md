# NCBI Taxonomy Database Setup Guide

## Overview

The `setup_taxonomy_db` module provides automated setup of NCBI taxonomy databases for virome_classifier. It handles downloading, extracting, and building SQLite databases from NCBI taxonomy dump files.

## Features

- **Auto-download**: Automatically downloads latest NCBI taxonomy dump from FTP
- **Progress tracking**: Visual progress bars for downloads
- **Smart caching**: Reuses existing files when possible
- **Verification**: Validates database after building
- **Flexible options**: Download-only mode, custom directories, etc.

## Quick Start

### 1. Auto-download and Build (Simplest)

```bash
cd <repo>/scripts

python -m virome_classifier.cli.setup_taxonomy_db \
    --output taxonomy.db
```

This will:
1. Download taxdump.tar.gz (~60 MB) from NCBI
2. Extract nodes.dmp and names.dmp
3. Build SQLite database
4. Verify the database

### 2. Specify Download Directory

```bash
python -m virome_classifier.cli.setup_taxonomy_db \
    --output taxonomy.db \
    --taxdump-dir ./ncbi_taxonomy
```

### 3. Use Existing Dump Files

If you already have NCBI taxonomy dump files:

```bash
python -m virome_classifier.cli.setup_taxonomy_db \
    --output taxonomy.db \
    --taxdump-dir /path/to/existing/taxdump
```

### 4. Download Only (Don't Build Database)

```bash
python -m virome_classifier.cli.setup_taxonomy_db \
    --download-only \
    --taxdump-dir ./taxdump
```

## Command-Line Options

```
--output, -o          Output SQLite database path (required unless --download-only)
--taxdump-dir, -d     Directory for taxdump files (default: temp dir next to output)
--download-only       Only download taxonomy files, don't build database
--batch-size, -b      Batch size for database inserts (default: 10000)
--force, -f           Force re-download and overwrite existing files
--no-keep-files       Remove downloaded files after building database
--verbose, -v         Enable verbose logging
```

## Complete Workflow

### Step 1: Setup Database

```bash
python -m virome_classifier.cli.setup_taxonomy_db \
    --output <repo>/resources/db/custom/taxonomy.db \
    --verbose
```

**Expected output:**
```
======================================================================
NCBI TAXONOMY DATABASE SETUP
======================================================================
======================================================================
STEP 1: DOWNLOAD NCBI TAXONOMY
======================================================================
📥 Downloading from: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
   Saving to: /home/share/.../taxdump_temp/taxdump.tar.gz
   [==================================================] 100.0% (60.5/60.5 MB)
✅ Download complete: taxdump.tar.gz

📦 Extracting: taxdump.tar.gz
   To: taxdump_temp/
✅ Extraction complete
   Extracted 2 files:
     - nodes.dmp (XXX.X MB)
     - names.dmp (XXX.X MB)

======================================================================
STEP 2: BUILD SQLITE DATABASE
======================================================================
📊 Building taxonomy database...
   Processing names.dmp...
   ✓ Loaded XXX,XXX taxonomy names

   Processing nodes.dmp...
   ✓ Loaded XXX,XXX taxonomy nodes

   Creating SQLite database: taxonomy.db
   ✓ Inserted XXX,XXX entries
   ✓ Database optimized

✅ Database built successfully

======================================================================
✅ SETUP COMPLETE
======================================================================
Database: /home/share/.../taxonomy.db
Size: XXX.XX MB

Usage:
  from virome_classifier import TaxonomyDB
  tax = TaxonomyDB.from_sqlite('/path/to/taxonomy.db')
```

### Step 2: Use in Python

```python
from virome_classifier import TaxonomyDB

# Load database
tax = TaxonomyDB.from_sqlite("/path/to/taxonomy.db")

# Get lineage
lineage = tax.get_lineage(10239)  # Viruses
print(f"Lineage: {lineage}")

# Get taxon name
name = tax.get_name(10239)
print(f"Name: {name}")

# Get rank
rank = tax.get_rank(10239)
print(f"Rank: {rank}")
```

### Step 3: Use in CLI

```bash
python -m virome_classifier.cli.classify \
    --taxonomy /path/to/taxonomy.db \
    --r1 sample_R1_mmseqs2_result.tsv \
    --r2 sample_R2_mmseqs2_result.tsv \
    --output results/
```

## Database Structure

The SQLite database contains a single table:

```sql
CREATE TABLE taxonomy (
    taxid INTEGER PRIMARY KEY,
    parent_taxid INTEGER,
    rank TEXT,
    scientific_name TEXT,
    FOREIGN KEY (parent_taxid) REFERENCES taxonomy(taxid)
);

CREATE INDEX idx_parent ON taxonomy(parent_taxid);
CREATE INDEX idx_rank ON taxonomy(rank);
CREATE INDEX idx_name ON taxonomy(scientific_name);
```

## File Locations

### Default Paths

```
.../resources/db/custom/
├── taxonomy.db              # SQLite database (output)
└── taxdump_temp/            # Temporary download directory
    ├── taxdump.tar.gz       # Downloaded archive
    ├── nodes.dmp            # Taxonomy tree structure
    ├── names.dmp            # Taxonomy names
    └── readme.txt           # NCBI documentation
```

### NCBI Source Files

- **Source URL**: `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`
- **nodes.dmp**: Taxonomic tree structure (taxid, parent_taxid, rank)
- **names.dmp**: Scientific and common names for each taxid

## Advanced Usage

### Force Rebuild

If the database already exists, use `--force` to rebuild:

```bash
python -m virome_classifier.cli.setup_taxonomy_db \
    --output taxonomy.db \
    --force \
    --verbose
```

### Custom Batch Size

For memory-constrained systems, reduce batch size:

```bash
python -m virome_classifier.cli.setup_taxonomy_db \
    --output taxonomy.db \
    --batch-size 5000
```

### Don't Keep Downloaded Files

To save disk space, remove taxdump files after building:

```bash
python -m virome_classifier.cli.setup_taxonomy_db \
    --output taxonomy.db \
    --no-keep-files
```

### Download to Specific Location

Keep downloaded files in a permanent location:

```bash
python -m virome_classifier.cli.setup_taxonomy_db \
    --output taxonomy.db \
    --taxdump-dir /data/ncbi_taxonomy \
    --verbose
```

## Troubleshooting

### Download Failed

```
❌ Download failed: [Error message]
```

**Solutions:**
1. Check internet connection
2. Check firewall settings
3. Try downloading manually from: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
4. Use `--taxdump-dir` with manually downloaded files

### Extraction Failed

```
❌ Extraction failed: [Error message]
```

**Solutions:**
1. Check disk space
2. Verify downloaded file is not corrupted
3. Try re-downloading with `--force`

### Database Build Failed

```
❌ Database build failed: [Error message]
```

**Solutions:**
1. Check disk space (needs ~200-500 MB)
2. Verify nodes.dmp and names.dmp exist
3. Check write permissions on output directory
4. Try with `--verbose` for detailed error info

### Permission Denied

```
❌ Permission denied: /path/to/output.db
```

**Solutions:**
1. Check write permissions on output directory
2. Use a different output path
3. Run with appropriate permissions

## Performance

### Download Time
- File size: ~60 MB compressed
- Time: 1-5 minutes (depends on connection speed)

### Build Time
- Processing: ~2-5 minutes
- Taxa processed: ~2.7 million entries
- Database size: ~200-500 MB

### System Requirements
- Disk space: ~500 MB free (download + database)
- RAM: ~500 MB for building
- Python: 3.8+

## Integration with Nextflow

To integrate into your Nextflow pipeline:

```groovy
process SETUP_TAXONOMY_DB {
    publishDir "${params.outdir}/databases", mode: 'copy'

    output:
    path "taxonomy.db", emit: taxonomy_db

    script:
    """
    python -m virome_classifier.cli.setup_taxonomy_db \\
        --output taxonomy.db \\
        --verbose
    """
}

workflow {
    SETUP_TAXONOMY_DB()

    // Use in classification
    VIROME_CLASSIFY(
        ch_reads,
        SETUP_TAXONOMY_DB.out.taxonomy_db,
        ch_mask_bed
    )
}
```

## Updating the Database

NCBI updates taxonomy regularly. To get the latest version:

```bash
# Re-download and rebuild
python -m virome_classifier.cli.setup_taxonomy_db \
    --output taxonomy.db \
    --force \
    --verbose
```

Recommended update frequency: Monthly or quarterly

## Related Documentation

- **Build script**: [build_taxonomy_db.py](build_taxonomy_db.py)
- **TaxonomyDB usage**: [../taxonomy/README.md](../taxonomy/README.md)
- **Classification CLI**: [classify.py](classify.py)
- **Coverage-based classification**: [classify_coverage.py](classify_coverage.py)

## Examples

### Example 1: Complete Setup for New Project

```bash
#!/bin/bash
# Setup taxonomy database for virome_classifier

DB_DIR="<repo>/resources/db/custom"
mkdir -p "$DB_DIR"

cd <repo>/scripts

python -m virome_classifier.cli.setup_taxonomy_db \
    --output "$DB_DIR/taxonomy.db" \
    --taxdump-dir "$DB_DIR/ncbi_taxdump" \
    --verbose

echo "Database ready at: $DB_DIR/taxonomy.db"
```

### Example 2: Quick Test

```bash
# Build in temp directory for testing
python -m virome_classifier.cli.setup_taxonomy_db \
    --output /tmp/test_taxonomy.db \
    --verbose

# Test the database
python3 -c "
from virome_classifier import TaxonomyDB
tax = TaxonomyDB.from_sqlite('/tmp/test_taxonomy.db')
print(f'Taxa count: {len(tax._taxonomy):,}')
print(f'Virus lineage: {tax.get_lineage(10239)}')
"
```

### Example 3: Download Only for Inspection

```bash
# Download files but don't build database
python -m virome_classifier.cli.setup_taxonomy_db \
    --download-only \
    --taxdump-dir ./inspect_taxonomy

# Inspect the files
head -20 ./inspect_taxonomy/nodes.dmp
head -20 ./inspect_taxonomy/names.dmp
```

## FAQs

**Q: How often should I update the taxonomy database?**
A: NCBI updates taxonomy regularly. Update monthly for active projects, quarterly for stable projects.

**Q: Can I use this with ICTV taxonomy?**
A: No, this tool specifically downloads NCBI taxonomy. For ICTV, see the notebook examples.

**Q: What's the difference between this and build_taxonomy_db.py?**
A: `setup_taxonomy_db.py` adds auto-download functionality. `build_taxonomy_db.py` requires manual download.

**Q: Can I run this in parallel?**
A: No, SQLite building is single-threaded. However, it's already quite fast (~2-5 minutes).

**Q: Does this work on Windows?**
A: Yes, but paths should use Windows format (e.g., `C:\path\to\taxonomy.db`)

**Q: Where can I find the original NCBI files?**
A: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/

## Support

For issues or questions:
1. Check this documentation
2. Run with `--verbose` for detailed logs
3. Check the [troubleshooting section](#troubleshooting)
4. Review the source code: [setup_taxonomy_db.py](setup_taxonomy_db.py)
