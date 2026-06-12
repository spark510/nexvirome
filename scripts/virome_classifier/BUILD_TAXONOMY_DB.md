# Building Taxonomy Database

Create a SQLite taxonomy database from NCBI taxonomy dump files.

## Quick Start

```bash
# Download NCBI taxonomy dump
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz

# Build database
python -m virome_classifier.cli.build_taxonomy_db \
    --nodes nodes.dmp \
    --names names.dmp \
    --output taxonomy.db \
    --verbose
```

## Prerequisites

### 1. Download NCBI Taxonomy Dump

The NCBI taxonomy dump is available at:
https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/

```bash
# Download and extract
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz

# Required files:
# - nodes.dmp: Taxonomy tree structure
# - names.dmp: Taxonomy names
```

### 2. Required Files

| File | Description | Size |
|------|-------------|------|
| `nodes.dmp` | Taxonomy tree (taxid, parent, rank) | ~200 MB |
| `names.dmp` | Taxonomy names (taxid, name, name_class) | ~300 MB |

## Usage

### Basic Build

```bash
python -m virome_classifier.cli.build_taxonomy_db \
    --nodes nodes.dmp \
    --names names.dmp \
    --output taxonomy.db
```

### With Options

```bash
python -m virome_classifier.cli.build_taxonomy_db \
    --nodes nodes.dmp \
    --names names.dmp \
    --output taxonomy.db \
    --batch-size 50000 \
    --force \
    --verbose
```

## Command-Line Options

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--nodes FILE` | `-n` | Required | Path to nodes.dmp |
| `--names FILE` | `-m` | Required | Path to names.dmp |
| `--output FILE` | `-o` | Required | Output database path |
| `--batch-size N` | `-b` | 10000 | Insert batch size |
| `--force` | `-f` | False | Overwrite existing DB |
| `--verbose` | `-v` | False | Verbose logging |

## File Formats

### nodes.dmp Format

Tab-delimited file with fields separated by ` | `:

```
taxid | parent_taxid | rank | embl_code | division_id | ...
1     | 1            | no rank | ...
2     | 131567       | superkingdom | ...
```

Fields used:
- **taxid**: Taxonomy ID
- **parent_taxid**: Parent taxonomy ID
- **rank**: Taxonomic rank (species, genus, family, etc.)

### names.dmp Format

Tab-delimited file with fields separated by ` | `:

```
taxid | name | unique_name | name_class
1     | root |             | scientific name
2     | Bacteria |         | scientific name
```

Fields used:
- **taxid**: Taxonomy ID
- **name**: Taxon name
- **name_class**: Type of name (only "scientific name" is used)

## Database Schema

### ncbi_taxonomy Table

```sql
CREATE TABLE ncbi_taxonomy (
    taxid           INTEGER PRIMARY KEY,
    parent_taxid    INTEGER NOT NULL,
    rank            TEXT NOT NULL,
    scientific_name TEXT NOT NULL
);

CREATE INDEX idx_parent ON ncbi_taxonomy(parent_taxid);
CREATE INDEX idx_rank ON ncbi_taxonomy(rank);
CREATE INDEX idx_name ON ncbi_taxonomy(scientific_name);
```

### metadata Table

```sql
CREATE TABLE metadata (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
```

Metadata includes:
- `source`: Data source (NCBI Taxonomy)
- `build_date`: Database build timestamp
- `nodes_file`: Source nodes.dmp filename
- `names_file`: Source names.dmp filename
- `record_count`: Total number of records

## Build Process

The build process consists of 6 steps:

### 1. Create Schema
```
🗄️  Creating database schema...
✅ Schema created successfully
```

Creates tables and indexes.

### 2. Parse Names
```
📄 Parsing scientific names from names.dmp...
✅ Found 2,419,564 scientific names from 3,850,123 lines
```

Extracts scientific names from names.dmp (filters out synonyms).

### 3. Parse Nodes & Insert
```
📄 Parsing taxonomy tree from nodes.dmp...
✅ Parsed 2,419,564 taxonomy nodes from 2,419,564 lines
💾 Inserting taxonomy data (batch size: 10,000)...
✅ Inserted 2,419,564 taxonomy records
```

Parses nodes.dmp and inserts data in batches for performance.

### 4. Add Metadata
```
📋 Adding metadata...
✅ Added 5 metadata entries
```

Stores build information in metadata table.

### 5. Optimize
```
⚡ Optimizing database...
  Running ANALYZE...
  Running VACUUM...
✅ Database optimized
```

Optimizes query planner and reclaims unused space.

### 6. Verify
```
🔍 Verifying database...
  Sample records:
    1: root (no rank)
    10239: Viruses (acellular root)
    9606: Homo sapiens (species)
✅ Verification complete: 2,419,564 records
```

Verifies record count and checks key taxonomy nodes.

## Output

After successful build:

```
✅ Database built successfully: taxonomy.db
   Records: 2,419,564
   Size: 456.78 MB
```

### Database File

The output is a SQLite database file that can be used with:
- `virome_classifier` classification pipeline
- Any SQLite-compatible tool
- Python `sqlite3` module

## Performance Tips

### 1. Batch Size

Larger batch sizes are faster but use more memory:

```bash
# Default (balanced)
--batch-size 10000

# Fast (more memory)
--batch-size 50000

# Slow (less memory)
--batch-size 1000
```

### 2. Storage

Use fast storage (SSD) for better performance:
- Input files: SSD preferred
- Output database: SSD strongly recommended

### 3. System Resources

Minimum requirements:
- **RAM**: 2 GB
- **Disk space**: 2 GB free (1 GB for output DB)
- **Time**: 2-5 minutes on modern hardware

## Verification

### Check Database

```bash
# Query database
sqlite3 taxonomy.db

# Check record count
SELECT COUNT(*) FROM ncbi_taxonomy;

# Check root nodes
SELECT taxid, scientific_name, rank
FROM ncbi_taxonomy
WHERE taxid IN (1, 10239, 9606);

# Check metadata
SELECT * FROM metadata;
```

### Test with virome_classifier

```python
from virome_classifier import TaxonomyDB

# Load database
tax = TaxonomyDB.from_sqlite("taxonomy.db", root_taxid=10239)

# Test queries
print(f"Loaded {len(tax):,} taxa")
print(tax.get_name(10359))  # Should print: Human betaherpesvirus 5
print(tax.get_rank(10359))  # Should print: no rank

# Test lineage
lineage = tax.get_lineage(10359)
for tid in lineage[:5]:
    print(f"  {tid}: {tax.get_name(tid)} ({tax.get_rank(tid)})")
```

## Troubleshooting

### FileNotFoundError: nodes.dmp not found

Make sure you've extracted the tar file:
```bash
tar -xzf taxdump.tar.gz
ls -lh nodes.dmp names.dmp
```

### FileExistsError: Database already exists

Use `--force` to overwrite:
```bash
python -m virome_classifier.cli.build_taxonomy_db \
    --nodes nodes.dmp --names names.dmp \
    --output taxonomy.db --force
```

### Warning: Root node (taxid=1) not found

The nodes.dmp file may be corrupted. Re-download:
```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz
```

### Database is very large

This is normal. Full NCBI taxonomy creates ~450-500 MB database.
To reduce size, consider filtering to specific clades.

## Advanced Usage

### Virus-Only Database

To create a smaller database with only viral taxonomy:

```python
# After building full database, extract virus subtree
from virome_classifier import TaxonomyDB
import sqlite3

# Load full taxonomy
full_tax = TaxonomyDB.from_sqlite("taxonomy_full.db", root_taxid=1)

# Create new database
conn = sqlite3.connect("taxonomy_viruses.db")
cursor = conn.cursor()

# Create schema
cursor.executescript("""
CREATE TABLE ncbi_taxonomy (
    taxid INTEGER PRIMARY KEY,
    parent_taxid INTEGER NOT NULL,
    rank TEXT NOT NULL,
    scientific_name TEXT NOT NULL
);
""")

# Get all virus taxids (starting from 10239)
virus_lineages = set()
conn_full = sqlite3.connect("taxonomy_full.db")
cursor_full = conn_full.cursor()

cursor_full.execute("""
WITH RECURSIVE virus_tree(taxid, parent_taxid, rank, scientific_name) AS (
  SELECT taxid, parent_taxid, rank, scientific_name
  FROM ncbi_taxonomy
  WHERE taxid = 10239
  UNION ALL
  SELECT t.taxid, t.parent_taxid, t.rank, t.scientific_name
  FROM ncbi_taxonomy t
  JOIN virus_tree v ON t.parent_taxid = v.taxid
)
SELECT * FROM virus_tree;
""")

# Insert into new database
for row in cursor_full.fetchall():
    cursor.execute(
        "INSERT INTO ncbi_taxonomy VALUES (?, ?, ?, ?)",
        row
    )

conn.commit()
conn.close()
conn_full.close()

print("✅ Virus-only database created")
```

## Updating Taxonomy

NCBI taxonomy is updated regularly. To update:

```bash
# Download latest taxonomy
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf taxdump.tar.gz

# Rebuild database
python -m virome_classifier.cli.build_taxonomy_db \
    --nodes nodes.dmp \
    --names names.dmp \
    --output taxonomy_latest.db \
    --force
```

## See Also

- [CLI Usage](CLI_USAGE.md) - Main classification pipeline
- [OTU Table](OTU_TABLE_USAGE.md) - OTU table generation
- [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) - Official documentation
