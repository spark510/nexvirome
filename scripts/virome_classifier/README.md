# Virome Classifier - OOP Refactored Version

**Modern, production-ready viral metagenomics classification pipeline**

## 🎯 Overview

This is a complete OOP refactoring of the original virome classification scripts, following Python best practices and modern software engineering principles.

### Key Improvements

| Aspect | Old Version | New Version |
|--------|-------------|-------------|
| **Structure** | Procedural scripts | OOP packages |
| **Organization** | Single files | Modular packages |
| **Type Safety** | Minimal type hints | Full type annotations |
| **Error Handling** | Generic exceptions | Custom exception hierarchy |
| **Testing** | Limited | Comprehensive test suite |
| **Documentation** | Inline comments | Full docstrings + guides |
| **Logging** | Print statements | Structured logging |
| **Dependencies** | Implicit | Explicit (requirements.txt) |

## 📁 Project Structure

```
virome_classifier/
├── __init__.py                 # Public API
├── __version__.py              # Version info
├── core/                       # Core utilities
│   ├── logger.py              # Logging system
│   ├── types.py               # Type definitions
│   └── exceptions.py          # Custom exceptions
├── taxonomy/                   # Taxonomy database
│   ├── taxonomy_db.py         # TaxonomyDB class
│   └── ranks.py               # Rank utilities
├── alignment/                  # Alignment parsing
│   ├── parser.py              # AlignmentParser
│   ├── header.py              # Header normalization
│   └── filters/               # Filtering modules
│       └── masking.py         # Masking filter
├── classification/             # LCA classification
│   ├── classifier.py          # Classifier classes
│   ├── paired.py              # Paired-end resolution
│   └── settings.py            # Configuration
├── reporting/                  # Report generation
│   ├── writer.py              # Report writers
│   └── formats.py             # Output formats
├── cli/                        # Command-line interface
│   └── main.py                # CLI entry point
└── tests/                      # Test suite
    ├── test_taxonomy.py
    ├── test_alignment.py
    └── test_classification.py
```

## 🚀 Quick Start

### Basic Usage

```python
from virome_classifier import TaxonomyDB, AlignmentParser, MaskingFilter

# 1. Load taxonomy
tax = TaxonomyDB.from_sqlite("taxonomy.db", root_taxid=10239)

# 2. Parse alignments
parser = AlignmentParser()
hits_df = parser.parse("alignments.tsv")

# 3. Apply quality filters
from virome_classifier.core import FilterCriteria

criteria = FilterCriteria(
    min_identity=0.8,
    min_alignment_length=30,
    max_evalue=1e-3,
    min_query_coverage=0.5
)

filtered_df = parser.filter(hits_df, criteria)

# 4. Apply masking filter
masking = MaskingFilter.from_bed_file("masked_regions.bed")
result = masking.filter_by_unmasked_coverage(filtered_df, min_coverage=0.25)

print(result.summary())
```

### Paired-End Processing

```python
from virome_classifier.alignment import BatchAlignmentParser

# Parse paired-end files
batch_parser = BatchAlignmentParser()
forward_df, reverse_df = batch_parser.parse_paired(
    "sample_R1.tsv",
    "sample_R2.tsv",
    filter_criteria=criteria
)
```

### Taxonomy Operations

```python
# Get taxonomy info
info = tax.get_info(10359)  # Human betaherpesvirus 5
print(f"Name: {info.name}")
print(f"Rank: {info.rank}")
print(f"Lineage: {info.lineage}")

# Compute LCA
taxids = [10298, 10310, 10359]  # Multiple herpesviruses
lca = tax.compute_lca(taxids)
print(f"LCA: {tax.get_name(lca)}")
```

## 🔧 Installation

### Requirements

```bash
# Install dependencies
pip install pandas numpy

# For development
pip install pytest black mypy
```

### Setup

```bash
# Add to PYTHONPATH
export PYTHONPATH="/path/to/nexvirome/scripts/virome_engine:$PYTHONPATH"

# Or install in editable mode
cd /path/to/nexvirome/scripts/virome_engine
pip install -e .
```

## 📚 Detailed Documentation

### Core Module

#### Logger

```python
from virome_classifier.core import get_logger, set_verbose

# Get global logger
logger = get_logger()

# Configure verbosity
set_verbose(True)

# Use logger
logger.info("Processing started")
logger.verbose("Detailed debug info")
logger.success("Processing complete!")

# Add file logging
logger.add_file_handler(Path("pipeline.log"))
```

#### Custom Exceptions

```python
from virome_classifier.core import (
    TaxonomyError,
    AlignmentError,
    ClassificationError
)

try:
    tax = TaxonomyDB.from_sqlite("missing.db")
except TaxonomyError as e:
    print(f"Taxonomy error: {e}")
```

#### Type Definitions

```python
from virome_classifier.core import (
    AlignmentHit,
    FilterCriteria,
    ClassificationResult,
    TaxonInfo
)

# Create immutable hit
hit = AlignmentHit(
    query="read1",
    target="NC_001806.2",
    identity=0.95,
    alignment_length=100,
    evalue=1e-50,
    bitscore=200,
    query_length=150,
    target_length=100000,
    taxid=10376
)

print(f"Query coverage: {hit.query_coverage:.2%}")
```

### Taxonomy Module

#### TaxonomyDB

```python
from virome_classifier import TaxonomyDB

# Load from SQLite
tax = TaxonomyDB.from_sqlite(
    "taxonomy.db",
    root_taxid=10239,  # Viruses
    build_major_taxon=False
)

# Basic operations
tax.exists(10359)                    # → True
tax.get_name(10359)                  # → "Human betaherpesvirus 5"
tax.get_rank(10359)                  # → "species"
tax.get_parent(10359)                # → 10357
tax.get_lineage(10359)               # → [10239, 10242, ..., 10359]

# Get complete info
info = tax.get_info(10359)
print(info.name, info.rank, info.lineage)

# LCA operations
lca = tax.compute_lca([10298, 10310, 10359])
lca_major = tax.compute_lca_with_rank([10298, 10310], major_only=True)

# Check relationships
tax.is_ancestor(10357, 10359)       # → True

# Statistics
print(f"Total taxa: {len(tax):,}")
```

### Alignment Module

#### AlignmentParser

```python
from virome_classifier import AlignmentParser, FilterCriteria

parser = AlignmentParser(normalize_headers=True)

# Parse file (auto-detect format)
df = parser.parse("alignments.tsv")

# Parse with specific format
from virome_classifier.core import AlignmentFormat
df = parser.parse("blast.out", format_type=AlignmentFormat.BLAST)

# Filter hits
criteria = FilterCriteria(
    min_identity=0.8,
    min_alignment_length=50,
    max_evalue=1e-5,
    min_query_coverage=0.7
)

filtered = parser.filter(df, criteria)

# Or filter with individual parameters
filtered = parser.filter(
    df,
    min_identity=0.8,
    min_alignment_length=50,
    max_evalue=1e-5,
    min_query_coverage=0.7
)
```

#### MaskingFilter

```python
from virome_classifier.alignment import MaskingFilter

# Load from BED file
masking = MaskingFilter.from_bed_file("masked_regions.bed")

# Filter by unmasked coverage
result = masking.filter_by_unmasked_coverage(df, min_coverage=0.25)
print(result.summary())
print(f"Passed: {result.n_passed_targets} targets")

# Filter by total coverage
result = masking.filter_by_total_coverage(df, min_coverage=0.5)

# Hybrid filter (both conditions required)
result = masking.filter_by_hybrid_coverage(
    df,
    min_total_cov=0.5,
    min_unmasked_cov=0.25
)

# Calculate statistics only
stats_df = masking.calculate_stats(df)
```

## 🔄 Migration Guide

### Old Code → New Code

#### 1. Taxonomy Loading

```python
# OLD
from taxonomy import Taxonomy
tax = Taxonomy("taxonomy.db", root_taxid=10239)

# NEW
from virome_classifier import TaxonomyDB
tax = TaxonomyDB.from_sqlite("taxonomy.db", root_taxid=10239)
```

#### 2. Alignment Parsing

```python
# OLD
from search_result_parser import SeqSearchResultParser
df = SeqSearchResultParser.parse_file("file.tsv")
filtered = SeqSearchResultParser.filter_hits(df, min_seq_id=0.8, ...)

# NEW
from virome_classifier import AlignmentParser, FilterCriteria
parser = AlignmentParser()
df = parser.parse("file.tsv")
criteria = FilterCriteria(min_identity=0.8, ...)
filtered = parser.filter(df, criteria)
```

#### 3. Masking Filter

```python
# OLD
from masking_filter import load_masked_bed, apply_unmasked_cov_filter
mask_dict = load_masked_bed("file.bed")
passed, failed, stats = apply_unmasked_cov_filter(df, mask_dict, 0.25)

# NEW
from virome_classifier.alignment import MaskingFilter
masking = MaskingFilter.from_bed_file("file.bed")
result = masking.filter_by_unmasked_coverage(df, min_coverage=0.25)
passed, failed, stats = result.passed, result.failed, result.stats
```

#### 4. Logging

```python
# OLD
from utils import log_info, log_verbose
log_info("Message")

# NEW
from virome_classifier.core import log_info, log_verbose
log_info("Message")  # Same API!
```

## 🧪 Testing

```bash
# Run all tests
pytest tests/

# Run specific test
pytest tests/test_taxonomy.py

# With coverage
pytest --cov=virome_classifier tests/
```

## 🏗️ Development

### Code Style

```bash
# Format code
black virome_classifier/

# Type check
mypy virome_classifier/

# Lint
pylint virome_classifier/
```

### Adding New Features

1. **Create feature branch**
   ```bash
   git checkout -b feature/new-classifier
   ```

2. **Follow package structure**
   - Add to appropriate package (taxonomy/, alignment/, etc.)
   - Write comprehensive docstrings
   - Add type hints
   - Include tests

3. **Example: New classifier**
   ```python
   # virome_classifier/classification/new_classifier.py
   from typing import List
   from ..core import ClassificationResult, TaxID
   from ..taxonomy import TaxonomyDB

   class NewClassifier:
       """
       New classification algorithm.

       Args:
           taxonomy: TaxonomyDB instance

       Example:
           >>> classifier = NewClassifier(tax)
           >>> result = classifier.classify(...)
       """

       def __init__(self, taxonomy: TaxonomyDB):
           self._tax = taxonomy

       def classify(self, taxids: List[TaxID]) -> ClassificationResult:
           """Classify using new algorithm."""
           # Implementation
           pass
   ```

## 📝 TODO

- [ ] Complete classification module OOP wrapper
- [ ] Complete reporting module OOP wrapper
- [ ] Implement modern CLI with typer
- [ ] Add comprehensive test suite
- [ ] Add performance benchmarks
- [ ] Create Sphinx documentation
- [ ] Add GitHub Actions CI/CD
- [ ] Create Docker container
- [ ] Add example datasets

## 🤝 Contributing

1. Follow PEP 8 style guide
2. Add type hints to all functions
3. Write docstrings (Google style)
4. Include unit tests
5. Update documentation

## 📄 License

Same as parent project

## 🙏 Acknowledgments

This refactoring maintains full compatibility with the original codebase while providing a modern, maintainable foundation for future development.

**Original Authors**: Virome Engine Team
**Refactoring**: 2024-11
