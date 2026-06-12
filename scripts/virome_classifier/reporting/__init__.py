"""
Reporting module for taxonomy classification results.

Provides Kraken-format report generation:
- Per-read classification output (.kraken)
- Hierarchical taxonomic report (.kreport)
- Abundance table (.tsv)

Quick Start:
    >>> from virome_classifier import TaxonomyDB
    >>> from virome_classifier.reporting import write_all_outputs
    >>>
    >>> # After LCA classification
    >>> tax = TaxonomyDB.from_sqlite("taxonomy.db")
    >>> files = write_all_outputs(
    ...     results_df=lca_df,
    ...     tax=tax,
    ...     output_dir="./output",
    ...     sample_name="sample1"
    ... )
"""

from .kraken_writer import (
    write_all_outputs,
    write_kraken_output,
    generate_kraken_report,
    write_abundance_table,
    build_abundance_from_results,
    attach_taxonomy_columns,
)
from .rank_utils import get_rank_code
from .otu_table import (
    build_otu_table,
    build_otu_table_at_rank,
    filter_otu_table,
    export_otu_table,
    create_otu_pipeline,
)

__all__ = [
    # High-level API
    "write_all_outputs",
    "create_otu_pipeline",

    # Individual writers
    "write_kraken_output",
    "generate_kraken_report",
    "write_abundance_table",

    # OTU table generation
    "build_otu_table",
    "build_otu_table_at_rank",
    "filter_otu_table",
    "export_otu_table",

    # Utilities
    "build_abundance_from_results",
    "attach_taxonomy_columns",
    "get_rank_code",
]
