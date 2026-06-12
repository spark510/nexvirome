#!/usr/bin/env python3
"""
Build Taxonomy Database from NCBI Taxonomy Dump Files

This script creates a SQLite database from NCBI taxonomy dump files:
- nodes.dmp: Contains taxonomy tree structure (taxid, parent_taxid, rank)
- names.dmp: Contains taxonomy names (taxid, name, name_class)

The resulting database is optimized for:
- Fast LCA (Lowest Common Ancestor) queries
- Efficient lineage traversal
- Minimal storage footprint

Usage:
    python -m virome_classifier.cli.build_taxonomy_db \\
        --nodes nodes.dmp \\
        --names names.dmp \\
        --output taxonomy.db
"""
from __future__ import annotations

import argparse
import sqlite3
import sys
from pathlib import Path
from typing import Generator, Dict, Tuple

from ..core import log_info, log_verbose, set_verbose


# ==========================================
# Database Schema
# ==========================================

SCHEMA_SQL = """
-- Main taxonomy table
CREATE TABLE IF NOT EXISTS ncbi_taxonomy (
    taxid           INTEGER PRIMARY KEY,
    parent_taxid    INTEGER NOT NULL,
    rank            TEXT NOT NULL,
    scientific_name TEXT NOT NULL
);

-- Indexes for fast queries
CREATE INDEX IF NOT EXISTS idx_parent ON ncbi_taxonomy(parent_taxid);
CREATE INDEX IF NOT EXISTS idx_rank ON ncbi_taxonomy(rank);
CREATE INDEX IF NOT EXISTS idx_name ON ncbi_taxonomy(scientific_name);

-- Metadata table
CREATE TABLE IF NOT EXISTS metadata (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
"""


# ==========================================
# Parsing Functions
# ==========================================

def parse_names_dmp(names_file: str) -> Dict[int, str]:
    """
    Parse names.dmp file to extract scientific names.

    Format of names.dmp:
        taxid | name | unique_name | name_class |

    Args:
        names_file: Path to names.dmp file

    Returns:
        Dictionary mapping taxid to scientific name
    """
    log_info(f"📄 Parsing scientific names from {names_file}...")

    names_map = {}
    line_count = 0

    with open(names_file, 'r', encoding='utf-8') as f:
        for line in f:
            line_count += 1

            # Split by delimiter '|'
            parts = [p.strip() for p in line.split('|')]

            if len(parts) < 4:
                continue

            # Only keep scientific names
            if parts[3] == 'scientific name':
                try:
                    taxid = int(parts[0])
                    name = parts[1]
                    names_map[taxid] = name
                except ValueError:
                    log_verbose(f"  Warning: Invalid taxid in line {line_count}")
                    continue

    log_info(f"✅ Found {len(names_map):,} scientific names from {line_count:,} lines")
    return names_map


def parse_nodes_dmp(
    nodes_file: str,
    names_map: Dict[int, str]
) -> Generator[Tuple[int, int, str, str], None, None]:
    """
    Parse nodes.dmp file and generate taxonomy records.

    Format of nodes.dmp:
        taxid | parent_taxid | rank | embl_code | division_id | ...

    Args:
        nodes_file: Path to nodes.dmp file
        names_map: Dictionary mapping taxid to scientific name

    Yields:
        Tuples of (taxid, parent_taxid, rank, scientific_name)
    """
    log_info(f"📄 Parsing taxonomy tree from {nodes_file}...")

    line_count = 0
    valid_count = 0

    with open(nodes_file, 'r', encoding='utf-8') as f:
        for line in f:
            line_count += 1

            # Split by tab-delimited '|'
            parts = [p.strip() for p in line.split('\t|\t')]

            try:
                taxid = int(parts[0])
                parent_taxid = int(parts[1])
                rank = parts[2]
            except (ValueError, IndexError):
                log_verbose(f"  Warning: Invalid format in line {line_count}")
                continue

            # Get scientific name from names_map
            name = names_map.get(taxid, "N/A")

            valid_count += 1
            yield (taxid, parent_taxid, rank, name)

    log_info(f"✅ Parsed {valid_count:,} taxonomy nodes from {line_count:,} lines")


# ==========================================
# Database Operations
# ==========================================

def create_database_schema(conn: sqlite3.Connection) -> None:
    """
    Create database schema with tables and indexes.

    Args:
        conn: SQLite database connection
    """
    log_info("🗄️  Creating database schema...")

    cursor = conn.cursor()
    cursor.executescript(SCHEMA_SQL)
    conn.commit()

    log_info("✅ Schema created successfully")


def bulk_insert_taxonomy(
    conn: sqlite3.Connection,
    data_generator: Generator[Tuple[int, int, str, str], None, None],
    batch_size: int = 10000
) -> int:
    """
    Bulk insert taxonomy data into database.

    Uses batched inserts for better performance.

    Args:
        conn: SQLite database connection
        data_generator: Generator yielding taxonomy records
        batch_size: Number of records per batch

    Returns:
        Total number of records inserted
    """
    log_info(f"💾 Inserting taxonomy data (batch size: {batch_size:,})...")

    cursor = conn.cursor()
    insert_query = """
        INSERT OR IGNORE INTO ncbi_taxonomy
        (taxid, parent_taxid, rank, scientific_name)
        VALUES (?, ?, ?, ?)
    """

    total_count = 0
    batch = []

    for record in data_generator:
        batch.append(record)

        if len(batch) >= batch_size:
            cursor.executemany(insert_query, batch)
            conn.commit()
            total_count += len(batch)
            log_verbose(f"  Inserted {total_count:,} records...")
            batch = []

    # Insert remaining records
    if batch:
        cursor.executemany(insert_query, batch)
        conn.commit()
        total_count += len(batch)

    log_info(f"✅ Inserted {total_count:,} taxonomy records")
    return total_count


def add_metadata(conn: sqlite3.Connection, metadata: Dict[str, str]) -> None:
    """
    Add metadata to database.

    Args:
        conn: SQLite database connection
        metadata: Dictionary of metadata key-value pairs
    """
    log_info("📋 Adding metadata...")

    cursor = conn.cursor()
    insert_query = "INSERT OR REPLACE INTO metadata (key, value) VALUES (?, ?)"

    for key, value in metadata.items():
        cursor.execute(insert_query, (key, value))

    conn.commit()
    log_info(f"✅ Added {len(metadata)} metadata entries")


def optimize_database(conn: sqlite3.Connection) -> None:
    """
    Optimize database for faster queries.

    Runs ANALYZE and VACUUM commands.

    Args:
        conn: SQLite database connection
    """
    log_info("⚡ Optimizing database...")

    cursor = conn.cursor()

    # Update query planner statistics
    log_verbose("  Running ANALYZE...")
    cursor.execute("ANALYZE")

    # Reclaim unused space and defragment
    log_verbose("  Running VACUUM...")
    cursor.execute("VACUUM")

    conn.commit()
    log_info("✅ Database optimized")


def verify_database(conn: sqlite3.Connection, expected_count: int) -> bool:
    """
    Verify database integrity and content.

    Args:
        conn: SQLite database connection
        expected_count: Expected number of records

    Returns:
        True if verification passed
    """
    log_info("🔍 Verifying database...")

    cursor = conn.cursor()

    # Check record count
    cursor.execute("SELECT COUNT(*) FROM ncbi_taxonomy")
    actual_count = cursor.fetchone()[0]

    if actual_count != expected_count:
        log_info(f"⚠️  Warning: Expected {expected_count:,} records, found {actual_count:,}")
        return False

    # Check for root node (taxid=1)
    cursor.execute("SELECT COUNT(*) FROM ncbi_taxonomy WHERE taxid = 1")
    has_root = cursor.fetchone()[0] > 0

    if not has_root:
        log_info("⚠️  Warning: Root node (taxid=1) not found")

    # Check for viruses root (taxid=10239)
    cursor.execute("SELECT COUNT(*) FROM ncbi_taxonomy WHERE taxid = 10239")
    has_viruses = cursor.fetchone()[0] > 0

    if not has_viruses:
        log_info("⚠️  Warning: Viruses root (taxid=10239) not found")

    # Sample query test
    cursor.execute("""
        SELECT taxid, scientific_name, rank
        FROM ncbi_taxonomy
        WHERE taxid IN (1, 10239, 9606)
        ORDER BY taxid
    """)
    sample_results = cursor.fetchall()

    log_info("  Sample records:")
    for taxid, name, rank in sample_results:
        log_info(f"    {taxid}: {name} ({rank})")

    log_info(f"✅ Verification complete: {actual_count:,} records")
    return True


# ==========================================
# Main Pipeline
# ==========================================

def build_taxonomy_database(
    nodes_file: str,
    names_file: str,
    output_file: str,
    batch_size: int = 10000,
    force: bool = False,
    verbose: bool = False,
) -> None:
    """
    Complete pipeline to build taxonomy database.

    Args:
        nodes_file: Path to nodes.dmp
        names_file: Path to names.dmp
        output_file: Output SQLite database path
        batch_size: Batch size for bulk inserts
        force: Overwrite existing database
        verbose: Enable verbose logging
    """
    import datetime

    if verbose:
        set_verbose(True)

    # Check input files
    nodes_path = Path(nodes_file)
    names_path = Path(names_file)
    output_path = Path(output_file)

    if not nodes_path.exists():
        raise FileNotFoundError(f"nodes.dmp not found: {nodes_file}")
    if not names_path.exists():
        raise FileNotFoundError(f"names.dmp not found: {names_file}")

    # Check if output exists
    if output_path.exists() and not force:
        raise FileExistsError(
            f"Database already exists: {output_file}\n"
            "Use --force to overwrite"
        )

    # Remove existing database if force
    if output_path.exists() and force:
        log_info(f"🗑️  Removing existing database: {output_file}")
        output_path.unlink()

    # Create output directory
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Connect to database
    log_info(f"🔗 Connecting to database: {output_file}")
    conn = sqlite3.connect(str(output_path))

    # Enable optimizations
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA cache_size=10000")
    conn.execute("PRAGMA temp_store=MEMORY")

    try:
        # Step 1: Create schema
        create_database_schema(conn)

        # Step 2: Parse names
        names_map = parse_names_dmp(str(names_path))

        # Step 3: Parse nodes and insert data
        data_generator = parse_nodes_dmp(str(nodes_path), names_map)
        record_count = bulk_insert_taxonomy(conn, data_generator, batch_size)

        # Step 4: Add metadata
        metadata = {
            "source": "NCBI Taxonomy",
            "build_date": datetime.datetime.now().isoformat(),
            "nodes_file": nodes_path.name,
            "names_file": names_path.name,
            "record_count": str(record_count),
        }
        add_metadata(conn, metadata)

        # Step 5: Optimize
        optimize_database(conn)

        # Step 6: Verify
        verify_database(conn, record_count)

        log_info(f"\n✅ Database built successfully: {output_file}")
        log_info(f"   Records: {record_count:,}")
        log_info(f"   Size: {output_path.stat().st_size / (1024**2):.2f} MB")

    except Exception as e:
        log_info(f"\n❌ Error building database: {e}")
        conn.close()
        if output_path.exists():
            output_path.unlink()
        raise

    finally:
        conn.close()


# ==========================================
# CLI
# ==========================================

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Build taxonomy database from NCBI dump files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  %(prog)s --nodes nodes.dmp --names names.dmp --output taxonomy.db

  # With custom batch size
  %(prog)s --nodes nodes.dmp --names names.dmp --output taxonomy.db \\
      --batch-size 50000

  # Force overwrite existing database
  %(prog)s --nodes nodes.dmp --names names.dmp --output taxonomy.db --force

Download NCBI Taxonomy:
  wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
  tar -xzf taxdump.tar.gz
        """,
    )

    parser.add_argument(
        "--nodes",
        "-n",
        type=str,
        required=True,
        help="Path to nodes.dmp file",
    )
    parser.add_argument(
        "--names",
        "-m",
        type=str,
        required=True,
        help="Path to names.dmp file",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output SQLite database path",
    )
    parser.add_argument(
        "--batch-size",
        "-b",
        type=int,
        default=10000,
        help="Batch size for inserts (default: 10000)",
    )
    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="Overwrite existing database",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )

    return parser.parse_args()


def main():
    """Main entry point."""
    args = parse_args()

    print("=" * 70)
    print("TAXONOMY DATABASE BUILDER")
    print("=" * 70)

    try:
        build_taxonomy_database(
            nodes_file=args.nodes,
            names_file=args.names,
            output_file=args.output,
            batch_size=args.batch_size,
            force=args.force,
            verbose=args.verbose,
        )
        return 0

    except Exception as e:
        log_info(f"\n❌ Fatal error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
