#!/usr/bin/env python3
"""
Setup Taxonomy Database - Complete Automation

This script automates the entire process:
1. Download NCBI taxonomy dump files
2. Extract the files
3. Build SQLite database
4. Verify the database

Usage:
    # Auto-download and build
    python -m virome_classifier.cli.setup_taxonomy_db \\
        --output taxonomy.db

    # Use existing dump files
    python -m virome_classifier.cli.setup_taxonomy_db \\
        --taxdump-dir /path/to/taxdump \\
        --output taxonomy.db

    # Download only
    python -m virome_classifier.cli.setup_taxonomy_db \\
        --download-only \\
        --taxdump-dir ./taxdump
"""

import argparse
import sys
import urllib.request
import tarfile
import shutil
from pathlib import Path
from typing import Optional

from ..core import log_info, log_verbose, set_verbose
from .build_taxonomy_db import build_taxonomy_database


# NCBI taxonomy dump URL
NCBI_TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"


def download_with_progress(url: str, output_file: Path) -> None:
    """
    Download file with progress bar.

    Args:
        url: URL to download
        output_file: Output file path
    """
    log_info(f"📥 Downloading from: {url}")
    log_info(f"   Saving to: {output_file}")

    def progress_hook(block_num, block_size, total_size):
        downloaded = block_num * block_size
        if total_size > 0:
            percent = min(100, (downloaded / total_size) * 100)
            mb_downloaded = downloaded / (1024 ** 2)
            mb_total = total_size / (1024 ** 2)

            # Progress bar
            bar_length = 50
            filled = int(bar_length * percent / 100)
            bar = '=' * filled + '-' * (bar_length - filled)

            print(f"\r   [{bar}] {percent:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)", end='')

            if downloaded >= total_size:
                print()  # New line when complete

    try:
        urllib.request.urlretrieve(url, str(output_file), progress_hook)
        log_info(f"✅ Download complete: {output_file}")
    except Exception as e:
        log_info(f"❌ Download failed: {e}")
        raise


def extract_taxdump(tar_file: Path, output_dir: Path) -> None:
    """
    Extract taxdump.tar.gz file.

    Args:
        tar_file: Path to taxdump.tar.gz
        output_dir: Directory to extract to
    """
    log_info(f"📦 Extracting: {tar_file}")
    log_info(f"   To: {output_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        with tarfile.open(tar_file, 'r:gz') as tar:
            # Extract only the files we need
            members = []
            for member in tar.getmembers():
                if member.name in ['nodes.dmp', 'names.dmp', 'readme.txt']:
                    members.append(member)
                    log_verbose(f"  Extracting: {member.name}")

            tar.extractall(path=output_dir, members=members)

        log_info(f"✅ Extraction complete")

        # List extracted files
        extracted = list(output_dir.glob('*.dmp'))
        log_info(f"   Extracted {len(extracted)} files:")
        for f in extracted:
            size_mb = f.stat().st_size / (1024 ** 2)
            log_info(f"     - {f.name} ({size_mb:.1f} MB)")

    except Exception as e:
        log_info(f"❌ Extraction failed: {e}")
        raise


def download_ncbi_taxonomy(output_dir: Path, force: bool = False) -> tuple[Path, Path]:
    """
    Download and extract NCBI taxonomy dump files.

    Args:
        output_dir: Directory to save files
        force: Re-download if files exist

    Returns:
        Tuple of (nodes.dmp path, names.dmp path)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    tar_file = output_dir / "taxdump.tar.gz"
    nodes_file = output_dir / "nodes.dmp"
    names_file = output_dir / "names.dmp"

    # Check if files already exist
    if nodes_file.exists() and names_file.exists() and not force:
        log_info(f"✅ Taxonomy files already exist:")
        log_info(f"   - {nodes_file}")
        log_info(f"   - {names_file}")
        log_info(f"   Use --force to re-download")
        return nodes_file, names_file

    # Download
    if not tar_file.exists() or force:
        download_with_progress(NCBI_TAXDUMP_URL, tar_file)
    else:
        log_info(f"✅ Archive already exists: {tar_file}")

    # Extract
    extract_taxdump(tar_file, output_dir)

    # Verify extracted files
    if not nodes_file.exists():
        raise FileNotFoundError(f"nodes.dmp not found after extraction: {nodes_file}")
    if not names_file.exists():
        raise FileNotFoundError(f"names.dmp not found after extraction: {names_file}")

    # Optionally remove tar file to save space
    # tar_file.unlink()

    return nodes_file, names_file


def setup_taxonomy_database(
    output_db: str,
    taxdump_dir: Optional[str] = None,
    download_only: bool = False,
    batch_size: int = 10000,
    force: bool = False,
    keep_files: bool = True,
    verbose: bool = False,
) -> None:
    """
    Complete setup of taxonomy database.

    Args:
        output_db: Output database path
        taxdump_dir: Directory with/for taxdump files
        download_only: Only download, don't build database
        batch_size: Batch size for database inserts
        force: Overwrite existing files
        keep_files: Keep downloaded files after building
        verbose: Enable verbose logging
    """
    if verbose:
        set_verbose(True)

    output_path = Path(output_db)

    # Determine taxdump directory
    if taxdump_dir:
        taxdump_path = Path(taxdump_dir)
    else:
        # Use temporary directory next to output
        taxdump_path = output_path.parent / "taxdump_temp"

    # Step 1: Get taxonomy files
    log_info("\n" + "="*70)
    log_info("STEP 1: DOWNLOAD NCBI TAXONOMY")
    log_info("="*70)

    nodes_file, names_file = download_ncbi_taxonomy(
        taxdump_path,
        force=force
    )

    if download_only:
        log_info("\n✅ Download complete (--download-only mode)")
        log_info(f"   Files saved to: {taxdump_path}")
        return

    # Step 2: Build database
    log_info("\n" + "="*70)
    log_info("STEP 2: BUILD SQLITE DATABASE")
    log_info("="*70)

    build_taxonomy_database(
        nodes_file=str(nodes_file),
        names_file=str(names_file),
        output_file=str(output_path),
        batch_size=batch_size,
        force=force,
        verbose=verbose,
    )

    # Step 3: Cleanup
    if not keep_files and taxdump_path.name == "taxdump_temp":
        log_info("\n" + "="*70)
        log_info("STEP 3: CLEANUP")
        log_info("="*70)
        log_info(f"🗑️  Removing temporary files: {taxdump_path}")
        shutil.rmtree(taxdump_path)
        log_info("✅ Cleanup complete")

    # Final summary
    log_info("\n" + "="*70)
    log_info("✅ SETUP COMPLETE")
    log_info("="*70)
    log_info(f"Database: {output_path}")
    log_info(f"Size: {output_path.stat().st_size / (1024**2):.2f} MB")

    if keep_files:
        log_info(f"\nTaxonomy files kept in: {taxdump_path}")

    log_info("\nUsage:")
    log_info(f"  from virome_classifier import TaxonomyDB")
    log_info(f"  tax = TaxonomyDB.from_sqlite('{output_path}')")


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Setup NCBI Taxonomy Database - Complete Automation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-download and build (simplest)
  %(prog)s --output taxonomy.db

  # Specify download directory
  %(prog)s --output taxonomy.db --taxdump-dir ./ncbi_taxonomy

  # Use existing dump files
  %(prog)s --output taxonomy.db --taxdump-dir /path/to/existing/taxdump

  # Download only (don't build database)
  %(prog)s --download-only --taxdump-dir ./taxdump

  # Force re-download and rebuild
  %(prog)s --output taxonomy.db --force

  # Don't keep downloaded files
  %(prog)s --output taxonomy.db --no-keep-files

Complete workflow for virome_classifier:
  1. Setup database:
     python -m virome_classifier.cli.setup_taxonomy_db \\
         --output /path/to/taxonomy.db

  2. Use in Python:
     from virome_classifier import TaxonomyDB
     tax = TaxonomyDB.from_sqlite("/path/to/taxonomy.db")

  3. Or use in CLI:
     python -m virome_classifier.cli.classify \\
         --taxonomy /path/to/taxonomy.db \\
         ...
        """,
    )

    parser.add_argument(
        "--output",
        "-o",
        type=str,
        help="Output SQLite database path (required unless --download-only)",
    )
    parser.add_argument(
        "--taxdump-dir",
        "-d",
        type=str,
        help="Directory for taxdump files (default: temp dir next to output)",
    )
    parser.add_argument(
        "--download-only",
        action="store_true",
        help="Only download taxonomy files, don't build database",
    )
    parser.add_argument(
        "--batch-size",
        "-b",
        type=int,
        default=10000,
        help="Batch size for database inserts (default: 10000)",
    )
    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="Force re-download and overwrite existing files",
    )
    parser.add_argument(
        "--no-keep-files",
        action="store_true",
        help="Remove downloaded files after building database",
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

    # Validate arguments
    if not args.download_only and not args.output:
        print("Error: --output is required unless --download-only is specified")
        return 1

    if args.download_only and not args.taxdump_dir:
        print("Error: --taxdump-dir is required with --download-only")
        return 1

    print("=" * 70)
    print("NCBI TAXONOMY DATABASE SETUP")
    print("=" * 70)

    try:
        setup_taxonomy_database(
            output_db=args.output or "taxonomy.db",
            taxdump_dir=args.taxdump_dir,
            download_only=args.download_only,
            batch_size=args.batch_size,
            force=args.force,
            keep_files=not args.no_keep_files,
            verbose=args.verbose,
        )
        return 0

    except KeyboardInterrupt:
        log_info("\n\n⚠️  Interrupted by user")
        return 130

    except Exception as e:
        log_info(f"\n❌ Fatal error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
