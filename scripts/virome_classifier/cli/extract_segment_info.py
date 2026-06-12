#!/usr/bin/env python3
"""
Extract segment information from RefSeq viral genomes.

This script parses a RefSeq viral genome FASTA file and extracts segment
information for multi-segmented viruses (e.g., influenza, bunyaviruses).

Usage:
    python -m virome_classifier.cli.extract_segment_info \
        --fasta viral.genomic.fna \
        --output segment_info.csv \
        --taxonomy taxonomy.db

Output CSV format:
    accession,species_taxid,virus_name,segment_name,genome_size,segment_size
"""

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Optional, Tuple
import csv


def parse_fasta_headers(fasta_file: Path) -> Dict[str, List[Dict]]:
    """
    Parse FASTA file and extract segment information.

    Args:
        fasta_file: Path to FASTA file

    Returns:
        Dictionary mapping virus name to list of segment info
    """
    segment_data = defaultdict(list)
    current_accession = None
    current_seq_length = 0

    print(f"📖 Parsing FASTA file: {fasta_file}")

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.rstrip()

            if line.startswith('>'):
                # Save previous sequence length
                if current_accession and current_seq_length > 0:
                    # Find and update the entry
                    for virus_name, segments in segment_data.items():
                        for seg in segments:
                            if seg['accession'] == current_accession:
                                seg['segment_size'] = current_seq_length
                                break

                # Parse new header
                header = line[1:].strip()
                accession = header.split()[0]
                current_accession = accession
                current_seq_length = 0

                # Check for segment information
                # Updated pattern to handle various formats:
                # - "segment 2, complete sequence"
                # - "segment 2 polymerase PB1 (PB1) gene, complete cds"
                # - "segment 7 nonstructural protein 2 (NS2) and nonstructural protein 1 (NS1), genes"
                segment_match = re.search(
                    r'segment\s+([A-Za-z0-9\-]+)',
                    header,
                    re.IGNORECASE
                )

                if segment_match:
                    segment_name = segment_match.group(1).strip()

                    # Extract virus name (before "segment")
                    virus_match = re.search(r'^(\w+\.\d+)\s+(.+?)\s+(?:isolate\s+.*?\s+)?segment', header)
                    if virus_match:
                        virus_name = virus_match.group(2).strip()
                    else:
                        # Alternative pattern
                        virus_match = re.search(r'^(\w+\.\d+)\s+(.+?)\s+segment', header)
                        if virus_match:
                            virus_name = virus_match.group(2).strip()
                        else:
                            continue

                    segment_data[virus_name].append({
                        'accession': accession,
                        'virus_name': virus_name,
                        'segment_name': segment_name,
                        'segment_size': 0,  # Will be filled when we see next header
                        'full_header': header
                    })
            else:
                # Count sequence length
                current_seq_length += len(line)

        # Handle last sequence
        if current_accession and current_seq_length > 0:
            for virus_name, segments in segment_data.items():
                for seg in segments:
                    if seg['accession'] == current_accession:
                        seg['segment_size'] = current_seq_length
                        break

    print(f"✅ Found {len(segment_data)} viruses with segment information")

    # Filter to only multi-segment viruses
    multi_segment = {k: v for k, v in segment_data.items() if len(v) > 1}
    print(f"✅ Found {len(multi_segment)} multi-segmented viruses")

    return multi_segment


def get_taxid_from_accession(accession: str, taxonomy_db: Optional[Path]) -> Optional[int]:
    """
    Get taxonomic ID from accession using taxonomy database.

    Args:
        accession: RefSeq accession
        taxonomy_db: Path to taxonomy database

    Returns:
        Taxonomic ID or None
    """
    # TODO: Implement actual lookup from taxonomy database or NCBI API
    # For now, return None - user can fill this in later
    return None


def write_segment_csv(
    segment_data: Dict[str, List[Dict]],
    output_file: Path,
    taxonomy_db: Optional[Path] = None,
) -> None:
    """
    Write segment information to CSV file.

    Args:
        segment_data: Dictionary of virus name to segment info
        output_file: Output CSV file
        taxonomy_db: Optional taxonomy database for taxid lookup
    """
    print(f"💾 Writing segment information to: {output_file}")

    total_segments = 0

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)

        # Header
        writer.writerow([
            'accession',
            'species_taxid',
            'virus_name',
            'segment_name',
            'genome_size',
            'segment_size',
        ])

        # Data
        for virus_name, segments in sorted(segment_data.items()):
            # Calculate total genome size
            genome_size = sum(seg['segment_size'] for seg in segments)

            # Get taxid (if possible)
            taxid = None
            if taxonomy_db and segments:
                taxid = get_taxid_from_accession(segments[0]['accession'], taxonomy_db)

            # Write each segment
            for seg in sorted(segments, key=lambda x: x['segment_name']):
                writer.writerow([
                    seg['accession'],
                    taxid if taxid else '',  # Empty if unknown
                    virus_name,
                    seg['segment_name'],
                    genome_size,
                    seg['segment_size'],
                ])
                total_segments += 1

    print(f"✅ Wrote {total_segments} segments from {len(segment_data)} viruses")

    # Print some statistics
    print("\n" + "=" * 70)
    print("STATISTICS")
    print("=" * 70)

    segment_counts = defaultdict(int)
    for virus_name, segments in segment_data.items():
        segment_counts[len(segments)] += 1

    print("Viruses by segment count:")
    for count in sorted(segment_counts.keys()):
        num_viruses = segment_counts[count]
        print(f"  {count} segments: {num_viruses} viruses")

    # Top viruses by segment count
    top_viruses = sorted(segment_data.items(), key=lambda x: len(x[1]), reverse=True)[:10]
    print("\nTop 10 viruses by segment count:")
    for i, (virus_name, segments) in enumerate(top_viruses, 1):
        print(f"  {i}. {virus_name}: {len(segments)} segments")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Extract segment information from RefSeq viral genomes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract from RefSeq viral genome file
  %(prog)s --fasta viral.genomic.fna --output segment_info.csv

  # With taxonomy database for taxid lookup
  %(prog)s --fasta viral.genomic.fna \\
           --output segment_info.csv \\
           --taxonomy taxonomy.db

Output format (CSV):
  accession,species_taxid,virus_name,segment_name,genome_size,segment_size
  NC_001477.1,11320,Dengue virus 1,,,10735
  NC_002640.1,11320,Dengue virus 2,,,10723
  ...

The segment_info.csv can be used with coverage-based classifier:
  python -m virome_classifier.cli.classify_coverage \\
      --segment-info segment_info.csv \\
      ...
        """,
    )

    parser.add_argument(
        '--fasta', '-f',
        type=str,
        required=True,
        help='Input FASTA file (RefSeq viral genomes)',
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        required=True,
        help='Output CSV file',
    )
    parser.add_argument(
        '--taxonomy', '-t',
        type=str,
        help='Taxonomy database for taxid lookup (optional)',
    )
    parser.add_argument(
        '--min-segments',
        type=int,
        default=2,
        help='Minimum number of segments (default: 2)',
    )

    args = parser.parse_args()

    # Validate input
    fasta_file = Path(args.fasta)
    if not fasta_file.exists():
        print(f"❌ Error: FASTA file not found: {fasta_file}")
        return 1

    output_file = Path(args.output)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    taxonomy_db = Path(args.taxonomy) if args.taxonomy else None
    if taxonomy_db and not taxonomy_db.exists():
        print(f"⚠️  Warning: Taxonomy database not found: {taxonomy_db}")
        print("   Continuing without taxid lookup...")
        taxonomy_db = None

    print("=" * 70)
    print("SEGMENT INFORMATION EXTRACTOR")
    print("=" * 70)
    print(f"Input FASTA: {fasta_file}")
    print(f"Output CSV: {output_file}")
    print(f"Min segments: {args.min_segments}")
    print("=" * 70 + "\n")

    try:
        # Parse FASTA
        segment_data = parse_fasta_headers(fasta_file)

        # Filter by minimum segments
        if args.min_segments > 1:
            before = len(segment_data)
            segment_data = {
                k: v for k, v in segment_data.items()
                if len(v) >= args.min_segments
            }
            after = len(segment_data)
            if before != after:
                print(f"🔍 Filtered to {after} viruses with ≥{args.min_segments} segments")

        if not segment_data:
            print("⚠️  No segmented viruses found")
            return 1

        # Write output
        write_segment_csv(segment_data, output_file, taxonomy_db)

        print("\n" + "=" * 70)
        print("✅ EXTRACTION COMPLETE")
        print("=" * 70)
        print(f"Output: {output_file}")
        print("\nNext steps:")
        print("1. Review the CSV file")
        print("2. Add species_taxid if needed")
        print("3. Use with coverage-based classifier:")
        print(f"   python -m virome_classifier.cli.classify_coverage \\")
        print(f"       --segment-info {output_file} \\")
        print(f"       ...")

        return 0

    except KeyboardInterrupt:
        print("\n\n⚠️  Interrupted by user")
        return 130

    except Exception as e:
        print(f"\n❌ Fatal error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
