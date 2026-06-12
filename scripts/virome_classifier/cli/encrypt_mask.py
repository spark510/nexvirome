#!/usr/bin/env python3
"""
Encrypt a plaintext masking BED into the distributable .enc container.

The pipeline loads the result transparently (MaskLoader auto-detects the magic
header and decrypts in memory), so distribute the .enc and point ``--mask`` at
it exactly as you would a plaintext BED.

Usage:
    python -m virome_classifier.cli.encrypt_mask \\
        --in  resources/.../viral_mask_20260525_combined.bed \\
        --out resources/.../viral_mask_20260525_combined.bed.enc

    # decrypt back (for your own audit; needs the same key)
    python -m virome_classifier.cli.encrypt_mask --decrypt \\
        --in  mask.bed.enc --out mask.roundtrip.bed

Key: derived from the code-embedded seed in mask_crypt, or overridden by the
NEXVIROME_MASK_KEY environment variable. Rotating the seed invalidates old .enc
files — re-run this tool after a rotation.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from ..alignment.filters import mask_crypt
from ..alignment.filters.mask_loader import MaskLoader


def _encrypt(in_path: Path, out_path: Path) -> int:
    bed_text = in_path.read_text(encoding="utf-8")
    blob = mask_crypt.encrypt_bed_text(bed_text)
    out_path.write_bytes(blob)

    # Round-trip self-check: decrypt + parse the just-written file and compare
    # target/interval counts against the source, so a silently corrupt or
    # key-mismatched output is caught here instead of in a pipeline run.
    src_dict = MaskLoader.from_bed_file(str(in_path))
    enc_dict = MaskLoader.from_bed_file(str(out_path))
    src_targets, enc_targets = len(src_dict), len(enc_dict)
    src_iv = sum(r.n_regions for r in src_dict.values())
    enc_iv = sum(r.n_regions for r in enc_dict.values())
    if (src_targets, src_iv) != (enc_targets, enc_iv):
        print(
            f"❌ Round-trip mismatch! source=({src_targets} targets, {src_iv} "
            f"intervals) vs encrypted=({enc_targets}, {enc_iv})",
            file=sys.stderr,
        )
        return 1

    print(
        f"🔒 Encrypted {in_path} -> {out_path}\n"
        f"   {src_targets} targets, {src_iv} intervals; round-trip verified.\n"
        f"   ({in_path.stat().st_size:,} B plaintext -> "
        f"{out_path.stat().st_size:,} B encrypted)"
    )
    return 0


def _decrypt(in_path: Path, out_path: Path) -> int:
    if not mask_crypt.is_encrypted(in_path):
        print(f"❌ {in_path} is not an encrypted mask file.", file=sys.stderr)
        return 1
    out_path.write_text(mask_crypt.decrypt_to_bed_text(in_path), encoding="utf-8")
    print(f"🔓 Decrypted {in_path} -> {out_path}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in", dest="in_path", required=True, help="Input file")
    ap.add_argument("--out", dest="out_path", required=True, help="Output file")
    ap.add_argument("--decrypt", action="store_true", help="Decrypt instead of encrypt")
    args = ap.parse_args()

    in_path, out_path = Path(args.in_path), Path(args.out_path)
    if not in_path.exists():
        print(f"❌ Input not found: {in_path}", file=sys.stderr)
        return 1
    try:
        return _decrypt(in_path, out_path) if args.decrypt else _encrypt(in_path, out_path)
    except Exception as e:  # noqa: BLE001 - surface a clean message to the CLI user
        print(f"❌ {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
