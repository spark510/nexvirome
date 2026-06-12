"""
Masking-BED obfuscation / format auto-detection (stdlib-only).

The masking BED encodes which reference coordinates the pipeline blanks out —
this is curated know-how we ship *with* the pipeline but do not want exposed as
plaintext coordinates. This module wraps a BED payload in an obfuscated
container so that opening the file in an editor shows binary noise, while the
pipeline transparently de-obfuscates it at load time (so the curated mask is
always applied, for us and for downstream users alike).

No third-party dependency: uses only gzip + hashlib + a SHA-256 keystream XOR,
so every environment that can run the pipeline can also read the .enc mask.

Honest threat model
-------------------
Because the pipeline MUST apply the mask to run, the de-obfuscation key ships
with the code (the seed below). Anyone who installs the pipeline therefore CAN,
with deliberate effort, intercept the recovered dict and read the coordinates.
This is the same limitation every DRM scheme has: to use the content you must
decode it, and decoding can be observed. What this buys us:

  * opening the .enc file reveals only binary noise (no casual inspection)
  * the coordinates are not greppable / diff-able / copy-pasteable
  * casually "just reusing our BED" requires reverse-engineering, not `cat`

It does NOT make recovery cryptographically impossible (and an XOR keystream is
weaker than AES — but under "key ships with the code" the practical bar is the
same: a motivated reverser wins either way; a casual viewer is stopped by
both). If you ever need true secrecy, the key must be WITHHELD from the
distribution (then the mask only works for the key holder) — a different
product decision than "everyone runs my mask".

Container format
----------------
    magic   : b"NVMSK1\\n"   (7 bytes) — distinguishes .enc from plaintext BED
    payload : SHA256-keystream XOR over gzip(bed_text.encode())

Plaintext BED files (which do NOT start with the magic) are passed through
untouched, so development uses plaintext and distribution uses .enc with the
exact same ``--mask`` argument and the exact same loader.
"""
from __future__ import annotations

import gzip
import hashlib
import os
from pathlib import Path
from typing import Union

# --- container magic -------------------------------------------------------
# Bytes that a plaintext BED can never start with (BED line 0 is an accession
# like "AC_000007.1\t..."). Used to auto-detect encrypted vs plaintext.
MAGIC = b"NVMSK1\n"


# --- key --------------------------------------------------------------------
# The seed is embedded so the shipped pipeline can de-obfuscate without any
# extra file. To rotate for a new distribution, change _SEED and re-run
# `python -m virome_classifier.cli.encrypt_mask`; old .enc files then stop
# decoding (the keystream no longer matches), which is usually what you want.
#
# NEXVIROME_MASK_KEY (env var), if set, overrides the embedded seed — this keeps
# a future "withhold the key" deployment on the same code path with no edit.
_SEED = b"nexvirome::viral-mask::v1::do-not-edit-without-re-encrypting"
_SALT = b"nexvirome-mask-salt-2026"


def _seed() -> bytes:
    env_key = os.environ.get("NEXVIROME_MASK_KEY")
    if env_key:
        return env_key.encode() if isinstance(env_key, str) else env_key
    return _SEED


def _keystream(n: int) -> bytes:
    """Deterministic SHA-256 keystream of length n (counter mode over the seed)."""
    seed = _seed()
    out = bytearray()
    counter = 0
    while len(out) < n:
        block = hashlib.sha256(seed + _SALT + counter.to_bytes(8, "big")).digest()
        out.extend(block)
        counter += 1
    return bytes(out[:n])


def _xor(data: bytes) -> bytes:
    ks = _keystream(len(data))
    return bytes(b ^ k for b, k in zip(data, ks))


# --- public API ------------------------------------------------------------
def is_encrypted(path: Union[str, Path]) -> bool:
    """True if the file begins with the NVMSK1 magic (i.e. an encrypted mask)."""
    try:
        with open(Path(path), "rb") as fh:
            return fh.read(len(MAGIC)) == MAGIC
    except OSError:
        return False


def encrypt_bed_text(bed_text: str) -> bytes:
    """Encode raw BED text into a self-describing container (magic + payload)."""
    return MAGIC + _xor(gzip.compress(bed_text.encode("utf-8")))


def decrypt_to_bed_text(path: Union[str, Path]) -> str:
    """Read an encrypted mask file and return the plaintext BED string.

    Raises ValueError on a missing magic header or a payload that does not
    de-obfuscate to valid gzip (e.g. wrong key / a rotated seed vs. an old file).
    """
    raw = Path(path).read_bytes()
    if not raw.startswith(MAGIC):
        raise ValueError(
            f"{path} is not an encrypted mask (missing {MAGIC!r} header)"
        )
    try:
        return gzip.decompress(_xor(raw[len(MAGIC):])).decode("utf-8")
    except (OSError, EOFError, UnicodeDecodeError) as e:
        # gzip.decompress raises BadGzipFile (subclass of OSError) / EOFError on a
        # wrong keystream; surface it as the same clean error as a bad key.
        raise ValueError(
            "Failed to decrypt mask: wrong key (NEXVIROME_MASK_KEY mismatch or "
            "the embedded seed was rotated after this file was encrypted)."
        ) from e
