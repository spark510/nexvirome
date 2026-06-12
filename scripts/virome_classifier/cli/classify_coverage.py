#!/usr/bin/env python3
"""
DEPRECATED entry point — kept for backward compatibility.

Coverage-based classification is now handled by the unified CLI:

    python -m virome_classifier.cli.classify --mode coverage ...

This shim translates the legacy `classify_coverage` flags to the unified CLI
and forwards the call, so existing scripts keep working. The previous
implementation imported modules (`virome_classifier.masking`,
`virome_classifier.kraken_output`, `coverage_based_classifier` v1) that no
longer exist; the canonical coverage engine is CoverageBasedClassifier3,
reached via `classify.py --mode coverage`.
"""
from __future__ import annotations

import argparse
import sys
from typing import List

from .classify import main as classify_main


# Legacy flags that the unified CLI does not (yet) accept. They are dropped
# with a warning rather than causing a hard failure. CoverageBasedClassifier3
# reads segment info from the taxonomy DB and uses CoverageThresholds defaults.
_DROPPED_FLAGS_WITH_VALUE = {
    "--min-segment-consistency",
    "--max-lca-rank",
    "--segment-info",
}
_DROPPED_FLAGS_BOOLEAN = {
    "--require-unique-reads",
    "--fractional-assignment",
    "--no-kraken-output",
}
# Legacy flag name -> unified CLI flag name
_RENAMED_FLAGS = {
    "--min-coverage": "--min-unmasked-coverage",
}


def _translate_argv(argv: List[str]) -> List[str]:
    """Map legacy classify_coverage flags onto the unified classify CLI."""
    out: List[str] = ["--mode", "coverage"]
    i = 0
    dropped = []
    while i < len(argv):
        tok = argv[i]
        if tok in _DROPPED_FLAGS_BOOLEAN:
            dropped.append(tok)
            i += 1
            continue
        if tok in _DROPPED_FLAGS_WITH_VALUE:
            dropped.append(tok)
            i += 2  # skip flag + its value
            continue
        if tok in _RENAMED_FLAGS:
            out.append(_RENAMED_FLAGS[tok])
            i += 1
            continue
        out.append(tok)
        i += 1

    if dropped:
        print(
            f"[classify_coverage] WARNING: ignoring legacy flags not supported "
            f"by the unified CLI: {', '.join(sorted(set(dropped)))}",
            file=sys.stderr,
        )
    return out


def main() -> int:
    print(
        "[classify_coverage] DEPRECATED: forwarding to "
        "`virome_classifier.cli.classify --mode coverage`.",
        file=sys.stderr,
    )
    sys.argv = [sys.argv[0]] + _translate_argv(sys.argv[1:])
    return classify_main()


if __name__ == "__main__":
    sys.exit(main())
