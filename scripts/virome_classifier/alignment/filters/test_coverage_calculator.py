#!/usr/bin/env python3
"""Regression tests for CoverageCalculator.calculate_merged_coverage.

The merge loop (covered-bp ratio of overlapping intervals) was vectorised with a
difference array for speed. These tests pin it byte-identical to the previous
Python merge loop, including the edge case where interval ends exceed the target
length (the old code did NOT clamp; it summed merged lengths and capped the final
ratio at 1.0). Run:
    python scripts/virome_classifier/alignment/filters/test_coverage_calculator.py
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
from virome_classifier.alignment.filters.coverage_calculator import CoverageCalculator


def _old_merge(intervals, L):
    """Reference implementation (the pre-vectorisation Python loop)."""
    if len(intervals) == 0:
        return 0.0
    si = intervals[intervals[:, 0].argsort()]
    merged = []
    cs, ce = si[0]
    for s, e in si[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            merged.append([cs, ce]); cs, ce = s, e
    merged.append([cs, ce])
    tot = sum(e - s for s, e in merged)
    return min(tot / L, 1.0)


def test_empty():
    cc = CoverageCalculator()
    assert cc.calculate_merged_coverage(np.empty((0, 2), int), 1000) == 0.0


def test_matches_old_merge_random():
    cc = CoverageCalculator()
    rng = np.random.default_rng(0)
    for _ in range(500):
        L = int(rng.integers(100, 5000))
        n = int(rng.integers(1, 50))
        s = rng.integers(0, L + 500, n)          # deliberately allow ends > L
        e = s + rng.integers(1, 300, n)
        iv = np.column_stack([s, e])
        a = _old_merge(iv.copy(), L)
        b = cc.calculate_merged_coverage(iv.copy(), L)
        assert abs(a - b) < 1e-9, (iv.tolist(), L, a, b)


def test_full_and_capped():
    cc = CoverageCalculator()
    # full single interval
    assert cc.calculate_merged_coverage(np.array([[0, 1000]]), 1000) == 1.0
    # overlapping intervals collapse to one
    assert abs(cc.calculate_merged_coverage(np.array([[0, 600], [400, 800]]), 1000)
               - 0.8) < 1e-9
    # ends beyond target are capped at ratio 1.0
    assert cc.calculate_merged_coverage(np.array([[0, 2000]]), 1000) == 1.0


# ───────────────────── masked-hit classification (Case 4-2 boundary) ──────────
from virome_classifier.alignment.filters.mask_loader import MaskedRegion
import pandas as pd


def _region(intervals):
    """MaskedRegion from list of (start,end)."""
    arr = np.array(intervals)
    return MaskedRegion("X", arr[:, 0], arr[:, 1])


def test_identify_masked_boundary_50_is_masked():
    """Case 4-2: a hit with EXACTLY 50% overlap is masked (>= threshold)."""
    cc = CoverageCalculator(overlap_threshold=0.5)
    region = _region([(400, 600)])
    # hit [300,500): len 200, overlap with [400,600) = 100 -> ratio 0.5 -> masked
    starts = np.array([300]); ends = np.array([500])
    out = cc.identify_masked_hits(starts, ends, region)
    assert out.tolist() == [True], out


def test_identify_masked_just_below_50_is_unmasked():
    """Case 1-3 / 4-2: < 50% overlap stays unmasked (< threshold)."""
    cc = CoverageCalculator(overlap_threshold=0.5)
    region = _region([(400, 600)])
    # hit [298,500): len 202, overlap [400,500)∩[400,600) = 100 -> 100/202=0.495 -> unmasked
    starts = np.array([298]); ends = np.array([500])
    out = cc.identify_masked_hits(starts, ends, region)
    assert out.tolist() == [False], out


def test_identify_masked_no_overlap_and_full():
    cc = CoverageCalculator(overlap_threshold=0.5)
    region = _region([(400, 600)])
    # no overlap, and fully inside
    starts = np.array([0,   450]); ends = np.array([200, 550])
    out = cc.identify_masked_hits(starts, ends, region)
    assert out.tolist() == [False, True], out


def test_per_hit_masked_labels_batch_and_index():
    """per_hit_masked_labels: all targets in one pass, robust to non-reset index,
    targets without a mask stay unmasked."""
    cc = CoverageCalculator(overlap_threshold=0.5)
    masks = {"A": _region([(100, 300)])}  # B has no mask
    df = pd.DataFrame({
        "target": ["A", "A", "B", "A"],
        "tstart": [120, 0,   50,  290],
        "tend":   [280, 80,  250, 400],
    }, index=[7, 8, 9, 10])  # deliberately non-contiguous index
    lbl = cc.per_hit_masked_labels(df, masks)
    # A[120,280): 160/160 masked; A[0,80): no overlap; B: no mask;
    # A[290,400): overlap [290,300)=10 / len 110 -> unmasked
    assert lbl.loc[7] == True
    assert lbl.loc[8] == False
    assert lbl.loc[9] == False
    assert lbl.loc[10] == False
    assert list(lbl.index) == [7, 8, 9, 10]


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(f"PASS {name}")
    print("all coverage_calculator tests passed")
