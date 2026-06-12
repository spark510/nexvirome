#!/usr/bin/env python3
"""
Shared helpers for the Fig4 per-sample "4-view" family (test/).

Data source: result_260603/fig4/tables/fig4_species_long.csv  (ALL 81 samples,
asthma + control; NexVirome method-B, built by fig4_build_species_long.py).

The composition is a SINGLE integrated stack per sample:
    phage species  -> rolled up to BACTERIAL/ARCHAEAL HOST GENUS ("<Genus> phage")
    non-phage virus-> kept at its own SPECIES or GENUS (per view)

Four views are produced by 4 thin driver scripts, each writing into its own
test/ subfolder:
    view_species_relabund/   species-level  x  relative abundance
    view_genus_relabund/     genus-level    x  relative abundance
    view_species_reads/      species-level  x  read count
    view_genus_reads/        genus-level    x  read count

Within every plot the four clinical groups are separated by a gap + a label:
    DNA·Control | DNA·Asthma   ||   RNA·Control | RNA·Asthma
"""
from __future__ import annotations
import os, sys, sqlite3
import numpy as np, pandas as pd

NX = "/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/scripts")
sys.path.insert(0, f"{NX}/result_260605")
from golden_rule import DB, GENOME_LENGTH_CSV, DNA_GROUPS, RNA_GROUPS, keep_samples
from virome_classifier.taxonomy import TaxonomyDB

SPECIES_LONG = f"{NX}/result_260605/fig4/tables/fig4_species_long.csv"

GROUPS = [("dna", "Control", "DNA · Ctrl"),
          ("dna", "Asthma",  "DNA · Asthma"),
          ("rna", "Control", "RNA · Ctrl"),
          ("rna", "Asthma",  "RNA · Asthma")]

_TAX = None
def _tax():
    global _TAX
    if _TAX is None:
        _TAX = TaxonomyDB.from_sqlite(DB)
    return _TAX


# ---------------------------------------------------------------- taxonomy maps
_ENV = ("metagenome", "sludge", "seawater", "sediment", "environment", "soil",
        "wastewater", "uncultured")
def _host_genus_name(h):
    if not h or not h.strip():
        return None
    low = h.lower()
    if any(k in low for k in _ENV):
        return "(metagenome)"
    # a phage can't truly have a human host: 'Homo sapiens' host metadata is an
    # error (taxid 38018 'Bacteriophage sp.'), so drop it -> 'Other-host phage'
    if low.startswith("homo sapiens") or low.startswith("homo "):
        return None
    return h.split()[0]


def build_phage_and_host(db=DB):
    """Return (phage_taxids:set, tax2hostgenus:dict). Phage set from the pipeline's
    own definition; host genus from refseq_metadata.host (+ title fallback)."""
    sys.path.insert(0, f"{NX}/scripts")
    from virome_classifier.classification.phage_host_rollup import build_phage_host_map
    _t, phage = build_phage_host_map(db)

    con = sqlite3.connect(db)
    acc2host = {a: h.strip() for a, h in con.execute(
        "SELECT accession, host FROM refseq_metadata WHERE host IS NOT NULL AND trim(host)!=''")}
    try:
        for a, h in con.execute("SELECT base_accession, host FROM phage_host_from_title"):
            acc2host.setdefault(a, h)
    except sqlite3.OperationalError:
        pass
    tax2hg = {}
    for acc, tx in con.execute("SELECT accession, taxid FROM refseq_sequences"):
        if tx in tax2hg:
            continue
        h = acc2host.get(acc.split(".")[0])
        if h:
            hg = _host_genus_name(h)
            if hg:
                tax2hg[tx] = hg
    con.close()
    return set(int(t) for t in phage), tax2hg


_GENUS = {}
def genus_name_of(taxid):
    """Genus NAME for a taxid; falls back to the species name if no genus rank."""
    t = int(taxid)
    if t not in _GENUS:
        tax = _tax()
        try:
            g = tax.get_taxid_at_rank(t, "genus")
        except Exception:
            g = None
        if g and int(g) > 0:
            _GENUS[t] = tax.get_name(int(g)) or f"taxid_{int(g)}"
        else:
            _GENUS[t] = tax.get_name(t) or f"taxid_{t}"
    return _GENUS[t]


# ---------------------------------------------------------------- the data
def load_labeled(level):
    """Load fig4_species_long and add:
       - `label`  : integrated taxon label for this `level` ('species' | 'genus'),
                    with phage -> '<HostGenus> phage' and non-phage -> species/genus
       - `is_phage`
    Returns the per-detection DataFrame (reads, rel_abund kept)."""
    df = pd.read_csv(SPECIES_LONG)
    df["taxid"] = df["taxid"].astype(int)
    phage, t2hg = build_phage_and_host()
    df["is_phage"] = df["taxid"].isin(phage)

    def lab(row):
        if row["is_phage"]:
            hg = t2hg.get(int(row["taxid"]))
            return (f"{hg} phage" if hg else "Other-host phage")
        # non-phage: species name as-is, or genus rollup
        if level == "genus":
            return genus_name_of(row["taxid"])
        return str(row["name"])

    df["label"] = df.apply(lab, axis=1)
    return df


def manifest():
    dna = pd.read_csv(DNA_GROUPS)[["sample", "group"]].assign(cohort="dna")
    rna = pd.read_csv(RNA_GROUPS)[["sample", "group"]].assign(cohort="rna")
    man = pd.concat([dna, rna], ignore_index=True)
    # result_260605: drop EXCLUDE_SAMPLES (vir17) so sample_layout never lays out a
    # vir17 column even if a stray row leaked into the species table.
    return man[man["sample"].isin(keep_samples(man["sample"]))].reset_index(drop=True)


# ---------------------------------------------------------------- aggregation
def build_matrix(df, metric, top_n=16, other="Other", subset=None):
    """Per-sample value matrix over `label`.
       metric='relabund' -> each sample column sums to 1 (renormalised over the
                            detections present in `d`, i.e. within the subset).
       metric='reads'     -> raw read counts.
       subset             -> None=all; 'phage'/'non-phage' restricts the rows.
    Returns (piv[label x sample], label_order incl. 'Other')."""
    d = df.copy()
    if subset == "phage":
        d = d[d["is_phage"]].copy()
    elif subset == "non-phage":
        d = d[~d["is_phage"]].copy()
    if metric == "relabund":
        # renormalise rel_abund within each sample over the full integrated set
        st = d.groupby("sample")["rel_abund"].transform("sum")
        d["val"] = np.where(st > 0, d["rel_abund"] / st, 0.0)
    else:
        d["val"] = d["reads"].astype(float)

    g = d.groupby(["sample", "label"])["val"].sum().reset_index()
    # global ranking to choose the top_n labels (by summed value), rest -> Other
    rank = g.groupby("label")["val"].sum().sort_values(ascending=False)
    keep = list(rank.head(top_n).index)
    g["label"] = g["label"].where(g["label"].isin(keep), other)
    order = [l for l in keep] + ([other] if (g["label"] == other).any() else [])
    piv = g.groupby(["label", "sample"])["val"].sum().unstack(fill_value=0.0)
    return piv, order


def phage_split_matrix(df, metric="reads"):
    """Per-sample phage vs non-phage composition.
       metric='reads' -> read fraction (each sample sums to 1).
    Returns (piv[{'phage','non-phage'} x sample], ['phage','non-phage'])."""
    d = df.copy()
    d["kind"] = np.where(d["is_phage"], "phage", "non-phage")
    val = d["reads"].astype(float) if metric == "reads" else d["rel_abund"].astype(float)
    d = d.assign(v=val)
    st = d.groupby("sample")["v"].transform("sum")
    d["v"] = np.where(st > 0, d["v"] / st, 0.0)
    piv = d.groupby(["kind", "sample"])["v"].sum().unstack(fill_value=0.0)
    return piv, ["phage", "non-phage"]


PHAGE_KIND_COLOR = {"phage": "#4C9F70", "non-phage": "#E6B800"}


def sample_layout(df, manifest_df):
    """x positions per sample, grouped with gaps; within a group richest-first by
    total read count (so the populated samples lead each block)."""
    tot = df.groupby("sample")["reads"].sum()
    xs, spans = [], []
    x = 0.0
    for cohort, grp, glab in GROUPS:
        sams = manifest_df[(manifest_df.cohort == cohort)
                           & (manifest_df.group == grp)]["sample"].tolist()
        sams = sorted(sams, key=lambda s: -float(tot.get(s, 0.0)))
        start = x
        for s in sams:
            xs.append((s, x)); x += 1.0
        spans.append((glab, start, x - 1.0, len(sams)))
        x += 1.8   # blank gap before next group
    return xs, spans


# ---------------------------------------------------------------- palette
# Panel-A-toned palette: built from the bright Baltimore-class hues
# (cyan/coral/mint/tan family of CLASS_COLOR) and hue-rotated so the MOST
# abundant taxon is not a heavy navy — fixes the previous all-blue look.
PALETTE = ["#00B2CA", "#FF595E", "#7DCFB6", "#FBD1A2", "#8AC926", "#FF924C",
           "#1982C4", "#F15BB5", "#6A4C93", "#52A675", "#E07A5F", "#118AB2",
           "#FFCA3A", "#C44E52", "#3CAEA3", "#9B5DE5", "#00BBF9", "#BC8A5F"]
# D (non-phage) uses a DIFFERENT palette from C (phage host) so the two
# composition panels read as distinct — warm/earthy + plum/teal register,
# led by a deep coral rather than C's cyan.
D_PALETTE = ["#E07A5F", "#3D8C8C", "#BC9C22", "#8E6C9B", "#6A994E", "#C44E52",
             "#4D7C8A", "#D08C60", "#A24936", "#7A9E7E", "#9C6644", "#5C5470",
             "#B5838D", "#3CAEA3", "#CB997E", "#6B705C"]
OTHER_COLOR = "#B0B0B0"


# --------- short display names for legends (keep data keys unchanged) ---------
_ABBR = {
    "Human lung-associated vientovirus FB": "Vientovirus FB",
    "Human betaherpesvirus 7": "HHV-7",
    "Human betaherpesvirus 6B": "HHV-6B",
    "human gammaherpesvirus 4": "HHV-4 (EBV)",
    "Human herpesvirus 4 type 2": "HHV-4 type 2",
    "Human alphaherpesvirus 1": "HHV-1 (HSV-1)",
    "Adelie penguin polyomavirus": "Adelie PyV",
    "Severe acute respiratory syndrome coronavirus 2": "SARS-CoV-2",
    "Tobacco mild green mosaic virus": "TMGMV",
    "Pepper mild mottle virus": "PMMoV",
    "Moloney murine leukemia virus": "MoMLV",
    "TTV-like mini virus": "TTV-like mini",
    "Torque teno virus 29": "TTV-29",
    "Torque teno virus 15": "TTV-15",
}


def short_name(label):
    """Shorten a long virus name for the legend. '<Genus> phage' kept (already
    short); known long species mapped to a compact alias; otherwise unchanged."""
    if label in _ABBR:
        return _ABBR[label]
    # generic fallback: 'Human <something>virus N' -> keep tail, drop 'Human'
    if label.startswith("Human ") and len(label) > 22:
        return label[len("Human "):]
    return label


def colour_map(order, other="Other", palette=None):
    palette = palette if palette is not None else PALETTE
    cols, i = {}, 0
    for lab in order:
        if lab == other:
            cols[lab] = OTHER_COLOR
        else:
            cols[lab] = palette[i % len(palette)]; i += 1
    return cols
