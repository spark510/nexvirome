#!/usr/bin/env python3
"""
Fig4 source-table builder — DNA vs RNA library protocol comparison (asthma cohort).

Question: how does the *library protocol* (DNA vs RNA) change the virome that
NexVirome detects? Three data-backed signals, one table each:

  (1) genome-type composition  — both cohorts are dsDNA(phage)-dominated, but
      only the RNA library captures ssRNA viruses (DNA protocol = DNase/RNase
      then DNA-only extraction, so RNA viruses are structurally absent).
  (2) detection depth/richness — DNA library detects far more (species/sample,
      viral reads/sample) because it reaches the dsDNA phageome deeply.
  (3) RNA-virus capture        — close-up of (1): ssRNA/dsRNA/RT read & species
      fraction, ~0 in DNA vs a real share in RNA.

Input
-----
- NexVirome lca detections (long): paper/figures/Fig3/source_data/realdata_fig3_mixed.csv
    columns: cohort, sample, tool, taxid, name, read_count, group   (tool filtered to "NexVirome lca")
- Cohort sample manifests (give the FULL denominator incl. zero-detection samples):
    paper/figures/Fig3/source_data/{dna,rna}_asthma_groups.csv
- genome-type labels from the taxonomy DB via a validated 3-tier fallback
  (VMR.genome_comp -> ictv_taxonomy.molecule by name -> NCBI lineage realm/phylum/kingdom),
  the same scheme as scripts/benchmark/analyze_molecule_distribution.py (99.4% non-Other coverage here).

IMPORTANT — denominator
-----------------------
RNA has 56 enrolled samples but only 52 have any detection; 4 are zero-detection.
Per-sample richness/depth tables include those zero-detection samples as 0 so the
RNA library's shallowness is not overstated. Cohort sizes: RNA=56, DNA=25.

Output (result_260603/fig4/tables/)
-----------------------------------
- fig4_genometype_per_sample.csv   cohort,sample,group,<molclass...>_reads,<...>_species,total_reads,total_species
- fig4_genometype_cohort.csv       cohort,molecule,reads,species,read_frac,species_frac,n_samples_with_class
- fig4_depth_per_sample.csv        cohort,sample,group,n_species,viral_reads   (incl. zeros)
- fig4_rnavirus_capture.csv        cohort,sample,group,rna_virus_reads,rna_virus_species,rna_read_frac,total_reads
- fig4_taxid_molecule_map.csv      taxid,name,molecule   (label provenance, for audit)

Run:  .venv/bin/python result_260603/fig4/code/fig4_build_tables.py
"""
from __future__ import annotations
import os, sqlite3, sys
import numpy as np
import pandas as pd

NX = "/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/result_260605")
from golden_rule import keep_samples
DB = f"{NX}/resources/db/custom/tax_seq_v20260526_MSL41.db"
# Detection source: REUSE Fig3's NexVirome species table directly. Fig3's
# nexvirome_counts (fig3_realdata_aggregate.py) IS method B on the HQ parquet
# cache with mask_v3_full + breadth>=0.01 — byte-identical to Fig2's
# nexvirome_counts. Verified: fig3_species_long (NexVirome) == an independent
# method-B rerun, 0/1612 rows differ. So Fig2/Fig3/Fig4 share ONE detection
# source; Fig4 only adds a genome-type label. (No separate reclassify needed.)
MIX = f"{NX}/result_260605/fig4/tables/fig4_species_long.csv"
DNA_GROUPS = f"{NX}/paper/figures/Fig3/source_data/dna_asthma_groups.csv"
RNA_GROUPS = f"{NX}/paper/figures/Fig3/source_data/rna_asthma_groups.csv"
OUT = f"{NX}/result_260605/fig4/tables"

# Molecule classes kept, in display order. RT = reverse-transcribing
# (ssRNA-RT retro + dsDNA-RT hepadna); kept separate because it is neither a
# clean DNA nor RNA library target.
MOL_ORDER = ["dsDNA", "ssDNA", "dsRNA", "ssRNA", "RT", "Other"]
# Which classes count as "RNA viruses" for the capture panel (panel 3).
RNA_VIRUS_CLASSES = {"ssRNA", "dsRNA"}


# ----------------------------- genome-type labels -----------------------------
def to_category(m) -> str:
    """Bucket an ICTV genome_comp / molecule string into a Baltimore-ish class."""
    if m is None or m == "" or (isinstance(m, float) and pd.isna(m)):
        return "Other"
    m = str(m)
    if m in ("Viroid", "Unassigned"):
        return "Other"
    if "RT" in m:
        return "RT"
    if "dsRNA" in m:
        return "dsRNA"
    if "ssRNA" in m:
        return "ssRNA"
    if "ssDNA" in m and "dsDNA" not in m:
        return "ssDNA"
    if "dsDNA" in m and "ssDNA" not in m:
        return "dsDNA"
    return "Other"


REALM_TO_CAT = {
    "Duplodnaviria": "dsDNA", "Varidnaviria": "dsDNA", "Adnaviria": "dsDNA",
    "Monodnaviria": "ssDNA", "Ribozyviria": "Other", "Riboviria": None,
}
PHYLUM_TO_CAT = {
    "Pisuviricota": "ssRNA", "Kitrinoviricota": "ssRNA", "Lenarviricota": "ssRNA",
    "Negarnaviricota": "ssRNA", "Duplornaviricota": "dsRNA", "Artverviricota": "RT",
}
KINGDOM_TO_CAT = {"Orthornavirae": "ssRNA", "Pararnavirae": "RT"}


def build_label_fns(db: str):
    """Return three lookups replicating analyze_molecule_distribution.py:
    taxid->cat (VMR), name->cat (ICTV taxonomy), taxid->cat (NCBI lineage trace)."""
    con = sqlite3.connect(db)
    vmr = pd.read_sql(
        "SELECT ncbi_taxid, genome_comp FROM ictv_vmr "
        "WHERE ncbi_taxid IS NOT NULL AND genome_comp IS NOT NULL AND genome_comp!='' "
        "GROUP BY ncbi_taxid", con)
    itx = pd.read_sql(
        "SELECT name, molecule FROM ictv_taxonomy "
        "WHERE rank='species' AND molecule IS NOT NULL GROUP BY name", con)
    nt = pd.read_sql(
        "SELECT taxid, parent_taxid, rank, scientific_name FROM ncbi_taxonomy", con)
    con.close()

    vmr_map = {int(t): to_category(g) for t, g in zip(vmr.ncbi_taxid, vmr.genome_comp)}
    name_map = {n: to_category(m) for n, m in zip(itx.name, itx.molecule)}

    parent = dict(zip(nt.taxid.astype(int), nt.parent_taxid.fillna(0).astype(int)))
    rank = dict(zip(nt.taxid.astype(int), nt["rank"].astype(str)))
    name = dict(zip(nt.taxid.astype(int), nt.scientific_name.astype(str)))

    def lineage_cat(t: int):
        cur = int(t)
        for _ in range(40):
            p = parent.get(cur)
            if not p or p == 0 or p == cur:
                break
            nm, r = name.get(cur, ""), rank.get(cur, "")
            if r == "phylum" and nm in PHYLUM_TO_CAT:
                return PHYLUM_TO_CAT[nm]
            if r == "kingdom" and nm in KINGDOM_TO_CAT:
                return KINGDOM_TO_CAT[nm]
            if r == "realm" and nm in REALM_TO_CAT and REALM_TO_CAT[nm]:
                return REALM_TO_CAT[nm]
            cur = p
        return None

    return vmr_map, name_map, lineage_cat


def label_detections(nv: pd.DataFrame, db: str) -> pd.DataFrame:
    """Add a `molecule` column to the detection table via the 3-tier fallback."""
    vmr_map, name_map, lineage_cat = build_label_fns(db)

    def lab(row):
        t, nm = int(row.taxid), str(row["name"]).strip()
        c = vmr_map.get(t)
        if c is None:
            c = name_map.get(nm)
        if c is None:
            c = lineage_cat(t)
        return c or "Other"

    nv = nv.copy()
    nv["molecule"] = nv.apply(lab, axis=1)
    cov = (nv["molecule"] != "Other").mean()
    print(f"  genome-type label coverage (non-Other): {cov:.1%}")
    return nv


# Genome-length table for TPM (same file Fig3 uses).
GLEN = f"{NX}/paper/figures/Fig5_extra/tables/species_genome_length.csv"


def _add_tpm_weight(nv):
    """Add a per-detection `tpm` weight = genome-length-normalised read count,
    renormalised so each sample's tpm sums to 1 (within that sample's detections).
    This is the SAME TPM Fig2/Fig3 use (reads / (genome_len/1000), then
    sample-normalised), so all Fig4 panels are TPM relative abundance, not raw
    reads. Detections with no genome length get weight 0 (rare; logged)."""
    glen = dict(pd.read_csv(GLEN).itertuples(index=False))
    nv = nv.copy()
    L = nv["taxid"].map(glen).fillna(0).astype(float)
    rpk = np.where(L > 0, nv["read_count"] / (L / 1000.0), 0.0)
    nv["rpk"] = rpk
    samp_tot = nv.groupby(["cohort", "sample"])["rpk"].transform("sum")
    nv["tpm"] = np.where(samp_tot > 0, nv["rpk"] / samp_tot, 0.0)
    missing = int((L <= 0).sum())
    if missing:
        print(f"  [TPM] {missing} detections lack genome length -> weight 0")
    return nv


# ----------------------------- inputs -----------------------------
def load_inputs():
    mix = pd.read_csv(MIX)
    # Fig3 species_long has a `tool` column (NexVirome / Ganon / ...) and a `reads`
    # column. Keep NexVirome only and normalise to read_count.
    if "tool" in mix.columns:
        nv = mix[mix["tool"].isin(["NexVirome", "NexVirome lca"])].copy()
    else:
        nv = mix.copy()
    if "read_count" not in nv.columns and "reads" in nv.columns:
        nv = nv.rename(columns={"reads": "read_count"})
    nv["taxid"] = nv["taxid"].astype(int)
    nv["read_count"] = nv["read_count"].astype(int)
    nv = _add_tpm_weight(nv)          # adds `tpm` (per-sample-normalised) + `rpk`

    dna_g = pd.read_csv(DNA_GROUPS)[["sample", "group"]].assign(cohort="dna")
    rna_g = pd.read_csv(RNA_GROUPS)[["sample", "group"]].assign(cohort="rna")
    manifest = pd.concat([dna_g, rna_g], ignore_index=True)  # FULL denominator
    # result_260605: drop EXCLUDE_SAMPLES (vir17) from the denominator + detections.
    manifest = manifest[manifest["sample"].isin(keep_samples(manifest["sample"]))] \
        .reset_index(drop=True)
    nv = nv[nv["sample"].isin(keep_samples(nv["sample"]))].copy()
    # fig4_species_long already carries a `group` column; drop it so the manifest
    # merge below does not create group_x/group_y.
    nv = nv.drop(columns=[c for c in ("group",) if c in nv.columns])
    # attach group to each detection (for cohort x group panels B'/C')
    nv = nv.merge(manifest[["cohort", "sample", "group"]], on=["cohort", "sample"],
                  how="left")
    return nv, manifest


# ----------------------------- TPM composition helper -----------------------------
def _mean_tpm(nv, label_col, by=("cohort",)):
    """Mean TPM composition over `label_col`, grouped by `by` (e.g. ("cohort",) or
    ("cohort","group")). Per sample the `tpm` weights sum to 1 over its detections;
    we sum tpm per label per sample, average across the n_samples of each `by`
    group, and renormalise so each `by` group sums to 1.

    n_samples denominator = distinct samples WITH detections in that `by` group.
    Returns long df: <by...>, <label_col>, frac."""
    by = list(by)
    g = nv.groupby(by + ["sample", label_col])["tpm"].sum().reset_index()
    nsamp = nv.groupby(by)["sample"].nunique()
    cm = g.groupby(by + [label_col])["tpm"].sum().reset_index()
    cm = cm.merge(nsamp.rename("nsamp").reset_index(), on=by)
    cm["frac"] = cm["tpm"] / cm["nsamp"]
    tot = cm.groupby(by)["frac"].transform("sum")
    cm["frac"] = cm["frac"] / tot
    return cm[by + [label_col, "frac"]]


def _cohort_mean_tpm(nv, label_col):
    """Back-compat wrapper: mean TPM composition by cohort only."""
    return _mean_tpm(nv, label_col, by=("cohort",))


# ----------------------------- table builders -----------------------------
def t_genometype_per_sample(nv, manifest):
    """Wide per-sample: reads & species per molecule class. Zero-detection
    samples (in manifest but absent from nv) appear as all-zero rows."""
    reads = (nv.groupby(["cohort", "sample", "molecule"])["read_count"].sum()
               .unstack(fill_value=0).reindex(columns=MOL_ORDER, fill_value=0))
    species = (nv.groupby(["cohort", "sample", "molecule"])["taxid"].nunique()
                 .unstack(fill_value=0).reindex(columns=MOL_ORDER, fill_value=0))
    reads.columns = [f"{c}_reads" for c in reads.columns]
    species.columns = [f"{c}_species" for c in species.columns]
    wide = reads.join(species).reset_index()

    out = manifest.merge(wide, on=["cohort", "sample"], how="left")
    val_cols = [c for c in out.columns if c.endswith("_reads") or c.endswith("_species")]
    out[val_cols] = out[val_cols].fillna(0).astype(int)
    out["total_reads"] = out[[f"{c}_reads" for c in MOL_ORDER]].sum(axis=1)
    out["total_species"] = out[[f"{c}_species" for c in MOL_ORDER]].sum(axis=1)
    return out.sort_values(["cohort", "sample"]).reset_index(drop=True)


def t_genometype_cohort(nv):
    """Cohort-level genome-type composition as cohort-mean TPM relative abundance
    (consistent with Fig2/Fig3). `tpm_frac` is the headline used by panel A; raw
    read counts/species are kept alongside for reference."""
    tpm = _cohort_mean_tpm(nv, "molecule").rename(columns={"frac": "tpm_frac"})
    rows = []
    for cohort, sub in nv.groupby("cohort"):
        tot_reads = sub["read_count"].sum()
        tot_sp = sub["taxid"].nunique()
        for mol in MOL_ORDER:
            m = sub[sub["molecule"] == mol]
            reads = int(m["read_count"].sum())
            sp = int(m["taxid"].nunique())
            tf = tpm[(tpm["cohort"] == cohort) & (tpm["molecule"] == mol)]["tpm_frac"]
            rows.append(dict(
                cohort=cohort, molecule=mol, reads=reads, species=sp,
                tpm_frac=round(float(tf.iloc[0]), 4) if len(tf) else 0.0,
                read_frac=round(reads / tot_reads, 4) if tot_reads else 0.0,
                species_frac=round(sp / tot_sp, 4) if tot_sp else 0.0,
                n_samples_with_class=int(m["sample"].nunique()),
            ))
    cat = pd.CategoricalDtype(MOL_ORDER, ordered=True)
    df = pd.DataFrame(rows)
    df["molecule"] = df["molecule"].astype(cat)
    return df.sort_values(["cohort", "molecule"]).reset_index(drop=True)


def t_phage_kind_cohort(nv, phage):
    """Panel-A companion: cohort-mean TPM relative abundance of phage vs non-phage
    (DNA, RNA). TPM-based so it is consistent with the genome-type bar."""
    df = nv.copy()
    df["kind"] = np.where(df["taxid"].isin(phage), "phage", "non-phage")
    cm = _cohort_mean_tpm(df, "kind").rename(columns={"frac": "tpm_frac"})
    return cm.sort_values(["cohort", "kind"]).reset_index(drop=True)


def t_depth_per_sample(per_sample):
    """Per-sample detection depth: n_species + viral_reads, zeros included.
    Derived from the per-sample genome-type table so denominators agree."""
    out = per_sample[["cohort", "sample", "group", "total_species", "total_reads"]].copy()
    out = out.rename(columns={"total_species": "n_species", "total_reads": "viral_reads"})
    return out.reset_index(drop=True)


def t_rnavirus_capture(per_sample):
    """Per-sample RNA-virus (ssRNA+dsRNA) capture: reads, species, read fraction.
    Zero-detection samples and DNA samples (≈0) included for honest contrast."""
    rna_read_cols = [f"{c}_reads" for c in RNA_VIRUS_CLASSES]
    rna_sp_cols = [f"{c}_species" for c in RNA_VIRUS_CLASSES]
    out = per_sample[["cohort", "sample", "group"]].copy()
    out["rna_virus_reads"] = per_sample[rna_read_cols].sum(axis=1)
    out["rna_virus_species"] = per_sample[rna_sp_cols].sum(axis=1)
    out["total_reads"] = per_sample["total_reads"]
    out["rna_read_frac"] = np.where(
        out["total_reads"] > 0, out["rna_virus_reads"] / out["total_reads"], 0.0).round(4)
    return out.reset_index(drop=True)


def t_taxid_molecule_map(nv):
    """Audit table: each detected taxid -> molecule label (one row per taxid)."""
    return (nv[["taxid", "name", "molecule"]]
            .drop_duplicates("taxid").sort_values("molecule").reset_index(drop=True))


def _phage_set(db: str):
    """Production phage-taxid set (Caudoviricetes/ssRNA-ssDNA-phage clades + the
    word 'phage' in the reference title), reused from the classifier so Fig4's
    phage/eukaryotic split matches the pipeline's own definition."""
    import sys as _sys
    _sys.path.insert(0, f"{NX}/scripts")
    from virome_classifier.classification.phage_host_rollup import build_phage_host_map
    _tax2host, phage = build_phage_host_map(db)
    return phage


def t_phage_eukaryotic(nv, phage):
    """phage vs eukaryotic-virus split per cohort (reads + species). The headline
    is the SPECIES axis: both cohorts are phage-dominated by reads, but the RNA
    library detects far more eukaryotic-virus SPECIES (DNA reaches the dsDNA
    phageome deeply; RNA opens the eukaryotic virome)."""
    df = nv.copy()
    df["kind"] = np.where(df["taxid"].isin(phage), "phage", "eukaryotic")
    rows = []
    for cohort, sub in df.groupby("cohort"):
        tot_r, tot_s = sub["read_count"].sum(), sub["taxid"].nunique()
        for kind in ("phage", "eukaryotic"):
            k = sub[sub["kind"] == kind]
            reads, sp = int(k["read_count"].sum()), int(k["taxid"].nunique())
            rows.append(dict(
                cohort=cohort, kind=kind, reads=reads, species=sp,
                read_frac=round(reads / tot_r, 4) if tot_r else 0.0,
                species_frac=round(sp / tot_s, 4) if tot_s else 0.0,
            ))
    return pd.DataFrame(rows)


def t_phage_eukaryotic_per_sample(nv, phage):
    """Per-sample eukaryotic-virus species count (for a box/strip overlay): the
    cleanest single number for 'RNA opens the eukaryotic virome'."""
    df = nv.copy()
    df["kind"] = np.where(df["taxid"].isin(phage), "phage", "eukaryotic")
    g = (df.groupby(["cohort", "sample", "kind"])["taxid"].nunique()
           .unstack(fill_value=0).reindex(columns=["phage", "eukaryotic"], fill_value=0)
           .reset_index())
    g.columns = ["cohort", "sample", "phage_species", "eukaryotic_species"]
    return g


def t_taxa_overlap(nv):
    """Detection-set overlap between cohorts: DNA-only / shared / RNA-only, at both
    species and genus level (the complementarity / Venn panel)."""
    rows = []
    dna_sp = set(nv[nv["cohort"] == "dna"]["taxid"])
    rna_sp = set(nv[nv["cohort"] == "rna"]["taxid"])
    rows.append(dict(rank="species", dna_only=len(dna_sp - rna_sp),
                     shared=len(dna_sp & rna_sp), rna_only=len(rna_sp - dna_sp)))
    return pd.DataFrame(rows)


def _hostgenus_map(db: str):
    """taxid -> bacterial/archaeal host genus (same map Fig3 builds:
    refseq_metadata.host + phage_host_from_title; env strings -> '(metagenome)')."""
    import sqlite3
    env = ("metagenome", "sludge", "seawater", "sediment", "environment", "soil",
           "wastewater", "uncultured")
    def hg(h):
        if not h or not h.strip():
            return None
        low = h.lower()
        if any(k in low for k in env):
            return "(metagenome)"
        # phage can't truly have a human host: 'Homo sapiens' is a RefSeq metadata
        # error (taxid 38018 'Bacteriophage sp.') -> drop, don't make 'Homo phage'
        if low.startswith("homo sapiens") or low.startswith("homo "):
            return None
        return h.split()[0]
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
        g = hg(acc2host.get(acc.split(".")[0]))
        if g:
            tax2hg[tx] = g
    con.close()
    return tax2hg


def t_hostgenus_cohort(nv, phage, top_n: int = 10, db: str = DB, by=("cohort",)):
    """Mean PHAGE host-genus composition as TPM relative abundance, grouped by
    `by` (("cohort",) for B; ("cohort","group") for B'). Phage detections only,
    rolled to bacterial/archaeal host genus; tpm re-normalised within phage per
    sample; top_n host genera kept, rest -> 'Other host'.

    THE DNA-vs-RNA contrast: Pseudomonas phage dominates DNA and is absent from RNA;
    Streptococcus/Bacteroides/Fusobacterium phage are shared (oral-gut)."""
    by = list(by)
    tax2hg = _hostgenus_map(db)
    df = nv[nv["taxid"].isin(phage)].copy()
    df["host_genus"] = df["taxid"].map(tax2hg)
    df = df[df["host_genus"].notna()]
    # re-normalise tpm within each sample's host-assignable phage
    samp_tot = df.groupby(["cohort", "sample"])["tpm"].transform("sum")
    df["tpm"] = np.where(samp_tot > 0, df["tpm"] / samp_tot, 0.0)
    cm = _mean_tpm(df, "host_genus", by=by)
    rank = cm.groupby("host_genus")["frac"].max().sort_values(ascending=False)
    keep = list(rank.head(top_n).index)
    cm["host_genus"] = cm["host_genus"].where(cm["host_genus"].isin(keep), "Other host")
    return cm.groupby(by + ["host_genus"])["frac"].sum().reset_index()


def _genus_map(taxids, db: str):
    """taxid -> genus NAME (via the production TaxonomyDB.get_taxid_at_rank +
    get_name). taxa with no genus ancestor keep their own name."""
    import sys as _sys
    _sys.path.insert(0, f"{NX}/scripts")
    from virome_classifier.taxonomy import TaxonomyDB
    tax = TaxonomyDB.from_sqlite(db)
    out = {}
    for t in set(int(x) for x in taxids):
        try:
            g = tax.get_taxid_at_rank(t, "genus")
        except Exception:
            g = None
        out[t] = (tax.get_name(int(g)) if g and g > 0 else None) or tax.get_name(t) or f"taxid_{t}"
    return out


def t_nonphage_cohort(nv, phage, level: str = "species", top_n: int = 10,
                      db: str = DB, by=("cohort",)):
    """Mean NON-PHAGE (eukaryotic-virus) composition as TPM relative abundance,
    at `level` ('species'/'genus'), grouped by `by` (("cohort",) for C;
    ("cohort","group") for C'). tpm re-normalised within each sample's non-phage
    detections; top_n labels kept, rest -> 'Other'.

    Contrast: DNA non-phage = dsDNA eukaryotic virus (herpes EBV/HHV-7, polyoma,
    anellovirus); RNA non-phage adds RNA/RT viruses (PMMoV, Mus provirus)."""
    by = list(by)
    df = nv[~nv["taxid"].isin(phage)].copy()
    if df.empty:
        return pd.DataFrame(columns=by + ["label", "frac"])
    if level == "genus":
        gmap = _genus_map(df["taxid"].unique(), db)
        df["label"] = df["taxid"].map(gmap)
    else:
        df["label"] = df["name"]
    # re-normalise tpm within each sample's non-phage detections, then mean by `by`
    samp_tot = df.groupby(["cohort", "sample"])["tpm"].transform("sum")
    df["tpm"] = np.where(samp_tot > 0, df["tpm"] / samp_tot, 0.0)
    g = _mean_tpm(df, "label", by=by)
    rank = g.groupby("label")["frac"].max().sort_values(ascending=False)
    keep = list(rank.head(top_n).index)
    g["label"] = g["label"].where(g["label"].isin(keep), "Other")
    return g.groupby(by + ["label"])["frac"].sum().reset_index()


# ----------------------------- main -----------------------------
def main():
    os.makedirs(OUT, exist_ok=True)
    print("Loading inputs...")
    nv, manifest = load_inputs()
    print(f"  NexVirome lca detections: {len(nv)} rows, "
          f"samples per cohort = {nv.groupby('cohort')['sample'].nunique().to_dict()}")
    print(f"  manifest (full denominator): {manifest.groupby('cohort')['sample'].nunique().to_dict()}")

    print("Labeling genome types...")
    nv = label_detections(nv, DB)

    print("Loading production phage set...")
    phage = _phage_set(DB)
    print(f"  phage taxids in DB: {len(phage):,}")

    print("Building tables...")
    per_sample = t_genometype_per_sample(nv, manifest)
    cohort = t_genometype_cohort(nv)
    depth = t_depth_per_sample(per_sample)
    rnacap = t_rnavirus_capture(per_sample)
    taxmap = t_taxid_molecule_map(nv)
    phg_cohort = t_phage_eukaryotic(nv, phage)
    phg_persample = t_phage_eukaryotic_per_sample(nv, phage)
    overlap = t_taxa_overlap(nv)
    phage_kind = t_phage_kind_cohort(nv, phage)
    hostgenus = t_hostgenus_cohort(nv, phage, top_n=10)
    nonphage_sp = t_nonphage_cohort(nv, phage, level="species", top_n=10)
    nonphage_gn = t_nonphage_cohort(nv, phage, level="genus", top_n=10)
    # NOTE: the cohort x clinical-group (DNA/RNA x Control/Asthma) expansion of
    # B'/C' lives in its own folder, result_260603/fig4_asthma/ (self-contained
    # build_tables.py + render.py). t_hostgenus_cohort / t_nonphage_cohort accept
    # by=("cohort","group") if you ever want those tables here too.

    tables = {
        "fig4_genometype_per_sample.csv": per_sample,
        "fig4_genometype_cohort.csv": cohort,
        "fig4_depth_per_sample.csv": depth,
        "fig4_rnavirus_capture.csv": rnacap,
        "fig4_taxid_molecule_map.csv": taxmap,
        "fig4_phage_eukaryotic_cohort.csv": phg_cohort,
        "fig4_phage_eukaryotic_per_sample.csv": phg_persample,
        "fig4_taxa_overlap.csv": overlap,
        "fig4_phage_kind_cohort.csv": phage_kind,
        "fig4_hostgenus_cohort.csv": hostgenus,
        "fig4_nonphage_species_cohort.csv": nonphage_sp,
        "fig4_nonphage_genus_cohort.csv": nonphage_gn,
    }
    for fn, df in tables.items():
        df.to_csv(f"{OUT}/{fn}", index=False)
        print(f"  -> {fn}  ({len(df)} rows)")

    # quick sanity print
    print("\n=== cohort genome-type composition (TPM%, read%, sp%) ===")
    show = cohort.copy()
    show["TPM%"] = (show["tpm_frac"] * 100).round(1)
    show["read%"] = (show["read_frac"] * 100).round(1)
    show["sp%"] = (show["species_frac"] * 100).round(1)
    print(show[["cohort", "molecule", "TPM%", "read%", "sp%"]].to_string(index=False))

    print("\n=== phage vs non-phage (cohort-mean TPM %) — panel A companion ===")
    pk = phage_kind.copy()
    pk["TPM%"] = (pk["tpm_frac"] * 100).round(1)
    print(pk[["cohort", "kind", "TPM%"]].to_string(index=False))

    print("\n=== detection depth (mean over FULL cohort incl. zero-detection) ===")
    d = depth.groupby("cohort").agg(
        n_samples=("sample", "nunique"),
        mean_species=("n_species", "mean"),
        median_species=("n_species", "median"),
        mean_reads=("viral_reads", "mean"),
        zero_detection_samples=("n_species", lambda s: int((s == 0).sum())),
    ).round(1)
    print(d.to_string())

    print("\n=== RNA-virus capture (ssRNA+dsRNA) ===")
    r = rnacap.groupby("cohort").agg(
        mean_read_frac=("rna_read_frac", "mean"),
        samples_with_rnavirus=("rna_virus_reads", lambda s: int((s > 0).sum())),
        total_rnavirus_reads=("rna_virus_reads", "sum"),
    ).round(4)
    print(r.to_string())

    print("\n=== phage vs eukaryotic virus (reads + species) ===")
    pe = phg_cohort.copy()
    pe["read%"] = (pe["read_frac"] * 100).round(1)
    pe["sp%"] = (pe["species_frac"] * 100).round(1)
    print(pe[["cohort", "kind", "reads", "species", "read%", "sp%"]].to_string(index=False))
    print("  (headline: eukaryotic-virus SPECIES — RNA detects far more)")

    print("\n=== taxa overlap (complementarity) ===")
    print(overlap.to_string(index=False))

    print("\n=== host-genus phage composition (cohort mean rel-abund %) ===")
    hg = hostgenus.copy()
    hg["pct"] = (hg["frac"] * 100).round(1)
    piv = hg.pivot_table(index="host_genus", columns="cohort", values="pct",
                         fill_value=0.0)
    piv["max"] = piv.max(axis=1)
    print(piv.sort_values("max", ascending=False).drop(columns="max").to_string())

    for lvl, t in (("species", nonphage_sp), ("genus", nonphage_gn)):
        print(f"\n=== non-phage {lvl} composition (cohort mean rel-abund %) ===")
        tt = t.copy()
        tt["pct"] = (tt["frac"] * 100).round(1)
        pv = tt.pivot_table(index="label", columns="cohort", values="pct", fill_value=0.0)
        pv["max"] = pv.max(axis=1)
        print(pv.sort_values("max", ascending=False).drop(columns="max").to_string())


if __name__ == "__main__":
    main()
