#!/usr/bin/env python3
"""
Fig3 real-data aggregation (FINALIZED 2026-06-03).

Real cohorts: RNA 56 (A*) + DNA 25 (vir*). NexVirome vs Ganon / Kraken2 /
Metabuli / Phanta. **External tools NATIVE DB only** (real data has no ground
truth; native = real-world use). NexVirome's viral-only DB is its native.

Per (cohort, tool, sample) we collect per-species {taxid: reads}, then build:
  - SPECIES level  : count + relative abundance
  - GENUS level    : species rolled to genus (taxonomy) — count + relative abundance
  - HOST-GENUS     : phage species rolled to their bacterial/archaeal host genus
                     — **NexVirome only** (per request), count + relative abundance

rel_abund definition (column is the same name for every tool, but tool-specific):
  - NexVirome    : TPM-style, genome-length normalized (its species taxids are 100%
                   covered by the viral genome_length table). reads/(genome_len/1000),
                   normalized to sum=1 per sample; genus/host roll up species weights.
  - external     : the tool's own read fraction among its viral calls (their NCBI
                   taxids are mostly absent from the viral genome_length table, so
                   TPM is undefined — read fraction is what they natively report).

Counts / abundance sources:
  NexVirome : HQ parquet (best-hit + unmasked breadth >= 0.01, mask_v3_full); reads
  Kraken2   : Bracken new_est_reads (result_260605/fig3/bracken)
  Metabuli  : Bracken new_est_reads
  Ganon     : .tre species cumulative reads (own estimator)
  Phanta    : bracken_species.filtered (own bracken)
All restricted to viral taxa (under 10239) for the species/genus tables.

Outputs -> result_260605/fig3/tables/:
  fig3_species_long.csv   cohort,tool,sample,taxid,name,reads,rel_abund
  fig3_genus_long.csv     cohort,tool,sample,genus_taxid,genus_name,reads,rel_abund
  fig3_hostgenus_long_nexvirome.csv  cohort,sample,host_genus,reads,rel_abund   (NexVirome only)
  fig3_richness_summary.csv  per (cohort,tool,sample): n_species, n_genus, total_reads
Run: conda run -n shotgun_virome python result_260605/fig3/code/fig3_realdata_aggregate.py
"""
from __future__ import annotations
import os, sys, glob, sqlite3
import numpy as np, pandas as pd

NX = "/home/share/programs/nexvirome"
RESULT_DIR = f"{NX}/result_260605"
sys.path.insert(0, RESULT_DIR)            # local golden_rule (vir17 exclusion)
sys.path.insert(0, f"{NX}/scripts"); sys.path.insert(0, f"{NX}/notebooks")
sys.path.insert(0, f"{NX}/scripts/benchmark")
from golden_rule import keep_samples
import three_strategy_breadth_lca as TS
from virome_classifier.alignment.filters.filter import MaskingFilter
from virome_classifier.taxonomy import TaxonomyDB
from benchmark_utils import (parse_kreport, kreport_to_species_counts,
                             parse_metabuli_report, metabuli_to_species_counts)

DB = f"{NX}/resources/db/custom/tax_seq_v20260526_MSL41.db"
MASK = f"{RESULT_DIR}/mask/mask_v3_full.bed"
GLEN = f"{NX}/paper/figures/Fig5_extra/tables/species_genome_length.csv"
CACHE = "/tmp/hq_cache"
NAT = f"{NX}/resources/db_20260525/external_native"
BRK = f"{NX}/result_260605/fig3/bracken"   # external-tool Bracken (Kraken2/Metabuli)
OUT = f"{RESULT_DIR}/fig3/tables"
# GOLDEN RULE (result_260605/GOLDEN_RULE.md): method B = best-hit + breadth>=CUT
# + per-taxon read floor n>=MIN_TAXON_READS, rel-abund OFF.
CUT = 0.01            # unmasked breadth (strict; was 0.005)
MIN_TAXON_READS = 3   # per-taxon read floor — drops 1-2 read cross-map FPs
os.makedirs(OUT, exist_ok=True)

# cohort -> sample list (from HQ parquet names)
RNA = keep_samples(sorted(os.path.basename(p)[:-8] for p in glob.glob(f"{CACHE}/A*.parquet")))
DNA = keep_samples(sorted(os.path.basename(p)[:-8] for p in glob.glob(f"{CACHE}/vir*.parquet")))
COHORTS = {"rna": RNA, "dna": DNA}

TS.DB = DB; TS._TAX_PATH = DB; TS._TLEN = TS._load_tlen()
TAX = TaxonomyDB.from_sqlite(DB)
GLEN_MAP = dict(pd.read_csv(GLEN).itertuples(index=False))
_MF = MaskingFilter.from_dataframe(
    pd.read_csv(MASK, sep="\t", header=None, usecols=[0, 1, 2], names=["target", "start", "end"]))

# caches
_name, _genus, _grank = {}, {}, {}
def nm(t):
    if t not in _name: _name[t] = TAX.get_name(int(t)) or f"taxid_{t}"
    return _name[t]
def genus_of(t):
    if t not in _genus:
        try: g = TAX.get_taxid_at_rank(int(t), "genus")
        except Exception: g = None
        _genus[t] = int(g) if g and g > 0 else int(t)
    return _genus[t]
_VIRAL_TAXIDS = None
def _load_viral_taxids():
    """All taxids present in our viral-only DB taxonomy (species + ancestors).
    An external-tool taxid is 'viral' iff it exists here. This is robust to the
    external tools' lineage formatting (e.g. ganon's abfv_rs_cgrg .tre omits the
    Viruses=10239 node, writing viral phyla directly under root as '1||phylum'),
    which a 10239-lineage regex wrongly scored as non-viral."""
    global _VIRAL_TAXIDS
    if _VIRAL_TAXIDS is None:
        con = sqlite3.connect(DB)
        ids = set()
        for (t,) in con.execute("SELECT DISTINCT taxid FROM refseq_sequences"):
            ids.add(int(t))
        con.close()
        _VIRAL_TAXIDS = ids
    return _VIRAL_TAXIDS

def is_viral(t):
    """viral iff the taxid (or its species) is named in our viral-only DB."""
    t = int(t)
    if t in _load_viral_taxids():
        return True
    return TAX.get_name(t) is not None


# ---- host-genus map (NexVirome host roll-up); phage species -> host genus ----
_ENV = ("metagenome", "sludge", "seawater", "sediment", "environment", "soil",
        "wastewater", "uncultured")
def _host_genus_name(h):
    if not h or not h.strip(): return None
    low = h.lower()
    if any(k in low for k in _ENV): return "(metagenome)"
    return h.split()[0]

def build_host_map():
    con = sqlite3.connect(DB)
    acc2host = {a: h.strip() for a, h in con.execute(
        "SELECT accession, host FROM refseq_metadata WHERE host IS NOT NULL AND trim(host)!=''")}
    try:
        for a, h in con.execute("SELECT base_accession, host FROM phage_host_from_title"):
            acc2host.setdefault(a, h)
    except sqlite3.OperationalError:
        pass
    tax2hg = {}
    for acc, tx in con.execute("SELECT accession, taxid FROM refseq_sequences"):
        if tx in tax2hg: continue
        h = acc2host.get(acc.split(".")[0])
        if h:
            hg = _host_genus_name(h)
            if hg: tax2hg[tx] = hg
    con.close()
    return tax2hg
TAX2HOSTGENUS = build_host_map()


# ---- per-tool species counts {taxid: reads} ----
def nexvirome_counts(sample):
    df = TS._add_tlen(pd.read_parquet(f"{CACHE}/{sample}.parquet"))
    best = TS._best_hit(df)
    refs = TS._refs_at(TS._breadth_by_ref(best, _MF), CUT)
    counts = TS._besthit_counts(best[best["target"].isin(refs)])
    return {int(t): int(n) for t, n in counts.items() if n >= MIN_TAXON_READS}  # GOLDEN: n>=3

def bracken_counts(path):
    """Bracken .bracken table: {taxonomy_id: new_est_reads}, kept if the taxid is
    viral. Viral test = is_viral (taxid present in our viral DB) — the lineage
    filter on the Bracken _report breaks because the Metabuli->kraken conversion
    writes the Viruses(10239) node with rank '-', so the indentation-subtree viral
    filter misses everything (e.g. vir22 had 118 viral species but returned 0)."""
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, sep="\t")
    out = {}
    for _, r in df.iterrows():
        tid = int(r["taxonomy_id"]); reads = int(r["new_est_reads"])
        if reads > 0 and is_viral(tid):
            out[tid] = reads
    return out or {}

def ganon_counts(cohort, s):
    # ganon abfv_rs_cgrg .tre: lineage omits the 10239 node, so filter viral by
    # cross-referencing taxids against our viral DB (is_viral) instead of lineage.
    p = f"{NAT}/{cohort}_ganon/{s}.tre"
    if not os.path.exists(p): return {}
    df = pd.read_csv(p, sep="\t", header=None,
                     names=["rank", "taxid", "lin", "name", "u", "sh", "ch", "cum", "pct"])
    sp = df[df["rank"] == "species"]
    return {int(r.taxid): int(r.cum) for r in sp.itertuples()
            if int(r.cum) > 0 and is_viral(int(r.taxid))}

def phanta_counts(cohort, s):
    for pat in [f"{NAT}/{cohort}_phanta/results/classification/{s}.krak.report_bracken_species.filtered",
                f"{NAT}/{cohort}_phanta/results/classification/{s}.krak.report.filtered"]:
        if os.path.exists(pat):
            return kreport_to_species_counts(parse_kreport(pat), filter_virus=True)
    return {}

def tool_counts(tool, cohort, s):
    if tool == "NexVirome": return nexvirome_counts(s)
    if tool == "Kraken2":   return bracken_counts(f"{BRK}/{cohort}_kraken2/{s}.bracken")
    if tool == "Metabuli":  return bracken_counts(f"{BRK}/{cohort}_metabuli/{s}.bracken")
    if tool == "Ganon":     return ganon_counts(cohort, s)
    if tool == "Phanta":    return phanta_counts(cohort, s)
    return {}


def rel_abund(counts):
    """taxid->reads -> taxid->relative abundance = each tool's own read fraction
    among its VIRAL calls (the tool's native percentage, re-normalized to the
    viral subset). Used for the EXTERNAL tools: their NCBI taxids are largely
    absent from our viral-only genome_length table so TPM is not defined for them,
    and read fraction is what every external tool natively reports."""
    s = sum(counts.values())
    return {t: r / s for t, r in counts.items()} if s else {}


def rpk(counts):
    """taxid->reads -> taxid->RPK (reads per kb of genome). NexVirome only:
    its species taxids are 100% covered by the viral genome_length table, so a
    proper genome-length-normalized (TPM-style) relative abundance is well defined.
    Taxa without a length (should be none for NexVirome) fall back to reads as-is."""
    out = {}
    for t, r in counts.items():
        glen = GLEN_MAP.get(int(t))
        out[t] = (r / (glen / 1000.0)) if glen and glen > 0 else float(r)
    return out


def abund_for(tool, counts):
    """NexVirome -> TPM-style (genome-length normalized) relative abundance.
       external tools -> their own read fraction. Both normalized to sum=1."""
    base = rpk(counts) if tool == "NexVirome" else dict(counts)
    s = sum(base.values())
    return {t: v / s for t, v in base.items()} if s else {}


def main():
    TOOLS = ["NexVirome", "Ganon", "Kraken2", "Metabuli", "Phanta"]
    sp_rows, gn_rows, host_rows, summ = [], [], [], []
    for cohort, samples in COHORTS.items():
        for tool in TOOLS:
            for s in samples:
                c = {int(t): int(r) for t, r in tool_counts(tool, cohort, s).items() if r > 0}
                if not c:
                    summ.append(dict(cohort=cohort, tool=tool, sample=s,
                                     n_species=0, n_genus=0, total_reads=0))
                    continue
                # relative abundance: NexVirome = TPM (genome-length normalized),
                # external tools = their own read fraction. Computed at SPECIES level
                # then rolled up so genus/host share the same per-species weights.
                sp_ra = abund_for(tool, c)
                for t, r in c.items():
                    sp_rows.append(dict(cohort=cohort, tool=tool, sample=s, taxid=t,
                                        name=nm(t), reads=r,
                                        rel_abund=round(sp_ra.get(t, 0.0), 6)))
                # genus rollup: reads summed for count, rel_abund summed from species weights
                gc, gra = {}, {}
                for t, r in c.items():
                    g = genus_of(t)
                    gc[g] = gc.get(g, 0) + r
                    gra[g] = gra.get(g, 0.0) + sp_ra.get(t, 0.0)
                for g, r in gc.items():
                    gn_rows.append(dict(cohort=cohort, tool=tool, sample=s, genus_taxid=g,
                                        genus_name=nm(g), reads=r,
                                        rel_abund=round(gra.get(g, 0.0), 6)))
                summ.append(dict(cohort=cohort, tool=tool, sample=s,
                                 n_species=len(c), n_genus=len(gc), total_reads=int(sum(c.values()))))
                # host-genus rollup — NexVirome only; rel_abund re-normalized over hosted taxa
                if tool == "NexVirome":
                    hg, hra = {}, {}
                    for t, r in c.items():
                        h = TAX2HOSTGENUS.get(int(t))
                        if h:
                            hg[h] = hg.get(h, 0) + r
                            hra[h] = hra.get(h, 0.0) + sp_ra.get(t, 0.0)
                    htot = sum(hra.values())
                    for h, r in sorted(hg.items(), key=lambda x: -x[1]):
                        host_rows.append(dict(cohort=cohort, sample=s, host_genus=h, reads=r,
                                              rel_abund=round(hra[h] / htot, 6) if htot else 0.0))
        print(f"  {cohort}: aggregated {len(samples)} samples x {len(TOOLS)} tools", flush=True)

    pd.DataFrame(sp_rows).to_csv(f"{OUT}/fig3_species_long.csv", index=False)
    pd.DataFrame(gn_rows).to_csv(f"{OUT}/fig3_genus_long.csv", index=False)
    pd.DataFrame(host_rows).to_csv(f"{OUT}/fig3_hostgenus_long_nexvirome.csv", index=False)
    pd.DataFrame(summ).to_csv(f"{OUT}/fig3_richness_summary.csv", index=False)

    # quick digest: mean richness per cohort x tool
    sm = pd.DataFrame(summ)
    dig = sm.groupby(["cohort", "tool"]).agg(
        mean_n_species=("n_species", "mean"), mean_n_genus=("n_genus", "mean"),
        mean_reads=("total_reads", "mean")).round(1)
    pd.set_option("display.width", 160)
    print("\n=== mean richness per cohort x tool ===")
    print(dig.to_string())
    print(f"\n-> {OUT}/fig3_{{species,genus}}_long.csv, fig3_hostgenus_long_nexvirome.csv, fig3_richness_summary.csv")


if __name__ == "__main__":
    main()
