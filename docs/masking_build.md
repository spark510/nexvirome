# Masking BED 빌드 절차 (NexVirome)

production viral masking BED를 처음부터 재현하는 전체 절차. 마스킹은 **여러 소스 BED를
모듈식으로 조립**한 것이며, 한 소스만 바꿔 재빌드하고 다시 조립할 수 있다 (전체를 처음부터
다시 만들 필요 없음).

작성: 2026-06-03 / DB: `tax_seq_v20260526_MSL41.db` / 도구: shotgun_virome (blastn/tblastx/dustmasker)

## 5 역할 (스크립트는 역할 prefix로 명명, `scripts/database/`)

| # | 역할 | 무엇을 마스킹 | 스크립트 |
|---|---|---|---|
| 1 | **CONTAMINATION** | 사람/벡터 **오염** (megablast ≥90%) | `build_mmseqs_refdb.py` (STEP 1) |
| 2 | **LOW-COMPLEXITY** | **저복잡도** 반복 (DustMasker) | `build_mmseqs_refdb.py` (STEP 1b) |
| 3 | **HOST-HOMOLOGY** | 바이러스↔**사람 진화적 상동** (blastn ≥75%) | `mask_host_homology_{retro,family}.py` |
| 4 | **CONSERVED-REGION** | **종간 보존부위** cross-map (tblastx/blastp) | `mask_conserved_{retro_gagpol,retro_env,reference,clustering,family}.py` |
| 5 | **ASSEMBLE** | 위 4역할 BED **조립** | `mask_assemble.py` |

각 스크립트 docstring 첫 줄에 `[ROLE n/5 ...]` 태그. 1·2는 게놈 서열 자체 기반(전체 바이러스),
3·4는 BLAST 기반(과별 표적). 1(오염 megablast 90%)과 3(상동 blastn 75%)은 **알고리즘·컷이
달라 겹치지 않는다** (retrovirus에서 base 90%는 0bp / host-homology 75%는 1786bp, 겹침 0).

---

## 0. 핵심 규약 (load-bearing)

- **좌표**: BED 출력은 **0-based half-open**. blastn/tblastx는 1-based inclusive →
  `core.write_bed(zero_based=True)`가 start-1 변환. (`maskbuild/core.py`)
- **모든 reference를 동일하게 처리**: 마스킹은 어느 종이 평가 대상인지 알지 못하며, 어떤
  reference도 사전 보호·제외하지 않는다. 한 ref의 보존부위가 진짜로 host-mimic(사람 상동)·
  저복잡도·과별 보존부위 cross-map이면 마스킹되고, robustness는 KIT(TP 6/6)로 사후 검증한다.
  - ⚠️ **침팬지 CMV(Panine betaHV2) 같은 cross-map FP는 마스킹으로 다루지 않는다.** 그것은
    사람 CMV read가 같은-속 ref로 흘러간 것이라 마스킹 대상이 아니고, **genus roll-up**으로
    같은 속(Cytomegalovirus)으로 합쳐 처리한다 ([[phage_host_processing.md]] / FP typology).
    마스킹은 host-mimic·저복잡도·과별 보존부위 cross-map에만 쓴다.
- **공통 리소스**: CHM13 사람 genome BLAST DB, OUT_DIR=`resources/db_20260525/mmseqs_refdb`.
- **머지**: 같은 ref 인접 구간 gap≤30bp 병합. breadth 계산은 중복 구간 자동 dedupe(차분배열)
  이라 조립 시 overlap 있어도 무해하나, 최종본은 깔끔히 머지 권장.

---

## 1. 소스 BED (production mask 구성)

| 역할 | 스크립트 | 도구 / cutoff | 내용 |
|---|---|---|---|
| **1 오염** | `build_mmseqs_refdb.py` STEP1 (`masking_genome.generate_all`) | **megablast** vs CHM13+UniVec, **id≥90%, len≥150bp** | viral genome 중 사람/벡터 오염 (양방향 4 BLAST). → `viral_mask_<ver>.bed` |
| **2 저복잡도** | `build_mmseqs_refdb.py` STEP1b (`dustmask_bed`) | **DustMasker** level 20, window 64, tract≥60bp | 글리신 반복 등 저복잡도 self-match FP (예: Xantho pIII `(DGGG)n`). 같은 base BED로 합쳐짐 |
| **3 host-homology (retro)** | `mask_host_homology_retro.py` | **blastn** vs CHM13, **pident≥75 & aln≥60** (§3) | retrovirus→사람 진화적 상동(ERV/oncogene). **retrovirus FP 제거의 핵심** |
| **3 host-homology (family)** | `mask_host_homology_family.py` | blastn vs CHM13 (일반화) | 임의 in-scope family→사람 상동 |
| **4 보존 (retro gag/pol)** | `mask_conserved_retro_gagpol.py` | tblastx, evalue 1e-4, min-aa 30 | Retroviridae gag-pol 보존부위 |
| **4 보존 (retro env)** | `mask_conserved_retro_env.py` | tblastx, evalue 1e-4, min-aa 30 | Retroviridae env 보존부위 |
| **4 보존 (reference)** | `mask_conserved_reference.py` | tblastx/blastp, gene-set | 알려진 gene-set (flavi NS3/NS5 등) |
| **4 보존 (clustering)** | `mask_conserved_clustering.py` | blastp 클러스터, CORE_FRAC 0.5, **medID≥50** | herpes 보존 코어 (§2-A) |
| **4 보존 (family 오케스트레이터)** | `mask_conserved_family.py` | reference+clustering+host-homology 묶음 | 과 단위 한번에 |
| **5 조립** | `mask_assemble.py` | concat+sort+dedup | 위 BED들 → production mask |

추가 (보존, polyprotein 군): `analysis/retro_flavi_conserved/step1_matpeptide.py` + `step2_finalize.py`
— GFF mat_peptide + tblastn(qcov≥0.6, pid≥30)로 flavi NS3/NS5·retro Pol/Env 정밀 추출 (§2-B).
※ 현재 analysis 폴더에 있고 maskbuild로 편입 예정.

✅ **저복잡도(DustMasker)는 이제 base 빌드(STEP 1b)에 정식 통합됨** — 과거엔 `scan_repeat_regions_db.py`로
별도 산출만 하고(7,395 ref / 82,627 구간) 어느 mask에도 안 들어가, Xantho pIII `(DGGG)n`
(NC_073753.1 2417-2596) 같은 FP가 처음부터 안 걸렸다. `build_mmseqs_refdb.py`에 `dustmask_bed`
단계를 넣어 base = 오염 + 저복잡도가 됐다 (`--no-dustmask`로 끌 수 있음).

### base 재현 명령 (오염 + 저복잡도)
```
conda run -n shotgun_virome python scripts/database/build_mmseqs_refdb.py \
    --target-fasta resources/db_20260525/ncbi/refseq/viral.1.1.genomic.fna \
    --chm13 /home/share/bowtie2_db/chm13v2.0/chm13v2.0.fa \
    --univec resources/db/ncbi/univec/UniVec \
    --min-length 150 --min-identity 90
# manifest: build_manifest_20260525.json / log: build_mmseqs_refdb.log
```

## 2. 군별 보존부위 제작법 (왜 군마다 다른가)

핵심 분기 = **유전자 경계가 GFF에 어떻게 있나**:

| 군 | 게놈 구조 | 경계 정보 | 보존 표적 |
|---|---|---|---|
| **herpes** | 유전자별 CDS 분리 (HSV-1=82 CDS) | CDS feature (바로 있음) | CD-HIT/medID로 데이터 선별 |
| **flavi** | 단일 polyprotein | mat_peptide feature 필요 | 알려진 NS3/NS5 |
| **retro** | gag-pol 융합 | mat_peptide feature 필요 | 알려진 Pol/Env |

**(A) Herpes — CDS blastp + medID≥50**: subfamily(Alpha/Beta/Gamma)별 CDS 추출 → all-vs-all
blastp → ortholog 클러스터 → **medID**(클러스터의 종간 단백질 일치도 중앙값) ≥50인 유전자만
(22유전자: pol/terminase/MCP/helicase/RNR/UNG/gB/nuclear egress/portal...). 게놈 덮임 중앙값
14%(보존 유전자만, 게놈 86%는 그대로).

**(B) Flavi/Retro — mat_peptide 2단계**: polyprotein이라 blastp 1덩어리 → medID 불가.
- **step1**: mat_peptide 있는 게놈(Flavi 46, Retro 15)은 GFF product로 직접 추출(RefSeq 절단점
  = 정확). Retro Pol은 PR/RT/IN 인접조각(gap≤90nt) 병합 → 완전 Pol.
- **step2**: 없는 게놈은 step1 query로 tblastn → HSP를 reference에 anchor 병합
  (qcov≥0.6 & pid≥30). 기존 tblastx 대비 정밀(Flavi 458→184, Retro 725→156).

**(C) Human-homolog (retrovirus 전용, 필수)**: retrovirus FP의 진짜 원인 = gag/pol/env가 아니라
**사람 oncogene 상동**(v-fos/src/myc). retrovirus genome→CHM13 직접 blastn. 2단계 Pol/Env와
41%만 겹쳐 **별도 필수**. (cutoff는 §3 — 2026-06-03 pident+aln 기반으로 변경)

---

## 3. human-homolog cutoff — evalue vs pident+aln (2026-06-03 변경)

**문제**: 기존 `--evalue 1e-10 --min-aln 60`이 ERV의 **짧고 발산한 끝 구간**(BaboonERV 5'/3',
79~90bp, pident 73~80%)을 놓쳤다. evalue는 길이-편향이라 짧은 상동을 불리하게 본다.

**해결**: evalue 대신 **pident + aln**으로 거른다 (길이 무관, 더 균일한 상동 기준).
sweep 결과 (BaboonERV 8507bp, 누락 4구간 = 77-156 / 2088-2178 / 7750-7825 / 8029-8108):

| pident | aln | 마스크% | 누락4구간 |
|---|---|---|---|
| 85 / 60 | 1% | 0/4 ❌ (기존류) |
| 80 / 60 | 5% | 3/4 |
| **75 / 60** | **21%** | **4/4 ✅ (선택)** |
| 70 / 50 | 33% | 4/4 (과마스킹↑) |

→ **`--min-pident 75 --min-aln 60 --evalue 1000`이 이제 스크립트 기본값** (evalue 사실상 해제,
pident/aln만). 과마스킹 21%로 억제하며 ERV 끝까지 커버. ⚠️ pident 75는 일부 외인성 retro
정상부위도 마스킹할 수 있으니 **재빌드 후 KIT TP 6/6 + 실데이터 진짜(EBV/Vientovirus) 보존
반드시 검증**.

빌드 명령 (기본값이 pident75/aln60이므로 인자 생략 가능):
```
conda run -n shotgun_virome python scripts/database/mask_host_homology_retro.py --date 2026-06-03
```

---

## 4. 최종 조립 (`mask_assemble.py`)

수동 `cat` 금지 — 이 스크립트가 concat → sort(acc,start) → exact-dup 제거 → write.
(col7 source label은 provenance로 보존. 모든 ref 동일 처리, 사전 제외 없음.)

```
conda run -n shotgun_virome python scripts/database/mask_assemble.py \
    --base resources/db_20260525/mmseqs_refdb/viral_mask_<ver>.bed \
    --bed retroviridae_gagpol_<...>.bed \
    --bed retroviridae_env_<...>.bed \
    --bed retroviridae_humanhomolog_CHM13_2026-06-03.bed \
    --bed flaviviridae_ns_<...>.bed \
    --bed <herpes_core>.bed \
    --date 2026-06-03
```

⚠️ **현재 `mask_v2_2step.bed`(2393)는 herpes 디렉토리에서 수동 cat+awk로 만든 것이다.**
production 승격 전 반드시 `mask_assemble.py`로 재조립할 것 (정렬·dedup·provenance).

---

## 5. 검증 (조립 후 필수)

1. **KIT 성능**: `kit_mock_uniformity.py` — TP 6/6 유지(어떤 ref가 마스킹됐든 진짜는 게놈
   나머지로 검출돼야 robust), FP 변화 확인.
2. **실데이터 진짜 보존**: EBV/HHV-7/Vientovirus 등 검출 유지 확인.
3. **과마스킹 점검**: family별 median masked fraction (`compare_conserved_masks.py`,
   >80%면 경고).

검증 이력 (FINAL_MASKING_SUMMARY, 2026-06-03, mask_v2_2step):
- KIT: baseline = 2step (마스킹이 진짜 양성 무해). 남은 FP=Tobamovirus(식물)·FN=Reovirus 둘 다 마스킹 무관.
- 실데이터: RNA 9개·DNA 114개 검출 사라짐 = **100% retrovirus oncogene cross-map FP**
  (Feline/Avian sarcoma, Abelson MLV...), **진짜 손실 0** (EBV 등 전부 보존).
- 핵심 통찰: cross-map FP는 보존 유전자에 붙고 진짜는 게놈 전체에 붙음 → 보존부위 마스킹해도
  진짜는 나머지로 breadth 확보해 생존. retrovirus FP의 진짜 원인은 human-homolog(oncogene).

---

## 관련 파일

- 빌드 (역할별): `scripts/database/build_mmseqs_refdb.py` (1 오염 + 2 저복잡도),
  `mask_host_homology_{retro,family}.py` (3 상동), `mask_conserved_*.py` (4 보존)
- 공통: `scripts/database/maskbuild/{core,blast}.py`
- 조립: `scripts/database/mask_assemble.py` (5)
- (과거 별도 DustMasker 도구: `scripts/benchmark/scan_repeat_regions_db.py` — 이제 base 빌드에 통합)
- 검증: `scripts/benchmark/{kit_mock_uniformity,compare_conserved_masks}.py`
- 산출: `resources/db_20260525/mmseqs_refdb/` (소스 BED), `analysis/herpes_conserved/bed/mask_v2_*.bed` (실험 조립)
- 관련 문서: [[phage_host_processing.md]], `paper/notes/2026-05-29_phage_dropout_root_cause_and_mask_strategy.md`
