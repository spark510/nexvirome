# Phage host 처리 (NexVirome)

phage 검출을 **host 박테리아/고세균 속(genus)** 단위로 정리하는 방법과 근거를 정리한 문서.
ICTV VMR과 **독립적인** host 소스를 사용하며, phage cross-map 분산을 host 단위로 통합하고
마이크로바이옴 해석을 직접 가능하게 한다.

작성일: 2026-06-03 / DB: `tax_seq_v20260526_MSL41.db`

---

## 1. 배경 — 왜 host로 정리하나

- DB에 phage가 **6,739종 / 6,821서열 (전체의 35%)** 존재. class 분포:
  Caudoviricetes 5,504 · Leviviricetes(ssRNA) 941 · Microviricetes(ssDNA) 61 · Cystoviridae 21 …
- 실데이터 phage 검출은 **cross-map 분산이 심하다**: 한 read가 같은 host의 여러 phage ref에
  multi-map되어 한 종이 수십 개 가짜 종으로 부풀려진다 (RNA phage median 19종, DNA 37종/샘플).
- **genus가 같으면 host가 같다 (93%)**, **genus 내 길이도 일정하다 (89% CV≤0.05)** — 검증됨.
  → phage를 **host 속**으로 roll-up하면 분산이 크게 통합되고, "어떤 박테리아의 phage인가" =
  마이크로바이옴 신호로 직접 해석된다. (genus roll-up보다 강력: vir16 91종→genus 30→host 9)

---

## 2. host 정보 소스 (3가지, VMR 독립이 핵심)

| 소스 | 채움 | 정밀도 | 비고 |
|---|---|---|---|
| **refseq_metadata.host** | 79% (taxid join) | 종 수준 ("Escherichia coli") | **NCBI Virus 포털 "Host" 필드. VMR 독립. 주 소스** |
| **phage_host_from_title** | +211종 (blank 보충) | 속 수준 | title 파싱으로 채운 보충 테이블 (아래 §3) |
| ictv_vmr.host_source | 100% (VMR 등록종) | 거침 (bacteria/plant) | MSL41 원본과 일치(정상). 거친 분류라 보조용 |

- **refseq_metadata 출처**: `ncbi_refseq_table.tsv` (NCBI Virus 포털 검색결과 수동 export,
  18,962행) → `import_ncbi_refseq_full.py`가 `Host` 컬럼을 `refseq_metadata.host`로 적재.
  **자동 fetch 목록(refseq_fasta/taxdump/acc2taxid/ICTV)에 없는 수동 파일.**
- **⚠️ accession 버전 함정**: `refseq_sequences.accession`은 버전 있음(`AC_000001.1`),
  `refseq_metadata.accession`은 버전 없음(`AC_000001`). **버전 제거 join 필수** (안 하면 0% 매칭).

---

## 3. title 기반 phage host 보충

`scripts/database/fill_phage_host_from_title.py` — metadata.host가 blank인 phage(728종)를
title/lineage에서 host를 파싱해 채운다. **결과는 `phage_host_from_title` 테이블에 기록**
(dry run은 테이블만 쓰고 metadata는 안 건드림; `--apply`로 blank만 채움, 기존값 절대 덮지 않음).

파싱 규칙 (우선순위):
1. **lineage clade host** (가장 신뢰): lineage에 `Crassvirales`/`Suoliviridae` → **Bacteroides**
   (crAssphage. "CrAssphage cr50_1"·"Uncultured phage cr35_1" 두 형식 모두 lineage로 잡힘)
2. **synthetic `<X>phage` 어간**: cyanophage→Cyanobacteria, mycobacteriophage→Mycobacterium …
3. **`<Host> phage` title 패턴**: 앞 단어를 host 후보로, 노이즈(ssRNA/uncultured/thermophilic/
   stx-converting/morphology) 제거 + 선행 수식어 strip + **박테리아/고세균 속명 검증**
   (ncbi_taxonomy Bacteria(2)+Archaea(2157) 자손 이름과 대조 — 서술어를 host로 만들지 않음)

결과: blank 728종 중 **211종(29%) host 파싱** (Bacteroides 63·Klebsiella 66·Escherichia 22 …).
나머지 517종은 `no_host_in_title` — **title에 host가 없는 게 정상** (대부분 ssRNA phage
"SsRNA phage Esthiorhiza.1_1", uncultured/MAG 환경 서열).

---

## 4. host 미상 phage (정상적으로 남는 것)

phage 6,739종 중 host 미상 **약 507종 (7.5%)** — 채울 수 없는 게 맞다:

| class | 수 | 정체 |
|---|---|---|
| Leviviricetes (ssRNA phage) | ~470 | 환경 메타게놈 TPA 서열. host 규명된 적 없음 |
| Microviricetes (ssDNA) | 34 | Chimpanzee faeces / Marine gokushovirus / Bog — 환경, host 미상 |
| virophage | 3 | Yellowstone Lake virophage (거대바이러스 기생, host가 박테리아 아님) |

→ host roll-up에서 이들은 **`(host_unknown)` 그룹**으로 둔다. 실데이터 영향 작음
(host_unknown reads: RNA 898 / DNA 4301 vs 상위 host Streptococcus 24996 / Pseudomonas 135317).

**주의**: 제 host-미상 조사 스크립트가 한때 `Duplodnaviria`를 phage clade로 넣어 herpesvirus
49종(Murine/Bovine herpesvirus 등)을 phage로 오분류했으나, **production `BIO.host_class()`는
Caudoviricetes/Microviridae/Inoviridae만 phage로 보고 Duplodnaviria는 안 봄 → 버그 없음**.
herpes를 phage로 세지 않는다 (검증: Porcine lymphotropic herpesvirus 3 → human/vertebrate).

---

## 5. 분절·중복 처리 (length 정규화의 전제)

- **phage 분절은 7종뿐** — 전부 Cystoviridae (Pseudomonas phage phi6/8/12/13/phiNN/phi2954/phiYY,
  L/M/S 3분절 dsRNA). phage 중 거의 유일한 분절 그룹 (분절은 본질적으로 진핵 바이러스 특징).
- **taxid 단위 합산**으로 분절·변이체를 한 종으로: 같은 taxid의 모든 accession length 합
  = 종 게놈 길이 (Cystoviridae L+M+S ≈ 13kb). best-hit 카운트(`_besthit_counts`)는 taxid
  value_counts라 read의 species 합산은 자동.
- **유전자조각/partial 제외**: `species_genome_length.csv` 빌드 시 `NG_` 접두사 + title
  `gene for`/`partial genome` 98개를 length 합에서 제외. **단 그게 유일한 ref인 종(250개)은
  보존**(length=0 방지). `expected_genome_length`는 자기 length 복사라 partial 판정에 무용 —
  **title이 유일한 판정 도구**.

---

## 6. roll-up + TPM (구현)

`scripts/benchmark/phage_host_rollup_analysis.py`:
1. 검출: best-hit + breadth≥0.005 + n≥3, rel-abund OFF (mask_v2_2step)
2. **phage → host 속 roll-up** (read 합산), **non-phage → species 유지**
   (herpes/EBV 등 같은속 진짜 종 보존; `host_class=="phage"`일 때만 roll-up)
3. TPM length 정규화: `species_genome_length.csv` 사용 — host 그룹 length = 멤버 phage 종들의
   평균 게놈 길이 (조각 제외돼 평균≈median)

검증 결과:
- **KIT TP 6/6 보존** (host roll-up이 진짜 6종 안 죽임; phage FP 거의 0)
- 실데이터 통합: RNA phage host 속/샘플 median **4** (종 19 → 5배↓), DNA median **5** (종 37 → 7배↓)
- 해석: RNA = Streptococcus/Actinomyces/Rothia/Fusobacterium(구강·기도 상재균),
  DNA = Pseudomonas(135k)/Streptococcus/Klebsiella(기도 박테리아) — 마이크로바이옴 직접 해석

---

## 관련 파일

- `scripts/database/fill_phage_host_from_title.py` — title→host 보충 (→ `phage_host_from_title` 테이블)
- `scripts/benchmark/phage_host_rollup_analysis.py` — host 속 roll-up + TPM 분석
- `paper/figures/Fig5_extra/tables/species_genome_length.csv` — 종 게놈 길이 (조각 제외, 분절 합산)
- `paper/figures/Fig5_extra/tables/phage_host_rollup_{kit,realdata}.csv` — 결과
