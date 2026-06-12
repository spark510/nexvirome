# EM Mode vs LCA Mode 비교 정리

Date: 2026-05-30

---

## 1. EM Mode 기본 파라미터

| 파라미터 | 기본값 | 설명 |
|---|---|---|
| `--em-iterations` | 20 | EM 반복 횟수 |
| `--em-convergence` | 1e-6 | 수렴 임계값 |

정의 위치: [scripts/virome_classifier/cli/classify.py:266-276](../scripts/virome_classifier/cli/classify.py)

---

## 2. EM 알고리즘 동작 방식

EM mode는 **breadth coverage threshold 기반의 "다음 후보로 보내는" 로직이 없다.**

작동 순서:
1. Quality filter (`fident`, `alnlen`, `qcov`, `evalue`) 통과
2. Masking filter 통과
3. **E-step**: multi-mapping reads를 현재 abundance × alignment score 비율로 fractional assignment
4. **M-step**: 분배된 counts로 species abundance 재추정 (genome-length normalization)
5. 수렴하거나 max_iterations 도달 시 종료
6. 최종: 각 read를 posterior probability 최대 species로 할당

LCA와의 핵심 차이: multi-mapping reads를 상위 분류로 올리는 대신, **전체 abundance 패턴 기반으로 species에 분배**

코드: [scripts/virome_classifier/classification/em_classifier.py](../scripts/virome_classifier/classification/em_classifier.py)

---

## 3. KIT Mock Community 벤치마크 결과

대상: MagNA_1, MagNA_2, Qiagen_1, Qiagen_2 (4샘플 평균)  
DB: tax_seq_v20260526_MSL41.db  
결과 위치: `/tmp/multimap_bench/`

| strategy | Precision | F1 | FP | TP |
|---|---|---|---|---|
| coverage_all | 0.581 | 0.682 | 3.8 | 5.0 |
| coverage_best_hit | **0.601** | 0.712 | **3.5** | 5.2 |
| coverage_local_depth | 0.581 | 0.682 | 3.8 | 5.0 |
| **em** | 0.576 | **0.716** | 4.5 | **5.8** |

KIT truth = 6 species

**해석:**
- EM은 TP가 가장 높고 F1도 가장 높지만, FP도 4.5로 가장 많음 (Precision 최저)
- `coverage_best_hit`이 Precision/FP 기준으로 가장 깔끔
- EM은 multi-mapping reads를 적극 할당 → TP↑ but FP↑ (trade-off)
- 이 벤치마크에는 LCA mode가 포함되지 않음

벤치마크 스크립트: [scripts/benchmark/bench_multimap_modes.py](../scripts/benchmark/bench_multimap_modes.py)

---

## 4. 두 모드 전체 파라미터 비교

### 공통 Quality Filter 파라미터 (LCA, EM 모두 동일)

| 파라미터 | LOCKED 기본값 | 설명 |
|---|---|---|
| `--min-identity` | 0.85 | fident |
| `--min-length` | 60 | alnlen |
| `--min-query-coverage` | 0.5 | qcov |
| `--max-evalue` | 1e-3 | e-value |
| `--breadth` | 0.05 | breadth threshold (LCA postfilter) |
| `--rel-abund` | 0.0005 | relative abundance threshold |

> LOCKED 2026-05-29 기준 (memory: `project_final_hitqual_params.md`)

### LCA mode 전용

| 파라미터 | 기본값 |
|---|---|
| `--lca-fix-rank` | none |
| `--mode` | lca |

### EM mode 전용

| 파라미터 | 기본값 |
|---|---|
| `--em-iterations` | 20 |
| `--em-convergence` | 1e-6 |

> **주의**: EM mode는 내부적으로 breadth filter를 사용하지 않음. postfilter(breadth ≥ 0.05)는 EM 이후 apply_hitqual_filter 단계에서 적용됨.

---

## 5. 실행 환경

```bash
conda activate shotgun_virome
# 또는
/usr/local/bin/miniconda3/envs/shotgun_virome/bin/python3
```
