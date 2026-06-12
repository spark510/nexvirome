# Masking Filter - OOP Refactored Version (v2.0)

깔끔하게 리팩토링된 OOP 기반 masking filter 패키지입니다.

## 📁 폴더 구조

```
masking_filter/
├── __init__.py              # Public API exports
├── models.py                # Data classes (CoverageStats, FilterResult, etc.)
├── mask_loader.py           # MaskedRegion, BED file loading
├── coverage_calculator.py   # Coverage calculation logic
├── filter.py                # MaskingFilter main class
├── legacy.py                # Backwards compatibility wrappers
└── README.md                # This file
```

## 🎯 주요 개선사항

### OOP 관점 개선

| 기존 문제점 | 해결 방법 |
|-----------|---------|
| ❌ 모든 함수가 전역 함수 | ✅ `MaskingFilter` 클래스로 캡슐화 |
| ❌ `mask_dict`를 매번 인자로 전달 | ✅ 인스턴스 변수로 관리 |
| ❌ `pd.Series`로 통계 반환 (타입 안정성 없음) | ✅ `CoverageStats` 데이터 클래스 사용 |
| ❌ 계산/필터링/로깅 로직 혼재 | ✅ 단일 책임 원칙(SRP) 준수 |
| ❌ 코드 중복 (80+ 라인) | ✅ 단일 진실 공급원(SSOT) |

### 코드 품질 개선

- **타입 힌트**: 모든 함수에 완전한 타입 힌트
- **불변 객체**: `@dataclass(frozen=True)` 사용
- **명확한 네이밍**: `calculate_merged_coverage` (기존: `calculate_total_coverage`)
- **Factory 메서드**: `MaskingFilter.from_bed_file()`
- **Fluent API**: 메서드 체이닝 가능

## 🚀 Quick Start

### 기본 사용법

```python
from masking_filter import MaskingFilter

# 1. BED 파일에서 필터 생성
filter = MaskingFilter.from_bed_file("masked_regions.bed")

# 2. Unmasked coverage로 필터링
result = filter.filter_by_unmasked_coverage(hits_df, min_coverage=0.25)

# 3. 결과 확인
print(result.summary())
print(f"Passed: {result.n_passed_targets} targets")

# 4. 통과한 hits 가져오기
passed_hits = result.passed
```

### 다양한 필터링 옵션

```python
# Total coverage 필터 (masked 포함)
result = filter.filter_by_total_coverage(hits_df, min_coverage=0.5)

# Hybrid 필터 (둘 다 만족해야 함)
result = filter.filter_by_hybrid_coverage(
    hits_df,
    min_total_cov=0.5,
    min_unmasked_cov=0.25
)

# Overlap threshold 조정
filter = MaskingFilter.from_bed_file("file.bed", overlap_threshold=0.7)
```

### 통계 계산만 하기

```python
# 필터링 없이 통계만 계산
stats_df = filter.calculate_stats(hits_df)

print(stats_df[['total_coverage_ratio', 'unmasked_coverage_ratio']])
```

### 필터 체이닝

```python
# 여러 필터 순차 적용
result1 = filter.filter_by_unmasked_coverage(df, min_coverage=0.02)
result2 = filter.filter_by_total_coverage(result1.passed, min_coverage=0.03)

final_hits = result2.passed
```

## 📊 데이터 모델

### CoverageStats

단일 타겟의 포괄적인 통계:

```python
stats = CoverageStats(
    target_name="NC_001806.2",
    target_length=100000,
    taxid=10376,
    has_mask=True,
    total=CoverageMetrics(...),      # 전체 hits
    masked=CoverageMetrics(...),     # Masked 영역의 hits
    unmasked=CoverageMetrics(...)    # Unmasked 영역의 hits
)

# 접근 방법
print(stats.total.coverage_ratio)
print(stats.unmasked.avg_depth)
```

### CoverageMetrics

Coverage 카테고리별 메트릭스:

```python
metrics = CoverageMetrics(
    hit_count=5,              # Hit 개수
    breadth_bp=4003,          # 커버된 bp 수
    coverage_ratio=0.04,      # Coverage 비율 (merged intervals)
    avg_depth=1.25            # 평균 depth
)
```

### FilterResult

필터링 결과:

```python
result = FilterResult(
    passed=passed_df,         # 통과한 hits
    failed=failed_df,         # 실패한 hits
    stats=stats_df,           # 통계 DataFrame
    filter_name="unmasked_coverage"
)

# 편리한 속성들
result.n_passed_targets
result.n_failed_targets
result.summary()  # 사람이 읽기 좋은 요약
```

## 🔄 마이그레이션 가이드

### 기존 코드 → 새 코드

#### 1. BED 파일 로딩

```python
# OLD
mask_dict = load_masked_bed("file.bed")

# NEW
filter = MaskingFilter.from_bed_file("file.bed")
```

#### 2. Unmasked coverage 필터

```python
# OLD
passed, failed, stats = apply_unmasked_cov_filter(
    df, mask_dict, min_unmasked_cov=0.25
)

# NEW
result = filter.filter_by_unmasked_coverage(df, min_coverage=0.25)
passed = result.passed
failed = result.failed
stats = result.stats
```

#### 3. Total coverage 필터

```python
# OLD
passed, failed, stats = apply_total_cov_filter(
    df, mask_dict, min_total_cov=0.5
)

# NEW
result = filter.filter_by_total_coverage(df, min_coverage=0.5)
```

#### 4. Hybrid 필터

```python
# OLD
passed, failed, stats = apply_hybrid_cov_filter(
    df, mask_dict, min_total_cov=0.5, min_unmasked_cov=0.2
)

# NEW
result = filter.filter_by_hybrid_coverage(
    df, min_total_cov=0.5, min_unmasked_cov=0.2
)
```

### Legacy API 사용 (권장하지 않음)

기존 코드를 수정하기 어려운 경우, legacy wrapper를 사용할 수 있습니다 (deprecation warning 표시됨):

```python
from masking_filter.legacy import (
    load_masked_bed,
    apply_unmasked_cov_filter,
    apply_total_cov_filter,
    apply_hybrid_cov_filter
)

# 기존 코드 그대로 작동 (warning 출력됨)
mask_dict = load_masked_bed("file.bed")
passed, failed, stats = apply_unmasked_cov_filter(df, mask_dict, 0.25)
```

## 🧪 테스트

```bash
cd <repo>/scripts/virome_engine/gpt
python test_new_masking_filter.py
```

테스트 커버리지:
- ✅ 기본 사용법
- ✅ Unmasked coverage 필터
- ✅ Total coverage 필터
- ✅ Hybrid 필터
- ✅ Overlap threshold 효과
- ✅ 필터 체이닝
- ✅ CoverageStats 객체 API
- ✅ Legacy API 호환성

## 📖 API Reference

### MaskingFilter

**생성자**

```python
MaskingFilter(
    mask_dict: Optional[Dict[str, MaskedRegion]] = None,
    overlap_threshold: float = 0.5
)
```

**Factory Methods**

```python
@classmethod
MaskingFilter.from_bed_file(bed_file: str, overlap_threshold: float = 0.5)

@classmethod
MaskingFilter.from_dataframe(df: pd.DataFrame, overlap_threshold: float = 0.5)
```

**Properties**

```python
filter.n_targets_with_masks  # int: 마스크 데이터가 있는 타겟 수
filter.has_masks             # bool: 마스크 데이터 존재 여부
```

**Methods**

```python
filter.calculate_stats(df: pd.DataFrame) -> pd.DataFrame

filter.filter_by_unmasked_coverage(
    df: pd.DataFrame,
    min_coverage: float = 0.25
) -> FilterResult

filter.filter_by_total_coverage(
    df: pd.DataFrame,
    min_coverage: float = 0.5
) -> FilterResult

filter.filter_by_hybrid_coverage(
    df: pd.DataFrame,
    min_total_cov: float = 0.5,
    min_unmasked_cov: float = 0.25
) -> FilterResult
```

## 💡 설계 철학

### SOLID Principles 준수

1. **Single Responsibility Principle (SRP)**
   - `CoverageCalculator`: 계산만
   - `MaskingFilter`: 필터링만
   - `MaskLoader`: 로딩만

2. **Open/Closed Principle (OCP)**
   - 새로운 필터 타입 추가 가능 (기존 코드 수정 없이)

3. **Liskov Substitution Principle (LSP)**
   - `FilterResult` 객체는 항상 동일한 인터페이스 제공

4. **Interface Segregation Principle (ISP)**
   - 작은 인터페이스로 분리 (`CoverageMetrics`, `CoverageStats`)

5. **Dependency Inversion Principle (DIP)**
   - 구체 클래스 대신 추상화(데이터 클래스)에 의존

### 불변성 (Immutability)

모든 데이터 클래스는 불변:

```python
@dataclass(frozen=True)
class CoverageMetrics:
    # 생성 후 수정 불가
    # Thread-safe, 예측 가능
```

### 타입 안정성

완전한 타입 힌트로 IDE 자동완성 및 정적 분석 지원:

```python
def filter_by_unmasked_coverage(
    self,
    df: pd.DataFrame,
    min_coverage: float = 0.25
) -> FilterResult:
    ...
```

## 🐛 알려진 제약사항

1. **메모리 사용**: `calculate_breadth_and_depth()`는 전체 target 길이만큼의 배열 생성
   - 매우 긴 genome (>100Mb)의 경우 메모리 주의

2. **Pandas 의존성**: 현재는 pandas DataFrame에 의존
   - 향후 polars 등 다른 DataFrame 라이브러리 지원 고려 가능

## 📝 버전 히스토리

- **v2.0.0** (2024-11): OOP 기반 완전 리팩토링
  - 클래스 기반 설계
  - 데이터 클래스 도입
  - Legacy API 호환성 유지

- **v1.x** (2024-08~10): 절차적 프로그래밍 버전
  - 기존 `masking_filter.py`

## 🤝 기여 가이드

새로운 필터 타입을 추가하려면:

1. `MaskingFilter` 클래스에 메서드 추가
2. `FilterResult` 반환
3. 테스트 작성

```python
def filter_by_custom_logic(
    self,
    df: pd.DataFrame,
    custom_param: float
) -> FilterResult:
    stats_df = self.calculate_stats(df)

    # Custom filtering logic
    valid_mask = (stats_df["unmasked_avg_depth"] > custom_param)
    valid_targets = stats_df[valid_mask].index

    passed, failed = self._split_by_targets(df, valid_targets)

    return FilterResult(
        passed=passed,
        failed=failed,
        stats=stats_df,
        filter_name="custom_logic"
    )
```

## 📧 문의

질문이나 버그 리포트는 이슈로 등록해주세요.
