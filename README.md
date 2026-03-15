# HLM_app - 다층모형분석 (Multilevel Modeling) 웹 애플리케이션

> **접속 주소**: [https://bernard-paintball-deaf-downloading.trycloudflare.com](https://bernard-paintball-deaf-downloading.trycloudflare.com)
>
> ※ Cloudflare Tunnel 기반 — 서버 재시작 시 URL이 변경될 수 있습니다.

R `lme4`/`lmerTest` 패키지 기반 다층모형(HLM/MLM) 분석을 브라우저에서 수행할 수 있는 웹 애플리케이션입니다. CSV 파일을 업로드하면 기저모형부터 교차수준 상호작용 모형까지 자동으로 적합하고, 결과를 시각적으로 확인하며, AI 보고서까지 생성할 수 있습니다.

---

## 주요 기능

### 1. 데이터 업로드 및 전처리
- CSV 파일 드래그 앤 드롭 또는 파일 선택 업로드
- 데이터 미리보기 테이블 (상위 5행)
- 수치형/범주형 변수 자동 감지

### 2. 기초 통계 분석 (5개 하위 탭)
| 탭 | 기능 |
|---|---|
| **기술통계** | N, 결측, 평균, SD, 최솟값, 최댓값, 왜도, 첨도 |
| **분포** | 히스토그램 + 정규분포 곡선 오버레이 (M, SD 표시) |
| **상관행렬** | Pearson 상관계수 히트맵, **변수 선택/제외 체크박스** 지원 |
| **집단별 비교** | 집단변수 × 수치변수 막대그래프 (평균 비교) |
| **시각화** | 산점도, 박스플롯, 히스토그램 등 다양한 차트 |

### 3. 분석 설정
- **종속변수 (Y)** / **집단변수 (Level-2 ID)** 선택
- **Level-1 독립변수**, **Level-2 공변량** 체크박스 선택
- **무선기울기 변수** 선택 (Level-1 변수 중)
- **공분산 구조 선택**:
  - **Unstructured (비구조적)**: 무선효과 간 분산과 공분산을 모두 추정 (lme4: `|`)
  - **Independent (독립)**: 분산만 추정, 공분산 = 0으로 제약 (lme4: `||`)
- **교차수준 상호작용** (L1 × L2) 설정
- **센터링 옵션**:
  - Level-1: 없음 / Grand Mean Centering (GMC) / Group Mean Centering (CWC)
  - Level-2: 없음 / Grand Mean Centering (GMC)
  - CWC 선택 시 집단평균이 자동으로 L2 공변량에 추가됨

### 4. 다층모형 분석 (4단계 모형)
| 모형 | 설명 | R 공식 예시 |
|---|---|---|
| **M0** 기저모형 (Null) | 무선절편만 포함, ICC 계산 | `Y ~ 1 + (1\|group)` |
| **M1** 무선절편 (Random Intercept) | 고정효과 추가, 기울기는 고정 | `Y ~ X1 + W1 + (1\|group)` |
| **M2** 무선기울기 (Random Slope) | 집단별 기울기 차이 허용 | `Y ~ X1 + W1 + (1 + X1\|group)` |
| **M3** 교차수준 상호작용 | L1×L2 상호작용 효과 | `Y ~ X1 + W1 + X1:W1 + (1 + X1\|group)` |

#### 추정 방식 (Stata mixed 방식)
- **EM 알고리즘 초기값** + **REML/ML 이중 추정**
  - 1단계: EM 알고리즘 (25회 반복)으로 분산 성분의 안정적인 초기값 추정
  - 2단계: lmer()의 Newton-Raphson/BOBYQA 옵티마이저로 최종 수렴
- **REML 추정** (보고용): 분산 성분의 비편향 추정치 (HLM7, SPSS, Stata 기본값)
- **ML 추정** (모형비교용): 중첩 모형 간 우도비 검정(LRT)에 사용

#### 각 모형 결과에 포함되는 정보:
- ICC (집단내 상관계수)
- R² marginal / conditional (Nakagawa & Schielzeth, 2013)
- AIC, BIC, Deviance, Log-Likelihood
- **고정효과 테이블**: 추정값(B), 표준오차(SE), t값, **p값** (Satterthwaite), 95% CI
- **무선효과 테이블 (Stata 스타일)**:
  - Estimate (분산/공분산 추정치)
  - Std. Err. (delta method 표준오차)
  - 95% CI 하한/상한 (profile likelihood 신뢰구간)
  - Unstructured 선택 시 공분산(cov) 행도 표시
- **LR test vs. linear model**: 혼합모형 vs OLS 우도비 검정 (Stata 스타일)
- **우도비 검정 (LRT)**: 이전 모형과의 χ² 비교

### 5. 모형 비교
- AIC, BIC, Deviance, △Deviance, χ², df, p값을 한 눈에 비교하는 표
- AIC/BIC 막대 차트
- 최적 모형 자동 표시 (★)

### 6. 모형 공식 탭 (수학적 표기)
- **Level 1** (개인 수준) 방정식: Y_ij = β₀ⱼ + β₁ⱼ·X₁ᵢⱼ + ... + eᵢⱼ
- **Level 2** (집단 수준) 방정식: β₀ⱼ = γ₀₀ + γ₀₁·W₁ⱼ + u₀ⱼ
- **통합 모형** (Combined): 모든 항을 하나의 방정식으로 결합
- 각 모형의 **R lmer() 코드** 표시
- 모형 선택 버튼으로 M0~M3 전환
- 범례: 고정효과(γ, 파란색), 무선효과(u, e, 빨간색) 구분

### 7. 진단 플롯 (9종)
| 차트 | 설명 |
|---|---|
| 잔차 히스토그램 | 표준화 잔차의 분포 + 정규곡선 |
| 잔차 vs 예측값 | 등분산성 확인 산점도 |
| 잔차 Q-Q | 잔차 정규성 검정 |
| 잔차 vs 독립변수 | 선형성 확인 (변수 선택) |
| 무선효과 히스토그램 | BLUP 분포 (u₀, u₁ 등 선택) |
| 무선효과 Q-Q | 무선효과 정규성 검정 |
| 회귀선 (전체) | 산점도 + OLS 회귀선 (X변수 선택) |
| 회귀선 (집단별) | 집단별 색상 구분 회귀선 (표시 집단 수 조절) |
| 변수 쌍 행렬 (Pairs Plot) | ggpairs 스타일: 대각선=히스토그램, 하삼각=산점도, 상삼각=상관계수 |

### 8. AI 보고서 생성
- **지원 AI**: OpenAI (GPT-4o), Anthropic (Claude), 로컬 LLM (Ollama 등)
- **보고서 스타일**: 학술 논문 / 요약 / 상세 분석
- **언어**: 한국어 / English
- **분석에 사용된 실제 R 코드**가 상세한 한국어 주석과 함께 AI 프롬프트에 포함
- 보고서 **Word(.docx) 다운로드** / 클립보드 복사

---

## 기술 스택

| 구성요소 | 기술 |
|---|---|
| **백엔드** | R 4.x + [Plumber](https://www.rplumber.io/) REST API |
| **통계 엔진** | [lme4](https://cran.r-project.org/package=lme4) + [lmerTest](https://cran.r-project.org/package=lmerTest) |
| **추정 방식** | EM 알고리즘 초기값 + REML/ML 이중 추정 (Stata mixed 방식) |
| **프론트엔드** | Vanilla HTML/CSS/JavaScript (프레임워크 없음) |
| **차트** | [Chart.js 4.x](https://www.chartjs.org/) |
| **CSV 파싱** | [PapaParse 5.x](https://www.papaparse.com/) |
| **Word 생성** | [docx.js 8.x](https://docx.js.org/) + FileSaver.js |

---

## 파일 구조

```
HLM_app/
├── api.R              # R Plumber API (분석 로직 전체)
│                       #   - em_lmer(): EM 알고리즘 초기값 + lmer 적합
│                       #   - calc_lrt_vs_lm(): LR test vs. linear model
│                       #   - calc_icc(), calc_r2(): ICC, R² 계산
│                       #   - extract_fixed(): 고정효과 추출 (B, SE, t, p, CI)
│                       #   - extract_random(): 무선효과 추출 (Stata 스타일)
│                       #   - extract_diagnostics(): 진단 플롯 데이터
│                       #   - /api/analyze: 4단계 모형 순차 적합
├── www/
│   └── index.html     # 프론트엔드 (HTML + CSS + JS 올인원)
│                       #   - 데이터 업로드/미리보기
│                       #   - 기초 통계 (기술통계, 분포, 상관, 시각화)
│                       #   - 분석 설정 (변수, 센터링, 공분산 구조)
│                       #   - 결과 테이블/차트/진단플롯
│                       #   - AI 보고서 생성
├── start_server.R     # R 서버 시작 스크립트 (작업 디렉토리 자동 감지)
├── run.bat            # Windows 더블클릭 실행 배치파일
├── install.R          # 필수 R 패키지 설치 스크립트
├── check_pkg.R        # 패키지 설치 확인
├── .gitignore
└── README.md
```

---

## 설치 및 실행

### 사전 요구사항
- **R 4.0 이상** ([다운로드](https://cran.r-project.org/))
- `Rscript`가 PATH에 등록되어 있어야 합니다
  - Windows: `C:\Program Files\R\R-x.x.x\bin\` 이 PATH에 포함되었는지 확인
  - macOS/Linux: 보통 자동 등록됨

### 1단계: 패키지 설치

```bash
Rscript install.R
```

필요 패키지: `plumber`, `lme4`, `lmerTest`

### 2단계: 서버 실행

**방법 A: Windows - 더블클릭**
```
run.bat 더블클릭
```

**방법 B: 터미널에서 직접 실행**
```bash
cd HLM_app
Rscript start_server.R
```

**방법 C: R 콘솔에서 실행**
```r
setwd("경로/HLM_app")
library(plumber)
pr <- plumb("api.R")
pr_run(pr, host = "127.0.0.1", port = 8001, docs = FALSE)
```

### 3단계: 브라우저 접속

```
http://localhost:8001
```

---

## API 엔드포인트

| Method | URL | 설명 |
|---|---|---|
| `GET` | `/api/health` | 서버 상태 확인 (R 버전, 패키지 버전) |
| `POST` | `/api/analyze` | 다층모형 분석 실행 |
| `GET` | `/*` | 정적 파일 제공 (`www/` 폴더) |

### POST /api/analyze 요청 형식
```json
{
  "data": { "col1": [1,2,3], "col2": [4,5,6] },
  "outcome": "종속변수명",
  "group_var": "집단변수명",
  "l1_preds": ["L1독립변수1", "L1독립변수2"],
  "l2_covs": ["L2공변량1"],
  "rand_slopes": ["무선기울기변수1"],
  "cross_interactions": [{"l1": "변수1", "l2": "변수2"}],
  "l1_centering": "none|gmc|cwc",
  "l2_centering": "none|gmc",
  "cov_struct": "un|ind"
}
```

### POST /api/analyze 응답 형식 (주요 필드)
```json
{
  "success": true,
  "data_info": { "n_obs": 7185, "n_groups": 160, "outcome": "MATHACH", ... },
  "null_model": {
    "icc": 0.1829,
    "fixed_effects": [{ "term": "(Intercept)", "estimate": 12.6370, "se": 0.2439, "t": 51.82, "p": 0.0000, "ci_lower": 12.1590, "ci_upper": 13.1150 }],
    "random_effects": [
      { "group": "school", "var1": "(Intercept)", "estimate": 8.6128, "se": 1.5820, "ci_lower": 6.0037, "ci_upper": 12.6017 },
      { "group": "Residual", "estimate": 39.1480, "se": 0.6632, "ci_lower": 37.8712, "ci_upper": 40.4655 }
    ],
    "r2": { "marginal": 0.0000, "conditional": 0.1804 },
    "fit": { "AIC": 47116.84, "BIC": 47137.50, "logLik": -23555.42, "deviance": 47110.84, "npar": 3 },
    "lrt_vs_lm": { "chi2": 1221.12, "df": 1, "p": 0.0000 }
  },
  "ri_model": { ... },
  "rs_model": { ... },
  "cl_model": { ... },
  "model_comparison": [ ... ],
  "diagnostics": { ... },
  "r_code": "# R 코드 로그..."
}
```

---

## 추정 방식 상세 설명

### EM 알고리즘 + REML/ML 이중 추정

이 앱은 Stata의 `mixed` 명령어와 동일한 추정 전략을 사용합니다:

```
[EM 알고리즘 (25회)]  →  [Newton-Raphson / BOBYQA]  →  최종 추정치
   (안정적 초기값)           (빠른 수렴)
```

1. **EM (Expectation-Maximization) 알고리즘**: 분산 성분(G, σ²)의 초기값을 안정적으로 추정
   - E-step: 각 집단의 조건부 무선효과 기대값과 분산 계산
   - M-step: 기대값으로부터 G(무선효과 공분산)와 σ²(잔차 분산) 갱신
2. **lmer() 최종 적합**: EM 초기값을 `start` 인자로 전달하여 빠르고 안정적으로 수렴

### REML vs ML

| | REML | ML |
|---|---|---|
| **용도** | 보고용 (분산 추정) | 모형 비교용 (LRT) |
| **특징** | 비편향 분산 추정치 | 고정효과가 다른 모형 간 LRT 가능 |
| **기본값** | HLM7, SPSS, Stata | - |

### 무선효과 SE 및 CI 계산

| 성분 | confint 스케일 | 변환 | SE 계산 |
|---|---|---|---|
| **분산** (var) | SD 스케일 [sd_lower, sd_upper] | CI: sd² | SE(σ²) ≈ 2σ × SE(σ) |
| **공분산** (cov) | 상관 스케일 [cor_lower, cor_upper] | CI: cor × sd₁ × sd₂ | SE(cov) ≈ SE(cor) × sd₁ × sd₂ |
| **잔차** (σ²) | SD 스케일 [sigma_lower, sigma_upper] | CI: sd² | SE(σ²) ≈ 2σ × SE(σ) |

---

## 다른 서버 배포 시 주의사항 및 예상 오류

### 1. R 환경 관련

#### `Rscript: command not found`
- **원인**: R이 설치되지 않았거나 PATH에 등록되지 않음
- **해결**:
  - Linux: `sudo apt install r-base` 또는 `sudo yum install R`
  - macOS: `brew install r`
  - Windows: R 설치 후 `C:\Program Files\R\R-x.x.x\bin\`을 시스템 PATH에 추가
  - Docker: `FROM r-base:4.x` 이미지 사용

#### 패키지 설치 실패 (`lme4`, `lmerTest`)
- **원인**: C/C++ 컴파일러 또는 시스템 라이브러리 누락
- **해결** (Linux):
  ```bash
  # Ubuntu/Debian
  sudo apt install build-essential gfortran liblapack-dev libblas-dev

  # CentOS/RHEL
  sudo yum install gcc gcc-c++ gcc-gfortran lapack-devel blas-devel
  ```
- **해결** (macOS):
  ```bash
  xcode-select --install
  ```

#### `lmerTest`가 없는 환경
- `lmerTest` 없이도 분석은 가능하지만 **p값이 계산되지 않습니다**
- 서버 시작 시 `[경고] lmerTest 없음` 메시지가 출력됩니다

### 2. 네트워크 / 포트 관련

#### `createTcpServer: address already in use`
- **원인**: 포트 8001이 이미 사용 중
- **해결**:
  ```bash
  # 사용 중인 프로세스 확인
  # Linux/macOS
  lsof -i :8001
  kill -9 <PID>

  # Windows
  netstat -ano | findstr :8001
  taskkill /F /PID <PID>
  ```
- **대안**: `start_server.R`에서 포트 번호를 변경하고, `www/index.html`의 `API_BASE` 변수도 수정
  ```javascript
  // index.html 상단 스크립트 영역
  const API_BASE = 'http://localhost:새포트번호';
  ```

#### CORS 오류 (프론트엔드와 백엔드가 다른 도메인/포트)
- `api.R`의 CORS 필터가 `Access-Control-Allow-Origin: *`으로 설정되어 있으므로 기본적으로 모든 출처 허용
- 프로덕션 환경에서는 특정 도메인만 허용하도록 수정 권장:
  ```r
  res$setHeader("Access-Control-Allow-Origin", "https://your-domain.com")
  ```

### 3. 프론트엔드 ↔ R 백엔드 연동 구조

```
┌───────────────────────┐      HTTP (JSON)      ┌─────────────────────┐
│   브라우저 (index.html) │ ◄──────────────────► │  R Plumber (api.R)  │
│                         │    POST /api/analyze  │                     │
│  - Chart.js 시각화      │    GET  /api/health   │  - em_lmer() (EM+NR)│
│  - PapaParse CSV 파싱   │                       │  - lme4::lmer()     │
│  - docx.js Word 생성    │                       │  - lmerTest (p값)   │
│  - AI API 호출          │                       │  - CORS 필터        │
└───────────────────────┘                        └─────────────────────┘
         │                                                │
         │  (프론트엔드에서 직접)                          │
         ▼                                                ▼
   OpenAI / Claude API                            R 통계 계산 엔진
```

- 프론트엔드(`index.html`)가 R 서버에서 정적으로 제공됨 (`@assets ./www /`)
- 분석 요청: 브라우저 → `POST /api/analyze` → R에서 lmer() 실행 → JSON 응답
- AI 보고서: 브라우저 → OpenAI/Claude API 직접 호출 (R 서버 경유하지 않음)
- AI API 키는 **프론트엔드에서 입력**하므로 서버에 저장되지 않음

### 4. JSON 직렬화 주의사항 (Plumber `auto_unbox`)

- Plumber의 `@serializer unboxedJSON` 옵션으로 인해 **길이가 1인 R 벡터가 스칼라로 변환**됩니다
- 예: `c("a")` → `"a"` (배열이 아닌 문자열)
- 프론트엔드의 `deepUnbox()` 함수와 `toArr()` 헬퍼가 이 문제를 처리합니다
- **새로운 API 필드를 추가할 때** 항상 `as.list()`로 감싸서 배열 형태를 보장하세요:
  ```r
  # 잘못된 예 (길이 1이면 스칼라로 변환됨)
  result$names <- some_vector

  # 올바른 예
  result$names <- as.list(some_vector)
  ```

### 5. 외부 접근 허용 (0.0.0.0 바인딩)

기본값은 `127.0.0.1` (로컬만 접근 가능). 외부 접근을 허용하려면:

```r
# start_server.R
pr_run(pr, host = "0.0.0.0", port = 8001, docs = FALSE)
```

그리고 방화벽에서 포트 8001을 허용해야 합니다:
```bash
# Linux (ufw)
sudo ufw allow 8001

# Linux (firewalld)
sudo firewall-cmd --add-port=8001/tcp --permanent
sudo firewall-cmd --reload
```

### 6. Docker 배포

```dockerfile
FROM r-base:4.4.0

RUN R -e "install.packages(c('plumber','lme4','lmerTest'), repos='https://cran.rstudio.com/')"

WORKDIR /app
COPY . /app

EXPOSE 8001

CMD ["Rscript", "-e", "library(plumber); pr <- plumb('api.R'); pr_run(pr, host='0.0.0.0', port=8001, docs=FALSE)"]
```

```bash
docker build -t hlm-app .
docker run -p 8001:8001 hlm-app
```

### 7. 리버스 프록시 (Nginx)

```nginx
server {
    listen 80;
    server_name your-domain.com;

    location / {
        proxy_pass http://127.0.0.1:8001;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_read_timeout 300s;  # 분석에 시간이 걸릴 수 있음
    }
}
```

### 8. 대용량 데이터 처리

- 기본 Plumber의 요청 크기 제한은 약 **100MB**
- 대용량 CSV 업로드 시 R 메모리 부족 가능:
  ```bash
  # R 메모리 제한 확인
  Rscript -e "gc(); cat(paste('Available:', round(as.numeric(system('free -m | grep Mem | awk \"{print $7}\"', intern=TRUE)), 0), 'MB\n'))"
  ```
- 10만 행 이상의 데이터에서 `lmer()` 수렴에 수 분 소요 가능
- 무선기울기 모형 (`m2`)이 특히 수렴이 어려움 → `bobyqa` 옵티마이저 + `maxfun = 200000` 적용됨

### 9. 흔한 통계 관련 경고 메시지

| 메시지 | 의미 | 대응 |
|---|---|---|
| `boundary (singular) fit` | 무선효과 분산이 0에 근접 (과모수화) | 무선기울기 제거 고려 |
| `Model failed to converge` | 최적화 알고리즘 수렴 실패 | 독립변수 축소, 센터링 적용 |
| `fixed-effect model matrix is rank deficient` | 다중공선성 문제 | 상관 높은 변수 제거 |

---

## 개발 히스토리

1. **기본 구조**: R Plumber API + HTML/JS 프론트엔드, 4단계 다층모형 분석
2. **기초 통계**: 기술통계, 분포(히스토그램), 상관행렬, 집단별 비교, 시각화
3. **교차수준 상호작용**: L1×L2 교차수준 상호작용 모형 지원
4. **진단 플롯**: 잔차 분석 6종 + 회귀선 시각화 2종 + Pairs Plot
5. **센터링**: GMC (Grand Mean Centering), CWC (Group Mean Centering) 옵션
6. **모형 공식 탭**: 수학적 표기(Level 1/2/Combined)로 모형 방정식 표시
7. **상관행렬 변수 선택**: 특정 행/열 포함·제외 체크박스
8. **AI 보고서 R 코드**: 분석에 사용된 실제 R 코드를 상세 주석과 함께 AI 보고서에 포함
9. **EM 알고리즘 + REML 이중 추정**: Stata mixed 방식 적용, HLM7/SPSS/Stata와 동일한 추정치
10. **Stata 스타일 무선효과 출력**: Estimate, Std. Err., 95% CI (profile likelihood)
11. **LR test vs. linear model**: Stata 스타일 무선효과 유의성 검정
12. **공분산 구조 선택**: Unstructured / Independent 옵션 추가
13. **공분산 행 표시**: Unstructured 선택 시 cov(var1, var2) 행 표시 (SE, CI 포함)

---

## 라이선스

MIT License
