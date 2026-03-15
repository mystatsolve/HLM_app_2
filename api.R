# ═══════════════════════════════════════════════════════════════════════════════
# api.R - 다층모형분석 (HLM/MLM) R Plumber REST API
# ═══════════════════════════════════════════════════════════════════════════════
#
# 이 파일은 다층모형(Hierarchical Linear Model / Multilevel Model) 분석을 위한
# R Plumber REST API 서버의 핵심 로직을 담고 있습니다.
#
# 주요 기능:
#   1. EM 알고리즘 + REML/ML 이중 추정 (Stata mixed 방식)
#   2. 4단계 모형 순차 적합 (M0 기저 → M1 무선절편 → M2 무선기울기 → M3 교차수준)
#   3. Stata 스타일 무선효과 출력 (Estimate, SE, 95% CI)
#   4. Profile likelihood 신뢰구간 및 delta method 표준오차
#   5. 공분산 구조 선택 (Unstructured / Independent)
#   6. LR test vs. linear model (Stata 스타일 무선효과 유의성 검정)
#   7. 센터링 옵션 (GMC, CWC)
#   8. 진단 플롯용 데이터 추출 (잔차, BLUP 등)
#
# 사용 패키지: plumber, lme4, lmerTest
# ═══════════════════════════════════════════════════════════════════════════════

library(plumber)
library(lme4)

# ── lmerTest 패키지 확인 ──────────────────────────────────────────────────────
# lmerTest는 Satterthwaite 근사 자유도로 고정효과의 p값을 계산합니다.
# lmerTest가 없으면 t값만 제공되고 p값은 NA가 됩니다.
# 또한 ranova() 함수로 개별 무선효과의 LRT 검정을 수행합니다.
has_lmerTest <- requireNamespace("lmerTest", quietly = TRUE)
if (has_lmerTest) {
  library(lmerTest)
  message("[OK] lmerTest 로드 - Satterthwaite p값 계산 가능")
} else {
  message("[경고] lmerTest 없음 - t값만 제공. install.R을 먼저 실행하세요.")
}

# ── 정적 파일 제공 ────────────────────────────────────────────────────────────
# www/ 폴더의 index.html 등을 루트 경로(/)에서 제공합니다.
# 브라우저에서 http://localhost:8001/ 접속 시 www/index.html이 로드됩니다.
#* @assets ./www /
list()

# ═══════════════════════════════════════════════════════════════════════════════
#  유틸리티 함수
# ═══════════════════════════════════════════════════════════════════════════════

# ── normalize_cross_interactions(): 교차상호작용 JSON 파싱 결과 정규화 ────────
#
# JSON 파싱 결과는 다양한 형태로 올 수 있으므로,
# 항상 list(list(l1="X", l2="Z"), ...) 형태로 정규화합니다.
#
# 가능한 입력:
#   data.frame, list-of-lists, flat named list, named character vector 등
#
# 반환: list of list(l1=character, l2=character)
normalize_cross_interactions <- function(x) {
  if (is.null(x) || length(x) == 0) return(list())

  # data.frame (jsonlite simplifyVector=TRUE로 배열 of 객체 파싱 시)
  if (is.data.frame(x)) {
    return(lapply(seq_len(nrow(x)), function(i) {
      list(l1 = as.character(x[["l1"]][i]),
           l2 = as.character(x[["l2"]][i]))
    }))
  }

  # 단일 상호작용이 flat named list로 온 경우: list(l1="A", l2="B")
  if (is.list(x) && all(c("l1", "l2") %in% names(x))) {
    return(list(list(
      l1 = as.character(x[["l1"]]),
      l2 = as.character(x[["l2"]])
    )))
  }

  # list of items (정상 케이스: list of lists, 또는 list of named vectors)
  if (is.list(x)) {
    result <- lapply(x, function(ci) {
      if (is.list(ci) && all(c("l1", "l2") %in% names(ci))) {
        list(l1 = as.character(ci[["l1"]]), l2 = as.character(ci[["l2"]]))
      } else if (is.character(ci) && all(c("l1", "l2") %in% names(ci))) {
        list(l1 = as.character(ci[["l1"]]), l2 = as.character(ci[["l2"]]))
      } else {
        NULL
      }
    })
    return(Filter(Negate(is.null), result))
  }

  # 그 외 (atomic vector 등): 안전하게 빈 리스트 반환
  list()
}

# ── em_lmer(): EM 알고리즘 초기값 + lmer 적합 ─────────────────────────────────
#
# Stata의 mixed 명령어와 동일한 2단계 추정 전략:
#   1단계: EM algorithm으로 분산 성분의 안정적인 초기값을 추정
#   2단계: EM 초기값을 lmer()에 전달하여 BOBYQA로 최종 수렴
#
# 장점: 복잡한 모형에서 수렴 안정성이 높음
# 단점: EM 반복으로 인해 계산 시간이 더 소요됨
#
# 인자:
#   formula  : lmer 공식 (예: Y ~ X + (1 + X | group))
#   data     : 데이터프레임
#   REML     : TRUE = REML 추정, FALSE = ML 추정
#   n_em     : EM 반복 횟수 (기본 25회)
#   control  : lmerControl 객체
em_lmer <- function(formula, data, REML = TRUE, n_em = 25,
                    control = lmerControl(), ...) {
  # lFormula()로 모형 행렬 계산, 실패 시 lmer() 직접 호출
  lmod <- tryCatch(lFormula(formula, data = data, REML = REML),
                   error = function(e) NULL)
  if (is.null(lmod)) {
    return(lmer(formula, data = data, REML = REML, control = control, ...))
  }

  # || (독립 공분산) 사용 시 EM 건너뛰기 (단일 G 행렬 가정 위배)
  formula_str <- paste(deparse(formula), collapse = "")
  if (grepl("\\|\\|", formula_str)) {
    return(lmer(formula, data = data, REML = REML, control = control, ...))
  }

  # 행렬 구성
  y <- as.numeric(lmod$fr[, 1])
  X <- as.matrix(lmod$X)
  Z <- as.matrix(t(lmod$reTrms$Zt))
  grp_idx <- as.integer(lmod$reTrms$flist[[1]])
  J <- nlevels(lmod$reTrms$flist[[1]])
  N <- length(y)
  p <- ncol(X)
  q <- length(lmod$reTrms$cnms[[1]])

  # EM 초기값
  sigma2 <- var(y) * 0.5
  G <- diag(q) * var(y) * 0.25

  for (iter in seq_len(n_em)) {
    # GLS로 beta 추정
    XtVX <- matrix(0, p, p)
    XtVY <- numeric(p)
    for (j in seq_len(J)) {
      idx  <- which(grp_idx == j)
      nj   <- length(idx)
      Xj   <- X[idx, , drop = FALSE]
      yj   <- y[idx]
      zcol <- ((j - 1) * q + 1):(j * q)
      Zj   <- Z[idx, zcol, drop = FALSE]
      Vj_inv <- solve(Zj %*% tcrossprod(G, Zj) + sigma2 * diag(nj))
      XtVX <- XtVX + crossprod(Xj, Vj_inv %*% Xj)
      XtVY <- XtVY + crossprod(Xj, Vj_inv %*% yj)
    }
    beta <- solve(XtVX, XtVY)

    # E-step + M-step
    ss_G <- matrix(0, q, q)
    ss_s <- 0
    for (j in seq_len(J)) {
      idx  <- which(grp_idx == j)
      nj   <- length(idx)
      Xj   <- X[idx, , drop = FALSE]
      yj   <- y[idx]
      zcol <- ((j - 1) * q + 1):(j * q)
      Zj   <- Z[idx, zcol, drop = FALSE]
      rj     <- yj - Xj %*% beta
      Vj_inv <- solve(Zj %*% tcrossprod(G, Zj) + sigma2 * diag(nj))
      GZt    <- tcrossprod(G, Zj)
      uj     <- GZt %*% (Vj_inv %*% rj)
      Cj     <- G - GZt %*% (Vj_inv %*% (Zj %*% G))
      ss_G <- ss_G + tcrossprod(uj) + Cj
      ej   <- rj - Zj %*% uj
      ss_s <- ss_s + sum(ej^2) + sum(Zj * (Zj %*% Cj))
    }

    G      <- ss_G / J
    sigma2 <- ss_s / N
    eig <- eigen(G, symmetric = TRUE)
    eig$values <- pmax(eig$values, 1e-10)
    G <- eig$vectors %*% diag(eig$values, nrow = q) %*% t(eig$vectors)
    sigma2 <- max(sigma2, 1e-10)
  }

  # EM 결과를 lmer theta 파라미터로 변환
  L <- tryCatch(t(chol(G)), error = function(e) {
    eig <- eigen(G, symmetric = TRUE)
    eig$values <- pmax(eig$values, 1e-10)
    t(chol(eig$vectors %*% diag(eig$values, nrow = q) %*% t(eig$vectors)))
  })
  theta_vec <- (L / sqrt(sigma2))[lower.tri(L, diag = TRUE)]

  lmer(formula, data = data, REML = REML,
       start = list(theta = theta_vec),
       control = control, ...)
}

# ── fit_model(): 추정 방법에 따라 lmer 또는 em_lmer 호출 ─────────────────────
#
# estimation_method 파라미터에 따라 적합 함수를 선택합니다:
#   "lmer"     : lme4::lmer() 직접 호출 (빠르고 표준적)
#   "em_lmer"  : EM 알고리즘 초기값 + lmer() (Stata 방식, 수렴 안정성 높음)
fit_model <- function(formula, data, REML, control, method = "lmer") {
  if (method == "em_lmer") {
    em_lmer(formula, data = data, REML = REML, control = control)
  } else {
    lmer(formula, data = data, REML = REML, control = control)
  }
}

# ── calc_lrt_vs_lm(): LR test vs. linear model (Stata 스타일) ─────────────────
#
# 혼합모형(random effects 포함)이 단순 선형모형(OLS)보다 유의하게 적합한지를
# 우도비 검정(Likelihood Ratio Test)으로 검정합니다.
#
# Stata의 "LR test vs. linear model: chi2(k) = ..., Prob > chi2 = ..." 출력과 동일합니다.
# 이 검정이 유의하면 → "집단 간 분산이 존재하므로 다층모형이 필요하다"는 의미입니다.
#
# 인자:
#   mixed_ml     : ML로 적합한 혼합모형 (REML=FALSE)
#   formula_fixed: 고정효과만 포함한 공식 (예: Y ~ X1 + X2)
#   data         : 데이터프레임
#
# 반환: list(chi2, df, p) 또는 오류 시 NULL
calc_lrt_vs_lm <- function(mixed_ml, formula_fixed, data) {
  tryCatch({
    m_lm    <- lm(formula_fixed, data = data)         # 단순 OLS 적합
    ll_lm   <- as.numeric(logLik(m_lm))               # OLS 로그우도
    ll_mix  <- as.numeric(logLik(mixed_ml))            # 혼합모형 로그우도
    chi2    <- -2 * (ll_lm - ll_mix)                   # LRT 통계량 = -2 × ΔLL
    df_test <- attr(logLik(mixed_ml), "df") - attr(logLik(m_lm), "df")  # 자유도 차이
    list(
      chi2 = round(max(chi2, 0), 4),
      df   = as.integer(df_test),
      p    = signif(pchisq(max(chi2, 0), df_test, lower.tail = FALSE), 4)
    )
  }, error = function(e) NULL)
}

# ── calc_icc(): ICC (집단내 상관계수) 계산 ────────────────────────────────────
#
# ICC = τ₀₀ / (τ₀₀ + σ²)
# τ₀₀: 집단 간 절편 분산 (between-group variance)
# σ² : 집단 내 잔차 분산 (within-group variance)
#
# ICC가 높을수록 종속변수의 분산 중 집단 간 차이가 차지하는 비율이 큽니다.
# 일반적으로 ICC > 0.05이면 다층모형 사용을 고려합니다.
calc_icc <- function(model) {
  tryCatch({
    vc <- as.data.frame(VarCorr(model))
    int_var <- vc$vcov[!is.na(vc$var1) & vc$var1 == "(Intercept)" & is.na(vc$var2)]
    res_var <- vc$vcov[vc$grp == "Residual"]
    if (length(int_var) == 0 || length(res_var) == 0) return(NA)
    round(int_var[1] / (int_var[1] + res_var[1]), 4)
  }, error = function(e) NA)
}

# ── calc_r2(): R² 계산 (Nakagawa & Schielzeth, 2013) ──────────────────────────
#
# R² marginal:    고정효과만으로 설명되는 분산 비율
#                 = Var(Xβ) / (Var(Xβ) + Var(무선효과) + σ²)
# R² conditional: 고정효과 + 무선효과로 설명되는 분산 비율
#                 = (Var(Xβ) + Var(무선효과)) / (Var(Xβ) + Var(무선효과) + σ²)
#
# 참고: Nakagawa, S., & Schielzeth, H. (2013). A general and simple method for
#       obtaining R² from generalized linear mixed-effects models.
#       Methods in Ecology and Evolution, 4(2), 133-142.
calc_r2 <- function(model) {
  tryCatch({
    vc    <- as.data.frame(VarCorr(model))
    sig2  <- sigma(model)^2
    X     <- model.matrix(model)
    beta  <- fixef(model)
    var_f <- var(as.vector(X %*% beta))
    var_r <- sum(vc$vcov[vc$grp != "Residual" & is.na(vc$var2)])
    tot   <- var_f + var_r + sig2
    list(
      marginal    = round(var_f / tot, 4),
      conditional = round((var_f + var_r) / tot, 4)
    )
  }, error = function(e) list(marginal = NA, conditional = NA))
}

# ── extract_fit(): 모형 적합도 지수 추출 ───────────────────────────────────────
# AIC, BIC: 작을수록 적합도 우수 (모형 비교에 사용)
# logLik: 로그우도 (클수록 적합도 우수)
# deviance: -2 × logLik (작을수록 적합도 우수, LRT에 사용)
# npar: 추정된 파라미터 수
extract_fit <- function(model) {
  list(
    AIC      = round(AIC(model), 2),
    BIC      = round(BIC(model), 2),
    logLik   = round(as.numeric(logLik(model)), 2),
    deviance = round(-2 * as.numeric(logLik(model)), 2),
    npar     = attr(logLik(model), "df")
  )
}

# ── extract_fixed(): 고정효과 추출 ─────────────────────────────────────────────
#
# 각 고정효과(절편, 독립변수)에 대해 다음을 추출합니다:
#   - estimate: 추정값 (γ 계수)
#   - se: 표준오차
#   - t: t통계량
#   - p: p값 (lmerTest의 Satterthwaite 자유도 기반)
#   - ci_lower, ci_upper: Wald 95% 신뢰구간
#
# 참고: summary()$coefficients 대신 coef(summary())를 사용합니다.
#       lmerTest가 로드되면 lmer 모형은 S4 클래스(lmerModLmerTest)가 되는데,
#       S4 객체에 $를 사용하면 오류가 발생하기 때문입니다.
extract_fixed <- function(model) {
  cm  <- as.data.frame(coef(summary(model)))
  cm$term <- rownames(cm)
  rownames(cm) <- NULL

  ci <- tryCatch(
    as.data.frame(confint(model, method = "Wald", level = 0.95)),
    error = function(e) NULL
  )

  lapply(seq_len(nrow(cm)), function(i) {
    row  <- cm[i, ]
    pcol <- grep("Pr\\(>\\|t\\|\\)", names(cm), value = TRUE)
    pval <- if (length(pcol) > 0) round(row[[pcol[1]]], 4) else NA

    ci_l <- NA; ci_u <- NA
    if (!is.null(ci) && row$term %in% rownames(ci)) {
      ci_l <- round(ci[row$term, 1], 4)
      ci_u <- round(ci[row$term, 2], 4)
    }
    list(
      term     = row$term,
      estimate = round(row$Estimate, 4),
      se       = round(row[["Std. Error"]], 4),
      t        = round(row[["t value"]], 4),
      p        = pval,
      ci_lower = ci_l,
      ci_upper = ci_u
    )
  })
}

# ── extract_random(): 무선효과 추출 (Stata 스타일) ─────────────────────────────
#
# Stata의 mixed 명령어 출력과 동일한 형식으로 무선효과를 추출합니다:
#   - estimate: 분산 추정치 (variance, σ² 또는 τ²)
#   - se: 표준오차 (delta method로 계산)
#   - ci_lower, ci_upper: 95% 신뢰구간 (Wald 근사)
#
# ── Wald 근사 CI 계산 방법 ──
#
# Profile likelihood CI는 가장 정확하지만 수분~수십분 소요되어
# 단일 스레드 Plumber 서버를 완전히 블로킹합니다.
# 대신 Wald 근사 CI를 사용합니다 (SPSS, HLM7과 동일 방식):
#
# 1) SD 파라미터: SE(sd) ≈ sd / sqrt(2 × n_groups)
#    CI = [max(0, sd - 1.96×SE), sd + 1.96×SE]
#    → 분산 CI: ci² (SD CI를 제곱)
#
# 2) 상관(correlation): Fisher z-transform
#    z = 0.5 × ln((1+r)/(1-r)), SE(z) = 1/sqrt(n_groups - 3)
#    CI(z) = z ± 1.96×SE(z), 역변환하여 CI(r) 계산
#    → 공분산 CI: r × sd₁ × sd₂
#
# 3) 잔차: SE(sigma) ≈ sigma / sqrt(2 × (n - p))
#
# 인자:
#   model: lmerMod 객체
#   data:  원본 데이터프레임 (ranova 실행에 필요)
#
# 반환: list of list (각 무선효과 성분별)
extract_random <- function(model, data = NULL) {
  vc <- as.data.frame(VarCorr(model))

  # ── Wald 근사 신뢰구간 계산 (즉시 완료, 서버 블로킹 없음) ──
  ci_prof <- tryCatch({
    ngrps_vec <- ngrps(model)
    n_obs <- nobs(model)
    n_fe  <- length(fixef(model))

    ci_rows <- list()
    for (i in seq_len(nrow(vc))) {
      grp  <- vc$grp[i]
      v1   <- as.character(vc$var1[i])
      v2   <- if (!is.na(vc$var2[i])) as.character(vc$var2[i]) else NA
      is_cov   <- !is.na(v2)
      is_resid <- grp == "Residual"
      sdcor    <- vc$sdcor[i]

      if (is_resid) {
        # 잔차 SD: chi-squared 근사
        df_r  <- max(n_obs - n_fe, 1)
        se_sd <- sdcor / sqrt(2 * df_r)
        ci_rows[["sigma"]] <- c(max(0, sdcor - 1.96 * se_sd), sdcor + 1.96 * se_sd)
      } else if (is_cov) {
        # 상관: Fisher z-transform CI
        r <- sdcor  # 공분산 행에서 sdcor = correlation
        m <- ngrps_vec[grp]
        if (!is.na(m) && m > 3 && abs(r) < 0.9999) {
          z     <- 0.5 * log((1 + r) / (1 - r))
          se_z  <- 1 / sqrt(max(m - 3, 1))
          z_lo  <- z - 1.96 * se_z
          z_hi  <- z + 1.96 * se_z
          r_lo  <- tanh(z_lo)
          r_hi  <- tanh(z_hi)
        } else {
          r_lo <- -1; r_hi <- 1
        }
        rn <- paste0("cor_", v1, ".", v2, "|", grp)
        ci_rows[[rn]] <- c(r_lo, r_hi)
      } else {
        # SD 파라미터: 점근 정규 근사
        m     <- ngrps_vec[grp]
        se_sd <- sdcor / sqrt(2 * max(m, 2))
        rn <- if (v1 == "(Intercept)") {
          paste0("sd_(Intercept)|", grp)
        } else {
          paste0("sd_", v1, "|", grp)
        }
        ci_rows[[rn]] <- c(max(0, sdcor - 1.96 * se_sd), sdcor + 1.96 * se_sd)
      }
    }

    ci_df <- do.call(rbind, ci_rows)
    colnames(ci_df) <- c("2.5 %", "97.5 %")
    as.data.frame(ci_df)
  }, error = function(e) {
    message("[Wald CI] 오류: ", e$message)
    NULL
  })

  # ── ranova: 개별 무선효과의 LRT 검정 ──
  # ranova()는 각 무선효과를 하나씩 제거한 축소 모형과 비교하여
  # 해당 무선효과가 유의한지 검정합니다.
  #
  # 주의: Plumber API 환경에서 ranova()는 내부적으로 update()를 호출하는데,
  # update()가 모형 호출에 사용된 데이터/control 변수명을 찾지 못합니다.
  # 해결: 데이터와 control을 전역 환경에 임시로 할당하고, on.exit()으로 정리합니다.
  ranova_df <- NULL
  if (has_lmerTest && !is.null(data)) {
    tryCatch({
      data_name <- deparse(getCall(model)$data)
      assign(data_name, data, envir = .GlobalEnv)
      on.exit(try(rm(list = data_name, envir = .GlobalEnv), silent = TRUE), add = TRUE)

      # control 변수도 전역 환경에 할당 (ranova의 update()가 참조)
      ctrl_sym <- getCall(model)$control
      if (!is.null(ctrl_sym)) {
        ctrl_name <- deparse(ctrl_sym)
        if (!exists(ctrl_name, envir = .GlobalEnv)) {
          assign(ctrl_name, lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)),
                 envir = .GlobalEnv)
          on.exit(try(rm(list = ctrl_name, envir = .GlobalEnv), silent = TRUE), add = TRUE)
        }
      }

      ra <- lmerTest::ranova(model)
      ranova_df <- as.data.frame(ra)
      ranova_df$rowname <- rownames(ra)
    }, error = function(e) {
      message("[ranova] 오류: ", e$message)
    })
  }

  lapply(seq_len(nrow(vc)), function(i) {
    is_cov   <- !is.na(vc$var2[i])
    is_resid <- vc$grp[i] == "Residual"
    var_name <- as.character(vc$var1[i])

    # --- Stata 스타일 SE 및 CI 계산 ---
    estimate <- vc$vcov[i]   # 분산 (Estimate)
    sd_val   <- vc$sdcor[i]  # 표준편차
    se_est   <- NA
    ci_lower <- NA
    ci_upper <- NA

    if (!is.null(ci_prof)) {
      ci_names <- rownames(ci_prof)

      if (is_resid) {
        # 잔차: "sigma" 파라미터
        idx <- which(ci_names == "sigma")
        if (length(idx) > 0) {
          sd_lower <- max(ci_prof[idx[1], 1], 0)
          sd_upper <- ci_prof[idx[1], 2]
          ci_lower <- sd_lower^2
          ci_upper <- sd_upper^2
          se_sd <- (sd_upper - sd_lower) / (2 * 1.96)
          se_est <- 2 * sd_val * se_sd
        }
      } else if (is_cov) {
        # 공분산: "cor_var1.var2|grp" 형태 → 상관 스케일 CI
        grp_name <- vc$grp[i]
        var2_name <- as.character(vc$var2[i])
        matched_ci <- NULL
        for (k in seq_along(ci_names)) {
          cn <- ci_names[k]
          if (grepl("^cor_", cn) && grepl(grp_name, cn, fixed = TRUE)) {
            v1_match <- if (var_name == "(Intercept)") grepl("Intercept", cn, fixed = TRUE) else grepl(var_name, cn, fixed = TRUE)
            v2_match <- if (var2_name == "(Intercept)") grepl("Intercept", cn, fixed = TRUE) else grepl(var2_name, cn, fixed = TRUE)
            if (v1_match && v2_match) { matched_ci <- k; break }
          }
        }
        if (!is.null(matched_ci)) {
          cor_lower <- ci_prof[matched_ci, 1]
          cor_upper <- ci_prof[matched_ci, 2]
          sd1 <- NA; sd2 <- NA
          for (r in seq_len(nrow(vc))) {
            if (vc$grp[r] == grp_name && is.na(vc$var2[r])) {
              if (as.character(vc$var1[r]) == var_name) sd1 <- vc$sdcor[r]
              if (as.character(vc$var1[r]) == var2_name) sd2 <- vc$sdcor[r]
            }
          }
          if (!is.na(sd1) && !is.na(sd2)) {
            ci_lower <- cor_lower * sd1 * sd2
            ci_upper <- cor_upper * sd1 * sd2
            se_cor <- (cor_upper - cor_lower) / (2 * 1.96)
            se_est <- se_cor * sd1 * sd2
          }
        }
      } else {
        # 분산: "sd_(Intercept)|grp" 또는 "sd_varname|grp" 형태
        grp_name <- vc$grp[i]
        matched_ci <- NULL
        # exact match로 confint rowname 찾기
        if (var_name == "(Intercept)") {
          exact_name <- paste0("sd_(Intercept)|", grp_name)
        } else {
          exact_name <- paste0("sd_", var_name, "|", grp_name)
        }
        idx <- which(ci_names == exact_name)
        if (length(idx) > 0) matched_ci <- idx[1]

        # 매칭 안 되면 부분 매칭 시도
        if (is.null(matched_ci)) {
          for (k in seq_along(ci_names)) {
            cn <- ci_names[k]
            if (grepl("^sd_", cn) && grepl(grp_name, cn, fixed = TRUE)) {
              if (var_name == "(Intercept)" && grepl("Intercept", cn, fixed = TRUE)) {
                matched_ci <- k; break
              } else if (var_name != "(Intercept)" && grepl(var_name, cn, fixed = TRUE)) {
                matched_ci <- k; break
              }
            }
          }
        }

        if (!is.null(matched_ci)) {
          sd_lower <- max(ci_prof[matched_ci, 1], 0)
          sd_upper <- ci_prof[matched_ci, 2]
          ci_lower <- sd_lower^2
          ci_upper <- sd_upper^2
          se_sd <- (sd_upper - sd_lower) / (2 * 1.96)
          se_est <- 2 * sd_val * se_sd
        }
      }
    }

    # --- ranova LRT ---
    lrt_df <- NA; lrt_chisq <- NA; lrt_p <- NA

    if (!is.null(ranova_df) && !is_cov && !is_resid) {
      for (j in seq_len(nrow(ranova_df))) {
        rn <- ranova_df$rowname[j]
        if (rn == "<none>") next
        matched <- FALSE
        if (var_name == "(Intercept)") {
          if (grepl("(1 | ", rn, fixed = TRUE) && !grepl(" in ", rn, fixed = TRUE)) matched <- TRUE
        } else {
          if (grepl(var_name, rn, fixed = TRUE)) matched <- TRUE
        }
        if (matched) {
          lrt_chisq <- ranova_df[j, "LRT"]
          lrt_df    <- ranova_df[j, "Df"]
          lrt_p     <- ranova_df[j, "Pr(>Chisq)"]
          break
        }
      }
    }

    list(
      group     = vc$grp[i],
      var1      = ifelse(is.na(vc$var1[i]), "", var_name),
      var2      = ifelse(is.na(vc$var2[i]), "", as.character(vc$var2[i])),
      estimate  = round(estimate, 4),
      se        = if (!is.na(se_est)) round(se_est, 4) else NA,
      ci_lower  = if (!is.na(ci_lower)) round(ci_lower, 4) else NA,
      ci_upper  = if (!is.na(ci_upper)) round(ci_upper, 4) else NA,
      vcov      = round(vc$vcov[i], 4),
      sdcor     = round(vc$sdcor[i], 4),
      lrt_df    = if (!is.na(lrt_df)) as.integer(lrt_df) else NA,
      lrt_chisq = if (!is.na(lrt_chisq)) round(lrt_chisq, 4) else NA,
      lrt_p     = if (!is.na(lrt_p)) signif(lrt_p, 4) else NA
    )
  })
}

# ── extract_diagnostics(): 모형 진단 데이터 추출 ──────────────────────────────
#
# 프론트엔드의 진단 플롯 (잔차 히스토그램, Q-Q plot, 잔차 vs 예측값 등)에
# 필요한 데이터를 추출합니다.
#
# 추출 항목:
#   - predicted: 모형 예측값
#   - residuals: 잔차 (관측값 - 예측값)
#   - zresid: 표준화 잔차 ((잔차 - 평균) / SD)
#   - re_data: 무선효과 BLUP (Best Linear Unbiased Prediction)
#   - predictor_values: 독립변수 값 (잔차 vs 독립변수 플롯용)
extract_diagnostics <- function(model, df, outcome, group_var, l1_preds, l2_covs, model_name) {
  tryCatch({
    predicted <- as.numeric(predict(model))
    obs       <- as.numeric(df[[outcome]])
    residuals <- obs - predicted
    res_mean  <- mean(residuals, na.rm = TRUE)
    res_sd    <- sd(residuals, na.rm = TRUE)
    zresid    <- if (res_sd > 0) (residuals - res_mean) / res_sd else rep(0, length(residuals))

    # 무선효과 BLUP
    re_list   <- ranef(model)
    re_data   <- list()
    re_names  <- character(0)
    if (group_var %in% names(re_list)) {
      re_raw   <- re_list[[group_var]]
      re_names <- colnames(re_raw)
      re_data  <- lapply(re_names, function(nm) {
        list(name = nm, values = round(as.numeric(re_raw[[nm]]), 6))
      })
    }

    # 독립변수 값 (잔차 vs 독립변수 플롯용)
    all_preds <- c(l1_preds, l2_covs)
    predictor_values <- list()
    if (length(all_preds) > 0) {
      for (v in all_preds) {
        if (v %in% names(df)) {
          predictor_values[[v]] <- as.numeric(df[[v]])
        }
      }
    }

    # re_names를 항상 배열로 보내기 위해 as.list()
    list(
      model_name       = model_name,
      predicted        = round(predicted, 6),
      residuals        = round(residuals, 6),
      zresid           = round(zresid, 6),
      re_data          = re_data,
      re_names         = as.list(re_names),
      predictor_names  = as.list(all_preds),
      predictor_values = predictor_values,
      outcome_values   = round(obs, 6),
      outcome_name     = outcome,
      group_values     = as.character(df[[group_var]]),
      group_name       = group_var
    )
  }, error = function(e) {
    list(model_name = model_name, error = e$message)
  })
}

# ═══════════════════════════════════════════════════════════════════════════════
#  Plumber 엔드포인트 (REST API)
# ═══════════════════════════════════════════════════════════════════════════════

# ── CORS 필터 ─────────────────────────────────────────────────────────────────
# 프론트엔드(브라우저)에서 API를 호출할 수 있도록 CORS 헤더를 설정합니다.
# 모든 출처(*)에서의 요청을 허용합니다.
# 프로덕션 환경에서는 특정 도메인만 허용하도록 수정하세요.
#* CORS 필터
#* @filter cors
function(req, res) {
  res$setHeader("Access-Control-Allow-Origin",  "*")
  res$setHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
  res$setHeader("Access-Control-Allow-Headers", "Content-Type, Accept")
  if (req$REQUEST_METHOD == "OPTIONS") { res$status <- 200; return(list()) }
  plumber::forward()
}

# ── /api/health: 서버 상태 확인 엔드포인트 ─────────────────────────────────────
# R 버전, lme4 버전, lmerTest 설치 여부를 반환합니다.
#* API 상태 확인
#* @get /api/health
#* @serializer unboxedJSON
function() {
  list(
    status   = "ok",
    R        = as.character(getRversion()),
    lme4     = as.character(packageVersion("lme4")),
    lmerTest = if (has_lmerTest) as.character(packageVersion("lmerTest")) else "미설치"
  )
}

# ── /api/analyze: 다층모형 분석 실행 (핵심 엔드포인트) ─────────────────────────
#
# 프론트엔드에서 CSV 데이터와 분석 설정을 JSON으로 전송하면,
# 이 엔드포인트에서 4단계 모형(M0~M3)을 순차적으로 적합하고 결과를 반환합니다.
#
# 요청 본문 (JSON):
#   data: { col1: [1,2,3], col2: [4,5,6] }  - 컬럼 형식 데이터
#   outcome: "종속변수명"
#   group_var: "집단변수명"
#   l1_preds: ["L1독립변수1", ...]
#   l2_covs: ["L2공변량1", ...]
#   rand_slopes: ["무선기울기변수1", ...]
#   cross_interactions: [{ l1: "변수1", l2: "변수2" }, ...]
#   l1_centering: "none" | "gmc" | "cwc"
#   l2_centering: "none" | "gmc"
#   cov_struct: "un" (비구조적) | "ind" (독립)
#
# 추정 방식:
#   - REML (보고용): 분산 성분의 비편향 추정치를 제공 (HLM7, SPSS, Stata 기본값)
#   - ML (모형비교용): 중첩 모형 간 LRT 비교에 사용 (REML은 고정효과가 다른 모형 비교 불가)
#   - 각 모형마다 REML + ML 두 번 적합합니다.
#* 다층모형 분석 실행
#* @post /api/analyze
#* @parser json
#* @serializer unboxedJSON
function(req) {
  tryCatch({
    body <- req$body

    # ── 1. 데이터 변환 (컬럼 형식 JSON → data.frame) ──────────────
    df <- as.data.frame(body$data, stringsAsFactors = FALSE)

    # ── 2. 파라미터 ────────────────────────────────────────────────
    outcome     <- as.character(body$outcome)
    group_var   <- as.character(body$group_var)
    l1_preds    <- if (!is.null(body$l1_preds))    as.character(unlist(body$l1_preds))    else character(0)
    l2_covs     <- if (!is.null(body$l2_covs))     as.character(unlist(body$l2_covs))     else character(0)
    rand_slopes <- if (!is.null(body$rand_slopes)) as.character(unlist(body$rand_slopes)) else character(0)
    l1_centering <- if (!is.null(body$l1_centering)) as.character(body$l1_centering) else "none"
    l2_centering <- if (!is.null(body$l2_centering)) as.character(body$l2_centering) else "none"
    # 공분산 구조: "un" = Unstructured (분산+공분산 모두 추정, lme4: | )
    #              "ind" = Independent (분산만 추정, 공분산=0, lme4: || )
    cov_struct   <- if (!is.null(body$cov_struct)) as.character(body$cov_struct) else "un"
    # 추정 방법: "lmer" = lme4 기본 (빠름), "em_lmer" = EM 초기값 + lmer (Stata 방식)
    est_method   <- if (!is.null(body$estimation_method)) as.character(body$estimation_method) else "lmer"
    l2_covs_original <- l2_covs   # CWC 처리 전 원본 보존 (rcode 로그용)
    cross_interactions <- normalize_cross_interactions(body$cross_interactions)

    # ── 3. 변수 존재 확인 ─────────────────────────────────────────
    needed  <- unique(c(outcome, group_var, l1_preds, l2_covs))
    missing <- setdiff(needed, names(df))
    if (length(missing) > 0) stop(paste0("데이터에 없는 변수: ", paste(missing, collapse = ", ")))

    # ── 4. 수치형 변환 ────────────────────────────────────────────
    df[[outcome]] <- suppressWarnings(as.numeric(df[[outcome]]))
    for (v in c(l1_preds, l2_covs)) {
      df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
    }

    # ── 5. 결측치 제거 ────────────────────────────────────────────
    df       <- df[complete.cases(df[needed]), ]
    n_obs    <- nrow(df)
    n_groups <- length(unique(df[[group_var]]))

    if (n_obs    <  10) stop("관측치가 너무 적습니다 (최소 10개 이상 필요).")
    if (n_groups <   3) stop("집단 수가 너무 적습니다 (최소 3개 이상 필요).")

    # ── 6. 센터링 (Centering) ──────────────────────────────────────
    # 센터링은 다층모형에서 계수의 해석을 용이하게 합니다.
    #
    # GMC (Grand Mean Centering): X_gmc = X - mean(X)
    #   - 절편(γ₀₀)이 "전체 평균에서의 종속변수 기대값"을 의미
    #   - L1, L2 변수 모두 적용 가능
    #
    # CWC (Centering Within Context = Group Mean Centering):
    #   - X_cwc = X - mean_j(X)  (각 관측치에서 해당 집단 평균을 뺌)
    #   - 순수한 개인 내 효과(within-group effect)를 분리
    #
    cwc_gm_vars <- character(0)

    if (l1_centering == "gmc" && length(l1_preds) > 0) {
      for (v in l1_preds) {
        df[[v]] <- df[[v]] - mean(df[[v]], na.rm = TRUE)
      }
    } else if (l1_centering == "cwc" && length(l1_preds) > 0) {
      for (v in l1_preds) {
        gm <- ave(df[[v]], df[[group_var]], FUN = function(x) mean(x, na.rm = TRUE))
        df[[v]] <- df[[v]] - gm                          # 집단평균 중심화 (within 효과)
      }
    }

    if (l2_centering == "gmc" && length(l2_covs) > 0) {
      for (v in l2_covs) {
        if (v %in% names(df)) {
          df[[v]] <- df[[v]] - mean(df[[v]], na.rm = TRUE)
        }
      }
    }

    result <- list(
      success   = TRUE,
      data_info = list(
        n_obs    = n_obs,
        n_groups = n_groups,
        outcome  = outcome,
        group_var = group_var,
        l1_centering = l1_centering,
        l2_centering = l2_centering,
        cov_struct   = cov_struct,
        cwc_gm_vars  = as.list(cwc_gm_vars),
        estimation_method = est_method
      )
    )

    # ── R 코드 로그 구축 ────────────────────────────────────────
    rcode <- character(0)
    rcode <- c(rcode,
      "# ═══════════════════════════════════════════════════",
      "# 다층모형(MLM/HLM) 분석 R 코드",
      "# 패키지: lme4, lmerTest",
      "# ═══════════════════════════════════════════════════",
      "",
      "# ── 1. 패키지 로드 ──",
      "library(lme4)",
      "library(lmerTest)  # Satterthwaite 자유도 기반 p값 계산",
      "",
      "# ── 2. 데이터 불러오기 및 전처리 ──",
      paste0('df <- read.csv("데이터파일.csv")'),
      paste0('# 종속변수: ', outcome),
      paste0('# 집단변수(Level-2 ID): ', group_var),
      if (length(l1_preds) > 0) paste0('# Level-1 독립변수: ', paste(l1_preds, collapse = ', ')) else '# Level-1 독립변수: (없음)',
      if (length(l2_covs_original) > 0) paste0('# Level-2 공변량: ', paste(l2_covs_original, collapse = ', ')) else '# Level-2 공변량: (없음)',
      paste0('# 관측치: ', n_obs, '개, 집단 수: ', n_groups, '개'),
      ""
    )

    # 센터링 코드 로그
    if (l1_centering == "gmc" && length(l1_preds) > 0) {
      rcode <- c(rcode,
        "# ── 3. 센터링: Grand Mean Centering (L1) ──",
        "# 각 L1 독립변수에서 전체 평균을 빼서 중심화",
        paste0("# X_gmc = X - mean(X)")
      )
      for (v in l1_preds) {
        rcode <- c(rcode, paste0('df$', v, ' <- df$', v, ' - mean(df$', v, ', na.rm = TRUE)'))
      }
      rcode <- c(rcode, "")
    } else if (l1_centering == "cwc" && length(l1_preds) > 0) {
      rcode <- c(rcode,
        "# ── 3. 센터링: Group Mean Centering (CWC, L1) ──",
        "# 각 L1 독립변수에서 해당 집단의 평균을 빼서 중심화",
        paste0("# X_cwc = X - mean_j(X)")
      )
      for (v in l1_preds) {
        rcode <- c(rcode,
          paste0('gm_', v, ' <- ave(df$', v, ', df$', group_var, ', FUN = function(x) mean(x, na.rm = TRUE))'),
          paste0('df$', v, ' <- df$', v, ' - gm_', v, '  # 집단평균 중심화 (CWC)')
        )
      }
      rcode <- c(rcode, "")
    } else {
      rcode <- c(rcode, "# ── 3. 센터링: 없음 (원점수 사용) ──", "")
    }

    if (l2_centering == "gmc" && length(l2_covs) > 0) {
      rcode <- c(rcode, "# L2 공변량 Grand Mean Centering")
      for (v in l2_covs) {
        rcode <- c(rcode, paste0('df$', v, ' <- df$', v, ' - mean(df$', v, ', na.rm = TRUE)'))
      }
      rcode <- c(rcode, "")
    }

    # ══════════════════════════════════════════════════════════════
    # Model 0: 기저모형 (Null Model / Empty Model)
    # ══════════════════════════════════════════════════════════════
    # 고정효과: 절편만 포함 (Y ~ 1)
    # 무선효과: 절편의 집단 간 변동만 허용 ((1|group))
    # 용도: ICC 계산 → 종속변수 분산 중 집단 간 차이 비율 확인
    #       ICC가 높으면 다층모형이 필요하다는 근거
    f0_str   <- paste0(outcome, " ~ 1 + (1 | ", group_var, ")")
    f0       <- as.formula(f0_str)
    f0_fixed <- as.formula(paste0(outcome, " ~ 1"))
    lmer_ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    m0    <- suppressWarnings(fit_model(f0, data = df, REML = TRUE, control = lmer_ctrl, method = est_method))
    m0_ml <- suppressWarnings(fit_model(f0, data = df, REML = FALSE, control = lmer_ctrl, method = est_method))

    rcode <- c(rcode,
      "# ── 4. Model 0: 기저모형 (Null Model) ──",
      "# 고정효과: 절편만, 무선효과: 절편의 집단 간 변동",
      paste0('m0 <- lmer(', f0_str, ','),
      "             data = df, REML = TRUE,",
      "             control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5)))",
      "summary(m0)",
      paste0("# ICC = ", round(calc_icc(m0), 4), " → 종속변수 분산의 ", round(calc_icc(m0) * 100, 1), "%가 집단 간 차이"),
      "",
      "# ML 추정 (모형 비교용 LRT에 필요)",
      paste0('m0_ml <- lmer(', f0_str, ','),
      "               data = df, REML = FALSE,",
      "               control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5)))",
      ""
    )

    result$null_model <- list(
      formula        = f0_str,
      icc            = calc_icc(m0),
      fixed_effects  = extract_fixed(m0),
      random_effects = extract_random(m0, df),
      r2             = calc_r2(m0),
      fit            = extract_fit(m0_ml),
      lrt_vs_lm      = calc_lrt_vs_lm(m0_ml, f0_fixed, df)
    )

    # ══════════════════════════════════════════════════════════════
    # Model 1: 무선절편 모형 (Random Intercept Model)
    # ══════════════════════════════════════════════════════════════
    # 고정효과: 절편 + L1 독립변수 + L2 공변량
    # 무선효과: 절편만 집단 간 변동 허용 (기울기는 고정)
    # 의미: 독립변수의 효과(기울기)는 모든 집단에서 동일하지만,
    #       종속변수의 평균(절편)은 집단마다 다를 수 있음
    all_fixed <- c(l1_preds, l2_covs)
    if (length(all_fixed) > 0) {
      f1_str <- paste0(outcome, " ~ ",
                       paste(all_fixed, collapse = " + "),
                       " + (1 | ", group_var, ")")
      f1_fixed <- as.formula(paste0(outcome, " ~ ", paste(all_fixed, collapse = " + ")))
      m1    <- suppressWarnings(fit_model(as.formula(f1_str), data = df, REML = TRUE, control = lmer_ctrl, method = est_method))
      m1_ml <- suppressWarnings(fit_model(as.formula(f1_str), data = df, REML = FALSE, control = lmer_ctrl, method = est_method))
      lrt01 <- suppressWarnings(anova(m0_ml, m1_ml))

      result$ri_model <- list(
        formula        = f1_str,
        fixed_effects  = extract_fixed(m1),
        random_effects = extract_random(m1, df),
        icc            = calc_icc(m1),
        r2             = calc_r2(m1),
        fit            = extract_fit(m1_ml),
        lrt_vs_null    = list(
          chi2 = round(lrt01[["Chisq"]][2], 4),
          df   = lrt01[["Df"]][2],
          p    = round(lrt01[["Pr(>Chisq)"]][2], 4)
        ),
        lrt_vs_lm      = calc_lrt_vs_lm(m1_ml, f1_fixed, df)
      )

      rcode <- c(rcode,
        "# ── 5. Model 1: 무선절편 모형 (Random Intercept) ──",
        "# 고정효과: L1 독립변수 + L2 공변량, 무선효과: 절편만",
        paste0('m1 <- lmer(', f1_str, ','),
        "             data = df, REML = TRUE,",
        "             control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5)))",
        "summary(m1)  # 고정효과 계수, 무선효과 분산 확인",
        "",
        "# ML 추정 (모형 비교용)",
        paste0('m1_ml <- lmer(', f1_str, ','),
        "               data = df, REML = FALSE,",
        "               control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5)))",
        "",
        "# Model 0 vs Model 1: 우도비 검정 (ML 기반)",
        "anova(m0_ml, m1_ml)  # 독립변수 투입 효과가 유의한지 검정",
        ""
      )

      # ══════════════════════════════════════════════════════════════
      # Model 2: 무선기울기 모형 (Random Slope Model)
      # ══════════════════════════════════════════════════════════════
      # 고정효과: M1과 동일
      # 무선효과: 절편 + 선택된 L1 변수의 기울기도 집단 간 변동 허용
      # 의미: 독립변수의 효과(기울기)가 집단마다 다를 수 있음
      #
      # 공분산 구조:
      #   un (Unstructured): (1 + X | group) → 절편-기울기 공분산도 추정
      #   ind (Independent):  (1 + X || group) → 공분산 = 0으로 제약
      valid_rs <- intersect(rand_slopes, l1_preds)
      if (length(valid_rs) > 0) {
        slope_part <- paste(c("1", valid_rs), collapse = " + ")
        re_bar <- if (cov_struct == "ind") " || " else " | "
        f2_str <- paste0(outcome, " ~ ",
                         paste(all_fixed, collapse = " + "),
                         " + (", slope_part, re_bar, group_var, ")")

        # 에러 시 문자열 반환 → is.character()로 판별 (S4 모형에 $ 사용 금지)
        m2 <- tryCatch(
          suppressWarnings(fit_model(as.formula(f2_str), data = df, REML = TRUE, control = lmer_ctrl, method = est_method)),
          error = function(e) e$message
        )

        if (is.character(m2)) {
          result$rs_model_error <- paste0("수렴 실패: ", m2)
        } else {
          m2_ml <- tryCatch(
            suppressWarnings(fit_model(as.formula(f2_str), data = df, REML = FALSE, control = lmer_ctrl, method = est_method)),
            error = function(e) e$message
          )
          lrt12 <- if (!is.character(m2_ml)) tryCatch(suppressWarnings(anova(m1_ml, m2_ml)), error = function(e) NULL) else NULL

          result$rs_model <- list(
            formula        = f2_str,
            fixed_effects  = extract_fixed(m2),
            random_effects = extract_random(m2, df),
            icc            = calc_icc(m2),
            r2             = calc_r2(m2),
            fit            = if (!is.character(m2_ml)) extract_fit(m2_ml) else extract_fit(m2),
            lrt_vs_ri      = if (!is.null(lrt12)) list(
              chi2 = round(lrt12[["Chisq"]][2], 4),
              df   = lrt12[["Df"]][2],
              p    = round(lrt12[["Pr(>Chisq)"]][2], 4)
            ) else NULL,
            lrt_vs_lm      = if (!is.character(m2_ml)) calc_lrt_vs_lm(m2_ml, f1_fixed, df) else NULL
          )

          rcode <- c(rcode,
            "# ── 6. Model 2: 무선기울기 모형 (Random Slope) ──",
            "# 무선절편에 무선기울기를 추가하여 집단별 기울기 차이를 허용",
            "# 즉, 독립변수의 효과가 집단마다 다를 수 있는지 검정",
            paste0('m2 <- lmer(', f2_str, ','),
            "             data = df, REML = TRUE,",
            "             control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5)))",
            "summary(m2)  # 고정효과 + 무선기울기 분산 확인",
            "",
            "# ML 추정 (모형 비교용)",
            paste0('m2_ml <- lmer(', f2_str, ','),
            "               data = df, REML = FALSE,",
            "               control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5)))",
            "",
            "# Model 1 vs Model 2: 우도비 검정 (ML 기반)",
            "anova(m1_ml, m2_ml)  # 무선기울기 추가가 유의한지 검정",
            ""
          )
        }
      }

      # ══════════════════════════════════════════════════════════════
      # Model 3: 교차수준 상호작용 모형 (Cross-Level Interaction)
      # ══════════════════════════════════════════════════════════════
      # 고정효과: M1 + L1×L2 상호작용항
      # 무선효과: M2와 동일 (상호작용에 포함된 L1 변수 자동 추가)
      # 의미: L2 변수(집단 특성)가 L1 변수의 기울기에 미치는 영향 검정
      #       예: 학교 크기(L2)가 SES→성취도 기울기(L1)에 영향을 주는가?
      if (length(cross_interactions) > 0) {
        int_terms   <- character(0)
        int_l1_vars <- character(0)
        for (i in seq_along(cross_interactions)) {
          ci <- cross_interactions[[i]]
          if (!is.list(ci)) next
          l1v <- if (!is.null(ci[["l1"]])) as.character(ci[["l1"]]) else NA
          l2v <- if (!is.null(ci[["l2"]])) as.character(ci[["l2"]]) else NA
          if (is.na(l1v) || is.na(l2v)) next
          if (l1v %in% l1_preds && l2v %in% l2_covs) {
            int_terms   <- c(int_terms, paste0(l1v, ":", l2v))
            int_l1_vars <- c(int_l1_vars, l1v)
          }
        }

        if (length(int_terms) > 0) {
          int_l1_vars <- unique(int_l1_vars)
          all_rs_m3   <- unique(c(valid_rs, int_l1_vars))
          fixed_part  <- paste(c(all_fixed, int_terms), collapse = " + ")
          slope_part3 <- paste(c("1", all_rs_m3), collapse = " + ")
          re_bar3 <- if (cov_struct == "ind") " || " else " | "
          f3_str <- paste0(outcome, " ~ ", fixed_part,
                           " + (", slope_part3, re_bar3, group_var, ")")

          m3 <- tryCatch(
            suppressWarnings(fit_model(as.formula(f3_str), data = df, REML = TRUE, control = lmer_ctrl, method = est_method)),
            error = function(e) e$message
          )

          if (is.character(m3)) {
            result$cl_model_error <- paste0("수렴 실패: ", m3)
          } else {
            m3_ml <- tryCatch(
              suppressWarnings(fit_model(as.formula(f3_str), data = df, REML = FALSE, control = lmer_ctrl, method = est_method)),
              error = function(e) e$message
            )
            prev_m_ml <- m1_ml
            if (length(valid_rs) > 0 && exists("m2_ml") && !is.character(m2_ml)) prev_m_ml <- m2_ml
            lrt_m3 <- if (!is.character(m3_ml)) tryCatch(suppressWarnings(anova(prev_m_ml, m3_ml)), error = function(e) NULL) else NULL
            f3_fixed <- as.formula(paste0(outcome, " ~ ", fixed_part))

            # lrt_vs_prev 안전하게 구축 ($ 대신 [[ 사용)
            lrt_vs_prev_val <- NULL
            if (!is.null(lrt_m3)) {
              tryCatch({
                lrt_vs_prev_val <- list(
                  chi2 = round(lrt_m3[["Chisq"]][2], 4),
                  df   = lrt_m3[["Df"]][2],
                  p    = round(lrt_m3[["Pr(>Chisq)"]][2], 4)
                )
              }, error = function(e) {
                message("[M3 lrt_vs_prev] 오류: ", e$message)
              })
            }

            result$cl_model <- tryCatch(list(
              formula           = f3_str,
              fixed_effects     = extract_fixed(m3),
              random_effects    = extract_random(m3, df),
              icc               = calc_icc(m3),
              r2                = calc_r2(m3),
              fit               = if (!is.character(m3_ml)) extract_fit(m3_ml) else extract_fit(m3),
              interaction_terms = as.list(int_terms),
              lrt_vs_prev       = lrt_vs_prev_val,
              lrt_vs_lm         = if (!is.character(m3_ml)) calc_lrt_vs_lm(m3_ml, f3_fixed, df) else NULL
            ), error = function(e) {
              message("[M3 result build] 오류: ", e$message)
              list(formula = f3_str, error = e$message)
            })

            prev_label <- if (length(valid_rs) > 0 && !is.null(result$rs_model)) "m2_ml" else "m1_ml"
            rcode <- c(rcode,
              "# ── 7. Model 3: 교차수준 상호작용 모형 (Cross-Level Interaction) ──",
              "# Level-1 변수와 Level-2 변수의 교차수준 상호작용 효과 검정",
              "# 집단 수준 변수가 개인 수준 변수의 기울기에 미치는 영향 분석",
              paste0('m3 <- lmer(', f3_str, ','),
              "             data = df, REML = TRUE,",
              "             control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5)))",
              "summary(m3)  # 교차수준 상호작용 계수 확인",
              "",
              "# ML 추정 (모형 비교용)",
              paste0('m3_ml <- lmer(', f3_str, ','),
              "               data = df, REML = FALSE,",
              "               control = lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 2e5)))",
              "",
              "# 이전 모형 vs Model 3: 우도비 검정 (ML 기반)",
              paste0("anova(", prev_label, ", m3_ml)  # 교차수준 상호작용 효과 검정"),
              ""
            )
          }
        }
      }
    }

    # ── 진단 데이터 추출 (최종 모형) ─────────────────────────────
    # result 에 이미 저장된 모형 키를 기준으로 최종 모형 결정 (exists 스코프 문제 회피)
    final_model <- m0
    final_name  <- "Model 0: 기저모형"
    if (!is.null(result$ri_model)) {
      final_model <- m1; final_name <- "Model 1: 무선절편"
    }
    if (!is.null(result$rs_model)) {
      final_model <- m2; final_name <- "Model 2: 무선기울기"
    }
    if (!is.null(result$cl_model)) {
      final_model <- m3; final_name <- "Model 3: 교차수준 상호작용"
    }
    result$diagnostics <- tryCatch(
      extract_diagnostics(final_model, df, outcome, group_var, l1_preds, l2_covs, final_name),
      error = function(e) list(model_name = final_name, error = e$message)
    )

    # ── 모형 비교표 ───────────────────────────────────────────────
    # 안전한 접근 헬퍼: 중첩 리스트에서 값 추출 ($ 대신 [[ 사용)
    safe_get <- function(lst, ...) {
      keys <- list(...)
      val <- lst
      for (k in keys) {
        if (is.null(val) || !is.list(val)) return(NA)
        val <- val[[k]]
      }
      if (is.null(val)) NA else val
    }

    comp <- list(list(
      name     = "Model 0: 기저모형",
      npar     = safe_get(result, "null_model", "fit", "npar"),
      AIC      = safe_get(result, "null_model", "fit", "AIC"),
      BIC      = safe_get(result, "null_model", "fit", "BIC"),
      logLik   = safe_get(result, "null_model", "fit", "logLik"),
      deviance = safe_get(result, "null_model", "fit", "deviance"),
      delta_deviance = NA, chi2 = NA, df_chi = NA, p_chi = NA
    ))

    if (!is.null(result[["ri_model"]])) {
      lrt <- result[["ri_model"]][["lrt_vs_null"]]
      ri_dev <- safe_get(result, "ri_model", "fit", "deviance")
      m0_dev <- safe_get(result, "null_model", "fit", "deviance")
      comp[[2]] <- list(
        name     = "Model 1: 무선절편",
        npar     = safe_get(result, "ri_model", "fit", "npar"),
        AIC      = safe_get(result, "ri_model", "fit", "AIC"),
        BIC      = safe_get(result, "ri_model", "fit", "BIC"),
        logLik   = safe_get(result, "ri_model", "fit", "logLik"),
        deviance = ri_dev,
        delta_deviance = if (!is.na(m0_dev) && !is.na(ri_dev)) round(m0_dev - ri_dev, 2) else NA,
        chi2     = safe_get(lrt, "chi2"),
        df_chi   = safe_get(lrt, "df"),
        p_chi    = safe_get(lrt, "p")
      )
    }

    if (!is.null(result[["rs_model"]])) {
      lrt <- result[["rs_model"]][["lrt_vs_ri"]]
      rs_dev <- safe_get(result, "rs_model", "fit", "deviance")
      ri_dev2 <- safe_get(result, "ri_model", "fit", "deviance")
      comp[[3]] <- list(
        name     = "Model 2: 무선기울기",
        npar     = safe_get(result, "rs_model", "fit", "npar"),
        AIC      = safe_get(result, "rs_model", "fit", "AIC"),
        BIC      = safe_get(result, "rs_model", "fit", "BIC"),
        logLik   = safe_get(result, "rs_model", "fit", "logLik"),
        deviance = rs_dev,
        delta_deviance = if (!is.null(lrt) && !is.na(ri_dev2) && !is.na(rs_dev)) round(ri_dev2 - rs_dev, 2) else NA,
        chi2     = if (!is.null(lrt)) safe_get(lrt, "chi2") else NA,
        df_chi   = if (!is.null(lrt)) safe_get(lrt, "df")   else NA,
        p_chi    = if (!is.null(lrt)) safe_get(lrt, "p")    else NA
      )
    }

    if (!is.null(result[["cl_model"]]) && is.null(result[["cl_model"]][["error"]])) {
      lrt <- result[["cl_model"]][["lrt_vs_prev"]]
      cl_dev <- safe_get(result, "cl_model", "fit", "deviance")
      prev_dev <- if (!is.null(result[["rs_model"]])) {
        safe_get(result, "rs_model", "fit", "deviance")
      } else {
        safe_get(result, "ri_model", "fit", "deviance")
      }
      comp[[length(comp) + 1]] <- list(
        name     = "Model 3: 교차수준 상호작용",
        npar     = safe_get(result, "cl_model", "fit", "npar"),
        AIC      = safe_get(result, "cl_model", "fit", "AIC"),
        BIC      = safe_get(result, "cl_model", "fit", "BIC"),
        logLik   = safe_get(result, "cl_model", "fit", "logLik"),
        deviance = cl_dev,
        delta_deviance = if (!is.null(lrt) && !is.na(prev_dev) && !is.na(cl_dev)) round(prev_dev - cl_dev, 2) else NA,
        chi2     = if (!is.null(lrt)) safe_get(lrt, "chi2") else NA,
        df_chi   = if (!is.null(lrt)) safe_get(lrt, "df")   else NA,
        p_chi    = if (!is.null(lrt)) safe_get(lrt, "p")    else NA
      )
    }

    result$model_comparison <- comp

    # ── R 코드 로그 마무리 ─────────────────────────────────────────
    rcode <- c(rcode,
      "# ── 모형 비교 (Model Comparison) ──",
      "# 중첩 모형 간 우도비 검정(LRT)으로 최적 모형 선택",
      "# LRT는 ML(REML=FALSE) 기반으로 수행해야 정확",
      "# AIC/BIC가 작을수록, LRT p값이 유의할수록 해당 모형이 우수"
    )
    if (!is.null(result$ri_model)) {
      rcode <- c(rcode, "anova(m0_ml, m1_ml)  # 기저모형 vs 무선절편 (ML 기반)")
    }
    if (!is.null(result$rs_model)) {
      rcode <- c(rcode, "anova(m1_ml, m2_ml)  # 무선절편 vs 무선기울기 (ML 기반)")
    }
    if (!is.null(result$cl_model)) {
      if (!is.null(result$rs_model)) {
        rcode <- c(rcode, "anova(m2_ml, m3_ml)  # 무선기울기 vs 교차수준 (ML 기반)")
      } else {
        rcode <- c(rcode, "anova(m1_ml, m3_ml)  # 무선절편 vs 교차수준 (ML 기반)")
      }
    }
    rcode <- c(rcode, "")
    result$r_code <- paste(rcode, collapse = "\n")

    result

  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
}
