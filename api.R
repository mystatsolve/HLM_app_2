# api.R - 다층모형분석 Plumber API
# 패키지: plumber, lme4, lmerTest

library(plumber)
library(lme4)

# lmerTest: Satterthwaite 자유도로 p값 계산
has_lmerTest <- requireNamespace("lmerTest", quietly = TRUE)
if (has_lmerTest) {
  library(lmerTest)
  message("[OK] lmerTest 로드 - p값 계산 가능")
} else {
  message("[경고] lmerTest 없음 - t값만 제공. install.R을 먼저 실행하세요.")
}

# 정적 파일 제공 (www/ 폴더)
#* @assets ./www /
list()

# ─────────────────────────────────────────
#  유틸리티 함수
# ─────────────────────────────────────────

# ICC 계산 (집단내 상관계수)
calc_icc <- function(model) {
  tryCatch({
    vc <- as.data.frame(VarCorr(model))
    int_var <- vc$vcov[!is.na(vc$var1) & vc$var1 == "(Intercept)" & is.na(vc$var2)]
    res_var <- vc$vcov[vc$grp == "Residual"]
    if (length(int_var) == 0 || length(res_var) == 0) return(NA)
    round(int_var[1] / (int_var[1] + res_var[1]), 4)
  }, error = function(e) NA)
}

# R² 계산 (Nakagawa & Schielzeth, 2013)
calc_r2 <- function(model) {
  tryCatch({
    vc    <- as.data.frame(VarCorr(model))
    sig2  <- sigma(model)^2
    X     <- model.matrix(model)
    beta  <- fixef(model)
    var_f <- var(as.vector(X %*% beta))
    var_r <- sum(vc$vcov[vc$grp != "Residual"])
    tot   <- var_f + var_r + sig2
    list(
      marginal    = round(var_f / tot, 4),
      conditional = round((var_f + var_r) / tot, 4)
    )
  }, error = function(e) list(marginal = NA, conditional = NA))
}

# 모형 적합도 지수
extract_fit <- function(model) {
  list(
    AIC      = round(AIC(model), 2),
    BIC      = round(BIC(model), 2),
    logLik   = round(as.numeric(logLik(model)), 2),
    deviance = round(-2 * as.numeric(logLik(model)), 2),
    npar     = attr(logLik(model), "df")
  )
}

# 고정효과 추출 (추정값, SE, t, p, 95% CI)
extract_fixed <- function(model) {
  # coef(summary()) → S3/S4 모두 안전 (summary()$coefficients 는 S4에서 $ 오류 발생)
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

# 무선효과 추출
extract_random <- function(model) {
  vc <- as.data.frame(VarCorr(model))
  lapply(seq_len(nrow(vc)), function(i) {
    list(
      group = vc$grp[i],
      var1  = ifelse(is.na(vc$var1[i]), "", as.character(vc$var1[i])),
      var2  = ifelse(is.na(vc$var2[i]), "", as.character(vc$var2[i])),
      vcov  = round(vc$vcov[i], 4),
      sdcor = round(vc$sdcor[i], 4)
    )
  })
}

# ─────────────────────────────────────────
#  Plumber 엔드포인트
# ─────────────────────────────────────────

#* CORS 필터
#* @filter cors
function(req, res) {
  res$setHeader("Access-Control-Allow-Origin",  "*")
  res$setHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
  res$setHeader("Access-Control-Allow-Headers", "Content-Type, Accept")
  if (req$REQUEST_METHOD == "OPTIONS") { res$status <- 200; return(list()) }
  plumber::forward()
}

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
    cross_interactions <- if (!is.null(body$cross_interactions)) body$cross_interactions else list()
    # jsonlite가 [{l1:"a",l2:"b"},...] 를 data.frame으로 변환하므로 list-of-list로 정규화
    if (is.data.frame(cross_interactions)) {
      cross_interactions <- lapply(seq_len(nrow(cross_interactions)), function(i) {
        list(l1 = as.character(cross_interactions$l1[i]),
             l2 = as.character(cross_interactions$l2[i]))
      })
    } else if (length(cross_interactions) > 0 && !is.null(names(cross_interactions))) {
      # 단일 상호작용이 flat named list로 온 경우
      cross_interactions <- list(cross_interactions)
    }

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

    result <- list(
      success   = TRUE,
      data_info = list(
        n_obs    = n_obs,
        n_groups = n_groups,
        outcome  = outcome,
        group_var = group_var
      )
    )

    # ── Model 0: 기저모형 (Null) ──────────────────────────────────
    f0 <- as.formula(paste0(outcome, " ~ 1 + (1|", group_var, ")"))
    m0 <- suppressWarnings(lmer(f0, data = df, REML = FALSE))

    result$null_model <- list(
      formula        = deparse(f0),
      icc            = calc_icc(m0),
      random_effects = extract_random(m0),
      fit            = extract_fit(m0)
    )

    # ── Model 1: 무선절편 모형 (Random Intercept) ─────────────────
    all_fixed <- c(l1_preds, l2_covs)
    if (length(all_fixed) > 0) {
      f1_str <- paste0(outcome, " ~ ",
                       paste(all_fixed, collapse = " + "),
                       " + (1|", group_var, ")")
      m1 <- suppressWarnings(lmer(as.formula(f1_str), data = df, REML = FALSE))
      lrt01 <- suppressWarnings(anova(m0, m1))

      result$ri_model <- list(
        formula        = f1_str,
        fixed_effects  = extract_fixed(m1),
        random_effects = extract_random(m1),
        icc            = calc_icc(m1),
        r2             = calc_r2(m1),
        fit            = extract_fit(m1),
        lrt_vs_null    = list(
          chi2 = round(lrt01[["Chisq"]][2], 4),
          df   = lrt01[["Df"]][2],
          p    = round(lrt01[["Pr(>Chisq)"]][2], 4)
        )
      )

      # ── Model 2: 무선기울기 모형 (Random Slope) ─────────────────
      valid_rs <- intersect(rand_slopes, l1_preds)
      if (length(valid_rs) > 0) {
        slope_part <- paste(c("1", valid_rs), collapse = " + ")
        f2_str <- paste0(outcome, " ~ ",
                         paste(all_fixed, collapse = " + "),
                         " + (", slope_part, "|", group_var, ")")

        # 에러 시 문자열 반환 → is.character()로 판별 (S4 모형에 $ 사용 금지)
        m2 <- tryCatch(
          suppressWarnings(
            lmer(as.formula(f2_str), data = df, REML = FALSE,
                 control = lmerControl(optimizer = "bobyqa",
                                       optCtrl   = list(maxfun = 2e5)))
          ),
          error = function(e) e$message
        )

        if (is.character(m2)) {
          result$rs_model_error <- paste0("수렴 실패: ", m2)
        } else {
          lrt12 <- tryCatch(suppressWarnings(anova(m1, m2)), error = function(e) NULL)

          result$rs_model <- list(
            formula        = f2_str,
            fixed_effects  = extract_fixed(m2),
            random_effects = extract_random(m2),
            icc            = calc_icc(m2),
            r2             = calc_r2(m2),
            fit            = extract_fit(m2),
            lrt_vs_ri      = if (!is.null(lrt12)) list(
              chi2 = round(lrt12[["Chisq"]][2], 4),
              df   = lrt12[["Df"]][2],
              p    = round(lrt12[["Pr(>Chisq)"]][2], 4)
            ) else NULL
          )
        }
      }

      # ── Model 3: 교차수준 상호작용 모형 (Cross-Level Interaction) ──
      if (length(cross_interactions) > 0) {
        int_terms   <- character(0)
        int_l1_vars <- character(0)
        for (ci in cross_interactions) {
          l1v <- as.character(ci$l1)
          l2v <- as.character(ci$l2)
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
          f3_str <- paste0(outcome, " ~ ", fixed_part,
                           " + (", slope_part3, "|", group_var, ")")

          m3 <- tryCatch(
            suppressWarnings(
              lmer(as.formula(f3_str), data = df, REML = FALSE,
                   control = lmerControl(optimizer = "bobyqa",
                                         optCtrl   = list(maxfun = 2e5)))
            ),
            error = function(e) e$message
          )

          if (is.character(m3)) {
            result$cl_model_error <- paste0("수렴 실패: ", m3)
          } else {
            prev_m <- m1
            if (length(valid_rs) > 0 && exists("m2") && !is.character(m2)) prev_m <- m2
            lrt_m3 <- tryCatch(suppressWarnings(anova(prev_m, m3)), error = function(e) NULL)

            result$cl_model <- list(
              formula           = f3_str,
              fixed_effects     = extract_fixed(m3),
              random_effects    = extract_random(m3),
              icc               = calc_icc(m3),
              r2                = calc_r2(m3),
              fit               = extract_fit(m3),
              interaction_terms = as.list(int_terms),
              lrt_vs_prev       = if (!is.null(lrt_m3)) list(
                chi2 = round(lrt_m3[["Chisq"]][2], 4),
                df   = lrt_m3[["Df"]][2],
                p    = round(lrt_m3[["Pr(>Chisq)"]][2], 4)
              ) else NULL
            )
          }
        }
      }
    }

    # ── 모형 비교표 ───────────────────────────────────────────────
    comp <- list(list(
      name     = "Model 0: 기저모형",
      npar     = result$null_model$fit$npar,
      AIC      = result$null_model$fit$AIC,
      BIC      = result$null_model$fit$BIC,
      logLik   = result$null_model$fit$logLik,
      deviance = result$null_model$fit$deviance,
      delta_deviance = NA, chi2 = NA, df_chi = NA, p_chi = NA
    ))

    if (!is.null(result$ri_model)) {
      lrt <- result$ri_model$lrt_vs_null
      comp[[2]] <- list(
        name     = "Model 1: 무선절편",
        npar     = result$ri_model$fit$npar,
        AIC      = result$ri_model$fit$AIC,
        BIC      = result$ri_model$fit$BIC,
        logLik   = result$ri_model$fit$logLik,
        deviance = result$ri_model$fit$deviance,
        delta_deviance = round(result$null_model$fit$deviance - result$ri_model$fit$deviance, 2),
        chi2 = lrt$chi2, df_chi = lrt$df, p_chi = lrt$p
      )
    }

    if (!is.null(result$rs_model)) {
      lrt <- result$rs_model$lrt_vs_ri
      comp[[3]] <- list(
        name     = "Model 2: 무선기울기",
        npar     = result$rs_model$fit$npar,
        AIC      = result$rs_model$fit$AIC,
        BIC      = result$rs_model$fit$BIC,
        logLik   = result$rs_model$fit$logLik,
        deviance = result$rs_model$fit$deviance,
        delta_deviance = if (!is.null(lrt)) round(result$ri_model$fit$deviance - result$rs_model$fit$deviance, 2) else NA,
        chi2     = if (!is.null(lrt)) lrt$chi2 else NA,
        df_chi   = if (!is.null(lrt)) lrt$df   else NA,
        p_chi    = if (!is.null(lrt)) lrt$p    else NA
      )
    }

    if (!is.null(result$cl_model)) {
      lrt <- result$cl_model$lrt_vs_prev
      prev_dev <- if (!is.null(result$rs_model)) result$rs_model$fit$deviance
                  else result$ri_model$fit$deviance
      comp[[length(comp) + 1]] <- list(
        name     = "Model 3: 교차수준 상호작용",
        npar     = result$cl_model$fit$npar,
        AIC      = result$cl_model$fit$AIC,
        BIC      = result$cl_model$fit$BIC,
        logLik   = result$cl_model$fit$logLik,
        deviance = result$cl_model$fit$deviance,
        delta_deviance = if (!is.null(lrt)) round(prev_dev - result$cl_model$fit$deviance, 2) else NA,
        chi2     = if (!is.null(lrt)) lrt$chi2 else NA,
        df_chi   = if (!is.null(lrt)) lrt$df   else NA,
        p_chi    = if (!is.null(lrt)) lrt$p    else NA
      )
    }

    result$model_comparison <- comp
    result

  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
}
