# ═══════════════════════════════════════════════════════════════════════════════
# start_server.R - R Plumber API 서버 시작 스크립트
# ═══════════════════════════════════════════════════════════════════════════════
#
# 사용법:
#   Rscript start_server.R          (터미널에서 직접 실행)
#   source("start_server.R")        (R 콘솔에서 실행)
#
# 실행 후 브라우저에서 http://localhost:8001 접속
#
# ── 작업 디렉토리 자동 설정 ──────────────────────────────────────────────────
# api.R과 www/ 폴더가 이 스크립트와 같은 디렉토리에 있어야 합니다.
# 실행 방식에 따라 스크립트 경로를 자동 감지합니다:
#   - source("start_server.R"): sys.frame(1)$ofile 에서 경로 추출
#   - Rscript start_server.R:   commandArgs()의 --file= 인자에서 경로 추출
#   - 둘 다 실패하면:          현재 작업 디렉토리(getwd()) 사용
# ═══════════════════════════════════════════════════════════════════════════════

script_path <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(sub("^--file=", "", file_arg))
    else getwd()
  }
)
setwd(script_path)

# Plumber API 서버 시작
# host = "0.0.0.0": 모든 네트워크 인터페이스에서 접근 가능 (외부 접속 허용)
# port = 8001: 포트 번호 (변경 시 index.html의 API_BASE도 수정 필요)
# docs = FALSE: Swagger UI 비활성화 (프론트엔드가 별도로 있으므로 불필요)
library(plumber)
pr <- plumb("api.R")
pr_run(pr, host = "0.0.0.0", port = 8001, docs = FALSE)
