#!/usr/bin/env bash
# ========================================
#   다층모형분석 웹앱 (MLM Web App)
#   외부 접속 가능 (포트 8001, sudo 불필요)
# ========================================

set -e

cd "$(dirname "$0")"

# conda 환경 활성화
CONDA_BASE="${CONDA_BASE:-$HOME/anaconda3}"
RSCRIPT="${CONDA_BASE}/envs/hlm_app/bin/Rscript"

if [ ! -f "$RSCRIPT" ]; then
  echo "[오류] hlm_app conda 환경을 찾을 수 없습니다."
  echo "  다음 명령으로 생성하세요:"
  echo "  conda create -n hlm_app -c conda-forge r-base r-plumber r-lme4 r-lmertest -y"
  exit 1
fi

echo "========================================"
echo "  다층모형분석 웹앱 (MLM Web App)"
echo "========================================"
echo ""

# [1단계] 패키지 확인
echo "[1단계] 패키지 확인 중..."
$RSCRIPT -e "
pkgs <- c('plumber','lme4','lmerTest')
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]
if(length(missing)>0){
  cat('누락 패키지:', paste(missing, collapse=', '), '\n')
  quit(status=1)
}
cat('모든 패키지 확인 완료\n')
"

# [2단계] 서버 시작
echo ""
echo "[2단계] R API 서버 시작 (포트: 8001, 외부 접속 허용)"
echo ""
echo "  접속 주소:"
echo "    로컬:  http://localhost:8001"

# 서버 IP 자동 감지
SERVER_IP=$(hostname -I 2>/dev/null | awk '{print $1}')
if [ -n "$SERVER_IP" ]; then
  echo "    외부:  http://${SERVER_IP}:8001"
fi

echo ""
echo "  종료: Ctrl+C"
echo ""

$RSCRIPT start_server.R
