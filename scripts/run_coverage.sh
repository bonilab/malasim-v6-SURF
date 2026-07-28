#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

BUILD_DIR="${REPO_ROOT}/build/coverage"
OUTPUT_ROOT="${REPO_ROOT}/coverage-results"
RUN_NAME=""
GTEST_FILTER="*"
BUILD_TYPE="Debug"
GENERATE_HTML=1
UPDATE_BASELINE=0
BASELINE_FILE="${REPO_ROOT}/tests/coverage-baseline.json"

usage() {
  cat <<'EOF'
Run the MalaSim unit tests with LLVM source coverage.

Usage:
  scripts/run_coverage.sh [options]

Options:
  --build-dir PATH    Instrumented CMake build directory.
                      Default: build/coverage
  --output-dir PATH   Root directory for saved coverage runs.
                      Default: coverage-results
  --run-name NAME     Name for this run. Default: current timestamp
  --filter PATTERN    GoogleTest filter. Default: *
  --build-type TYPE   CMake build type. Default: Debug
  --no-html           Skip the detailed HTML report.
  --update-baseline   Refresh tests/coverage-baseline.json from a full run.
  -h, --help          Show this help.

Environment:
  VCPKG_ROOT          Used by the repository CMake configuration when set.
  CC, CXX             Set these to clang and clang++ when Clang is not the
                      platform default.

Artifacts:
  <output-dir>/runs/<run-name>/coverage-summary.txt
  <output-dir>/runs/<run-name>/coverage-summary.json
  <output-dir>/runs/<run-name>/test-output.log
  <output-dir>/runs/<run-name>/coverage-warnings.log
  <output-dir>/runs/<run-name>/html/index.html
  <output-dir>/history.csv
  <output-dir>/latest-run.txt
  <output-dir>/latest-summary.txt
  <output-dir>/latest-summary.json
  tests/coverage-baseline.json (only with --update-baseline)
EOF
}

require_value() {
  if [[ $# -lt 2 || -z "$2" ]]; then
    echo "error: $1 requires a value" >&2
    exit 2
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --build-dir)
      require_value "$@"
      BUILD_DIR="$2"
      shift 2
      ;;
    --output-dir)
      require_value "$@"
      OUTPUT_ROOT="$2"
      shift 2
      ;;
    --run-name)
      require_value "$@"
      RUN_NAME="$2"
      shift 2
      ;;
    --filter)
      require_value "$@"
      GTEST_FILTER="$2"
      shift 2
      ;;
    --build-type)
      require_value "$@"
      BUILD_TYPE="$2"
      shift 2
      ;;
    --no-html)
      GENERATE_HTML=0
      shift
      ;;
    --update-baseline)
      UPDATE_BASELINE=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "error: unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ "${UPDATE_BASELINE}" -eq 1 && "${GTEST_FILTER}" != "*" ]]; then
  echo "error: --update-baseline requires a full run without --filter" >&2
  exit 2
fi

if [[ -z "${RUN_NAME}" ]]; then
  RUN_NAME="$(date '+%Y%m%d-%H%M%S')"
fi

if [[ "${RUN_NAME}" == */* || "${RUN_NAME}" == "." || "${RUN_NAME}" == ".." ]]; then
  echo "error: --run-name must be a single directory name" >&2
  exit 2
fi

case "${BUILD_DIR}" in
  /*) ;;
  *) BUILD_DIR="${REPO_ROOT}/${BUILD_DIR}" ;;
esac

case "${OUTPUT_ROOT}" in
  /*) ;;
  *) OUTPUT_ROOT="${REPO_ROOT}/${OUTPUT_ROOT}" ;;
esac

RUNS_DIR="${OUTPUT_ROOT}/runs"
RUN_DIR="${RUNS_DIR}/${RUN_NAME}"
if [[ -e "${RUN_DIR}" ]]; then
  RUN_DIR="${RUN_DIR}-$$"
  RUN_NAME="$(basename "${RUN_DIR}")"
fi

RAW_DIR="${RUN_DIR}/raw"
HTML_DIR="${RUN_DIR}/html"
TEST_LOG="${RUN_DIR}/test-output.log"
SUMMARY_FILE="${RUN_DIR}/coverage-summary.txt"
JSON_FILE="${RUN_DIR}/coverage-summary.json"
PROFDATA_FILE="${RUN_DIR}/coverage.profdata"
METADATA_FILE="${RUN_DIR}/metadata.txt"
COVERAGE_TOOL_LOG="${RUN_DIR}/coverage-warnings.log"
HISTORY_FILE="${OUTPUT_ROOT}/history.csv"
LATEST_FILE="${OUTPUT_ROOT}/latest-run.txt"
LATEST_SUMMARY_FILE="${OUTPUT_ROOT}/latest-summary.txt"
LATEST_JSON_FILE="${OUTPUT_ROOT}/latest-summary.json"

find_llvm_tool() {
  local tool_name="$1"
  if command -v "${tool_name}" >/dev/null 2>&1; then
    command -v "${tool_name}"
    return
  fi
  if command -v xcrun >/dev/null 2>&1; then
    xcrun --find "${tool_name}" 2>/dev/null && return
  fi
  echo "error: ${tool_name} was not found; install LLVM or Xcode command-line tools" >&2
  exit 1
}

command -v cmake >/dev/null 2>&1 || {
  echo "error: cmake was not found" >&2
  exit 1
}

LLVM_COV="$(find_llvm_tool llvm-cov)"
LLVM_PROFDATA="$(find_llvm_tool llvm-profdata)"

mkdir -p "${RAW_DIR}"

CMAKE_GENERATOR=""
if [[ ! -f "${BUILD_DIR}/CMakeCache.txt" ]] && command -v ninja >/dev/null 2>&1; then
  CMAKE_GENERATOR="Ninja"
fi

echo "Configuring coverage build: ${BUILD_DIR}"
CMAKE_CONFIGURE_ARGS=(
  -S "${REPO_ROOT}"
  -B "${BUILD_DIR}"
  -DCMAKE_BUILD_TYPE="${BUILD_TYPE}"
  -DENABLE_COVERAGE=ON
)
if [[ -n "${CMAKE_GENERATOR}" ]]; then
  cmake -G "${CMAKE_GENERATOR}" "${CMAKE_CONFIGURE_ARGS[@]}"
else
  cmake "${CMAKE_CONFIGURE_ARGS[@]}"
fi

echo "Building instrumented test executable"
cmake --build "${BUILD_DIR}" --target malasim_test --config "${BUILD_TYPE}" -j 6

TEST_BINARY="${BUILD_DIR}/bin/malasim_test"
if [[ ! -x "${TEST_BINARY}" ]]; then
  TEST_BINARY="${BUILD_DIR}/bin/${BUILD_TYPE}/malasim_test"
fi
if [[ ! -x "${TEST_BINARY}" ]]; then
  echo "error: instrumented test executable was not found under ${BUILD_DIR}/bin" >&2
  exit 1
fi

PROFILE_PATTERN="${RAW_DIR}/malasim-%p.profraw"
echo "Running tests with filter: ${GTEST_FILTER}"
if (
  cd "$(dirname "${TEST_BINARY}")"
  LLVM_PROFILE_FILE="${PROFILE_PATTERN}" \
    "./$(basename "${TEST_BINARY}")" --gtest_filter="${GTEST_FILTER}"
) > "${TEST_LOG}" 2>&1; then
  grep -E '^\[(==========|  PASSED  |  SKIPPED |  FAILED  )\]' "${TEST_LOG}" | tail -n 6 || true
else
  echo "error: coverage test run failed; complete output follows" >&2
  cat "${TEST_LOG}" >&2
  exit 1
fi

RAW_PROFILES=("${RAW_DIR}"/*.profraw)
if [[ ! -e "${RAW_PROFILES[0]}" ]]; then
  echo "error: the test run did not produce LLVM raw profile data" >&2
  exit 1
fi

echo "Merging raw profiles"
"${LLVM_PROFDATA}" merge -sparse "${RAW_PROFILES[@]}" -o "${PROFDATA_FILE}"

# Report only MalaSim production code. Tests and dependency sources are not part
# of the product coverage target.
IGNORE_REGEX='/(tests|vcpkg|vcpkg_installed|build)/|/usr/|/Applications/Xcode'
COVERAGE_ARGS=(
  "${TEST_BINARY}"
  "-instr-profile=${PROFDATA_FILE}"
  "-ignore-filename-regex=${IGNORE_REGEX}"
)

echo "Writing coverage summaries"
: > "${COVERAGE_TOOL_LOG}"
"${LLVM_COV}" report \
  "${COVERAGE_ARGS[@]}" \
  --show-branch-summary \
  --show-region-summary > "${SUMMARY_FILE}" 2>> "${COVERAGE_TOOL_LOG}"

"${LLVM_COV}" export \
  "${COVERAGE_ARGS[@]}" \
  --summary-only > "${JSON_FILE}" 2>> "${COVERAGE_TOOL_LOG}"

if [[ "${GENERATE_HTML}" -eq 1 ]]; then
  echo "Writing HTML coverage report"
  "${LLVM_COV}" show \
    "${COVERAGE_ARGS[@]}" \
    -format=html \
    "-output-dir=${HTML_DIR}" \
    --show-branches=count \
    --show-line-counts-or-regions 2>> "${COVERAGE_TOOL_LOG}"
fi

GIT_COMMIT="$(git -C "${REPO_ROOT}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
GIT_BRANCH="$(git -C "${REPO_ROOT}" branch --show-current 2>/dev/null || echo unknown)"
if [[ -z "${GIT_BRANCH}" ]]; then
  GIT_BRANCH="detached"
fi
if [[ -z "$(git -C "${REPO_ROOT}" status --porcelain --untracked-files=normal 2>/dev/null)" ]]; then
  GIT_DIRTY="false"
else
  GIT_DIRTY="true"
fi
GENERATED_AT="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"

{
  printf 'generated_at=%s\n' "${GENERATED_AT}"
  printf 'git_commit=%s\n' "${GIT_COMMIT}"
  printf 'git_branch=%s\n' "${GIT_BRANCH}"
  printf 'git_dirty=%s\n' "${GIT_DIRTY}"
  printf 'gtest_filter=%s\n' "${GTEST_FILTER}"
  printf 'build_type=%s\n' "${BUILD_TYPE}"
  printf 'test_binary=%s\n' "${TEST_BINARY}"
  printf 'llvm_cov=%s\n' "${LLVM_COV}"
  printf 'llvm_profdata=%s\n' "${LLVM_PROFDATA}"
} > "${METADATA_FILE}"

TOTAL_ROW="$(awk '$1 == "TOTAL" { print; exit }' "${SUMMARY_FILE}")"
if [[ -z "${TOTAL_ROW}" ]]; then
  echo "error: LLVM coverage summary did not contain a TOTAL row" >&2
  exit 1
fi

read -r _ REGION_COUNT REGION_MISSED REGION_COVERAGE \
  FUNCTION_COUNT FUNCTION_MISSED FUNCTION_COVERAGE \
  LINE_COUNT LINE_MISSED LINE_COVERAGE \
  BRANCH_COUNT BRANCH_MISSED BRANCH_COVERAGE <<< "${TOTAL_ROW}"

REGION_COVERED=$((REGION_COUNT - REGION_MISSED))
FUNCTION_COVERED=$((FUNCTION_COUNT - FUNCTION_MISSED))
LINE_COVERED=$((LINE_COUNT - LINE_MISSED))
BRANCH_COVERED=$((BRANCH_COUNT - BRANCH_MISSED))

csv_escape() {
  local value="${1//\"/\"\"}"
  printf '"%s"' "${value}"
}

RECORDED_AS_FULL_RUN=0
if [[ "${GTEST_FILTER}" == "*" ]]; then
  if [[ ! -f "${HISTORY_FILE}" ]]; then
    printf '%s\n' \
      'generated_at,git_commit,git_branch,git_dirty,gtest_filter,region_coverage,function_coverage,line_coverage,branch_coverage,run_directory' \
      > "${HISTORY_FILE}"
  fi

  {
    csv_escape "${GENERATED_AT}"
    printf ','
    csv_escape "${GIT_COMMIT}"
    printf ','
    csv_escape "${GIT_BRANCH}"
    printf ','
    csv_escape "${GIT_DIRTY}"
    printf ','
    csv_escape "${GTEST_FILTER}"
    printf ','
    csv_escape "${REGION_COVERAGE}"
    printf ','
    csv_escape "${FUNCTION_COVERAGE}"
    printf ','
    csv_escape "${LINE_COVERAGE}"
    printf ','
    csv_escape "${BRANCH_COVERAGE}"
    printf ','
    csv_escape "${RUN_DIR}"
    printf '\n'
  } >> "${HISTORY_FILE}"

  printf '%s\n' "${RUN_DIR}" > "${LATEST_FILE}"
  cp "${SUMMARY_FILE}" "${LATEST_SUMMARY_FILE}"
  cp "${JSON_FILE}" "${LATEST_JSON_FILE}"
  RECORDED_AS_FULL_RUN=1
fi

if [[ "${UPDATE_BASELINE}" -eq 1 ]]; then
  BASELINE_TEMP="${BASELINE_FILE}.tmp.$$"
  {
    printf '{\n'
    printf '  "schema_version": 1,\n'
    printf '  "description": "MalaSim full unit-test coverage baseline",\n'
    printf '  "generated_at": "%s",\n' "${GENERATED_AT}"
    printf '  "source_commit": "%s",\n' "${GIT_COMMIT}"
    printf '  "source_branch": "%s",\n' "${GIT_BRANCH}"
    printf '  "working_tree_dirty": %s,\n' "${GIT_DIRTY}"
    printf '  "build_type": "%s",\n' "${BUILD_TYPE}"
    printf '  "coverage": {\n'
    printf '    "regions": {"covered": %s, "total": %s, "percent": %s},\n' \
      "${REGION_COVERED}" "${REGION_COUNT}" "${REGION_COVERAGE%\%}"
    printf '    "functions": {"covered": %s, "total": %s, "percent": %s},\n' \
      "${FUNCTION_COVERED}" "${FUNCTION_COUNT}" "${FUNCTION_COVERAGE%\%}"
    printf '    "lines": {"covered": %s, "total": %s, "percent": %s},\n' \
      "${LINE_COVERED}" "${LINE_COUNT}" "${LINE_COVERAGE%\%}"
    printf '    "branches": {"covered": %s, "total": %s, "percent": %s}\n' \
      "${BRANCH_COVERED}" "${BRANCH_COUNT}" "${BRANCH_COVERAGE%\%}"
    printf '  }\n'
    printf '}\n'
  } > "${BASELINE_TEMP}"
  mv "${BASELINE_TEMP}" "${BASELINE_FILE}"
fi

echo
head -n 2 "${SUMMARY_FILE}"
echo "..."
tail -n 2 "${SUMMARY_FILE}"
echo
echo "Coverage run saved to: ${RUN_DIR}"
echo "Text summary:          ${SUMMARY_FILE}"
echo "JSON summary:          ${JSON_FILE}"
echo "Test output:           ${TEST_LOG}"
if [[ -s "${COVERAGE_TOOL_LOG}" ]]; then
  echo "Coverage warnings:     ${COVERAGE_TOOL_LOG}"
fi
if [[ "${GENERATE_HTML}" -eq 1 ]]; then
  echo "HTML report:           ${HTML_DIR}/index.html"
fi
if [[ "${RECORDED_AS_FULL_RUN}" -eq 1 ]]; then
  echo "Coverage history:      ${HISTORY_FILE}"
  echo "Latest text summary:   ${LATEST_SUMMARY_FILE}"
  echo "Latest JSON summary:   ${LATEST_JSON_FILE}"
else
  echo "Filtered run: full-suite history and latest summaries were not updated."
fi
if [[ "${UPDATE_BASELINE}" -eq 1 ]]; then
  echo "Tracked baseline:       ${BASELINE_FILE}"
fi
