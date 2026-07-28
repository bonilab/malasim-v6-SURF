# Automated Test Coverage

MalaSim uses Clang/LLVM source-based coverage. The repository script
`scripts/run_coverage.sh` configures an isolated instrumented build, runs the tests,
and saves the results so coverage work can be repeated and compared.

## Quick start

Run the complete suite and generate all reports:

```sh
make coverage
```

This is equivalent to:

```sh
./scripts/run_coverage.sh
```

The normal development build is not modified. Coverage compilation uses
`build/coverage/`, while reports are stored under `coverage-results/`.

## Saved results

Each run creates `coverage-results/runs/<timestamp>/` containing:

- `coverage-summary.txt` — human-readable totals and per-file results;
- `coverage-summary.json` — machine-readable LLVM summary data;
- `test-output.log` — complete GoogleTest output;
- `coverage-warnings.log` — diagnostics emitted while generating LLVM reports;
- `metadata.txt` — commit, branch, dirty state, filter, tools, and build type;
- `coverage.profdata` and `raw/` — merged and raw LLVM profile data;
- `html/index.html` — browsable line and branch coverage.

Two files at the output root make repeated runs easier:

- `coverage-results/latest-run.txt` contains the path of the newest full-suite run.
- `coverage-results/latest-summary.txt` and `latest-summary.json` are stable copies of
  the newest full-suite main results.
- `coverage-results/history.csv` appends the region, function, line, and branch totals
  from every full-suite run.

The output directory is ignored by Git. Copy a result elsewhere if it must be retained
outside the local workspace or attached to a review.

## Tracked baseline

The repository keeps one compact baseline in `tests/coverage-baseline.json`. It records
the covered and total region, function, line, and branch counts, plus the source commit
and build metadata. Detailed per-file reports remain local because they are large and
produce noisy diffs.

Refresh the tracked baseline only after a successful full-suite run:

```sh
make coverage-baseline
```

You can still name the saved detailed run:

```sh
make coverage-baseline COVERAGE_ARGS='--run-name coverage-after-mosquito-tests'
```

Review the baseline diff before committing it. Coverage percentages should normally
stay level or increase; if one decreases, document why the new behavior or generated
code makes the decrease intentional. A filtered run cannot update the tracked
baseline.

## Common runs

Name a run for the feature or coverage iteration:

```sh
make coverage COVERAGE_ARGS='--run-name mosquito-baseline'
```

Run only one GoogleTest suite while developing:

```sh
make coverage COVERAGE_ARGS="--run-name mosquito-unit --filter 'MosquitoTest.*'"
```

Skip HTML generation for a faster feedback run:

```sh
make coverage COVERAGE_ARGS='--run-name quick-check --no-html'
```

Choose different build or result directories:

```sh
./scripts/run_coverage.sh \
  --build-dir /tmp/malasim-coverage-build \
  --output-dir /tmp/malasim-coverage-results
```

See every option:

```sh
./scripts/run_coverage.sh --help
```

Filtered runs are useful during development, but their totals are not comparable with
full-suite totals. They are saved in their named run directory but do not update
`history.csv` or the stable `latest-*` files. Use a complete `make coverage` run for
baselines and final results.

## Prerequisites

The script requires:

- CMake;
- a Clang or AppleClang coverage build;
- `llvm-cov` and `llvm-profdata`;
- the project dependencies available through the normal CMake/vcpkg setup.

On macOS, the script locates LLVM tools through `xcrun` when they are not directly on
`PATH`. On other platforms, install LLVM and ensure `llvm-cov` and `llvm-profdata` are
on `PATH`. If Clang is not the default compiler, configure it explicitly:

```sh
CC=clang CXX=clang++ make coverage
```

Set `VCPKG_ROOT` as for a normal MalaSim build.

## Interpreting results

The primary improvement metric is line coverage, supported by branch coverage for
decision-heavy code. Function and region coverage help identify large untouched
areas, but a higher percentage alone does not demonstrate meaningful behavior.

When adding coverage:

- start from a named full-suite baseline;
- prioritize important untested behavior and error paths;
- keep deterministic unit tests separate from generated-file integration tests;
- rerun the relevant filtered suite during development;
- finish with a named full-suite run and compare its `coverage-summary.txt` or the
  corresponding rows in `history.csv`;
- update `tests/coverage-baseline.json` when the improved full-suite result is ready to
  become the new repository reference.

Do not weaken assertions or add implementation-only calls merely to execute lines.

## Troubleshooting

### No raw profile was generated

Confirm the test binary was built with `ENABLE_COVERAGE=ON` and that the selected
compiler supports `-fprofile-instr-generate` and `-fcoverage-mapping`. The script exits
instead of recording an empty result.

### LLVM reports a profile-version mismatch

`llvm-cov`, `llvm-profdata`, and the compiler should come from the same LLVM toolchain.
Remove `build/coverage/` after changing compilers, then rerun the script.

### A newly added test is missing

The script reconfigures CMake on every run, so new `*.cpp` test files are normally
discovered automatically. Check that the source is under `tests/` and matches the
current test-source glob.

### Third-party or test files appear in totals

The runner excludes `tests/`, build directories, vcpkg sources, and system/Xcode
sources. Update `IGNORE_REGEX` in `scripts/run_coverage.sh` if a new dependency path
must be excluded.
