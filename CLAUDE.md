# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```sh
# Full build (requires VCPKG_ROOT set)
make generate BUILD_TESTS=ON
make build

# Or use the build script
./scripts/build.sh

# Run tests
make test                    # ctest in verbose mode
./build/bin/malasim_test --gtest_filter=*TestName*  # Single test

# Run with verbose logging
MALASIM_LOG_LEVEL=info ./build/bin/malasim_test

# Lint and format
make lint    # clang-tidy
make format  # clang-format

# Coverage
make generate-coverage
make coverage
```

## Project Overview

Malaria simulation (C++20) using an event-driven architecture. Key subsystems:
- **Core/Scheduler** — Time management, event scheduling via `EventManager<T>`
- **Population** — `Person`, immune system, parasite density tracking
- **Events** — Population and environment event builders
- **Spatial** — Location-based and grid-based spatial models
- **Treatment** — Drug databases, therapy strategies, treatment coverage models
- **Reporters** — Output reporting via SQLite and Specialist reporters
- **Mosquito** — Mosquito population dynamics

## Architecture

**Singleton Model**: `Model::get_instance()` is the central access point for all subsystems (config, scheduler, population, random, etc.). Components access shared state through `Model::get_*()` static methods.

**Typed IDs with Sentinels** (in `Core/types.h`):
- `core::LocationId`, `core::AgeClass`, `core::MovingLevel`, `core::PersonId`, etc.
- Invalid sentinels: `K_INVALID_LOCATION_ID`, `K_INVALID_AGE_CLASS`, etc.
- Prefer comparing against sentinels over `>= 0` checks.

**Event System**: `EventManager<T>` is templated; `Scheduler` holds `EventManager<WorldEvent>`. Events use `notify_change()` callbacks for property changes.

**Person Indexing**: `Person` inherits from multiple index handlers (`PersonIndexByLocationStateAgeClassHandler`, `PersonIndexByLocationMovingLevelHandler`). Indexes track persons by location/state/age and location/moving_level for O(1) lookup.

## Important Patterns

- **Mock factories** in `tests/fixtures/MockFactories.h` for isolated unit tests
- **TestFileGenerators.h** for integration tests that need full config — generates raster `.asc` files programmatically
- Use `test_fixtures::setup_model_with_mocks()` for mock-only tests
- RAII is preferred; avoid manual `release()` calls where possible
- Prefer `std::unique_ptr` over raw pointers; use `gsl::observer_ptr<T>` for non-owning raw pointers

## Coding Style

- C++20, Google-style with `-2` access modifier offset
- `snake_case` for functions/variables, `PascalCase` for types
- `ALL_CAPS` for constants, `kCamelCase` for enums
- Members have `_` suffix (e.g., `config_`)
- Use `[[nodiscard]]` on functions returning values that shouldn't be ignored
- Tests follow Arrange-Act-Assert; test variables: `inputX`, `mockX`, `actualX`, `expectedX`

## Debugging Segmentation Faults

```sh
# Build in Debug mode
cmake -Bbuild -DCMAKE_BUILD_TYPE=Debug -DBUILD_TESTS=ON
make -C build -j6

# Run under lldb
lldb -- ./build/bin/malasim_test --gtest_filter=*TargetTest*

# Common lldb commands
(lldb) run [args]                    # Start execution
(lldb) bt                            # Backtrace on crash
(lldb) fr va                         # Frame variables
(lldb) p variable_name               # Print value
(lldb) break set -n function_name    # Set breakpoint
(lldb) n                             # Next line
(lldb) s                             # Step in
(lldb) c                             # Continue

# For the main executable
lldb -- ./build/bin/malasim -- path/to/input.yml
```

## Known Issues (from code review)

- `Model::get_instance()` is globally accessible — makes isolated unit testing difficult
- `Person` inherits from multiple index handlers (violates single responsibility)
- `EventManager::has_event<T>()` uses O(n) dynamic_cast — potential bottleneck
- 28+ unresolved TODOs in the codebase
