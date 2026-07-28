# Testing and Test Development Guide

This is the canonical guide for adding and maintaining MalaSim tests. It is based on the
current `tests/` tree, CMake configuration, CI workflow, shared fixtures, and repository
coding rules.

Related documents:

- [`tests/README.md`](../tests/README.md) describes the fixture infrastructure and the
  migration away from external input files. Its migration counts and status notes are
  historical.
- [`docs/test_coverage.md`](test_coverage.md) contains the macOS LLVM coverage workflow.
- [`docs/test_tasks.md`](test_tasks.md) is a historical test-suite review and improvement
  backlog; it is not the current contribution standard.

## Test stack and layout

MalaSim tests use C++20, GoogleTest, and GoogleMock. All test sources are linked into the
single `malasim_test` executable and registered with CTest as one test target.

Place a test next to its domain under `tests/`, mirroring the production tree:

```text
src/Population/Person/...         -> tests/Population/Person/...Test.cpp
src/Core/Scheduler/...            -> tests/Core/Scheduler/...Test.cpp
src/Treatment/Strategies/...      -> tests/Treatment/Strategies/...Test.cpp
```

Use these naming rules for new tests:

- Name source files `<Subject>Test.cpp`.
- Name fixtures `<Subject>Test` or `<Scenario>Test`.
- Name cases for observable behavior, for example
  `TEST_F(PersonTest, SetAgeUpdatesAgeClass)`.
- Use `TEST` when no shared state is needed and `TEST_F` when setup or cleanup is shared.
- Keep reusable test-only support in `tests/fixtures/`, a domain-specific `TestHelpers.h`,
  or `tests/helpers/` when it is shared across domains.
- Do not add backup, generated, or disabled test sources to `tests/`.

`tests/CMakeLists.txt` discovers `*.cpp` files when CMake configures the build. After
adding or moving a test source, run `make generate-test` before building so the new file
is included.

## Choose the smallest test environment

Use the least expensive setup that exercises the behavior.

### 1. Plain unit test

Prefer a real object with explicit inputs when the subject has no global-model
dependency. These tests should not read files or initialize `Model`.

```cpp
TEST(DrugTypeTest, CalculateConcentrationAppliesHalfLife) {
  DrugType drug_type;
  // Arrange inputs.

  const double actual_concentration = /* Act */;

  EXPECT_NEAR(actual_concentration, /* expected */, /* justified tolerance */);
}
```

### 2. Unit test with model mocks

For code coupled to the singleton model, reuse
`tests/fixtures/MockFactories.h` and
`test_fixtures::setup_model_with_mocks()`. Keep ownership in `Model` or a
`std::unique_ptr`; pointers returned in `ModelMocks` are non-owning observation
pointers.

Release the singleton in `TearDown()` so state cannot leak into the next test:

```cpp
class MyComponentTest : public ::testing::Test {
protected:
  void SetUp() override {
    model_ = Model::get_instance();
    const auto mocks = test_fixtures::setup_model_with_mocks(model_);
    mock_config_ = mocks.config;
  }

  void TearDown() override { model_->release(); }

  Model* model_{nullptr};
  MockConfig* mock_config_{nullptr};
};
```

Before creating a new fixture, look for an existing domain base such as
`PersonTestBase`, `RandomTestBase`, or `EventManagerTestCommon`.

### 3. In-memory configuration test

For YAML encode/decode and validation behavior, use `YAML::Load` and the snippets in
`tests/fixtures/InMemoryYamlConfig.h`. Do not create a file merely to test a
`YAML::convert<T>` specialization.

Test at least:

- valid decode and the resulting typed values;
- missing required fields;
- invalid values and boundary values;
- encode/decode round trips when encoding is supported.

### 4. Generated-file integration test

Use `tests/fixtures/TestFileGenerators.h` only when the behavior genuinely requires a
complete configuration, raster/CSV input, or multi-subsystem model initialization.

```cpp
class ModelIntegrationTest : public ::testing::Test {
protected:
  void SetUp() override {
    test_fixtures::setup_test_environment(
        "test_input.yml", [](YAML::Node& config) {
          config["model_settings"]["initial_seed_number"] = 42;
        });

    utils::Cli::MaSimAppInput cli_input;
    cli_input.input_path = "test_input.yml";
    Model::set_cli_input(cli_input);

    ASSERT_TRUE(Model::get_instance()->initialize());
  }

  void TearDown() override {
    Model::get_instance()->release();
    test_fixtures::cleanup_test_files();
  }
};
```

Rules for generated-file tests:

- Generate inputs from `test_input_template.yml`; do not depend on files in
  `sample_inputs/` or on developer-machine paths.
- Change only the configuration fields required by the scenario.
- Call `cleanup_test_files()` after releasing objects that may still use those files.
- If setup can fail after files are created, make cleanup unconditional in
  `TearDown()` or an RAII helper.
- Use the existing two-location and district raster helpers when their topology matches
  the scenario.
- Generated tests currently use shared fixed filenames and run inside one test
  executable. Do not introduce concurrent access to those files.

## Test design rules

Follow Arrange-Act-Assert and keep each test focused on one behavior or one tightly
related state transition.

- Use descriptive local names such as `input_age`, `actual_class`, and
  `expected_class`.
- Verify public behavior. Inspect private implementation details only when an existing
  test seam intentionally exposes them.
- Cover the normal case, boundaries, invalid inputs, and ownership/lifecycle behavior
  relevant to the change.
- Prefer typed-ID invalid sentinels such as `core::K_INVALID_LOCATION_ID` over numeric
  assumptions such as `value >= 0`.
- Use `std::unique_ptr` for ownership and `gsl::observer_ptr<T>` or a clearly
  non-owning raw pointer where that is the production convention.
- Do not add sleeps, depend on test execution order, or leave mutable singleton/global
  state behind.
- Do not use a commented-out test, `DISABLED_` name, or TODO as a substitute for a
  runnable regression test.
- Use `GTEST_SKIP()` only for a genuinely unsupported optional environment. Missing
  data that the fixture is responsible for creating should fail setup, not silently
  skip the test.

### Assertions

- Use `ASSERT_*` for prerequisites whose failure makes later statements unsafe, such
  as a pointer that will be dereferenced or a collection size required for indexing.
- Use `EXPECT_*` for independent outcomes so one failure does not hide the rest.
- Use `EXPECT_DOUBLE_EQ` for computations expected to be exactly the same floating
  representation.
- Use `EXPECT_NEAR` with a domain-justified tolerance for numerical algorithms; avoid
  arbitrary wide tolerances.
- Use `EXPECT_THROW`/`ASSERT_THROW` with the most specific expected exception type.
- Add `SCOPED_TRACE` or assertion messages when a loop or parameterized scenario would
  otherwise be difficult to diagnose.

### Mocks

- Mock subsystem boundaries and observable collaboration, not every internal call.
- Set required interactions with `EXPECT_CALL`.
- Use `ON_CALL` for harmless defaults that are not the behavior under test.
- Use `NiceMock` only when unrelated calls are intentionally irrelevant.
- Avoid overspecifying call order unless order is part of the contract.

### Random and statistical tests

- Set an explicit seed when asserting a reproducible sequence or sampled outcome.
- For distribution properties, use the helpers in `tests/helpers/test_helpers.h`,
  a sufficient sample size, and a statistically defensible tolerance.
- Avoid assertions that can fail at an appreciable rate for a correct implementation.
- Keep statistical tests bounded so they remain suitable for the normal test suite.

## Build and run

The Make targets are the shortest local workflow:

```sh
# Configure after cloning or after adding/moving a test source.
make generate-test

# Build and run the complete CTest suite.
make test

# Run one suite or case. Quote wildcard filters so the shell does not expand them.
make gtest filter='PersonBasicTest.Age*'

# List registered GoogleTest cases.
./build/bin/malasim_test --gtest_list_tests
```

Useful direct commands:

```sh
cd build
GTEST_COLOR=1 ctest -V

cd bin
./malasim_test --gtest_filter='PersonBasicTest.AgeSetInitialClass'
```

The test target’s CTest working directory is `build/bin`; generated fixture files will
normally appear there. CI uses the `ci` configure preset followed by the
`ci-release` build and test presets. To reproduce that structure locally with
`VCPKG_ROOT` set:

```sh
cmake --preset local
cmake --build --preset local-release
ctest --preset local-release
```

For a repeatable coverage run with saved text, JSON, HTML, test-log, and history
artifacts, use `make coverage`. Options and result layout are documented in
[`docs/test_coverage.md`](test_coverage.md).

The committed reference totals live in `tests/coverage-baseline.json`. Refresh them
after a successful full-suite improvement with `make coverage-baseline`, then review
the baseline diff before committing it.

## Before submitting

For a change that adds or modifies tests:

1. Reconfigure if a source file was added, removed, or moved.
2. Run the narrow filter while developing.
3. Format the touched C++ files with the repository `.clang-format`.
4. Run `make test`.
5. Run `make lint` when practical; warnings are configured as errors.
6. Confirm the new test fails without the production fix when it is a regression test.
7. Confirm no generated `test_*`, database, profile, or coverage artifacts are staged.
8. Update shared fixtures or this guide when introducing a reusable testing pattern.
