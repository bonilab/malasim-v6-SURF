#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>

#include "Utils/Cli.h"

namespace {

class CliParseTest : public ::testing::Test {
protected:
  void TearDown() override {
    // Reset CLI state between tests
  }
};

TEST_F(CliParseTest, ParseDefaultValues) {
  const char* argv[] = {"malasim"};
  int argc = 1;

  auto cli_input = utils::Cli::parse_args(argc, const_cast<char**>(argv));

  EXPECT_EQ(cli_input.input_path, "input.yml");
  EXPECT_EQ(cli_input.output_path, "");
  EXPECT_EQ(cli_input.reporter, "");
  EXPECT_EQ(cli_input.verbosity, 0);
  EXPECT_EQ(cli_input.job_number, 0);
  EXPECT_EQ(cli_input.replicate, 1);
  EXPECT_FALSE(cli_input.dump_movement_matrix);
  EXPECT_FALSE(cli_input.record_individual_movement);
  EXPECT_FALSE(cli_input.record_cell_movement);
  EXPECT_FALSE(cli_input.record_district_movement);
}

TEST_F(CliParseTest, ParseCustomInputPath) {
  const char* argv[] = {"malasim", "-i", "/custom/path.yml"};
  int argc = 3;

  auto cli_input = utils::Cli::parse_args(argc, const_cast<char**>(argv));

  EXPECT_EQ(cli_input.input_path, "/custom/path.yml");
}

TEST_F(CliParseTest, ParseCustomOutputPath) {
  const char* argv[] = {"malasim", "-o", "/output/dir"};
  int argc = 3;

  auto cli_input = utils::Cli::parse_args(argc, const_cast<char**>(argv));

  EXPECT_EQ(cli_input.output_path, "/output/dir");
}

TEST_F(CliParseTest, ParseVerbosity) {
  const char* argv[] = {"malasim", "-v", "2"};
  int argc = 3;

  auto cli_input = utils::Cli::parse_args(argc, const_cast<char**>(argv));

  EXPECT_EQ(cli_input.verbosity, 2);
}

TEST_F(CliParseTest, ParseMultipleOptions) {
  const char* argv[] = {"malasim", "-i", "/input.yml", "-o", "/output", "-v", "1", "-j", "42"};
  int argc = 9;

  auto cli_input = utils::Cli::parse_args(argc, const_cast<char**>(argv));

  EXPECT_EQ(cli_input.input_path, "/input.yml");
  EXPECT_EQ(cli_input.output_path, "/output");
  EXPECT_EQ(cli_input.verbosity, 1);
  EXPECT_EQ(cli_input.job_number, 42);
}

TEST_F(CliParseTest, ParseMovementOptions) {
  const char* argv[] = {"malasim", "--im", "--mc"};
  int argc = 3;

  auto cli_input = utils::Cli::parse_args(argc, const_cast<char**>(argv));

  EXPECT_TRUE(cli_input.record_individual_movement);
  EXPECT_TRUE(cli_input.record_cell_movement);
  EXPECT_FALSE(cli_input.record_district_movement);
}

TEST_F(CliParseTest, ParseDistrictMovement) {
  const char* argv[] = {"malasim", "--md"};
  int argc = 2;

  auto cli_input = utils::Cli::parse_args(argc, const_cast<char**>(argv));

  EXPECT_TRUE(cli_input.record_district_movement);
}

TEST_F(CliParseTest, ParseMemoryStatsFlag) {
  const char* argv[] = {"malasim", "--memory-stats"};
  int argc = 2;

  auto cli_input = utils::Cli::parse_args(argc, const_cast<char**>(argv));

  EXPECT_TRUE(cli_input.print_memory_stats);
}

TEST_F(CliParseTest, MemoryStatsFlagDefaultFalse) {
  const char* argv[] = {"malasim"};
  int argc = 1;

  auto cli_input = utils::Cli::parse_args(argc, const_cast<char**>(argv));

  EXPECT_FALSE(cli_input.print_memory_stats);
}

// validate_config exercises the input-file existence check, the mutual-exclusion
// rule for movement flags, the aggregate movement flag, and the verbosity
// switch. None of these touch the Model singleton.
class CliValidateTest : public ::testing::Test {
protected:
  static constexpr const char* kTempInput = "cli_validate_input.yml";

  void create_temp_input() {
    std::ofstream out(kTempInput);
    out << "# temp\n";
  }

  void TearDown() override { std::remove(kTempInput); }
};

TEST_F(CliValidateTest, FailsWhenInputFileMissing) {
  utils::Cli::MaSimAppInput input;
  input.input_path = "definitely_missing_cli_input.yml";
  EXPECT_FALSE(utils::Cli::validate_config(input));
}

TEST_F(CliValidateTest, RejectsCellAndDistrictMovementTogether) {
  create_temp_input();
  utils::Cli::MaSimAppInput input;
  input.input_path = kTempInput;
  input.record_cell_movement = true;
  input.record_district_movement = true;
  EXPECT_FALSE(utils::Cli::validate_config(input));
}

TEST_F(CliValidateTest, AggregatesMovementFlagFromIndividual) {
  create_temp_input();
  utils::Cli::MaSimAppInput input;
  input.input_path = kTempInput;
  input.record_individual_movement = true;
  ASSERT_TRUE(utils::Cli::validate_config(input));
  EXPECT_TRUE(input.record_movement);
}

TEST_F(CliValidateTest, NoMovementLeavesAggregateFalse) {
  create_temp_input();
  utils::Cli::MaSimAppInput input;
  input.input_path = kTempInput;
  ASSERT_TRUE(utils::Cli::validate_config(input));
  EXPECT_FALSE(input.record_movement);
}

TEST_F(CliValidateTest, AcceptsEachVerbosityBranch) {
  create_temp_input();
  for (int verbosity : {0, 1, 2, 99}) {
    utils::Cli::MaSimAppInput input;
    input.input_path = kTempInput;
    input.verbosity = verbosity;
    EXPECT_TRUE(utils::Cli::validate_config(input)) << "verbosity=" << verbosity;
  }
}

TEST_F(CliValidateTest, CreateCliOptionsParsesFullCommandLine) {
  CLI::App app;
  utils::Cli::MaSimAppInput input;
  utils::Cli::create_cli_options(app, input);

  ASSERT_NO_THROW(
      app.parse(std::string{"-i in.yml -o /tmp/o -r MMC -v 1 -j 4 --replicate 2"}));

  EXPECT_EQ(input.input_path, "in.yml");
  EXPECT_EQ(input.output_path, "/tmp/o");
  EXPECT_EQ(input.reporter, "MMC");
  EXPECT_EQ(input.verbosity, 1);
  EXPECT_EQ(input.job_number, 4);
  EXPECT_EQ(input.replicate, 2);
}

}  // namespace
