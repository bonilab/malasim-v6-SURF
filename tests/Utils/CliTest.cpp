#include <gtest/gtest.h>

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

}  // namespace