#include <gtest/gtest.h>

#include <cstdio>

#include <spdlog/spdlog.h>

#include "Reporters/ValidationReporter.h"
#include "Simulation/Model.h"
#include "Utils/Cli.h"
#include "fixtures/TestFileGenerators.h"

class ValidationReporterTest : public ::testing::Test {
 protected:
  void SetUp() override {
    test_fixtures::setup_test_environment();
    utils::Cli::MaSimAppInput cli_input;
    cli_input.input_path = "test_input.yml";
    Model::set_cli_input(cli_input);
    ASSERT_TRUE(Model::get_instance()->initialize());
  }

  void TearDown() override {
    Model::get_instance()->release();
    spdlog::drop_all();
    test_fixtures::cleanup_test_files();
    for (const auto* name : {"validation_monthly_data_0.txt", "validation_summary_0.txt",
                             "validation_gene_freq_0.txt", "validation_gene_db_0.txt",
                             "validation_monthly_mutation_0.txt",
                             "validation_mosquito_res_count_0.txt"}) {
      std::remove(name);
    }
  }
};

TEST_F(ValidationReporterTest, RunsDirectReporterLifecycle) {
  ValidationReporter reporter;
  EXPECT_NO_THROW(reporter.initialize(0, "."));
  EXPECT_NO_THROW(reporter.before_run());
  EXPECT_NO_THROW(reporter.begin_time_step());
  EXPECT_NO_THROW(reporter.monthly_report());
  EXPECT_NO_THROW(reporter.after_run());
}
