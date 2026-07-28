#include <gtest/gtest.h>

#include <cstdio>
#include <sstream>

#include "Reporters/Utility/ReporterUtils.h"
#include "Simulation/Model.h"
#include "Utils/Index/PersonIndexByLocationStateAgeClass.h"
#include "Utils/Cli.h"
#include "fixtures/TestFileGenerators.h"

TEST(ReporterUtilsTest, EmptyPopulationProducesZeroFrequencyRecords) {
  PersonIndexByLocationStateAgeClass person_index(1, Person::NUMBER_OF_STATE, 1);
  std::stringstream genotype_frequency_1;
  std::stringstream genotype_frequency_2;
  std::stringstream genotype_frequency_3;
  std::stringstream combined_frequency;

  ReporterUtils::output_genotype_frequency1(genotype_frequency_1, 2, &person_index);
  ReporterUtils::output_genotype_frequency2(genotype_frequency_2, 2, &person_index);
  ReporterUtils::output_genotype_frequency3(genotype_frequency_3, 2, &person_index);
  ReporterUtils::output_3_genotype_frequency(combined_frequency, 2, &person_index);

  EXPECT_NE(genotype_frequency_1.str().find("0"), std::string::npos);
  EXPECT_NE(genotype_frequency_2.str().find("0"), std::string::npos);
  EXPECT_NE(genotype_frequency_3.str().find("0"), std::string::npos);
  EXPECT_NE(combined_frequency.str().find("0"), std::string::npos);
}

class ReporterUtilsModelTest : public ::testing::Test {
protected:
  void SetUp() override {
    test_fixtures::setup_test_environment();
    utils::Cli::MaSimAppInput cli_input;
    cli_input.input_path = "test_input.yml";
    cli_input.job_number = 0;
    Model::set_cli_input(cli_input);
    ASSERT_TRUE(Model::get_instance()->initialize());
    ReporterUtils::initialize_moi_file_logger();
  }

  void TearDown() override {
    Model::get_instance()->release();
    std::remove("moi_0.txt");
    test_fixtures::cleanup_test_files();
  }
};

TEST_F(ReporterUtilsModelTest, EmptyPopulationSupportsMosquitoAndMoiReports) {
  PersonIndexByLocationStateAgeClass person_index(1, Person::NUMBER_OF_STATE, 1);
  std::stringstream genotype_frequency_4;
  std::stringstream prmc_frequency;
  std::stringstream moi;

  EXPECT_NO_THROW(ReporterUtils::output_genotype_frequency4(
      genotype_frequency_4, prmc_frequency, 1, &person_index));
  EXPECT_NO_THROW(ReporterUtils::output_moi(moi, &person_index));
  EXPECT_TRUE(moi.str().empty());
}
