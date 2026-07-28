#include <gtest/gtest.h>

#include "Simulation/Model.h"
#include "Utils/Cli.h"
#include "fixtures/TestFileGenerators.h"

class ModelCalibrationVerificationTest : public ::testing::Test {
protected:
  void TearDown() override {
    Model::get_instance()->release();
    test_fixtures::cleanup_test_files();
  }
};

TEST_F(ModelCalibrationVerificationTest, BeforeRunVerifiesAllCalibrationOverridePaths) {
  test_fixtures::setup_test_environment("test_input.yml", [](YAML::Node &config) {
    config["version6_pfpr_incidence_calibrations"] = YAML::Load(R"(
      chosen_calibration_id: 0
      random_selection: false
      calibration_ids:
        0:
          "immune_system_parameters:immune_effect_on_progression_to_clinical": 1.2
          "immune_system_parameters:factor_effect_age_mature_immunity": 0.4
          "immune_system_parameters:midpoint": 0.8
          "epidemiological_parameters:allow_new_coinfection_to_cause_symptoms:probability": 0.2
          "epidemiological_parameters:age_based_probability_of_seeking_treatment:power:base": 0.7
          "genotype_parameters:mutation_probability_per_locus": 0.0005
          "genotype_parameters:default_cnv_reversion_multiplier": 0.2
          "unknown:field": 7.0
    )");
  });

  utils::Cli::MaSimAppInput cli_input;
  cli_input.input_path = "test_input.yml";
  cli_input.reporter = "Console";
  Model::set_cli_input(cli_input);
  ASSERT_TRUE(Model::get_instance()->initialize());
  EXPECT_NO_THROW(Model::get_instance()->before_run());
}
