#include <gtest/gtest.h>


#define private public
#include "Configuration/Config.h"
#undef private

#include "Simulation/Model.h"
#include "Utils/Cli.h"
#include "fixtures/TestFileGenerators.h"

class ConfigValidationEdgeCasesTest : public ::testing::Test {
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
    test_fixtures::cleanup_test_files();
  }

  Config* config() { return Model::get_config(); }
};

TEST_F(ConfigValidationEdgeCasesTest, RejectsAdditionalTimeframeAndDemographicValues) {
  config()->simulation_timeframe_.start_collect_data_day_ = -1;
  EXPECT_THROW(config()->validate_simulation_timeframe(), std::invalid_argument);

  config()->simulation_timeframe_.start_collect_data_day_ = 0;
  config()->population_demographic_.number_of_age_classes_ = 0;
  EXPECT_THROW(config()->validate_population_demographic(), std::invalid_argument);

  config()->population_demographic_.number_of_age_classes_ = 6;
  config()->population_demographic_.initial_age_structure_ = {-1};
  EXPECT_THROW(config()->validate_population_demographic(), std::invalid_argument);

  config()->population_demographic_.initial_age_structure_ = {};
  config()->population_demographic_.artificial_rescaling_of_population_size_ = 0;
  EXPECT_THROW(config()->validate_population_demographic(), std::invalid_argument);
}

TEST_F(ConfigValidationEdgeCasesTest, RejectsGenotypeAndTherapyParameterRanges) {
  config()->genotype_parameters_.mutation_mask_.clear();
  EXPECT_THROW(config()->validate_genotype_parameters(), std::invalid_argument);

  config()->genotype_parameters_.mutation_mask_ = {true};
  config()->genotype_parameters_.mutation_probability_per_locus_ = 2.0;
  EXPECT_THROW(config()->validate_genotype_parameters(), std::invalid_argument);

  config()->therapy_parameters_.tf_testing_day_ = -1;
  EXPECT_THROW(config()->validate_therapy_parameters(), std::invalid_argument);
  config()->therapy_parameters_.tf_testing_day_ = 0;
  config()->therapy_parameters_.tf_rate_ = 2.0;
  EXPECT_THROW(config()->validate_therapy_parameters(), std::invalid_argument);
}

TEST_F(ConfigValidationEdgeCasesTest, RejectsAdditionalEpidemiologicalRanges) {
  config()->epidemiological_parameters_.number_of_tracking_days_ = 0;
  EXPECT_THROW(config()->validate_epidemiological_treatment_parameters(10), std::invalid_argument);

  config()->epidemiological_parameters_.number_of_tracking_days_ = 1;
  config()->epidemiological_parameters_.days_to_clinical_under_five_ = -1;
  EXPECT_THROW(config()->validate_epidemiological_treatment_parameters(10), std::invalid_argument);

  config()->epidemiological_parameters_.days_to_clinical_under_five_ = 1;
  config()->epidemiological_parameters_.gametocyte_level_full_ = 2;
  EXPECT_THROW(config()->validate_epidemiological_treatment_parameters(10), std::invalid_argument);

  config()->epidemiological_parameters_.gametocyte_level_full_ = 0.5;
  config()->epidemiological_parameters_.relapse_duration_ = -1;
  EXPECT_THROW(config()->validate_epidemiological_transmission_parameters(), std::invalid_argument);
  config()->epidemiological_parameters_.relapse_duration_ = 1;
  config()->epidemiological_parameters_.tf_window_size_ = -1;
  EXPECT_THROW(config()->validate_epidemiological_transmission_parameters(), std::invalid_argument);
}

TEST_F(ConfigValidationEdgeCasesTest, HandlesMissingConfigurationFile) {
  Config unloaded_config;
  EXPECT_FALSE(unloaded_config.load("config-file-that-does-not-exist.yml"));
}

TEST_F(ConfigValidationEdgeCasesTest, SelectsRandomCalibrationAndHandlesEmptyCandidates) {
  config()->version6_pfpr_incidence_calibrations_.set_calibration_ids({
      {4, ImmuneSystemParametercalibration_id{}},
      {9, ImmuneSystemParametercalibration_id{}},
  });
  config()->version6_pfpr_incidence_calibrations_.set_chosen_calibration_id(-1);
  config()->select_random_immune_system_parameter_calibration_id();
  EXPECT_TRUE(config()->version6_pfpr_incidence_calibrations_.get_chosen_calibration_id() == 4 ||
              config()->version6_pfpr_incidence_calibrations_.get_chosen_calibration_id() == 9);

  config()->version6_pfpr_incidence_calibrations_.set_calibration_ids({});
  EXPECT_NO_THROW(config()->select_random_immune_system_parameter_calibration_id());
}

TEST_F(ConfigValidationEdgeCasesTest, AppliesCalibrationOverridesAndSkipsNegativeValues) {
  namespace P = ImmuneSystemOverridePaths;
  ImmuneSystemParametercalibration_id overrides;
  overrides.overrides[P::K_Z] = 0.25;
  overrides.overrides[P::K_KAPPA] = 0.75;
  overrides.overrides[P::K_MIDPOINT] = 2.0;
  overrides.overrides[P::K_P_CI_SYMP] = 0.4;
  overrides.overrides[P::K_P_SEEK_BASE] = 0.3;
  overrides.overrides[P::K_MUTATION_PROB] = -1.0;
  overrides.overrides[P::K_DEFAULT_CNV_REVERSION_MULTIPLIER] = -1.0;
  config()->version6_pfpr_incidence_calibrations_.set_chosen_calibration_id(7);
  config()->version6_pfpr_incidence_calibrations_.set_calibration_ids({{7, overrides}});

  EXPECT_NO_THROW(config()->apply_selected_immune_system_parameter_calibration_id());
  EXPECT_DOUBLE_EQ(config()->immune_system_parameters_.get_immune_effect_on_progression_to_clinical(), 0.25);
  EXPECT_DOUBLE_EQ(config()->immune_system_parameters_.get_factor_effect_age_mature_immunity(), 0.75);
  EXPECT_DOUBLE_EQ(config()->immune_system_parameters_.get_midpoint(), 2.0);
}
